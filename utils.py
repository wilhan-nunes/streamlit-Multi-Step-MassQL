import os
from dataclasses import dataclass
from io import StringIO
from typing import List

import pandas as pd
import streamlit as st
import yaml
from gnpsdata import workflow_fbmn, taskinfo, taskresult

import massql_launch


def get_git_short_rev():
    try:
        with open('.git/logs/HEAD', 'r') as f:
            last_line = f.readlines()[-1]
            hash_val = last_line.split()[1]
        return hash_val[:7]
    except Exception:
        return ".git/ not found"


@dataclass
class MassQLQueries:
    with open('massql_queries.yaml', 'r') as file:
        data = yaml.safe_load(file)
    ALL_MASSQL_QUERIES = data['ALL_MASSQL_QUERIES']
    stage1 = {key: value for key, value in ALL_MASSQL_QUERIES.items() if "stage1" in key.lower()}
    stage2 = {key: value for key, value in ALL_MASSQL_QUERIES.items() if "stage2" in key.lower()}
    mono_queries = {key: value for key, value in ALL_MASSQL_QUERIES.items() if "mono" in key.lower()}
    di_queries = {key: value for key, value in ALL_MASSQL_QUERIES.items() if "di" in key.lower()}
    tri_queries = {key: value for key, value in ALL_MASSQL_QUERIES.items() if "tri" in key.lower()}


with open('bile_acid_tree.yaml', 'r') as file:
    bile_acid_tree = yaml.safe_load(file)


@st.cache_data
def download_and_filter_mgf(task_id: str) -> (str, str):
    os.makedirs("temp_mgf", exist_ok=True)
    mgf_file_path = f"temp_mgf/{task_id}_mgf_all.mgf"
    cleaned_mgf = f"temp_mgf/{task_id}_mgf_cleaned.mgf"

    # Skip if cleaned file already exists
    if os.path.exists(cleaned_mgf):
        print(f"Skipping download, using existing file: {cleaned_mgf}")
        scans_list = []
        with open(cleaned_mgf, "r") as mgf_file:
            for line in mgf_file:
                if line.startswith("SCANS="):
                    scans_list.append(line.strip().split("=")[1])
        return cleaned_mgf, scans_list

    print("Downloading mgf...")
    task_info = taskinfo.get_task_information(task_id)
    workflowname = task_info.get('workflowname')
    if workflowname == 'feature_based_molecular_networking_workflow':
        fbmn_download_mgf_wrapper(mgf_file_path, task_id)
    elif workflowname == 'classical_networking_workflow':
        gnps2_download_resultfile_wrapper(mgf_file_path, task_id)
    else:
        print(f"Unsupported workflow: {workflowname}. Cannot download MGF.")
        raise ValueError(f"Unsupported workflow: {workflowname}. Cannot download MGF.")

    print(f"MGF saved to {mgf_file_path}")

    print("Starting MGF filtering...")
    ## Extract all scan numbers from the MGF file
    scans_list = []
    with open(mgf_file_path, "r") as mgf_file:
        lines = mgf_file.readlines()
    cleaned_mgf_lines = []
    inside_scan = False
    current_scan = []
    for line in lines:
        if line.startswith("BEGIN IONS"):
            inside_scan = True
            current_scan = [line]  # Start a new scan block
        elif line.startswith("END IONS"):
            current_scan.append(line)
            if any(
                    len(peak.split()) == 2
                    and all(part.replace(".", "", 1).isdigit() for part in peak.split())
                    for peak in current_scan
            ):
                cleaned_mgf_lines.extend(current_scan)
            inside_scan = False
        elif inside_scan:
            current_scan.append(line)
        else:
            cleaned_mgf_lines.append(line)

    with open(cleaned_mgf, "w") as fout:
        fout.writelines(cleaned_mgf_lines)
    print(f"Cleaned MGF saved to {cleaned_mgf}")

    # Extract all scan numbers from the cleaned MGF file
    with open(cleaned_mgf, "r") as mgf_file:
        for line in mgf_file:
            if line.startswith("SCANS="):
                scans_list.append(line.strip().split("=")[1])

    return cleaned_mgf, scans_list


@st.cache_data
def gnps2_download_resultfile_wrapper(mgf_file_path, task_id):
    return taskresult.download_gnps2_task_resultfile(task_id, "nf_output/clustering/specs_ms.mgf", mgf_file_path)


@st.cache_data
def gnps2_get_library_match_dataframe(task_id):
    return taskresult.get_gnps2_task_resultfile_dataframe(task_id, "nf_output/library/merged_results_with_gnps.tsv")


@st.cache_data
def fbmn_download_mgf_wrapper(mgf_file_path, task_id):
    return workflow_fbmn.download_mgf(task_id, mgf_file_path)


@st.cache_data
def filter_mgf_by_scans(input_mgf_path, output_mgf_path, scans_to_keep):
    """
    Write a new MGF file containing only the scans in scans_to_keep.
    :param input_mgf_path: Path to the input MGF file.
    :param output_mgf_path: Path to the output filtered MGF file.
    :param scans_to_keep: List of scan numbers (as strings or ints) to keep.
    """
    scans_to_keep = set(str(s) for s in scans_to_keep)
    total_scans = 0
    kept_scans = 0
    with open(input_mgf_path, "r") as infile, open(output_mgf_path, "w") as outfile:
        write_block = False
        block_lines = []
        for line in infile:
            if line.strip() == "BEGIN IONS":
                block_lines = [line]
                write_block = False
            elif line.startswith("SCANS="):
                total_scans += 1
                scan_num = line.strip().split("=")[1]
                if scan_num in scans_to_keep:
                    write_block = True
                block_lines.append(line)
            elif line.strip() == "END IONS":
                block_lines.append(line)
                if write_block:
                    outfile.writelines(block_lines)
                    kept_scans += 1
            else:
                block_lines.append(line)
    print(f"Total Scans: {total_scans} ** Kept: {kept_scans} scans ** Excluded: {total_scans - kept_scans}")
    return output_mgf_path


def add_df_and_filtering(df, key_prefix: str, default_cols: List = None) -> pd.DataFrame:
    # Session state for tracking number of filters
    if f"{key_prefix}_filter_count" not in st.session_state:
        st.session_state[f"{key_prefix}_filter_count"] = 1

    add_col, remove_col, _, _ = st.columns(4)
    with add_col:
        # Button to add more filter fields
        if st.button("➕ Add Filter Field", use_container_width=True, key=f"{key_prefix}_add_btn"):
            st.session_state[f"{key_prefix}_filter_count"] += 1
    with remove_col:
        if st.button("➖ Remove Filter Field", use_container_width=True, key=f"{key_prefix}_rmv_btn"):
            st.session_state[f"{key_prefix}_filter_count"] -= 1

    filtered_df = df.copy()
    cols = st.columns([1, 2])  # for headers
    cols[0].markdown("**Filter Column**")
    cols[1].markdown("**Search String**")

    # Generate filter fields
    for i in range(st.session_state[f"{key_prefix}_filter_count"]):
        col1, col2 = st.columns([1, 2])
        with col1:
            selected_col = st.selectbox(
                f"Column {i + 1}", df.columns, key=f"{key_prefix}_col_select_{i}"
            )
        with col2:
            search_term = st.text_input(
                f"Contains (Column {i + 1})", key=f"{key_prefix}_search_input_{i}"
            )

        if selected_col and search_term:
            filtered_df = filtered_df[filtered_df[selected_col].str.contains(search_term, case=False, na=False)]

    # Show result
    st.markdown("### 🔎 Filtered Results")
    st.write(f"Total results: {len(filtered_df)}")
    all_cols = df.columns
    if default_cols:
        with st.expander('Cols to show'):
            cols_to_show = st.multiselect("Columns to show", options=all_cols, default=default_cols,
                                          label_visibility='collapsed')
    else:
        cols_to_show = all_cols

    return filtered_df[cols_to_show]


def highlight_hydroxy(s):
    styles = []
    for v in s:
        v_str = str(v)
        if 'Trihydroxy' in v_str:
            styles.append('color: #4287f5')  # blue
        elif 'Dihydroxy' in v_str:
            styles.append('color: #ae0775')  # purple
        elif 'Monohydroxy' in v_str:
            styles.append('color: #18b760')  # green
        else:
            styles.append('')
    return styles


def insert_mgf_info(task: str, input_mgf: str, validation_df: pd.DataFrame) -> StringIO:
    print(f"Inserting MGF info for task {task}...")

    mask = validation_df["query_validation"] != "Did not pass stage1 filtering"
    valid_scans = set(
        pd.to_numeric(validation_df.loc[mask, "#Scan#"], errors="coerce")
        .dropna().astype(int).tolist()
    )
    scan_to_validation = {
        int(k): v for k, v in zip(
            pd.to_numeric(validation_df["#Scan#"], errors="coerce").fillna(-1).astype(int),
            validation_df["classification"]
        ) if k != -1
    }

    buffer = StringIO()
    spectrum_lines = []
    skip_spectrum = False
    print(f"Processing MGF file: {input_mgf}")
    print(f"Filtering to {len(valid_scans)} scans that passed validation (out of {len(validation_df)} total scans)")

    file_contents = open(input_mgf, "r").readlines()
    for line in file_contents:
        if line.startswith("BEGIN IONS"):
            spectrum_lines = [line]
            skip_spectrum = False
        elif line.startswith("SCANS"):
            scan_number = int(line.split("=")[1].strip())
            spectrum_lines.append(line)

            if scan_number not in valid_scans:
                skip_spectrum = True
                continue

            validation_status = scan_to_validation.get(scan_number, "Unknown")

            insert_string = f"MASSQL_VALIDATION={validation_status}\n"

            for prev_line in spectrum_lines[:-1]:
                buffer.write(prev_line)
            buffer.write(insert_string)
            buffer.write(line)
            spectrum_lines = []

        elif line.startswith("END IONS"):
            if not skip_spectrum:
                spectrum_lines.append(line)
                for spectrum_line in spectrum_lines:
                    buffer.write(spectrum_line)
            spectrum_lines = []
        else:
            if not skip_spectrum:
                if spectrum_lines:
                    spectrum_lines.append(line)
                else:
                    buffer.write(line)
    print(f"Processed {input_mgf}")
    buffer.seek(0)
    return buffer


if __name__ == "__main__":
    task_id = "4e5f76ebc4c6481aba4461356f20bc35"
    cleaned_mgf, scans_list = download_and_filter_mgf(task_id)
    mgf_path = cleaned_mgf

    with open('massql_queries.yaml', 'r') as file:
        data = yaml.safe_load(file)
        ALL_MASSQL_QUERIES = data['ALL_MASSQL_QUERIES']

    only_stage1 = {key: value for key, value in ALL_MASSQL_QUERIES.items() if "stage1" in key}

    stage1_all_results = massql_launch.run_massql(mgf_path, queries_dict=only_stage1)
    stage1_results_df = pd.DataFrame(stage1_all_results)
    scans_to_keep = set(sum(stage1_results_df['scan_list'], []))

    stage1_passed_mgf = filter_mgf_by_scans(mgf_path, f"temp_mgf/{task_id}_scans_passed_stg1.mgf",
                                            scans_to_keep)


def insert_plot_download_button(identifier, svg_data, key_prefix):
    col1, col2, _, _ = st.columns(4)
    col1.download_button(
        label=":material/download: Download Plot as SVG",
        data=svg_data,
        file_name=f"{identifier}.svg",
        mime="image/svg+xml",  # Set the MIME type to SVG
        key=f'{key_prefix}_svg_download',
        type="primary"
    )
