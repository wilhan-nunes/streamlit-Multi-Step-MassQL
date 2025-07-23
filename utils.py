import os
import subprocess
import uuid
import logging
from typing import List

from streamlit import cache_data
from dataclasses import dataclass

import pandas as pd
import yaml
from gnpsdata import workflow_fbmn

import massql_launch
import streamlit as st

logging.basicConfig(
    level=logging.DEBUG,
    format="%(asctime)s - %(levelname)s - %(message)s",
    stream=logging.sys.stdout,
)


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

@cache_data
def download_and_filter_mgf(task_id: str) -> (str, str):
    os.makedirs("temp_mgf", exist_ok=True)
    unique_uuid = str(uuid.uuid4())
    mgf_file_path = f"temp_mgf/{unique_uuid}_mgf_all.mgf"

    logging.info("Downloading mgf...")
    workflow_fbmn.download_mgf(task_id, mgf_file_path)
    logging.info(f"MGF saved to {mgf_file_path}")

    logging.info("Starting MGF filtering...")
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
    # Save the cleaned MGF file
    cleaned_mgf = f"temp_mgf/{unique_uuid}_mgf_cleaned.mgf"
    with open(cleaned_mgf, "w") as fout:
        fout.writelines(cleaned_mgf_lines)
    logging.info(f"Cleaned MGF saved to {cleaned_mgf}")

    # Extract all scan numbers from the cleaned MGF file
    with open(cleaned_mgf, "r") as mgf_file:
        for line in mgf_file:
            if line.startswith("SCANS="):
                scans_list.append(line.strip().split("=")[1])

    return cleaned_mgf, scans_list


@cache_data
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
    logging.info(f"Total Scans: {total_scans} ** Kept: {kept_scans} scans ** Excluded: {total_scans - kept_scans}")
    return output_mgf_path


def add_df_and_filtering(df, key_prefix:str, default_cols: List = None) -> pd.DataFrame:
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
                f"Column {i+1}", df.columns, key=f"{key_prefix}_col_select_{i}"
            )
        with col2:
            search_term = st.text_input(
                f"Contains (Column {i+1})", key=f"{key_prefix}_search_input_{i}"
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

    stage1_passed_mgf = filter_mgf_by_scans(mgf_path, f"temp_mgf/{uuid.uuid4()}_scans_passed_stg1.mgf",
                                            scans_to_keep)
