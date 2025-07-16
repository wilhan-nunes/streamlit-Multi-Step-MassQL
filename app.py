import ast
import glob
import os
import uuid
from typing import List

import gnpsdata
import pandas as pd
from gnpsdata import workflow_fbmn

import massql_launch
from utils import (
    download_and_filter_mgf,
    filter_mgf_by_scans,
    MassQLQueries,
    bile_acid_tree,
    add_df_and_filtering
)
from tree_plotter import create_custom_tree
from tree_classifier import check_classification_paths
import streamlit as st


massql_queries = MassQLQueries()
ALL_MASSQL_QUERIES = massql_queries.ALL_MASSQL_QUERIES
stage1 = massql_queries.stage1
stage2 = massql_queries.stage2
mono_queries = massql_queries.mono_queries
di_queries = massql_queries.di_queries
tri_queries = massql_queries.tri_queries

# Set page configuration
page_title = "Multi-Step MassQL Bile Acid Isomer Annotation"

# TODO: Bump version
app_version = "2025-07-16"

st.set_page_config(
    page_title=page_title,
    layout="wide",
    # page_icon=":pill:",
    menu_items={"About": ("**App version: %s**" % app_version)},
)


def process_results(massql_results_df: List, library_matches: pd.DataFrame, all_scans: List[str]):
    """
    Process results and include scans without library matches in full_table output.

    Args:
        massql_results_df: List of MassQL results
        library_matches: DataFrame with library matches
        all_scans: List of all scan numbers as strings
    """

    # Process MassQL results
    all_query_results_df = pd.DataFrame(massql_results_df)
    all_query_results_df["scan_list"] = all_query_results_df["scan_list"].apply(
        lambda x: ast.literal_eval(x) if isinstance(x, str) else x
    )
    all_query_results_df = all_query_results_df.explode("scan_list")
    all_query_results_df = all_query_results_df.rename(
        columns={"scan_list": "#Scan#", "query": "query_validation"}
    )

    with st.spinner("Merging results..."):
        # Ensure consistent data types
        all_query_results_df["#Scan#"] = all_query_results_df["#Scan#"].astype(str)
        library_matches["#Scan#"] = library_matches["#Scan#"].astype(str)

        # Create complete scan list DataFrame
        all_scans_df = pd.DataFrame({"#Scan#": [str(scan) for scan in all_scans]})

        # Merge everything: all_scans -> library_matches -> query_results
        full_table = (all_scans_df
                      .merge(library_matches, on="#Scan#", how="left")
                      .merge(all_query_results_df, on="#Scan#", how="left"))

        # Fill missing values
        full_table["query_validation"] = full_table["query_validation"].fillna("Did not pass stage1 filtering")

        # Library matches only (existing functionality)
        library_matches_only = library_matches.merge(all_query_results_df, on="#Scan#", how="left")

        # Reorder columns and aggregate
        cols = ["query_validation", "Compound_Name"] + [col for col in full_table.columns
                                                        if col not in ["query_validation", "Compound_Name"]]
        full_table = full_table[cols]

        # Group by scan and aggregate
        full_table = full_table.groupby("#Scan#", as_index=False).agg({
            "query_validation": lambda x: ";".join(set(x.dropna())),
            **{col: "first" for col in full_table.columns
               if col not in ["#Scan#", "query_validation"]}
        })

    return library_matches_only, full_table, all_query_results_df


def cleanup_massql_files():
    feather_files = glob.glob("temp_mgf/*.feather")
    for file in feather_files:
        try:
            os.remove(file)
        except Exception as e:
            st.warning(f"Could not delete {file}: {e}")


def load_example_data():
    st.session_state.library_matches = pd.read_csv(
        "examples/example_library_matches.csv", dtype=str
    )
    st.session_state.all_query_results_df = pd.read_csv(
        "examples/example_all_results.csv", dtype=str
    )
    st.session_state.full_table = pd.read_csv(
        "examples/example_lib_and_query_results.csv", dtype=str
    ).fillna("No match")
    with open("examples/example_massql_results_after_stg1.txt", "r") as f:
        st.session_state.massql_results_df = ast.literal_eval(f.read())
    with open("examples/example_all_scans.txt", "r") as f:
        st.session_state.all_scans = [line.strip() for line in f]



def get_bile_acids_classifications(results_df, exclude_string:str):
    passed_queries = results_df[
        ~results_df["query_validation"].str.contains(exclude_string, case=False)
    ]
    passed_queries["classification"] = passed_queries["query_validation"].apply(
        lambda x: check_classification_paths(str(x).split(";"), bile_acid_tree)[
            "satisfied_paths"
        ]
    )
    filtered_classifications = passed_queries[
        passed_queries["classification"].apply(lambda x: bool(x))
    ]

    return filtered_classifications


with st.sidebar:
    st.subheader("Analysis configuration")
    load_example = st.checkbox("Load query example", value=False, key="load_example_checkbox")
    if load_example:
        task_id_value = "4e5f76ebc4c6481aba4461356f20bc35"
    else:
        task_id_value = ""
    task_id = st.text_input("FBMN task ID:", value=task_id_value)

    col1, col2 = st.columns(2)
    with col1:
        run_query = st.button("Run Query", key="run_query", use_container_width=True)

    with col2:
        if st.button(
            "Restart Session",
            icon="♻️",
            key="restart_session",
            type="primary",
        ):
            # Reset the session state
            st.session_state.clear()
            st.session_state.load_example_checkbox = False
            st.rerun()

    st.subheader("Contributors")
    st.markdown(
        """
    - [Ipsita Mohanty PhD](https://scholar.google.com/citations?user=iHJ3vgsAAAAJ) - UC San Diego
    - [Wilhan Nunes PhD](https://scholar.google.com/citations?user=4cPVoeIAAAAJ) - UC San Diego
    - [Helena Russo PhD](https://sites.google.com/view/helenamrusso/home) - UC San Diego
    - [Mingxun Wang PhD](https://www.cs.ucr.edu/~mingxunw/) - UC Riverside
    """
    )

    st.subheader("Documentations and Resources")
    st.markdown("""
    [Feature Based Molecular Networking](https://wang-bioinformatics-lab.github.io/GNPS2_Documentation/fbmn/)
    """
    )

if not run_query and "run_query_done" not in st.session_state:
    from welcome import welcome_page

    welcome_page()

if run_query:
    st.session_state["run_query_done"] = True
    if not load_example:
        with st.spinner("Downloading files..."):
            library_matches = workflow_fbmn.get_library_match_dataframe(task_id)
            cleaned_mgf_path, all_mgf_scans = download_and_filter_mgf(task_id)
            mgf_path = cleaned_mgf_path

        with st.spinner("Running Stage 1 queries..."):
            stage1_all_results = massql_launch.run_massql(
                mgf_path, queries_dict=stage1
            )
            stage1_results_df = pd.DataFrame(stage1_all_results)
            # create a new mgf filtering to just maintain the scans that passed stage1
            scans_to_keep = set(sum(stage1_results_df["scan_list"], []))
            stage1_passed_mgf = filter_mgf_by_scans(
                mgf_path,
                f"temp_mgf/{task_id}_stg1_passed.mgf",
                scans_to_keep,
            )

        with st.spinner(
            "Running MassQL for filtered scans... This may take a while, please be patient!"
        ):
            container = st.empty()
            with container:
                # Run all queries for the filtered data
                massql_results_df = massql_launch.run_massql(
                    stage1_passed_mgf, ALL_MASSQL_QUERIES
                )

        cleanup_massql_files()

    else:
        # this function stores the static result file dataframes in st.session_state, just as above.
        load_example_data()
        all_mgf_scans = st.session_state.get('all_scans')
        massql_results_df =st.session_state.get('massql_results_df')
        library_matches =st.session_state.get('library_matches')

    with st.spinner("Processing tables..."):
        (
            only_library_matches,
            full_table,
            all_query_results_df,
        ) = process_results(massql_results_df, library_matches, all_mgf_scans)

        # #TODO:remove later (saving example files)
        # pd.DataFrame(massql_results_df).to_csv('examples/example_massql_results_after_stg1.csv', index=False)
        # full_table.to_csv('example_lib_and_query_results.csv', index=False)
        # only_library_matches.to_csv('example_library_matches.csv', index=False)
        # all_query_results_df.to_csv('example_all_results.csv', index=False)


        # Store results in session state
        st.session_state["only_library_matches"] = only_library_matches
        st.session_state["full_table"] = full_table
        st.session_state["all_query_results_df"] = all_query_results_df


if run_query or st.session_state.get("run_query_done"):
    st.title("Multi-step MassQL Results")
    only_library_matches = st.session_state["only_library_matches"]
    full_table = st.session_state["full_table"]
    full_table["Compound_Name"] = full_table[
        "Compound_Name"
    ].fillna("No match")

    filtered_classifications = get_bile_acids_classifications(full_table, exclude_string="did not pass")
    if len(filtered_classifications) > 0:
        feature_ids_dict = filtered_classifications[["#Scan#", "Compound_Name"]].astype(str)
        feature_ids_dict = feature_ids_dict.set_index("#Scan#")["Compound_Name"].to_dict()
        feature_ids_dict = dict(sorted(feature_ids_dict.items(), key=lambda item: item[1]))
    else:
        st.warning("No classifications retrieved for this task ID. Inspect the full table below for details")
        st.write(full_table)
        st.stop()

    viz_tab, class_tab, lib_tab, full_tab = st.tabs(
        ["👓 Visualizations", "🗂️ Classified", "📚 Library Matches", "📋 Full Table",]
    )

    default_cols = ["#Scan#", "Compound_Name", 'classification', 'query_validation']

    with viz_tab:
        st.subheader("Feature Classification")
        selected_feature = st.selectbox(
            f"Select a feature : :blue-badge[{len(feature_ids_dict)} of {len(full_table)}]",
            [f"{v}: {k}" for v, k in feature_ids_dict.items()],
            index=0,
        )
        fid = selected_feature.split(":")[0]

        validation_lists = filtered_classifications[
            filtered_classifications["#Scan#"] == fid
        ]["classification"].values[0]

        if isinstance(validation_lists, list):
            selected_classification = st.selectbox(
                "Select the classification to see:", validation_lists
            )
        else:
            selected_classification = validation_lists

        if len(validation_lists) >=2 and all([lst[-1].endswith("stage2") for lst in validation_lists]):
            st.warning("This is potentially a chimeric spectrum since it was classified in more than one Stage2 query", icon='❗️')


        if selected_classification:
            ba_tree_fig = create_custom_tree(selected_classification, selected_feature)
            st.plotly_chart(ba_tree_fig)

    with class_tab:
        add_df_and_filtering(filtered_classifications, key_prefix="class_table", default_cols=default_cols)
        with st.expander("How to interpret this table"):
            st.markdown("""
            The **"classification"** column displays the queries that support the compound's annotation as the most likely isomer.  
            The **"query_validation"** column lists all queries that matched a given spectra (scan).
            """)

    with lib_tab:
        only_library_matches = only_library_matches.merge(filtered_classifications[['#Scan#', 'classification']], on='#Scan#', how='left')
        only_library_matches = only_library_matches[default_cols + [col for col in only_library_matches.columns if col not in default_cols]]
        add_df_and_filtering(only_library_matches, key_prefix="lib_matches")

    with full_tab:
        full_table = full_table.merge(filtered_classifications[['#Scan#', 'classification']], on='#Scan#', how='left')
        full_table = full_table[default_cols + [col for col in full_table.columns if col not in default_cols]]
        add_df_and_filtering(full_table, key_prefix="full")
