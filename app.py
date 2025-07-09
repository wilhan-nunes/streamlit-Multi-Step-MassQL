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
only_stage1 = massql_queries.only_stage1
stage2 = massql_queries.stage2
mono_queries = massql_queries.mono_queries
di_queries = massql_queries.di_queries
tri_queries = massql_queries.tri_queries

# Set page configuration
page_title = "Multi-Step MassQL Bile Acid Isomer Annotation"

# TODO: Bump version
app_version = "2025-07-07"

st.set_page_config(
    page_title=page_title,
    layout="wide",
    # page_icon=":pill:",
    menu_items={"About": ("**App version: %s**" % app_version)},
)


def process_results(massql_results_df: List):
    all_query_results_df = pd.DataFrame(massql_results_df)
    all_query_results_df["scan_list"] = all_query_results_df["scan_list"].apply(
        lambda x: ast.literal_eval(x) if isinstance(x, str) else x
    )
    all_query_results_df = all_query_results_df.explode("scan_list")
    all_query_results_df = all_query_results_df.rename(
        columns={"scan_list": "#Scan#", "query": "query_validation"}
    )

    with st.spinner("Merging results..."):
        all_query_results_df["#Scan#"] = all_query_results_df["#Scan#"].astype(str)
        library_matches["#Scan#"] = library_matches["#Scan#"].astype(str)

        library_and_query_results = pd.merge(
            library_matches, all_query_results_df, on="#Scan#", how="outer"
        )

        library_matches_only = pd.merge(
            library_matches, all_query_results_df, on="#Scan#", how="left"
        )
        fallback_label = "Did not pass any selected query"
        library_and_query_results["query_validation"] = library_and_query_results[
            "query_validation"
        ].fillna(fallback_label)

        library_and_query_results = library_and_query_results[
            ["query_validation", "Compound_Name"]
            + [
                col
                for col in library_and_query_results.columns
                if col not in ["query_validation", "Compound_Name"]
            ]
        ]

        library_and_query_results = library_and_query_results.groupby(
            "#Scan#", as_index=False
        ).agg(
            {
                "query_validation": lambda x: ";".join(set(x)),
                **{
                    col: "first"
                    for col in library_and_query_results.columns
                    if col not in ["#Scan#", "query_validation"]
                },
            }
        )

    return library_matches_only, library_and_query_results, all_query_results_df


def cleanup_massql_files():
    feather_files = glob.glob("temp_mgf/*.feather")
    for file in feather_files:
        try:
            os.remove(file)
        except Exception as e:
            st.warning(f"Could not delete {file}: {e}")


def load_example_data():
    if "only_library_matches" not in st.session_state:
        st.session_state.only_library_matches = pd.read_csv(
            "examples/example_library_matches.csv", dtype=str
        )
    if "all_query_results_df" not in st.session_state:
        st.session_state.all_query_results_df = pd.read_csv("examples/example_all_results.csv", dtype=str)
    if "library_and_query_results" not in st.session_state:
        st.session_state.library_and_query_results = pd.read_csv(
            "examples/example_lib_and_query_results.csv", dtype=str
        ).fillna("No match")


def get_bile_acids_classifications(results_df):
    passed_queries = results_df[
        results_df["query_validation"] != "Did not pass any selected query"
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
    load_example = st.checkbox("Load query example", value=False)
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
            st.rerun()

    st.subheader("Contributors")
    st.markdown(
        """
    - [Helena Russo PhD](https://sites.google.com/view/helenamrusso/home) - UC San Diego
    - [Wilhan Nunes PhD](https://scholar.google.com/citations?user=4cPVoeIAAAAJ) - UC San Diego
    - [Ipsita Mohanty PhD](https://scholar.google.com/citations?user=iHJ3vgsAAAAJ) - UC San Diego
    - [Mingxun Wang PhD](https://www.cs.ucr.edu/~mingxunw/) - UC Riverside
    """
    )

    st.subheader("Documentations and Resources")
    st.write(
        """<a href="https://cmmc.gnps2.org/network_enrichment/">CMMC Enrichment Workflow</a><br>
                        <a href="https://wang-bioinformatics-lab.github.io/GNPS2_Documentation/fbmn/">Feature Based Molecular Networking</a><br>
                        <a href="https://cmmc-kb.gnps2.org" target="_blank">CMMC Knowledge Base</a>""",
        unsafe_allow_html=True,
    )

if not run_query and "run_query_done" not in st.session_state:
    from welcome import welcome_page

    welcome_page()

if run_query:
    st.session_state["run_query_done"] = True
    st.title("Multi-step MassQL Results")
    if not load_example:
        with st.spinner("Downloading files and running queries..."):
            library_matches = workflow_fbmn.get_library_match_dataframe(task_id)
            cleaned_mgf_path, all_scans = download_and_filter_mgf(task_id)
            mgf_path = cleaned_mgf_path

        with st.spinner("Running Stage 1 queries..."):
            stage1_all_results = massql_launch.run_massql(
                mgf_path, queries_dict=only_stage1
            )
            stage1_results_df = pd.DataFrame(stage1_all_results)
            # create a new mgf filtering to just maintain the scans that passed stage1
            scans_to_keep = set(sum(stage1_results_df["scan_list"], []))
            stage1_passed_mgf = filter_mgf_by_scans(
                mgf_path,
                f"temp_mgf/{uuid.uuid4()}_stg1_passed.mgf",
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
                (
                    only_library_matches,
                    library_and_query_results,
                    all_query_results_df,
                ) = process_results(massql_results_df)

                # #TODO:remove later (saving example files)
                # pd.DataFrame(massql_results_df).to_csv('example_massql_results_after_stg1.csv', index=False)
                # library_and_query_results.to_csv('example_lib_and_query_resuts.csv', index=False)
                # only_library_matches.to_csv('example_library_matches.csv', index=False)
                # all_query_results_df.to_csv('example_all_results.csv', index=False)

            cleanup_massql_files()
        # Store results in session state
        st.session_state["only_library_matches"] = only_library_matches
        st.session_state["library_and_query_results"] = library_and_query_results
        st.session_state["all_query_results_df"] = all_query_results_df

    else:
        # this function stores the static result file dataframes in st.session_state, just as above.
        load_example_data()

if run_query or st.session_state.get("run_query_done"):
    only_library_matches = st.session_state["only_library_matches"]
    library_and_query_results = st.session_state["library_and_query_results"]
    library_and_query_results["Compound_Name"] = library_and_query_results[
        "Compound_Name"
    ].fillna("No match")

    filtered_classifications = get_bile_acids_classifications(library_and_query_results)
    feature_ids_dict = filtered_classifications[["#Scan#", "Compound_Name"]].astype(str)

    feature_ids_dict = feature_ids_dict.set_index("#Scan#")["Compound_Name"].to_dict()
    feature_ids_dict = dict(sorted(feature_ids_dict.items(), key=lambda item: item[1]))

    viz_tab, class_tab, lib_tab, full_tab = st.tabs(
        ["👓 Visualizations", "🗂️ Classified", "📚 Library Matches", "📋 Full Table",]
    )

    with viz_tab:
        st.subheader("Feature Classification")
        selected_feature = st.selectbox(
            "Select a feature:",
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

        if selected_classification:
            ba_tree_fig = create_custom_tree(selected_classification, selected_feature)
            st.plotly_chart(ba_tree_fig)

    with class_tab:
        default_cols = ["#Scan#", "Compound_Name", 'classification', 'query_validation']
        add_df_and_filtering(filtered_classifications, key_prefix="class_table", default_cols=default_cols)
        with st.expander("How to interpret this table"):
            st.markdown("""
            The **"classification"** column displays the queries that support the compound's annotation as the most likely isomer.  
            The **"query_validation"** column lists all queries that matched a given compound.
            """)

    with lib_tab:
        add_df_and_filtering(only_library_matches, key_prefix="lib_matches")

    with full_tab:
        add_df_and_filtering(library_and_query_results, key_prefix="full")
