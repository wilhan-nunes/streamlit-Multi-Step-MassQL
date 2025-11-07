import ast
import glob
import os
import uuid
from typing import List

import gnpsdata
import pandas as pd
from gnpsdata import workflow_fbmn, taskinfo
from streamlit.components.v1 import html

import massql_launch
from utils import (
    download_and_filter_mgf,
    filter_mgf_by_scans,
    highlight_hydroxy,
    MassQLQueries,
    bile_acid_tree,
    add_df_and_filtering,
    get_git_short_rev,
    gnps2_get_library_match_dataframe,
    insert_mgf_info,
    insert_plot_download_button,
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

EXAMPLE_CONFIG = {
    "task_id": "a322acf7936c4f91a41fd2f267d9b613",
    "description": "HNRC cohort samples of 10 cognitively impaired, 10 non impaired patients, all from the HIV+ group.",
}

# Set page configuration
page_title = "Multi-Step MassQL Bile Acid Isomer Annotation"

# TODO: Bump version
app_version = "2025-11-07"
git_hash = get_git_short_rev()
repo_link = "https://github.com/wilhan-nunes/streamlit-Multi-Step-MassQL"

st.set_page_config(
    page_title=page_title,
    layout="wide",
    page_icon="🪜",
    initial_sidebar_state="expanded",
    menu_items={
        "About": (
            f"**App version**: {app_version} | "
            f"[**Git Hash**: {git_hash}]({repo_link}/commit/{git_hash})"
        )
    },
)

# Add a tracking token
html(
    '<script async defer data-website-id="<your_website_id>" src="https://analytics.gnps2.org/umami.js"></script>',
    width=0,
    height=0,
)


def process_results(
    massql_results_df: List, library_matches: pd.DataFrame, all_scans: List[str]
):
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
        full_table = all_scans_df.merge(library_matches, on="#Scan#", how="left").merge(
            all_query_results_df, on="#Scan#", how="left"
        )

        # Fill missing values
        full_table["query_validation"] = full_table["query_validation"].fillna(
            "Did not pass stage1 filtering"
        )

        # Library matches only (existing functionality)
        library_matches_only = library_matches.merge(
            all_query_results_df, on="#Scan#", how="left"
        )

        # Reorder columns and aggregate
        cols = ["query_validation", "Compound_Name"] + [
            col
            for col in full_table.columns
            if col not in ["query_validation", "Compound_Name"]
        ]
        full_table = full_table[cols]

        # Group by scan and aggregate
        full_table = full_table.groupby("#Scan#", as_index=False).agg(
            {
                "query_validation": lambda x: ";".join(set(x.dropna())),
                **{
                    col: "first"
                    for col in full_table.columns
                    if col not in ["#Scan#", "query_validation"]
                },
            }
        )

    return library_matches_only, full_table, all_query_results_df


def cleanup_massql_files():
    feather_files = glob.glob("temp_mgf/*.feather")
    for file in feather_files:
        try:
            os.remove(file)
        except Exception as e:
            st.warning(f"Could not delete {file}: {e}")


# ============================================================================
# Main Processing Pipeline
# This section contains the core processing logic that can be cached or not
# based on whether the task_id matches the example configuration.
# ============================================================================


def process_task_pipeline(task_id: str):
    """
    Main processing pipeline that executes all steps for a given task_id.
    Returns all necessary results for downstream processing.

    Returns:
        tuple: (only_library_matches, full_table, all_query_results_df,
                filtered_classifications, feature_ids_dict)
    """
    # Step 1: Get library matches
    task_info = taskinfo.get_task_information(task_id)
    workflowname = task_info.get("workflowname")
    if workflowname == "feature_based_molecular_networking_workflow":
        library_matches = workflow_fbmn.get_library_match_dataframe(task_id)
    elif workflowname == "classical_networking_workflow":
        library_matches = gnps2_get_library_match_dataframe(task_id)

    # Step 2: Download and filter MGF
    cleaned_mgf_path, all_mgf_scans = download_and_filter_mgf(task_id)
    mgf_path = cleaned_mgf_path

    # Step 3: Run Stage 1 MassQL queries
    stage1_all_results = massql_launch.run_massql(mgf_path, queries_dict=stage1)
    stage1_results_df = pd.DataFrame(stage1_all_results)
    scans_to_keep = set(sum(stage1_results_df["scan_list"], []))
    stage1_passed_mgf = filter_mgf_by_scans(
        mgf_path,
        f"temp_mgf/{task_id}_stg1_passed.mgf",
        scans_to_keep,
    )

    # Step 4: Run all MassQL queries on filtered data
    massql_results_df = massql_launch.run_massql(stage1_passed_mgf, ALL_MASSQL_QUERIES)

    # Step 5: Process results into tables
    (
        only_library_matches,
        full_table,
        all_query_results_df,
    ) = process_results(massql_results_df, library_matches, all_mgf_scans)

    # Step 6: Fill missing compound names
    full_table["Compound_Name"] = full_table["Compound_Name"].fillna("No match")

    # Step 7: Extract bile acid classifications
    filtered_classifications = get_bile_acids_classifications(
        full_table, exclude_string="did not pass"
    )

    # Step 8: Create feature IDs dictionary for visualization
    feature_ids_dict = {}
    if len(filtered_classifications) > 0:
        feature_ids_dict_df = filtered_classifications[
            ["#Scan#", "Compound_Name"]
        ].astype(str)
        feature_ids_dict = feature_ids_dict_df.set_index("#Scan#")[
            "Compound_Name"
        ].to_dict()
        feature_ids_dict = dict(
            sorted(feature_ids_dict.items(), key=lambda item: item[1])
        )

    return (
        only_library_matches,
        full_table,
        all_query_results_df,
        filtered_classifications,
        feature_ids_dict,
    )


@st.cache_data
def _cached_process_task_pipeline(task_id: str):
    """Cached version of the main processing pipeline."""
    return process_task_pipeline(task_id)


def run_task_analysis(task_id: str):
    """
    Wrapper function that decides whether to use cached processing or not
    based on whether task_id matches the example configuration.
    """
    if task_id == EXAMPLE_CONFIG["task_id"]:
        return _cached_process_task_pipeline(task_id)
    else:
        return process_task_pipeline(task_id)


def get_bile_acids_classifications(results_df, exclude_string: str):
    """Extract bile acid classifications from results dataframe."""
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


task_id_value = st.query_params.get("task_id", "")

with st.sidebar:
    st.subheader("Analysis configuration")
    load_example = st.checkbox(
        "Load query example", value=False, key="load_example_checkbox"
    )
    if load_example:
        task_id_value = EXAMPLE_CONFIG["task_id"]
        st.info(
            f"[FBMN Job](https://gnps2.org/status?task={task_id_value}): {EXAMPLE_CONFIG['description']}",
            icon="ℹ️",
        )

    task_id = st.text_input(
        "Task ID:",
        value=task_id_value,
        disabled=load_example,
        help="Enter a valid GNPS task ID (Feature Based or Classical MN) to run the analysis.",
    )

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
    st.markdown(
        """
    [Feature Based Molecular Networking](https://wang-bioinformatics-lab.github.io/GNPS2_Documentation/fbmn/)<br>
    [MassQL documentation](https://mwang87.github.io/MassQueryLanguage_Documentation/)
    """,
        unsafe_allow_html=True,
    )

if not run_query and "run_query_done" not in st.session_state:
    from welcome import welcome_page

    welcome_page()

if run_query:
    st.session_state["run_query_done"] = True

    # Execute the entire processing pipeline (includes all table processing)
    with st.spinner("Processing task... This may take a while for new tasks."):
        (
            only_library_matches,
            full_table,
            all_query_results_df,
            filtered_classifications,
            feature_ids_dict,
        ) = run_task_analysis(task_id)

    cleanup_massql_files()

    # Store results in session state
    st.session_state["only_library_matches"] = only_library_matches
    st.session_state["full_table"] = full_table
    st.session_state["all_query_results_df"] = all_query_results_df
    st.session_state["filtered_classifications"] = filtered_classifications
    st.session_state["feature_ids_dict"] = feature_ids_dict


if run_query or st.session_state.get("run_query_done"):
    st.title("🔢 Multi-step MassQL Results")
    only_library_matches = st.session_state["only_library_matches"]
    full_table = st.session_state["full_table"]
    filtered_classifications = st.session_state["filtered_classifications"]
    feature_ids_dict = st.session_state["feature_ids_dict"]

    # Check if we have classifications to display
    if len(filtered_classifications) == 0:
        st.warning(
            "No classifications retrieved for this task ID. Inspect the full table below for details"
        )
        st.write(full_table)
        st.stop()

    viz_tab, class_tab, lib_tab, full_tab = st.tabs(
        [
            "👓 Visualizations",
            "🗂️ Classified",
            "📚 Library Matches",
            "📋 Full Table",
        ]
    )

    default_cols = ["#Scan#", "Compound_Name", "classification"]

    with viz_tab:
        st.subheader("Feature Classification")

        @st.fragment
        def render_tree():
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
                if len(validation_lists) >= 2:
                    st.warning(
                        "This is potentially a chimeric spectrum since it was classified in more than one Stage2 query",
                        icon="❗️",
                    )
                    selected_classification = st.selectbox(
                        "Select the classification to see:", validation_lists
                    )

                else:
                    selected_classification = validation_lists[0]

            if selected_classification:

                ba_tree_fig = create_custom_tree(
                    selected_classification, selected_feature[:20]
                )
                st.plotly_chart(
                    ba_tree_fig,
                    config={
                        "toImageButtonOptions": {
                            "height": None,
                            "width": None,
                            "format": "svg",
                        }
                    },
                )

                if st.button(
                    "Generate SVG",
                    type="secondary",
                    icon=":material/manufacturing:",
                    key="download_tree_svg",
                ):

                    with st.spinner("Generating SVG..."):
                        svg_bytes = ba_tree_fig.to_image(
                            format="svg", width=1300, height=600, scale=1
                        )
                        insert_plot_download_button("bile_acid_tree", svg_bytes, "tree")

        render_tree()

    with class_tab:

        @st.fragment
        def render_classified_table():
            class_df = add_df_and_filtering(
                filtered_classifications,
                key_prefix="class_table",
                default_cols=default_cols,
            )
            st.dataframe(
                class_df.style.apply(highlight_hydroxy, subset=["classification"])
            )

            # Download button for classified table
            csv_classified = class_df.to_csv(index=False).encode("utf-8")
            st.download_button(
                label="Download Classified Table as CSV",
                data=csv_classified,
                file_name=f"{task_id}_classified_results.csv",
                mime="text/csv",
                key="download_classified_csv",
                icon=":material/download:",
            )

            with st.expander("How to interpret this table"):
                st.markdown(
                    """
                The **"classification"** column displays the queries that support the compound's annotation as the most likely isomer.  
                The **"query_validation"** column lists all queries that matched a given spectra (scan).
                """
                )

        render_classified_table()

    with lib_tab:
        only_library_matches = only_library_matches.merge(
            filtered_classifications[["#Scan#", "classification"]],
            on="#Scan#",
            how="left",
        )
        only_library_matches = only_library_matches[
            default_cols
            + [col for col in only_library_matches.columns if col not in default_cols]
        ]

        @st.fragment
        def render_library_matches_table():
            library_df = add_df_and_filtering(
                only_library_matches, key_prefix="lib_matches"
            )
            st.dataframe(library_df)

            # Download button for library matches table
            csv_library = library_df.to_csv(index=False).encode("utf-8")
            st.download_button(
                label="Download Library Matches as CSV",
                data=csv_library,
                file_name=f"{task_id}_library_matches.csv",
                mime="text/csv",
                key="download_library_csv",
                icon=":material/download:",
            )

        render_library_matches_table()

        with full_tab:
            full_table = full_table.merge(
                filtered_classifications[["#Scan#", "classification"]],
                on="#Scan#",
                how="left",
            )
            full_table = full_table[
                default_cols
                + [col for col in full_table.columns if col not in default_cols]
            ]

            @st.fragment
            def render_full_table():
                full_df = add_df_and_filtering(full_table, key_prefix="full")
                st.dataframe(full_df)

                # Download button for full table
                csv_full = full_df.to_csv(index=False).encode("utf-8")
                st.download_button(
                    label="Download Full Table as CSV",
                    data=csv_full,
                    file_name=f"{task_id}_full_results.csv",
                    mime="text/csv",
                    key="download_full_csv",
                    icon=":material/download:",
                )

            render_full_table()

    st.markdown("---")
    st.subheader("Download MGF with validated scans")

    @st.fragment
    def render_mgf_download():
        if st.button(
            "Generate MGF with validated scans",
            type="secondary",
            icon=":material/manufacturing:",
        ):
            input_mgf = f"./temp_mgf/{task_id}_stg1_passed.mgf"
            buf = insert_mgf_info(
                task_id,
                input_mgf,
                full_table[["#Scan#", "query_validation", "classification"]].astype(
                    str
                ),
            )
            st.download_button(
                label="Download validated MGF",
                data=buf.getvalue(),
                file_name=f"{task_id}_validated_scans.mgf",
                mime="txt/plain",
                icon=":material/download:",
                type="primary",
            )

    render_mgf_download()
