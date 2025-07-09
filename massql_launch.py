import pandas as pd
import yaml
from massql import msql_engine
import logging



def run_massql(mgf_path: str, queries_dict: dict):
    logger = logging.getLogger(__name__)
    executed_queries = []
    all_query_results_list = []
    for query_name, query_string in queries_dict.items():
        logger.info(f"Running query: {query_name}")
        executed_queries.append(query_name)
        try:
            results_df = msql_engine.process_query(query_string, mgf_path, parallel=True)
        except KeyError:
            logger.error(f"KeyError encountered for query: {query_name}")
            results_df = pd.DataFrame()

        if len(results_df) == 0:
            all_query_results_list.append({"query": query_name, "scan_list": []})
        else:
            passed_scan_ls = results_df["scan"].values.tolist()
            passed_scan_ls = [int(x) for x in passed_scan_ls]
            all_query_results_list.append({"query": query_name, "scan_list": passed_scan_ls})

    return all_query_results_list


### In case we want to redirect the stdout to streamlit
# # Run the script file
# result = subprocess.Popen(['bash', 'run.sh'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
# stdout, stderr = result.communicate()
#
# # Display the terminal output
# st.text('\n'.join(stdout.decode().split('\n')[1:][:-1]))

##Maybe also this: https://discuss.streamlit.io/t/cannot-print-the-terminal-output-in-streamlit/6602/9
