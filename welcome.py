import streamlit as st


def welcome_page():
    st.title("Multi-step MassQL Bile Acid Isomer Annotation App")
    st.markdown("""
    ### 📖 Citation
    1. Ipsita Mohanty, Shipei Xing, Vanessa Castillo, Julius Agongo, Abubaker Patan, Yasin El Abiead, et al. 
    **MS/MS Mass Spectrometry Filtering Tree for Bile Acid Isomer Annotation**
    bioRxiv 2025.03.04.641505; doi: https://doi.org/10.1101/2025.03.04.641505
    2. Damiani, T., Jarmusch, A.K., Aron, A.T. et al. 
    **A universal language for finding mass spectrometry data patterns**. 
    Nat Methods 22, 1247–1254 (2025). https://doi.org/10.1038/s41592-025-02660-z
    
    ### 🧭 Purpose
    This application is designed to support the annotation of bile acid isomers in mass spectrometry data using **sequential MassQL queries**. It extends the functionality of standard MassQL searches by enabling **multi-step spectral pattern queries**, which are particularly useful for distinguishing structurally similar compounds such as bile acids.
    
    ### 📘 How It Works
    1. **Input a GNPS FBMN Task ID**: Enter your Feature-Based Molecular Networking (FBMN) task ID from GNPS. The MGF file will be downloaded and used to run the queries.
    2. **MassQL Query Execution**:
       - The original MGF will be pre-processed (remove blank scans)
       - The app runs an initial filter on the MGF file using stage 1 queries to remove spectra that do not pass this stage.
       - It then performs a set of secondary, more specific MassQL queries only on the remaining spectra.
    3. **Library Match Integration**: The results from MassQL are merged with GNPS library matches from the FBMN job.
    4. **Classification**: A bile acid classification tree is used to assign query results to specific bile acid subclasses/isomers.
    
    ### 🧩 Features
    - 👓 Visual classification of features using an interactive bile acid ontology tree
    - 📚 Explore GNPS Library matches and full annotated results
    - 🔍 Rapid exploration of GNPS jobs using MassQL query patterns
    
    ### 🧪 Example Dataset
    To explore the functionality without inputting your own data, use the **"Load query example"** checkbox in the sidebar.
    
    ---    
    """)
    st.info("""
    - This application is part of the GNPS downstream analysis ecosystem known as **MetaboApps**.
    - If you encounter any issues or have suggestions, please reach out to the app maintainers.
    - [Checkout other tools](https://wang-bioinformatics-lab.github.io/GNPS2_Documentation/toolindex/#gnps2-web-tools)
    """)