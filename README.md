# Multi-step MassQL

Annotate bile acid isomers in a GNPS2 molecular networking job with a fixed
sequence of MassQL queries. Each spectrum goes through three steps, and the
deepest step it passes is its annotation:

1. A stage 1 filter: is it a mono-, di- or trihydroxy bile acid?
2. Stage 2 confirmation queries, run only on the spectra that passed stage 1.
3. Isomer-specific queries.

The queries come from the bile acid MS/MS filtering tree in Mohanty et al.
(see Citation). The figures `ba_tree_mono_tri.png` and `ba_tree_di.png` show
the tree.

## Use it online

- **Streamlit app:** <https://multistep-massql.gnps2.org/>
- **Browser version:** <https://apps.gnps2.org/multistep_massql>. It's the same
  pipeline in a single HTML file, with the queries running in your browser, and
  it adds filters on the results table and a search box for scans. See
  [`webapp/`](webapp/) and its [changelog](webapp/CHANGELOG.md).

Both take the task ID of a Feature-Based or Classical Molecular Networking job,
and both include an example task.

## Output

- The classification of each spectrum, from bile acid class down to isomer.
  Spectra that match more than one stage 2 class are flagged as possibly
  chimeric.
- The classification tree, with the selected spectrum's path highlighted.
- Tables that merge the classifications with the GNPS library matches,
  downloadable as TSV, plus an MGF with only the validated spectra.

## Files

| File | Contents |
|---|---|
| `massql_queries.yaml` | The MassQL queries, in the order they run |
| `bile_acid_tree.yaml` | The classification tree that the query names map onto |
| `tree_classifier.py` | Assigns each spectrum its deepest matching path in the tree |
| `app.py`, `welcome.py` | The Streamlit app |
| `webapp/` | The browser version |

The browser version contains copies of the two YAML files and the two tree
figures. If you change any of them, the browser version needs to be rebuilt
(see [webapp/CHANGELOG.md](webapp/CHANGELOG.md)).

## Run locally

```bash
git clone https://github.com/wilhan-nunes/streamlit-Multi-Step-MassQL.git
cd streamlit-Multi-Step-MassQL
pip install -r requirements.txt
streamlit run app.py
```

With Docker, run `make server-compose`.

The browser version needs no installation: download
[`webapp/multistep_massql.html`](webapp/multistep_massql.html) and open it in
your browser. It works offline with an `.mgf` file from your computer (and
optionally the library-matches `.tsv`). Loading a job by task ID, including
the example, only works on the online version, because GNPS2 only lets
`apps.gnps2.org` request its result files.

## Citation

- Mohanty I., Xing S., Castillo V. et al. *MS/MS Mass Spectrometry Filtering
  Tree for Bile Acid Isomer Annotation.* bioRxiv (2025).
  <https://doi.org/10.1101/2025.03.04.641505>
- Damiani T. et al. *A universal language for finding mass spectrometry data
  patterns.* Nature Methods 22, 1247–1254 (2025).
  <https://doi.org/10.1038/s41592-025-02660-z>

Part of the GNPS2 MetaboApps. See also
[PostMN MassQL](https://massqlpostmn.gnps2.org/) for running any MassQL query on
a molecular networking job.
