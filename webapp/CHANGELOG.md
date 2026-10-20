# Changelog — Multi-step MassQL web app

All notable changes to `multistep_massql.html` are listed here. This is the
client-side version of the Multi-step MassQL bile-acid isomer app: one
self-contained HTML file that runs the multi-step MassQL pipeline in the
browser, without Streamlit or a server.

- **Live:** <https://apps.gnps2.org/multistep_massql>
- **Releases** are named by the date shown as **App version** in the page footer.
- **Source:** the file is generated from `build/metaboapps/` in the
  [NextGenWebApps](https://github.com/Wang-Bioinformatics-Lab/NextGenWebApps)
  repository (`multistep.head.html` and `multistep.app.js`, assembled by
  `build.py`). It embeds this repository's `massql_queries.yaml`,
  `bile_acid_tree.yaml`, `ba_tree_mono_tri.png` and `ba_tree_di.png`, copied
  into `build/metaboapps/data/`. When any of those change here, update the
  copies there, rebuild, and copy the new file back. If you edit this copy by
  hand, the next build overwrites your change.
- **Running it locally:** you can open the file directly or serve it from
  anywhere, and dropping an `.mgf` works offline. Fetching a GNPS2 task by ID
  only works from `apps.gnps2.org`, because GNPS2 allows only that origin to
  request its result files.

## [2026-10-20]

### Added
- A **Load example** button next to **Fetch task**. It loads the public GNPS2
  task `a322acf7936c4f91a41fd2f267d9b613` from `gnps2.org`.
- The two bile-acid ontology figures from the Streamlit app's welcome page
  (mono- and trihydroxy, and dihydroxy), shown in a collapsible card under the
  introduction. Click a figure to see it large. In the enlarged view, the arrow
  buttons or the ← / → keys switch figures, and Esc or a click outside closes
  it.
- **Classified** tab: a `most_specific_class` column, the deepest tree node the
  scan reached (for example `Di-3,7-OH|Di-3,12b-OH|Di-7,12b-OH`). It comes first,
  and the full path stays next to it. A chimeric scan shows one class per path.
  The classified TSV export gains the same column, and its `classification`
  column (the full path) is unchanged.
- **Classified** tab: a filter panel, which replaces the single text filter. It
  filters by:
  - most specific class;
  - tree node, which keeps scans whose path passes through that node;
  - text in the full path or the query validation;
  - scan and precursor m/z, which accept comparisons such as `>=400` or `<300`;
  - chimeric or not.

  Active fields are highlighted, and a **Clear filters** button and a "Showing X
  of Y scans" line appear as in the GNPS2-AnnEx app.
- **Visualizations** tab: the feature dropdown is now a search box, like the
  Streamlit app's select box. Type a scan number (scans that start with it are
  listed first), a compound name or a class. Pick a result with ↑ / ↓ and
  Enter, or by clicking. ← / → buttons step to the previous or next feature,
  and a counter shows the position (for example "33 of 143 classified").

### Changed
- Clicking a node in the classification tree now sets the **Tree node** filter
  on the **Classified** tab. Before, it typed the node's name into the text
  filter, which could also match other nodes whose names contain it.
- The selected feature's path is now drawn in its level colors, with a dark
  outline and dark connecting lines, and the rest of the tree is faded. Before,
  the path was filled light blue.

### Fixed
- The tree's level colors never appeared. Every node was drawn near-white
  because a stylesheet rule overrode the color set for its level, so the legend
  matched nothing in the tree. The nodes now use the palette from
  `tree_plotter.py`, with the light blue darkened a little so white text stays
  readable. The legend now lists all six levels (category, stage 1, stage 2 and
  isomer levels 1–3). Before, it listed four.

## [2026-10-06]

### Added
- The **App version** date in the page footer.

### Fixed
- Results were shifted onto the wrong scans when the MGF contained empty
  entries. FBMN writes an entry without peaks for every feature that has no
  MS2, and the MassQL engine carried that entry's `SCANS=` and `PEPMASS=` into
  the next spectrum. So every real spectrum after an empty entry was reported
  under the wrong scan number, with a precursor m/z of 0. Empty entries are now
  skipped before the pipeline runs, and the file summary says how many were
  skipped.

## [2026-08-17]

The first release. The page did not show a version yet.

### Added
- A client-side port of this Streamlit app. Give it a GNPS2 task ID and a server
  (`gnps2.org` or `beta.gnps2.org`), and it fetches the task's consensus MGF and
  library matches. Or drop an `.mgf` (and optionally a library-matches `.tsv`)
  from your computer.
- The same fixed pipeline as the Streamlit app, run in the browser:
  1. the stage-1 hydroxylation filter;
  2. an MGF with only the scans that passed stage 1;
  3. the stage-2 confirmation and isomer queries;
  4. each scan's classification down the bile-acid tree.
- MassQL runs with the same Rust engine, compiled to WebAssembly, as the GNPS2
  MassQL Playground. The engine is embedded in the file, so nothing is uploaded.
- The query set and the classification tree are copied unchanged from
  `massql_queries.yaml` and `bile_acid_tree.yaml`. The classifier is a
  JavaScript port of `tree_classifier.py`. It gave the same result as the
  Python code for all 196 classified scans of GNPS2 task
  `58e0e2959ec748049cb2c5f8bb8b87dc`.
- **Visualizations** (the classification tree, with the selected feature's
  path highlighted, and an SVG download), **Classified**, **Library Matches**,
  **Full Table** and **Executed Queries** tabs.
- TSV exports, and a **validated MGF** export.
- Scans classified on more than one stage-2 path are marked as potentially
  chimeric spectra.
