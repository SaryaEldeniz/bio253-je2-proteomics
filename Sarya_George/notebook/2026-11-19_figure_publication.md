# Publish the updated network and overlap figures

The student requested updating GitHub with the new figure. Prepared the complete
75-protein poster network in PNG, PDF and SVG, the earlier 49/67 comparison
figures, and the three-protein phosphosite/parent sample-value figure. Included
the corresponding R scripts, selected derived evidence tables, identifier and
annotation provenance, input lists and dated notebook records.

The 75-protein figure retains all eight restored candidates, Stk1/Stp1 at the
centre, fresh STRING KEGG outlines and unchanged source interaction halos.
The network contains 75 nodes and 412 functional associations at confidence
>=0.4. The source/output reconciliation script passed again before publication.
The overlap figure concerns IsaA, Pyk and QoxB; it preserves missing values and
the distinction between total protein and phosphosite signal.

Publication scope excludes raw experimental data, raw API downloads/caches,
full sequence files, machine configuration, the abstract draft and older
untracked design/slide drafts. The published README and R index are scoped to
the files included in this update; local abstract changes are retained.
An independent read-only agent checked the selected report/figure folders for
raw inputs, private paths, credentials and broken local Markdown links.

GitHub main was checked against the local base commit before preparing the
update. No source DEA effect, FDR or historical enrichment result was changed.
AI assistance: Codex prepared the publication scope and ran source/output
checks; a second agent independently audited publication contents.

Automated publication was blocked in this environment: Git's Windows shell
credential helper could not start, the connected GitHub app rejected the blob
upload with HTTP403 (Resource not accessible by integration), and the native
Windows credential-manager route was unavailable. The complete scoped update
was committed locally and retained for a normal authenticated Git push. No
partial figure upload or remote branch change was reported as successful.
