# Bio253 – *S. aureus* JE2 proteome & phosphoproteome: Sarya Repo

Starter code for the Bio253 block course (UZH/ETHZ, HS2026). You work in small groups
(2–3 people) on a **fresh, unpublished** proteomics data set of *Staphylococcus aureus*
JE2. You'll formulate hypotheses, test them in R, and use an AI coding agent
along the way. Use it carefully (see below).

Lecturer: Jonas Grossmann and José Luis Ruiz Rodriguez (Functional Genomics Center Zurich, FGCZ).

> **The raw data are not in this repository and are never committed** (see
> [data/README.md](data/README.md)). Your code and rendered reports are — you hand in
> your work by pushing your clone of this repo to GitHub (see *Handing in*).

## Biology in 10 lines or even more

1. *S. aureus* and *Pseudomonas aeruginosa* often co-infect, e.g. in cystic fibrosis lungs. *P. aeruginosa* usually suppresses *S. aureus* with secreted molecules (HQNO, PQS, pyocyanin, siderophores, proteases …).
2. In this experiment, *S. aureus* **never meets living *P. aeruginosa***. It is grown in medium containing **30 % cell-free, sterile-filtered spent medium of *P. aeruginosa* PAO1** ("PA supernatant", **PASN**). It sees only the secreted molecules, in a medium that is partly depleted of nutrients.
3. Niggli *et al.* 2023 (*mBio*): *S. aureus* quickly evolves resistance to PASN. Schwyter *et al.* 2025 (*mSystems*): the resistance comes mostly from a large **plastic** proteome response. One evolved line also had a ***tcyA*** mutation, which destroyed the cystine transporter TcyABC.
4. Bacteria sense their environment through phosphorylation. *S. aureus* has one main Ser/Thr kinase, **Stk1 (PknB)**, and its partner phosphatase **Stp1**. The genes sit next to each other: *stp1* first, then *stk1*.
5. Huemer *et al.* 2021 (*Sci Signal*): PknB/Stp control quiescence and antibiotic tolerance. Prust *et al.* 2021 (*MCP*): a deep *S. aureus* phosphoproteome, which also reports **phospho-histidine**.
6. Here: the JE2 wild type (**anc**) and two transposon mutants, **stk** (*stk1::Tn*) and **stp** (*stp1::Tn*), each grown in three media.
7. **TSB** = control. **SN** = ~15 h in 70 % TSB + 30 % PASN (long-term response). **SNpulse** = grown in TSB, then 15 min in PASN (immediate response).
8. Every sample was measured twice: **total proteome** (protein amounts) and **phospho-enriched** (Fe-IMAC; phospho-sites).
9. The phospho data were searched twice, with phospho allowed on **S/T/Y** and on **S/T/Y/H**.
10. Open questions: what changes in PASN, and how fast? What depends on Stk1/Stp1? Is there pHis? Do plastic responses mirror evolved ones?

## Design

| | TSB | SNpulse (15 min PASN) | SN (~15 h PASN) |
|---|---|---|---|
| **anc** (JE2 wild type) | 3 | 3 | 3 |
| **stk** (*stk1*::Tn, kinase) | 3 | 3 | 3 |
| **stp** (*stp1*::Tn, phosphatase) | 3 | 3 | 3 |

27 samples × 2 measurements (total proteome, phospho-enriched).
File names write `pulseSN`; the tables and contrasts write `SNpulse`.

### The one rule you must not get wrong

A contrast named **`A_vs_B`** has **`diff = log2(A) − log2(B)`**.
So `TSB_vs_SN_at_anc > 0` means *higher in TSB*, which is *lower in PASN*.
If in doubt, run `explain_contrast("TSB_vs_SN_at_anc")`.

## Getting the data

The data are **not** in this repository. You get them from the lecturer as **one zip file**.

1. Extract it into the `data/` folder of your clone (you then have `data/raw/`, `data/derived/`, `data/reference/`).
2. Verify all files against the MD5 checksums in `data/checksums.md5`:
   `bash data/check_data.sh` (macOS/Linux) or `Rscript data/check_data.R` (Windows). Every file must say `OK`.

Details and troubleshooting: [data/README.md](data/README.md). Everything in `data/` except
the README, the checksum file and the check scripts is git-ignored and never pushed.

## Quick start

```bash
git clone <this repo>   # or unzip the folder you got
cd bio253-je2-proteomics
Rscript install_packages.R              # once
unzip /path/to/bio253_je2_data.zip -d data   # extract the course data zip into data/
bash data/check_data.sh                 # all files must be OK
Rscript tests/run_tests.R               # should end with "DONE" and no failures
cd analyses && quarto render            # renders all notebooks to analyses/_output/
```

In R / RStudio (open `bio253-je2-proteomics.Rproj`):

```r
source("R/setup.R")
find_protein("cystine")                     # search by locus tag, product, protein id
plot_protein("B7H15_13375")                 # per-sample abundances in all 9 groups
dea <- read_dea("total")                    # all contrasts, long format
volcano(dea, "TSB_vs_SN_at_anc")
sites <- read_sites("STYH")                 # FragPipe phospho-sites, long format
```

## Repository layout

| path | what |
|---|---|
| `R/` | helper functions, loaded with `source("R/setup.R")` |
| `scripts/00_fetch_data.R` | lecturer only: rebuilds `data/` from the original course folder |
| `data/checksums.md5`, `data/check_data.sh/.R` | MD5 check of everything extracted from the course data zip |
| `analyses/` | exercise notebooks (Quarto, R) |
| `notebook/` | **your** lab notebook: one dated entry per session |
| `tests/` | `testthat` checks of the helpers and of the sign convention |
| `data/` | local data only (git-ignored), see `data/README.md` |
| `CLAUDE.md` / `AGENTS.md` | instructions for your AI coding agent |

## Exercises

| notebook | question |
|---|---|
| `01_explore_data.qmd` | What is in the data? PCA, correlation, missing values, hits per contrast |
| `02_hypothesis_TcyABC.qmd` | **Flagship.** Is the cystine transporter TcyABC down-regulated in PASN? |
| `03_mutant_validation.qmd` | Are the *stk* and *stp* mutants what we think they are? |
| `04_phospho_two_roadmaps.qmd` | Sparse phospho data: reproducible vs exploratory analysis |
| `05_histidine_phospho.qmd` | Is there phospho-histidine, and how sure can we be? |
| `06_template_your_question.qmd` | Template for your own question |

Each notebook renders on its own in under a minute. Start with 01 and 02. Then pick
a topic as a group and copy 06.

### Ideas for your own question

* Which pathways change in SN? Use the STRING enrichment terms in `data/raw/resources/`.
* Do virulence factors go down and ROS scavengers (AhpC, KatA) go up, as in *mSystems*?
* Which phospho-sites depend on Stk1? Compare `stk` with `anc` and use the
  phospho-sites from Prust *et al.* as a reference.
* Is there a fast phospho-response in SNpulse when the proteome hardly changes?
* Correct phospho-site changes for changes in protein amount.

## Handing in: your repo on GitHub

Each group works in its own clone of this starter repo and hands in by pushing it to GitHub.

```bash
git clone <starter-repo-url> bio253-<group>      # or use "Use this template" / fork on GitHub
cd bio253-<group>
# ... work, commit small steps ...
git remote set-url origin https://github.com/<you>/bio253-<group>.git   # your own repo
git push -u origin main
```

What goes in: all code (`R/`, `scripts/`, `analyses/*.qmd`, your own `reports/`), the rendered
HTML reports (`analyses/_output/`, `reports/*.html`), and `notebook/`.
What stays out: everything in `data/` (raw Excel/TSV/parquet/FASTA — `.gitignore` takes care of it).
Before the final push run `git status` and check that no data file is staged.

## Working with an AI coding agent – responsibly

You're encouraged to use an agent (Claude Code, Copilot, Cursor …). Some rules:

1. **Use whatever AI tools help you** (Claude Code, Codex, Cursor, Copilot, chat
   assistants …), directly on the course data. Just never commit the raw files in
   `data/`. `CLAUDE.md` / `AGENTS.md` tell the agents the same.
2. **You are responsible for every line.** Ask the agent to explain code you don't
   understand. Run it, look at the output, and check one number by hand.
3. **Plot the raw values first.** Before you believe a fold change, look at
   `plot_protein()` for that protein.
4. **Ask for the sign.** Agents (and people) often get `A_vs_B` backwards.
5. **Verify claims about biology.** Agents make up references and protein
   functions. Check them in UniProt, PubMed or the papers, and cite what you used.
6. **Keep a lab notebook** in `notebook/`: what you asked, what you ran, what you
   found, and what was wrong.
7. **Small steps.** Ask for one plot or one table at a time, not "analyse everything".

## Literature

* Schwyter L. *et al.* (2025) *mSystems*, doi:10.1128/msystems.00262-25 – TcyABC, plastic vs genetic response to PASN.
* Niggli S. *et al.* (2023) *mBio* – rapid, strain-specific resistance evolution of *S. aureus* against PASN.
* Prust N. *et al.* (2021) *Mol Cell Proteomics* 20:100034 – *S. aureus* phosphoproteome, Stk1 targets, pHis.
* Huemer M. *et al.* (2021) *Science Signaling*, doi:10.1126/scisignal.abj8194 – PknB/Stp, quiescence and antibiotic tolerance.

## License

The code is MIT-licensed (see `LICENSE`). The **data are not part of this repository
and are not covered by the license**: they are unpublished and remain property of
the authors/FGCZ.
