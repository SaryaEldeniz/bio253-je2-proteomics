"""Regenerate the English Step 9 phosphosite figure PDFs."""
import argparse
import csv
import json
import os
import re
import subprocess
import tempfile
import time
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode, urlparse
from urllib.request import Request, urlopen

import linear_models
import sequence_logos
import structure_models

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
DEFAULT_OUTPUT = REPO / "Sarya_George" / "reports" / "step09_phosphosite_figures"


def get_public_data(url):
    request = Request(url, headers={"User-Agent": "JE2-phosphosite-figures/1.0"})
    for attempt in range(3):
        try:
            with urlopen(request, timeout=30) as response:
                return response.read()
        except (URLError, TimeoutError) as error:
            temporary = not isinstance(error, HTTPError) or error.code in {429, 500, 502, 503, 504}
            if not temporary or attempt == 2:
                raise
            time.sleep(2 ** (attempt + 1))


def cached_json(path, url):
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8-sig"))
    data = json.loads(get_public_data(url))
    path.write_text(json.dumps(data), encoding="utf-8")
    return data


def public_references(accessions, cache):
    query = " OR ".join("accession:" + accession for accession in accessions)
    params = urlencode({"query": query, "format": "json", "size": len(accessions)})
    uniprot = cached_json(cache / "uniprot-proteins.json",
                          f"https://rest.uniprot.org/uniprotkb/search?{params}")
    entries = {entry["primaryAccession"]: entry for entry in uniprot["results"]}
    if set(entries) != set(accessions):
        raise ValueError("The UniProt cache must contain exactly the Step 9 accessions.")
    for accession, entry in entries.items():
        if entry["organism"]["taxonId"] != 93061:
            raise ValueError(f"{accession}: expected the NCTC 8325 UniProt reference.")
    return entries


def download_annotations(accessions, cache):
    structures = cache / "models"
    structures.mkdir(exist_ok=True)
    for accession in accessions:
        pfam = cached_json(cache / f"pfam-{accession}.json",
                            f"https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot/{accession}/")
        if pfam.get("next"):
            raise ValueError(f"{accession}: Pfam result spans multiple pages.")

        pdb_path = structures / f"{accession}.pdb"
        info_path = structures / f"{accession}.json"
        if info_path.exists():
            model = json.loads(info_path.read_text(encoding="utf-8-sig"))
        else:
            payload = get_public_data(f"https://alphafold.ebi.ac.uk/api/prediction/{accession}")
            entries = json.loads(payload)
            canonical = [entry for entry in entries if
                         (entry.get("modelEntityId") or entry.get("entryId")) == f"AF-{accession}-F1"]
            if len(canonical) != 1:
                raise ValueError(f"{accession}: expected one canonical AlphaFold F1 model.")
            model = canonical[0]
            info_path.write_text(json.dumps(model), encoding="utf-8")

        if not pdb_path.exists():
            url = model["pdbUrl"]
            parsed = urlparse(url)
            if parsed.scheme != "https" or parsed.hostname not in {
                "alphafold.ebi.ac.uk", "www.alphafold.ebi.ac.uk", "storage.googleapis.com",
            }:
                raise ValueError(f"{accession}: unexpected AlphaFold model URL.")
            pdb = get_public_data(url)
            if b"\nATOM  " not in pdb:
                raise ValueError(f"{accession}: the AlphaFold file contains no atom records.")
            pdb_path.write_bytes(pdb)


def reference_map(reference, je2):
    """Align the complete sequences and map 1-based UniProt positions to JE2."""
    if reference == je2:
        return {position: position for position in range(1, len(je2) + 1)}

    gap = -5
    width = len(je2) + 1
    previous = [j * gap for j in range(width)]
    moves = [bytearray(width) for _ in range(len(reference) + 1)]
    for i, a in enumerate(reference, 1):
        current = [i * gap]
        for j, b in enumerate(je2, 1):
            values = (previous[j - 1] + (4 if a == b else -3),
                      previous[j] + gap, current[-1] + gap)
            move = max(range(3), key=values.__getitem__)
            moves[i][j] = move
            current.append(values[move])
        previous = current
    i, j = len(reference), len(je2)
    mapping = {}
    while i and j:
        move = moves[i][j]
        if move == 0:
            mapping[i] = j
            i -= 1
            j -= 1
        elif move == 1:
            i -= 1
        else:
            j -= 1
    return mapping


def make_models(rows, references):
    grouped = {}
    for row in rows:
        accession = row["accession"]
        key = (accession, row["locus_tag"])
        if key not in grouped:
            reference = references[accession]["sequence"]["value"]
            grouped[key] = {
                "accession": accession,
                "name": row["name"],
                "locus_tag": row["locus_tag"],
                "sequence": row["sequence"],
                "reference_sequence": reference,
                "reference_map": reference_map(reference, row["sequence"]),
                "sites": {},
            }
        protein = grouped[key]
        position = int(row["position"])
        site_key = (position, row["residue"])
        protein["sites"].setdefault(site_key, []).append([row["medium"], row["mutant"]])

    models = list(grouped.values())
    for model in models:
        position_in_reference = {je2: ref for ref, je2 in model["reference_map"].items()}
        sites = []
        for (position, residue), groups in sorted(model["sites"].items()):
            ref_position = position_in_reference.get(position)
            reference = model["reference_sequence"]
            if ref_position is None or reference[ref_position - 1] != residue:
                raise ValueError(f'{model["accession"]}: JE2 {residue}{position} has no matching reference residue.')
            context_start = max(position - 5, 1)
            sites.append({
                "position": position,
                "residue": residue,
                "groups": groups,
                "ref_position": ref_position,
                "ref_residue": reference[ref_position - 1],
                "context_start": context_start,
                "context": model["sequence"][context_start - 1:position + 5],
            })
        model["sites"] = sites
    return models


def default_data_dir():
    if os.environ.get("JE2_DATA_DIR"):
        return Path(os.environ["JE2_DATA_DIR"])
    config = REPO / "Sarya_George" / "config" / "data_dir.txt"
    if config.is_file():
        data_dir = config.read_text(encoding="utf-8-sig").splitlines()
        if data_dir and data_dir[0].strip():
            directory = Path(data_dir[0].strip()).expanduser()
            return config.parent.parent / directory
    return REPO / "data"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=default_data_dir())
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--cache-dir", type=Path, help="Save public reference files and structure renders for reuse.")
    parser.add_argument("--rscript", default="Rscript", help="Rscript executable or absolute path.")
    args = parser.parse_args()

    with tempfile.TemporaryDirectory(prefix="je2-step09-") as work:
        work = Path(work)
        hits = work / "hits.tsv"
        subprocess.run(
            [args.rscript, "--vanilla", str(HERE / "extract_proteins.R"), str(args.data_dir), str(hits)],
            check=True
        )
        with hits.open(encoding="utf-8-sig", newline="") as stream:
            rows = list(csv.DictReader(stream, delimiter="\t"))
        accessions = list(dict.fromkeys(row["accession"] for row in rows))
        if len(accessions) != 12 or len(rows) != 14 or any(
            not re.fullmatch(r"[A-Z0-9]{6,10}", accession) for accession in accessions
        ):
            raise ValueError("Expected the 14 original Step 9 hits across 12 UniProt accessions.")
        cache = args.cache_dir or work / "references"
        cache.mkdir(parents=True, exist_ok=True)
        models = make_models(rows, public_references(accessions, cache))
        download_annotations(accessions, cache)

        output = args.output_dir
        output.mkdir(parents=True, exist_ok=True)
        sequence_logos.render(models, output / "JE2_Sequence_Logos_Plus_Minus_10_AA.pdf")
        linear_models.render(models, cache, output / "JE2_Phosphosite_Protein_Models.pdf")
        structure_models.render(models, cache, output / "JE2_Phosphosites_3D_Structures.pdf")
        print(f"Three Step 9 figures saved in {output.resolve()}")


if __name__ == "__main__":
    main()
