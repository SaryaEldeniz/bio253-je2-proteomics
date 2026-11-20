"""Render linear maps of the JE2 phosphosite results."""

import json
import math
from pathlib import Path

from reportlab.lib.colors import HexColor
from reportlab.lib.pagesizes import A4, landscape
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfgen.canvas import Canvas

W, H = landscape(A4)
BODY = "Helvetica"
BOLD = "Helvetica-Bold"
MONO = "Courier"
MONO_BOLD = "Courier-Bold"

INK = "#172637"
MUTED = "#4C6171"
PAPER = "#F4F7F9"
RULE = "#CFDCE3"
AXIS = "#79929F"
TEAL = "#006B73"
ORANGE = "#9F492B"
MUTANT_COLOR = {"stk": TEAL, "stp": ORANGE}

FEATURE_STYLES = {
    "Signal": ("Signal peptide", "#F2D58D", "#69480D"),
    "Chain": ("Mature chain", "#DAE4EB", "#41576B"),
    "Domain": ("Domain", "#CFE8D9", "#1F5C3F"),
    "Transmembrane": ("TM", "#D8CCF0", "#52367D"),
    "Region": ("Region", "#D6E4ED", "#33576F"),
    "Coiled coil": ("Coiled coil", "#E1D8EC", "#55426D"),
    "Helix": ("Alpha helix", "#BEE1EA", "#15566A"),
    "Beta strand": ("Beta strand", "#D8E6A7", "#45581C"),
    "Turn": ("Turn", "#D5E7E9", "#33576F"),
}
POSITIONAL_SITES = {"Binding site", "Site", "Active site", "Modified residue", "Cross-link", "Disulfide bond"}


def color(value):
    return HexColor(value)


def text(draw, x, y, value, size=9, tone=INK, face=BODY, align="left"):
    draw.setFillColor(color(tone))
    draw.setFont(face, size)
    add = {"left": draw.drawString, "right": draw.drawRightString, "center": draw.drawCentredString}[align]
    add(x, y, value)


def linked_text(draw, x, y, value, url, size=9, align="left"):
    text(draw, x, y, value, size, TEAL, BOLD, align)
    width = pdfmetrics.stringWidth(value, BOLD, size)
    left = x - width if align == "right" else x
    draw.linkURL(url, (left, y - 2, left + width, y + size), relative=0)


def segment(draw, x0, x1, y, fill, label, tone, height=8):
    draw.setFillColor(color(fill))
    draw.roundRect(x0, y, max(1.5, x1 - x0), height, min(2, height / 3), stroke=0, fill=1)
    if label and pdfmetrics.stringWidth(label, BOLD, 6.1) + 6 <= x1 - x0:
        text(draw, (x0 + x1) / 2, y + (height - 6.1) / 2 + 0.7, label, 6.1, tone, BOLD, "center")


def project(ref_to_je2, start, end):
    matched = [ref_to_je2[pos] for pos in range(start, end + 1) if pos in ref_to_je2]
    if matched:
        return matched[0], matched[-1]
    return None


def feature_tracks(protein, ref, cache):
    ref_to_je2 = {int(key): value for key, value in protein["reference_map"].items()}
    structural = []
    key_sites = set()
    for feature in ref.get("features", []):
        kind = feature["type"]
        if kind not in FEATURE_STYLES and kind not in POSITIONAL_SITES:
            continue
        start = feature["location"]["start"].get("value")
        stop = feature["location"]["end"].get("value")
        if start is None or stop is None or (kind == "Chain" and start == 1):
            continue
        mapped = project(ref_to_je2, start, stop)
        if not mapped:
            continue
        if kind in POSITIONAL_SITES:
            key_sites.update(mapped)
            continue
        base_label, fill, tone = FEATURE_STYLES[kind]
        label = feature["description"] if kind == "Domain" and feature.get("description") else base_label
        if kind == "Region" and feature.get("description") == "Disordered":
            label = "Disordered"
        structural.append((*mapped, label, fill, tone))
    domains = []
    pfam_names = {
        match["id"]: next(
            (item["value"] for item in match.get("properties", []) if item["key"] == "EntryName"),
            ""
        )
        for match in ref.get("uniProtKBCrossReferences", [])
        if match["database"] == "Pfam"
    }
    domain_file = cache / f'pfam-{protein["accession"]}.json'
    if domain_file.exists() and domain_file.stat().st_size:
        payload = json.loads(domain_file.read_text(encoding="utf-8-sig"))
        assert payload["next"] is None
        for match in payload.get("results", []):
            accession = match["metadata"]["accession"]
            for entry in match["proteins"]:
                assert entry["accession"].upper() == protein["accession"]
                assert entry["protein_length"] == len(protein["reference_sequence"])
                for location in entry["entry_protein_locations"]:
                    for fragment in location["fragments"]:
                        mapped = project(ref_to_je2, fragment["start"], fragment["end"])
                        if mapped:
                            domains.append((*mapped, accession, pfam_names.get(accession, "")))
    return sorted(structural), sorted(domains), sorted(key_sites)


def ticks(length):
    step = next(item for item in (10, 20, 25, 50, 100, 200) if item >= length / 7)
    inner = [value for value in range(step, length, step) if length - value > step / 4]
    return [1, *inner, length]


def group_text(groups):
    media = " + ".join(medium for medium, _ in groups)
    mutants = {mutant for _, mutant in groups}
    assert len(mutants) == 1
    return f"{media} · {groups[0][1]}"


def arrow(draw, position_x, lower_y, box_y, x_start, x_end, site):
    label = f'{site["residue"]}{site["position"]}'
    group = group_text(site["groups"])
    tone = MUTANT_COLOR[site["groups"][0][1]]
    label_width = pdfmetrics.stringWidth(label, BOLD, 9)
    group_width = pdfmetrics.stringWidth(group, BODY, 8.5)
    width = 22 + label_width + group_width
    left = max(x_start - 3, min(position_x - width / 2, x_end - width + 3))
    draw.setFillColor(color("#E0F0F1" if tone == TEAL else "#F7E9DF"))
    draw.roundRect(left, box_y, width, 15, 4, stroke=0, fill=1)
    text(draw, left + 7, box_y + 4, label, 9, tone, BOLD)
    text(draw, left + 14 + label_width, box_y + 4, group, 8.5, tone)
    draw.setStrokeColor(color(tone))
    draw.setLineWidth(1.18)
    draw.line(left + width / 2, box_y, position_x, lower_y + 5)
    tip = draw.beginPath()
    tip.moveTo(position_x, lower_y)
    tip.lineTo(position_x - 3.1, lower_y + 5.5)
    tip.lineTo(position_x + 3.1, lower_y + 5.5)
    tip.close()
    draw.setFillColor(color(tone))
    draw.drawPath(tip, stroke=0, fill=1)


def context(draw, x, y, site):
    sequence = site["context"]
    index = site["position"] - site["context_start"]
    assert sequence[index] == site["residue"]
    label = f'{site["residue"]}{site["position"]}'
    tone = MUTANT_COLOR[site["groups"][0][1]]
    text(draw, x, y, label, 8.2, tone, BOLD)
    cursor = x + 31
    glyph = pdfmetrics.stringWidth("A", MONO, 9)
    for i, aa in enumerate(sequence):
        if i == index:
            draw.setFillColor(color(tone))
            draw.roundRect(cursor - 1.3, y - 2, glyph + 2.6, 12, 2, fill=1, stroke=0)
        text(draw, cursor, y, aa, 9, "#FFFFFF" if i == index else INK, MONO_BOLD if i == index else MONO)
        cursor += glyph + 1.1
    ref_label = f'UniProt {site["ref_residue"]}{site["ref_position"]}'
    differs = site["position"] != site["ref_position"] or site["residue"] != site["ref_residue"]
    text(draw, cursor + 10, y, ref_label, 8.2, INK if differs else MUTED, BOLD if differs else BODY)


def reference_features(draw, x_start, x_end, y, x_at, structural, key_sites):
    draw.setStrokeColor(color(RULE))
    draw.setLineWidth(0.5)
    draw.line(x_start, y + 3.8, x_end, y + 3.8)
    end_by_lane = []
    for start, end, label, fill, tone in structural:
        lane = next((i for i, previous in enumerate(end_by_lane) if previous < start), len(end_by_lane))
        if lane == len(end_by_lane):
            end_by_lane.append(end)
        else:
            end_by_lane[lane] = end
        lane_y = y + 5 if lane == 0 and len(structural) > 1 else y
        if lane > 0:
            lane_y = y - (lane - 1) * 3.5
        segment(draw, max(x_start, x_at(start) - 0.7), min(x_end, x_at(end) + 0.7), lane_y, fill, label, tone, 8 if lane == 0 else 4.5)
    draw.setStrokeColor(color("#8198AC"))
    draw.setLineWidth(0.72)
    for pos in key_sites:
        x = x_at(pos)
        draw.line(x, y + 1, x, y + 5.6)


def protein_panel(draw, top, protein, ref, cache):
    left, right = 38, W - 38
    draw.setFillColor(color("#FFFFFF"))
    draw.setStrokeColor(color(RULE))
    draw.setLineWidth(0.7)
    draw.roundRect(left, top - 144, right - left, 144, 8, stroke=1, fill=1)

    accession = protein["accession"]
    sequence = protein["sequence"]
    n = len(sequence)
    name = protein["name"].replace("HslU--HslV", "HslU/HslV")
    name = name[0].upper() + name[1:]
    text(draw, left + 15, top - 20, name, 14, face=BOLD)
    text(draw, left + 15, top - 34, f'{protein["locus_tag"]}  ·  {n} AA', 9, MUTED)
    linked_text(draw, right - 15, top - 20, f"UniProt {accession}", f"https://www.uniprot.org/uniprotkb/{accession}/feature-viewer", 10, "right")
    text(draw, right - 15, top - 34, "NCTC 8325 reference strain", 8.5, MUTED, align="right")

    x_start, x_end = left + 120, right - 25

    def x_at(pos):
        return x_start + (pos - 1) / (n - 1) * (x_end - x_start)

    axis_y = top - 68
    text(draw, left + 15, axis_y - 3, "JE2 position", 8.3, MUTED, BOLD)
    draw.setStrokeColor(color(AXIS))
    draw.setLineCap(1)
    draw.setLineWidth(2.6)
    draw.line(x_start, axis_y, x_end, axis_y)
    draw.setLineWidth(0.65)
    for pos in ticks(n):
        x = x_at(pos)
        draw.line(x, axis_y, x, axis_y - 5)
        text(draw, x, axis_y - 14, str(pos), 7.4, MUTED, align="center")
    for site in protein["sites"]:
        arrow(draw, x_at(site["position"]), axis_y + 0.8, top - 55, x_start, x_end, site)

    structural, domains, other_sites = feature_tracks(protein, ref, cache)
    text(draw, left + 15, top - 97, "UniProt features", 8.3, MUTED, BOLD)
    reference_features(draw, x_start, x_end, top - 104, x_at, structural, other_sites)

    text(draw, left + 15, top - 116, "Pfam", 8.3, MUTED, BOLD)
    if domains:
        for start, end, domain, name in domains:
            start_x, end_x = x_at(start), x_at(end)
            full_label = f"{name} · {domain}" if name else domain
            label = full_label if pdfmetrics.stringWidth(full_label, BOLD, 6.1) + 6 <= end_x - start_x else domain
            segment(draw, start_x, end_x, top - 119, "#CFE8D9", label, "#1F5C3F")
            draw.linkURL(f"https://www.ebi.ac.uk/interpro/entry/pfam/{domain}/", (start_x, top - 119, end_x, top - 111), relative=0)
    else:
        pfam = [
            item["id"] for item in ref.get("uniProtKBCrossReferences", [])
            if item["database"] == "Pfam"
        ]
        if pfam:
            for index, domain in enumerate(pfam):
                linked_text(draw, x_start + index * 53, top - 116, domain, f"https://www.ebi.ac.uk/interpro/entry/pfam/{domain}/", 8)
            text(draw, x_start + len(pfam) * 53 + 5, top - 116, "Database cross-references", 8.2, MUTED)

    text(draw, left + 15, top - 132, "JE2 sequence", 8.3, MUTED, BOLD)
    for index, site in enumerate(protein["sites"]):
        context(draw, x_start + 262 * index, top - 132, site)
    return bool(domains)


def page_header(draw):
    draw.setFillColor(color(PAPER))
    draw.rect(0, 0, W, H, fill=1, stroke=0)
    text(draw, 38, H - 38, "Phosphosites in S. aureus JE2", 24, face=BOLD)
    text(draw, 38, H - 58, "12 linear protein models · 13 phosphosites · 14 hits across media and mutants", 10.5, MUTED)

    legend_y = H - 78
    draw.setFillColor(color(TEAL))
    draw.circle(43, legend_y + 2.4, 3.2, stroke=0, fill=1)
    text(draw, 52, legend_y, "stk vs. anc", 8.8, INK, BOLD)
    draw.setFillColor(color(ORANGE))
    draw.circle(118, legend_y + 2.4, 3.2, stroke=0, fill=1)
    text(draw, 127, legend_y, "stp vs. anc", 8.8, INK, BOLD)
    text(draw, 199, legend_y, "Arrows mark phosphosites on the JE2 sequence.", 8.8, MUTED)
    text(draw, W - 38, legend_y, "Thin ticks mark other UniProt sites.", 8.8, MUTED, align="right")


def page_footer(draw, page, page_count):
    draw.setStrokeColor(color(RULE))
    draw.setLineWidth(0.6)
    draw.line(38, 34, W - 38, 34)
    text(draw, 38, 22, "UniProt positions and features were mapped to the JE2 sequence.", 7.7, MUTED)
    linked_text(draw, 367, 22, "STYH analysis", "https://github.com/SaryaEldeniz/bio253-je2-proteomics/blob/ff9a6546437974b2834cd182ba2c376403019bb4/Sarya_George/reports/step09_exclusive_phosphosite_proteins.md", 7.7)
    linked_text(draw, 450, 22, "UniProt feature types", "https://www.uniprot.org/help/sequence_annotation", 7.7)
    linked_text(draw, 546, 22, "Pfam / InterPro", "https://www.ebi.ac.uk/interpro/entry/pfam/", 7.7)
    text(draw, W - 38, 21, f"{page} / {page_count}", 8, MUTED, BOLD, "right")


def render(models, cache: Path, output: Path):
    """Write the model atlas using the cached UniProt and Pfam annotations."""
    database = json.loads((cache / "uniprot-proteins.json").read_text(encoding="utf-8-sig"))
    references = {entry["primaryAccession"]: entry for entry in database["results"]}
    if len(models) != 12 or {protein["accession"] for protein in models} != references.keys():
        raise ValueError("Expected one UniProt reference for each of the 12 JE2 protein models.")

    site_count = sum(len(protein["sites"]) for protein in models)
    group_count = sum(len(site["groups"]) for protein in models for site in protein["sites"])
    if (site_count, group_count) != (13, 14):
        raise ValueError("Expected 13 distinct phosphosites and 14 occurrences by medium and mutant.")

    for protein in models:
        accession = protein["accession"]
        reference_sequence = references[accession]["sequence"]["value"]
        if protein["reference_sequence"] != reference_sequence:
            raise ValueError(f"The alignment and cached reference sequence differ for {accession}.")
        for site in protein["sites"]:
            if protein["sequence"][site["position"] - 1] != site["residue"]:
                raise ValueError(f"JE2 residue mismatch at {accession} position {site['position']}.")
            if reference_sequence[site["ref_position"] - 1] != site["ref_residue"]:
                raise ValueError(f"UniProt residue mismatch at {accession} position {site['ref_position']}.")

    output.parent.mkdir(exist_ok=True, parents=True)
    draw = Canvas(str(output), pagesize=(W, H), pageCompression=1, invariant=1)
    draw.setTitle("Phosphosites in Staphylococcus aureus JE2")
    draw.setSubject("Linear models of twelve proteins with exclusive observed STYH hits")
    draw.setCreator("BIO253 JE2 proteomics figure pipeline")
    with_domains = 0
    page_count = math.ceil(len(models) / 3)
    for offset in range(0, len(models), 3):
        page_header(draw)
        for row, protein in enumerate(models[offset:offset + 3]):
            reference = references[protein["accession"]]
            with_domains += protein_panel(draw, H - 94 - 151 * row, protein, reference, cache)
        page_footer(draw, offset // 3 + 1, page_count)
        draw.showPage()
    draw.save()
    print(
        f"{page_count} pages, {len(models)} proteins, {site_count} sites, "
        f"{group_count} group occurrences. Pfam coordinate tracks: {with_domains}."
    )
