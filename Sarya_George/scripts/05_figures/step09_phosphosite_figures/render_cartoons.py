import hashlib
import json
import math
import os
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
from PIL import Image, ImageDraw, ImageFont

SIDE = 1100
PAPER = (244, 247, 249, 255)
DARK = (21, 35, 53, 255)
COLORS = [
    ("af_dark_blue", "0053D6"), ("af_light_blue", "65CBF3"),
    ("af_yellow", "FFDB13"), ("af_orange", "FF7D45"),
    ("selected_yellow", "FFF000"),
]
AA3 = dict(zip(
    "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split(),
    "ARNDCQEGHILKMFPSTWYV"
))


def label_font(size, bold=False):
    windows_fonts = Path(os.environ.get("WINDIR", "C:/Windows")) / "Fonts"
    arial = "arialbd.ttf" if bold else "arial.ttf"
    mac_arial = "Arial Bold.ttf" if bold else "Arial.ttf"
    dejavu = "DejaVuSans-Bold.ttf" if bold else "DejaVuSans.ttf"
    for candidate in (windows_fonts / arial, Path("/System/Library/Fonts/Supplemental") / mac_arial, dejavu):
        try:
            return ImageFont.truetype(str(candidate), size)
        except OSError:
            continue
    return ImageFont.load_default(size=size)


def read_model(protein, pdb_bytes):
    accession = protein["accession"]
    pdb = pdb_bytes.decode("ascii").splitlines()
    atoms = [line for line in pdb if line.startswith("ATOM  ") and line[12:16].strip() == "CA"]
    residues = {int(line[22:26]): line for line in atoms}
    positions = list(range(1, len(protein["reference_sequence"]) + 1))
    if len(atoms) != len(residues) or sorted(residues) != positions or {line[21] for line in atoms} != {"A"}:
        raise ValueError(f"{accession}: expected one alpha carbon per residue in chain A.")
    records = [residues[position] for position in positions]
    reference = "".join(AA3[line[17:20]] for line in records)
    if reference != protein["reference_sequence"]:
        raise ValueError(f"{accession}: PDB and UniProt sequences differ.")

    ca = np.array([[float(line[i:i + 8]) for i in (30, 38, 46)] for line in records])
    site_indices = [site["ref_position"] - 1 for site in protein["sites"]]
    for site, index in zip(protein["sites"], site_indices, strict=True):
        residue = site["residue"]
        je2_index = site["position"] - 1
        if (
            not 0 <= index < len(reference)
            or reference[index] != site["ref_residue"]
            or not 0 <= je2_index < len(protein["sequence"])
            or protein["sequence"][je2_index] != residue
        ):
            raise ValueError(f'{accession}: sequence mismatch at JE2 {residue}{site["position"]}.')
    scores = [float(records[index][60:66]) for index in site_indices]
    return pdb, ca, site_indices, scores


def choose_views(ca, indices):
    coords = ca - ca.mean(axis=0)
    views = []
    ref_up = np.array([0.0, 1.0, 0.0])
    for sample in range(220):
        z = 1 - 2 * (sample + 0.5) / 220
        angle = sample * math.pi * (3 - math.sqrt(5))
        radial = math.sqrt(1 - z * z)
        facing = np.array([radial * math.cos(angle), radial * math.sin(angle), z])
        horizontal = np.cross(ref_up, facing)
        horizontal /= np.linalg.norm(horizontal)
        vertical = np.cross(facing, horizontal)
        view = np.stack([horizontal, vertical, facing])
        projected = coords @ view.T
        covariance = np.cov(projected[:, :2].T)
        _, axes = np.linalg.eigh(covariance)
        u, v = axes[:, 1]
        twist = np.array([[u, v, 0], [-v, u, 0], [0, 0, 1]])
        view = twist @ view
        projected = coords @ view.T
        size = np.ptp(projected[:, :2], axis=0)
        sites = projected[indices]
        visibility = 0.0
        for target_index, target in zip(indices, sites):
            delta = projected - target
            near = np.linalg.norm(delta[:, :2], axis=1) < 4.0
            in_front = delta[:, 2] > 1.7
            window = np.abs(np.arange(ca.shape[0]) - target_index) > 1
            visibility -= 2.2 * np.count_nonzero(near & in_front & window)
            visibility += 0.9 * np.linalg.norm(target[:2] / size)
        if len(sites) > 1:
            visibility += min(0.8, 3.4 * np.linalg.norm((sites[0, :2] - sites[1, :2]) / size))
        visibility -= 0.4 * abs(math.log(size[0] / size[1]))
        views.append((visibility, view))
    views.sort(key=lambda item: item[0], reverse=True)
    return [matrix for _, matrix in views]


def rotated_pdb(lines, origin, matrix):
    output = []
    for line in lines:
        if line.startswith("ATOM  "):
            xyz = np.array([float(line[i:i + 8]) for i in (30, 38, 46)])
            xyz = matrix @ (xyz - origin)
            line = line[:30] + "".join(f"{value:8.3f}" for value in xyz) + line[54:]
        output.append(line)
    return "\n".join(output) + "\n"


def png(cmd, path):
    cmd.png(str(path), width=SIDE, height=SIDE, ray=1, quiet=1)
    with Image.open(path) as rendered:
        return rendered.convert("RGBA")


def install_style(cmd):
    for name, hex_color in COLORS:
        cmd.set_color(name, [int(hex_color[i:i + 2], 16) / 255 for i in (0, 2, 4)])
    cmd.set("orthoscopic", 1)
    cmd.set("ray_opaque_background", 0)
    cmd.set("ray_shadows", 0)
    cmd.set("ambient", 0.48)
    cmd.set("direct", 0.57)
    cmd.set("specular", 0.20)
    cmd.set("shininess", 35)
    cmd.set("depth_cue", 0)
    cmd.set("antialias", 2)
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("cartoon_smooth_loops", 1)
    cmd.set("cartoon_loop_radius", 0.26)
    cmd.set("cartoon_tube_radius", 0.30)
    cmd.viewport(SIDE, SIDE)


def render_view(cmd, pdb, ca, indices, matrix, accession, scratch):
    cmd.delete("all")
    cmd.read_pdbstr(rotated_pdb(pdb, ca.mean(axis=0), matrix), "reference", quiet=1)
    cmd.hide("everything", "all")
    cmd.dss("reference", quiet=1)
    cmd.show("cartoon", "reference")
    cmd.color("af_orange", "reference")
    cmd.color("af_yellow", "reference and b >= 50")
    cmd.color("af_light_blue", "reference and b >= 70")
    cmd.color("af_dark_blue", "reference and b > 90")
    markers = []
    for i, sequence_index in enumerate(indices):
        resi = sequence_index + 1
        selector = f"reference and chain A and resi {resi}"
        cmd.color("selected_yellow", selector)
        cmd.set("cartoon_color", "selected_yellow", selector)
        cmd.show("sticks", selector)
        cmd.set("stick_radius", 0.20, selector)
        coord = cmd.get_atom_coords(f"{selector} and name CA")
        marker = f"site_marker_{i}"
        cmd.pseudoatom(marker, pos=coord, elem="C", vdw=0.90, color="selected_yellow")
        cmd.show("spheres", marker)
        markers.append(marker)

    cmd.reset()
    extent = np.ptp((ca - ca.mean(axis=0)) @ matrix.T, axis=0)
    cmd.zoom("reference", buffer=max(4.0, 0.11 * max(extent[:2])), complete=1)
    view = cmd.get_view()
    raw = png(cmd, scratch / f"{accession}-raw.png")
    rgb = np.asarray(raw)[:, :, :3]
    visibly_yellow = (rgb[:, :, 0] > 155) & (rgb[:, :, 1] > 140) & (rgb[:, :, 2] < 75)

    centres = []
    quality = []
    for i, marker in enumerate(markers):
        cmd.disable("all")
        cmd.enable(marker)
        cmd.set_view(view)
        points = png(cmd, scratch / f"{accession}-point-{i}.png")
        alpha = np.asarray(points)[:, :, 3] > 150
        y, x = np.nonzero(alpha)
        if not len(x):
            raise ValueError(f"{accession}: isolated site marker missing in ray render")
        centres.append((float(x.mean()), float(y.mean())))
        quality.append(float(visibly_yellow[alpha].mean()))
    return raw, centres, min(quality)


def closest(box, anchor):
    a, b, c, d = box
    x, y = anchor
    return np.array([min(c, max(a, x)), min(d, max(b, y))])


def label_options(draw, anchor, texts, mask, anchors, fonts):
    title, detail = texts
    bold, regular = fonts
    width = max(draw.textlength(title, bold), draw.textlength(detail, regular)) + 42
    height = 97
    x, y = anchor
    options = []
    for radius in (150, 205, 255, 320):
        for theta in np.linspace(0, 2 * math.pi, 24, endpoint=False):
            cx = x + math.cos(theta) * radius
            cy = y + math.sin(theta) * radius
            left = int(max(24, min(SIDE - width - 24, cx - width / 2)))
            top = int(max(24, min(SIDE - height - 24, cy - height / 2)))
            box = (left, top, int(left + width), top + height)
            start = closest(box, anchor)
            distance = np.linalg.norm(start - anchor)
            if distance < 56:
                continue
            if any(box[0] - 15 < px < box[2] + 15 and box[1] - 15 < py < box[3] + 15 for px, py in anchors):
                continue
            cutout = mask[box[1]:box[3], box[0]:box[2]]
            line_probe = np.linspace(start, anchor, 30)[:-4].astype(int)
            crosses = mask[line_probe[:, 1], line_probe[:, 0]].mean()
            score = 1500 * cutout.mean() + 250 * crosses + distance * 0.35
            options.append((score, box))
    if not options:
        raise RuntimeError("No clear label placement around the phosphosite")
    return sorted(options)


def boxes_conflict(first, second, centres):
    a, b, c, d = first
    e, f, g, h = second
    if a < g + 18 and c > e - 18 and b < h + 18 and d > f - 18:
        return True
    lines = [np.linspace(closest(box, anchor), anchor, 27)
             for box, anchor in zip((first, second), centres)]
    for line, box in zip(lines, (second, first)):
        x, y = line.T
        if np.any((box[0] - 10 < x) & (x < box[2] + 10) & (box[1] - 10 < y) & (y < box[3] + 10)):
            return True
    spacing = np.linalg.norm(lines[0][:, None, :] - lines[1][None, :, :], axis=2).min()
    return spacing < 24


def choose_boxes(draw, centres, labels, mask, fonts):
    candidates = [label_options(draw, coord, label, mask, centres, fonts)
                  for coord, label in zip(centres, labels)]
    if len(candidates) == 1:
        return [candidates[0][0][1]]
    best = (math.inf, None)
    for first_cost, first in candidates[0]:
        if first_cost + candidates[1][0][0] >= best[0]:
            break
        for second_cost, second in candidates[1]:
            if first_cost + second_cost >= best[0]:
                break
            if not boxes_conflict(first, second, centres):
                best = (first_cost + second_cost, [first, second])
                break
    if best[1] is None:
        raise RuntimeError("The two residue arrows could not be placed without overlap")
    return best[1]


def add_labels(raw, centres, protein, fonts):
    bold, regular = fonts
    alpha = np.asarray(raw)[:, :, 3] > 190
    figure = Image.new("RGBA", raw.size, PAPER)
    figure.alpha_composite(raw)
    pen = ImageDraw.Draw(figure)

    labels = []
    for site in protein["sites"]:
        title = f'{site["residue"]}{site["position"]}'
        groups = site["groups"]
        if len({mutant for _, mutant in groups}) == 1:
            media = " + ".join(medium for medium, _ in groups)
            detail = f"{media} · {groups[0][1]}"
        else:
            detail = " / ".join(f"{medium} · {mutant}" for medium, mutant in groups)
        labels.append((title, detail))
    boxes = choose_boxes(pen, centres, labels, alpha, fonts)
    for centre, box, (title, detail) in zip(centres, boxes, labels):
        anchor = np.array(centre)
        start = closest(box, anchor)
        vector = (anchor - start) / np.linalg.norm(anchor - start)
        edge = anchor - 17 * vector
        perpendicular = np.array([-vector[1], vector[0]])
        head = [tuple(edge), tuple(edge - 21 * vector + 10 * perpendicular),
                tuple(edge - 21 * vector - 10 * perpendicular)]

        x, y = centre
        pen.ellipse((x - 15, y - 15, x + 15, y + 15), outline=(255, 255, 255, 255), width=8)
        pen.ellipse((x - 15, y - 15, x + 15, y + 15), outline=DARK, width=3)
        pen.line([tuple(start), tuple(edge - 13 * vector)], fill=(255, 255, 255, 255), width=11)
        pen.line([tuple(start), tuple(edge - 13 * vector)], fill=DARK, width=5)
        pen.polygon(head, fill=DARK)
        pen.rounded_rectangle(box, radius=15, fill=(255, 255, 255, 255), outline=(158, 177, 190, 255), width=2)
        pen.text((box[0] + 19, box[1] + 12), title, font=bold, fill=DARK)
        pen.text((box[0] + 19, box[1] + 57), detail, font=regular, fill=(68, 88, 104, 255))

    outline = [raw.getchannel("A").getbbox(), *boxes]
    frame = (
        max(0, min(box[0] for box in outline) - 42),
        max(0, min(box[1] for box in outline) - 42),
        min(SIDE, max(box[2] for box in outline) + 42),
        min(SIDE, max(box[3] for box in outline) + 42),
    )
    return figure.crop(frame).convert("RGB"), boxes, frame


def render_protein(cmd, protein, pdb_bytes, metadata, scratch, fonts):
    accession = protein["accession"]
    model_id = metadata.get("modelEntityId") or metadata.get("entryId")
    if not model_id or accession not in model_id:
        raise ValueError(f"{accession}: AlphaFold metadata does not match the UniProt accession.")
    if not 1 <= len(protein["sites"]) <= 2:
        raise ValueError(f"{accession}: the cartoon layout requires one or two sites.")

    pdb, ca, indices, scores = read_model(protein, pdb_bytes)
    best_quality = 0.0
    best_view = None
    for matrix in choose_views(ca, indices)[:5]:
        image, centres, quality = render_view(cmd, pdb, ca, indices, matrix, accession, scratch)
        if quality > best_quality:
            best_view = image, centres
            best_quality = quality
        if quality >= 0.58:
            break
    if best_quality < 0.40:
        raise RuntimeError(f"{accession}: marker visibility too low ({best_quality:.2f}).")

    figure, boxes, frame = add_labels(*best_view, protein, fonts)
    origin_x, origin_y = frame[:2]
    return figure, {
        "model_id": model_id,
        "structure_url": metadata["pdbUrl"],
        "sites": [
            {
                "position": site["position"],
                "ref_position": site["ref_position"],
                "plddt": score,
                "pixel": [pixel[0] - origin_x, pixel[1] - origin_y],
                "label_box": [box[0] - origin_x, box[1] - origin_y, box[2] - origin_x, box[3] - origin_y],
            }
            for site, score, pixel, box in zip(protein["sites"], scores, best_view[1], boxes, strict=True)
        ],
    }


def render_cartoons(models, cache):
    """Render from cached PDBs and reuse figures only if all inputs match.

    Render entries in cache/renders include a hash of the model and site data,
    PDB and metadata files, and this renderer source. The entry is published
    only once both the PNG and its accompanying metadata have been generated.
    """
    model_dir = cache / "models"
    output = cache / "renders"
    output.mkdir(parents=True, exist_ok=True)
    renderer_source = Path(__file__).read_bytes()
    manifest = {}
    pending = []
    for protein in models:
        accession = protein["accession"]
        pdb_bytes = (model_dir / f"{accession}.pdb").read_bytes()
        metadata_bytes = (model_dir / f"{accession}.json").read_bytes()
        digest = hashlib.sha256(renderer_source)
        digest.update(pdb_bytes)
        digest.update(metadata_bytes)
        digest.update(json.dumps(protein, sort_keys=True).encode("utf-8"))
        entry = output / f"{accession}-{digest.hexdigest()}"
        image = entry.with_suffix(".png")
        info = entry.with_suffix(".json")
        if image.exists() and info.exists():
            manifest[accession] = json.loads(info.read_text(encoding="utf-8"))
            manifest[accession]["image"] = image
        else:
            pending.append((protein, pdb_bytes, json.loads(metadata_bytes), image, info))

    if pending:
        from pymol2 import PyMOL

        fonts = (label_font(42, bold=True), label_font(29))
        with TemporaryDirectory(dir=cache) as scratch_dir, PyMOL() as pymol:
            scratch = Path(scratch_dir)
            cmd = pymol.cmd
            install_style(cmd)
            for protein, pdb_bytes, metadata, image, info in pending:
                figure, details = render_protein(cmd, protein, pdb_bytes, metadata, scratch, fonts)
                temporary_image = scratch / image.name
                temporary_info = scratch / info.name
                figure.save(temporary_image, optimize=True)
                temporary_info.write_text(json.dumps(details, indent=2) + "\n", encoding="utf-8")
                temporary_image.replace(image)
                temporary_info.replace(info)
                manifest[protein["accession"]] = {**details, "image": image}

    return manifest
