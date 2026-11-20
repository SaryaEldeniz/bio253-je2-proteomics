from pathlib import Path
from tempfile import TemporaryDirectory

from reportlab.lib.colors import HexColor
from reportlab.lib.pagesizes import A4, landscape
from reportlab.lib.utils import ImageReader
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfgen.canvas import Canvas


from render_cartoons import render_cartoons


WIDTH, HEIGHT = landscape(A4)
MARGIN = 31
GAP = 16
PANEL_WIDTH = (WIDTH - 2 * MARGIN - GAP) / 2

INK = "#182939"
MUTED = "#4A6070"
LINK = "#126271"
BORDER = "#CBD8DF"
PAPER = "#F5F8FA"
WHITE = "#FFFFFF"

CONFIDENCE = (
    (">90", "#0053D6"),
    ("70 to 90", "#65CBF3"),
    ("50 to <70", "#FFDB13"),
    ("<50", "#FF7D45"),
)
HIGHLIGHT = "#FFE549"


def draw_text(pdf, x, y, label, size=9, color=INK, bold=False, right=False):
    font = "Helvetica-Bold" if bold else "Helvetica"
    pdf.setFillColor(HexColor(color))
    pdf.setFont(font, size)
    draw = pdf.drawRightString if right else pdf.drawString
    draw(x, y, label)
    return pdfmetrics.stringWidth(label, font, size)


def draw_link(pdf, x, y, label, destination, size=8, right=False):
    width = draw_text(pdf, x, y, label, size, LINK, bold=True, right=right)
    start = x - width if right else x
    pdf.linkURL(destination, (start, y - 2, start + width, y + size + 1), thickness=0)


def draw_header(pdf, protein_count, site_count):
    pdf.setFillColor(HexColor(PAPER))
    pdf.rect(0, 0, WIDTH, HEIGHT, stroke=0, fill=1)
    draw_text(pdf, MARGIN, HEIGHT - 35, "Phosphosites in 3D", 25, bold=True)
    draw_text(
        pdf,
        MARGIN,
        HEIGHT - 54,
        "AlphaFold reference models of strain NCTC 8325 showing phosphosites observed in JE2",
        10,
        MUTED,
    )
    draw_text(
        pdf,
        WIDTH - MARGIN,
        HEIGHT - 35,
        f"{protein_count} proteins  ·  {site_count} phosphosites",
        9,
        MUTED,
        bold=True,
        right=True,
    )
    y = HEIGHT - 79
    x = MARGIN + draw_text(pdf, MARGIN, y, "Model confidence (pLDDT)", 8.4, bold=True) + 15
    for label, color in CONFIDENCE:
        pdf.setFillColor(HexColor(color))
        pdf.roundRect(x, y - 1, 11, 7.5, 1.6, stroke=0, fill=1)
        x += 16 + draw_text(pdf, x + 15, y, label, 8.1, MUTED) + 10

    caption = "Yellow highlights mark the phosphosites"
    caption_w = pdfmetrics.stringWidth(caption, "Helvetica", 8.5)
    pdf.setFillColor(HexColor(HIGHLIGHT))
    pdf.setStrokeColor(HexColor(INK))
    pdf.setLineWidth(0.65)
    pdf.circle(WIDTH - MARGIN - caption_w - 8, y + 3, 3.1, stroke=1, fill=1)
    draw_text(pdf, WIDTH - MARGIN, y, caption, 8.5, MUTED, right=True)


def draw_model_image(pdf, x, y, width, height, image):
    """Fit the full render including its existing arrow, label and margins."""
    native_width, native_height = image.getSize()
    scale = min(width / native_width, height / native_height)
    image_width = native_width * scale
    image_height = native_height * scale
    pdf.drawImage(
        image,
        x + (width - image_width) / 2,
        y + (height - image_height) / 2,
        image_width,
        image_height,
        mask="auto",
    )


def draw_panel(pdf, x, model, manifest):
    y = 49
    top = HEIGHT - 96
    right = x + PANEL_WIDTH
    pdf.setFillColor(HexColor(WHITE))
    pdf.setStrokeColor(HexColor(BORDER))
    pdf.setLineWidth(0.7)
    pdf.roundRect(x, y, PANEL_WIDTH, top - y, 7, stroke=1, fill=1)

    accession = model["accession"]
    name = model["name"].replace("HslU--HslV", "HslU/HslV")
    name = name[0].upper() + name[1:]
    draw_text(pdf, x + 12, top - 21, name, 13.1, bold=True)
    draw_text(pdf, x + 12, top - 37, model["locus_tag"], 8.5, MUTED)
    draw_link(
        pdf,
        right - 12,
        top - 37,
        f"AlphaFold DB · {accession}",
        f"https://alphafold.ebi.ac.uk/entry/{accession}",
        8.6,
        right=True,
    )
    draw_model_image(
        pdf,
        x + 5,
        y + 34,
        PANEL_WIDTH - 10,
        top - y - 78,
        ImageReader(str(manifest["image"])),
    )

    sites = model["sites"]
    rendered_sites = manifest["sites"]
    site_info = []
    position_notes = []
    for site, rendered in zip(sites, rendered_sites, strict=True):
        label = f'{site["residue"]}{site["position"]}'
        reference = f'{site["ref_residue"]}{site["ref_position"]}'
        if rendered["position"] != site["position"] or rendered["ref_position"] != site["ref_position"]:
            raise ValueError(f"Render residue numbers do not match the site data for {accession} {label}.")
        score = float(rendered["plddt"])
        if not 0.0 <= score <= 100.0:
            raise ValueError(f"Invalid pLDDT for {accession} {label}")
        site_info.append(f"{label}  {score:.1f}")
        if label != reference:
            position_notes.append(f"JE2 {label} · UniProt {reference}")

    model_id = manifest["model_id"]
    draw_text(pdf, x + 12, y + 10, model_id, 7.6, MUTED)
    scores = "Local pLDDT   " + "   /   ".join(site_info)
    draw_text(pdf, right - 12, y + 10, scores, 8.0, MUTED, right=True)
    if position_notes:
        draw_text(pdf, x + 12, y + 24, "   /   ".join(position_notes), 9.5, LINK, bold=True)


def draw_footer(pdf, page, total):
    pdf.setStrokeColor(HexColor(BORDER))
    pdf.setLineWidth(0.7)
    pdf.line(MARGIN, 35, WIDTH - MARGIN, 35)
    x = MARGIN
    x += draw_text(
        pdf,
        x,
        21,
        "Arrow labels use JE2 residue numbers.",
        7.8,
        MUTED,
    ) + 12
    draw_link(
        pdf,
        x,
        21,
        "AlphaFold DB",
        "https://alphafold.ebi.ac.uk/",
        7.8,
    )
    draw_text(pdf, WIDTH - MARGIN, 21, f"{page} / {total}", 8, MUTED, bold=True, right=True)


def render(models, cache, output):
    """Render the selected proteins and assemble their 3D structure report.

    The cache's models directory must contain each UniProt accession's PDB
    coordinates and its AlphaFold metadata JSON. The report is written at output.
    """
    cache = Path(cache)
    output = Path(output)
    if not models or len({model["accession"] for model in models}) != len(models):
        raise ValueError("Expected a non-empty list of proteins with unique UniProt accessions.")
    manifest = render_cartoons(models, cache)
    site_count = sum(len(model["sites"]) for model in models)
    pages = (len(models) + 1) // 2

    output.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(dir=output.parent) as temporary_dir:
        working_copy = Path(temporary_dir) / output.name
        pdf = Canvas(str(working_copy), pagesize=(WIDTH, HEIGHT), pageCompression=1)
        pdf.setTitle("Phosphosites in 3D | Staphylococcus aureus JE2")
        pdf.setSubject("AlphaFold reference models with annotated JE2 phosphosites")
        for page in range(pages):
            draw_header(pdf, len(models), site_count)
            for column, model in enumerate(models[2 * page:2 * page + 2]):
                x = MARGIN + column * (PANEL_WIDTH + GAP)
                draw_panel(pdf, x, model, manifest[model["accession"]])
            draw_footer(pdf, page + 1, pages)
            pdf.showPage()
        pdf.save()
        working_copy.replace(output)
