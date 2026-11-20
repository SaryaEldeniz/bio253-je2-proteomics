from collections import Counter
from pathlib import Path

from reportlab.lib.colors import HexColor
from reportlab.lib.pagesizes import A4, landscape
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfgen.canvas import Canvas

W, H = landscape(A4)
MARGIN = 32
PAPER = "#F3F7F8"
INK = "#152839"
MUTED = "#496173"
RULE = "#CFDEE4"
COLUMN = "#EDF1FB"
FOCUS = "#313D8C"
TEAL = "#00646D"
GAP_COLOR = "#965231"
COLORS = {
    **dict.fromkeys("GSTYC", "#18704E"),
    **dict.fromkeys("NQ", "#764897"),
    **dict.fromkeys("KRH", "#1D5293"),
    **dict.fromkeys("DE", "#A7394D"),
    **dict.fromkeys("AVLIPWMF", "#34485A"),
}
OFFSETS = tuple(range(-10, 11))
GROUPS = (("TSB", "stk"), ("TSB", "stp"), ("SN", "stk"), ("SN", "stp"))
GROUP_COUNTS = (3, 4, 3, 4)
ANALYSIS_URL = (
    "https://github.com/SaryaEldeniz/bio253-je2-proteomics/blob/"
    "ff9a6546437974b2834cd182ba2c376403019bb4/"
    "Sarya_George/reports/step09_exclusive_phosphosite_proteins.md"
)
METHOD_URL = "https://weblogo.berkeley.edu/manual.html"


def text(draw, x, y, label, size=9, color=INK, font="Helvetica", anchor="left"):
    draw.setFont(font, size)
    draw.setFillColor(HexColor(color))
    positions = {"left": draw.drawString, "middle": draw.drawCentredString, "right": draw.drawRightString}
    positions[anchor](x, y, label)


def link(draw, x, y, label, url, size=8.1):
    text(draw, x, y, label, size, TEAL, "Helvetica-Bold")
    width = pdfmetrics.stringWidth(label, "Helvetica-Bold", size)
    draw.linkURL(url, (x, y - 2, x + width, y + size + 1), thickness=0)
    return width


def geometry(draw, x, y, width, height, color, radius=0):
    draw.setFillColor(HexColor(color))
    if radius:
        draw.roundRect(x, y, width, height, radius, fill=1, stroke=0)
    else:
        draw.rect(x, y, width, height, fill=1, stroke=0)


def line(draw, x0, x1, y, ink=RULE, weight=0.6):
    draw.setStrokeColor(HexColor(ink))
    draw.setLineWidth(weight)
    draw.line(x0, y, x1, y)


def protein_name(name):
    name = name.replace("HslU--HslV", "HslU/HslV")
    return name[0].upper() + name[1:]


def extract_windows(models):
    rows = []
    for protein in models:
        sequence = protein["sequence"]
        for site in protein["sites"]:
            index = site["position"] - 1
            sequence_window = "".join(
                sequence[index + offset] if 0 <= index + offset < len(sequence) else "-"
                for offset in OFFSETS
            )
            assert sequence_window[10] == site["residue"]
            assert all(letter == "-" or letter in COLORS for letter in sequence_window)
            rows.append({
                "name": protein_name(protein["name"]),
                "id": protein["accession"],
                "site": f'{site["residue"]}{site["position"]}',
                "groups": [tuple(group) for group in site["groups"]],
                "window": sequence_window,
            })

    assert len(rows) == len({(row["id"], row["site"]) for row in rows}) == 13
    assert [sum(group in row["groups"] for row in rows) for group in GROUPS] == list(GROUP_COUNTS)
    missing = {(row["id"], row["site"]): [i for i, aa in zip(OFFSETS, row["window"]) if aa == "-"]
               for row in rows if "-" in row["window"]}
    assert missing == {("Q2FZY9", "S5"): [-10, -9, -8, -7, -6, -5], ("Q2FXL6", "T157"): [10]}
    assert next(row["window"] for row in rows if row["id"] == "Q2FV52") == "SNGQVNAYNPSGASGLFQTMP"
    return rows


def letter(draw, aa, bottom, centre, width, height):
    # Courier Bold gives isoleucine visible top and bottom bars.
    font, cap = ("Courier-Bold", 562) if aa == "I" else ("Helvetica-Bold", 718)
    advance = pdfmetrics.stringWidth(aa, font, 1000)
    draw.saveState()
    draw.translate(centre - width / 2, bottom + 0.12)
    draw.scale(width / advance, max(0.1, height - 0.55) / cap)
    text(draw, 0, 0, aa, 1000, COLORS[aa], font)
    draw.restoreState()


def draw_logo(draw, x, y, width, height, label, rows, large=False):
    geometry(draw, x, y, width, height, "#FFFFFF", 7)
    text(draw, x + 12, y + height - 20, label, 13.5 if large else 12.5, INK, "Helvetica-Bold")
    text(draw, x + width - 12, y + height - 20, f"{len(rows)} sequence windows",
         8.9, MUTED, "Helvetica", "right")
    plot_left = x + 47
    plot_right = x + width - 13
    col = (plot_right - plot_left) / len(OFFSETS)
    bottom = y + 42
    top = y + height - 33
    plot_height = top - bottom

    geometry(draw, plot_left + 10 * col + 0.55, y + 9, col - 1.1, top - y - 9, COLUMN, 2)
    for amount, number in ((0, "0%"), (0.5, "50%"), (1, "100%")):
        yy = bottom + amount * plot_height
        line(draw, plot_left, plot_right, yy, RULE, 0.72 if not amount else 0.4)
        text(draw, plot_left - 8, yy - 2.2, number, 7, MUTED, "Helvetica", "right")
    text(draw, plot_left - 8, y + 28, "Pos.", 7.2, MUTED, "Helvetica-Bold", "right")
    text(draw, plot_left - 8, y + 13, "n", 7.2, MUTED, "Helvetica-Bold", "right")

    for j, offset in enumerate(OFFSETS):
        horizontal = plot_left + (j + 0.5) * col
        counts = Counter(row["window"][j] for row in rows if row["window"][j] != "-")
        total = sum(counts.values())
        cursor = bottom
        for aa, count in sorted(counts.items(), key=lambda pair: (pair[1], pair[0])):
            column_height = plot_height * count / total
            letter(draw, aa, cursor, horizontal, col * (0.77 if large else 0.81), column_height)
            cursor += column_height
        display = f"{offset:+d}" if offset else "0"
        text(draw, horizontal, y + 28, display, 7.55 if large else 6.8,
             FOCUS if not offset else MUTED, "Helvetica-Bold" if not offset else "Helvetica", "middle")
        text(draw, horizontal, y + 13, str(total), 7.25,
             GAP_COLOR if total < len(rows) else MUTED, "Helvetica-Bold" if total < len(rows) else "Helvetica", "middle")


def draw_header(draw, title, subtitle, page):
    geometry(draw, 0, 0, W, H, PAPER)
    text(draw, MARGIN, H - 39, title, 25, INK, "Helvetica-Bold")
    text(draw, MARGIN, H - 58, subtitle, 10.2, MUTED)
    text(draw, W - MARGIN, H - 38, "S. aureus JE2", 9.5, TEAL, "Helvetica-Bold", "right")
    line(draw, MARGIN, W - MARGIN, 32)
    x = MARGIN
    x += link(draw, x, 19, "JE2 phosphosites and experimental groups", ANALYSIS_URL) + 17
    link(draw, x, 19, "Frequency-based sequence logos", METHOD_URL)
    text(draw, W - MARGIN, 19, f"{page} / 2", 8, MUTED, "Helvetica-Bold", "right")


def page_logos(draw, rows):
    draw_header(draw, "Sequence logos around JE2 phosphosites",
                "Amino acid frequencies from -10 to +10 relative to the phosphosite · Position 0 is highlighted", 1)
    content_width = W - 2 * MARGIN
    draw_logo(draw, MARGIN, 378, content_width, 139, "Pooled · 13 distinct phosphosites", rows, True)
    gap = 12
    half_width = (content_width - gap) / 2
    for i, (medium, mutant) in enumerate(GROUPS):
        selection = [row for row in rows if (medium, mutant) in row["groups"]]
        draw_logo(draw, MARGIN + (i % 2) * (half_width + gap), 218 - (i // 2) * 155,
                  half_width, 149, f"{medium} · {mutant} vs. anc", selection)
    text(draw, MARGIN, 45,
         "Letter height shows the frequency among available residues. n is the number of sequences at each position. "
         "Each group contains 3 to 4 hits.", 8.2, MUTED)
    draw.showPage()


def page_windows(draw, rows):
    draw_header(draw, "Aligned sequence windows",
                "13 distinct JE2 phosphosites · 10 amino acids on either side of position 0", 2)
    table_x, table_y = MARGIN, 137
    width, top = W - 2 * MARGIN, 516
    row_height = 25
    header_y = top - 22
    seq_x = table_x + 329
    col = 10.65
    seq_width = len(OFFSETS) * col
    group_x = seq_x + seq_width + 17
    geometry(draw, table_x, table_y, width, top - table_y, "#FFFFFF", 7)
    text(draw, table_x + 12, header_y, "JE2 protein", 9.2, INK, "Helvetica-Bold")
    text(draw, table_x + 219, header_y, "UniProt", 9.2, INK, "Helvetica-Bold")
    text(draw, table_x + 287, header_y, "Site", 9.2, INK, "Helvetica-Bold")
    text(draw, group_x, header_y, "Experimental group", 9.2, INK, "Helvetica-Bold")

    first_y = header_y - 21
    for i in range(len(rows)):
        base_y = first_y - i * row_height
        if i % 2 == 0:
            geometry(draw, table_x + 1, base_y - 9, width - 2, row_height, "#F7FAFB")

    center_x = seq_x + 10.5 * col
    geometry(draw, center_x - col / 2, first_y - (len(rows) - 1) * row_height - 10,
             col, 21 + 10 + (len(rows) - 1) * row_height + 11, COLUMN, 2)
    for j, value in enumerate(OFFSETS):
        if value in (-10, -5, 0, 5, 10):
            display = f"{value:+d}" if value else "0"
            text(draw, seq_x + (j + 0.5) * col, header_y, display, 7.7,
                 FOCUS if not value else MUTED, "Helvetica-Bold" if not value else "Helvetica", "middle")
    for i, row in enumerate(rows):
        y = first_y - i * row_height
        text(draw, table_x + 12, y, row["name"], 9.3)
        text(draw, table_x + 219, y, row["id"], 9.1, MUTED, "Courier")
        text(draw, table_x + 287, y, row["site"], 9.3, FOCUS, "Helvetica-Bold")
        group_label = "  /  ".join(f"{medium} · {mutant}" for medium, mutant in row["groups"])
        text(draw, group_x, y, group_label, 9.3)
        for j, aa in enumerate(row["window"]):
            tone = FOCUS if j == 10 else (GAP_COLOR if aa == "-" else INK)
            text(draw, seq_x + (j + 0.5) * col, y, aa, 10.2, tone,
                 "Courier-Bold" if j == 10 else "Courier", "middle")

    text(draw, MARGIN, 108, "Gaps at protein termini", 10.2, INK, "Helvetica-Bold")
    text(draw, MARGIN, 91, "CsbD S5 has gaps at -10 to -5 and UspA T157 has a gap at +10. "
         "Frequencies use the residues available at each position.", 9.0, MUTED)
    text(draw, MARGIN, 71, "CsbD S5 contributes to both stp groups and appears once in the pooled logo. "
         "Windows use the sequences of the experimental strain JE2.", 9.0, MUTED)
    draw.showPage()


def render(models, output):
    """Render the frequency logos and aligned JE2 sequence windows to a PDF."""
    rows = extract_windows(models)
    output = Path(output)
    output.parent.mkdir(exist_ok=True, parents=True)
    draw = Canvas(str(output), pagesize=(W, H), pageCompression=1)
    draw.setTitle("JE2 phosphosite sequence logos (±10 amino acids)")
    draw.setSubject("Four experimental groups and a pooled logo for the observed STYH hits")
    page_logos(draw, rows)
    page_windows(draw, rows)
    draw.save()
