"""Optional: rebuild the five-page Chinese report and four readable panels.

Install requirements.txt and requirements-report.txt; run after run_all.py.
Supports only the Markdown subset used in docs/REPORT.md, not arbitrary Markdown.
"""
from pathlib import Path
import argparse
import csv
import html
import json
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.cidfonts import UnicodeCIDFont
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, Table, TableStyle, PageBreak
from reportlab.lib.utils import ImageReader

ROOT = Path(__file__).resolve().parents[1]


def panels():
    figs, docs = ROOT/"figures", ROOT/"docs"
    data = np.genfromtxt(figs/"cv.csv", delimiter=",", names=True)
    metrics = json.loads((figs/"validation.json").read_text())
    cox = metrics["device"]["Cox_F_m2"]
    plt.rcParams.update({"font.size": 10, "axes.titlesize": 11, "legend.fontsize": 8})
    def finish(fig, ax, name):
        ax.grid(alpha=.25)
        fig.savefig(docs/name, dpi=180, bbox_inches="tight")
        plt.close(fig)
    fig, ax = plt.subplots(figsize=(7, 2.8), constrained_layout=True)
    ax.plot(data["Vg_V"], data["Cqs_F_m2"]/cox, label="Quasistatic Poisson")
    ax.plot(data["Vg_V"], data["Chf_frozen_F_m2"]/cox, "--", label="HF: frozen electron AC response")
    ax.set(xlabel="Gate voltage (V)", ylabel="C / Cox")
    ax.legend()
    finish(fig, ax, "report_cv.png")
    fig, ax = plt.subplots(figsize=(7, 2.6), constrained_layout=True)
    traces = json.loads((figs/"newton_traces.json").read_text())
    for name, rows in traces.items():
        ax.semilogy([r["iteration"] for r in rows], [max(r["residual_inf"], 1e-16) for r in rows], "o-", ms=3, label=name)
    ax.set(xlabel="Damped Newton iteration", ylabel="Scaled residual")
    ax.legend()
    finish(fig, ax, "report_newton.png")
    mesh = np.genfromtxt(figs/"mesh_convergence.csv", delimiter=",", names=True)
    fig, ax = plt.subplots(figsize=(7, 2.7), constrained_layout=True)
    ax.loglog(mesh["nodes"]-1, mesh["max_capacitance_error_over_Cox"], "o-", label="Max capacitance error / Cox")
    ax.loglog(mesh["nodes"]-1, mesh["max_capacitance_error_over_Cox"][0]*((mesh["nodes"][0]-1)/(mesh["nodes"]-1))**2, "--", label="Second-order guide")
    ax.set(xlabel="Mesh intervals", ylabel="Error / Cox")
    ax.legend()
    finish(fig, ax, "report_mesh.png")
    data = np.genfromtxt(figs/"data/synthetic_cv.csv", delimiter=",", names=True)
    fig, ax = plt.subplots(figsize=(7, 2.8), constrained_layout=True)
    ax.plot(data["Vg_V"], data["clean_reference_C_F_m2"]*100, label="Synthetic first-integral truth")
    ax.plot(data["Vg_V"], data["synthetic_noisy_C_F_m2"]*100, ".", label="Synthetic noisy data")
    ax.plot(data["Vg_V"], data["fitted_Poisson_C_F_m2"]*100, "--", label="Fitted Poisson model")
    ax.set(xlabel="Gate voltage (V)", ylabel="QS C (uF/cm^2)")
    ax.legend()
    finish(fig, ax, "report_fit.png")


def build(font_path: Path | None = None):
    docs = ROOT/"docs"
    panels()
    candidate = Path("/System/Library/Fonts/Supplemental/Arial Unicode.ttf")
    if font_path is None and candidate.exists():
        font_path = candidate
    if font_path:
        pdfmetrics.registerFont(TTFont("ReportCJK", str(font_path)))
        font = "ReportCJK"
    else:
        pdfmetrics.registerFont(UnicodeCIDFont("STSong-Light"))
        font = "STSong-Light"
    body = ParagraphStyle("Body", fontName=font, fontSize=10.2, leading=16.2,
                          textColor=colors.HexColor("#253343"), wordWrap="CJK", spaceAfter=7)
    h1 = ParagraphStyle("H1", parent=body, fontSize=19, leading=25,
                        textColor=colors.HexColor("#12556a"), spaceAfter=13)
    h2 = ParagraphStyle("H2", parent=body, fontSize=12.5, leading=18,
                        textColor=colors.HexColor("#12556a"), spaceBefore=8, spaceAfter=7)
    small = ParagraphStyle("Small", parent=body, fontSize=9, leading=13.5, spaceAfter=2)
    caption = ParagraphStyle("Caption", parent=small, textColor=colors.HexColor("#526678"), spaceAfter=9)
    width = A4[0]-88
    story, paragraph = [], []
    def ptext(text):
        return html.escape(text.replace("`", "").replace("**", ""))
    def flush():
        if paragraph:
            story.append(Paragraph(ptext(" ".join(paragraph)), body))
            paragraph.clear()
    lines = (docs/"REPORT.md").read_text().splitlines()
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if not line:
            flush()
        elif line == "<!-- page -->":
            flush()
            story.append(PageBreak())
        elif line.startswith("# "):
            flush()
            story.append(Paragraph(ptext(line[2:]), h1))
        elif line.startswith("## "):
            flush()
            story.append(Paragraph(ptext(line[3:]), h2))
        elif line.startswith("!["):
            flush()
            match = re.fullmatch(r"!\[(.*?)\]\((.*?)\)", line)
            if not match:
                raise ValueError(f"Unsupported image line: {line}")
            file = (docs/match[2]).resolve()
            iw, ih = ImageReader(str(file)).getSize()
            scale = min(width/iw, 215/ih)
            story.append(Image(str(file), width=iw*scale, height=ih*scale))
            story.append(Paragraph(ptext(match[1]), caption))
        elif line.startswith("|"):
            flush()
            rows = []
            while i < len(lines) and lines[i].strip().startswith("|"):
                cells = [s.strip() for s in lines[i].strip().strip("|").split("|")]
                if not all(re.fullmatch(r":?-+:?", cell) for cell in cells):
                    rows.append([Paragraph(ptext(cell), small) for cell in cells])
                i += 1
            i -= 1
            count = len(rows[0])
            fractions = {2: [.48, .52], 3: [.36, .32, .32], 4: [.26, .25, .27, .22]}.get(count, [1/count]*count)
            table = Table(rows, colWidths=[width*f for f in fractions], hAlign="LEFT")
            table.setStyle(TableStyle([
                ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#e7f1f4")),
                ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#f5f7f9")]),
                ("LINEBELOW", (0, 0), (-1, 0), .6, colors.HexColor("#b6cbd3")),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("LEFTPADDING", (0, 0), (-1, -1), 8),
                ("RIGHTPADDING", (0, 0), (-1, -1), 8),
                ("TOPPADDING", (0, 0), (-1, -1), 5),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 5)]))
            story.extend((table, Spacer(1, 10)))
        else:
            paragraph.append(line)
        i += 1
    flush()
    def decoration(canvas, document):
        canvas.setStrokeColor(colors.HexColor("#b6cbd3"))
        canvas.line(44, A4[1]-40, A4[0]-44, A4[1]-40)
        canvas.setFont(font, 8)
        canvas.setFillColor(colors.HexColor("#526678"))
        canvas.drawString(44, A4[1]-30, "mini-tcad-moscap | Technical report | Baseline v1")
        canvas.drawString(44, 25, "2026-10-19 | Classical MOS electrostatics and synthetic inverse validation")
        canvas.drawRightString(A4[0]-44, 25, f"{document.page} / 5")
    output = docs/"REPORT.pdf"
    document = SimpleDocTemplate(str(output), pagesize=A4, leftMargin=44, rightMargin=44,
                                 topMargin=61, bottomMargin=43, title="MOS capacitor self-consistent simulation",
                                 author="mini-tcad-moscap educational project")
    document.build(story, onFirstPage=decoration, onLaterPages=decoration)
    print(f"Built {output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--font", type=Path, help="Optional CJK-capable TrueType font")
    build(parser.parse_args().font)
