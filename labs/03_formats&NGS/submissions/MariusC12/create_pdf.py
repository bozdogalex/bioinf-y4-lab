from pathlib import Path
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.styles import ParagraphStyle
from xml.sax.saxutils import escape

folder = Path(__file__).parent
input_file = folder / "notes.txt"
output_file = folder / "notes.pdf"

styles = getSampleStyleSheet()

title_style = ParagraphStyle(
    "TitleCustom",
    parent=styles["Title"],
    alignment=TA_CENTER,
    fontSize=16
)

story = []

with open(input_file, "r", encoding="utf-8") as file:
    for line in file:
        line = line.strip()

        if not line:
            story.append(Spacer(1, 8))
        elif line.startswith("# "):
            story.append(Paragraph(escape(line[2:]), title_style))
        elif line.startswith("## "):
            story.append(Paragraph(escape(line[3:]), styles["Heading2"]))
        else:
            story.append(Paragraph(escape(line), styles["BodyText"]))

        story.append(Spacer(1, 5))

document = SimpleDocTemplate(
    str(output_file),
    pagesize=A4,
    rightMargin=50,
    leftMargin=50,
    topMargin=45,
    bottomMargin=45
)

document.build(story)

print(f"PDF creat: {output_file}")