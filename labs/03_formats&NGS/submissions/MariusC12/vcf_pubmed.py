
import gzip
import urllib.request
import urllib.parse
import xml.etree.ElementTree as ET
from pathlib import Path

handle = "MariusC12"

vcf_file = Path(f"data/work/{handle}/lab03/clinvar.vcf.gz")
output_file = Path(
    f"labs/03_formats&NGS/submissions/{handle}/variants_{handle}.txt"
)

base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

variants = []

with gzip.open(vcf_file, "rt") as vcf:
    for line in vcf:
        if line.startswith("#"):
            continue

        columns = line.strip().split("\t")

        if len(columns) < 8:
            continue

        chrom, pos, variant_id, ref, alt = columns[:5]
        info = columns[7]

        if "GENEINFO=TP53:" not in info:
            continue

        rs_id = None

        for item in info.split(";"):
            if item.startswith("RS="):
                rs_id = item.split("=", 1)[1]
                break

        if rs_id:
            variants.append({
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "rs_id": rs_id
            })

        if len(variants) == 2:
            break

if len(variants) < 2:
    raise RuntimeError("Nu am gasit doua variante TP53 cu ID RS.")

with open(output_file, "w", encoding="utf-8") as out:
    for index, variant in enumerate(variants, 1):
        query = f"rs{variant['rs_id']} AND TP53"

        params = urllib.parse.urlencode({
            "db": "pubmed",
            "term": query,
            "retmax": 5,
            "retmode": "xml"
        })

        url = base_url + "esearch.fcgi?" + params

        request = urllib.request.Request(
            url,
            headers={"User-Agent": "MariusC12-Bioinformatics-Student"}
        )

        with urllib.request.urlopen(request) as response:
            search_xml = ET.fromstring(response.read())

        ids = [
            element.text
            for element in search_xml.findall(".//Id")
        ]

        out.write(f"VARIANTA {index}\n")
        out.write(f"Cromozom: {variant['chrom']}\n")
        out.write(f"Pozitie: {variant['pos']}\n")
        out.write(f"REF: {variant['ref']}\n")
        out.write(f"ALT: {variant['alt']}\n")
        out.write(f"ID: rs{variant['rs_id']}\n")

        if ids:
            fetch_params = urllib.parse.urlencode({
                "db": "pubmed",
                "id": ",".join(ids),
                "retmode": "xml"
            })

            fetch_url = base_url + "efetch.fcgi?" + fetch_params

            request = urllib.request.Request(
                fetch_url,
                headers={"User-Agent": "MariusC12-Bioinformatics-Student"}
            )

            with urllib.request.urlopen(request) as response:
                articles_xml = ET.fromstring(response.read())

            out.write("Articole PubMed:\n")

            for article in articles_xml.findall(".//PubmedArticle"):
                title = article.findtext(
                    ".//ArticleTitle",
                    default="Titlu indisponibil"
                )
                out.write(f"- {title}\n")
        else:
            out.write("Nu au fost gasite articole PubMed.\n")

        out.write("-" * 60 + "\n")

print(f"[OK] Rezultate salvate in: {output_file}")
