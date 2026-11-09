
import urllib.request
import urllib.parse
import xml.etree.ElementTree as ET
from pathlib import Path

handle = "MariusC12"
query = "TP53 AND cancer"

output_file = Path(
    f"labs/03_formats&NGS/submissions/{handle}/pubmed_{handle}.txt"
)

base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

search_params = urllib.parse.urlencode({
    "db": "pubmed",
    "term": query,
    "retmax": 5,
    "retmode": "xml"
})

search_url = base_url + "esearch.fcgi?" + search_params

request = urllib.request.Request(
    search_url,
    headers={"User-Agent": "MariusC12-Bioinformatics-Student"}
)

with urllib.request.urlopen(request) as response:
    search_xml = ET.fromstring(response.read())

ids = [
    element.text
    for element in search_xml.findall(".//Id")
]

if not ids:
    raise RuntimeError("Nu au fost gasite articole.")

fetch_params = urllib.parse.urlencode({
    "db": "pubmed",
    "id": ",".join(ids),
    "rettype": "abstract",
    "retmode": "xml"
})

fetch_url = base_url + "efetch.fcgi?" + fetch_params

request = urllib.request.Request(
    fetch_url,
    headers={"User-Agent": "MariusC12-Bioinformatics-Student"}
)

with urllib.request.urlopen(request) as response:
    articles_xml = ET.fromstring(response.read())

with open(output_file, "w", encoding="utf-8") as out:
    for index, article in enumerate(
        articles_xml.findall(".//PubmedArticle"), 1
    ):
        title = article.findtext(
            ".//ArticleTitle", default="No title"
        )

        authors = []
        for author in article.findall(".//Author"):
            last = author.findtext("LastName", default="")
            fore = author.findtext("ForeName", default="")
            collective = author.findtext("CollectiveName")

            if collective:
                authors.append(collective)
            else:
                name = f"{fore} {last}".strip()
                if name:
                    authors.append(name)

        abstract_parts = []
        for part in article.findall(".//Abstract/AbstractText"):
            label = part.get("Label", "")
            text = "".join(part.itertext())
            if label:
                abstract_parts.append(f"{label}: {text}")
            else:
                abstract_parts.append(text)

        abstract = "\n".join(abstract_parts)

        out.write(f"ARTICLE {index}\n")
        out.write(f"Title: {title}\n")
        out.write(f"Authors: {', '.join(authors) or 'Not available'}\n")
        out.write(f"Abstract:\n{abstract or 'Not available'}\n")
        out.write("-" * 70 + "\n\n")

print(f"[OK] PubMed results -> {output_file.resolve()}")
