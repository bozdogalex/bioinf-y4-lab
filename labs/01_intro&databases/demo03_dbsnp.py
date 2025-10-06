from Bio import Entrez

Entrez.email = "student@example.com"

res = Entrez.esearch(
    db="snp",
    term="BRCA1[gene] AND Homo sapiens[organism]",
    retmax=5,
)
ids = Entrez.read(res)["IdList"]
print(f"Am găsit {len(ids)} SNP IDs.")

if ids:
    summ = Entrez.esummary(db="snp", id=",".join(ids), retmode="xml")
    records = Entrez.read(summ)
    docsums = records["DocumentSummarySet"]["DocumentSummary"]
    for d in docsums:
        snp_id = d.get("SNP_ID")
        doc = d.get("DOCSUM")
        print("SNP_ID:", snp_id, "| DOCSUM:", doc)
