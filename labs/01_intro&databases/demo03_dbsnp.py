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
    handle = Entrez.esummary(db="snp", id=",".join(ids), retmode="xml")
    docsums = Entrez.read(handle)
    handle.close()

    for d in docsums['DocumentSummarySet']['DocumentSummary']:
        snp_id = d.attributes.get('uid')
        chrpos = d.get('CHRPOS', 'N/A')
        fxn = d.get('FXN_CLASS', 'N/A')
        print(f"SNP_ID: {snp_id} | CHRPOS: {chrpos} | Funcție: {fxn}")
