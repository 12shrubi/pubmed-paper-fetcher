from typing import List, Dict
from Bio import Entrez

Entrez.email = "shrubikumari2001@gmail.com"

def fetch_pubmed_data(query: str, debug: bool = False) -> List[Dict]:
    if debug:
        print(f"🔍 Fetching papers for query: {query}")

    handle = Entrez.esearch(db="pubmed", term=query, retmax=10)
    record = Entrez.read(handle)
    ids = record["IdList"]
    handle.close()

    if debug:
        print(f"📄 Found {len(ids)} paper IDs: {ids}")

    papers = []
    if ids:
        fetch_handle = Entrez.efetch(db="pubmed", id=",".join(ids), rettype="medline", retmode="text")
        data = fetch_handle.read()
        fetch_handle.close()
        for id_ in ids:
            papers.append({
                "PubmedID": id_,
                "Title": f"Title for {id_}",
                "Publication Date": "2023-01-01",
                "Authors": [
                    {"name": "John Doe", "affiliation": "ABC Biotech Inc.", "email": "john@abcbiotech.com"},
                    {"name": "Dr. Smith", "affiliation": "Harvard University", "email": "smith@harvard.edu"}
                ],
                "Corresponding Author Email": "john@abcbiotech.com"
            })
    return papers
