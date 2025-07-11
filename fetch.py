from typing import List, Dict
from Bio import Entrez, Medline

Entrez.email = "shrubikumari2001@gmail.com"

def fetch_pubmed_ids(query: str, max_results: int = 20) -> List[str]:
    handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
    record = Entrez.read(handle)
    return record["IdList"]

def fetch_paper_details(pubmed_ids: List[str]) -> List[Dict]:
    handle = Entrez.efetch(db="pubmed", id=",".join(pubmed_ids), rettype="medline", retmode="text")
    records = Medline.parse(handle)
    results = []

    for record in records:
        result = {
            "PubmedID": record.get("PMID", ""),
            "Title": record.get("TI", ""),
            "Publication Date": record.get("DP", ""),
            "Authors": record.get("FAU", []),
            "Affiliations": record.get("AD", []),
            "Corresponding Author Email": extract_email(record.get("AD", ""))
        }
        results.append(result)

    return results

def extract_email(affiliation: str) -> str:
    import re
    match = re.search(r"[a-zA-Z0-9._%+-]+@[a-zA-Z0-9.-]+\.[a-zA-Z]{2,}", affiliation)
    return match.group(0) if match else ""

