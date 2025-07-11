from typing import List, Dict

# A sample list of known pharma/biotech indicators
COMPANY_KEYWORDS = [
    "pharma", "biotech", "Pfizer", "Roche", "Novartis", "Merck", "AstraZeneca", "Sanofi", "GSK", "Bayer"
]

def is_company_affiliation(affiliation: str) -> bool:
    affiliation_lower = affiliation.lower()
    return any(keyword.lower() in affiliation_lower for keyword in COMPANY_KEYWORDS)

def filter_papers(papers: List[Dict]) -> List[Dict]:
    filtered = []

    for paper in papers:
        affiliations = paper.get("Affiliations", [])
        if isinstance(affiliations, str):
            affiliations = [affiliations]

        company_affils = [a for a in affiliations if is_company_affiliation(a)]
        if company_affils:
            paper["Company Affiliation(s)"] = company_affils
            paper["Non-academic Author(s)"] = paper.get("Authors", [])
            filtered.append(paper)

    return filtered

