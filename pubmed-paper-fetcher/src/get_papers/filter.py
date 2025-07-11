from typing import List, Dict

def extract_non_academic_authors(papers: List[Dict], debug: bool = False) -> List[Dict]:
    results = []

    academic_keywords = ["university", "college", ".edu", "institute", "school", ".ac."]
    
    for paper in papers:
        non_academic_authors = []
        companies = []

        for author in paper.get("Authors", []):
            affil = author.get("affiliation", "").lower()
            if not any(keyword in affil for keyword in academic_keywords):
                non_academic_authors.append(author["name"])
                companies.append(author["affiliation"])

        if non_academic_authors:
            results.append({
                "PubmedID": paper["PubmedID"],
                "Title": paper["Title"],
                "Publication Date": paper["Publication Date"],
                "Non-academic Author(s)": ", ".join(non_academic_authors),
                "Company Affiliation(s)": ", ".join(companies),
                "Corresponding Author Email": paper["Corresponding Author Email"]
            })

            if debug:
                print(f"✅ Paper {paper['PubmedID']} has non-academic authors: {non_academic_authors}")

    return results
