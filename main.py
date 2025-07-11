import argparse
import csv
from get_papers.fetch import fetch_pubmed_ids, fetch_paper_details
from get_papers.filter import filter_papers

def save_to_csv(papers, filename):
    with open(filename, mode="w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "PubmedID", "Title", "Publication Date", "Non-academic Author(s)",
            "Company Affiliation(s)", "Corresponding Author Email"
        ])
        writer.writeheader()
        for paper in papers:
            writer.writerow({
                "PubmedID": paper.get("PubmedID", ""),
                "Title": paper.get("Title", ""),
                "Publication Date": paper.get("Publication Date", ""),
                "Non-academic Author(s)": ", ".join(paper.get("Non-academic Author(s)", [])),
                "Company Affiliation(s)": ", ".join(paper.get("Company Affiliation(s)", [])),
                "Corresponding Author Email": paper.get("Corresponding Author Email", ""),
            })

def main():
    parser = argparse.ArgumentParser(description="Fetch PubMed papers with pharma/biotech authors.")
    parser.add_argument("query", help="Search query")
    parser.add_argument("-f", "--file", help="Filename to save results (CSV)")
    parser.add_argument("-d", "--debug", action="store_true", help="Enable debug output")
    args = parser.parse_args()

    if args.debug:
        print(f"🔍 Searching PubMed for: {args.query}")

    ids = fetch_pubmed_ids(args.query)
    papers = fetch_paper_details(ids)
    filtered = filter_papers(papers)

    if args.file:
        save_to_csv(filtered, args.file)
        print(f"✅ Results saved to {args.file}")
    else:
        for paper in filtered:
            print("\n---")
            print(f"📄 Title: {paper['Title']}")
            print(f"🏢 Company: {', '.join(paper['Company Affiliation(s)'])}")
            print(f"📧 Email: {paper['Corresponding Author Email']}")

