import argparse
from get_papers.fetch import fetch_pubmed_data
from get_papers.filter import extract_non_academic_authors
from get_papers.utils import save_to_csv

def main():
    parser = argparse.ArgumentParser(description="PubMed paper fetcher CLI")
    parser.add_argument("query", help="PubMed search query")
    parser.add_argument("-f", "--file", help="Output CSV filename")
    parser.add_argument("-d", "--debug", action="store_true", help="Enable debug output")
    
    args = parser.parse_args()

    print("🚀 PubMed paper fetcher CLI is working!")
    print(f"🔍 Searching PubMed for: {args.query}")

    # Step 1: Fetch papers
    papers = fetch_pubmed_data(args.query, debug=args.debug)

    # Step 2: Filter results
    filtered_results = extract_non_academic_authors(papers, debug=args.debug)

    # Step 3: Output
    if args.file:
        save_to_csv(filtered_results, args.file)
        print(f"✅ Results saved to {args.file}")
    else:
        for row in filtered_results:
            print(row)
