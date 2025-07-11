# 🧪 PubMed Paper Fetcher (get-papers-list)

A Python command-line tool to fetch PubMed research papers based on a user-defined query. The program filters papers that have at least one author affiliated with a **pharmaceutical** or **biotech** company and outputs results to a CSV file or the console.

---

## ✅ Features

* 🔎 Supports full **PubMed query syntax**
* 🏢 Identifies **non-academic authors** (e.g., pharma, biotech companies)
* 📬 Extracts **corresponding author emails**
* 📄 Outputs to CSV or displays results in the terminal
* 💻 Built with **Poetry**, uses **typed Python**, and structured for extensibility

---

## 📂 Project Structure

```
pubmed-paper-fetcher/
│
├── src/
│   └── get_papers/
│       ├── __init__.py
│       ├── main.py         # CLI entry point (argparse logic)
│       ├── fetch.py        # Functions to fetch paper metadata from PubMed
│       └── filter.py       # Heuristics for identifying non-academic/company authors
│
├── pyproject.toml          # Poetry configuration
├── README.md               # Project documentation
```

---

## ⚙️ Installation

### 1. Clone the repository

```bash
git clone git@github.com:12shrubi/pubmed-paper-fetcher.git
cd pubmed-paper-fetcher
```

### 2. Install Poetry (if not installed)

```bash
curl -sSL https://install.python-poetry.org | python3 -
```

### 3. Install dependencies

```bash
poetry install
```

---

## 🚀 Usage

Run the CLI tool with Poetry:

```bash
poetry run get-papers-list "your pubmed search query"
```

### Options

| Flag          | Description                                            |
| ------------- | ------------------------------------------------------ |
| `-h, --help`  | Show usage instructions                                |
| `-d, --debug` | Print debug output (IDs, fetch status, etc.)           |
| `-f, --file`  | Save output to CSV file instead of printing to console |

### Examples

```bash
# Print to console
poetry run get-papers-list "cancer immunotherapy"

# Save results to CSV
poetry run get-papers-list "covid vaccine" -f results.csv

# Debug mode
poetry run get-papers-list "gene therapy" -d
```

---

## 🏢 How Non-Academic Authors Are Identified

We use basic heuristics to detect **company affiliations**:

* Emails not ending in `.edu` or `.ac.in`
* Affiliation strings containing:

  * Company names (e.g., Pfizer, Roche, Merck)
  * Keywords: `biotech`, `pharma`, `inc`, `corp`

These can be extended in `filter.py`.

---

## 🧲 Sample CSV Output Format

| PubmedID | Title | Publication Date | Non-academic Author(s) | Company Affiliation(s) | Corresponding Author Email |
| -------- | ----- | ---------------- | ---------------------- | ---------------------- | -------------------------- |

---

## 📡 Tools & Libraries Used

* [`Biopython`](https://biopython.org/) – to access PubMed API
* [`Poetry`](https://python-poetry.org/) – dependency management and packaging
* [`argparse`](https://docs.python.org/3/library/argparse.html) – for CLI interface
* LLM: [ChatGPT](https://openai.com/chatgpt) used to help structure, debug, and finalize functionality

---

## 🛆 Optional: Publish to TestPyPI (Bonus)

To prepare the package:

```bash
poetry config repositories.test-pypi https://test.pypi.org/legacy/
poetry build
poetry publish -r test-pypi
```

---

## 📌 Evaluation Criteria Checklist

* ✅ Functional: Fetch + filter papers correctly
* ✅ CLI: Works with flags (`-f`, `-d`, etc.)
* ✅ Typed Python throughout
* ✅ Good code organization: modular + documented
* ✅ Robust to errors: invalid query, empty results
* ✅ Hosted on GitHub

---

## 🧐 Future Improvements

* Smarter NLP-based author filtering
* Unit tests for fetch & filter logic
* Save JSON output alongside CSV

---

## 👤 Author

**Shrubi Kumari**
[GitHub: 12shrubi](https://github.com/12shrubi)
📧 [shrubikumari2001@gmail.com](mailto:shrubikumari2001@gmail.com)

---

## 🧪 License

MIT License – Free to use and modify.
