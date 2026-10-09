"""Independent check of the paper transcription used by Code B.

Extracts Tables 2, 3, 15-19 directly from the article PDF (PyMuPDF text layer),
parses every model row (8 loss values + 8 MCS p-values) and compares each cell
with docs/paper_values/code_b_forecast_tables.csv.
Usage (from project root):
  python3 -I docs/review_B/verify_paper_cells.py [path/to/article.pdf]

If no argument is supplied, PAPER_PDF is used when set, otherwise the historical
project-root filename is tried.
Writes docs/review_B/paper_cells_check.csv and docs/review_B/paper_text.txt.
"""
import csv, re, sys, os
try:
    import pymupdf as fitz
except ImportError:  # PyMuPDF versions before the pymupdf import name
    import fitz

ROOT = os.getcwd()
pdf = (sys.argv[1] if len(sys.argv) > 1 else
       os.environ.get("PAPER_PDF", os.path.join(ROOT, "04_Bibinger_Yu_Zhang_2026_JBES.pdf")))
if not os.path.isfile(pdf):
    raise SystemExit(f"Article PDF not found: {pdf}")
doc = fitz.open(pdf)
text = "".join(f"\n===== PAGE {i + 1}\n" + p.get_text() for i, p in enumerate(doc))
open(os.path.join(ROOT, "docs/review_B/paper_text.txt"), "w").write(text)

H = [1, 2, 3, 4, 5, 10, 15, 20]
TABLES = {"2", "3", "15", "16", "17", "18", "19"}
MODELS = ["fBm", "bfBm", "mfBm3", "mfBm4", "mfBm5", "HAR", "LHAR"] + \
    [p + str(k) for p in ("VHAR", "VHARF", "VLHAR") for k in range(2, 6)]
num = re.compile(r"^-?(\d*\.\d+|0)$")  # p-values are sometimes printed as a bare 0

def panel_of(table, model, seen):
    if table == "2":
        return ["full", "period1", "period2"][seen]
    if table in {"3", "15", "16"}:
        return "full"
    if model.startswith("VHARF"): return "vharf"
    if model.startswith("VLHAR") or model == "LHAR": return "vlhar"
    if model.startswith("VHAR"): return "vhar"
    if model == "HAR": return ["vhar", "vharf"][seen]
    return "mfbm"

cells, cur, seen, page = {}, None, {}, None
lines = [l.strip() for l in text.split("\n")]
i = 0
while i < len(lines):
    l = lines[i]
    if l.startswith("===== PAGE"):
        page = int(l.split()[-1])
    m = re.match(r"^Table (\d+):", l)
    if m:
        cur = m.group(1) if m.group(1) in TABLES else None
        seen = {}
    elif cur and l in MODELS:
        vals, j = [], i + 1
        while len(vals) < 16 and j < len(lines):
            if num.match(lines[j]): vals.append(float(lines[j]))
            elif lines[j] in MODELS: break
            j += 1
        if len(vals) == 16:
            k = seen.get(l, 0); seen[l] = k + 1
            pan = panel_of(cur, l, k)
            for h, v, p in zip(H, vals[:8], vals[8:]):
                cells[(f"Table {cur}", pan, l, h)] = (v, p, page)
            i = j - 1
    i += 1

ref = list(csv.DictReader(open(os.path.join(ROOT, "docs/paper_values/code_b_forecast_tables.csv"))))
out, bad = [], 0
for r in ref:
    key = (r["table"], r["panel"], r["model"], int(r["h"]))
    v, p, pg = cells.get(key, (None, None, None))
    ok_v = v is not None and abs(v - float(r["value"])) < 5e-5
    ok_p = p is not None and abs(p - float(r["mcs_p"])) < 5e-5
    bad += not (ok_v and ok_p)
    out.append(dict(table=key[0], panel=key[1], model=key[2], h=key[3], pdf_page=pg,
                    csv_value=r["value"], pdf_value=v, csv_mcs_p=r["mcs_p"], pdf_mcs_p=p,
                    value_match=ok_v, mcs_p_match=ok_p))
w = csv.DictWriter(open(os.path.join(ROOT, "docs/review_B/paper_cells_check.csv"), "w", newline=""),
                   fieldnames=list(out[0]))
w.writeheader(); w.writerows(out)
print(f"parsed {len(cells)} cells from PDF; reference rows {len(ref)}; mismatches {bad}")
sys.exit(1 if bad else 0)
