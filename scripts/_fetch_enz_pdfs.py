import subprocess
from pathlib import Path

downloads = {
    "shenai.pdf": "https://www.jbc.org/article/S0021-9258(19)61288-3/pdf",
    "sijwali.pdf": "https://portlandpress.com/biochemj/article-pdf/360/2/481/710218/bj3600481.pdf",
    "wyatt.pdf": "https://febs.onlinelibrary.wiley.com/doi/pdfdirect/10.1016/S0014-5793%2802%2902241-X",
    "hill.pdf": "https://febs.onlinelibrary.wiley.com/doi/pdfdirect/10.1016/0014-5793%2894%2900940-6",
    "li.pdf": "https://www.jbc.org/article/S0021-9258(20)74591-8/pdf",
    "luker.pdf": "https://www.sciencedirect.com/science/article/pii/0166685196026515/pdfft?md5=&pid=1-s2.0-0166685196026515-main.pdf",
}

outdir = Path("/tmp/enz_pdfs")
outdir.mkdir(exist_ok=True)
for name, url in downloads.items():
    path = outdir / name
    cmd = [
        "curl", "-sL", "-A", "Mozilla/5.0",
        "-o", str(path), "-w", f"{name} %{{http_code}} %{{size_download}}\\n",
        url,
    ]
    print(subprocess.check_output(cmd, text=True).strip())

# extract text where possible
for pdf in outdir.glob("*.pdf"):
    if pdf.stat().st_size < 1000:
        print(pdf.name, "too small")
        continue
    try:
        txt = subprocess.check_output(["pdftotext", "-layout", str(pdf), "-"], text=True, errors="ignore")
        out = outdir / (pdf.stem + ".txt")
        out.write_text(txt, encoding="utf-8")
        print(pdf.name, "->", len(txt), "chars")
    except Exception as e:
        print(pdf.name, "pdftotext fail", e)
