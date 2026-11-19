#!/usr/bin/env python3
"""Reproduce the repository's theory documents as LaTeX appendices.

The appendices of the technical report are *verbatim conversions* (pandoc) of
the Markdown sources in docs/extension, at the commit the report is pinned to.
Nothing in the text is rewritten: proofs, tags and correction notes appear as in
the source.  The script only adds typesetting aids:
  * result-level anchors (\\label) before each bold result label, so that the
    main text can link to a result's proof;
  * a per-appendix proof index (generated/idx_<name>.tex);
  * \\path{} for file names, so long paths can break;
  * scaling of single-line displays that are wider than the page.
Run from docs/report:  python3 -I tools/make_appendices.py
"""
import pathlib, re, subprocess

HERE = pathlib.Path(__file__).resolve().parent.parent      # docs/report
SRC = HERE.parent / "extension"
OUT = HERE / "generated"

# (output name, source file, anchor prefix, optional list of '## ' headings to keep)
JOBS = [
    ("app_phase1", "PHASE1_foundations.md", "pA", None),
    ("app_phase2a", "PHASE2a_rxc.md", "pB", None),
    ("app_phase2b", "PHASE2b_districts.md", "pC", None),
    ("app_phase2c", "PHASE2c_monotone.md", "pD", None),
    ("app_phase3", "PHASE3_inference.md", "pE", None),
    ("app_phase5a", "PHASE5a_sign_module.md", "pF", None),
    ("app_pkgdesign", "PACKAGE_DESIGN.md", "pG",
     ["## 2. Audit of the GPT conversation", "## 5. New results for the extension modules",
      "## 7. Extension to difference-in-differences", "## 8. Proof register"]),
]

READER = "markdown+tex_math_dollars+lists_without_preceding_blankline-implicit_figures"
KINDS = r"Lemma|Proposition|Theorem|Corollary|Conjecture|Remark|Definition"
RESULT = re.compile(r"^(\*\*)(" + KINDS + r")\s+([^\s*(]*(?:[0-9′]|-[A-Z])[^\s*(]*?(?:\([a-z0-9]+\))?)[.:]?(\s|\*\*)")


def slug(kind, num):
    """Anchor slug for a result label, shared with the main text (tools/link_proofs.py)."""
    num = num.replace("′", "p").replace("Δ", "Delta").rstrip(".:")
    return (kind.lower() + "-" + re.sub(r"[^A-Za-z0-9.]+", "-", num)).strip("-")


def select_sections(text, keep):
    if keep is None:
        return text
    lines = text.splitlines()
    out, on = [lines[0]], False
    for ln in lines[1:]:
        if ln.startswith("## "):
            on = any(ln.startswith(k) for k in keep)
        if on:
            out.append(ln)
    return "\n".join(out) + "\n"


def add_anchors(text, prefix):
    """Insert \\phantomsection\\label before each bold result label; return text and index."""
    index, seen, out = [], set(), []
    for ln in text.splitlines():
        m = RESULT.match(ln)
        if m:
            kind, num = m.group(2), m.group(3)
            lab = f"{prefix}-res-{slug(kind, num)}"
            if lab not in seen:
                seen.add(lab)
                index.append((f"{kind} {num.rstrip('.:')}", lab))
                ln = f"\\phantomsection\\label{{{lab}}}" + ln
        out.append(ln)
    return "\n".join(out) + "\n", index


PATHLIKE = re.compile(r"\\texttt\{((?:[A-Za-z0-9.\-/]|\\_)+(?:\.(?:md|py|R|txt|tex|json|csv|Rmd)|/))\}")


def postprocess(tex, prefix):
    tex = tex.replace("\u0303", "").replace("\u0304", "")
    # unique anchors across appendices (pandoc's --id-prefix does not reach \label)
    for cmd in ("\\hypertarget{", "\\label{", "\\hyperlink{"):
        tex = tex.replace(cmd + prefix + "-res-", "\x00")          # protect our anchors
        tex = tex.replace(cmd, cmd + prefix + "-")
        tex = tex.replace("\x00", cmd + prefix + "-res-")
    # file names: \path allows line breaks at / . _
    # (not in headings, which are moving arguments)
    tex = "\n".join(ln if ("section{" in ln or "paragraph{" in ln) else
                    PATHLIKE.sub(lambda m: "\\path{" + m.group(1).replace("\\_", "_") + "}", ln)
                    for ln in tex.split("\n"))
    # narrow table cells: allow the first word of a cell to hyphenate
    tex = tex.replace("\\raggedright\\arraybackslash}p{", "\\raggedright\\arraybackslash\\hspace{0pt}}p{")
    tex = tex.replace("\\begin{minipage}[b]{\\linewidth}\\raggedright\n", "\\begin{minipage}[b]{\\linewidth}\\raggedright\\hspace{0pt}")
    # single-line displays wider than the page: scale to the line width
    def fit(m):
        body = m.group(1)
        if len(body) > 170 and ";\\qquad" in body and "\\\\" not in body:
            return "\\[\\begin{gathered}" + body.replace(";\\qquad", ";\\\\") + "\\end{gathered}\\]"
        if len(body) > 170 and "\\tag" not in body and "\\\\" not in body and "&" not in body:
            return "\\[\\resizebox{\\linewidth}{!}{$\\displaystyle " + body.strip() + "$}\\]"
        return m.group(0)
    tex = re.sub(r"\\\[(.*?)\\\]", fit, tex, flags=re.S)
    return tex


def convert(name, src, prefix, keep):
    text = (SRC / src).read_text(encoding="utf-8")
    text = select_sections(text, keep)
    text = re.sub(r"\A# .*\n", "", text)            # the appendix supplies its own heading
    text, index = add_anchors(text, prefix)
    res = subprocess.run(["pandoc", "-f", READER, "-t", "latex", "--wrap=preserve", "--columns=50"],
                         input=text, capture_output=True, text=True, check=True)
    tex = postprocess(res.stdout, prefix)
    header = (f"% GENERATED by tools/make_appendices.py from docs/extension/{src}\n"
              f"% Do not edit by hand; edit the source and regenerate.\n")
    (OUT / f"{name}.tex").write_text(header + tex, encoding="utf-8")
    # proof index
    rows = "\n".join(f"{lab_txt} & \\hyperref[{lab}]{{p.~\\pageref{{{lab}}}}} \\\\"
                     for lab_txt, lab in index)
    idx = ("% GENERATED proof index\n\\begin{center}\\footnotesize\n"
           "\\begin{longtable}{@{}L{0.62\\linewidth}l@{}}\n"
           "\\toprule Result (source label) & Page \\\\ \\midrule\\endfirsthead\n"
           "\\toprule Result (source label) & Page \\\\ \\midrule\\endhead\n"
           f"{rows}\n\\bottomrule\n\\end{{longtable}}\n\\end{{center}}\n")
    (OUT / f"idx_{name}.tex").write_text(idx if index else "% no labelled results\n", encoding="utf-8")
    (OUT / f"labels_{name}.txt").write_text("\n".join(f"{t}\t{l}" for t, l in index) + "\n", encoding="utf-8")
    print(f"{name}: {len(tex.splitlines())} lines, {len(index)} anchored results")


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True)
    for job in JOBS:
        convert(*job)
