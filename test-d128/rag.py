#!/usr/bin/env python3
"""Indice/RAG minimo del repo (sin embeddings, sin deps).

Uso:
  python test-d128/rag.py                          # reconstruir indice
  python test-d128/rag.py "donde se define BG_0 del plus"
  python test-d128/rag.py -k 10 "parche serial cero en verify"
  python test-d128/rag.py --def BG_0
  python test-d128/rag.py --fn verify
  python test-d128/rag.py --file auditor/verify.c "hash"

El indice guarda cada linea de codigo con archivo:linea en
test-d128/.indice.json. La busqueda puntua por tokens (coincidencia de
palabras), no solo regex, para poder preguntar en lenguaje natural.
"""
import argparse
import json
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
OUT = pathlib.Path(__file__).resolve().parent / ".indice.json"
SKIP_DIRS = {".git", "artifacts", "__pycache__", "auditor_d128",
             "auditor_d128_audit", ".indice.json"}
EXTS = (".c", ".h", ".py", ".yml", ".md", ".sh")

DEF_RE = re.compile(r"^\s*#\s*define\s+([A-Za-z_]\w*)\s*(.*)$")
FN_RE = re.compile(
    r"^\s*(?:static\s+)?(?:inline\s+)?(?:const\s+)?[A-Za-z_][\w \*]*?"
    r"\b([A-Za-z_]\w*)\s*\([^;{]*\)\s*\{?\s*$")
STRUCT_RE = re.compile(r"^\s*(?:typedef\s+)?struct\s+([A-Za-z_]\w*)?\s*\{")


def scan():
    idx = {"defines": {}, "functions": {}, "structs": {}, "lines": [],
           "files": []}
    for p in ROOT.rglob("*"):
        if p.suffix not in EXTS or any(s in p.parts for s in SKIP_DIRS):
            continue
        rel = p.relative_to(ROOT).as_posix()
        idx["files"].append(rel)
        try:
            lines = p.read_text(encoding="utf-8", errors="replace").splitlines()
        except OSError:
            continue
        for n, raw in enumerate(lines, 1):
            s = raw.strip()
            idx["lines"].append([rel, n, s])
            m = DEF_RE.match(raw)
            if m:
                idx["defines"].setdefault(m.group(1), []).append(
                    {"v": m.group(2).strip(), "f": rel, "l": n})
            m = STRUCT_RE.match(raw)
            if m and m.group(1):
                idx["structs"].setdefault(m.group(1), []).append(
                    {"f": rel, "l": n})
            m = FN_RE.match(raw)
            if m and m.group(1) not in ("if", "for", "while", "switch",
                                        "return", "sizeof"):
                idx["functions"].setdefault(m.group(1), []).append(
                    {"s": s, "f": rel, "l": n})
    OUT.write_text(json.dumps(idx, ensure_ascii=False), encoding="utf-8")
    return idx


def load(rebuild=False):
    if rebuild or not OUT.exists():
        return scan()
    return json.loads(OUT.read_text(encoding="utf-8"))


STOP = {"donde", "se", "define", "del", "de", "la", "el", "los", "las",
        "en", "un", "una", "por", "con", "para", "al", "y", "o", "que",
        "como", "esta", "este", "eso", "si", "no", "es", "son", "the",
        "where", "is", "are", "of", "to", "in", "for", "and", "or"}


def tokens(q):
    return [t.lower() for t in re.findall(r"[A-Za-z_]\w*|\d+", q)
            if t.lower() not in STOP]


def search(idx, query, k, only_file=None):
    toks = tokens(query)
    scored = []
    for f, n, s in idx["lines"]:
        if only_file and only_file not in f:
            continue
        if not s or len(s) > 240:
            continue
        low = s.lower()
        score = 0
        for t in toks:
            if t in low:
                score += 3 if re.search(r"\b" + re.escape(t) + r"\b", low) else 1
        if score:
            scored.append((score, f, n, s))
    scored.sort(key=lambda x: (-x[0], x[1], x[2]))
    for score, f, n, s in scored[:k]:
        print(f"{f}:{n}  [s={score}]  {s}")
    if not scored:
        print("(sin resultados)")
    return len(scored)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("query", nargs="?")
    ap.add_argument("-k", type=int, default=20, help="max resultados")
    ap.add_argument("--def", dest="define", help="valor de un #define")
    ap.add_argument("--fn", dest="fn", help="firma de una funcion")
    ap.add_argument("--struct", dest="struct", help="donde esta un struct")
    ap.add_argument("--file", help="limitar busqueda a archivos que contengan")
    ap.add_argument("--rebuild", action="store_true")
    a = ap.parse_args()

    idx = load(a.rebuild)
    if not (a.query or a.define or a.fn or a.struct):
        print(f"indice: {len(idx['files'])} archivos, {len(idx['lines'])} lineas, "
              f"{len(idx['defines'])} defines, {len(idx['functions'])} funciones")
        return 0
    if a.define:
        for name in a.define.split(","):
            for e in idx["defines"].get(name, []):
                print(f"{name} = {e['v']}   ({e['f']}:{e['l']})")
        return 0
    if a.fn:
        for name in a.fn.split(","):
            for e in idx["functions"].get(name, []):
                print(f"{e['s']}   ({e['f']}:{e['l']})")
        return 0
    if a.struct:
        for name in a.struct.split(","):
            for e in idx["structs"].get(name, []):
                print(f"struct {name}   ({e['f']}:{e['l']})")
        return 0
    search(idx, a.query, a.k, a.file)
    return 0


if __name__ == "__main__":
    sys.exit(main())
