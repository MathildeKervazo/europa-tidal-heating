# -*- coding: utf-8 -*-
"""Extrait la geometrie vectorielle et le texte de figure1.pdf (manuscrit) en JSON (necessite PyMuPDF).
Usage : python3 extract_geometry.py published/figure1.pdf figure1_geometry     (ecrit figure1_geometry.json)"""
import sys, json
try:
    import pymupdf as fitz
except ImportError:
    import fitz

src, out = sys.argv[1], sys.argv[2]
doc = fitz.open(src)
pg = doc[0]
W, H = pg.rect.width, pg.rect.height

def pt(p):
    return [round(float(p.x), 4), round(float(p.y), 4)]

def conv_items(items):
    ops = []
    for it in items:
        k = it[0]
        if k == "l":
            ops.append(["l", pt(it[1]), pt(it[2])])
        elif k == "c":
            ops.append(["c", pt(it[1]), pt(it[2]), pt(it[3]), pt(it[4])])
        elif k == "re":
            r = it[1]
            ops.append(["re", [round(r.x0, 4), round(r.y0, 4), round(r.x1, 4), round(r.y1, 4)]])
        elif k == "qu":
            q = it[1]
            ops.append(["qu", [pt(q.ul), pt(q.ur), pt(q.lr), pt(q.ll)]])
        else:
            raise ValueError("item inconnu %s" % k)
    return ops

def rgb(c):
    return None if c is None else [round(float(v), 6) for v in c]

ext = pg.get_drawings(extended=True)
stack = []
prims = []
for el in ext:
    L = el["level"]
    while stack and stack[-1][0] >= L:
        stack.pop()
    t = el["type"]
    if t == "clip":
        stack.append((L, "clip", conv_items(el["items"])))
    elif t == "group":
        stack.append((L, "group", float(el["opacity"])))
    elif t in ("f", "s", "fs"):
        alpha_group = 1.0
        clips = []
        for _, kind, data in stack:
            if kind == "group":
                alpha_group *= data
            else:
                clips.append(data)
        prims.append({
            "type": t,
            "fill": rgb(el.get("fill")),
            "color": rgb(el.get("color")),
            "width": float(el.get("width") or 0.0),
            "dashes": el.get("dashes"),
            "lineCap": list(el.get("lineCap") or [0]),
            "lineJoin": float(el.get("lineJoin") or 0),
            "fill_opacity": float(el["fill_opacity"]) if el.get("fill_opacity") is not None else 1.0,
            "stroke_opacity": float(el["stroke_opacity"]) if el.get("stroke_opacity") is not None else 1.0,
            "group_alpha": alpha_group,
            "even_odd": bool(el.get("even_odd", False)),
            "closePath": bool(el.get("closePath", False)),
            "clips": clips,
            "items": conv_items(el["items"]),
            "rect": [round(el["rect"].x0, 3), round(el["rect"].y0, 3), round(el["rect"].x1, 3), round(el["rect"].y1, 3)],
        })

spans = []
td = pg.get_text("dict")
for b in td["blocks"]:
    if b["type"] != 0:
        continue
    for l in b["lines"]:
        for s in l["spans"]:
            c = s["color"]
            spans.append({
                "text": s["text"], "font": s["font"], "size": round(s["size"], 3),
                "origin": [round(s["origin"][0], 3), round(s["origin"][1], 3)],
                "color": [((c >> 16) & 255) / 255., ((c >> 8) & 255) / 255., (c & 255) / 255.],
                "dir": [round(l["dir"][0], 4), round(l["dir"][1], 4)],
                "bbox": [round(v, 3) for v in s["bbox"]],
            })

json.dump({"page": [W, H], "prims": prims, "spans": spans}, open(out + ".json", "w"))
print("page %.2f x %.2f pt | %d primitives | %d textes" % (W, H, len(prims), len(spans)))
n_alpha = sum(1 for p in prims if p["group_alpha"] != 1.0 or p["fill_opacity"] != 1.0 or p["stroke_opacity"] != 1.0)
print("primitives avec opacite :", n_alpha)
