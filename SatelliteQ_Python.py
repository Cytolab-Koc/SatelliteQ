from __future__ import annotations

"""SatelliteQ analysis and report generator."""

import sys
import warnings
import html as html_lib
import re
import math
from pathlib import Path
from typing import Optional, List, Tuple
from datetime import datetime

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

SATELLITEQ_VERSION = "V1.0.0"
PCI_INNER_RADIUS_UM = 3.0
PCI_OUTER_RADIUS_UM = 12.0
PCI_STD_INNER_COL = "Std_Cell_Pericentrosomal_Inner_0_3um_Intensity"
PCI_STD_OUTER_COL = "Std_Cell_Pericentrosomal_Outer_3_12um_Intensity"
PCI_STD_TOTAL_COL = "Std_Cell_Pericentrosomal_Total_0_12um_Intensity"
PCI_CELL_INNER_COL = "Pericentrosomal_inner_0_3um_intensity"
PCI_CELL_OUTER_COL = "Pericentrosomal_outer_3_12um_intensity"
PCI_CELL_TOTAL_COL = "Pericentrosomal_total_0_12um_intensity"
PERICENTROSOMAL_CELL_METRICS = [
    "Cell_pericentrosomal_clustering_index",
    PCI_CELL_INNER_COL,
    PCI_CELL_OUTER_COL,
    PCI_CELL_TOTAL_COL,
]

SATQ_PURPLE = "#8A05FD"
SATQ_LIGHT_PURPLE = "#DCC6FF"
SATQ_PURPLE_DARK = "#4B216E"
SATQ_PURPLE_PALETTE = [
    "#7E3FB5", "#A66BDA", "#C39BEA", "#6B4FA3",
    "#B57EDC", "#8E5CC2", "#D6B7FF", "#5D3A91"
]

def choose_input_gui() -> Optional[Path]:
    try:
        import tkinter as tk
        from tkinter import filedialog, messagebox

        root = tk.Tk()
        root.withdraw()
        root.update()

        choose_folder = messagebox.askyesno(
            "SatelliteQ Python report",
            "Select a folder?\n\n"
            "Yes = select a folder containing one or more SatelliteQ outputs\n"
            "No = select one CSV/Excel file"
        )

        if choose_folder:
            selected = filedialog.askdirectory(title="Choose folder containing SatelliteQ outputs")
        else:
            selected = filedialog.askopenfilename(
                title="Choose one SatelliteQ CSV or Excel report",
                filetypes=[
                    ("SatelliteQ CSV/Excel", "*.csv *.xlsx *.xlsm"),
                    ("CSV files", "*.csv"),
                    ("Excel files", "*.xlsx *.xlsm"),
                    ("All files", "*.*"),
                ],
            )

        root.destroy()
        if not selected:
            return None
        return Path(selected)
    except Exception:
        return None

def ask_multiple_mode(n_files: int) -> str:
    """
    Return:
        combine
        separate
        both
    """
    try:
        import tkinter as tk

        choice = {"value": "combine"}

        root = tk.Tk()
        root.title("SatelliteQ multiple files")
        root.geometry("860x735")
        root.resizable(False, False)

        label = tk.Label(
            root,
            text=(
                f"Found {n_files} compatible SatelliteQ files.\n\n"
                "What should Python do?"
            ),
            justify="center",
            font=("Arial", 11),
        )
        label.pack(pady=(14, 10))

        info = tk.Label(
            root,
            text="Recommended folder structure for combined analysis",
            justify="center",
            font=("Arial", 10, "bold"),
            fg="#5A2D82",
        )
        info.pack(pady=(2, 6))

        diagram_frame = tk.Frame(root, bd=1, relief="solid", bg="white")
        diagram_frame.pack(padx=18, pady=(0, 12), fill="x")

        canvas = tk.Canvas(diagram_frame, width=790, height=385, bg="white", highlightthickness=0)
        canvas.pack(padx=10, pady=10)

        purple = "#7E57C2"
        green = "#4CAF50"
        dark = "#333333"
        light_purple = "#F7F2FF"
        light_green = "#F3FFF4"

        def rect(x1, y1, x2, y2, text_lines, outline, fill, font=("Arial", 11), bold_first=False):
            canvas.create_rectangle(x1, y1, x2, y2, outline=outline, width=2, fill=fill)
            if isinstance(text_lines, str):
                text_lines = [text_lines]
            cy = (y1 + y2) / 2
            if len(text_lines) == 1:
                canvas.create_text((x1 + x2) / 2, cy, text=text_lines[0], font=font, fill=dark)
            else:
                offset = 8 if len(text_lines) == 2 else 12
                for i, line in enumerate(text_lines):
                    if bold_first and i == 0:
                        line_font = (font[0], font[1], "bold")
                    else:
                        line_font = font
                    canvas.create_text((x1 + x2) / 2, cy - offset + i * 16, text=line, font=line_font, fill=dark)

        canvas.create_text(
            395, 18,
            text="Experiment / Condition / Replicate / Fiji_output.csv",
            font=("Arial", 12, "bold"),
            fill=dark,
        )

        # Experiment node
        exp = (25, 162, 185, 202)
        rect(*exp, "Experiment", purple, light_purple, font=("Arial", 12), bold_first=True)

        # Main vertical trunk from experiment
        trunk_x = 105
        canvas.create_line(trunk_x, 182, trunk_x, 315, width=3, fill=purple)

        # Conditions A, B, and C
        condA = (220, 60, 390, 100)
        condB = (220, 170, 390, 210)
        condC = (220, 280, 390, 320)
        canvas.create_line(trunk_x, 80, condA[0], 80, width=3, fill=purple)
        canvas.create_line(trunk_x, 190, condB[0], 190, width=3, fill=purple)
        canvas.create_line(trunk_x, 300, condC[0], 300, width=3, fill=purple)
        rect(*condA, "Condition_A", purple, light_purple, font=("Arial", 12))
        rect(*condB, "Condition_B", purple, light_purple, font=("Arial", 12))
        rect(*condC, "Condition_C", purple, light_purple, font=("Arial", 12))

        # Replicates under each condition
        repA1 = (430, 38, 590, 78)
        repA2 = (430, 83, 590, 123)
        repB1 = (430, 148, 590, 188)
        repB2 = (430, 193, 590, 233)
        repC1 = (430, 258, 590, 298)
        repC2 = (430, 303, 590, 343)

        def branch_replicates(cond, rep1, rep2):
            branch_x = 410
            center_y = (cond[1] + cond[3]) / 2
            rep1_y = (rep1[1] + rep1[3]) / 2
            rep2_y = (rep2[1] + rep2[3]) / 2
            canvas.create_line(cond[2], center_y, branch_x, center_y, width=3, fill=purple)
            canvas.create_line(branch_x, rep1_y, branch_x, rep2_y, width=3, fill=purple)
            canvas.create_line(branch_x, rep1_y, rep1[0], rep1_y, width=3, fill=purple)
            canvas.create_line(branch_x, rep2_y, rep2[0], rep2_y, width=3, fill=purple)

        branch_replicates(condA, repA1, repA2)
        branch_replicates(condB, repB1, repB2)
        branch_replicates(condC, repC1, repC2)

        for rep in [repA1, repB1, repC1]:
            rect(*rep, "Replicate_1", green, light_green, font=("Arial", 12))
        for rep in [repA2, repB2, repC2]:
            rect(*rep, "Replicate_2", green, light_green, font=("Arial", 12))

        # Fiji output files
        outA1 = (660, 38, 755, 78)
        outA2 = (660, 83, 755, 123)
        outB1 = (660, 148, 755, 188)
        outB2 = (660, 193, 755, 233)
        outC1 = (660, 258, 755, 298)
        outC2 = (660, 303, 755, 343)

        for rep, outnode in [(repA1, outA1), (repA2, outA2), (repB1, outB1), (repB2, outB2), (repC1, outC1), (repC2, outC2)]:
            y = (rep[1] + rep[3]) / 2
            canvas.create_line(rep[2], y, outnode[0], y, width=3, fill=green)
            rect(*outnode, ["Fiji_output", ".csv/.xlsx"], green, light_green, font=("Arial", 9), bold_first=True)

        canvas.create_text(
            395, 365,
            text="Tip: put each Fiji output inside its replicate folder, inside its condition folder.",
            font=("Arial", 10),
            fill="#555555",
        )

        def set_choice(v: str):
            choice["value"] = v
            root.destroy()

        btn_frame = tk.Frame(root)
        btn_frame.pack(pady=8)

        tk.Button(btn_frame, text="Combine into one report", width=28, command=lambda: set_choice("combine")).pack(pady=4)
        tk.Button(btn_frame, text="Analyze separately", width=28, command=lambda: set_choice("separate")).pack(pady=4)
        tk.Button(btn_frame, text="Do both", width=28, command=lambda: set_choice("both")).pack(pady=4)

        root.mainloop()
        return choice["value"]
    except Exception:
        print(f"Found {n_files} compatible files.")
        print("Choose mode:")
        print("  1 = combine into one report")
        print("  2 = analyze separately")
        print("  3 = do both")
        print("Recommended structure: Experiment/Condition/Replicate/Fiji_output.csv")
        try:
            ans = input("Enter 1, 2, or 3: ").strip()
        except EOFError:
            ans = "1"
        if ans == "2":
            return "separate"
        if ans == "3":
            return "both"
        return "combine"


def remove_date_time_tokens(text: str) -> str:
    """Remove common microscope/Fiji date-time stamps from generated file names.

    This only affects output file/folder names. Full original file names remain
    preserved inside the Excel report and source columns.
    """
    s = str(text)

    patterns = [

        r"(?<!\d)(20\d{2}|19\d{2})[._ -]?\d{1,2}[._ -]?\d{1,2}[._ -]?\d{1,2}[._ -]?\d{1,2}[._ -]?\d{1,2}(?!\d)",

        r"(?<!\d)(20\d{2}|19\d{2})[._ -]?\d{1,2}[._ -]?\d{1,2}(?!\d)",

        r"(?<!\d)\d{1,2}[._ -]?\d{1,2}[._ -]?(20\d{2}|19\d{2})(?!\d)",

        r"(?:(?<=_)|(?<=-)|(?<=\s)|^)\d{1,2}[._-]\d{2}[._-]\d{2}(?=(_|-|\s|$))",

        r"(?<!\d)(20\d{2}|19\d{2})_\d{1,2}_\d{1,2}_\d{1,2}_\d{1,2}_\d{1,2}(?!\d)",
    ]

    for pat in patterns:
        s = re.sub(pat, "_", s)

    s = re.sub(r"(?i)(date|time|timestamp|runstamp|run_timestamp)", "_", s)

    return s

def safe_name(text: str) -> str:
    out = remove_date_time_tokens(str(text))
    for ch in ['\\', '/', ':', '*', '?', '"', '<', '>', '|', ' ', '.', '(', ')', '[', ']']:
        out = out.replace(ch, "_")
    while "__" in out:
        out = out.replace("__", "_")
    return out.strip("_")[:120] or "SatelliteQ"


def satelliteq_date_stamp() -> str:
    return datetime.now().strftime("%Y_%m_%d")

def satelliteq_report_base_name(outdir: Path, inputs: Optional[List[Path]] = None, combined: bool = False) -> str:
    """Return a clean experiment-like base name for output files."""
    if combined:
        base = outdir.parent.name if outdir.name == "SatelliteQ_Combined_Report" else outdir.name
    elif inputs:
        base = Path(inputs[0]).stem
    else:
        base = outdir.name
        for suffix in ["_Python_Report", "_SatelliteQ_Combined_Report", "_Combined_Report"]:
            if base.endswith(suffix):
                base = base[: -len(suffix)]
    return safe_name(base) or "SatelliteQ"

def satelliteq_dated_name(outdir: Path, label: str, inputs: Optional[List[Path]] = None, combined: bool = False, ext: str = "") -> str:
    return f"{satelliteq_report_base_name(outdir, inputs, combined)}_{satelliteq_date_stamp()}_{label}{ext}"

def find_candidate_files(path: Path) -> List[Path]:
    if path.is_file():
        return [path]

    patterns = [
        # 2D Fiji outputs
        "**/*_SatelliteQuantify_Fiji_Output.csv",
        "**/*_SatelliteQuantify_Fiji_Output.xlsx",
        "**/*_SatelliteQuantify_Fiji_Output.xlsm",

        # Old 3D Fiji outputs
        "**/satellite_3D_results_clean.csv",
        "**/satellite_3D_results_clean.xlsx",
        "**/satellite_3D_results_clean.xlsm",
        "**/satellite_3D_results_clean_excel_friendly.tsv",
        "**/satellite_3D_results_clean_excel_friendly.xlsx",
        "**/satellite_3D_results_clean_excel_friendly.xlsm",

        # 3D Fiji outputs, for example:
        # ImageName_3D_2026_7_8_16_42_10_results.csv/xlsx
        "**/*_3D_*_results.csv",
        "**/*_3D_*_results.tsv",
        "**/*_3D_*_results.xlsx",
        "**/*_3D_*_results.xlsm",
        "**/*_3D_*_results_excel_friendly.tsv",
        "**/*_3D_*_results_excel_friendly.xlsx",
        "**/*_3D_*_results_excel_friendly.xlsm",

        # Manually renamed 3D result Excel files
        "**/*3D*results*.csv",
        "**/*3D*results*.tsv",
        "**/*3D*results*.xlsx",
        "**/*3D*results*.xlsm",

        # Time-lapse outputs
        "**/*timelapse_satellite_rows_combined.csv",
        "**/*timelapse_satellite_rows_combined.xlsx",

    ]

    out: List[Path] = []
    seen = set()

    for pat in patterns:
        for f in sorted(path.glob(pat)):
            if not f.is_file():
                continue
            low = str(f).lower()
            if any(skip in low for skip in [
                "~$",
                "satelliteq_combined_report",
                "combined_quantification",
                "combined_python_report",
                "quantification_report",
                "_python_report",
                "graph_files",
                "_csv_note",
                "_selected_options_summary",
                "review_label_note",
                "ciliumquantification",
                "ciliumcellsummary",
                "ciliumskeletonpoints",
            ]):
                continue
            if f not in seen:
                out.append(f)
                seen.add(f)

    return out



def _meta_first_existing(columns: List[str], candidates: List[str]) -> Optional[str]:
    lower = {str(c).strip().lower(): c for c in columns}
    for c in candidates:
        if c in columns:
            return c
        if c.lower() in lower:
            return lower[c.lower()]
    return None

def find_experiment_metadata_file(root: Path) -> Optional[Path]:
    """Find a user-provided condition/replicate annotation table."""
    roots = []
    if root.is_file():
        roots.append(root.parent)
    else:
        roots.append(root)
    if roots[0].parent not in roots:
        roots.append(roots[0].parent)

    names = [
        "SatelliteQ_experiment_metadata.csv",
        "satelliteq_experiment_metadata.csv",
        "experiment_metadata.csv",
        "conditions.csv",
    ]
    for r in roots:
        for name in names:
            p = r / name
            if p.exists():
                return p

    try:
        for p in roots[0].glob("*metadata*.csv"):
            if "experiment" in p.name.lower() or "condition" in p.name.lower():
                return p
    except Exception:
        pass
    return None

def looks_like_replicate_folder(name: str) -> bool:
    """Detect common replicate folder names such as Rep1, Replicate 1, R2."""
    s = str(name).strip().lower().replace("_", " ").replace("-", " ")
    compact = s.replace(" ", "")
    if not s:
        return False
    if "replicate" in s or "biological replicate" in s:
        return True
    if compact.startswith("rep") and any(ch.isdigit() for ch in compact):
        return True
    if compact.startswith("r") and len(compact) <= 4 and any(ch.isdigit() for ch in compact):
        return True
    return False

def looks_like_frame_token(name: str) -> bool:
    """Detect common time-lapse frame/timepoint folder or file tokens."""
    txt = str(name).strip().lower()
    if not txt:
        return False
    stem = Path(txt).stem if any(txt.endswith(ext) for ext in [".csv", ".xlsx", ".xls", ".xlsm"]) else txt
    stem = stem.replace("-", "_").replace(" ", "_")
    patterns = [
        r"^(frame|frames|frm|fr|f)[_]*\d+$",
        r"^(timepoint|time_point|tp|t)[_]*\d+$",
        r"^.*[_](frame|frm|fr|f)[_]*\d+.*$",
        r"^.*[_](timepoint|time_point|tp|t)[_]*\d+.*$",
        r"^t\d{2,4}$",
    ]
    return any(re.match(pat, stem) for pat in patterns)

def infer_frame_from_path(file_path: Path) -> Optional[int]:
    """Infer a frame number from a time-lapse file or parent folder name when the CSV lacks a Frame column."""
    file_path = Path(file_path)
    candidates = [file_path.stem] + [part for part in reversed(file_path.parent.parts[-4:])]
    patterns = [
        r"(?:^|[_\-\s])frame[_\-\s]*(\d+)(?:$|[_\-\s])",
        r"(?:^|[_\-\s])frm[_\-\s]*(\d+)(?:$|[_\-\s])",
        r"(?:^|[_\-\s])fr[_\-\s]*(\d+)(?:$|[_\-\s])",
        r"(?:^|[_\-\s])timepoint[_\-\s]*(\d+)(?:$|[_\-\s])",
        r"(?:^|[_\-\s])time_point[_\-\s]*(\d+)(?:$|[_\-\s])",
        r"(?:^|[_\-\s])tp[_\-\s]*(\d+)(?:$|[_\-\s])",
        r"(?:^|[_\-\s])t[_\-\s]*(\d{1,4})(?:$|[_\-\s])",
        r"^t(\d{1,4})$",
        r"^f(\d{1,4})$",
    ]
    for cand in candidates:
        c = str(cand).strip().lower().replace(".", "_")
        for pat in patterns:
            m = re.search(pat, c)
            if m:
                try:
                    return int(m.group(1))
                except Exception:
                    pass
    return None

def remove_frame_tokens_from_parts(parts: List[str]) -> List[str]:
    """Remove frame/timepoint folder levels so frame folders are not treated as conditions or replicates."""
    return [p for p in parts if not looks_like_frame_token(p)]

def infer_condition_replicate_from_path(file_path: Path, root: Optional[Path] = None) -> Tuple[str, str]:
    """Infer Condition and Replicate from folder structure.

    Supported structures:
        Experiment/WT/Replicate 1/output.csv
        Experiment/WT/Replicate 2/output.csv
        Experiment/Condition 1/Replicate 1/output.csv

    Also supported:
        Experiment/Replicate 1/WT/output.csv
        Experiment/Replicate 1/Condition 1/output.csv

    Fallback:
        Condition = parent folder
        Replicate = file stem
    """
    file_path = Path(file_path)
    parent = file_path.parent.name if file_path.parent else ""
    grandparent = file_path.parent.parent.name if file_path.parent and file_path.parent.parent else ""
    fallback_condition = parent if parent else "Condition 1"
    fallback_replicate = file_path.stem if file_path.stem else "Replicate 1"

    parts = []
    if root is not None:
        try:
            root = Path(root).resolve()
            rel_parent = file_path.resolve().parent.relative_to(root)
            parts = [p for p in rel_parent.parts if p not in ["", "."]]
        except Exception:
            parts = []

    if parts:
        parts = remove_frame_tokens_from_parts(parts)

    if len(parts) >= 2:
        first = parts[0]
        last = parts[-1]
        second_last = parts[-2]

        if looks_like_replicate_folder(first) and not looks_like_replicate_folder(last):
            return last, first

        if looks_like_replicate_folder(last):
            return second_last, last

        return second_last, last

    if len(parts) == 1:
        return parts[0], fallback_replicate

    if looks_like_frame_token(parent) and grandparent:

        gg = file_path.parent.parent.parent.name if file_path.parent and file_path.parent.parent and file_path.parent.parent.parent else ""
        if looks_like_replicate_folder(grandparent) and gg:
            return gg, grandparent
        return grandparent, fallback_replicate
    if looks_like_replicate_folder(parent) and grandparent:
        return grandparent, parent
    if looks_like_replicate_folder(grandparent) and parent:
        return parent, grandparent

    if looks_like_frame_token(file_path.stem) and parent:
        return parent, grandparent if looks_like_replicate_folder(grandparent) else "Replicate 1"

    return fallback_condition, fallback_replicate

def experiment_template_df(inputs: List[Path], df: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    rows = []
    common_root = None
    if inputs:
        try:
            common_root = Path(__import__("os").path.commonpath([str(p.parent.resolve()) for p in inputs]))
        except Exception:
            common_root = inputs[0].parent

    for f in inputs:
        condition, replicate = infer_condition_replicate_from_path(f, common_root)
        rows.append({
            "Input_file": f.name,
            "File_stem": f.stem,
            "Condition": condition,
            "Replicate": replicate,
            "Notes": "Folder-based labels were inferred automatically.",
        })
    if not rows and df is not None and "Source_Input" in df.columns:
        for src_name, g in df.groupby("Source_Input", dropna=False):
            rows.append({
                "Input_file": src_name,
                "File_stem": Path(str(src_name)).stem,
                "Condition": g["Condition"].iloc[0] if "Condition" in g.columns else "",
                "Replicate": g["Replicate"].iloc[0] if "Replicate" in g.columns else "",
                "Notes": "Folder-based labels are inferred automatically.",
            })
    return pd.DataFrame(rows)

def _metadata_match_keys(path_or_name) -> set:
    p = Path(str(path_or_name))
    s = str(path_or_name)
    keys = {
        s,
        s.replace("\\", "/"),
        p.name,
        p.stem,
        str(p),
        str(p).replace("\\", "/"),
    }
    if p.parent.name:
        keys.add(p.parent.name + "/" + p.name)
        keys.add(p.parent.name + "/" + p.stem)
    return {k.strip().lower() for k in keys if str(k).strip()}

def apply_experiment_metadata(df: pd.DataFrame, inputs: List[Path], root: Path) -> Tuple[pd.DataFrame, Optional[Path]]:
    """Apply Condition and Replicate columns from a user metadata CSV if present.

    Metadata file format:
        Input_file,Condition,Replicate
        my_file_1.csv,WT,rep1
        my_file_2.csv,WT,rep2
        my_file_3.csv,Treatment,rep1
    """
    out = df.copy()
    if "Condition" not in out.columns:
        out["Condition"] = out["Source_Folder"] if "Source_Folder" in out.columns else "Condition 1"
    if "Replicate" not in out.columns:
        out["Replicate"] = out["Source_Input"].astype(str).map(lambda x: Path(x).stem) if "Source_Input" in out.columns else "Replicate 1"

    meta_path = find_experiment_metadata_file(root)
    if meta_path is None:

        inferred_conditions = []
        inferred_replicates = []
        for _, r in out.iterrows():
            src_path = str(r.get("Source_Path", ""))
            if src_path and Path(src_path).exists():
                cond, rep = infer_condition_replicate_from_path(Path(src_path), root)
            else:
                cond = r.get("Condition", "")
                rep = r.get("Replicate", "")
            inferred_conditions.append(cond)
            inferred_replicates.append(rep)
        if inferred_conditions:
            out["Condition"] = inferred_conditions
            out["Replicate"] = inferred_replicates
        return out, None

    try:
        meta = read_csv_robust(meta_path)
    except Exception as e:
        print(f"Could not read experiment metadata file {meta_path}: {e}")
        return out, None

    cols = list(meta.columns)
    file_col = _meta_first_existing(cols, ["Input_file", "File", "Filename", "Source_Input", "Source_File", "Full_path", "Path"])
    cond_col = _meta_first_existing(cols, ["Condition", "Group", "Treatment", "Genotype", "Experiment_condition"])
    rep_col = _meta_first_existing(cols, ["Replicate", "Rep", "Biological_replicate", "Experiment_replicate"])

    if file_col is None or cond_col is None:
        print(f"Metadata file found but missing Input_file/File and Condition columns: {meta_path}")
        return out, meta_path

    mapping = {}
    for _, row in meta.iterrows():
        key_raw = str(row.get(file_col, "")).strip()
        if not key_raw:
            continue
        cond = row.get(cond_col, "")
        rep = row.get(rep_col, "") if rep_col else ""
        for k in _metadata_match_keys(key_raw):
            mapping[k] = (cond, rep)

    input_by_name = {}
    for f in inputs:
        for k in _metadata_match_keys(f):
            input_by_name[k] = f

    def row_lookup(src_input: str, src_path: str):
        keys = _metadata_match_keys(src_input)
        if src_path:
            keys |= _metadata_match_keys(src_path)
        for k in keys:
            if k in mapping:
                return mapping[k]
        return None

    conds = []
    reps = []
    for _, r in out.iterrows():
        match = row_lookup(str(r.get("Source_Input", "")), str(r.get("Source_Path", "")))
        if match is not None:
            cond, rep = match
            conds.append(cond if str(cond).strip() else r.get("Condition", ""))
            reps.append(rep if str(rep).strip() else r.get("Replicate", ""))
        else:
            conds.append(r.get("Condition", ""))
            reps.append(r.get("Replicate", ""))

    out["Condition"] = conds
    out["Replicate"] = reps
    return out, meta_path

def read_csv_robust(path: Path) -> pd.DataFrame:
    try:
        return pd.read_csv(path, sep=None, engine="python")
    except Exception:
        try:
            return pd.read_csv(path)
        except Exception:
            return pd.read_csv(path, sep="\t")

def first_existing(columns: List[str], candidates: List[str]) -> Optional[str]:
    lower_map = {str(c).lower(): str(c) for c in columns}
    for c in candidates:
        if c in columns:
            return c
        if c.lower() in lower_map:
            return lower_map[c.lower()]
    return None


def _normalized_column_key(name: str) -> str:
    """Normalize common Fiji/Excel unit spelling differences for dimensional columns."""
    s = str(name).strip().lower()
    s = s.replace("µ", "u").replace("μ", "u")
    s = s.replace("²", "2").replace("³", "3")
    s = s.replace("micrometers", "um").replace("micrometer", "um")
    s = s.replace("microns", "um").replace("micron", "um")
    s = re.sub(r"[\s_().\[\]{}\-]+", "", s)
    s = s.replace("^", "")
    return s


def first_existing_normalized(columns: List[str], candidates: List[str]) -> Optional[str]:
    """Find a column while tolerating spaces, underscores, µ/um and exponent formatting."""
    exact = first_existing(columns, candidates)
    if exact is not None:
        return exact

    normalized = {}
    for col in columns:
        normalized.setdefault(_normalized_column_key(col), col)

    for candidate in candidates:
        key = _normalized_column_key(candidate)
        if key in normalized:
            return normalized[key]
    return None


def _numeric_dimension_series(df: pd.DataFrame, candidates: List[str]) -> pd.Series:
    """Return first usable numeric dimensional measurement across known aliases.

    Values are filled row-by-row, rather than selecting one global source column.
    This is important when an already-combined table contains both 2D and 3D rows,
    or when a standardized column exists but is blank while the original Fiji
    measurement column is still populated.
    """
    result = pd.Series(np.nan, index=df.index, dtype=float)

    # Build candidate order from exact/normalized matches, without duplicates.
    matched = []
    seen = set()
    cols = list(df.columns)
    for candidate in candidates:
        col = first_existing_normalized(cols, [candidate])
        if col is not None and col not in seen:
            matched.append(col)
            seen.add(col)

    for col in matched:
        raw = df[col]
        numeric_values = pd.to_numeric(raw, errors="coerce")

        # Limited decimal-comma fallback for dimensional values only.
        if numeric_values.notna().sum() == 0 and raw.dtype == object:
            cleaned = (
                raw.astype(str)
                .str.strip()
                .str.replace("\u00a0", "", regex=False)
                .str.replace(",", ".", regex=False)
            )
            numeric_values = pd.to_numeric(cleaned, errors="coerce")

        fill_mask = result.isna() & numeric_values.notna()
        result.loc[fill_mask] = numeric_values.loc[fill_mask]

    return result

def to_numeric_if_possible(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for col in out.columns:
        if out[col].dtype == object:
            converted = pd.to_numeric(out[col], errors="coerce")
            if converted.notna().sum() > 0 and converted.notna().sum() >= max(3, int(0.5 * len(out))):
                out[col] = converted
    return out

def detect_mode(df: pd.DataFrame, path: Path) -> str:
    cols = list(df.columns)

    if "Mode" in cols:
        vals = df["Mode"].dropna().astype(str).str.strip()
        vals = vals[vals.isin(["2D", "3D"])]
        if not vals.empty:
            unique_modes = vals.unique()
            if len(unique_modes) == 1:
                return str(unique_modes[0])

    # Prefer actual dimensional measurements over file/folder names.
    volume_probe = _numeric_dimension_series(
        df,
        ["Std_Volume_3D", "Volume_um3", "Volume (micron^3)", "Volume (µm³)",
         "Volume (µm^3)", "Volume_um^3", "Volume"]
    )
    area_probe = _numeric_dimension_series(
        df,
        ["Std_Area_2D", "Satellite_Area_um2", "Area_um2", "Area (micron^2)",
         "Area (µm²)", "Area (µm^2)", "Area_um^2", "Area"]
    )

    z_col = first_existing_normalized(cols, ["Std_Z", "Z_um", "Z", "ZM"])
    z_has_values = False
    if z_col is not None:
        z_has_values = pd.to_numeric(df[z_col], errors="coerce").notna().any()

    if volume_probe.notna().any() or z_has_values:
        return "3D"
    if area_probe.notna().any():
        return "2D"

    # Filename/path is now only a last-resort fallback.
    if "3d" in Path(path).name.lower():
        return "3D"
    return "2D"

def detect_analysis_type(df: pd.DataFrame, path: Path) -> str:
    cols = list(df.columns)
    if "Analysis_Type" in cols and df["Analysis_Type"].dropna().size:
        return str(df["Analysis_Type"].dropna().iloc[0])
    low = str(path).lower()
    if "timelapse" in low or "time_lapse" in low or infer_frame_from_path(path) is not None or first_existing(cols, ["Frame", "Std_Frame"]):
        return "timelapse_rows"
    return "single"

def standardize_rows(df_raw: pd.DataFrame, source_path: Path) -> Tuple[pd.DataFrame, str, str]:
    df = df_raw.copy()

    # Preserve image/calibration metadata from Fiji Setting rows before they are
    # removed from the satellite-level table. These values are used by the
    # cell/image-center distance calculation.
    setting_values = {}
    raw_cols = list(df.columns)
    raw_row_col = first_existing(raw_cols, ["Row_Type"])
    raw_setting_col = first_existing(raw_cols, ["Setting"])
    raw_value_col = first_existing(raw_cols, ["Value"])
    if raw_row_col is not None and raw_setting_col is not None and raw_value_col is not None:
        try:
            smask = df[raw_row_col].astype(str).str.strip().str.lower().eq("setting")
            for _, sr in df.loc[smask, [raw_setting_col, raw_value_col]].iterrows():
                skey = str(sr.get(raw_setting_col, "")).strip().lower()
                if skey:
                    setting_values[skey] = sr.get(raw_value_col, np.nan)
        except Exception:
            setting_values = {}

    row_col = first_existing(list(df.columns), ["Row_Type"])
    if row_col is not None:
        df = df[df[row_col].astype(str).str.lower().eq("satellite")].copy()

    df = df.dropna(how="all").copy()
    df = to_numeric_if_possible(df)

    mode = detect_mode(df, source_path)
    analysis_type = detect_analysis_type(df, source_path)

    if "Source_File" not in df.columns:
        df.insert(0, "Source_File", source_path.name)

    df["Source_Input"] = source_path.name
    df["Source_Path"] = str(source_path)

    df["Source_UID"] = str(source_path.resolve()) if hasattr(source_path, "resolve") else str(source_path)
    df["Source_Folder"] = source_path.parent.name

    def _setting_number(*keys):
        for key in keys:
            value = setting_values.get(str(key).strip().lower(), np.nan)
            try:
                num = float(str(value).strip().replace(",", "."))
                if np.isfinite(num):
                    return num
            except Exception:
                pass
        return np.nan

    # baseline; coordinate extent is used when these values are absent.
    for _col, _keys in {
        "Std_Image_Width_Pixels": ("Image width pixels", "Image width", "Width pixels"),
        "Std_Image_Height_Pixels": ("Image height pixels", "Image height", "Height pixels"),
        "Std_Pixel_Width_um": ("Pixel width um", "Pixel width µm", "Pixel width"),
        "Std_Pixel_Height_um": ("Pixel height um", "Pixel height µm", "Pixel height"),
        "Std_Voxel_Depth_um": ("Voxel depth um", "Z-step size", "Z step", "Voxel depth"),
    }.items():
        if _col not in df.columns or pd.to_numeric(df[_col], errors="coerce").notna().sum() == 0:
            df[_col] = _setting_number(*_keys)

    # 3D outputs historically stored calibration in the sibling selected-options
    # text file rather than Setting rows. Read it when needed.
    if any(pd.to_numeric(df[c], errors="coerce").notna().sum() == 0 for c in ["Std_Image_Width_Pixels", "Std_Image_Height_Pixels", "Std_Pixel_Width_um", "Std_Pixel_Height_um", "Std_Voxel_Depth_um"]):
        try:
            option_files = sorted(Path(source_path).parent.glob("*_SELECTED_OPTIONS_SUMMARY.txt"))
            if option_files:
                opt_text = option_files[0].read_text(encoding="utf-8", errors="replace")
                opt_map = {}
                for line in opt_text.splitlines():
                    if ":" in line:
                        k, v = line.split(":", 1)
                        opt_map[k.strip().lower()] = v.strip()
                key_map = {
                    "Std_Image_Width_Pixels": ["image width pixels", "image width"],
                    "Std_Image_Height_Pixels": ["image height pixels", "image height"],
                    "Std_Pixel_Width_um": ["pixel width um", "pixel width µm"],
                    "Std_Pixel_Height_um": ["pixel height um", "pixel height µm"],
                    "Std_Voxel_Depth_um": ["voxel depth um", "z-step size", "z step", "voxel depth"],
                }
                for c, keys in key_map.items():
                    if pd.to_numeric(df[c], errors="coerce").notna().sum() > 0:
                        continue
                    for k in keys:
                        if k in opt_map:
                            try:
                                val = float(opt_map[k].replace(",", "."))
                                if np.isfinite(val):
                                    df[c] = val
                                    break
                            except Exception:
                                pass
        except Exception:
            pass

    if "Condition" not in df.columns:
        df["Condition"] = source_path.parent.name
    if "Replicate" not in df.columns:
        df["Replicate"] = source_path.stem

    if "Mode" in df.columns:
        existing_mode = df["Mode"].astype(str).str.strip()
        valid_mode = existing_mode.isin(["2D", "3D"])
        df["Mode"] = existing_mode.where(valid_mode, mode)
    else:
        df["Mode"] = mode

    df["Analysis_Type"] = analysis_type

    cols = list(df.columns)

    sat_id = first_existing(cols, ["Std_Satellite_ID", "Satellite_ID", "Object", "Obj", "ID"])
    cell_id = first_existing(cols, ["Std_Cell_ID", "Cell_ID"])
    cent_id = first_existing(cols, ["Std_Centrosome_ID", "Assigned_Centrosome_ID", "Nearest_Centrosome_ID", "Centrosome_ID"])
    frame = first_existing(cols, ["Std_Frame", "Frame"])
    time_min = first_existing(cols, ["Std_Time_min", "Time_min"])

    # Keep dimensional measurements in separate standardized columns.
    # Values are recovered independently from all known Fiji aliases. This means
    # area and volume remain available even if Mode was inferred incorrectly,
    # a parent folder contains "3D", or a previous standardized column is blank.
    area_2d_candidates = [
        "Satellite_Area_um2", "Area_um2", "Area (micron^2)", "Area (µm²)",
        "Area (µm^2)", "Area_um^2", "Area", "Std_Area_2D"
    ]
    volume_3d_candidates = [
        "Volume_um3", "Volume (micron^3)", "Volume (µm³)", "Volume (µm^3)",
        "Volume_um^3", "Volume", "Std_Volume_3D"
    ]
    cell_area_2d_candidates = [
        "Cell_Area_um2", "Cell_Area", "Cell Area (micron^2)", "Cell Area (µm²)",
        "Cell_Area_um^2", "Std_Cell_Area", "Std_Cell_Area_2D"
    ]
    cell_volume_3d_candidates = [
        "Cell_Volume_um3", "Cell_Volume_3D", "Cell_Volume",
        "Cell Volume (micron^3)", "Cell Volume (µm³)", "Cell_Volume_um^3",
        "Std_Cell_Volume", "Std_Cell_Volume_3D"
    ]

    if mode == "3D":
        size = first_existing_normalized(cols, ["Std_Size", "Volume_um3", "Volume (micron^3)", "Volume (µm³)", "Volume (µm^3)", "Volume_um^3", "Volume"])
        surface = first_existing(cols, ["Std_Surface", "Surface_um2", "Surface (micron^2)", "Surface (µm²)", "Surface_Area_um2", "Surface Area", "Surface"])
        distance = first_existing(cols, ["Std_Distance", "Distance_to_Assigned_Centrosome_3D_um", "Distance_to_Assigned_Centrosome_3D", "Distance_to_nearest_centrosome_3D_um", "Distance_to_nearest_centrosome_3D"])
        diameter = first_existing(cols, ["Std_Diameter", "Equivalent_Sphere_Diameter_um"])

        # Prefer the dedicated morphology sphericity reported by 3D Objects Counter+.
        # Sphericity_3D values are not preferred because they may
        # have been calculated from an ambiguous generic "Surface" column.
        sphericity = first_existing(cols, [
            "Morph_Sphericity", "Morph Sphericity", "Morph.Sphericity",
            "MorphSphericity", "Std_Sphericity"
        ])

        # Formula fallback is allowed only when a column explicitly represents
        # calibrated surface AREA. Generic legacy "Surface" is excluded.
        sphericity_surface_area = first_existing(cols, [
            "Surface_um2", "Surface (micron^2)", "Surface (µm²)",
            "Surface_Area_um2", "Surface area", "Surface Area",
            "Morph_SurfaceArea_um2", "Morph_Surface_Area_um2",
            "Morph_SurfaceArea", "Morph_Surface_Area"
        ])

        circularity_2d = first_existing(cols, ["Std_Circularity_2D", "Circularity_2D", "Circ.", "Circularity"])
        cell_measure = first_existing_normalized(cols, ["Std_Cell_Measure", "Std_Cell_Volume", "Cell_Volume_um3", "Cell_Volume_3D", "Cell_Volume"])
        cell_measure_name = "Cell volume (µm³)"
    else:
        sphericity_surface_area = None
        size = first_existing_normalized(cols, ["Std_Size", "Satellite_Area_um2", "Area_um2", "Area (micron^2)", "Area (µm²)", "Area (µm^2)", "Area_um^2", "Area"])
        surface = first_existing(cols, ["Std_Surface"])
        distance = first_existing(cols, ["Std_Distance", "Distance_to_Nearest_Centrosome_um", "Nearest_Centrosome_Distance_um", "Nearest_Centrosome_Distance", "Distance_to_Assigned_Centrosome", "Distance_to_nearest_centrosome"])
        diameter = first_existing(cols, ["Std_Diameter", "Equivalent_Circle_Diameter_um", "Equivalent_Circle_Diameter"])
        sphericity = first_existing(cols, ["Std_Sphericity", "Sphericity_3D", "Sphericity"])
        circularity_2d = first_existing(cols, ["Std_Circularity_2D", "Circularity_2D", "Circ.", "Circularity"])
        cell_measure = first_existing_normalized(cols, ["Std_Cell_Measure", "Std_Cell_Area", "Cell_Area_um2", "Cell_Area"])
        cell_measure_name = "Cell area (µm²)"

    intensity = first_existing(cols, ["Std_Intensity", "Integrated_Intensity_AU", "Integrated_Intensity", "IntDen", "Integrated density", "Integrated_Density", "Mean_Intensity_AU", "Mean_Intensity", "Mean"])
    mean_intensity = first_existing(cols, ["Std_Mean_Intensity", "Mean_Intensity_AU", "Mean_Intensity", "Mean"])
    x = first_existing(cols, ["Std_X", "X_um", "X", "XM"])
    y = first_existing(cols, ["Std_Y", "Y_um", "Y", "YM"])
    z = first_existing(cols, ["Std_Z", "Z_um", "Z", "ZM"])
    cent_x = first_existing(cols, ["Std_Centrosome_X", "Assigned_Centrosome_X_um", "Nearest_Centrosome_X_um", "Centrosome_X_um", "Assigned_Centrosome_X", "Nearest_Centrosome_X", "Centrosome_X", "Centroid_X_um"])
    cent_y = first_existing(cols, ["Std_Centrosome_Y", "Assigned_Centrosome_Y_um", "Nearest_Centrosome_Y_um", "Centrosome_Y_um", "Assigned_Centrosome_Y", "Nearest_Centrosome_Y", "Centrosome_Y", "Centroid_Y_um"])
    cent_z = first_existing(cols, ["Std_Centrosome_Z", "Assigned_Centrosome_Z_um", "Nearest_Centrosome_Z_um", "Centrosome_Z_um", "Assigned_Centrosome_Z", "Nearest_Centrosome_Z", "Centrosome_Z", "Centroid_Z_um"])

    def copy_std(std_col: str, src: Optional[str], numeric: bool = True):
        if src and src in df.columns:
            if numeric:
                df[std_col] = pd.to_numeric(df[src], errors="coerce")
            else:
                df[std_col] = df[src]
        else:
            df[std_col] = np.nan

    copy_std("Std_Satellite_ID", sat_id, numeric=False)
    if df["Std_Satellite_ID"].isna().all():
        df["Std_Satellite_ID"] = np.arange(1, len(df) + 1)

    copy_std("Std_Cell_ID", cell_id, numeric=False)
    copy_std("Std_Centrosome_ID", cent_id, numeric=False)
    copy_std("Std_Frame", frame, numeric=True)
    inferred_frame = infer_frame_from_path(source_path)
    if inferred_frame is not None and df["Std_Frame"].dropna().empty:
        df["Std_Frame"] = float(inferred_frame)
        df["Analysis_Type"] = "timelapse_rows"
    copy_std("Std_Time_min", time_min, numeric=True)

    # Recover dimensional measurements directly from Fiji output columns.
    # Do NOT blank one measurement merely because the file was classified as
    # 2D or 3D. The presence of the actual Fiji measurement is the source of truth.
    df["Std_Area_2D"] = _numeric_dimension_series(df, area_2d_candidates)
    df["Std_Volume_3D"] = _numeric_dimension_series(df, volume_3d_candidates)

    # Std_Size is the unified area/volume column.
    # Prefer the correct dimensional measurement row-by-row.
    legacy_size = pd.to_numeric(df[size], errors="coerce") if size and size in df.columns else pd.Series(np.nan, index=df.index)
    if "Mode" in df.columns:
        mode_text = df["Mode"].astype(str).str.strip()
    else:
        mode_text = pd.Series(mode, index=df.index)
    df["Std_Size"] = np.where(
        mode_text.eq("3D"),
        df["Std_Volume_3D"].combine_first(legacy_size),
        df["Std_Area_2D"].combine_first(legacy_size),
    )
    df["Std_Size"] = pd.to_numeric(df["Std_Size"], errors="coerce")

    copy_std("Std_Surface", surface, numeric=True)
    copy_std("Std_Intensity", intensity, numeric=True)

    df["Std_Total_Intensity"] = pd.to_numeric(df["Std_Intensity"], errors="coerce")
    copy_std("Std_Mean_Intensity", mean_intensity, numeric=True)
    copy_std("Std_Distance", distance, numeric=True)
    copy_std("Std_Diameter", diameter, numeric=True)
    copy_std("Std_Sphericity", sphericity, numeric=True)
    copy_std("Std_Circularity_2D", circularity_2d, numeric=True)

    # Keep a transparent record of where 3D sphericity came from.
    df["Std_Sphericity_Source"] = ""
    if mode == "3D" and sphericity and sphericity in df.columns:
        has_source = pd.to_numeric(df["Std_Sphericity"], errors="coerce").notna()
        df.loc[has_source, "Std_Sphericity_Source"] = str(sphericity)

    if mode == "3D":
        # Valid sphericity is dimensionless and bounded from 0 to 1.
        # Invalid values are marked missing rather than silently clipped.
        sph = pd.to_numeric(df["Std_Sphericity"], errors="coerce")
        invalid = sph.notna() & ((sph < 0) | (sph > 1))
        df.loc[invalid, "Std_Sphericity"] = np.nan
        df.loc[invalid, "Std_Sphericity_Source"] = "Invalid source value (outside 0-1)"

        # Safe fallback: conventional sphericity may be calculated only from
        # volume plus an explicitly named calibrated surface-area column.
        # Never derive it from the ambiguous generic legacy "Surface" field.
        if sphericity_surface_area and sphericity_surface_area in df.columns:
            missing = pd.to_numeric(df["Std_Sphericity"], errors="coerce").isna()
            volume_values = pd.to_numeric(df["Std_Volume_3D"], errors="coerce")
            surface_area_values = pd.to_numeric(df[sphericity_surface_area], errors="coerce")
            valid = missing & (volume_values > 0) & (surface_area_values > 0)
            calculated = (
                np.pi ** (1.0 / 3.0)
                * (6.0 * volume_values.loc[valid]) ** (2.0 / 3.0)
                / surface_area_values.loc[valid]
            )
            calculated = calculated[(calculated >= 0) & (calculated <= 1)]
            df.loc[calculated.index, "Std_Sphericity"] = calculated
            df.loc[calculated.index, "Std_Sphericity_Source"] = (
                "Calculated from Std_Volume_3D + " + str(sphericity_surface_area)
            )

        still_missing = pd.to_numeric(df["Std_Sphericity"], errors="coerce").isna()
        blank_source = df["Std_Sphericity_Source"].astype(str).str.strip().eq("")
        df.loc[still_missing & blank_source, "Std_Sphericity_Source"] = (
            "Unavailable: no Morph_Sphericity or explicit calibrated surface-area column"
        )

    df["Std_Shape_Compactness"] = np.where(
        df["Mode"].astype(str).eq("3D"),
        pd.to_numeric(df["Std_Sphericity"], errors="coerce"),
        pd.to_numeric(df["Std_Circularity_2D"], errors="coerce"),
    )

    copy_std("Std_X", x, numeric=True)
    copy_std("Std_Y", y, numeric=True)
    copy_std("Std_Z", z, numeric=True)

    # Remember the coordinate space used by the source columns.
    # Fiji 2D outputs normally provide X_um/Y_um, whereas 3D Objects Counter
    # commonly exports centroid X/Y in pixels. The entropy boundary must use
    # the same coordinate space as the satellite centroids.
    x_name = str(x) if x is not None else ""
    y_name = str(y) if y is not None else ""
    if ("um" in x_name.lower() or "µm" in x_name.lower()) and ("um" in y_name.lower() or "µm" in y_name.lower()):
        df["Std_XY_Coordinate_Space"] = "um"
    elif x_name in ["X", "XM"] or y_name in ["Y", "YM"]:
        df["Std_XY_Coordinate_Space"] = "pixel"
    else:
        df["Std_XY_Coordinate_Space"] = "um"

    # Keep an explicitly calibrated XYZ representation for true-3D spatial metrics.
    # 3D Objects Counter+ commonly reports X/Y as pixels and Z as a slice index,
    # whereas 2D SatelliteQ output normally already contains calibrated µm values.
    raw_x = pd.to_numeric(df["Std_X"], errors="coerce")
    raw_y = pd.to_numeric(df["Std_Y"], errors="coerce")
    raw_z = pd.to_numeric(df["Std_Z"], errors="coerce")
    pw = pd.to_numeric(df.get("Std_Pixel_Width_um", np.nan), errors="coerce")
    ph = pd.to_numeric(df.get("Std_Pixel_Height_um", np.nan), errors="coerce")
    pz = pd.to_numeric(df.get("Std_Voxel_Depth_um", np.nan), errors="coerce")
    if not isinstance(pw, pd.Series): pw = pd.Series(pw, index=df.index)
    if not isinstance(ph, pd.Series): ph = pd.Series(ph, index=df.index)
    if not isinstance(pz, pd.Series): pz = pd.Series(pz, index=df.index)

    if df["Std_XY_Coordinate_Space"].astype(str).eq("pixel").any():
        df["Std_X_um_3D"] = raw_x * pw
        df["Std_Y_um_3D"] = raw_y * ph
    else:
        df["Std_X_um_3D"] = raw_x
        df["Std_Y_um_3D"] = raw_y

    z_name = str(z) if z is not None else ""
    if "um" in z_name.lower() or "µm" in z_name.lower():
        df["Std_Z_um_3D"] = raw_z
    elif mode == "3D":
        # Fiji's current 3D SatelliteQ macro treats the 3D Objects Counter+ Z
        # coordinate as a 1-based slice index when calculating distances.
        df["Std_Z_um_3D"] = (raw_z - 1.0) * pz
    else:
        df["Std_Z_um_3D"] = raw_z
    copy_std("Std_Centrosome_X", cent_x, numeric=True)
    copy_std("Std_Centrosome_Y", cent_y, numeric=True)
    copy_std("Std_Centrosome_Z", cent_z, numeric=True)

    # Robust 2D distance fallback. If a Fiji table still contains calibrated
    # satellite and centrosome coordinates but its distance header was renamed,
    # reconstruct the missing XY distance so distance-based metrics remain usable.
    # In 3D, Fiji's calibrated 3D distance remains the source of truth.
    _dist = pd.to_numeric(df["Std_Distance"], errors="coerce")
    _mode_series = df["Mode"].astype(str).str.strip() if "Mode" in df.columns else pd.Series(mode, index=df.index)
    _sx = pd.to_numeric(df["Std_X"], errors="coerce")
    _sy = pd.to_numeric(df["Std_Y"], errors="coerce")
    _cx = pd.to_numeric(df["Std_Centrosome_X"], errors="coerce")
    _cy = pd.to_numeric(df["Std_Centrosome_Y"], errors="coerce")
    _fallback_mask = _dist.isna() & _mode_series.eq("2D") & _sx.notna() & _sy.notna() & _cx.notna() & _cy.notna()
    if _fallback_mask.any():
        df.loc[_fallback_mask, "Std_Distance"] = np.sqrt(
            (_sx.loc[_fallback_mask] - _cx.loc[_fallback_mask]) ** 2 +
            (_sy.loc[_fallback_mask] - _cy.loc[_fallback_mask]) ** 2
        )

    # Recover cell area and cell volume independently too.
    df["Std_Cell_Area_2D"] = _numeric_dimension_series(df, cell_area_2d_candidates)
    df["Std_Cell_Volume_3D"] = _numeric_dimension_series(df, cell_volume_3d_candidates)

    legacy_cell_measure = (
        pd.to_numeric(df[cell_measure], errors="coerce")
        if cell_measure and cell_measure in df.columns
        else pd.Series(np.nan, index=df.index)
    )
    if "Mode" in df.columns:
        mode_text = df["Mode"].astype(str).str.strip()
    else:
        mode_text = pd.Series(mode, index=df.index)
    df["Std_Cell_Measure"] = np.where(
        mode_text.eq("3D"),
        df["Std_Cell_Volume_3D"].combine_first(legacy_cell_measure),
        df["Std_Cell_Area_2D"].combine_first(legacy_cell_measure),
    )
    df["Std_Cell_Measure"] = pd.to_numeric(df["Std_Cell_Measure"], errors="coerce")

    df["Std_Cell_Measure_Name"] = cell_measure_name if pd.to_numeric(df["Std_Cell_Measure"], errors="coerce").notna().any() else ""

    df = add_cell_spatial_entropy_metrics(df, source_path)
    df = add_cell_pericentrosomal_clustering_index(df, source_path)
    df = add_center_distance_metrics(df, source_path)
    return df, mode, analysis_type

def _standardized_if_usable(raw: pd.DataFrame, path: Path) -> Optional[pd.DataFrame]:
    """Return standardized rows only when the table looks like a SatelliteQ quantification table."""
    if raw is None or raw.empty:
        return None
    try:
        df, _, _ = standardize_rows(raw, path)
    except Exception:
        return None
    if df is None or df.empty:
        return None

    # Require at least one SatelliteQ-like measurement/ID column after standardization.
    useful_cols = [
        "Satellite_ID", "Std_Size", "Std_Intensity", "Std_X", "Std_Y",
        "Volume", "Volume_um3", "Area", "Satellite_Area_um2",
        "Integrated_Intensity_AU", "Mean_Intensity_AU",
    ]
    if any(c in df.columns for c in useful_cols):
        return df
    return df if len(df.columns) > 5 else None


def read_input(path: Path) -> Optional[pd.DataFrame]:
    try:
        suffix = path.suffix.lower()

        if suffix in [".csv", ".tsv", ".txt"]:
            raw = read_csv_robust(path)
            return _standardized_if_usable(raw, path)

        if suffix in [".xlsx", ".xlsm", ".xls"]:
            # First try the official SatelliteQ report sheet names.
            preferred_sheets = [
                "Satellite_rows_clean",
                "Combined_satellite_rows",
                "Sheet1",
            ]

            tried = set()

            for sheet in preferred_sheets:
                try:
                    raw = pd.read_excel(path, sheet_name=sheet)
                    tried.add(sheet)
                    df = _standardized_if_usable(raw, path)
                    if df is not None and not df.empty:
                        return df
                except Exception:
                    pass

            # If the user converted a CSV/TSV to a normal Excel file, the sheet
            # may have any name. Try every sheet and use the first compatible one.
            try:
                xls = pd.ExcelFile(path)
                for sheet in xls.sheet_names:
                    if sheet in tried:
                        continue
                    try:
                        raw = pd.read_excel(path, sheet_name=sheet)
                        df = _standardized_if_usable(raw, path)
                        if df is not None and not df.empty:
                            return df
                    except Exception:
                        continue
            except Exception:
                return None

    except Exception as e:
        print(f"Skipped {path.name}: {e}")
        return None
    return None

def numeric(df: pd.DataFrame, col: str) -> pd.Series:
    if col not in df.columns:
        return pd.Series(dtype=float)
    return pd.to_numeric(df[col], errors="coerce")


def _read_cilium_csv(path: Path) -> pd.DataFrame:
    try:
        return read_csv_robust(Path(path)).dropna(how="all")
    except Exception:
        return pd.DataFrame()


def _cilium_companion_candidates(source_path: Path, kind: str) -> List[Path]:
    parent=Path(source_path).parent
    patterns = {
        "cells": ["*CiliumCellSummary.csv", "*Cilium_Cell_Summary.csv", "*cilium*cell*summary*.csv"],
        "rows": ["*CiliumQuantification.csv", "*Cilium_Quantification.csv", "*cilium*quantification*.csv"],
        "skeleton": ["*CiliumSkeletonPoints.csv", "*cilium*skeleton*points*.csv"],
    }[kind]
    found=[]
    for pat in patterns:
        found.extend(parent.glob(pat))
    return sorted(set(p for p in found if p.is_file()))


def _best_cilium_companion(source_path: Path, df: pd.DataFrame, kind: str) -> Optional[Path]:
    files=_cilium_companion_candidates(source_path,kind)
    if not files: return None
    stem=Path(source_path).stem.lower()
    prefix=stem
    for token in ["_satellitequantify_fiji_output","_timelapse_satellite_rows_combined","_results"]:
        if token in prefix: prefix=prefix.split(token)[0]
    direct=[p for p in files if prefix and prefix in p.stem.lower()]
    if len(direct)==1: return direct[0]
    # Match the Fiji Image field if several companion files share a folder.
    wanted=set(df.get("Image",pd.Series(dtype=object)).dropna().astype(str).str.strip())
    if wanted:
        scored=[]
        for p in files:
            t=_read_cilium_csv(p)
            score=0
            if "Image" in t.columns:
                score=len(wanted.intersection(set(t["Image"].dropna().astype(str).str.strip())))
            scored.append((score,p))
        scored.sort(key=lambda z:z[0],reverse=True)
        if scored and scored[0][0]>0: return scored[0][1]
    return files[0] if len(files)==1 else None


def load_cilium_companions_for_source(source_path: Path, df: pd.DataFrame) -> Tuple[pd.DataFrame,pd.DataFrame,pd.DataFrame]:
    """Return (cilium rows, per-cell cilium summary, skeleton path points)."""
    row_path=_best_cilium_companion(source_path,df,"rows")
    cell_path=_best_cilium_companion(source_path,df,"cells")
    skel_path=_best_cilium_companion(source_path,df,"skeleton")
    rows=_read_cilium_csv(row_path) if row_path else pd.DataFrame()
    cells=_read_cilium_csv(cell_path) if cell_path else pd.DataFrame()
    skel=_read_cilium_csv(skel_path) if skel_path else pd.DataFrame()
    return rows,cells,skel


def _standardize_cilium_cells(cells: pd.DataFrame) -> pd.DataFrame:
    if cells is None or cells.empty: return pd.DataFrame()
    c=cells.copy()
    aliases={
        "Cell_ID":["Cell_ID","Std_Cell_ID"], "Frame":["Frame","Std_Frame"], "Image":["Image"],
        "Cilium_Count":["Cilium_Count","Cilia_Count"],
        "Cilium_Length_um":["Longest_Cilium_Length_um","Cilium_Length_um","Max_Cilium_Length_um"],
        "Mean_Cilium_Length_um":["Mean_Cilium_Length_um"],
        "Ciliated":["Ciliated","Is_Ciliated"],
        "Ciliation_Length_Threshold_um":["Ciliation_Length_Threshold_um","Ciliated_Threshold_um"],
    }
    out=pd.DataFrame(index=c.index)
    for dst,cands in aliases.items():
        src=first_existing(list(c.columns),cands)
        out[dst]=c[src] if src else np.nan
    for col in ["Cilium_Count","Cilium_Length_um","Mean_Cilium_Length_um","Ciliated","Ciliation_Length_Threshold_um","Frame"]:
        out[col]=pd.to_numeric(out[col],errors="coerce")
    # Cells without a cilium should not carry a 0 µm cilium length into graphs/tables.
    # Keep ciliation status/count information, but blank cilium-length metrics unless the cell is ciliated.
    for col in ["Cilium_Length_um", "Mean_Cilium_Length_um"]:
        out.loc[pd.to_numeric(out[col], errors="coerce") <= 0, col] = np.nan
    non_ciliated_mask = out["Ciliated"].ne(1)
    count_zero_mask = pd.to_numeric(out["Cilium_Count"], errors="coerce").fillna(0).le(0)
    out.loc[non_ciliated_mask | count_zero_mask, ["Cilium_Length_um", "Mean_Cilium_Length_um"]] = np.nan
    out["Cell_ID_key"]=out["Cell_ID"].map(_normalize_id_text)
    out["Frame_key"]=out["Frame"].map(_normalize_id_text)
    out["Image_key"]=out["Image"].fillna("").astype(str).str.strip()
    out["Ciliation_status"]=np.where(out["Ciliated"].eq(1),"Ciliated",np.where(out["Ciliated"].eq(0),"Non-ciliated","Unknown"))
    return out


def attach_cilium_measurements(df: pd.DataFrame, source_path: Path) -> pd.DataFrame:
    """Attach optional Fiji cilium companion results to every satellite row."""
    out=df.copy()
    for col in ["Std_Cilium_Count","Std_Cilium_Length_um","Std_Cilium_Mean_Length_um","Std_Ciliated","Std_Ciliation_Length_Threshold_um"]:
        if col not in out.columns: out[col]=np.nan
    if "Ciliation_status" not in out.columns: out["Ciliation_status"]=""
    rows,cells_raw,skel=load_cilium_companions_for_source(source_path,out)
    cells=_standardize_cilium_cells(cells_raw)
    if cells.empty or "Std_Cell_ID" not in out.columns: return out

    for idx,row in out.iterrows():
        ck=_normalize_id_text(row.get("Std_Cell_ID",np.nan))
        if not ck: continue
        candidates=cells[cells["Cell_ID_key"].eq(ck)]
        fk=_normalize_id_text(row.get("Std_Frame",np.nan))
        if fk and candidates["Frame_key"].ne("").any():
            same=candidates[candidates["Frame_key"].eq(fk)]
            if not same.empty: candidates=same
        img=str(row.get("Image","")).strip()
        if img and candidates["Image_key"].ne("").any():
            same=candidates[candidates["Image_key"].eq(img)]
            if not same.empty: candidates=same
        if candidates.empty: continue
        cr=candidates.iloc[0]
        mapping={"Std_Cilium_Count":"Cilium_Count","Std_Cilium_Length_um":"Cilium_Length_um","Std_Cilium_Mean_Length_um":"Mean_Cilium_Length_um","Std_Ciliated":"Ciliated","Std_Ciliation_Length_Threshold_um":"Ciliation_Length_Threshold_um"}
        for dst,src_col in mapping.items(): out.at[idx,dst]=cr.get(src_col,np.nan)
        out.at[idx,"Ciliation_status"]=cr.get("Ciliation_status","")
    return out


def collect_cilium_companions(df: pd.DataFrame) -> Tuple[pd.DataFrame,pd.DataFrame,pd.DataFrame]:
    """Collect unique cilium rows/cell summaries/skeleton points for report export/viewers."""
    rows_all=[]; cells_all=[]; skel_all=[]; seen=set()
    if "Source_Path" not in df.columns: return pd.DataFrame(),pd.DataFrame(),pd.DataFrame()
    for sp in df["Source_Path"].dropna().astype(str).unique():
        p=Path(sp)
        subset=df[df["Source_Path"].astype(str).eq(sp)]
        rows,cells,skel=load_cilium_companions_for_source(p,subset)
        key=str(p.resolve()) if p.exists() else str(p)
        if key in seen: continue
        seen.add(key)
        meta={
            "Source_UID": subset["Source_UID"].iloc[0] if "Source_UID" in subset.columns and not subset.empty else key,
            "Source_Input": subset["Source_Input"].iloc[0] if "Source_Input" in subset.columns and not subset.empty else p.name,
            "Condition": subset["Condition"].iloc[0] if "Condition" in subset.columns and not subset.empty else "",
            "Replicate": subset["Replicate"].iloc[0] if "Replicate" in subset.columns and not subset.empty else "",
        }
        for table,bucket in [(rows,rows_all),(cells,cells_all),(skel,skel_all)]:
            if table is not None and not table.empty:
                t=table.copy()
                for k,v in meta.items():
                    if k not in t.columns: t[k]=v
                bucket.append(t)
    cat=lambda xs: pd.concat(xs,ignore_index=True,sort=False).drop_duplicates() if xs else pd.DataFrame()
    return cat(rows_all),cat(cells_all),cat(skel_all)


def sanitize_cilium_length_metrics(table: pd.DataFrame) -> pd.DataFrame:
    """Treat absence of a cilium as missing length, never as a 0-µm cilium.

    This is intentionally applied late in the reporting pipeline as well as during
    cilium import, because outer merges with Fiji's per-cell cilium companion table
    can otherwise reintroduce zero placeholders. Ciliation/count fields are retained.
    """
    if table is None:
        return table
    out = table.copy()
    if out.empty:
        return out

    length_cols = [c for c in [
        "Cilium_length_um", "Mean_cilium_length_um",
        "Std_Cilium_Length_um", "Std_Cilium_Mean_Length_um",
        "Cilium_Length_um", "Mean_Cilium_Length_um",
    ] if c in out.columns]
    if not length_cols:
        return out

    non_ciliated = pd.Series(False, index=out.index)

    # Explicit ciliated yes/no columns.
    for c in ["Ciliated", "Std_Ciliated"]:
        if c in out.columns:
            cv = pd.to_numeric(out[c], errors="coerce")
            non_ciliated |= cv.eq(0)

    # Text status is useful when numeric ciliation is absent.
    if "Ciliation_status" in out.columns:
        st = out["Ciliation_status"].fillna("").astype(str).str.strip().str.lower()
        non_ciliated |= st.isin(["non-ciliated", "non ciliated", "nonciliated", "no", "false"])

    # A confirmed count of zero also means there is no cilium length to plot.
    for c in ["Cilium_count", "Cilium_Count", "Std_Cilium_Count"]:
        if c in out.columns:
            cnt = pd.to_numeric(out[c], errors="coerce")
            non_ciliated |= cnt.eq(0)

    for c in length_cols:
        vals = pd.to_numeric(out[c], errors="coerce")
        # Zero/negative placeholders are never valid biological cilium lengths.
        vals = vals.mask(vals <= 0, np.nan)
        vals = vals.mask(non_ciliated, np.nan)
        out[c] = vals

    return out


def merge_cilium_only_cells_into_summary(cell_summary: pd.DataFrame, cilium_cells_raw: pd.DataFrame) -> pd.DataFrame:
    """Add segmented cells that have cilium data but zero detected satellites.

    Fiji writes one row for every segmented cell to CiliumCellSummary, including
    cells in which no accepted cilium was found. A satellite-derived per-cell
    table normally cannot contain a cell with zero satellites. This outer merge
    keeps those cells in ciliation comparisons and assigns Satellite_count=0.
    """
    if cilium_cells_raw is None or cilium_cells_raw.empty:
        return cell_summary
    base = cell_summary.copy() if cell_summary is not None else pd.DataFrame()
    raw = cilium_cells_raw.copy()
    std = _standardize_cilium_cells(raw)
    if std.empty:
        return base
    # Reattach report metadata that _standardize_cilium_cells intentionally ignores.
    for col in ["Source_UID","Source_Input","Source_Path","Source_Folder","Condition","Replicate","Mode"]:
        if col in raw.columns:
            std[col] = raw[col].values
        elif col not in std.columns:
            std[col] = ""

    if base.empty:
        base = pd.DataFrame(columns=[
            "Source_UID","Source_Input","Source_Path","Source_Folder","Condition","Replicate",
            "Unique_Cell","Cell_ID","Mode","Satellite_count","Cilium_length_um","Ciliated",
            "Ciliation_status","Ciliation_length_threshold_um"
        ])
    # Ensure cilium result columns exist in the main cell table.
    required = {
        "Cilium_length_um": np.nan,
        "Ciliated": np.nan,
        "Ciliation_status": "",
        "Ciliation_length_threshold_um": np.nan,
    }
    for col, default in required.items():
        if col not in base.columns:
            base[col] = default

    def key_from_row(row):
        uid=str(row.get("Source_UID","")).strip()
        cid=_normalize_id_text(row.get("Cell_ID",np.nan))
        frame=_normalize_id_text(row.get("Frame",np.nan))
        return uid,cid,frame

    existing={}
    for idx,row in base.iterrows():
        existing[key_from_row(row)] = idx

    additions=[]
    for _,r in std.iterrows():
        cid=_normalize_id_text(r.get("Cell_ID",np.nan))
        if not cid:
            continue
        uid=str(r.get("Source_UID","")).strip()
        frame=_normalize_id_text(r.get("Frame",np.nan))
        key=(uid,cid,frame)
        # When Frame is absent, match on source and cell.
        match_idx=existing.get(key)
        if match_idx is None and frame:
            match_idx=existing.get((uid,cid,""))
        vals={
            "Cilium_length_um": r.get("Cilium_Length_um",np.nan),
            "Ciliated": r.get("Ciliated",np.nan),
            "Ciliation_status": r.get("Ciliation_status",""),
            "Ciliation_length_threshold_um": r.get("Ciliation_Length_Threshold_um",np.nan),
        }
        if match_idx is not None:
            for col,val in vals.items():
                if pd.notna(val) and (col != "Ciliation_status" or str(val).strip()):
                    base.at[match_idx,col]=val
            continue

        row={c:np.nan for c in base.columns}
        for c in ["Source_UID","Source_Input","Source_Path","Source_Folder","Condition","Replicate","Mode"]:
            if c in row:
                row[c]=r.get(c,"")
        row["Cell_ID"]=r.get("Cell_ID",cid)
        source_label=str(r.get("Source_Input","")).strip() or uid
        row["Unique_Cell"]=source_label + " | Cell " + cid
        row["Satellite_count"]=0
        row.update(vals)
        if "Frame" in base.columns:
            row["Frame"]=r.get("Frame",np.nan)
        additions.append(row)
    if additions:
        base=pd.concat([base,pd.DataFrame(additions)],ignore_index=True,sort=False)
    return sanitize_cilium_length_metrics(base)


def update_image_ciliation_from_cell_summary(image_summary: pd.DataFrame, cell_summary: pd.DataFrame) -> pd.DataFrame:
    """Recalculate ciliation fields using all cell rows, including zero-satellite cells."""
    if image_summary is None or image_summary.empty or cell_summary is None or cell_summary.empty:
        return image_summary
    out=image_summary.copy()
    if "Source_UID" not in out.columns or "Source_UID" not in cell_summary.columns:
        return out
    for uid,g in cell_summary.groupby("Source_UID",dropna=False):
        idx=out.index[out["Source_UID"].astype(str).eq(str(uid))]
        if len(idx)==0: continue
        ci=pd.to_numeric(g.get("Ciliated",pd.Series(dtype=float)),errors="coerce").dropna()
        lens=pd.to_numeric(g.get("Cilium_length_um",pd.Series(dtype=float)),errors="coerce").dropna()
        if not ci.empty:
            out.loc[idx,"Ciliated_cell_count"]=float((ci==1).sum())
            out.loc[idx,"Cell_count_with_ciliation_data"]=float(len(ci))
            out.loc[idx,"Ciliation_fraction"]=float(ci.mean())
        if not lens.empty:
            out.loc[idx,"Mean_cilium_length_um"]=float(lens.mean())
    return out

def _normalize_id_text(value) -> str:
    """Normalize Cell_ID / Frame-like values so 1 and 1.0 match."""
    if pd.isna(value):
        return ""
    s = str(value).strip()
    if s.lower() in ["", "nan", "none", "na"]:
        return ""
    try:
        f = float(s)
        if np.isfinite(f) and float(f).is_integer():
            return str(int(f))
    except Exception:
        pass
    return s


def _source_run_id_candidates(source_path: Path, df: Optional[pd.DataFrame] = None) -> List[str]:
    """Return plausible Fiji Run_ID values for one quantification source."""
    found: List[str] = []

    def add(value):
        if value is None:
            return
        s = str(value).strip()
        if not s or s.lower() in {"nan", "none", "na"}:
            return
        if s not in found:
            found.append(s)

    if df is not None and not df.empty:
        for col in ["Run_ID", "Run ID"]:
            if col in df.columns:
                for value in df[col].dropna().astype(str).unique():
                    add(value)

    p = Path(source_path)
    names = [p.stem, p.name, p.parent.name]
    if df is not None and not df.empty:
        for col in ["Source_Input", "Source_File"]:
            if col in df.columns:
                names.extend(df[col].dropna().astype(str).unique().tolist())

    for raw_name in names:
        name = Path(str(raw_name)).stem.strip()
        if not name:
            continue

        # Run-specific Fiji output folders.
        for suffix in [
            "_SatelliteQuantify_3D_satellites_only",
            "_SatelliteQuantify_results",
            "_SQ3D",
            "_SQ2D",
        ]:
            if name.endswith(suffix):
                add(name[:-len(suffix)])

        # 2D single-image output: <Run_ID>_SatelliteQuantify_Fiji_Output.csv
        suffix = "_SatelliteQuantify_Fiji_Output"
        if name.endswith(suffix):
            add(name[:-len(suffix)])

        # Combined time-lapse output generally keeps the base/run ID.
        suffix = "_timelapse_satellite_rows_combined"
        if name.endswith(suffix):
            add(name[:-len(suffix)])

        # 3D single/TL frame output:
        #   <base>_3D_YYYY_M_D_H_M_S_results.csv
        # while the boundary file is:
        #   <base>_YYYY_M_D_H_M_S_cell_boundaries_for_python.csv
        m = re.match(
            r"^(?P<base>.+)_3D_(?P<stamp>\d{4}_\d{1,2}_\d{1,2}_\d{1,2}_\d{1,2}_\d{1,2})_results(?:_excel_friendly)?$",
            name,
            flags=re.IGNORECASE,
        )
        if m:
            add(f"{m.group('base')}_{m.group('stamp')}")

        # Already looks like a Fiji Run_ID ending in a timestamp.
        if re.search(r"_\d{4}_\d{1,2}_\d{1,2}_\d{1,2}_\d{1,2}_\d{1,2}$", name):
            add(name)

    return found


def _load_boundary_polygons_for_entropy(source_path: Path, df: pd.DataFrame) -> dict:
    """Load source-specific Fiji cell-boundary polygons using Run_ID, image, frame and cell identity."""
    parent = Path(source_path).parent
    all_files = sorted(parent.glob("*cell_boundaries_for_python.csv"))
    if not all_files:
        return {}

    wanted_runs = set(_source_run_id_candidates(source_path, df))
    wanted_images = set()
    if "Image" in df.columns:
        wanted_images = {
            str(v).strip() for v in df["Image"].dropna().astype(str).unique()
            if str(v).strip() and str(v).strip().lower() not in {"nan", "none", "na"}
        }

    exact_files = []
    for run_id in wanted_runs:
        candidate = parent / f"{run_id}_cell_boundaries_for_python.csv"
        if candidate.exists() and candidate not in exact_files:
            exact_files.append(candidate)

    parsed = []
    for f in all_files:
        try:
            b = pd.read_csv(f)
        except Exception:
            continue
        if b.empty or "Cell_ID" not in b.columns:
            continue
        file_run = f.name[:-len("_cell_boundaries_for_python.csv")] if f.name.endswith("_cell_boundaries_for_python.csv") else f.stem
        run_values = set()
        if "Run_ID" in b.columns:
            run_values = {
                str(v).strip() for v in b["Run_ID"].dropna().astype(str).unique()
                if str(v).strip() and str(v).strip().lower() not in {"nan", "none", "na"}
            }
        image_values = set()
        if "Image" in b.columns:
            image_values = {
                str(v).strip() for v in b["Image"].dropna().astype(str).unique()
                if str(v).strip() and str(v).strip().lower() not in {"nan", "none", "na"}
            }
        score = 0
        if f in exact_files:
            score += 1000
        if wanted_runs and (file_run in wanted_runs or bool(run_values & wanted_runs)):
            score += 500
        if wanted_images and bool(image_values & wanted_images):
            score += 100
        parsed.append((f, b, file_run, run_values, image_values, score))

    if not parsed:
        return {}

    if exact_files:
        selected = [item for item in parsed if item[0] in exact_files]
    else:
        best_score = max(item[5] for item in parsed)
        if best_score > 0:
            selected = [item for item in parsed if item[5] == best_score]
        elif len(parsed) == 1:
            selected = parsed
        else:
            return {}

    polygons = {}
    for f, b, file_run, run_values, image_values, _score in selected:
        coordinate_pairs = []
        if "X_um" in b.columns and "Y_um" in b.columns:
            coordinate_pairs.append(("um", "X_um", "Y_um"))
        if "X_pixel" in b.columns and "Y_pixel" in b.columns:
            coordinate_pairs.append(("pixel", "X_pixel", "Y_pixel"))
        if not coordinate_pairs:
            continue

        # Once a source was selected by Run_ID, retain only rows from that run
        # if the CSV itself happens to contain more than one run.
        if wanted_runs and "Run_ID" in b.columns:
            run_mask = b["Run_ID"].astype(str).str.strip().isin(wanted_runs)
            if run_mask.any():
                b = b.loc[run_mask].copy()

        if wanted_images and "Image" in b.columns:
            image_mask = b["Image"].astype(str).str.strip().isin(wanted_images)
            if image_mask.any():
                b = b.loc[image_mask].copy()

        for space, x_col, y_col in coordinate_pairs:
            work = pd.DataFrame({
                "_run": b["Run_ID"].astype(str) if "Run_ID" in b.columns else file_run,
                "_image": b["Image"].astype(str) if "Image" in b.columns else "",
                "_frame": b["Frame"] if "Frame" in b.columns else "",
                "_cell": b["Cell_ID"],
                "_x": pd.to_numeric(b[x_col], errors="coerce"),
                "_y": pd.to_numeric(b[y_col], errors="coerce"),
                "_point": pd.to_numeric(b["Point_Index"], errors="coerce") if "Point_Index" in b.columns else np.arange(len(b)),
            }).dropna(subset=["_x", "_y"])

            if work.empty:
                continue

            work["_run_key"] = work["_run"].astype(str).str.strip()
            work["_image_key"] = work["_image"].astype(str).str.strip()
            work["_frame_key"] = work["_frame"].map(_normalize_id_text)
            work["_cell_key"] = work["_cell"].map(_normalize_id_text)

            for key, g in work.groupby(["_run_key", "_image_key", "_frame_key", "_cell_key"], dropna=False, sort=False):
                if len(g) < 3:
                    continue
                g = g.sort_values("_point", kind="stable")
                coords = g[["_x", "_y"]].to_numpy(dtype=float)
                if np.isfinite(coords).all(axis=1).sum() >= 3:
                    polygons[(space,) + key] = coords

    return polygons


def _choose_boundary_polygon(boundaries: dict, group: pd.DataFrame):
    if not boundaries:
        return None

    run_keys = []
    if "Run_ID" in group.columns:
        for value in group["Run_ID"].dropna().astype(str).unique():
            value = str(value).strip()
            if value and value.lower() not in {"nan", "none", "na"} and value not in run_keys:
                run_keys.append(value)

    if not run_keys:
        source_value = None
        for col in ["Source_Path", "Source_Input", "Source_File"]:
            if col in group.columns and group[col].dropna().size:
                source_value = str(group[col].dropna().iloc[0])
                break
        if source_value:
            run_keys.extend(_source_run_id_candidates(Path(source_value), group))

    image_key = ""
    if "Image" in group.columns:
        vals = group["Image"].dropna().astype(str)
        if not vals.empty:
            image_key = str(vals.iloc[0]).strip()

    frame_key = ""
    if "Std_Frame" in group.columns:
        vals = pd.to_numeric(group["Std_Frame"], errors="coerce").dropna()
        if not vals.empty:
            frame_key = _normalize_id_text(vals.iloc[0])

    cell_key = ""
    if "Std_Cell_ID" in group.columns:
        vals = group["Std_Cell_ID"].dropna()
        if not vals.empty:
            cell_key = _normalize_id_text(vals.iloc[0])

    space = "um"
    if "Std_XY_Coordinate_Space" in group.columns:
        vals = group["Std_XY_Coordinate_Space"].dropna().astype(str).str.strip().str.lower()
        if not vals.empty and vals.iloc[0] in ["um", "pixel"]:
            space = vals.iloc[0]

    for run_key in run_keys or [""]:
        exact = boundaries.get((space, run_key, image_key, frame_key, cell_key))
        if exact is not None:
            return exact

    candidates = []
    for (poly_space, run_id, img, fr, cell), poly in boundaries.items():
        if poly_space != space or cell != cell_key:
            continue
        if run_keys and run_id and run_id not in run_keys:
            continue
        if frame_key and fr and fr != frame_key:
            continue
        if image_key and img and img != image_key:
            continue
        candidates.append(poly)
    if len(candidates) == 1:
        return candidates[0]

    fallback = []
    for (_poly_space, run_id, img, fr, cell), poly in boundaries.items():
        if cell != cell_key:
            continue
        if run_keys and run_id and run_id not in run_keys:
            continue
        if frame_key and fr and fr != frame_key:
            continue
        if image_key and img and img != image_key:
            continue
        fallback.append(poly)
    return fallback[0] if len(fallback) == 1 else None

def _normalized_shannon_entropy_from_grid(
    x: np.ndarray,
    y: np.ndarray,
    weights=None,
    polygon=None,
    cell_area_um2=None,
    grid_size: int = 10,
) -> float:
    """Return normalized spatial Shannon entropy for one cell at one grid scale.

    Satellite centroid positions define where signal is located. ``weights``
    define how much satellite material each centroid contributes. For the main
    SatelliteQ metric, weights are per-satellite integrated intensities.

    0 = almost all satellite signal concentrated in one spatial region.
    1 = the satellite signal is maximally evenly distributed across the
        available regions for the observed number of detected satellites.

    The denominator uses log(min(K, N)), where K is the number of available
    cell regions and N is the number of detected satellites. This prevents a
    cell with fewer satellites than grid bins from being penalized simply
    because it cannot occupy every bin.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    if weights is None:
        weights = np.ones(len(x), dtype=float)
    else:
        weights = np.asarray(weights, dtype=float)

    valid = np.isfinite(x) & np.isfinite(y) & np.isfinite(weights) & (weights >= 0)
    x, y, weights = x[valid], y[valid], weights[valid]
    n = len(x)
    if n < 2:
        return np.nan

    # If all usable weights are zero, fall back to equal weights rather than
    # returning a misleading zero entropy.
    if float(np.sum(weights)) <= 0:
        weights = np.ones(n, dtype=float)

    g = max(2, int(grid_size))
    poly = None if polygon is None else np.asarray(polygon, dtype=float)

    if poly is not None:
        poly = poly[np.isfinite(poly).all(axis=1)]
        if len(poly) >= 3:
            xmin, ymin = np.min(poly[:, 0]), np.min(poly[:, 1])
            xmax, ymax = np.max(poly[:, 0]), np.max(poly[:, 1])
            if not (np.isfinite([xmin, xmax, ymin, ymax]).all() and xmax > xmin and ymax > ymin):
                poly = None
        else:
            poly = None

    if poly is None:
        try:
            area = float(cell_area_um2)
        except Exception:
            area = np.nan
        if not np.isfinite(area) or area <= 0:
            return np.nan
        radius = float(np.sqrt(area / np.pi))
        cx, cy = float(np.nanmean(x)), float(np.nanmean(y))
        xmin, xmax = cx - radius, cx + radius
        ymin, ymax = cy - radius, cy + radius
        domain = "circle"
    else:
        domain = "polygon"

    ix = np.floor((x - xmin) / (xmax - xmin) * g).astype(int)
    iy = np.floor((y - ymin) / (ymax - ymin) * g).astype(int)
    ix, iy = np.clip(ix, 0, g - 1), np.clip(iy, 0, g - 1)

    weighted_bins = np.zeros((g, g), dtype=float)
    np.add.at(weighted_bins, (iy, ix), weights)

    xc = xmin + (np.arange(g) + 0.5) * (xmax - xmin) / g
    yc = ymin + (np.arange(g) + 0.5) * (ymax - ymin) / g
    xx, yy = np.meshgrid(xc, yc)
    centers = np.column_stack([xx.ravel(), yy.ravel()])

    if domain == "polygon":
        try:
            from matplotlib.path import Path as MplPath
            active = MplPath(poly).contains_points(centers, radius=1e-12).reshape(g, g)
        except Exception:
            active = np.ones((g, g), dtype=bool)

        vx = np.floor((poly[:, 0] - xmin) / (xmax - xmin) * g).astype(int)
        vy = np.floor((poly[:, 1] - ymin) / (ymax - ymin) * g).astype(int)
        vx, vy = np.clip(vx, 0, g - 1), np.clip(vy, 0, g - 1)
        active[vy, vx] = True
    else:
        cx, cy = (xmin + xmax) / 2.0, (ymin + ymax) / 2.0
        radius = (xmax - xmin) / 2.0
        active = ((xx - cx) ** 2 + (yy - cy) ** 2) <= radius ** 2

    # Always include bins containing detected satellites, including edge bins.
    active[iy, ix] = True

    vals = weighted_bins[active]
    k_available = int(active.sum())
    total = float(vals.sum())
    if total <= 0 or k_available < 2:
        return np.nan

    p = vals[vals > 0] / total
    h = float(-np.sum(p * np.log(p)))

    # Granule-number-aware normalization.
    max_states = min(k_available, n)
    if max_states <= 1:
        return np.nan
    hmax = float(np.log(max_states))
    return float(np.clip(h / hmax, 0.0, 1.0)) if hmax > 0 else np.nan


def _multiscale_spatial_entropy(
    x: np.ndarray,
    y: np.ndarray,
    weights,
    polygon=None,
    cell_area_um2=None,
    grid_sizes=(5, 10, 20),
):
    """Average normalized entropy across coarse, medium and fine spatial scales."""
    values = []
    for grid_size in grid_sizes:
        value = _normalized_shannon_entropy_from_grid(
            x=x,
            y=y,
            weights=weights,
            polygon=polygon,
            cell_area_um2=cell_area_um2,
            grid_size=grid_size,
        )
        values.append(value)

    finite = [v for v in values if np.isfinite(v)]
    mean_value = float(np.mean(finite)) if finite else np.nan
    return mean_value, values


def add_cell_spatial_entropy_metrics(df: pd.DataFrame, source_path: Path) -> pd.DataFrame:
    """Add collective cell spatial entropy and repeat it on all satellite rows.

    Main metric:
      Std_Cell_Spatial_Entropy
        = mean of 5x5, 10x10 and 20x20 normalized Shannon entropies,
          weighted by each satellite's integrated intensity.

    QC metric:
      Std_Cell_Spatial_Entropy_Count_QC
        = the same multiscale calculation with equal weight per detected
          satellite. It is retained in Excel for troubleshooting segmentation
          but is not used as the primary graphing metric.

    Entropy is calculated separately for every cell and every time-lapse frame.
    """
    out = df.copy()
    out["Std_Cell_Spatial_Entropy"] = np.nan
    out["Std_Cell_Spatial_Entropy_5x5"] = np.nan
    out["Std_Cell_Spatial_Entropy_10x10"] = np.nan
    out["Std_Cell_Spatial_Entropy_20x20"] = np.nan
    out["Std_Cell_Spatial_Entropy_Count_QC"] = np.nan

    if out.empty or "Std_X" not in out.columns or "Std_Y" not in out.columns or "Std_Cell_ID" not in out.columns:
        return out

    boundaries = _load_boundary_polygons_for_entropy(Path(source_path), out)

    group_source = "Source_UID" if "Source_UID" in out.columns else "Source_Input"
    group_cols = [group_source]
    if "Std_Frame" in out.columns and pd.to_numeric(out["Std_Frame"], errors="coerce").notna().any():
        group_cols.append("Std_Frame")
    group_cols.append("Std_Cell_ID")

    cell_text = out["Std_Cell_ID"].astype(str).str.strip()
    valid_cell = out["Std_Cell_ID"].notna() & ~cell_text.str.lower().isin(["", "nan", "none", "na"])

    for _, group in out.loc[valid_cell].groupby(group_cols, dropna=False, sort=False):
        x = pd.to_numeric(group["Std_X"], errors="coerce").to_numpy(dtype=float)
        y = pd.to_numeric(group["Std_Y"], errors="coerce").to_numpy(dtype=float)
        valid_xy = np.isfinite(x) & np.isfinite(y)
        if valid_xy.sum() < 2:
            continue

        # Integrated intensity is the preferred weight because a large merged
        # aggregate should contribute according to how much satellite signal it
        # contains, not as only one centroid.
        if "Std_Intensity" in group.columns:
            intensity = pd.to_numeric(group["Std_Intensity"], errors="coerce").to_numpy(dtype=float)
        elif "Std_Total_Intensity" in group.columns:
            intensity = pd.to_numeric(group["Std_Total_Intensity"], errors="coerce").to_numpy(dtype=float)
        else:
            intensity = np.ones(len(group), dtype=float)

        # Replace missing/negative intensity values with zero. If the entire
        # cell has unusable intensities, the entropy helper falls back to equal
        # object weights so the result is still defined.
        intensity = np.where(np.isfinite(intensity) & (intensity >= 0), intensity, 0.0)

        polygon = _choose_boundary_polygon(boundaries, group)

        coordinate_space = "um"
        if "Std_XY_Coordinate_Space" in group.columns:
            vals = group["Std_XY_Coordinate_Space"].dropna().astype(str).str.strip().str.lower()
            if not vals.empty and vals.iloc[0] in ["um", "pixel"]:
                coordinate_space = vals.iloc[0]

        area = pd.to_numeric(
            group.get("Std_Cell_Area_2D", pd.Series(np.nan, index=group.index)),
            errors="coerce",
        )
        area = area[(area > 0) & np.isfinite(area)]
        cell_area = float(area.median()) if not area.empty else np.nan

        if coordinate_space == "um" and not np.isfinite(cell_area):
            volume = pd.to_numeric(
                group.get("Std_Cell_Volume_3D", pd.Series(np.nan, index=group.index)),
                errors="coerce",
            )
            volume = volume[(volume > 0) & np.isfinite(volume)]
            if not volume.empty:
                r = (3.0 * float(volume.median()) / (4.0 * np.pi)) ** (1.0 / 3.0)
                cell_area = np.pi * r ** 2

        if polygon is None and coordinate_space == "pixel":
            xv = x[np.isfinite(x)]
            yv = y[np.isfinite(y)]
            if len(xv) >= 2 and len(yv) >= 2:
                xmin, xmax = float(np.min(xv)), float(np.max(xv))
                ymin, ymax = float(np.min(yv)), float(np.max(yv))
                dx = max(xmax - xmin, 1.0)
                dy = max(ymax - ymin, 1.0)
                pad_x = max(dx * 0.10, 1.0)
                pad_y = max(dy * 0.10, 1.0)
                polygon = np.array([
                    [xmin - pad_x, ymin - pad_y],
                    [xmax + pad_x, ymin - pad_y],
                    [xmax + pad_x, ymax + pad_y],
                    [xmin - pad_x, ymax + pad_y],
                ], dtype=float)

        domain_area = cell_area if coordinate_space == "um" and np.isfinite(cell_area) else None

        weighted_entropy, scale_values = _multiscale_spatial_entropy(
            x=x,
            y=y,
            weights=intensity,
            polygon=polygon,
            cell_area_um2=domain_area,
            grid_sizes=(5, 10, 20),
        )

        count_entropy, _ = _multiscale_spatial_entropy(
            x=x,
            y=y,
            weights=np.ones(len(group), dtype=float),
            polygon=polygon,
            cell_area_um2=domain_area,
            grid_sizes=(5, 10, 20),
        )

        out.loc[group.index, "Std_Cell_Spatial_Entropy"] = weighted_entropy
        out.loc[group.index, "Std_Cell_Spatial_Entropy_5x5"] = scale_values[0]
        out.loc[group.index, "Std_Cell_Spatial_Entropy_10x10"] = scale_values[1]
        out.loc[group.index, "Std_Cell_Spatial_Entropy_20x20"] = scale_values[2]
        out.loc[group.index, "Std_Cell_Spatial_Entropy_Count_QC"] = count_entropy

    return out












def add_cell_pericentrosomal_clustering_index(df: pd.DataFrame, source_path: Path) -> pd.DataFrame:
    """Calculate the pericentrosomal clustering index and component intensities per cell."""
    out = df.copy()
    out[PCI_STD_INNER_COL] = np.nan
    out[PCI_STD_OUTER_COL] = np.nan
    out[PCI_STD_TOTAL_COL] = np.nan
    out["Std_Cell_Pericentrosomal_Clustering_Index"] = np.nan
    out["Std_Cell_Pericentrosomal_Clustering_Status"] = ""
    out["Std_Cell_Pericentrosomal_Geometry"] = ""

    required = {"Std_Cell_ID", "Std_Distance"}
    if out.empty or not required.issubset(out.columns):
        return out

    intensity_col = None
    for candidate in ["Std_Total_Intensity", "Std_Intensity"]:
        if candidate in out.columns and pd.to_numeric(out[candidate], errors="coerce").notna().any():
            intensity_col = candidate
            break
    if intensity_col is None:
        return out

    group_source = "Source_UID" if "Source_UID" in out.columns else "Source_Input"
    group_cols = [group_source]
    if "Std_Frame" in out.columns and pd.to_numeric(out["Std_Frame"], errors="coerce").notna().any():
        group_cols.append("Std_Frame")
    group_cols.append("Std_Cell_ID")

    cell_text = out["Std_Cell_ID"].astype(str).str.strip()
    valid_cell = out["Std_Cell_ID"].notna() & ~cell_text.str.lower().isin(["", "nan", "none", "na"])


    for _, group in out.loc[valid_cell].groupby(group_cols, dropna=False, sort=False):
        distances = pd.to_numeric(group["Std_Distance"], errors="coerce").to_numpy(dtype=float)
        intensities = pd.to_numeric(group[intensity_col], errors="coerce").to_numpy(dtype=float)
        valid = np.isfinite(distances) & (distances >= 0) & np.isfinite(intensities) & (intensities >= 0)
        if not np.any(valid):
            continue

        d = distances[valid]
        inten = intensities[valid]
        inner = float(np.sum(inten[d <= 3.0]))
        total_outer = float(np.sum(inten[d <= 12.0]))
        outer = float(total_outer - inner)
        if abs(outer) < 1e-12:
            outer = 0.0

        if outer > 0:
            ratio = float(inner / outer)
            status = "Calculated"
        elif inner > 0:
            ratio = np.nan
            status = "Outer 3-12 um fluorescence is zero; index undefined"
        else:
            ratio = np.nan
            status = "No fluorescence within 12 um"

        mode_vals = group.get("Mode", pd.Series(dtype=object)).dropna().astype(str).str.strip()
        geometry = "3D radial distance" if (not mode_vals.empty and mode_vals.iloc[0] == "3D") else "2D radial distance"

        out.loc[group.index, PCI_STD_INNER_COL] = inner
        out.loc[group.index, PCI_STD_OUTER_COL] = outer
        out.loc[group.index, PCI_STD_TOTAL_COL] = total_outer
        out.loc[group.index, "Std_Cell_Pericentrosomal_Clustering_Index"] = ratio
        out.loc[group.index, "Std_Cell_Pericentrosomal_Clustering_Status"] = status
        out.loc[group.index, "Std_Cell_Pericentrosomal_Geometry"] = geometry

    return out


def _polygon_centroid_xy(polygon) -> tuple:
    """Return geometric XY centroid of a polygon, or (nan, nan) if unavailable."""
    if polygon is None:
        return np.nan, np.nan
    p = np.asarray(polygon, dtype=float)
    p = p[np.isfinite(p).all(axis=1)]
    if len(p) < 3:
        return np.nan, np.nan
    x, y = p[:, 0], p[:, 1]
    x2 = np.r_[x, x[0]]
    y2 = np.r_[y, y[0]]
    cross = x2[:-1] * y2[1:] - x2[1:] * y2[:-1]
    area2 = float(np.sum(cross))
    if np.isfinite(area2) and abs(area2) > 1e-12:
        cx = float(np.sum((x2[:-1] + x2[1:]) * cross) / (3.0 * area2))
        cy = float(np.sum((y2[:-1] + y2[1:]) * cross) / (3.0 * area2))
        if np.isfinite(cx) and np.isfinite(cy):
            return cx, cy
    return float(np.mean(x)), float(np.mean(y))


def _median_positive(group: pd.DataFrame, column: str) -> float:
    if column not in group.columns:
        return np.nan
    s = pd.to_numeric(group[column], errors="coerce")
    s = s[np.isfinite(s) & (s > 0)]
    return float(s.median()) if not s.empty else np.nan


def _image_center_xy(group: pd.DataFrame, coordinate_space: str):
    """Return image center in the same coordinate system as Std_X/Std_Y."""
    width = _median_positive(group, "Std_Image_Width_Pixels")
    height = _median_positive(group, "Std_Image_Height_Pixels")
    pw = _median_positive(group, "Std_Pixel_Width_um")
    ph = _median_positive(group, "Std_Pixel_Height_um")
    if np.isfinite(width) and np.isfinite(height):
        if coordinate_space == "pixel":
            return width / 2.0, height / 2.0, "Image center (exact image dimensions)"
        if np.isfinite(pw) and np.isfinite(ph):
            return width * pw / 2.0, height * ph / 2.0, "Image center (exact image dimensions)"

    # If image width and height are unavailable,
    # height. This is explicitly labelled as estimated because a one-sided
    # satellite distribution can bias the detected-coordinate extent.
    xv = pd.to_numeric(group.get("Std_X", pd.Series(dtype=float)), errors="coerce")
    yv = pd.to_numeric(group.get("Std_Y", pd.Series(dtype=float)), errors="coerce")
    xv = xv[np.isfinite(xv)]
    yv = yv[np.isfinite(yv)]
    if not xv.empty and not yv.empty:
        return float((xv.min() + xv.max()) / 2.0), float((yv.min() + yv.max()) / 2.0), "Estimated image center (detected-coordinate extent fallback)"
    return np.nan, np.nan, "Unavailable"


def add_center_distance_metrics(df: pd.DataFrame, source_path: Path) -> pd.DataFrame:
    """Add XY distance from every satellite to the center of its cell or image.

    Reference hierarchy:
      1. Exact Fiji cell-boundary polygon centroid when a real cell ROI exists.
      2. The center used by centrosome-radius / centrosome-centered-sphere cell
         assignment, because that ROI is explicitly defined around that center.
      3. Image center when there is no cell assignment or when a cell center is
         unavailable. Exact image dimensions are preferred; otherwise use a
         clearly labelled coordinate-extent estimate.

    The per-satellite value is repeated in the standardized satellite table.
    Per-cell/image/frame summaries report the mean of these distances.
    """
    out = df.copy()
    out["Std_Distance_to_Cell_or_Image_Center_XY"] = np.nan
    out["Std_Center_Reference_X"] = np.nan
    out["Std_Center_Reference_Y"] = np.nan
    out["Std_Center_Reference_Type"] = ""
    out["Std_Center_Distance_Unit"] = ""
    if out.empty or "Std_X" not in out.columns or "Std_Y" not in out.columns:
        return out

    boundaries = _load_boundary_polygons_for_entropy(Path(source_path), out)
    group_source = "Source_UID" if "Source_UID" in out.columns else "Source_Input"
    group_cols = [group_source]
    if "Std_Frame" in out.columns and pd.to_numeric(out["Std_Frame"], errors="coerce").notna().any():
        group_cols.append("Std_Frame")

    for _, frame_group in out.groupby(group_cols, dropna=False, sort=False):
        coord_space = "um"
        if "Std_XY_Coordinate_Space" in frame_group.columns:
            vals = frame_group["Std_XY_Coordinate_Space"].dropna().astype(str).str.strip().str.lower()
            if not vals.empty and vals.iloc[0] in ["um", "pixel"]:
                coord_space = vals.iloc[0]
        unit = "µm" if coord_space == "um" else "pixels"
        image_cx, image_cy, image_source = _image_center_xy(frame_group, coord_space)

        # Treat each valid cell separately. Unassigned rows are processed as an
        # image-level group below.
        if "Std_Cell_ID" in frame_group.columns:
            ctext = frame_group["Std_Cell_ID"].astype(str).str.strip()
            valid_cell = frame_group["Std_Cell_ID"].notna() & ~ctext.str.lower().isin(["", "nan", "none", "na", "0", "-1"])
        else:
            valid_cell = pd.Series(False, index=frame_group.index)

        for cell_value, cell_group in frame_group.loc[valid_cell].groupby("Std_Cell_ID", dropna=False, sort=False):
            center_x = center_y = np.nan
            center_source = ""
            assignment = ""
            if "Cell_Assignment_Mode" in cell_group.columns:
                av = cell_group["Cell_Assignment_Mode"].dropna().astype(str)
                if not av.empty:
                    assignment = av.iloc[0].strip().lower()

            # "No cell assignment" often uses Cell_ID=1 as a whole-image group.
            no_cell_mode = any(k in assignment for k in ["no cell assignment", "skip cell assignment", "whole image"])
            if not no_cell_mode:
                poly = _choose_boundary_polygon(boundaries, cell_group)
                center_x, center_y = _polygon_centroid_xy(poly)
                if np.isfinite(center_x) and np.isfinite(center_y):
                    center_source = "Cell center (exact Fiji boundary centroid)"

            if (not np.isfinite(center_x) or not np.isfinite(center_y)) and any(k in assignment for k in ["centrosome radius", "centrosome-centered sphere", "centrosome centered sphere"]):
                cx = pd.to_numeric(cell_group.get("Std_Centrosome_X", pd.Series(np.nan, index=cell_group.index)), errors="coerce").dropna()
                cy = pd.to_numeric(cell_group.get("Std_Centrosome_Y", pd.Series(np.nan, index=cell_group.index)), errors="coerce").dropna()
                if not cx.empty and not cy.empty:
                    center_x = float(cx.median())
                    center_y = float(cy.median())
                    if coord_space == "pixel":
                        pw = _median_positive(cell_group, "Std_Pixel_Width_um")
                        ph = _median_positive(cell_group, "Std_Pixel_Height_um")
                        if np.isfinite(pw) and np.isfinite(ph):
                            center_x /= pw
                            center_y /= ph
                        else:
                            center_x = center_y = np.nan
                    if np.isfinite(center_x) and np.isfinite(center_y):
                        center_source = "Cell-assignment center (centrosome radius/sphere)"

            if not np.isfinite(center_x) or not np.isfinite(center_y):
                center_x, center_y, center_source = image_cx, image_cy, image_source
                if not no_cell_mode and center_source != "Unavailable":
                    center_source += " — fallback because cell center was unavailable"

            gx = pd.to_numeric(cell_group["Std_X"], errors="coerce")
            gy = pd.to_numeric(cell_group["Std_Y"], errors="coerce")
            dist = np.sqrt((gx - center_x) ** 2 + (gy - center_y) ** 2) if np.isfinite(center_x) and np.isfinite(center_y) else pd.Series(np.nan, index=cell_group.index)
            out.loc[cell_group.index, "Std_Distance_to_Cell_or_Image_Center_XY"] = dist
            out.loc[cell_group.index, "Std_Center_Reference_X"] = center_x
            out.loc[cell_group.index, "Std_Center_Reference_Y"] = center_y
            out.loc[cell_group.index, "Std_Center_Reference_Type"] = center_source
            out.loc[cell_group.index, "Std_Center_Distance_Unit"] = unit

        unassigned_idx = frame_group.index.difference(frame_group.loc[valid_cell].index)
        if len(unassigned_idx):
            ug = frame_group.loc[unassigned_idx]
            gx = pd.to_numeric(ug["Std_X"], errors="coerce")
            gy = pd.to_numeric(ug["Std_Y"], errors="coerce")
            dist = np.sqrt((gx - image_cx) ** 2 + (gy - image_cy) ** 2) if np.isfinite(image_cx) and np.isfinite(image_cy) else pd.Series(np.nan, index=ug.index)
            out.loc[ug.index, "Std_Distance_to_Cell_or_Image_Center_XY"] = dist
            out.loc[ug.index, "Std_Center_Reference_X"] = image_cx
            out.loc[ug.index, "Std_Center_Reference_Y"] = image_cy
            out.loc[ug.index, "Std_Center_Reference_Type"] = image_source
            out.loc[ug.index, "Std_Center_Distance_Unit"] = unit

    return out


def mean_unique_cell_metric(g: pd.DataFrame, column: str) -> float:
    """Average a repeated cell-level metric using each biological cell once.

    Cell_ID values restart at 1 for every Fiji input, so source identity must be
    part of the key. Frame is included when present for time-lapse data.
    """
    if column not in g.columns or "Std_Cell_ID" not in g.columns:
        return np.nan
    tmp = pd.DataFrame({
        "cell": g["Std_Cell_ID"].astype(str),
        "value": pd.to_numeric(g[column], errors="coerce"),
    }, index=g.index)
    subset = []
    if "Source_UID" in g.columns:
        tmp["source"] = g["Source_UID"].astype(str)
        subset.append("source")
    elif "Source_Input" in g.columns:
        tmp["source"] = g["Source_Input"].astype(str)
        subset.append("source")
    subset.append("cell")
    if "Std_Frame" in g.columns:
        tmp["frame"] = pd.to_numeric(g["Std_Frame"], errors="coerce")
        subset.append("frame")
    tmp = tmp.dropna(subset=["value"]).drop_duplicates(subset=subset)
    return float(tmp["value"].mean()) if not tmp.empty else np.nan


def summarize_overall(df: pd.DataFrame, input_count: int, graph_count: int) -> pd.DataFrame:
    rows = [
        ("Input files", input_count),
        ("Satellite rows", len(df)),
        ("Images / inputs", df["Source_Input"].nunique() if "Source_Input" in df.columns else 1),
        ("Modes", ", ".join(sorted(df["Mode"].dropna().astype(str).unique())) if "Mode" in df.columns else ""),
("Graphs and viewers created", graph_count),
    ]
    for label, col in [
        ("Mean 2D satellite area (µm²)", "Std_Area_2D"),
        ("Median 2D satellite area (µm²)", "Std_Area_2D"),
        ("Total 2D satellite area (µm²)", "Std_Area_2D"),
        ("Mean 3D satellite volume (µm³)", "Std_Volume_3D"),
        ("Median 3D satellite volume (µm³)", "Std_Volume_3D"),
        ("Total 3D satellite volume (µm³)", "Std_Volume_3D"),
        ("Mean integrated intensity (A.U.)", "Std_Intensity"),
        ("Total integrated intensity (A.U.)", "Std_Intensity"),
        ("Mean distance (µm)", "Std_Distance"),
        ("Mean 3D sphericity", "Std_Sphericity"),
        ("Mean 2D circularity", "Std_Circularity_2D"),
        ("Mean shape compactness (2D circularity / 3D sphericity)", "Std_Shape_Compactness"),
        ("Mean cell spatial entropy", "Std_Cell_Spatial_Entropy"),
        ("Mean pericentrosomal clustering index", "Std_Cell_Pericentrosomal_Clustering_Index"),
        ("Mean distance to cell/image center (XY)", "Std_Distance_to_Cell_or_Image_Center_XY"),
        ("Mean 2D cell area (µm²)", "Std_Cell_Area_2D"),
        ("Mean 3D cell volume (µm³)", "Std_Cell_Volume_3D"),
    ]:
        s = numeric(df, col).dropna()
        if s.empty:
            val = np.nan
        elif col in ["Std_Cell_Spatial_Entropy", "Std_Cell_Pericentrosomal_Clustering_Index"]:
            # These are repeated on every satellite row; summarize each cell/frame once.
            val = mean_unique_cell_metric(df, col)
        elif "Median" in label:
            val = s.median()
        elif "Total" in label:
            val = s.sum()
        else:
            val = s.mean()
        rows.append((label, val))
    return pd.DataFrame(rows, columns=["Metric", "Value"])

def summarize_by_image(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    group_key = "Source_UID" if "Source_UID" in df.columns else "Source_Input"
    for src_uid, g in df.groupby(group_key, dropna=False):
        display_name = g["Source_Input"].iloc[0] if "Source_Input" in g.columns else str(src_uid)
        rows.append({
            "Source_UID": src_uid,
            "Source_Input": display_name,
            "Source_Path": g["Source_Path"].iloc[0] if "Source_Path" in g.columns else "",
            "Source_Folder": g["Source_Folder"].iloc[0] if "Source_Folder" in g.columns else "",
            "Condition": g["Condition"].iloc[0] if "Condition" in g.columns else "",
            "Replicate": g["Replicate"].iloc[0] if "Replicate" in g.columns else "",
            "Mode": g["Mode"].iloc[0] if "Mode" in g.columns else "",
            "Satellite_count": len(g),
            "Assigned_satellites": g["Std_Cell_ID"].replace("", np.nan).dropna().shape[0] if "Std_Cell_ID" in g.columns else np.nan,
            "Cells": g["Std_Cell_ID"].replace("", np.nan).dropna().nunique() if "Std_Cell_ID" in g.columns else np.nan,
            "Mean_size": numeric(g, "Std_Size").mean(),
            "Median_size": numeric(g, "Std_Size").median(),
            "Total_size": numeric(g, "Std_Size").sum(min_count=1),
            "Mean_area_2D": numeric(g, "Std_Area_2D").mean(),
            "Median_area_2D": numeric(g, "Std_Area_2D").median(),
            "Total_area_2D": numeric(g, "Std_Area_2D").sum(min_count=1),
            "Mean_volume_3D": numeric(g, "Std_Volume_3D").mean(),
            "Median_volume_3D": numeric(g, "Std_Volume_3D").median(),
            "Total_volume_3D": numeric(g, "Std_Volume_3D").sum(min_count=1),
            "Mean_intensity": numeric(g, "Std_Intensity").mean(),
            "Total_intensity": numeric(g, "Std_Intensity").sum(min_count=1),
            "Mean_distance": numeric(g, "Std_Distance").mean(),
            "Mean_sphericity": numeric(g, "Std_Sphericity").mean(),
            "Mean_circularity_2D": numeric(g, "Std_Circularity_2D").mean(),
            "Mean_shape_compactness": numeric(g, "Std_Shape_Compactness").mean(),
            "Mean_cell_spatial_entropy": mean_unique_cell_metric(g, "Std_Cell_Spatial_Entropy"),
            "Mean_cell_pericentrosomal_clustering_index": mean_unique_cell_metric(g, "Std_Cell_Pericentrosomal_Clustering_Index"),
            "Mean_distance_to_cell_or_image_center": numeric(g, "Std_Distance_to_Cell_or_Image_Center_XY").mean(),
            "Ciliated_cell_count": float((pd.DataFrame({"cell":g.get("Std_Cell_ID",pd.Series(index=g.index,dtype=object)).map(_normalize_id_text),"v":numeric(g,"Std_Ciliated")}).dropna(subset=["v"]).drop_duplicates("cell")["v"]==1).sum()) if "Std_Ciliated" in g.columns else np.nan,
            "Cell_count_with_ciliation_data": float(pd.DataFrame({"cell":g.get("Std_Cell_ID",pd.Series(index=g.index,dtype=object)).map(_normalize_id_text),"v":numeric(g,"Std_Ciliated")}).dropna(subset=["v"]).drop_duplicates("cell").shape[0]) if "Std_Ciliated" in g.columns else np.nan,
            "Ciliation_fraction": float(pd.DataFrame({"cell":g.get("Std_Cell_ID",pd.Series(index=g.index,dtype=object)).map(_normalize_id_text),"v":numeric(g,"Std_Ciliated")}).dropna(subset=["v"]).drop_duplicates("cell")["v"].mean()) if "Std_Ciliated" in g.columns and pd.DataFrame({"cell":g.get("Std_Cell_ID",pd.Series(index=g.index,dtype=object)).map(_normalize_id_text),"v":numeric(g,"Std_Ciliated")}).dropna(subset=["v"]).drop_duplicates("cell").shape[0] else np.nan,
            "Mean_cilium_length_um": mean_unique_cell_metric(g, "Std_Cilium_Length_um") if "Std_Cilium_Length_um" in g.columns else np.nan,
        })
    return pd.DataFrame(rows)

def count_centrosomes_in_cell(g: pd.DataFrame) -> float:
    if "Std_Centrosome_ID" in g.columns:
        vals = g["Std_Centrosome_ID"].replace("", np.nan).dropna().astype(str)
        vals = vals[~vals.str.lower().isin(["nan", "none", "na"])]
        if not vals.empty:
            return float(vals.nunique())
    coord_cols = [c for c in ["Std_Centrosome_X", "Std_Centrosome_Y", "Std_Centrosome_Z"] if c in g.columns]
    if len(coord_cols) >= 2:
        tmp = g[coord_cols].apply(pd.to_numeric, errors="coerce").dropna(how="all")
        if not tmp.empty:
            tmp = tmp.round(3).drop_duplicates()
            return float(len(tmp))
    return np.nan

def summarize_by_cell(df: pd.DataFrame) -> pd.DataFrame:
    if "Std_Cell_ID" not in df.columns:
        return pd.DataFrame()
    d = df.copy()
    d["Std_Cell_ID"] = d["Std_Cell_ID"].replace("", np.nan)
    d = d.dropna(subset=["Std_Cell_ID"])
    if d.empty:
        return pd.DataFrame()

    group_key = "Source_UID" if "Source_UID" in d.columns else "Source_Input"

    def _short_source_label(g: pd.DataFrame) -> str:
        src = str(g["Source_Input"].iloc[0]) if "Source_Input" in g.columns else "Input"
        for suffix in [
            "_SatelliteQuantify_Fiji_Output.csv",
            "_timelapse_satellite_rows_combined.csv",
            "satellite_3D_results_clean.csv",
            ".csv", ".xlsx", ".xlsm"
        ]:
            src = src.replace(suffix, "")
        cond = str(g["Condition"].iloc[0]) if "Condition" in g.columns else ""
        rep = str(g["Replicate"].iloc[0]) if "Replicate" in g.columns else ""
        parts = []
        for part in [cond, rep, src]:
            if part and part.lower() not in ["nan", "none"] and part not in parts:
                parts.append(part)
        return " / ".join(parts) if parts else src

    rows = []
    group_cols = [group_key, "Std_Cell_ID"]
    frame_key = "Std_Frame" if "Std_Frame" in d.columns else None
    if frame_key is not None:
        group_cols.append(frame_key)
    for keys, g in d.groupby(group_cols, dropna=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        src_uid = keys[0]
        cell = keys[1]
        frame = keys[2] if frame_key is not None and len(keys) > 2 else np.nan
        vals = numeric(g, "Std_Cell_Measure").dropna()
        measure_label = "Cell measure"
        if "Std_Cell_Measure_Name" in g.columns and g["Std_Cell_Measure_Name"].replace("", np.nan).dropna().size:
            measure_label = str(g["Std_Cell_Measure_Name"].replace("", np.nan).dropna().iloc[0])
        source_label = _short_source_label(g)
        row = {
            "Source_UID": src_uid,
            "Source_Input": g["Source_Input"].iloc[0] if "Source_Input" in g.columns else str(src_uid),
            "Source_Path": g["Source_Path"].iloc[0] if "Source_Path" in g.columns else "",
            "Source_Folder": g["Source_Folder"].iloc[0] if "Source_Folder" in g.columns else "",
            "Condition": g["Condition"].iloc[0] if "Condition" in g.columns else "",
            "Replicate": g["Replicate"].iloc[0] if "Replicate" in g.columns else "",
            "Unique_Cell": source_label + " | Cell " + str(cell) + ((" | Frame " + _normalize_id_text(frame)) if frame_key is not None and _normalize_id_text(frame) else ""),
            "Cell_ID": cell,
            "Frame": frame,
            "Mode": g["Mode"].iloc[0] if "Mode" in g.columns else "",
            "Cell_measure": vals.iloc[0] if vals.size else np.nan,
            "Cell_measure_name": measure_label,
            "Cell_area_2D": numeric(g, "Std_Cell_Area_2D").dropna().iloc[0] if numeric(g, "Std_Cell_Area_2D").dropna().size else np.nan,
            "Cell_volume_3D": numeric(g, "Std_Cell_Volume_3D").dropna().iloc[0] if numeric(g, "Std_Cell_Volume_3D").dropna().size else np.nan,
            "Satellite_count": len(g),
            "Centrosome_number_per_cell": count_centrosomes_in_cell(g),
            "Mean_size": numeric(g, "Std_Size").mean(),
            "Median_size": numeric(g, "Std_Size").median(),
            "Total_size": numeric(g, "Std_Size").sum(min_count=1),
            "Mean_area_2D": numeric(g, "Std_Area_2D").mean(),
            "Median_area_2D": numeric(g, "Std_Area_2D").median(),
            "Total_area_2D": numeric(g, "Std_Area_2D").sum(min_count=1),
            "Mean_volume_3D": numeric(g, "Std_Volume_3D").mean(),
            "Median_volume_3D": numeric(g, "Std_Volume_3D").median(),
            "Total_volume_3D": numeric(g, "Std_Volume_3D").sum(min_count=1),
            "Mean_intensity": numeric(g, "Std_Intensity").mean(),
            "Total_intensity": numeric(g, "Std_Intensity").sum(min_count=1),
            "Mean_distance": numeric(g, "Std_Distance").mean(),
            "Mean_sphericity": numeric(g, "Std_Sphericity").mean(),
            "Mean_circularity_2D": numeric(g, "Std_Circularity_2D").mean(),
            "Mean_shape_compactness": numeric(g, "Std_Shape_Compactness").mean(),
            "Cell_spatial_entropy": numeric(g, "Std_Cell_Spatial_Entropy").dropna().iloc[0] if numeric(g, "Std_Cell_Spatial_Entropy").dropna().size else np.nan,
            "Cell_pericentrosomal_clustering_index": numeric(g, "Std_Cell_Pericentrosomal_Clustering_Index").dropna().iloc[0] if numeric(g, "Std_Cell_Pericentrosomal_Clustering_Index").dropna().size else np.nan,
            PCI_CELL_INNER_COL: numeric(g, PCI_STD_INNER_COL).dropna().iloc[0] if numeric(g, PCI_STD_INNER_COL).dropna().size else np.nan,
            PCI_CELL_OUTER_COL: numeric(g, PCI_STD_OUTER_COL).dropna().iloc[0] if numeric(g, PCI_STD_OUTER_COL).dropna().size else np.nan,
            PCI_CELL_TOTAL_COL: numeric(g, PCI_STD_TOTAL_COL).dropna().iloc[0] if numeric(g, PCI_STD_TOTAL_COL).dropna().size else np.nan,
            "Pericentrosomal_geometry": g["Std_Cell_Pericentrosomal_Geometry"].replace("",np.nan).dropna().astype(str).iloc[0] if "Std_Cell_Pericentrosomal_Geometry" in g.columns and g["Std_Cell_Pericentrosomal_Geometry"].replace("",np.nan).dropna().size else "",
            "Pericentrosomal_clustering_status": g["Std_Cell_Pericentrosomal_Clustering_Status"].replace("",np.nan).dropna().astype(str).iloc[0] if "Std_Cell_Pericentrosomal_Clustering_Status" in g.columns and g["Std_Cell_Pericentrosomal_Clustering_Status"].replace("",np.nan).dropna().size else "",
            "Cilium_count": numeric(g,"Std_Cilium_Count").dropna().iloc[0] if numeric(g,"Std_Cilium_Count").dropna().size else np.nan,
            "Cilium_length_um": numeric(g,"Std_Cilium_Length_um").dropna().iloc[0] if numeric(g,"Std_Cilium_Length_um").dropna().size else np.nan,
            "Mean_cilium_length_um": numeric(g,"Std_Cilium_Mean_Length_um").dropna().iloc[0] if numeric(g,"Std_Cilium_Mean_Length_um").dropna().size else np.nan,
            "Ciliated": numeric(g,"Std_Ciliated").dropna().iloc[0] if numeric(g,"Std_Ciliated").dropna().size else np.nan,
            "Ciliation_status": g["Ciliation_status"].replace("",np.nan).dropna().astype(str).iloc[0] if "Ciliation_status" in g.columns and g["Ciliation_status"].replace("",np.nan).dropna().size else "",
            "Ciliation_length_threshold_um": numeric(g,"Std_Ciliation_Length_Threshold_um").dropna().iloc[0] if numeric(g,"Std_Ciliation_Length_Threshold_um").dropna().size else np.nan,
            "Mean_distance_to_cell_or_image_center": numeric(g, "Std_Distance_to_Cell_or_Image_Center_XY").mean(),
            "Center_reference_type": g["Std_Center_Reference_Type"].dropna().astype(str).iloc[0] if "Std_Center_Reference_Type" in g.columns and g["Std_Center_Reference_Type"].dropna().size else "",
        }
        if measure_label in ["Cell area", "Cell volume"]:
            row[measure_label] = vals.iloc[0] if vals.size else np.nan
        rows.append(row)
    return pd.DataFrame(rows)

def make_spatial_metrics_table(cell_summary: pd.DataFrame) -> pd.DataFrame:
    """One explicit row per cell with SatelliteQ's two primary spatial metrics."""
    wanted = [
        "Source_Input", "Source_Path", "Source_Folder", "Condition", "Replicate",
        "Unique_Cell", "Cell_ID", "Mode",
        "Cell_spatial_entropy", "Cell_pericentrosomal_clustering_index",
        PCI_CELL_INNER_COL, PCI_CELL_OUTER_COL, PCI_CELL_TOTAL_COL,
        "Pericentrosomal_geometry", "Pericentrosomal_clustering_status",
    ]
    if cell_summary is None or cell_summary.empty:
        return pd.DataFrame(columns=wanted)
    out = cell_summary.copy()
    for col in wanted:
        if col not in out.columns:
            out[col] = np.nan
    return out[wanted].copy()


def mean_unique_cell_measure(g: pd.DataFrame) -> float:
    """Mean cell area/volume for one grouped frame/image, using each cell once."""
    if "Std_Cell_Measure" not in g.columns:
        return np.nan
    d = g.copy()
    vals = pd.to_numeric(d["Std_Cell_Measure"], errors="coerce")
    if vals.notna().sum() == 0:
        return np.nan
    if "Std_Cell_ID" in d.columns:
        tmp = pd.DataFrame({"cell": d["Std_Cell_ID"].astype(str), "measure": vals}, index=d.index)
        key_cols = []
        if "Source_UID" in d.columns:
            tmp["source"] = d["Source_UID"].astype(str); key_cols.append("source")
        elif "Source_Input" in d.columns:
            tmp["source"] = d["Source_Input"].astype(str); key_cols.append("source")
        key_cols.append("cell")
        if "Std_Frame" in d.columns:
            tmp["frame"] = pd.to_numeric(d["Std_Frame"], errors="coerce"); key_cols.append("frame")
        tmp = tmp.dropna(subset=["measure"]).drop_duplicates(subset=key_cols)
        return float(tmp["measure"].mean()) if not tmp.empty else np.nan
    return float(vals.mean())

def total_unique_cell_measure(g: pd.DataFrame) -> float:
    """Total cell area/volume for one grouped frame/image, using each cell once."""
    if "Std_Cell_Measure" not in g.columns:
        return np.nan
    d = g.copy()
    vals = pd.to_numeric(d["Std_Cell_Measure"], errors="coerce")
    if vals.notna().sum() == 0:
        return np.nan
    if "Std_Cell_ID" in d.columns:
        tmp = pd.DataFrame({"cell": d["Std_Cell_ID"].astype(str), "measure": vals}, index=d.index)
        key_cols = []
        if "Source_UID" in d.columns:
            tmp["source"] = d["Source_UID"].astype(str); key_cols.append("source")
        elif "Source_Input" in d.columns:
            tmp["source"] = d["Source_Input"].astype(str); key_cols.append("source")
        key_cols.append("cell")
        if "Std_Frame" in d.columns:
            tmp["frame"] = pd.to_numeric(d["Std_Frame"], errors="coerce"); key_cols.append("frame")
        tmp = tmp.dropna(subset=["measure"]).drop_duplicates(subset=key_cols)
        return float(tmp["measure"].sum()) if not tmp.empty else np.nan
    return float(vals.sum())

def mean_unique_cell_measure_for_column(g: pd.DataFrame, column: str) -> float:
    """Mean an explicit cell area/volume column using each cell once."""
    if column not in g.columns:
        return np.nan
    d = g.copy()
    vals = pd.to_numeric(d[column], errors="coerce")
    if vals.notna().sum() == 0:
        return np.nan
    if "Std_Cell_ID" in d.columns:
        tmp = pd.DataFrame({"cell": d["Std_Cell_ID"].astype(str), "measure": vals}, index=d.index)
        key_cols = []
        if "Source_UID" in d.columns:
            tmp["source"] = d["Source_UID"].astype(str); key_cols.append("source")
        elif "Source_Input" in d.columns:
            tmp["source"] = d["Source_Input"].astype(str); key_cols.append("source")
        key_cols.append("cell")
        if "Std_Frame" in d.columns:
            tmp["frame"] = pd.to_numeric(d["Std_Frame"], errors="coerce"); key_cols.append("frame")
        tmp = tmp.dropna(subset=["measure"]).drop_duplicates(subset=key_cols)
        return float(tmp["measure"].mean()) if not tmp.empty else np.nan
    return float(vals.mean())


def total_unique_cell_measure_for_column(g: pd.DataFrame, column: str) -> float:
    """Total an explicit cell area/volume column using each cell once."""
    if column not in g.columns:
        return np.nan
    d = g.copy()
    vals = pd.to_numeric(d[column], errors="coerce")
    if vals.notna().sum() == 0:
        return np.nan
    if "Std_Cell_ID" in d.columns:
        tmp = pd.DataFrame({"cell": d["Std_Cell_ID"].astype(str), "measure": vals}, index=d.index)
        key_cols = []
        if "Source_UID" in d.columns:
            tmp["source"] = d["Source_UID"].astype(str); key_cols.append("source")
        elif "Source_Input" in d.columns:
            tmp["source"] = d["Source_Input"].astype(str); key_cols.append("source")
        key_cols.append("cell")
        if "Std_Frame" in d.columns:
            tmp["frame"] = pd.to_numeric(d["Std_Frame"], errors="coerce"); key_cols.append("frame")
        tmp = tmp.dropna(subset=["measure"]).drop_duplicates(subset=key_cols)
        return float(tmp["measure"].sum()) if not tmp.empty else np.nan
    return float(vals.sum())


def mean_centrosome_number_per_cell(g: pd.DataFrame) -> float:
    if "Std_Cell_ID" not in g.columns:
        return np.nan
    d = g.copy()
    d["Std_Cell_ID"] = d["Std_Cell_ID"].replace("", np.nan)
    d = d.dropna(subset=["Std_Cell_ID"])
    if d.empty:
        return np.nan
    vals = []
    group_cols = []
    if "Source_UID" in d.columns:
        group_cols.append("Source_UID")
    elif "Source_Input" in d.columns:
        group_cols.append("Source_Input")
    group_cols.append("Std_Cell_ID")
    if "Std_Frame" in d.columns:
        group_cols.append("Std_Frame")
    for _, gg in d.groupby(group_cols, dropna=False):
        v = count_centrosomes_in_cell(gg)
        if pd.notna(v):
            vals.append(v)
    return float(np.mean(vals)) if vals else np.nan

def summarize_by_frame(df: pd.DataFrame) -> pd.DataFrame:
    """Summarize time-lapse satellite rows by experiment/replicate/frame.

    For time-lapse analyses, Fiji may export each frame as a separate file.
    This groups those frame files back into one continuous experiment using
    Condition + Replicate + Frame, instead of treating each frame file as a
    separate condition.
    """
    if "Std_Frame" not in df.columns or numeric(df, "Std_Frame").dropna().empty:
        return pd.DataFrame()
    rows = []
    d = df.dropna(subset=["Std_Frame"]).copy()
    group_cols = [c for c in ["Condition", "Replicate", "Std_Frame"] if c in d.columns]
    if "Std_Frame" not in group_cols:
        group_cols = ["Std_Frame"]
    for keys, g in d.groupby(group_cols, dropna=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        key_map = dict(zip(group_cols, keys))
        frame = key_map.get("Std_Frame", np.nan)
        cond = key_map.get("Condition", g["Condition"].iloc[0] if "Condition" in g.columns else "")
        rep = key_map.get("Replicate", g["Replicate"].iloc[0] if "Replicate" in g.columns else "")
        rows.append({
            "Condition": cond,
            "Replicate": rep,
            "Frame": frame,
            "Time_min": numeric(g, "Std_Time_min").dropna().iloc[0] if numeric(g, "Std_Time_min").dropna().size else np.nan,
            "Source_Input": "; ".join(sorted(set(g["Source_Input"].astype(str)))) if "Source_Input" in g.columns else "",
            "Source_Path": "; ".join(sorted(set(g["Source_Path"].astype(str)))) if "Source_Path" in g.columns else "",
            "Frame_file_count": g["Source_UID"].nunique() if "Source_UID" in g.columns else g["Source_Input"].nunique() if "Source_Input" in g.columns else np.nan,
            "Satellite_count": len(g),
            "Mean_size": numeric(g, "Std_Size").mean(),
            "Median_size": numeric(g, "Std_Size").median(),
            "Total_size": numeric(g, "Std_Size").sum(min_count=1),
            "Mean_area_2D": numeric(g, "Std_Area_2D").mean(),
            "Median_area_2D": numeric(g, "Std_Area_2D").median(),
            "Total_area_2D": numeric(g, "Std_Area_2D").sum(min_count=1),
            "Mean_volume_3D": numeric(g, "Std_Volume_3D").mean(),
            "Median_volume_3D": numeric(g, "Std_Volume_3D").median(),
            "Total_volume_3D": numeric(g, "Std_Volume_3D").sum(min_count=1),
            "Mean_intensity": numeric(g, "Std_Intensity").mean(),
            "Total_intensity": numeric(g, "Std_Intensity").sum(min_count=1),
            "Mean_distance": numeric(g, "Std_Distance").mean(),
            "Mean_sphericity": numeric(g, "Std_Sphericity").mean(),
            "Mean_circularity_2D": numeric(g, "Std_Circularity_2D").mean(),
            "Mean_shape_compactness": numeric(g, "Std_Shape_Compactness").mean(),
            "Mean_cell_spatial_entropy": mean_unique_cell_metric(g, "Std_Cell_Spatial_Entropy"),
            "Mean_cell_pericentrosomal_clustering_index": mean_unique_cell_metric(g, "Std_Cell_Pericentrosomal_Clustering_Index"),
            "Mean_distance_to_cell_or_image_center": numeric(g, "Std_Distance_to_Cell_or_Image_Center_XY").mean(),
            "Mean_centrosome_number_per_cell": mean_centrosome_number_per_cell(g),
            "Mean_cell_size": mean_unique_cell_measure(g),
            "Total_cell_size": total_unique_cell_measure(g),
            "Mean_cell_area_2D": mean_unique_cell_measure_for_column(g, "Std_Cell_Area_2D"),
            "Total_cell_area_2D": total_unique_cell_measure_for_column(g, "Std_Cell_Area_2D"),
            "Mean_cell_volume_3D": mean_unique_cell_measure_for_column(g, "Std_Cell_Volume_3D"),
            "Total_cell_volume_3D": total_unique_cell_measure_for_column(g, "Std_Cell_Volume_3D"),
        })
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    return out.sort_values([c for c in ["Condition", "Replicate", "Frame"] if c in out.columns])

def mode_label(df: pd.DataFrame) -> str:
    if "Mode" not in df.columns:
        return ""
    vals = sorted(df["Mode"].dropna().astype(str).unique())
    if len(vals) == 1:
        return vals[0]
    if len(vals) > 1:
        return "mixed"
    return ""

def cell_measure_label(df: pd.DataFrame) -> str:

    for col in ["Std_Cell_Measure_Name", "Cell_measure_name"]:
        if col in df.columns:
            vals = df[col].replace("", np.nan).dropna().astype(str).unique()
            vals = [v for v in vals if v.lower() not in ["nan", "none"]]
            if len(vals) == 1:
                if vals[0] == "Cell area":
                    return "Cell area (µm²)"
                if vals[0] == "Cell volume":
                    return "Cell volume (µm³)"
                return vals[0]
            if len(vals) > 1:
                return "Cell area/volume (µm² or µm³)"
    m = mode_label(df)
    if m == "2D":
        return "Cell area (µm²)"
    if m == "3D":
        return "Cell volume (µm³)"
    return "Cell area/volume (µm² or µm³)"

def size_label(df: pd.DataFrame) -> str:
    m = mode_label(df)
    if m == "2D":
        return "Satellite area (µm²)"
    if m == "3D":
        return "Satellite volume (µm³)"
    return "Satellite area/volume (µm² or µm³)"

def label_for(col: str, df: pd.DataFrame, context: str = "satellite") -> str:
    labels = {
        "Source_Input": "Image / Fiji result file",
        "Source_Folder": "Input folder",
        "Source_short": "Fiji result",
        "Plot_Image": "Fiji result / image",
        "Unique_Cell": "Unique cell ID (file + cell)",
        "Plot_Cell": "Unique cell ID (file + cell)",
        "Unique_Satellite": "Unique satellite ID (file + satellite)",
        "Unique_Frame": "Unique frame ID (file + frame)",
        "Frame": "Frame number",
        "Std_Frame": "Frame number",
        "Time_min": "Time (min)",
        "Std_Time_min": "Time (min)",
        "Std_Satellite_ID": "Satellite ID",
        "Satellite_ID": "Satellite ID",
        "Std_Cell_ID": "Cell ID",
        "Cell_ID": "Cell ID",
        "Std_Centrosome_ID": "Centrosome ID",
        "Centrosome_ID": "Centrosome ID",
        "Satellite_count": "Number of satellites",
        "Assigned_satellites": "Number of assigned satellites",
        "Cells": "Number of cells",
        "Centrosome_number_per_cell": "Centrosome number per cell",
        "Mean_centrosome_number_per_cell": "Mean centrosome number per cell",
        "Std_Intensity": "Satellite integrated intensity (a.u.)",
        "Std_Total_Intensity": "Satellite total intensity (a.u.)",
        "Integrated_Intensity": "Satellite integrated intensity (a.u.)",
        "IntDen": "Satellite integrated intensity (a.u.)",
        "Integrated_Density": "Satellite integrated intensity (a.u.)",
        "Std_Mean_Intensity": "Satellite mean intensity (a.u./pixel or voxel)",
        "Mean_Intensity": "Satellite mean intensity (a.u./pixel or voxel)",
        "Mean": "Mean intensity (a.u./pixel or voxel)",
        "Std_Distance": "Satellite-to-centrosome distance (µm)",
        "Nearest_Centrosome_Distance": "Satellite-to-centrosome distance (µm)",
        "Distance_to_Assigned_Centrosome_3D": "Satellite-to-centrosome distance (µm)",
        "Distance_to_nearest_centrosome_3D": "Satellite-to-centrosome distance (µm)",
        "Std_Diameter": "Satellite equivalent diameter (µm)",
        "Equivalent_Circle_Diameter": "Satellite equivalent diameter (µm)",
        "Equivalent_Sphere_Diameter_um": "Satellite equivalent diameter (µm)",
        "Std_Area_2D": "Satellite area (2D; µm²)",
        "Std_Volume_3D": "Satellite volume (3D; µm³)",
        "Std_Cell_Area_2D": "Cell area (2D; µm²)",
        "Std_Cell_Volume_3D": "Cell volume (3D; µm³)",
        "Std_Surface": "Satellite surface/surface-area value reported by Fiji (3D)",
        "Std_Sphericity": "Satellite sphericity (3D; Morph_Sphericity preferred; 1 = sphere)",
        "Std_Sphericity_Source": "Source used for standardized 3D sphericity",
        "Std_Circularity_2D": "Satellite circularity (2D; 1 = circle)",
        "Std_Shape_Compactness": "Satellite shape compactness (2D circularity or 3D sphericity)",
        "Std_Cell_Spatial_Entropy": "Cell spatial entropy (intensity-weighted, multiscale; 0 concentrated → 1 spread)",
        "Std_Cell_Pericentrosomal_Clustering_Index": "Pericentrosomal clustering index (0–3 µm / 3–12 µm satellite fluorescence)",
        PCI_STD_INNER_COL: "Satellite fluorescence within 0–3 µm of assigned centrosome(s) (a.u.)",
        PCI_STD_OUTER_COL: "Satellite fluorescence within 3–12 µm of assigned centrosome(s) (a.u.)",
        PCI_STD_TOTAL_COL: "Satellite fluorescence within 0–12 µm of assigned centrosome(s) (a.u.)",
        "Cell_pericentrosomal_clustering_index": "Pericentrosomal clustering index",
        PCI_CELL_INNER_COL: "Pericentrosomal satellite fluorescence, 0–3 µm (a.u.)",
        PCI_CELL_OUTER_COL: "Outer satellite fluorescence, 3–12 µm (a.u.)",
        PCI_CELL_TOTAL_COL: "Total satellite fluorescence, 0–12 µm (a.u.)",
        "Pericentrosomal_geometry": "Pericentrosomal clustering distance geometry",
        "Pericentrosomal_clustering_status": "Pericentrosomal clustering index status",
        "Mean_cell_pericentrosomal_clustering_index": "Mean pericentrosomal clustering index",
        "Std_Distance_to_Cell_or_Image_Center_XY": "Distance to cell/image center (XY)",
        "Mean_distance_to_cell_or_image_center": "Mean distance of satellites to cell/image center (XY)",
        "Mean_sphericity": "Mean satellite sphericity (3D)",
        "Mean_circularity_2D": "Mean satellite circularity (2D)",
        "Mean_shape_compactness": "Mean satellite shape compactness (2D circularity / 3D sphericity)",
        "Cell_spatial_entropy": "Cell spatial entropy (0–1)",
        "Mean_cell_spatial_entropy": "Mean cell spatial entropy (0–1)",
        "Surface (micron^2)": "Satellite surface area (µm²)",
        "Surface": "Satellite surface area (µm²)",
        "Std_X": "X position (µm or pixels)",
        "Std_Y": "Y position (µm or pixels)",
        "Std_Z": "Z position (slice or µm)",
        "Mean_size": "Mean satellite size",
        "Median_size": "Median satellite size",
        "Total_size": "Total satellite size",
        "Mean_area_2D": "Mean satellite area (2D; µm²)",
        "Median_area_2D": "Median satellite area (2D; µm²)",
        "Total_area_2D": "Total satellite area (2D; µm²)",
        "Mean_volume_3D": "Mean satellite volume (3D; µm³)",
        "Median_volume_3D": "Median satellite volume (3D; µm³)",
        "Total_volume_3D": "Total satellite volume (3D; µm³)",
        "Cell_area_2D": "Cell area (2D; µm²)",
        "Cell_volume_3D": "Cell volume (3D; µm³)",
        "Mean_cell_area_2D": "Mean cell area (2D; µm²)",
        "Total_cell_area_2D": "Total cell area (2D; µm²)",
        "Mean_cell_volume_3D": "Mean cell volume (3D; µm³)",
        "Total_cell_volume_3D": "Total cell volume (3D; µm³)",
        "Mean_intensity": "Mean satellite integrated intensity (a.u.)",
        "Total_intensity": "Total satellite integrated intensity (a.u.)",
        "Std_Total_Intensity": "Satellite total intensity (a.u.)",
        "Total_intensity_A": "Total channel A intensity inside cell/ROI (a.u.)",
        "Total_intensity_B": "Total channel B intensity inside cell/ROI (a.u.)",
        "Integrated_A": "Total channel A intensity inside cell/ROI (a.u.)",
        "Integrated_B": "Total channel B intensity inside cell/ROI (a.u.)",
        "ROI_area_um2": "Cell/ROI area (µm²)",
        "ROI_volume_um3": "Cell/ROI volume (µm³)",
        "Cell_size": "Cell/ROI size (area or volume)",
        "Mean_cell_size": "Mean cell area/volume per frame",
        "Total_cell_size": "Total cell area/volume per frame",
        "Mean_distance": "Mean satellite-to-centrosome distance (µm)",
        "Std_Cilium_Count": "Cilium count per cell",
        "Std_Cilium_Length_um": "Cilium length (longest cilium per cell; µm)",
        "Std_Cilium_Mean_Length_um": "Mean cilium length per cell (µm)",
        "Std_Ciliated": "Ciliated (1=yes, 0=no)",
        "Ciliated": "Ciliated (1=yes, 0=no; mean = ciliated fraction)",
        "Ciliation_status": "Ciliation status",
        "Cilium_count": "Cilium count per cell",
        "Cilium_length_um": "Cilium length (µm)",
        "Mean_cilium_length_um": "Mean cilium length (µm)",
        "Ciliation_fraction": "Fraction of cells classified as ciliated",
        "Outlier_Cell": "Outlier cell",
        "Outlier_Score": "Outlier score",
        "Outlier_Reason": "Outlier reason",
        "Outlier_Label": "Outlier label",
    }
    if col == "Std_Area_2D":
        return "Satellite area (2D; µm²)"
    if col == "Std_Volume_3D":
        return "Satellite volume (3D; µm³)"
    if col == "Std_Cell_Area_2D":
        return "Cell area (2D; µm²)"
    if col == "Std_Cell_Volume_3D":
        return "Cell volume (3D; µm³)"
    if col in ["Std_Size", "Area", "Volume", "Volume (micron^3)"]:
        return size_label(df)
    if col in ["Std_Cell_Measure", "Cell_measure", "Cell_Area", "Cell_Volume"]:
        return cell_measure_label(df)
    if col in ["Mean_size", "Median_size", "Total_size"]:
        return col.replace("_", " ").capitalize() + " (µm² for 2D / µm³ for 3D)"
    return labels.get(col, str(col).replace("_", " "))

def linear_regression_r2(x_values, y_values):
    """Return slope, intercept, R², and n for usable numeric x/y values.

    This uses a centered least-squares calculation instead of np.polyfit,
    so Python does not print distracting "Polyfit may be poorly conditioned"
    warnings when x values are large, nearly repeated, or categorical-like.
    """
    x_arr = np.asarray(x_values, dtype=float)
    y_arr = np.asarray(y_values, dtype=float)
    mask = np.isfinite(x_arr) & np.isfinite(y_arr)
    x_arr = x_arr[mask]
    y_arr = y_arr[mask]

    if len(x_arr) < 3:
        return None

    if np.nanstd(x_arr) == 0 or np.nanstd(y_arr) == 0:
        return None

    x_mean = float(np.mean(x_arr))
    y_mean = float(np.mean(y_arr))
    x_centered = x_arr - x_mean
    y_centered = y_arr - y_mean
    denom = float(np.sum(x_centered ** 2))

    if denom == 0 or not np.isfinite(denom):
        return None

    slope = float(np.sum(x_centered * y_centered) / denom)
    intercept = float(y_mean - slope * x_mean)

    y_pred = slope * x_arr + intercept
    ss_res = float(np.sum((y_arr - y_pred) ** 2))
    ss_tot = float(np.sum((y_arr - y_mean) ** 2))

    if ss_tot == 0 or not np.isfinite(ss_tot):
        return None

    r2 = 1.0 - (ss_res / ss_tot)

    if not np.isfinite(r2):
        return None

    return slope, intercept, r2, len(x_arr)

def add_r2_to_scatter(ax, x_values, y_values):
    """Add a linear trendline and R² label to relationship scatter plots."""
    fit = linear_regression_r2(x_values, y_values)
    if fit is None:
        return None
    slope, intercept, r2, n = fit
    x_arr = np.asarray(x_values, dtype=float)
    mask = np.isfinite(x_arr)
    x_min = float(np.nanmin(x_arr[mask]))
    x_max = float(np.nanmax(x_arr[mask]))
    if x_min == x_max:
        return None
    line_x = np.array([x_min, x_max])
    line_y = slope * line_x + intercept
    ax.plot(line_x, line_y, linestyle="--", linewidth=1.4, color=SATQ_PURPLE_DARK)
    ax.text(
        0.04, 0.96,
        f"R² = {r2:.3f}\nn = {n}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        color=SATQ_PURPLE_DARK,
        bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor=SATQ_LIGHT_PURPLE, alpha=0.9),
    )
    return r2

def graphpad_axes(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#444444")
    ax.spines["bottom"].set_color("#444444")
    ax.grid(False)
    ax.tick_params(axis="both", labelsize=9, colors="#333333")
    ax.xaxis.label.set_size(10)
    ax.yaxis.label.set_size(10)
    ax.title.set_size(11)
    ax.title.set_color(SATQ_PURPLE)
    ax.title.set_weight("bold")

def save_fig(fig, path: Path) -> Optional[Path]:
    path.parent.mkdir(parents=True, exist_ok=True)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        try:
            fig.tight_layout()
        except Exception:
            pass

    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path

def scatter(df: pd.DataFrame, x: str, y: str, outdir: Path, title: str, xlabel: Optional[str] = None, ylabel: Optional[str] = None) -> Optional[Path]:
    if x not in df.columns or y not in df.columns:
        return None
    plot_df = sanitize_cilium_length_metrics(df)
    xv = pd.to_numeric(plot_df[x], errors="coerce")
    yv = pd.to_numeric(plot_df[y], errors="coerce")
    mask = xv.notna() & yv.notna()
    cilium_length_names = {
        "Cilium_length_um", "Mean_cilium_length_um",
        "Std_Cilium_Length_um", "Std_Cilium_Mean_Length_um",
        "Cilium_Length_um", "Mean_Cilium_Length_um",
    }
    if x in cilium_length_names:
        mask &= xv > 0
    if y in cilium_length_names:
        mask &= yv > 0
    if mask.sum() < 2:
        return None
    fig, ax = plt.subplots(figsize=(5.8, 4.7))
    ax.scatter(xv[mask], yv[mask], s=38, alpha=0.82, color=SATQ_LIGHT_PURPLE, edgecolors=SATQ_PURPLE, linewidths=0.7)
    add_r2_to_scatter(ax, xv[mask].to_numpy(), yv[mask].to_numpy())
    ax.set_title(title)
    ax.set_xlabel(xlabel or label_for(x, df))
    ax.set_ylabel(ylabel or label_for(y, df))
    graphpad_axes(ax)
    return save_fig(fig, outdir / f"{safe_name(title)}.png")

def bar(df: pd.DataFrame, x: str, y: str, outdir: Path, title: str, xlabel: Optional[str] = None, ylabel: Optional[str] = None) -> Optional[Path]:
    if x not in df.columns or y not in df.columns:
        return None
    d = df[[x, y]].copy()
    d[y] = pd.to_numeric(d[y], errors="coerce")
    d = d.dropna(subset=[x, y])
    if d.empty:
        return None
    fig_w = max(6.5, min(18, 0.60 * len(d) + 4))
    fig, ax = plt.subplots(figsize=(fig_w, 4.9))
    ax.bar(np.arange(len(d)), d[y], color=SATQ_LIGHT_PURPLE, edgecolor=SATQ_PURPLE, linewidth=0.8)
    ax.set_xticks(np.arange(len(d)))
    ax.set_xticklabels(d[x].astype(str), rotation=45, ha="right")
    ax.set_title(title)
    ax.set_xlabel(xlabel or label_for(x, df))
    ax.set_ylabel(ylabel or label_for(y, df))
    graphpad_axes(ax)
    return save_fig(fig, outdir / f"{safe_name(title)}.png")

def hist(df: pd.DataFrame, col: str, outdir: Path, title: str, xlabel: Optional[str] = None) -> Optional[Path]:
    if col not in df.columns:
        return None
    v = pd.to_numeric(df[col], errors="coerce").dropna()
    if v.empty:
        return None
    fig, ax = plt.subplots(figsize=(5.8, 4.7))
    ax.hist(v, bins=25, color=SATQ_LIGHT_PURPLE, edgecolor=SATQ_PURPLE, linewidth=0.8)
    ax.set_title(title)
    ax.set_xlabel(xlabel or label_for(col, df))
    ax.set_ylabel("Number of satellites")
    graphpad_axes(ax)
    return save_fig(fig, outdir / f"{safe_name(title)}.png")

def line(df: pd.DataFrame, x: str, y: str, group: Optional[str], outdir: Path, title: str, xlabel: Optional[str] = None, ylabel: Optional[str] = None) -> Optional[Path]:
    if x not in df.columns or y not in df.columns:
        return None
    d = df.copy()
    d[x] = pd.to_numeric(d[x], errors="coerce")
    d[y] = pd.to_numeric(d[y], errors="coerce")
    d = d.dropna(subset=[x, y])
    if d.empty:
        return None
    fig, ax = plt.subplots(figsize=(6.8, 4.8))
    if group and group in d.columns:
        for i, (name, g) in enumerate(d.groupby(group)):
            g = g.sort_values(x)
            c = SATQ_PURPLE_PALETTE[i % len(SATQ_PURPLE_PALETTE)]
            ax.plot(
                g[x], g[y],
                marker="o",
                linewidth=1.7,
                label=str(name),
                color=c,
                markerfacecolor=SATQ_LIGHT_PURPLE,
                markeredgecolor=c,
                markeredgewidth=0.8,
            )
        if d[group].nunique() <= 12:
            ax.legend(frameon=False, fontsize=8)
    else:
        d = d.sort_values(x)
        ax.plot(
            d[x], d[y],
            marker="o",
            linewidth=1.7,
            color=SATQ_PURPLE,
            markerfacecolor=SATQ_LIGHT_PURPLE,
            markeredgecolor=SATQ_PURPLE,
            markeredgewidth=0.8,
        )
    ax.set_title(title)
    ax.set_xlabel(xlabel or label_for(x, df))
    ax.set_ylabel(ylabel or label_for(y, df))
    graphpad_axes(ax)
    return save_fig(fig, outdir / f"{safe_name(title)}.png")


def _satq_numeric_series(df: pd.DataFrame, col: str) -> pd.Series:
    """Numeric helper used by the interactive outlier overlay."""
    if df is None or col not in df.columns:
        return pd.Series(dtype=float)
    return pd.to_numeric(df[col], errors="coerce")


def _satq_robust_zscore(series: pd.Series) -> pd.Series:
    """Median/MAD robust z-score; falls back to standard z-score when MAD is zero."""
    x = pd.to_numeric(series, errors="coerce")
    out = pd.Series(np.nan, index=x.index, dtype=float)
    clean = x.dropna()
    if clean.size < 5:
        return out
    med = float(clean.median())
    mad = float(np.median(np.abs(clean - med)))
    if mad and np.isfinite(mad):
        out.loc[clean.index] = 0.6745 * (clean - med) / mad
        return out
    sd = float(clean.std(ddof=1))
    if sd and np.isfinite(sd):
        out.loc[clean.index] = (clean - float(clean.mean())) / sd
    return out


def _satq_iqr_bounds(series: pd.Series, k: float = 1.5) -> Tuple[float, float]:
    """IQR outlier bounds for a numeric series."""
    clean = pd.to_numeric(series, errors="coerce").dropna()
    if clean.size < 5:
        return np.nan, np.nan
    q1 = float(clean.quantile(0.25))
    q3 = float(clean.quantile(0.75))
    iqr = q3 - q1
    if not np.isfinite(iqr) or iqr == 0:
        return np.nan, np.nan
    return q1 - k * iqr, q3 + k * iqr


def add_satelliteq_outlier_flags(
    df: pd.DataFrame,
    cell_summary: Optional[pd.DataFrame],
    image_summary: Optional[pd.DataFrame] = None,
    frame_summary: Optional[pd.DataFrame] = None,
) -> Tuple[pd.DataFrame, Optional[pd.DataFrame], Optional[pd.DataFrame], Optional[pd.DataFrame]]:
    """Add cell-level outlier annotations used by the interactive graph viewer.

    The detector is intentionally conservative and cell-level:
      - flags unusual cells based on satellite count, mean distance, mean/total size, and mean/total intensity when available
      - uses robust z-score and IQR rules
      - propagates outlier status from per-cell summary rows to satellite rows from the same image/cell

    It only adds columns. It does not remove or transform any quantitative data.
    """
    df_out = df.copy() if df is not None else pd.DataFrame()
    cell_out = cell_summary.copy() if cell_summary is not None else None

    if cell_out is None or cell_out.empty:
        return df_out, cell_out, image_summary, frame_summary

    metrics = [
        ("Satellite_count", "satellite count"),
        ("Mean_distance", "mean distance to centrosome"),
        ("Mean_size", "mean satellite size"),
        ("Total_size", "total satellite size"),
        ("Mean_intensity", "mean satellite intensity"),
        ("Total_intensity", "total satellite intensity"),
        ("Cell_measure", "cell area/volume"),
        ("Centrosome_number_per_cell", "centrosome number"),
    ]
    metrics = [(c, label) for c, label in metrics if c in cell_out.columns and _satq_numeric_series(cell_out, c).notna().sum() >= 5]

    cell_out["Outlier_Cell"] = False
    cell_out["Outlier_Score"] = np.nan
    cell_out["Outlier_Reason"] = ""
    cell_out["Outlier_Label"] = ""

    if metrics:
        score = pd.Series(0.0, index=cell_out.index, dtype=float)
        reasons = {idx: [] for idx in cell_out.index}

        for col, readable in metrics:
            vals = _satq_numeric_series(cell_out, col)
            rz = _satq_robust_zscore(vals)
            low, high = _satq_iqr_bounds(vals, k=1.5)

            for idx in cell_out.index:
                v = vals.loc[idx]
                z = rz.loc[idx] if idx in rz.index else np.nan
                if pd.isna(v):
                    continue

                if pd.notna(z) and np.isfinite(z):
                    score.loc[idx] = max(float(score.loc[idx]), abs(float(z)))

                high_flag = False
                low_flag = False
                if pd.notna(z) and np.isfinite(z) and abs(float(z)) >= 3.5:
                    high_flag = float(z) > 0
                    low_flag = float(z) < 0
                if np.isfinite(high) and v > high:
                    high_flag = True
                if np.isfinite(low) and v < low:
                    low_flag = True

                if high_flag:
                    if pd.notna(z) and np.isfinite(z):
                        reasons[idx].append(f"high {readable} (z={float(z):.2f})")
                    else:
                        reasons[idx].append(f"high {readable}")
                elif low_flag:
                    if pd.notna(z) and np.isfinite(z):
                        reasons[idx].append(f"low {readable} (z={float(z):.2f})")
                    else:
                        reasons[idx].append(f"low {readable}")

        cell_out["Outlier_Score"] = score.replace(0, np.nan)
        cell_out["Outlier_Reason"] = ["; ".join(reasons[idx]) for idx in cell_out.index]
        cell_out["Outlier_Cell"] = cell_out["Outlier_Reason"].astype(str).str.strip() != ""

    def _make_label(row) -> str:
        if "Unique_Cell" in row.index and str(row.get("Unique_Cell", "")).strip():
            return str(row.get("Unique_Cell"))
        src = str(row.get("Source_Input", "")) if "Source_Input" in row.index else ""
        cell = str(row.get("Cell_ID", "")) if "Cell_ID" in row.index else ""
        parts = []
        for c in ["Condition", "Replicate"]:
            if c in row.index and str(row.get(c, "")).strip() and str(row.get(c, "")).lower() not in ["nan", "none"]:
                parts.append(str(row.get(c)))
        if src and cell:
            parts.append(Path(src).stem + " | Cell " + cell)
        elif cell:
            parts.append("Cell " + cell)
        return " / ".join([p for p in parts if p]) or "Outlier cell"

    cell_out["Outlier_Label"] = cell_out.apply(_make_label, axis=1)

    # Propagate cell-level outlier flags to satellite-level rows.
    if not df_out.empty and "Std_Cell_ID" in df_out.columns:
        map_cols = [c for c in ["Source_UID", "Source_Input", "Cell_ID", "Outlier_Cell", "Outlier_Score", "Outlier_Reason", "Outlier_Label"] if c in cell_out.columns]
        mapper = cell_out[map_cols].copy()
        if "Cell_ID" in mapper.columns:
            mapper["__cell_key"] = mapper["Cell_ID"].astype(str)
        df_out["__cell_key"] = df_out["Std_Cell_ID"].astype(str)

        if "Source_UID" in mapper.columns and "Source_UID" in df_out.columns:
            merge_keys = ["Source_UID", "__cell_key"]
        elif "Source_Input" in mapper.columns and "Source_Input" in df_out.columns:
            merge_keys = ["Source_Input", "__cell_key"]
        else:
            merge_keys = ["__cell_key"]

        mapper = mapper.drop_duplicates(subset=merge_keys)
        df_out = df_out.merge(
            mapper[[*merge_keys, "Outlier_Cell", "Outlier_Score", "Outlier_Reason", "Outlier_Label"]],
            on=merge_keys,
            how="left",
        )
        df_out = df_out.drop(columns=["__cell_key"], errors="ignore")
        df_out["Outlier_Cell"] = df_out["Outlier_Cell"].fillna(False)
        df_out["Outlier_Reason"] = df_out["Outlier_Reason"].fillna("")
        df_out["Outlier_Label"] = df_out["Outlier_Label"].fillna("")
    elif not df_out.empty:
        df_out["Outlier_Cell"] = False
        df_out["Outlier_Score"] = np.nan
        df_out["Outlier_Reason"] = ""
        df_out["Outlier_Label"] = ""

    return df_out, cell_out, image_summary, frame_summary


def make_interactive_html(df: pd.DataFrame, image_summary: pd.DataFrame, cell_summary: pd.DataFrame, frame_summary: pd.DataFrame, outdir: Path) -> Optional[Path]:
    """Create a self-contained Plotly HTML graph explorer.

    SatelliteQ creates
    the interactive HTML:
        - clear analysis levels
        - unique cell IDs in combined reports
        - categorical x-axis for cell/image IDs
        - point plots for satellite-level and per-cell summaries
    """
    try:
        import json
        import plotly.offline as po
    except Exception:
        return None

    # Final viewer-level protection: a non-ciliated cell must not appear as a
    # biological 0-µm cilium in any Graph Explorer plot.
    df = sanitize_cilium_length_metrics(df)
    cell_summary = sanitize_cilium_length_metrics(cell_summary)

    def _short_source(value: object, max_len: int = 34) -> str:
        s = "" if value is None else str(value)
        for suffix in [
            "_SatelliteQuantify_Fiji_Output.csv",
            "_timelapse_satellite_rows_combined.csv",
            "_SatelliteQ_Python_Report.xlsx",
            "satellite_3D_results_clean.csv",
            ".csv", ".xlsx", ".xlsm"
        ]:
            s = s.replace(suffix, "")
        return (s[:max_len - 3] + "...") if len(s) > max_len else s

    def _add_display_columns(table: pd.DataFrame) -> pd.DataFrame:
        t = table.copy()

        if "Source_Input" in t.columns:
            base_source = t["Source_Input"].map(_short_source)
        else:
            base_source = pd.Series(["Input"] * len(t), index=t.index)

        def _clean_label_value(v):
            if v is None:
                return ""
            s = str(v)
            if s.lower() in ["nan", "none", ""]:
                return ""
            return s

        def _combine_source_label(row):
            cond = _clean_label_value(row.get("Condition", ""))
            rep = _clean_label_value(row.get("Replicate", ""))
            src = _clean_label_value(row.get("_base_source_for_label", "Input"))
            parts = []
            for part in [cond, rep, src]:
                if part and part not in parts:
                    parts.append(part)
            return " / ".join(parts) if parts else "Input"

        t["_base_source_for_label"] = base_source
        t["Source_short"] = t.apply(_combine_source_label, axis=1)
        t = t.drop(columns=["_base_source_for_label"], errors="ignore")

        cell_col = None
        for c in ["Std_Cell_ID", "Cell_ID"]:
            if c in t.columns:
                cell_col = c
                break
        if cell_col is not None:
            cell_values = t[cell_col].replace("", np.nan)
            t["Unique_Cell"] = np.where(
                cell_values.notna(),
                t["Source_short"].astype(str) + " | Cell " + cell_values.astype(str),
                ""
            )
            t["Plot_Cell"] = t["Unique_Cell"]
        else:
            t["Unique_Cell"] = ""
            t["Plot_Cell"] = ""

        sat_col = None
        for c in ["Std_Satellite_ID", "Satellite_ID", "Object", "ID"]:
            if c in t.columns:
                sat_col = c
                break
        if sat_col is not None:
            sat_values = t[sat_col].replace("", np.nan)
            t["Unique_Satellite"] = np.where(
                sat_values.notna(),
                t["Source_short"].astype(str) + " | Sat " + sat_values.astype(str),
                ""
            )
        else:
            t["Unique_Satellite"] = ""

        frame_col = None
        for c in ["Std_Frame", "Frame"]:
            if c in t.columns:
                frame_col = c
                break
        if frame_col is not None:
            frame_values = t[frame_col].replace("", np.nan)
            t["Unique_Frame"] = np.where(
                frame_values.notna(),
                t["Source_short"].astype(str) + " | Frame " + frame_values.astype(str),
                ""
            )
        else:
            t["Unique_Frame"] = ""

        t["Plot_Image"] = t["Source_short"]
        return t

    def _visible_columns(table_name: str, t: pd.DataFrame) -> List[str]:
        if table_name.startswith("Satellite-level"):
            preferred = [
                "Condition", "Replicate", "Plot_Cell", "Plot_Image", "Unique_Satellite",
                "Std_Area_2D", "Std_Volume_3D", "Std_Size", "Std_Surface", "Std_Sphericity", "Std_Sphericity_Source", "Std_Circularity_2D", "Std_Shape_Compactness",
                "Std_Total_Intensity", "Std_Intensity", "Std_Mean_Intensity", "Std_Distance",
                "Std_Cell_Spatial_Entropy", "Std_Cell_Pericentrosomal_Clustering_Index", "Std_Distance_to_Cell_or_Image_Center_XY",
                "Std_Center_Reference_X", "Std_Center_Reference_Y", "Std_Center_Reference_Type", "Std_Center_Distance_Unit",
                "Std_Cell_Area_2D", "Std_Cell_Volume_3D", "Std_Cell_Measure", "Std_X", "Std_Y", "Std_Z",
                "Std_Frame", "Std_Time_min", "Mode", "Analysis_Type"
            ]
        elif table_name.startswith("Per-cell"):
            preferred = [
                "Condition", "Replicate", "Plot_Cell", "Unique_Cell", "Plot_Image", "Cell_ID", "Mode",
                "Satellite_count", "Centrosome_number_per_cell",
                "Mean_area_2D", "Median_area_2D", "Total_area_2D",
                "Mean_volume_3D", "Median_volume_3D", "Total_volume_3D",
                "Mean_size", "Median_size", "Total_size",
                "Mean_intensity", "Total_intensity", "Mean_distance", "Mean_sphericity", "Mean_circularity_2D", "Mean_shape_compactness",
                "Cell_spatial_entropy",
                "Cell_pericentrosomal_clustering_index",
                PCI_CELL_INNER_COL, PCI_CELL_OUTER_COL, PCI_CELL_TOTAL_COL,
                "Mean_distance_to_cell_or_image_center",
                "Cell_area_2D", "Cell_volume_3D", "Cell_measure", "Cell area", "Cell volume"
            ]
        elif table_name.startswith("Per-image"):
            preferred = [
                "Condition", "Replicate", "Plot_Image", "Source_Input", "Source_Folder", "Mode",
                "Satellite_count", "Assigned_satellites", "Cells",
                "Mean_area_2D", "Median_area_2D", "Total_area_2D",
                "Mean_volume_3D", "Median_volume_3D", "Total_volume_3D",
                "Mean_size", "Median_size", "Total_size",
                "Mean_intensity", "Total_intensity", "Mean_distance", "Mean_sphericity", "Mean_circularity_2D", "Mean_shape_compactness",
                "Mean_cell_spatial_entropy", "Mean_cell_pericentrosomal_clustering_index", "Mean_distance_to_cell_or_image_center"
            ]
        else:
            preferred = [
                "Condition", "Replicate", "Unique_Frame", "Plot_Image", "Frame", "Time_min",
                "Satellite_count",
                "Mean_area_2D", "Median_area_2D", "Total_area_2D",
                "Mean_volume_3D", "Median_volume_3D", "Total_volume_3D",
                "Mean_size", "Median_size", "Total_size", "Mean_intensity", "Total_intensity",
                "Mean_distance", "Mean_sphericity", "Mean_circularity_2D", "Mean_shape_compactness",
                "Mean_cell_spatial_entropy", "Mean_cell_pericentrosomal_clustering_index", "Mean_distance_to_cell_or_image_center",
                "Mean_cell_area_2D", "Total_cell_area_2D", "Mean_cell_volume_3D", "Total_cell_volume_3D",
                "Mean_cell_size", "Total_cell_size"
            ]

        cols = [c for c in preferred if c in t.columns]
        for oc in ["Outlier_Cell", "Outlier_Score", "Outlier_Reason", "Outlier_Label"]:
            if oc in t.columns and oc not in cols:
                cols.append(oc)

        metadata_like = {
            "Source_Input", "Source_File", "Source_Path", "Source_UID", "Source_Folder", "Source_short",
            "Run_ID", "Image", "Mode", "Condition", "Replicate", "Cell_ID", "Unique_Cell", "Plot_Cell", "Plot_Image",
            "Threshold_method", "Threshold_scope", "Gaussian_smoothing_enabled"
        }
        for c in t.columns:
            if c in cols:
                continue
            if c.startswith("Std_") or c in ["Area", "Volume", "Integrated_Intensity", "IntDen", "Mean", "Distance_to_nearest_centrosome_3D"]:
                cols.append(c)
                continue
            if table_name.startswith("Per-cell"):
                if c not in metadata_like and pd.to_numeric(t[c], errors="coerce").notna().sum() > 0:
                    cols.append(c)

        return cols if cols else list(t.columns)

    df, cell_summary, image_summary, frame_summary = add_satelliteq_outlier_flags(df, cell_summary, image_summary, frame_summary)

    tables = {
        "Satellite-level data (one row = one satellite)": df,
    }
    if cell_summary is not None and not cell_summary.empty:
        tables["Per-cell summary (one row = one cell)"] = cell_summary
    if frame_summary is not None and not frame_summary.empty:
        tables["Time-lapse frame summary (one row = one frame)"] = frame_summary

    prepared = {}
    for name, table in tables.items():
        if table is None or table.empty:
            continue

        t = _add_display_columns(table)

        if len(t) > 50000:
            t = t.sample(50000, random_state=1).copy()

        cols = _visible_columns(name, t)
        t = t[cols].copy()
        t = t.replace({np.nan: None})

        cols = list(t.columns)
        numeric_cols = []
        for c in cols:
            vals = pd.to_numeric(t[c], errors="coerce")
            if vals.notna().sum() > 0:
                numeric_cols.append(c)

        if name.startswith("Per-cell"):
            for c in PERICENTROSOMAL_CELL_METRICS:
                if c in cols and c not in numeric_cols:
                    numeric_cols.append(c)

        # IMPORTANT: JavaScript Number(null) == 0. For cilium-length metrics,
        # missing means "no cilium", not a biological zero-length cilium.
        # Encode missing/non-positive length values as a nonnumeric token in the
        # interactive Graph Explorer so they are excluded rather than plotted at 0.
        for c in [
            "Cilium_length_um", "Mean_cilium_length_um",
            "Std_Cilium_Length_um", "Std_Cilium_Mean_Length_um",
            "Cilium_Length_um", "Mean_Cilium_Length_um",
        ]:
            if c in t.columns:
                vals = pd.to_numeric(t[c], errors="coerce")
                invalid = vals.isna() | (vals <= 0)
                t[c] = t[c].astype(object)
                t.loc[invalid, c] = "NA"

        categorical_axis_cols = []
        for c in cols:
            lc = c.lower()
            if (
                c not in numeric_cols
                or "id" in lc
                or c.startswith("Plot_")
                or c.startswith("Unique_")
                or c in ["Source_Input", "Source_short", "Cell_ID", "Mode", "Analysis_Type"]
            ):
                categorical_axis_cols.append(c)

        labels = {c: label_for(c, table) for c in cols}

        if name.startswith("Satellite-level"):
            default_x = "Condition" if "Condition" in cols else ("Plot_Cell" if "Plot_Cell" in cols else ("Plot_Image" if "Plot_Image" in cols else cols[0]))
            default_y = (
                "Std_Area_2D" if "Std_Area_2D" in cols and pd.to_numeric(t["Std_Area_2D"], errors="coerce").notna().any()
                else "Std_Volume_3D" if "Std_Volume_3D" in cols and pd.to_numeric(t["Std_Volume_3D"], errors="coerce").notna().any()
                else "Std_Size" if "Std_Size" in cols
                else (numeric_cols[0] if numeric_cols else cols[0])
            )
            default_type = "point"
        elif name.startswith("Per-cell"):
            if any(("Pearson" in str(c) or "Manders" in str(c) or str(c) == "Overlap_R") for c in cols):
                default_x = "Condition" if "Condition" in cols else (cols[0] if cols else "")
                for candidate in ["Pearson_R_recommended_for_puncta", "Pearson_R_gaussian_smoothed", "Overlap_R", "Pearson_R_raw_all_pixels", "Manders_M1_A_in_B"]:
                    if candidate in cols:
                        default_y = candidate
                        break
                else:
                    default_y = numeric_cols[0] if numeric_cols else (cols[0] if cols else "")
            else:
                default_x = "Satellite_count" if "Satellite_count" in cols else (numeric_cols[0] if numeric_cols else cols[0])
                default_y = "Mean_size" if "Mean_size" in cols else ("Mean_intensity" if "Mean_intensity" in cols else (numeric_cols[1] if len(numeric_cols) > 1 else (numeric_cols[0] if numeric_cols else cols[0])))
            default_type = "point"
        elif name.startswith("Per-image"):
            default_x = "Plot_Image" if "Plot_Image" in cols else cols[0]
            default_y = "Satellite_count" if "Satellite_count" in cols else (numeric_cols[0] if numeric_cols else cols[0])
            default_type = "bar"
        else:
            default_x = "Time_min" if "Time_min" in cols else ("Frame" if "Frame" in cols else cols[0])
            default_y = "Satellite_count" if "Satellite_count" in cols else (numeric_cols[0] if numeric_cols else cols[0])
            default_type = "line"

        prepared[name] = {
            "rows": t.to_dict(orient="records"),
            "columns": cols,
            "numeric": numeric_cols,
            "categoricalAxis": categorical_axis_cols,
            "labels": labels,
            "defaultX": default_x,
            "defaultY": default_y,
            "defaultType": default_type,
        }

    if not prepared:
        return None

    plotlyjs = po.get_plotlyjs()
    data_json = json.dumps(prepared, ensure_ascii=False)

    html_template = """<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>SatelliteQ Interactive Graph Explorer</title>
<style>
body { font-family: Arial, sans-serif; margin: 0 auto; padding: 18px 22px; max-width: 1680px; background: #fbf8ff; color: #222; }
h1 { margin: 0 0 10px 0; color: #7E3FB5; font-size: 22px; }
.controls { display: grid; grid-template-columns: repeat(4, minmax(170px, 1fr)); gap: 10px; background: white; border: 1px solid #DCC6FF; padding: 12px; border-radius: 10px; box-shadow: 0 2px 10px rgba(126,63,181,0.08); }
label { font-weight: bold; font-size: 13px; display: block; margin-bottom: 4px; color: #4B216E; }
select { width: 100%; padding: 6px; }
.celltools { margin-top: 12px; background: #ffffff; border: 1px solid #DCC6FF; padding: 12px; border-radius: 10px; box-shadow: 0 2px 10px rgba(126,63,181,0.08); }
.celltools_header { display:flex; justify-content:space-between; gap:12px; align-items:flex-start; margin-bottom:10px; }
.celltools .hint { font-size: 12px; color: #555; margin-top: 4px; line-height: 1.3; }
.cellbuttons { display:flex; gap:8px; flex-wrap:wrap; justify-content:flex-end; }
.cellchips { display:flex; gap:8px; flex-wrap:wrap; max-height:140px; overflow:auto; padding:6px; background:#fbf8ff; border:1px solid #eadcff; border-radius:8px; }
.cellchip { display:inline-flex; align-items:center; gap:5px; border:1px solid #DCC6FF; background:#F6EFFF; border-radius:999px; padding:5px 9px; cursor:pointer; color:#4B216E; font-size:13px; }
.cellchip input { margin:0; }
.stattools { margin-top:10px; display:flex; gap:18px; flex-wrap:wrap; align-items:flex-end; color:#4B216E; }
.stattools label { font-weight:normal; color:#4B216E; margin:0; }
.stattools select { min-width:260px; }
.statbar { margin-top: 10px; display:flex; gap:16px; flex-wrap:wrap; align-items:center; background:white; border:1px solid #DCC6FF; padding:10px 12px; border-radius:10px; box-shadow:0 2px 10px rgba(126,63,181,0.08); color:#4B216E; }
.statbar label { display:flex; align-items:center; gap:6px; margin:0; font-weight:bold; color:#4B216E; }
.statbar select { min-width:220px; padding:6px; }
button { padding: 7px 10px; border-radius: 7px; border: 1px solid #7E3FB5; background: #F1E7FF; color: #4B216E; cursor: pointer; }
button:hover { background: #E6D4FF; }
#plot { width: 100%; max-width: none; background: white; border: 1px solid #DCC6FF; border-radius: 10px; margin: 14px auto 0 auto; padding: 8px; box-shadow: 0 2px 10px rgba(126,63,181,0.08); }
.graphwrap { display:flex; gap:12px; align-items:flex-start; }
.graphwrap #plot { flex:1 1 auto; min-width:0; }
.filterpanel { width:240px; flex:0 0 240px; margin-top:14px; background:white; border:1px solid #DCC6FF; border-radius:10px; padding:10px; box-shadow:0 2px 10px rgba(126,63,181,0.08); }
.filterpanel h3 { margin:0 0 8px 0; color:#7E3FB5; font-size:15px; }
.filtersection { margin-bottom:12px; }
.filtertitle { font-weight:bold; color:#4B216E; font-size:13px; margin-bottom:5px; }
.filterlist { max-height:150px; overflow:auto; border:1px solid #eadcff; border-radius:8px; padding:6px; background:#fbf8ff; }
.filteritem { display:flex; align-items:center; gap:6px; font-size:12px; color:#333; margin:4px 0; cursor:pointer; }
.filterbuttons { display:flex; gap:6px; margin-top:6px; }
.outlierlist { max-height:520px; overflow:auto; font-size:12px; line-height:1.35; color:#333; }
.outlieritem { border-left:4px solid #D62728; padding:7px 8px; margin:7px 0; background:#fff6f6; border-radius:7px; }
.outlieritem b { color:#A00000; }
.outlierhint { color:#555; font-size:12px; margin-bottom:8px; }
.outlierbottom { margin-top:12px; background:white; border:1px solid #DCC6FF; border-radius:10px; padding:12px; box-shadow:0 2px 10px rgba(126,63,181,0.08); }
.outlierbottom h3 { margin:0 0 8px 0; color:#7E3FB5; font-size:15px; }
.outlierbottom .outlierlist { max-height:270px; overflow:auto; display:grid; grid-template-columns:repeat(auto-fit, minmax(260px, 1fr)); gap:8px; }
.outlierbottom .outlieritem { margin:0; }
.legendbottom { margin-top:12px; background:white; border:1px solid #DCC6FF; border-radius:10px; padding:12px; box-shadow:0 2px 10px rgba(126,63,181,0.08); }
.legendbottom h3 { margin:0 0 8px 0; color:#7E3FB5; font-size:15px; }
.legendhint { color:#555; font-size:12px; margin-bottom:8px; }
.legenditems { display:grid; grid-template-columns:repeat(auto-fit, minmax(260px, 1fr)); gap:8px; }
.legenditem { display:flex; align-items:flex-start; gap:8px; padding:7px 8px; background:#fbf8ff; border:1px solid #eadcff; border-radius:7px; font-size:12px; line-height:1.25; min-width:0; }
.legenddot { width:13px; height:13px; border-radius:50%; margin-top:2px; flex:0 0 13px; border:1px solid rgba(0,0,0,0.25); }
.legendtext { overflow-wrap:anywhere; word-break:break-word; }

.filterbuttons button { font-size:11px; padding:4px 7px; }
.presetbar { margin-top: 12px; display: flex; flex-wrap: wrap; gap: 8px; }
.satstats { background:white; border:1px solid #DCC6FF; border-radius:10px; padding:9px 12px; box-shadow:0 2px 10px rgba(126,63,181,0.08); display:none; align-items:center; gap:10px; }
.satnormalize { background:white; border:1px solid #DCC6FF; border-radius:10px; padding:9px 12px; box-shadow:0 2px 10px rgba(126,63,181,0.08); display:none; align-items:center; gap:8px; }
.satnormalize label { margin:0; font-weight:bold; color:#4B216E; }
.satnormalize select { min-width:150px; }
.satstats label { margin:0; display:flex; align-items:center; gap:6px; font-weight:bold; color:#4B216E; }
.namedisplay { background:white; border:1px solid #DCC6FF; border-radius:10px; padding:9px 12px; box-shadow:0 2px 10px rgba(126,63,181,0.08); display:flex; align-items:center; gap:8px; }
.namedisplay label { margin:0; font-weight:bold; color:#4B216E; }
.namedisplay select { min-width:140px; }
</style>
<script type="text/javascript">__PLOTLYJS__</script>
</head>
<body>
<h1>SatelliteQ Interactive Graph Explorer</h1>
<div class="hint">SatelliteQ version __SATELLITEQ_VERSION__</div>
<div class="controls">
  <div id="tableControl"><label>Analysis level</label><select id="tableSelect"></select></div>
  <div id="graphTypeControl"><label>Graph type</label><select id="typeSelect"><option value="point">point plot</option><option value="cell_panels">cell panels</option><option value="bar">bar</option><option value="box">box</option><option value="line">line</option><option value="histogram">histogram</option></select></div>
  <div id="xAxisControl"><label>X-axis</label><select id="xSelect"></select></div>
  <div id="yAxisControl"><label>Y-axis</label><select id="ySelect"></select></div>
  <div id="groupControl"><label>Group / color</label><select id="groupSelect"></select></div>
  <div id="nameDisplayControl" class="namedisplay"><label>Long names</label><select id="nameDisplayMode"><option value="short" selected>Short</option><option value="wrapped">Wrapped</option><option value="full">Full</option></select></div>
  <div id="satelliteStatsControl" class="satstats">
    <label><input type="checkbox" id="satelliteShowMeanSD" checked> Show mean ± SD and n</label>
    <label><input type="checkbox" id="satelliteShowOutliers"> Show outlier cells/points</label>
  </div>
  <div id="satelliteNormalizeControl" class="satnormalize"><label>Normalize to</label><select id="satelliteNormalizeCondition"></select></div>
  <div id="perCellStatsControl" class="satstats">
    <label><input type="checkbox" id="perCellShowMeanSD" checked> Show mean ± SD and n</label>
    <label><input type="checkbox" id="perCellShowOutliers"> Show outlier cells/points</label>
    <label><input type="checkbox" id="perCellShowTrendline" checked> Show trendline</label>
    <label><input type="checkbox" id="perCellShowRValue" checked> Show R value</label>
  </div>
</div>


<div class="celltools" id="cellPanelTools" style="display:none;">
  <div class="celltools_header">
    <div>
      <label>Select cells for cell panels</label>
      <div class="hint">Click cells normally. No Ctrl/Command selection needed. Each selected cell becomes its own small graph.</div>
    </div>
    <div class="cellbuttons">
      <button onclick="selectAllCells()">Select all cells</button>
      <button onclick="clearSelectedCells()">Clear cells</button>
    </div>
  </div>
  <div id="cellPanelSelect" class="cellchips"></div>
  <div class="stattools">
    <div>
      <label>Per-cell plot style</label>
      <select id="cellPlotType">
        <option value="point">Points: one point per satellite</option>
        <option value="profile">Plot profile: line + points</option>
        <option value="bar">Bar plot: one bar per satellite</option>
        <option value="box">Box plot: distribution per cell</option>
        <option value="scatterxy">X vs Y scatter per cell</option>
      </select>
    </div>
    <div>
      <label>Panels per row</label>
      <select id="panelsPerRow">
        <option value="2">2</option>
        <option value="3" selected>3</option>
        <option value="4">4</option>
        <option value="5">5</option>
      </select>
    </div>
    <label><input type="checkbox" id="showMean" checked> Show mean</label>
    <label><input type="checkbox" id="showSD"> Show mean ± SD</label>
  </div>
</div>

<div class="presetbar" id="presetbar" style="display:none;"></div>

<div class="graphwrap">
  <div id="plot" style="height:620px;"></div>
  <aside id="satelliteFilterPanel" class="filterpanel" style="display:none;">
    <h3 id="filterPanelTitle">Filters</h3>
    <div class="filtersection">
      <div class="filtertitle">Condition</div>
      <div id="conditionFilters" class="filterlist"></div>
      <div class="filterbuttons"><button onclick="setAllFilter('conditionFilters', true)">All</button><button onclick="setAllFilter('conditionFilters', false)">None</button></div>
    </div>
    <div class="filtersection">
      <div class="filtertitle">Replicate</div>
      <div id="replicateFilters" class="filterlist"></div>
      <div class="filterbuttons"><button onclick="setAllFilter('replicateFilters', true)">All</button><button onclick="setAllFilter('replicateFilters', false)">None</button></div>
    </div>
    <div class="filtersection" id="frameFilterSection" style="display:none;">
      <div class="filtertitle">Frame</div>
      <div id="frameFilters" class="filterlist"></div>
      <div class="filterbuttons"><button onclick="setAllFilter('frameFilters', true)">All</button><button onclick="setAllFilter('frameFilters', false)">None</button></div>
    </div>
  </aside>
</div>
<section id="legendPanel" class="legendbottom" style="display:none;">
  <h3>Graph legend</h3>
  <div id="legendSummary"></div>
</section>
<section id="outlierPanel" class="outlierbottom" style="display:none;">
  <h3>Outlier cells</h3>
  <div id="outlierSummary" class="outlierlist"></div>
</section>

<script>
const DATA = __DATA_JSON__;
const SATQ_PURPLE = '#7E3FB5';
const SATQ_LIGHT_PURPLE = '#DCC6FF';
const SATQ_PURPLE_DARK = '#4B216E';
const SATQ_PURPLE_PALETTE = ['#7E3FB5','#A66BDA','#C39BEA','#6B4FA3','#B57EDC','#8E5CC2','#D6B7FF','#5D3A91'];

const tableSelect = document.getElementById('tableSelect');
const typeSelect = document.getElementById('typeSelect');
const xSelect = document.getElementById('xSelect');
const ySelect = document.getElementById('ySelect');
const groupSelect = document.getElementById('groupSelect');
const graphTypeControl = document.getElementById('graphTypeControl');
const xAxisControl = document.getElementById('xAxisControl');
const groupControl = document.getElementById('groupControl');
const nameDisplayMode = document.getElementById('nameDisplayMode');
const presetbar = document.getElementById('presetbar');
const satelliteFilterPanel = document.getElementById('satelliteFilterPanel');
const conditionFilters = document.getElementById('conditionFilters');
const replicateFilters = document.getElementById('replicateFilters');
const frameFilters = document.getElementById('frameFilters');
const frameFilterSection = document.getElementById('frameFilterSection');
const satelliteStatsControl = document.getElementById('satelliteStatsControl');
const satelliteShowMeanSD = document.getElementById('satelliteShowMeanSD');
const satelliteShowOutliers = document.getElementById('satelliteShowOutliers');
const satelliteNormalizeControl = document.getElementById('satelliteNormalizeControl');
const satelliteNormalizeCondition = document.getElementById('satelliteNormalizeCondition');
const perCellStatsControl = document.getElementById('perCellStatsControl');
const perCellShowMeanSD = document.getElementById('perCellShowMeanSD');
const perCellShowOutliers = document.getElementById('perCellShowOutliers');
const perCellShowTrendline = document.getElementById('perCellShowTrendline');
const perCellShowRValue = document.getElementById('perCellShowRValue');
const outlierPanel = document.getElementById('outlierPanel');
const outlierSummary = document.getElementById('outlierSummary');
const legendPanel = document.getElementById('legendPanel');
const legendSummary = document.getElementById('legendSummary');
const filterPanelTitle = document.getElementById('filterPanelTitle');
const cellPanelTools = document.getElementById('cellPanelTools');
const cellPanelSelect = document.getElementById('cellPanelSelect');
const cellPlotType = document.getElementById('cellPlotType');
const panelsPerRow = document.getElementById('panelsPerRow');
const showMean = document.getElementById('showMean');
const showSD = document.getElementById('showSD');

function label(table, col) { return DATA[table].labels[col] || col; }
function isNumericCol(table, col) { return DATA[table].numeric.includes(col); }
function isCategoricalAxis(table, col) { return DATA[table].categoricalAxis.includes(col); }

function cleanDisplayName(value) {
  let s = String(value === null || value === undefined ? '' : value);
  const suffixes = [
    '_SatelliteQuantify_Fiji_Output.csv',
    '_timelapse_satellite_rows_combined.csv',
    '_SatelliteQ_Python_Report.xlsx',
    'satellite_3D_results_clean.csv',
    'satellite_3D_results_clean',
    '_SatelliteQuantify_3D_satellites_only',
    'SatelliteQuantify_3D_satellites_only',
    '.csv', '.xlsx', '.xlsm'
  ];
  suffixes.forEach(suf => { s = s.split(suf).join(''); });
  // Remove common date/time stamps from display labels only. Original labels remain in hover/data.
  s = s.replace(/(?:^|[_\\s-])(?:20\\d{2}|19\\d{2})[_\\s-]?\\d{1,2}[_\\s-]?\\d{1,2}[_\\s-]?\\d{1,2}[_\\s-]?\\d{1,2}[_\\s-]?\\d{1,2}(?=$|[_\\s-])/g, '_');
  s = s.replace(/_{2,}/g, '_').replace(/\\s{2,}/g, ' ').trim();
  return s || String(value || '');
}

function shortenDisplayPart(value, maxLen=28) {
  let s = cleanDisplayName(value);
  if (s.length <= maxLen) return s;
  const head = Math.max(8, Math.floor(maxLen * 0.58));
  const tail = Math.max(6, maxLen - head - 1);
  return s.slice(0, head) + '…' + s.slice(s.length - tail);
}

function wrapOneDisplayLine(line, maxLen=30) {
  const s = String(line || '');
  if (s.length <= maxLen) return [s];
  const tokens = s.split(/([_\\s-]+)/);
  const out = [];
  let cur = '';
  tokens.forEach(tok => {
    if (!tok) return;
    if ((cur + tok).length <= maxLen) {
      cur += tok;
    } else {
      if (cur.trim()) out.push(cur.trim());
      if (tok.length > maxLen) {
        for (let i=0; i<tok.length; i+=maxLen) out.push(tok.slice(i, i+maxLen));
        cur = '';
      } else {
        cur = tok;
      }
    }
  });
  if (cur.trim()) out.push(cur.trim());
  return out.length ? out : [s];
}

function wrapDisplayName(value, maxLen=30, maxLines=5) {
  let s = cleanDisplayName(value);
  const roughLines = s.split(/\\s+\\/\\s+|\\s+—\\s+/);
  let lines = [];
  roughLines.forEach(part => { lines = lines.concat(wrapOneDisplayLine(part, maxLen)); });
  if (lines.length > maxLines) {
    lines = lines.slice(0, maxLines);
    lines[maxLines - 1] = lines[maxLines - 1].replace(/…?$/, '') + '…';
  }
  return lines.join('<br>');
}

function displayName(value, mode=null) {
  const m = mode || (nameDisplayMode ? nameDisplayMode.value : 'short');
  if (m === 'full') return String(value === null || value === undefined ? '' : value);
  if (m === 'wrapped') return wrapDisplayName(value, 30, 5);
  const s = cleanDisplayName(value);
  const parts = s.split(/\\s+\\/\\s+/);
  if (parts.length > 1) return parts.map(p => shortenDisplayPart(p, 24)).join(' / ');
  return shortenDisplayPart(s, 46);
}

function displayTickLabel(value) {
  const m = nameDisplayMode ? nameDisplayMode.value : 'short';
  if (m === 'full') return String(value === null || value === undefined ? '' : value);
  if (m === 'wrapped') return wrapDisplayName(value, 22, 4);
  return shortenDisplayPart(value, 28);
}

function patchAxisTickLabels(axisObj) {
  if (!axisObj || !Array.isArray(axisObj.ticktext)) return axisObj;
  const copy = Object.assign({}, axisObj);
  copy.ticktext = axisObj.ticktext.map(displayTickLabel);
  return copy;
}

function patchLayoutDisplayNames(layout) {
  if (!layout || !nameDisplayMode || nameDisplayMode.value === 'full') return layout;
  const out = Object.assign({}, layout);
  ['xaxis','yaxis','xaxis2','yaxis2','xaxis3','yaxis3','xaxis4','yaxis4','xaxis5','yaxis5','xaxis6','yaxis6'].forEach(k => {
    if (out[k]) out[k] = patchAxisTickLabels(out[k]);
  });
  return out;
}

function patchTraceDisplayNames(traces) {
  if (!Array.isArray(traces) || !nameDisplayMode || nameDisplayMode.value === 'full') return traces;
  return traces.map(tr => {
    if (!tr || tr.name === undefined || tr.name === null) return tr;
    const t = Object.assign({}, tr);
    t.name = displayName(tr.name);
    return t;
  });
}


function legendEscapeHtml(value) {
  return String(value === null || value === undefined ? '' : value)
    .replace(/&/g, '&amp;')
    .replace(/</g, '&lt;')
    .replace(/>/g, '&gt;')
    .replace(/"/g, '&quot;');
}
function legendTraceColor(trace, index) {
  let c = null;
  if (trace && trace.marker && trace.marker.color !== undefined) c = trace.marker.color;
  if ((c === null || c === undefined) && trace && trace.line && trace.line.color !== undefined) c = trace.line.color;
  if (Array.isArray(c)) c = c.find(v => typeof v === 'string') || null;
  if (typeof c !== 'string' || c.length === 0) c = SATQ_PURPLE_PALETTE[index % SATQ_PURPLE_PALETTE.length] || SATQ_PURPLE;
  return c;
}
function updateCustomLegend(traces, layout) {
  if (!legendPanel || !legendSummary) return;
  const entries = [];
  const seen = new Set();
  (traces || []).forEach((tr, i) => {
    if (!tr || tr.showlegend === false) return;
    if (tr.name === undefined || tr.name === null || String(tr.name).trim() === '') return;
    const name = String(tr.name);
    if (name.toLowerCase() === 'undefined' || name.toLowerCase() === 'null') return;
    const color = legendTraceColor(tr, i);
    const key = name + '|' + color;
    if (seen.has(key)) return;
    seen.add(key);
    entries.push({name:name, color:color});
  });
  if (!entries.length) {
    legendPanel.style.display = 'none';
    legendSummary.innerHTML = '';
    return;
  }
  let title = 'Graph legend';
  try {
    if (layout && layout.legend && layout.legend.title) {
      title = typeof layout.legend.title === 'string' ? layout.legend.title : (layout.legend.title.text || title);
    }
  } catch(e) {}
  legendPanel.style.display = '';
  let html = '<div class="legendhint"><b>' + legendEscapeHtml(title) + '</b> — colors and dot labels for the current graph.</div>';
  html += '<div class="legenditems">';
  entries.forEach(e => {
    html += '<div class="legenditem"><span class="legenddot" style="background:' + legendEscapeHtml(e.color) + '"></span><span class="legendtext">' + legendEscapeHtml(e.name) + '</span></div>';
  });
  html += '</div>';
  legendSummary.innerHTML = html;
}

const SatelliteQOriginalNewPlot = Plotly.newPlot.bind(Plotly);
Plotly.newPlot = function(div, traces, layout, config) {
  const patchedTraces = patchTraceDisplayNames(traces);
  const patchedLayout = patchLayoutDisplayNames(layout);
  updateCustomLegend(patchedTraces, patchedLayout);
  let layoutNoBuiltInLegend = patchedLayout;
  if (legendPanel && patchedLayout) layoutNoBuiltInLegend = Object.assign({}, patchedLayout, {showlegend:false});
  return SatelliteQOriginalNewPlot(div, patchedTraces, layoutNoBuiltInLegend, config);
};

function addOptions(sel, opts, labels=null, includeNone=false) {
  sel.innerHTML = '';
  if (includeNone) { let o=document.createElement('option'); o.value=''; o.textContent='None'; sel.appendChild(o); }
  opts.forEach(c => { let o=document.createElement('option'); o.value=c; o.textContent=labels ? labels[c] || c : c; sel.appendChild(o); });
}

function isSatelliteLevelTable(tableName=null) {
  const t = tableName || tableSelect.value;
  return String(t).startsWith('Satellite-level');
}
function isPerCellTable(tableName=null) {
  const t = tableName || tableSelect.value;
  return String(t).startsWith('Per-cell');
}
function isTimeLapseTable(tableName=null) {
  const t = tableName || tableSelect.value;
  return String(t).startsWith('Time-lapse');
}
function isFilterableLevelTable(tableName=null) {
  return isSatelliteLevelTable(tableName) || isPerCellTable(tableName) || isTimeLapseTable(tableName);
}

function selectedFilterValues(containerId) {
  const el = document.getElementById(containerId);
  if (!el) return [];
  return Array.from(el.querySelectorAll('input[type="checkbox"]:checked')).map(cb => cb.value);
}

function buildFilterCheckboxes(container, valuesArr, previousSelected) {
  container.innerHTML = '';
  const valuesClean = uniqueNonEmpty(valuesArr.map(v => String(v))).sort((a,b) => a.localeCompare(b, undefined, {numeric:true}));
  const prev = new Set(previousSelected || []);
  const usePrev = previousSelected && previousSelected.length > 0;
  valuesClean.forEach(v => {
    const lab = document.createElement('label');
    lab.className = 'filteritem';
    const cb = document.createElement('input');
    cb.type = 'checkbox';
    cb.value = v;
    cb.checked = usePrev ? prev.has(v) : true;
    cb.addEventListener('change', drawPlot);
    const span = document.createElement('span');
    span.textContent = v;
    lab.appendChild(cb);
    lab.appendChild(span);
    container.appendChild(lab);
  });
}

function isNormalizedY(col) {
  return String(col || '').startsWith('NORM::');
}
function rawYColumn(col) {
  return isNormalizedY(col) ? String(col).replace('NORM::', '') : col;
}
function satelliteAllowedYColumns(info) {
  const dropCols = new Set([
    'Std_X','Std_Y','X','Y','X_um','Y_um','XM','YM',
    'Std_Satellite_ID','Satellite_ID','Object','Obj','ID',
    'Std_Cell_ID','Cell_ID',
    'Std_Centrosome_ID','Assigned_Centrosome_ID','Nearest_Centrosome_ID','Centrosome_ID',
    'Std_Centrosome_X','Std_Centrosome_Y','Assigned_Centrosome_X_um','Assigned_Centrosome_Y_um',
    'Nearest_Centrosome_X_um','Nearest_Centrosome_Y_um','Centrosome_X_um','Centrosome_Y_um',
    // Hide raw/duplicate integrated-intensity names; keep the clean Std_Total_Intensity option.
    'Std_Intensity','Integrated_Intensity_AU','Integrated_Intensity','IntDen','Integrated density','Integrated_Density'
  ]);
  return (info.numeric || []).filter(c => {
    const col = String(c);
    const lab = String((info.labels && info.labels[c]) || c).toLowerCase();
    const lc = col.toLowerCase();
    if (dropCols.has(col)) return false;
    if (lc.includes('centrosome_x') || lc.includes('centrosome_y')) return false;
    if (lc === 'std_x' || lc === 'std_y') return false;
    if (lc.includes('satellite_id') || lc.includes('cell_id') || lc.includes('centrosome_id')) return false;
    if (lab.includes('x position') || lab.includes('y position')) return false;
    if (lab.includes('satellite id') || lab.includes('cell id') || lab.includes('centrosome id')) return false;
    if (lab.includes('integrated intensity')) return false;
    return true;
  });
}
function populateSatelliteYOptions(info) {
  const prev = ySelect.value;
  ySelect.innerHTML = '';
  const allowed = satelliteAllowedYColumns(info);
  allowed.forEach(c => {
    let o = document.createElement('option');
    o.value = c;
    o.textContent = label(tableSelect.value, c);
    ySelect.appendChild(o);
  });
  if (allowed.length > 0) {
    const group = document.createElement('optgroup');
    group.label = 'Normalized value';
    allowed.forEach(c => {
      let o = document.createElement('option');
      o.value = 'NORM::' + c;
      o.textContent = label(tableSelect.value, c) + ' — normalized';
      group.appendChild(o);
    });
    ySelect.appendChild(group);
  }
  if ([...ySelect.options].some(o => o.value === prev)) ySelect.value = prev;
  else if (allowed.includes('Std_Area_2D') && (info.rows || []).some(r => r.Std_Area_2D !== null && r.Std_Area_2D !== undefined && r.Std_Area_2D !== '')) ySelect.value = 'Std_Area_2D';
  else if (allowed.includes('Std_Volume_3D') && (info.rows || []).some(r => r.Std_Volume_3D !== null && r.Std_Volume_3D !== undefined && r.Std_Volume_3D !== '')) ySelect.value = 'Std_Volume_3D';
  else if (allowed.includes('Std_Size')) ySelect.value = 'Std_Size';
  else if (allowed.length) ySelect.value = allowed[0];
  updateSatelliteNormalizeControl();
}
function updateSatelliteNormalizeControl() {
  if (!satelliteNormalizeControl || !satelliteNormalizeCondition) return;
  if (!isSatelliteLevelTable() || !isNormalizedY(ySelect.value)) {
    satelliteNormalizeControl.style.display = 'none';
    return;
  }
  const info = DATA[tableSelect.value];
  const prev = satelliteNormalizeCondition.value;
  const conds = uniqueNonEmpty(values(info.rows, 'Condition').map(v => String(v))).sort((a,b)=>a.localeCompare(b, undefined, {numeric:true}));
  satelliteNormalizeCondition.innerHTML = '';
  conds.forEach(c => {
    let o = document.createElement('option');
    o.value = c;
    o.textContent = c;
    satelliteNormalizeCondition.appendChild(o);
  });
  if (conds.includes(prev)) satelliteNormalizeCondition.value = prev;
  else if (conds.includes('WT')) satelliteNormalizeCondition.value = 'WT';
  else if (conds.includes('Wild Type')) satelliteNormalizeCondition.value = 'Wild Type';
  else if (conds.length) satelliteNormalizeCondition.value = conds[0];
  satelliteNormalizeControl.style.display = 'flex';
}
function makeNormalizationDenominators(allRows, rawCol, refCondition) {
  const byRep = {};
  const refRows = (allRows || []).filter(r => String(r.Condition || '') === String(refCondition));
  const reps = uniqueNonEmpty(refRows.map(r => String(r.Replicate || '')));
  reps.forEach(rep => {
    const vals = refRows.filter(r => String(r.Replicate || '') === rep).map(r => safeNumber(r[rawCol])).filter(v => !isNaN(v) && isFinite(v));
    if (vals.length > 0) byRep[rep] = vals.reduce((a,b)=>a+b,0) / vals.length;
  });
  const allVals = refRows.map(r => safeNumber(r[rawCol])).filter(v => !isNaN(v) && isFinite(v));
  const global = allVals.length ? allVals.reduce((a,b)=>a+b,0) / allVals.length : null;
  return {byRep:byRep, global:global, refCondition:refCondition};
}
function satelliteYValue(row, yCol, denom) {
  const rawCol = rawYColumn(yCol);
  const v = safeNumber(row[rawCol]);
  if (isNaN(v) || !isFinite(v)) return null;
  if (!isNormalizedY(yCol)) return v;
  const rep = String(row.Replicate || '');
  const d = (denom && denom.byRep && denom.byRep[rep]) ? denom.byRep[rep] : (denom ? denom.global : null);
  if (!d || isNaN(d) || !isFinite(d) || d === 0) return null;
  return v / d;
}
function satelliteYLabel(table, yCol) {
  const rawCol = rawYColumn(yCol);
  const base = label(table, rawCol);
  if (!isNormalizedY(yCol)) return base;
  const ref = satelliteNormalizeCondition ? satelliteNormalizeCondition.value : 'reference';
  return 'Normalized ' + base + ' (reference: ' + ref + ' = 1)';
}

function perCellAllowedAxisColumns(info) {
  const dropCols = new Set([
    'Mode','Condition','Replicate','Unique_Cell','Plot_Cell','Plot_Image','Source_Input','Source_UID','Source_Path','Source_Folder','Source_short',
    'Source_File','Run_ID','Image','Threshold_method','Threshold_scope','Gaussian_smoothing_enabled'
  ]);
  const out = [];
  (info.numeric || []).forEach(c => {
    const col = String(c);
    const lab = String((info.labels && info.labels[c]) || c).toLowerCase();
    const lc = col.toLowerCase();
    const isCellId = col === 'Cell_ID' || col === 'Std_Cell_ID';
    if (dropCols.has(col) && !isCellId) return;
    if (!isCellId && (lc.includes('unique_cell') || lc.includes('condition') || lc.includes('replicate') || lc === 'mode')) return;
    if (!isCellId && (lab.includes('mode') || lab.includes('condition') || lab.includes('replicate') || lab.includes('unique cell') || lab.includes('fiji result'))) return;
    if (!out.includes(col)) out.push(col);
  });
  ['Ciliation_status','Cell_ID','Std_Cell_ID'].forEach(c => { if ((info.columns || []).includes(c) && !out.includes(c)) out.push(c); });
  return out;
}

function timeLapseAllowedYColumns(info) {
  const preferred = ['Satellite_count','Mean_area_2D','Median_area_2D','Total_area_2D','Mean_volume_3D','Median_volume_3D','Total_volume_3D','Mean_size','Median_size','Total_size','Mean_intensity','Total_intensity','Mean_distance','Mean_sphericity','Mean_circularity_2D','Mean_shape_compactness','Mean_cell_spatial_entropy','Mean_cell_pericentrosomal_clustering_index','Mean_distance_to_cell_or_image_center','Mean_centrosome_number_per_cell','Mean_cell_area_2D','Total_cell_area_2D','Mean_cell_volume_3D','Total_cell_volume_3D','Mean_cell_size','Total_cell_size'];
  const dropCols = new Set(['Condition','Replicate','Frame','Time_min','Frame_file_count','Source_Input','Source_Path']);
  const available = (info.numeric || []).filter(c => !dropCols.has(String(c)));
  const out = [];
  preferred.forEach(c => { if (available.includes(c) && !out.includes(c)) out.push(c); });
  available.forEach(c => { if (!out.includes(c)) out.push(c); });
  return out;
}

function updateSatelliteFilters() {
  if (!isFilterableLevelTable()) {
    if (satelliteFilterPanel) satelliteFilterPanel.style.display = 'none';
    return;
  }
  const info = DATA[tableSelect.value];
  const prevConditions = selectedFilterValues('conditionFilters');
  const prevReplicates = selectedFilterValues('replicateFilters');
  const prevFrames = selectedFilterValues('frameFilters');
  buildFilterCheckboxes(conditionFilters, values(info.rows, 'Condition'), prevConditions);
  buildFilterCheckboxes(replicateFilters, values(info.rows, 'Replicate'), prevReplicates);
  if (isTimeLapseTable()) {
    buildFilterCheckboxes(frameFilters, values(info.rows, 'Frame').map(frameStr), prevFrames);
    if (frameFilterSection) frameFilterSection.style.display = 'block';
  } else {
    if (frameFilters) frameFilters.innerHTML = '';
    if (frameFilterSection) frameFilterSection.style.display = 'none';
  }
  if (filterPanelTitle) {
    if (isSatelliteLevelTable()) filterPanelTitle.textContent = 'Satellite filters';
    else if (isTimeLapseTable()) filterPanelTitle.textContent = 'Time-lapse filters';
    else filterPanelTitle.textContent = 'Cell filters';
  }
  if (satelliteFilterPanel) satelliteFilterPanel.style.display = 'block';
}

function setAllFilter(containerId, checked) {
  const el = document.getElementById(containerId);
  if (!el) return;
  el.querySelectorAll('input[type="checkbox"]').forEach(cb => cb.checked = checked);
  drawPlot();
}

function rowsAfterSatelliteFilters(rows) {
  if (!isFilterableLevelTable()) return rows;
  const condSelected = new Set(selectedFilterValues('conditionFilters'));
  const repSelected = new Set(selectedFilterValues('replicateFilters'));
  const frameSelected = new Set(selectedFilterValues('frameFilters'));
  return rows.filter(r => {
    const c = String(r.Condition || '');
    const rep = String(r.Replicate || '');
    const frame = frameStr(r.Frame);
    const keepC = condSelected.size === 0 ? false : condSelected.has(c);
    const keepR = repSelected.size === 0 ? false : repSelected.has(rep);
    const keepF = isTimeLapseTable() ? (frameSelected.size === 0 ? false : frameSelected.has(frame)) : true;
    return keepC && keepR && keepF;
  });
}

function setPerCellAdvancedStatsVisibility(showTrendControls) {
  if (perCellShowTrendline && perCellShowTrendline.parentElement) perCellShowTrendline.parentElement.style.display = showTrendControls ? '' : 'none';
  if (perCellShowRValue && perCellShowRValue.parentElement) perCellShowRValue.parentElement.style.display = showTrendControls ? '' : 'none';
}

function configureControlsForTable() {
  if (isSatelliteLevelTable()) {
    typeSelect.innerHTML = '<option value="point">point graph</option>';
    typeSelect.value = 'point';
    if (graphTypeControl) graphTypeControl.style.display = 'none';
    if (xAxisControl) xAxisControl.style.display = 'none';
    if (groupControl) groupControl.style.display = 'none';
    if (satelliteStatsControl) satelliteStatsControl.style.display = 'flex';
    if (perCellStatsControl) perCellStatsControl.style.display = 'none';
    setPerCellAdvancedStatsVisibility(true);
    updateSatelliteNormalizeControl();
    if (presetbar) presetbar.style.display = 'none';
    if ([...xSelect.options].some(o => o.value === 'Condition')) xSelect.value = 'Condition';
    if (cellPanelTools) cellPanelTools.style.display = 'none';
  } else if (isPerCellTable()) {
    typeSelect.innerHTML = '<option value="point">point graph</option>';
    typeSelect.value = 'point';
    if (graphTypeControl) graphTypeControl.style.display = 'none';
    if (xAxisControl) xAxisControl.style.display = '';
    if (groupControl) groupControl.style.display = 'none';
    if (satelliteStatsControl) satelliteStatsControl.style.display = 'none';
    if (satelliteNormalizeControl) satelliteNormalizeControl.style.display = 'none';
    if (perCellStatsControl) perCellStatsControl.style.display = 'flex';
    setPerCellAdvancedStatsVisibility(true);
    if (presetbar) presetbar.style.display = 'none';
    if (cellPanelTools) cellPanelTools.style.display = 'none';
  } else if (isTimeLapseTable()) {
    typeSelect.innerHTML = '<option value="line">line plot</option><option value="point">point graph</option>';
    if (![...typeSelect.options].some(o => o.value === typeSelect.value)) typeSelect.value = 'line';
    if (graphTypeControl) graphTypeControl.style.display = '';
    if (xAxisControl) xAxisControl.style.display = '';
    if (groupControl) groupControl.style.display = 'none';
    if (satelliteStatsControl) satelliteStatsControl.style.display = 'none';
    if (satelliteNormalizeControl) satelliteNormalizeControl.style.display = 'none';
    if (perCellStatsControl) perCellStatsControl.style.display = 'none';
    if (presetbar) presetbar.style.display = 'none';
    if (cellPanelTools) cellPanelTools.style.display = 'none';
  } else {
    typeSelect.innerHTML = '<option value="point">point plot</option>';
    typeSelect.value = 'point';
    if (graphTypeControl) graphTypeControl.style.display = 'none';
    if (xAxisControl) xAxisControl.style.display = '';
    if (groupControl) groupControl.style.display = 'none';
    if (satelliteStatsControl) satelliteStatsControl.style.display = 'none';
    if (satelliteNormalizeControl) satelliteNormalizeControl.style.display = 'none';
    if (perCellStatsControl) perCellStatsControl.style.display = 'none';
    setPerCellAdvancedStatsVisibility(true);
    if (presetbar) presetbar.style.display = 'none';
  }
}

function seededJitter(i) {
  const v = ((i * 9301 + 49297) % 233280) / 233280;
  return (v - 0.5) * 0.46;
}

const CONDITION_BASE_COLORS = ['#7E3FB5', '#1F8A84', '#D17A22', '#2F6FDB', '#C44E8F', '#5C6BC0', '#78909C', '#8A6D3B', '#2E7D32', '#C62828'];
function hexToRgb(hex) {
  const h = String(hex).replace('#','');
  return {r:parseInt(h.substring(0,2),16), g:parseInt(h.substring(2,4),16), b:parseInt(h.substring(4,6),16)};
}
function rgbToHex(r,g,b) {
  const f = v => Math.max(0, Math.min(255, Math.round(v))).toString(16).padStart(2,'0');
  return '#' + f(r) + f(g) + f(b);
}
function blendWithWhite(hex, amount) {
  const c = hexToRgb(hex);
  return rgbToHex(c.r + (255-c.r)*amount, c.g + (255-c.g)*amount, c.b + (255-c.b)*amount);
}
function conditionBaseColor(condition, conditions) {
  const idx = Math.max(0, conditions.indexOf(String(condition)));
  return CONDITION_BASE_COLORS[idx % CONDITION_BASE_COLORS.length];
}
function replicateShadeColor(condition, replicate, conditions, replicates) {
  const base = conditionBaseColor(condition, conditions);
  const reps = replicates.length ? replicates : [String(replicate)];
  const ri = Math.max(0, reps.indexOf(String(replicate)));
  const amount = reps.length <= 1 ? 0.10 : (0.48 - (0.38 * ri / Math.max(1, reps.length - 1)));
  return blendWithWhite(base, amount);
}
function meanAndSd(vals) {
  const clean = vals.map(v => safeNumber(v)).filter(v => !isNaN(v) && isFinite(v));
  if (clean.length === 0) return {n:0, mean:null, sd:null};
  const mean = clean.reduce((a,b)=>a+b,0) / clean.length;
  const sd = clean.length > 1 ? Math.sqrt(clean.reduce((a,b)=>a + Math.pow(b-mean,2),0) / (clean.length - 1)) : 0;
  return {n:clean.length, mean:mean, sd:sd};
}


function satqStarLabel(p) {
  if (p === null || isNaN(p) || !isFinite(p)) return '';
  if (p < 0.0001) return '****';
  if (p < 0.001) return '***';
  if (p < 0.01) return '**';
  if (p < 0.05) return '*';
  return 'ns';
}

function satqPLabel(p) {
  if (p === null || isNaN(p) || !isFinite(p)) return 'P = NA';
  if (p < 0.0001) return 'P < 0.0001';
  return 'P = ' + p.toPrecision(3);
}

function satqLogGamma(z) {
  const g = 7;
  const p = [0.99999999999980993,676.5203681218851,-1259.1392167224028,771.32342877765313,-176.61502916214059,12.507343278686905,-0.13857109526572012,9.9843695780195716e-6,1.5056327351493116e-7];
  if (z < 0.5) return Math.log(Math.PI) - Math.log(Math.sin(Math.PI * z)) - satqLogGamma(1 - z);
  z -= 1;
  let x = p[0];
  for (let i = 1; i < p.length; i++) x += p[i] / (z + i);
  const t = z + g + 0.5;
  return 0.5 * Math.log(2 * Math.PI) + (z + 0.5) * Math.log(t) - t + Math.log(x);
}

function satqBetaCf(a, b, x) {
  const MAXIT = 100;
  const EPS = 3e-7;
  const FPMIN = 1e-30;
  let qab = a + b, qap = a + 1, qam = a - 1;
  let c = 1, d = 1 - qab * x / qap;
  if (Math.abs(d) < FPMIN) d = FPMIN;
  d = 1 / d;
  let h = d;
  for (let m = 1; m <= MAXIT; m++) {
    const m2 = 2 * m;
    let aa = m * (b - m) * x / ((qam + m2) * (a + m2));
    d = 1 + aa * d;
    if (Math.abs(d) < FPMIN) d = FPMIN;
    c = 1 + aa / c;
    if (Math.abs(c) < FPMIN) c = FPMIN;
    d = 1 / d;
    h *= d * c;
    aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
    d = 1 + aa * d;
    if (Math.abs(d) < FPMIN) d = FPMIN;
    c = 1 + aa / c;
    if (Math.abs(c) < FPMIN) c = FPMIN;
    d = 1 / d;
    const del = d * c;
    h *= del;
    if (Math.abs(del - 1.0) < EPS) break;
  }
  return h;
}

function satqIBeta(x, a, b) {
  if (x <= 0) return 0;
  if (x >= 1) return 1;
  const bt = Math.exp(satqLogGamma(a + b) - satqLogGamma(a) - satqLogGamma(b) + a * Math.log(x) + b * Math.log(1 - x));
  if (x < (a + 1) / (a + b + 2)) return bt * satqBetaCf(a, b, x) / a;
  return 1 - bt * satqBetaCf(b, a, 1 - x) / b;
}

function satqStudentTCdf(t, df) {
  if (df <= 0 || !isFinite(t)) return NaN;
  const x = df / (df + t * t);
  const ib = satqIBeta(x, df / 2, 0.5);
  return t >= 0 ? 1 - 0.5 * ib : 0.5 * ib;
}

function satqTPValue(t, df) {
  const cdf = satqStudentTCdf(Math.abs(t), df);
  if (isNaN(cdf)) return NaN;
  return Math.max(0, Math.min(1, 2 * (1 - cdf)));
}

function satqFPValue(f, df1, df2) {
  if (!isFinite(f) || f < 0 || df1 <= 0 || df2 <= 0) return NaN;
  const x = (df1 * f) / (df1 * f + df2);
  const cdf = satqIBeta(x, df1 / 2, df2 / 2);
  return Math.max(0, Math.min(1, 1 - cdf));
}

function satqMeanVar(vals) {
  const a = vals.map(Number).filter(v => !isNaN(v) && isFinite(v));
  const n = a.length;
  if (n === 0) return {n:0, mean:NaN, variance:NaN};
  const mean = a.reduce((s,v)=>s+v,0) / n;
  const variance = n > 1 ? a.reduce((s,v)=>s + Math.pow(v - mean, 2), 0) / (n - 1) : 0;
  return {n:n, mean:mean, variance:variance};
}

function satqGroupValues(rows, categories, catGetter, yGetter) {
  const out = {};
  categories.forEach(c => out[String(c)] = []);
  rows.forEach(r => {
    const c = String(catGetter(r));
    if (!(c in out)) return;
    const y = safeNumber(yGetter(r));
    if (!isNaN(y) && isFinite(y)) out[c].push(y);
  });
  return categories.map(c => out[String(c)] || []);
}

function satqGroupReplicateMeans(rows, categories, catGetter, yGetter) {
  const maps = {};
  categories.forEach(c => maps[String(c)] = {});
  rows.forEach((r, idx) => {
    const c = String(catGetter(r));
    if (!(c in maps)) return;
    const y = safeNumber(yGetter(r));
    if (isNaN(y) || !isFinite(y)) return;
    let rep = r.Replicate !== undefined && r.Replicate !== null && String(r.Replicate).trim() !== '' ? String(r.Replicate) : ('row_' + idx);
    if (!maps[c][rep]) maps[c][rep] = [];
    maps[c][rep].push(y);
  });
  return categories.map(c => {
    const m = maps[String(c)] || {};
    return Object.keys(m).map(k => m[k].reduce((a,b)=>a+b,0) / m[k].length);
  });
}

function satqSamples(rows, categories, catGetter, yGetter) {
  const level = statSampleLevel ? statSampleLevel.value : 'replicate_means';
  return level === 'each_value' ? satqGroupValues(rows, categories, catGetter, yGetter) : satqGroupReplicateMeans(rows, categories, catGetter, yGetter);
}

function satqUnpairedT(a, b, welch=true) {
  const A = satqMeanVar(a), B = satqMeanVar(b);
  if (A.n < 2 || B.n < 2) return {p:NaN, test:welch ? 'Welch t-test' : 'Student t-test'};
  const diff = A.mean - B.mean;
  if (welch) {
    const se2 = A.variance / A.n + B.variance / B.n;
    if (se2 <= 0) return {p:NaN, test:'Welch t-test'};
    const t = diff / Math.sqrt(se2);
    const df = Math.pow(se2, 2) / (Math.pow(A.variance / A.n, 2) / (A.n - 1) + Math.pow(B.variance / B.n, 2) / (B.n - 1));
    return {p:satqTPValue(t, df), test:'Welch t-test'};
  }
  const pooled = ((A.n - 1) * A.variance + (B.n - 1) * B.variance) / (A.n + B.n - 2);
  if (pooled <= 0) return {p:NaN, test:'Student t-test'};
  const t = diff / Math.sqrt(pooled * (1 / A.n + 1 / B.n));
  return {p:satqTPValue(t, A.n + B.n - 2), test:'Student t-test'};
}

function satqPairedT(rows, categories, catGetter, yGetter) {
  if (categories.length !== 2) return {p:NaN, test:'Paired t-test'};
  const maps = {};
  categories.forEach(c => maps[String(c)] = {});
  rows.forEach((r, idx) => {
    const c = String(catGetter(r));
    if (!(c in maps)) return;
    const y = safeNumber(yGetter(r));
    if (isNaN(y) || !isFinite(y)) return;
    const rep = r.Replicate !== undefined && r.Replicate !== null && String(r.Replicate).trim() !== '' ? String(r.Replicate) : ('row_' + idx);
    if (!maps[c][rep]) maps[c][rep] = [];
    maps[c][rep].push(y);
  });
  const reps = Object.keys(maps[String(categories[0])] || {}).filter(k => maps[String(categories[1])] && maps[String(categories[1])][k]);
  const diffs = reps.map(k => {
    const a = maps[String(categories[0])][k].reduce((s,v)=>s+v,0) / maps[String(categories[0])][k].length;
    const b = maps[String(categories[1])][k].reduce((s,v)=>s+v,0) / maps[String(categories[1])][k].length;
    return a - b;
  });
  const D = satqMeanVar(diffs);
  if (D.n < 2 || D.variance <= 0) return {p:NaN, test:'Paired t-test'};
  const t = D.mean / Math.sqrt(D.variance / D.n);
  return {p:satqTPValue(t, D.n - 1), test:'Paired t-test'};
}

function satqAnova(samples) {
  const groups = samples.map(s => s.map(Number).filter(v => !isNaN(v) && isFinite(v))).filter(s => s.length > 0);
  const k = groups.length;
  const nTotal = groups.reduce((s,g)=>s+g.length,0);
  if (k < 2 || nTotal <= k) return {p:NaN, test:'One-way ANOVA'};
  const all = [].concat(...groups);
  const grand = all.reduce((s,v)=>s+v,0) / all.length;
  let ssb = 0, ssw = 0;
  groups.forEach(g => {
    const mean = g.reduce((s,v)=>s+v,0) / g.length;
    ssb += g.length * Math.pow(mean - grand, 2);
    g.forEach(v => ssw += Math.pow(v - mean, 2));
  });
  const df1 = k - 1, df2 = nTotal - k;
  const msb = ssb / df1, msw = ssw / df2;
  if (msw <= 0) return {p:NaN, test:'One-way ANOVA'};
  const f = msb / msw;
  return {p:satqFPValue(f, df1, df2), test:'One-way ANOVA'};
}

function satqAddSignificance(shapes, annotations, categories, catIndex, rows, catGetter, yGetter, yMin, yMax, ySpan) {
  if (!showSignificance || !showSignificance.checked) return null;
  if (!categories || categories.length < 2) return null;
  let test = statTestSelect ? statTestSelect.value : 'auto';
  if (test === 'auto') test = categories.length === 2 ? 'welch_t' : 'anova';

  let result = null;
  const samples = satqSamples(rows, categories, catGetter, yGetter);
  if (test === 'paired_t') {
    result = satqPairedT(rows, categories, catGetter, yGetter);
  } else if (test === 'student_t') {
    if (categories.length !== 2) return null;
    result = satqUnpairedT(samples[0], samples[1], false);
  } else if (test === 'welch_t') {
    if (categories.length !== 2) return null;
    result = satqUnpairedT(samples[0], samples[1], true);
  } else if (test === 'anova') {
    result = satqAnova(samples);
  }
  if (!result) return null;

  const x0 = 1;
  const x1 = categories.length;
  const y = yMax + ySpan * 0.16;
  const h = ySpan * 0.045;
  shapes.push({type:'line', xref:'x', yref:'y', x0:x0, x1:x0, y0:y, y1:y+h, line:{color:'#111111', width:1.7}});
  shapes.push({type:'line', xref:'x', yref:'y', x0:x0, x1:x1, y0:y+h, y1:y+h, line:{color:'#111111', width:1.7}});
  shapes.push({type:'line', xref:'x', yref:'y', x0:x1, x1:x1, y0:y+h, y1:y, line:{color:'#111111', width:1.7}});

  const stars = satqStarLabel(result.p);
  const pTxt = satqPLabel(result.p);
  const sampleTxt = statSampleLevel && statSampleLevel.value === 'each_value' ? 'each value' : 'replicate averages';
  annotations.push({
    x:(x0+x1)/2,
    y:y + h * 1.65,
    text:(stars || 'NA') + '<br><span style="font-size:10px">' + pTxt + '<br>' + result.test + ', ' + sampleTxt + '</span>',
    showarrow:false,
    align:'center',
    font:{color:'#111111', size:13},
    bgcolor:'rgba(255,255,255,0.88)',
    bordercolor:'#111111',
    borderpad:3
  });
  return y + h * 2.4;
}

function drawSatellitePointPlot(table, info, rows, x, y, allRows) {
  // Satellite-level graphs are intentionally fixed: X-axis = condition, one point = one satellite.
  x = 'Condition';
  if (!rows || rows.length === 0) {
    document.getElementById('plot').style.height = '680px';
    Plotly.newPlot('plot', [], {title:{text:'No satellites match the selected filters', font:{color:SATQ_PURPLE}}, paper_bgcolor:'white', plot_bgcolor:'white'});
    return;
  }

  updateSatelliteNormalizeControl();
  const rawY = rawYColumn(y);
  const refCond = satelliteNormalizeCondition ? satelliteNormalizeCondition.value : '';
  const denom = isNormalizedY(y) ? makeNormalizationDenominators(allRows || rows, rawY, refCond) : null;
  const yValsAll = rows.map(r => satelliteYValue(r, y, denom)).filter(v => v !== null && !isNaN(Number(v)) && isFinite(Number(v)));
  const conditionsVisible = uniqueNonEmpty(values(rows, 'Condition').map(v => String(v))).sort((a,b)=>a.localeCompare(b, undefined, {numeric:true}));
  const replicatesVisible = uniqueNonEmpty(values(rows, 'Replicate').map(v => String(v))).sort((a,b)=>a.localeCompare(b, undefined, {numeric:true}));
  const catIndex = {};
  conditionsVisible.forEach((c, i) => { catIndex[c] = i + 1; });
  const repCount = Math.max(1, replicatesVisible.length);
  const showStats = satelliteShowMeanSD ? satelliteShowMeanSD.checked : true;

  let traces = [];
  let outlierX = [];
  let outlierY = [];
  let outlierText = [];
  let outlierRows = [];
  const showOutliersNow = outlierToggleOn();
  if (showOutliersNow) {
    computeDynamicOutliers(
      rows,
      satelliteYLabel(table, y),
      r => satelliteYValue(r, y, denom),
      r => 'condition ' + String(r.Condition || 'visible') + (r.Mode ? ' / mode ' + r.Mode : '')
    );
  } else {
    clearDynamicOutliers(rows);
  }

  conditionsVisible.forEach((cond, ci) => {
    replicatesVisible.forEach((rep, ri) => {
      const rr = rows.filter(r => String(r.Condition || '') === String(cond) && String(r.Replicate || '') === String(rep));
      if (rr.length === 0) return;
      const yy = rr.map(r => satelliteYValue(r, y, denom));
      const repOffset = repCount <= 1 ? 0 : (ri - (repCount - 1) / 2) * Math.min(0.20, 0.52 / repCount);
      const xx = rr.map((r, i) => (catIndex[String(cond)] || 0) + repOffset + seededJitter(i + ci * 997 + ri * 389) * 0.28);
      const color = replicateShadeColor(cond, rep, conditionsVisible, replicatesVisible);
      const hover = rr.map((r, i) => {
        let txt = [];
        txt.push('Condition: ' + (r.Condition || ''));
        txt.push('Replicate: ' + (r.Replicate || ''));
        if (r.Plot_Cell) txt.push(prettyCellLabel(r.Plot_Cell));
        if (r.Unique_Satellite) txt.push(String(r.Unique_Satellite).replace(/^.*\\|\\s*/, ''));
        txt.push(satelliteYLabel(table, y) + ': ' + fmt(satelliteYValue(r, y, denom)));
        return txt.join('<br>');
      });
      if (showOutliersNow) {
        rr.forEach((r, i) => {
          const yv = yy[i];
          if (isOutlierRow(r) && yv !== null && !isNaN(Number(yv)) && isFinite(Number(yv))) {
            outlierX.push(xx[i]);
            outlierY.push(yv);
            outlierText.push(hover[i] + outlierHoverSuffix(r));
            outlierRows.push(r);
          }
        });
      }
      traces.push({
        type:'scatter',
        mode:'markers',
        x:xx,
        y:yy,
        name:String(cond) + ' / ' + String(rep),
        legendgroup:String(cond),
        marker:{size:7, opacity:0.84, color:color, line:{color:'#333333', width:0.25}},
        text:hover,
        hovertemplate:'%{text}<extra></extra>'
      });
    });
  });

  addOutlierTrace(traces, outlierX, outlierY, outlierText, 'satellite-level');
  updateOutlierPanel(outlierRows, 'Outliers are recalculated from the current satellite-level metric within each visible condition.');

  let annotations = [];
  let shapes = [];
  let yMin = yValsAll.length ? Math.min(...yValsAll) : 0;
  let yMax = yValsAll.length ? Math.max(...yValsAll) : 1;
  let ySpan = Math.max(1e-9, yMax - yMin);

  if (showStats) {
    conditionsVisible.forEach((cond, ci) => {
      const vals = rows.filter(r => String(r.Condition || '') === String(cond)).map(r => satelliteYValue(r, y, denom)).filter(v => v !== null && !isNaN(v) && isFinite(v));
      const st = meanAndSd(vals);
      const xpos = catIndex[cond];
      annotations.push({
        x:xpos,
        y:yMax + ySpan * 0.10,
        text:'n = ' + st.n,
        showarrow:false,
        font:{color:SATQ_PURPLE_DARK, size:12},
        bgcolor:'rgba(255,255,255,0.82)',
        bordercolor:SATQ_LIGHT_PURPLE,
        borderpad:3
      });
      if (st.n > 0 && st.mean !== null) {
        const sd = st.sd || 0;
        const x0 = xpos - 0.26;
        const x1 = xpos + 0.26;
        const cap0 = xpos - 0.12;
        const cap1 = xpos + 0.12;
        shapes.push({type:'line', xref:'x', yref:'y', x0:x0, x1:x1, y0:st.mean, y1:st.mean, line:{color:'#111111', width:3}});
        shapes.push({type:'line', xref:'x', yref:'y', x0:xpos, x1:xpos, y0:st.mean - sd, y1:st.mean + sd, line:{color:'#111111', width:2}});
        shapes.push({type:'line', xref:'x', yref:'y', x0:cap0, x1:cap1, y0:st.mean - sd, y1:st.mean - sd, line:{color:'#111111', width:2}});
        shapes.push({type:'line', xref:'x', yref:'y', x0:cap0, x1:cap1, y0:st.mean + sd, y1:st.mean + sd, line:{color:'#111111', width:2}});
        yMin = Math.min(yMin, st.mean - sd);
        yMax = Math.max(yMax, st.mean + sd);
      }
    });
  }

  ySpan = Math.max(1e-9, yMax - yMin);
  const yRange = [yMin - ySpan * 0.12, yMax + ySpan * (showStats ? 0.24 : 0.14)];
  const title = satelliteYLabel(table, y) + ' by condition';
  document.getElementById('plot').style.height = '680px';
  Plotly.newPlot('plot', traces, {
    title:{text:title, font:{color:SATQ_PURPLE, size:19}},
    xaxis:{
      title:{text:'Condition', font:{color:'#333333'}},
      tickmode:'array',
      tickvals:conditionsVisible.map((_, i) => i + 1),
      ticktext:conditionsVisible,
      range:[0.45, conditionsVisible.length + 0.55],
      tickangle:conditionsVisible.length > 4 ? -35 : 0,
      zeroline:false,
      linecolor:'#444444',
      tickfont:{color:'#333333'},
      automargin:true
    },
    yaxis:{title:{text:satelliteYLabel(table, y), font:{color:'#333333'}}, range:yRange, zeroline:false, linecolor:'#444444', tickfont:{color:'#333333'}, automargin:true},
    shapes:shapes,
    annotations:annotations,
    plot_bgcolor:'white',
    paper_bgcolor:'white',
    hovermode:'closest',
    legend:{title:{text:'Condition / replicate'}},
    margin:{l:86, r:45, t:84, b:115}
  }, {responsive:true, toImageButtonOptions:{format:'png', filename:'SatelliteQ_satellite_level_points', scale:3}});
}

function drawPerCellPointPlot(table, info, rows, x, y) {
  // Per-cell graphs: one point = one cell measurement. Conditions use different colors;
  // replicates use different shades of the same condition color.
  if (!rows || rows.length === 0) {
    document.getElementById('plot').style.height = '680px';
    Plotly.newPlot('plot', [], {title:{text:'No cells match the selected filters', font:{color:SATQ_PURPLE}}, paper_bgcolor:'white', plot_bgcolor:'white'});
    return;
  }

  const xCategorical = isCategoricalAxis(table, x) || !isNumericCol(table, x);
  const conditionsVisible = uniqueNonEmpty(values(rows, 'Condition').map(v => String(v))).sort((a,b)=>a.localeCompare(b, undefined, {numeric:true}));
  const replicatesVisible = uniqueNonEmpty(values(rows, 'Replicate').map(v => String(v))).sort((a,b)=>a.localeCompare(b, undefined, {numeric:true}));
  const xTitle = label(table, x);
  const yTitle = label(table, y);
  const repCount = Math.max(1, replicatesVisible.length);
  let traces = [];
  let outlierX = [];
  let outlierY = [];
  let outlierText = [];
  let outlierRows = [];
  const showOutliersNow = outlierToggleOn();
  if (showOutliersNow) {
    computeDynamicOutliers(
      rows,
      label(table, y),
      r => safeNumber(r[y]),
      r => dynamicGroupForRow(r, x, table, 'point')
    );
  } else {
    clearDynamicOutliers(rows);
  }
  let allXNum = [];
  let allYNum = [];

  let categories = [];
  let catIndex = {};
  if (xCategorical) {
    categories = uniqueNonEmpty(values(rows, x).map(v => String(v))).sort((a,b)=>a.localeCompare(b, undefined, {numeric:true}));
    categories.forEach((c, i) => { catIndex[c] = i + 1; });
  }

  conditionsVisible.forEach((cond, ci) => {
    replicatesVisible.forEach((rep, ri) => {
      const rrRaw = rows.filter(r => String(r.Condition || '') === String(cond) && String(r.Replicate || '') === String(rep));
      const rr = rrRaw.filter(r => {
        const yv = safeNumber(r[y]);
        if (isNaN(yv) || !isFinite(yv)) return false;
        if (!xCategorical) {
          const xv = safeNumber(r[x]);
          return !isNaN(xv) && isFinite(xv);
        }
        return String(r[x] || '') !== '';
      });
      if (rr.length === 0) return;
      const color = replicateShadeColor(cond, rep, conditionsVisible, replicatesVisible);
      const yy = rr.map(r => safeNumber(r[y]));
      let xx;
      if (xCategorical) {
        const repOffset = repCount <= 1 ? 0 : (ri - (repCount - 1) / 2) * Math.min(0.20, 0.52 / repCount);
        xx = rr.map((r, i) => (catIndex[String(r[x])] || 0) + repOffset + seededJitter(i + ci * 997 + ri * 389) * 0.28);
      } else {
        xx = rr.map(r => safeNumber(r[x]));
      }
      xx.forEach((v, i) => { if (!isNaN(v) && isFinite(v) && !isNaN(yy[i]) && isFinite(yy[i])) { allXNum.push(v); allYNum.push(yy[i]); } });
      const hover = rr.map((r, i) => {
        let txt = [];
        txt.push('Condition: ' + (r.Condition || ''));
        txt.push('Replicate: ' + (r.Replicate || ''));
        if (r.Plot_Cell) txt.push(prettyCellLabel(r.Plot_Cell));
        txt.push(xTitle + ': ' + r[x]);
        txt.push(yTitle + ': ' + fmt(safeNumber(r[y])));
        return txt.join('<br>');
      });
      if (showOutliersNow) {
        rr.forEach((r, i) => {
          if (isOutlierRow(r) && !isNaN(Number(yy[i])) && isFinite(Number(yy[i]))) {
            outlierX.push(xx[i]);
            outlierY.push(yy[i]);
            outlierText.push(hover[i] + outlierHoverSuffix(r));
            outlierRows.push(r);
          }
        });
      }
      traces.push({
        type:'scatter',
        mode:'markers',
        x:xx,
        y:yy,
        name:String(cond) + ' / ' + String(rep),
        legendgroup:String(cond),
        marker:{size:8, opacity:0.86, color:color, line:{color:'#333333', width:0.35}},
        text:hover,
        hovertemplate:'%{text}<extra></extra>'
      });
    });
  });

  addOutlierTrace(traces, outlierX, outlierY, outlierText, 'per-cell');
  updateOutlierPanel(outlierRows, 'Outliers are recalculated from the current per-cell metric and graph set.');

  let shapes = [];
  let annotations = [];
  let yClean = allYNum.filter(v => !isNaN(v) && isFinite(v));
  let yMin = yClean.length ? Math.min(...yClean) : 0;
  let yMax = yClean.length ? Math.max(...yClean) : 1;
  let ySpan = Math.max(1e-9, yMax - yMin);

  if (perCellShowMeanSD && perCellShowMeanSD.checked) {
    if (xCategorical) {
      categories.forEach(cat => {
        const vals = rows.filter(r => String(r[x]) === String(cat)).map(r => safeNumber(r[y])).filter(v => !isNaN(v) && isFinite(v));
        const st = meanAndSd(vals);
        const xpos = catIndex[cat];
        if (st.n > 0) {
          annotations.push({x:xpos, y:yMax + ySpan * 0.10, text:'n = ' + st.n, showarrow:false, font:{color:SATQ_PURPLE_DARK, size:12}, bgcolor:'rgba(255,255,255,0.82)', bordercolor:SATQ_LIGHT_PURPLE, borderpad:3});
          shapes.push({type:'line', xref:'x', yref:'y', x0:xpos-0.24, x1:xpos+0.24, y0:st.mean, y1:st.mean, line:{color:'#111111', width:3}});
          const sd = st.sd || 0;
          shapes.push({type:'line', xref:'x', yref:'y', x0:xpos, x1:xpos, y0:st.mean-sd, y1:st.mean+sd, line:{color:'#111111', width:2}});
          shapes.push({type:'line', xref:'x', yref:'y', x0:xpos-0.10, x1:xpos+0.10, y0:st.mean-sd, y1:st.mean-sd, line:{color:'#111111', width:2}});
          shapes.push({type:'line', xref:'x', yref:'y', x0:xpos-0.10, x1:xpos+0.10, y0:st.mean+sd, y1:st.mean+sd, line:{color:'#111111', width:2}});
          yMin = Math.min(yMin, st.mean - sd);
          yMax = Math.max(yMax, st.mean + sd);
        }
      });
    } else {
      const st = meanAndSd(yClean);
      const xClean = allXNum.filter(v => !isNaN(v) && isFinite(v));
      const xMin = xClean.length ? Math.min(...xClean) : 0;
      const xMax = xClean.length ? Math.max(...xClean) : 1;
      const xSpan = Math.max(1e-9, xMax - xMin);
      if (st.n > 0) {
        shapes.push({type:'line', xref:'x', yref:'y', x0:xMin-xSpan*0.03, x1:xMax+xSpan*0.03, y0:st.mean, y1:st.mean, line:{color:'#111111', width:2}});
        if (st.sd !== null) {
          shapes.push({type:'line', xref:'x', yref:'y', x0:xMin-xSpan*0.03, x1:xMax+xSpan*0.03, y0:st.mean+st.sd, y1:st.mean+st.sd, line:{color:'#111111', width:1.4, dash:'dot'}});
          shapes.push({type:'line', xref:'x', yref:'y', x0:xMin-xSpan*0.03, x1:xMax+xSpan*0.03, y0:st.mean-st.sd, y1:st.mean-st.sd, line:{color:'#111111', width:1.4, dash:'dot'}});
        }
        annotations.push({text:'n = ' + st.n + '<br>mean = ' + fmt(st.mean) + (st.sd !== null ? '<br>SD = ' + fmt(st.sd) : ''), xref:'paper', yref:'paper', x:0.98, y:0.98, xanchor:'right', showarrow:false, align:'right', bgcolor:'rgba(255,255,255,0.9)', bordercolor:SATQ_LIGHT_PURPLE, font:{color:SATQ_PURPLE_DARK, size:12}});
      }
    }
  }

  if (!xCategorical && perCellShowTrendline && perCellShowTrendline.checked) {
    const fit = linearRegressionR2(allXNum, allYNum);
    if (fit) {
      const rVal = correlationRFromFit(fit);
      traces.push({type:'scatter', mode:'lines', x:[fit.x1, fit.x2], y:[fit.y1, fit.y2], name:'Trendline' + (perCellShowRValue && perCellShowRValue.checked ? ' (R=' + (rVal === null ? 'NA' : rVal.toFixed(3)) + ')' : ''), line:{color:'#111111', width:2, dash:'dash'}, hovertemplate:'R=' + (rVal === null ? 'NA' : rVal.toFixed(3)) + '<br>R²=' + fit.r2.toFixed(3) + '<br>n=' + fit.n + '<extra></extra>'});
      if (perCellShowRValue && perCellShowRValue.checked) {
        annotations.push({text:'R = ' + (rVal === null ? 'NA' : rVal.toFixed(3)) + '<br>R² = ' + fit.r2.toFixed(3) + '<br>n = ' + fit.n, xref:'paper', yref:'paper', x:0.02, y:0.98, showarrow:false, align:'left', bgcolor:'rgba(255,255,255,0.9)', bordercolor:SATQ_LIGHT_PURPLE, font:{color:SATQ_PURPLE_DARK, size:13}});
      }
    }
  }

  ySpan = Math.max(1e-9, yMax - yMin);
  let yRange = [yMin - ySpan * 0.12, yMax + ySpan * 0.20];
  let xAxisConfig = {title:{text:xTitle, font:{color:'#333333'}}, zeroline:false, linecolor:'#444444', tickfont:{color:'#333333'}, automargin:true};
  if (xCategorical) {
    xAxisConfig.tickmode = 'array';
    xAxisConfig.tickvals = categories.map((_, i) => i + 1);
    xAxisConfig.ticktext = categories;
    xAxisConfig.range = [0.45, categories.length + 0.55];
    xAxisConfig.tickangle = categories.length > 4 ? -35 : 0;
  }

  document.getElementById('plot').style.height = '680px';
  Plotly.newPlot('plot', traces, {
    title:{text:yTitle + ' vs ' + xTitle, font:{color:SATQ_PURPLE, size:19}},
    xaxis:xAxisConfig,
    yaxis:{title:{text:yTitle, font:{color:'#333333'}}, range:yRange, zeroline:false, linecolor:'#444444', tickfont:{color:'#333333'}, automargin:true},
    shapes:shapes,
    annotations:annotations,
    plot_bgcolor:'white',
    paper_bgcolor:'white',
    hovermode:'closest',
    legend:{title:{text:'Condition / replicate'}},
    margin:{l:86, r:45, t:84, b:115}
  }, {responsive:true, toImageButtonOptions:{format:'png', filename:'SatelliteQ_per_cell_points', scale:3}});
}

function getCellColumnForCurrentTable() {
  const table = tableSelect.value;
  const cols = DATA[table].columns || [];
  if (cols.includes('Plot_Cell')) return 'Plot_Cell';
  if (cols.includes('Unique_Cell')) return 'Unique_Cell';
  if (cols.includes('Cell_ID')) return 'Cell_ID';
  if (cols.includes('Std_Cell_ID')) return 'Std_Cell_ID';
  return null;
}

function prettyCellLabel(value) {
  let s = String(value);
  if (s.includes('|')) return s;
  const m = s.match(/Cell\\s+([0-9]+(?:\\.0+)?|[0-9]*\\.[0-9]+)/i);
  if (m) {
    let n = m[1];
    if (n.endsWith('.0')) n = n.slice(0, -2);
    return 'Cell ' + n;
  }
  if (/^[0-9]+(?:\\.0+)?$/.test(s)) {
    if (s.endsWith('.0')) s = s.slice(0, -2);
    return 'Cell ' + s;
  }
  return s;
}

function cellSortKey(value) {
  const label = prettyCellLabel(value);
  const m = label.match(/Cell\\s+([0-9.]+)/i);
  return m ? Number(m[1]) : 999999;
}

function availableCellsForCurrentTable() {
  const table = tableSelect.value;
  const info = DATA[table];
  const cellCol = getCellColumnForCurrentTable();
  if (!cellCol) return [];
  return uniqueNonEmpty(info.rows.map(r => r[cellCol]).map(v => String(v))).sort((a,b) => cellSortKey(a) - cellSortKey(b) || String(a).localeCompare(String(b)));
}

function updateCellPanelVisibility() {
  if (cellPanelTools) cellPanelTools.style.display = 'none';
}

function updateCellSelector() {
  const cells = availableCellsForCurrentTable();
  const previous = selectedCells();
  cellPanelSelect.innerHTML = '';
  cells.forEach(c => {
    const lab = document.createElement('label');
    lab.className = 'cellchip';
    const cb = document.createElement('input');
    cb.type = 'checkbox';
    cb.value = c;
    cb.checked = previous.includes(c);
    cb.addEventListener('change', drawPlot);
    const span = document.createElement('span');
    span.textContent = prettyCellLabel(c);
    lab.appendChild(cb);
    lab.appendChild(span);
    cellPanelSelect.appendChild(lab);
  });
  updateCellPanelVisibility();
}

function selectedCells() {
  return Array.from(cellPanelSelect.querySelectorAll('input[type="checkbox"]:checked')).map(o => o.value);
}

function selectAllCells() {
  cellPanelSelect.querySelectorAll('input[type="checkbox"]').forEach(o => o.checked = true);
  drawPlot();
}

function clearSelectedCells() {
  cellPanelSelect.querySelectorAll('input[type="checkbox"]').forEach(o => o.checked = false);
  drawPlot();
}

function rowsAfterCellSelection(rows) {
  // Only filter the main graph by selected cells in cell-panel mode.
  // Other graph types stay clean and uncluttered.
  if (typeSelect.value !== 'cell_panels') return rows;
  const cellCol = getCellColumnForCurrentTable();
  const cells = selectedCells();
  if (!cellCol || cells.length === 0) return rows;
  const keep = new Set(cells.map(String));
  return rows.filter(r => keep.has(String(r[cellCol])));
}

function presetCellPanels() {
  switchToTable('Satellite-level');
  typeSelect.value = 'cell_panels';
  if (cellPlotType) cellPlotType.value = 'point';
  setIfExists(ySelect, 'Std_Size');
  setIfExists(groupSelect, '');
  updateCellSelector();
  updateCellPanelVisibility();
  drawPlot();
}

function presetAreaVsIntensityPanels() {
  switchToTable('Satellite-level');
  typeSelect.value = 'cell_panels';
  if (cellPlotType) cellPlotType.value = 'scatterxy';
  setIfExists(xSelect, 'Std_Size');
  setIfExists(ySelect, 'Std_Intensity');
  setIfExists(groupSelect, '');
  updateCellSelector();
  updateCellPanelVisibility();
  drawPlot();
}

function init() {
  addOptions(tableSelect, Object.keys(DATA));

  updateColumns();
}

function updateColumns() {
  const table = tableSelect.value;
  const info = DATA[table];
  if (isSatelliteLevelTable(table)) {
    addOptions(xSelect, info.columns, info.labels);
    populateSatelliteYOptions(info);
  } else if (isPerCellTable(table)) {
    const allowed = perCellAllowedAxisColumns(info);
    addOptions(xSelect, allowed.length ? allowed : (info.numeric.length ? info.numeric : info.columns), info.labels);
    addOptions(ySelect, allowed.length ? allowed : (info.numeric.length ? info.numeric : info.columns), info.labels);
  } else if (isTimeLapseTable(table)) {
    const xOpts = info.columns.filter(c => ['Frame','Time_min'].includes(c));
    addOptions(xSelect, xOpts.length ? xOpts : info.columns, info.labels);
    const allowed = timeLapseAllowedYColumns(info);
    addOptions(ySelect, allowed.length ? allowed : (info.numeric.length ? info.numeric : info.columns), info.labels);
  } else {
    addOptions(xSelect, info.columns, info.labels);
    addOptions(ySelect, info.numeric.length ? info.numeric : info.columns, info.labels);
  }
  addOptions(groupSelect, [], {}, true);
  if (info.defaultX && [...xSelect.options].some(o => o.value === info.defaultX)) xSelect.value = info.defaultX;
  if (!isSatelliteLevelTable(table) && info.defaultY && [...ySelect.options].some(o => o.value === info.defaultY)) ySelect.value = info.defaultY;
  configureControlsForTable();
  updateSatelliteFilters();
  updateCellSelector();
  updateCellPanelVisibility();
  drawPlot();
}

function values(rows, col) { return rows.map(r => r[col]); }
function frameStr(v) { return (v === null || v === undefined) ? '' : String(v); }
function safeNumber(v) {
  if (v === null || v === undefined) return NaN;
  const t = String(v).trim();
  if (t === '' || ['NA','NAN','NONE','NULL','N/A'].includes(t.toUpperCase())) return NaN;
  const n = Number(v);
  return Number.isFinite(n) ? n : NaN;
}
function numericValues(rows, col) { return rows.map(r => { const n = safeNumber(r[col]); return Number.isFinite(n) ? n : null; }); }
function uniqueNonEmpty(arr) { return Array.from(new Set(arr.filter(v => v !== null && v !== undefined && v !== ''))); }
function rowsForGroup(rows, col, val) { return rows.filter(r => String(r[col]) === String(val)); }
function setIfExists(sel, value) { if ([...sel.options].some(o => o.value === value)) sel.value = value; }

function switchToTable(prefix) {
  const key = Object.keys(DATA).find(k => k.startsWith(prefix));
  if (key) { tableSelect.value = key; updateColumns(); }
}

function presetSatelliteAreaByCell() {
  switchToTable('Satellite-level');
  typeSelect.value = 'point';
  setIfExists(xSelect, 'Plot_Cell');
  setIfExists(ySelect, 'Std_Size');
  setIfExists(groupSelect, '');
  drawPlot();
}
function presetSatelliteIntensityByCell() {
  switchToTable('Satellite-level');
  typeSelect.value = 'point';
  setIfExists(xSelect, 'Plot_Cell');
  setIfExists(ySelect, 'Std_Intensity');
  setIfExists(groupSelect, '');
  drawPlot();
}
function presetPerCellCounts() {
  switchToTable('Per-cell');
  typeSelect.value = 'bar';
  setIfExists(xSelect, 'Plot_Cell');
  setIfExists(ySelect, 'Satellite_count');
  setIfExists(groupSelect, '');
  drawPlot();
}
function presetPerCellMeanSize() {
  switchToTable('Per-cell');
  typeSelect.value = 'point';
  setIfExists(xSelect, 'Plot_Cell');
  setIfExists(ySelect, 'Mean_size');
  setIfExists(groupSelect, '');
  drawPlot();
}
function presetPerImageCounts() {
  switchToTable('Per-image');
  typeSelect.value = 'bar';
  setIfExists(xSelect, 'Plot_Image');
  setIfExists(ySelect, 'Satellite_count');
  setIfExists(groupSelect, '');
  drawPlot();
}

function linearRegressionR2(xArr, yArr) {
  let pts = [];
  for (let i=0; i<xArr.length; i++) {
    let xVal = safeNumber(xArr[i]);
    let yVal = safeNumber(yArr[i]);
    if (!isNaN(xVal) && !isNaN(yVal) && isFinite(xVal) && isFinite(yVal)) pts.push([xVal, yVal, xArr[i]]);
  }
  if (pts.length < 3) return null;
  let n = pts.length;
  let sx = pts.reduce((a,p)=>a+p[0],0);
  let sy = pts.reduce((a,p)=>a+p[1],0);
  let mx = sx/n;
  let my = sy/n;
  let sxx = pts.reduce((a,p)=>a+Math.pow(p[0]-mx,2),0);
  let syy = pts.reduce((a,p)=>a+Math.pow(p[1]-my,2),0);
  if (sxx === 0 || syy === 0) return null;
  let sxy = pts.reduce((a,p)=>a+(p[0]-mx)*(p[1]-my),0);
  let slope = sxy/sxx;
  let intercept = my - slope*mx;
  let ssRes = pts.reduce((a,p)=>a+Math.pow(p[1]-(slope*p[0]+intercept),2),0);
  let r2 = 1 - ssRes/syy;
  let sorted = pts.slice().sort((a,b)=>a[0]-b[0]);
  let x1 = sorted[0][0];
  let x2 = sorted[sorted.length-1][0];
  return {slope:slope, intercept:intercept, r2:r2, n:n, x1:x1, x2:x2, y1:slope*x1+intercept, y2:slope*x2+intercept};
}

function meanAndSD(valuesArr) {
  const vals = valuesArr.map(Number).filter(v => !isNaN(v) && isFinite(v));
  if (vals.length === 0) return {mean:null, sd:null, n:0};
  const mean = vals.reduce((a,b)=>a+b,0) / vals.length;
  if (vals.length < 2) return {mean:mean, sd:null, n:vals.length};
  const variance = vals.reduce((a,b)=>a+Math.pow(b-mean,2),0) / (vals.length - 1);
  return {mean:mean, sd:Math.sqrt(variance), n:vals.length};
}

function correlationRFromFit(fit) {
  if (!fit) return null;
  const r = Math.sqrt(Math.max(0, fit.r2));
  return fit.slope < 0 ? -r : r;
}

function fmt(v) {
  if (v === null || v === undefined || isNaN(Number(v))) return 'NA';
  const n = Number(v);
  if (Math.abs(n) >= 1000 || Math.abs(n) < 0.01) return n.toExponential(2);
  return n.toFixed(3);
}


function outlierToggleOn() {
  if (isSatelliteLevelTable()) return satelliteShowOutliers && satelliteShowOutliers.checked;
  return perCellShowOutliers && perCellShowOutliers.checked;
}
function boolish(v) {
  return v === true || v === 1 || String(v).toLowerCase() === 'true' || String(v).toLowerCase() === 'yes';
}
function clearDynamicOutliers(rows) {
  (rows || []).forEach(r => {
    if (!r) return;
    r._Dynamic_Outlier_Cell = false;
    r._Dynamic_Outlier_Score = null;
    r._Dynamic_Outlier_Reason = '';
    r._Dynamic_Outlier_Label = '';
  });
}
function percentileSorted(sorted, p) {
  if (!sorted || sorted.length === 0) return NaN;
  if (sorted.length === 1) return sorted[0];
  const idx = (sorted.length - 1) * p;
  const lo = Math.floor(idx);
  const hi = Math.ceil(idx);
  if (lo === hi) return sorted[lo];
  return sorted[lo] + (sorted[hi] - sorted[lo]) * (idx - lo);
}
function outlierBaseLabel(row) {
  return (row && (row.Unique_Cell || row.Plot_Cell || row.Cell_ID || row.Std_Cell_ID)) ? String(row.Unique_Cell || row.Plot_Cell || ('Cell ' + (row.Cell_ID || row.Std_Cell_ID))) : 'Outlier point';
}
function computeDynamicOutliers(rows, metricName, valueGetter, groupGetter) {
  clearDynamicOutliers(rows);
  const groups = new Map();
  (rows || []).forEach(r => {
    let v = null;
    try { v = safeNumber(valueGetter(r)); } catch(e) { v = NaN; }
    if (isNaN(v) || !isFinite(v)) return;
    let g = 'visible data';
    try { g = String(groupGetter ? groupGetter(r) : 'visible data'); } catch(e) { g = 'visible data'; }
    if (!groups.has(g)) groups.set(g, []);
    groups.get(g).push({row:r, value:v});
  });

  groups.forEach((items, groupName) => {
    if (items.length < 5) return;
    const vals = items.map(d => d.value).sort((a,b)=>a-b);
    const q1 = percentileSorted(vals, 0.25);
    const med = percentileSorted(vals, 0.50);
    const q3 = percentileSorted(vals, 0.75);
    const iqr = q3 - q1;
    const absDev = vals.map(v => Math.abs(v - med)).sort((a,b)=>a-b);
    const mad = percentileSorted(absDev, 0.50);
    const meanVal = vals.reduce((a,b)=>a+b, 0) / vals.length;
    const sdVal = Math.sqrt(vals.reduce((a,b)=>a + Math.pow(b - meanVal, 2), 0) / Math.max(1, vals.length - 1));
    const low = q1 - 1.5 * iqr;
    const high = q3 + 1.5 * iqr;

    items.forEach(d => {
      const v = d.value;
      let reasons = [];
      let scores = [];
      if (iqr > 0 && (v < low || v > high)) {
        const direction = v > high ? 'high' : 'low';
        const iqrScore = v > high ? (v - high) / iqr : (low - v) / iqr;
        scores.push(Math.abs(iqrScore) + 1.5);
        reasons.push(direction + ' ' + metricName + ' for current graph set');
      }
      let rz = 0;
      if (mad > 0) rz = 0.6745 * (v - med) / mad;
      else if (sdVal > 0) rz = (v - meanVal) / sdVal;
      if (Math.abs(rz) >= 3) {
        scores.push(Math.abs(rz));
        reasons.push('robust z-score ' + rz.toFixed(2) + ' for ' + metricName);
      }
      if (reasons.length === 0) return;
      const score = scores.length ? Math.max(...scores) : Math.abs(rz);
      d.row._Dynamic_Outlier_Cell = true;
      d.row._Dynamic_Outlier_Score = score;
      d.row._Dynamic_Outlier_Label = outlierBaseLabel(d.row);
      d.row._Dynamic_Outlier_Reason = reasons.join('; ') + '. Set: ' + groupName + '. Value=' + fmt(v) + ', median=' + fmt(med) + ', IQR=' + fmt(iqr) + ', n=' + vals.length + '.';
    });
  });
  return (rows || []).filter(r => isOutlierRow(r));
}
function dynamicGroupForRow(row, xCol, table, graphType) {
  const parts = [];
  if (row && row.Mode) parts.push('mode ' + row.Mode);
  const xLooksCategorical = xCol && isCategoricalAxis(table, xCol);
  const xName = xCol ? String(xCol).toLowerCase() : '';
  const xIsIdLike = xName.includes('cell') || xName.includes('satellite') || xName.includes('source') || xName.includes('input') || xName.includes('uid');
  if (xLooksCategorical && !xIsIdLike && graphType !== 'histogram') {
    parts.push(label(table, xCol) + ' ' + String(row[xCol] || ''));
  } else if (row && row.Condition) {
    parts.push('condition ' + row.Condition);
  }
  if (!parts.length && row && row.Replicate) parts.push('replicate ' + row.Replicate);
  return parts.join(' / ') || 'current visible graph';
}
function isOutlierRow(row) {
  return row && boolish(row._Dynamic_Outlier_Cell);
}
function outlierLabel(row) {
  return (row && (row._Dynamic_Outlier_Label || row.Unique_Cell || row.Plot_Cell || row.Cell_ID || row.Std_Cell_ID)) ? String(row._Dynamic_Outlier_Label || row.Unique_Cell || row.Plot_Cell || ('Cell ' + (row.Cell_ID || row.Std_Cell_ID))) : 'Outlier point';
}
function outlierReason(row) {
  return row && row._Dynamic_Outlier_Reason ? String(row._Dynamic_Outlier_Reason) : 'Flagged dynamically for the current graph, current metric, and visible condition/set.';
}
function outlierScore(row) {
  const s = row ? Number(row._Dynamic_Outlier_Score) : NaN;
  return (!isNaN(s) && isFinite(s)) ? s.toFixed(2) : 'NA';
}
function outlierHoverSuffix(row) {
  return '<br><b style="color:#D62728">Current-graph outlier</b><br>Reason: ' + outlierReason(row) + '<br>Outlier score: ' + outlierScore(row);
}
function addOutlierTrace(traces, xVals, yVals, hoverText, nameSuffix) {
  if (!xVals || xVals.length === 0) return;
  traces.push({
    type:'scatter',
    mode:'markers',
    x:xVals,
    y:yVals,
    name:'Current graph outliers' + (nameSuffix ? ' — ' + nameSuffix : ''),
    marker:{size:16, color:'#D62728', opacity:0.96, symbol:'circle', line:{color:'#7A0000', width:1.4}},
    text:hoverText,
    hovertemplate:'%{text}<extra></extra>'
  });
}
function updateOutlierPanel(rows, contextText) {
  if (!outlierPanel || !outlierSummary) return;
  if (!outlierToggleOn()) {
    outlierPanel.style.display = 'none';
    outlierSummary.innerHTML = '';
    return;
  }
  const seen = new Set();
  const uniqueRows = [];
  (rows || []).forEach(r => {
    if (!isOutlierRow(r)) return;
    const key = outlierLabel(r) + '|' + outlierReason(r);
    if (seen.has(key)) return;
    seen.add(key);
    uniqueRows.push(r);
  });
  outlierPanel.style.display = '';
  if (uniqueRows.length === 0) {
    outlierSummary.innerHTML = '<div class="outlierhint">No outliers are visible for the current graph, metric, filters, and condition/set.</div>';
    return;
  }
  uniqueRows.sort((a,b) => {
    const sa = Number(a._Dynamic_Outlier_Score);
    const sb = Number(b._Dynamic_Outlier_Score);
    return (isNaN(sb) ? -1 : sb) - (isNaN(sa) ? -1 : sa);
  });
  const shown = uniqueRows.slice(0, 80);
  let html = '<div class="outlierhint">Large red dots are recalculated from the current graph only. Changing the metric, graph type, condition/filter, or analysis level changes which points are outliers.</div>';
  html += '<div class="outlierhint"><b>Visible current-graph outliers:</b> ' + uniqueRows.length + (contextText ? '<br>' + contextText : '') + '</div>';
  shown.forEach(r => {
    html += '<div class="outlieritem"><b>' + escapeHtml(outlierLabel(r)) + '</b><br>Score: ' + escapeHtml(outlierScore(r)) + '<br>Reason: ' + escapeHtml(outlierReason(r)) + '</div>';
  });
  if (uniqueRows.length > shown.length) {
    html += '<div class="outlierhint">Showing first ' + shown.length + ' outliers.</div>';
  }
  outlierSummary.innerHTML = html;
}

function escapeHtml(value) {
  return String(value === null || value === undefined ? '' : value)
    .replace(/&/g, '&amp;')
    .replace(/</g, '&lt;')
    .replace(/>/g, '&gt;')
    .replace(/"/g, '&quot;');
}


function borderColorForCell(i) {
  const palette = ['#7E3FB5','#3A86FF','#06A77D','#FF9F1C','#E76F51','#2A9D8F','#F15BB5','#00BBF9','#9B5DE5','#70C1B3','#F28482','#6D597A'];
  return palette[i % palette.length];
}

function pastelColorForCell(i) {
  const palette = ['#DCC6FF','#BDE0FE','#CDEAC0','#FFE5B4','#FFD6C2','#BEE3DB','#FBC4D8','#CAF0F8','#E0BBE4','#D8F3DC','#FAD2E1','#E9D8FD'];
  return palette[i % palette.length];
}

function purpleForCell(i) {
  return borderColorForCell(i);
}

function buildHoverForRow(row, xCol, yCol, indexInCell) {
  let txt = [];
  txt.push('Point: ' + indexInCell);
  if (row.Unique_Satellite) txt.push(String(row.Unique_Satellite).replace(/^.*\\|\\s*/, ''));
  if (row.Plot_Cell) txt.push(prettyCellLabel(row.Plot_Cell));
  if (row.Std_Frame !== null && row.Std_Frame !== undefined && row.Std_Frame !== '') txt.push('Frame: ' + row.Std_Frame);
  if (row.Std_Time_min !== null && row.Std_Time_min !== undefined && row.Std_Time_min !== '') txt.push('Time min: ' + row.Std_Time_min);
  if (xCol && row[xCol] !== null && row[xCol] !== undefined && row[xCol] !== '') txt.push(label(tableSelect.value, xCol) + ': ' + row[xCol]);
  txt.push(label(tableSelect.value, yCol) + ': ' + row[yCol]);
  return txt.join('<br>');
}

function addMeanAndSDTraces(traces, xref, yref, xStart, xEnd, stats, panelName, color, showLegend) {
  if (showMean.checked && stats.mean !== null) {
    traces.push({
      type:'scatter',
      mode:'lines',
      x:[xStart, xEnd],
      y:[stats.mean, stats.mean],
      xaxis:xref,
      yaxis:yref,
      name:'Mean',
      showlegend:showLegend,
      line:{color:color, width:2.5},
      hovertemplate:panelName + '<br>Mean=' + fmt(stats.mean) + '<extra></extra>'
    });
  }
  if (showSD.checked && stats.mean !== null && stats.sd !== null) {
    traces.push({
      type:'scatter',
      mode:'lines',
      x:[xStart, xEnd],
      y:[stats.mean + stats.sd, stats.mean + stats.sd],
      xaxis:xref,
      yaxis:yref,
      name:'Mean + SD',
      showlegend:showLegend,
      line:{color:color, width:1.7, dash:'dash'},
      hovertemplate:panelName + '<br>Mean+SD=' + fmt(stats.mean + stats.sd) + '<extra></extra>'
    });
    traces.push({
      type:'scatter',
      mode:'lines',
      x:[xStart, xEnd],
      y:[stats.mean - stats.sd, stats.mean - stats.sd],
      xaxis:xref,
      yaxis:yref,
      name:'Mean - SD',
      showlegend:showLegend,
      line:{color:color, width:1.7, dash:'dash'},
      hovertemplate:panelName + '<br>Mean-SD=' + fmt(stats.mean - stats.sd) + '<extra></extra>'
    });
  }
}

function drawCellPanels(table, info, rows, x, y) {
  const cellCol = getCellColumnForCurrentTable();
  if (!cellCol) {
    Plotly.newPlot('plot', [], {title:{text:'No cell column found for this table', font:{color:SATQ_PURPLE}}});
    return;
  }

  let cells = selectedCells();
  const allCells = uniqueNonEmpty(rows.map(r => r[cellCol]).map(v => String(v))).sort((a,b) => cellSortKey(a) - cellSortKey(b) || String(a).localeCompare(String(b)));
  if (cells.length === 0) cells = allCells;
  cells = cells.filter(c => allCells.includes(String(c)));

  if (cells.length === 0) {
    Plotly.newPlot('plot', [], {title:{text:'No cells selected / no cell data available', font:{color:SATQ_PURPLE}}});
    return;
  }

  const style = cellPlotType ? cellPlotType.value : 'point';
  const requestedCols = panelsPerRow ? Number(panelsPerRow.value) : 3;
  const maxCols = Math.max(1, Math.min(5, requestedCols));
  const nCols = Math.min(maxCols, cells.length);
  const nRows = Math.ceil(cells.length / nCols);
  const xGap = 0.045;
  const yGap = 0.15;
  let traces = [];
  let layout = {
    title: {text: (style === 'scatterxy' ? (label(table, y) + ' vs ' + label(table, x)) : label(table, y)) + ' — selected cells', font:{color:SATQ_PURPLE, size:20}},
    paper_bgcolor:'white',
    plot_bgcolor:'white',
    hovermode:'closest',
    showlegend:true,
    margin:{l:78, r:30, t:105, b:75},
    annotations:[]
  };

  cells.forEach((cell, i) => {
    const rr = rows.filter(r => String(r[cellCol]) === String(cell));
    const panelRow = Math.floor(i / nCols);
    const panelCol = i % nCols;
    const x0 = panelCol / nCols + xGap;
    const x1 = (panelCol + 1) / nCols - xGap;
    const y1 = 1 - panelRow / nRows - yGap / 2;
    const y0 = 1 - (panelRow + 1) / nRows + yGap / 2;

    const axisIndex = i + 1;
    const xa = axisIndex === 1 ? 'xaxis' : 'xaxis' + axisIndex;
    const ya = axisIndex === 1 ? 'yaxis' : 'yaxis' + axisIndex;
    const xref = axisIndex === 1 ? 'x' : 'x' + axisIndex;
    const yref = axisIndex === 1 ? 'y' : 'y' + axisIndex;

    const yy = rr.map(r => satelliteYValue(r, y, denom));
    const xxIndex = rr.map((r, j) => j + 1);
    const xxNum = numericValues(rr, x);
    const hover = rr.map((r, j) => buildHoverForRow(r, x, y, j + 1));
    const stats = meanAndSD(yy);
    const panelName = prettyCellLabel(cell);
    const cellColor = borderColorForCell(i);
    const lightColor = pastelColorForCell(i);

    const isScatterXY = style === 'scatterxy';
    const xTitle = isScatterXY ? label(table, x) : (style === 'box' ? '' : 'Satellite index');
    const isBottomRow = panelRow === nRows - 1;

    layout[xa] = {
      domain:[x0, x1],
      title:{text:isBottomRow ? xTitle : '', font:{size:10}},
      tickfont:{size:9},
      showticklabels:isBottomRow,
      zeroline:false,
      linecolor:'#444444',
      automargin:true
    };
    layout[ya] = {
      domain:[y0, y1],
      title:{text: panelCol === 0 ? label(table, y) : '', font:{size:10}},
      tickfont:{size:9},
      showticklabels:panelCol === 0,
      zeroline:false,
      linecolor:'#444444',
      automargin:true
    };

    if (style === 'bar') {
      traces.push({
        type:'bar',
        x:xxIndex,
        y:yy,
        xaxis:xref,
        yaxis:yref,
        name:panelName,
        showlegend:false,
        marker:{color:lightColor, line:{color:cellColor, width:1}},
        text:hover,
        hovertemplate:'%{text}<extra></extra>'
      });
      addMeanAndSDTraces(traces, xref, yref, 0.5, xxIndex.length + 0.5, stats, panelName, cellColor, i===0);
    } else if (style === 'profile') {
      traces.push({
        type:'scatter',
        mode:'lines+markers',
        x:xxIndex,
        y:yy,
        xaxis:xref,
        yaxis:yref,
        name:panelName,
        showlegend:false,
        line:{color:cellColor, width:2},
        marker:{size:7, opacity:0.88, color:lightColor, line:{color:cellColor, width:1}},
        text:hover,
        hovertemplate:'%{text}<extra></extra>'
      });
      addMeanAndSDTraces(traces, xref, yref, 0.5, xxIndex.length + 0.5, stats, panelName, cellColor, i===0);
    } else if (style === 'box') {
      const bx = yy.map(_ => panelName);
      traces.push({
        type:'box',
        x:bx,
        y:yy,
        xaxis:xref,
        yaxis:yref,
        name:panelName,
        showlegend:false,
        boxpoints:'all',
        jitter:0.35,
        pointpos:0,
        marker:{size:5, color:lightColor, line:{color:cellColor, width:0.8}},
        line:{color:cellColor},
        fillcolor:'rgba(166,107,218,0.28)',
        hovertemplate:panelName + '<br>' + label(table, y) + '=%{y}<extra></extra>'
      });
      if (showMean.checked && stats.mean !== null) {
        traces.push({
          type:'scatter',
          mode:'markers',
          x:[panelName],
          y:[stats.mean],
          xaxis:xref,
          yaxis:yref,
          name:'Mean',
          showlegend:i===0,
          marker:{symbol:'diamond', size:11, color:cellColor, line:{color:'#111111', width:1}},
          hovertemplate:panelName + '<br>Mean=' + fmt(stats.mean) + '<extra></extra>'
        });
      }
      if (showSD.checked && stats.mean !== null && stats.sd !== null) {
        traces.push({
          type:'scatter',
          mode:'markers',
          x:[panelName],
          y:[stats.mean],
          xaxis:xref,
          yaxis:yref,
          name:'Mean ± SD',
          showlegend:i===0,
          marker:{symbol:'diamond', size:9, color:cellColor},
          error_y:{type:'data', array:[stats.sd], visible:true, color:cellColor, thickness:1.5, width:8},
          hovertemplate:panelName + '<br>Mean=' + fmt(stats.mean) + '<br>SD=' + fmt(stats.sd) + '<extra></extra>'
        });
      }
    } else if (style === 'scatterxy') {
      traces.push({
        type:'scatter',
        mode:'markers',
        x:xxNum,
        y:yy,
        xaxis:xref,
        yaxis:yref,
        name:panelName,
        showlegend:false,
        marker:{size:8, opacity:0.86, color:lightColor, line:{color:cellColor, width:1}},
        text:hover,
        hovertemplate:'%{text}<extra></extra>'
      });
      const xValsClean = xxNum.map(Number).filter(v => !isNaN(v) && isFinite(v));
      if (xValsClean.length > 0) {
        addMeanAndSDTraces(traces, xref, yref, Math.min(...xValsClean), Math.max(...xValsClean), stats, panelName, cellColor, i===0);
      }
    } else {
      traces.push({
        type:'scatter',
        mode:'markers',
        x:xxIndex,
        y:yy,
        xaxis:xref,
        yaxis:yref,
        name:panelName,
        showlegend:false,
        marker:{size:8, opacity:0.86, color:lightColor, line:{color:cellColor, width:1}},
        text:hover,
        hovertemplate:'%{text}<extra></extra>'
      });
      addMeanAndSDTraces(traces, xref, yref, 0.5, xxIndex.length + 0.5, stats, panelName, cellColor, i===0);
    }

    let statText = '<b>' + panelName + '</b><br>n=' + rr.length;
    if (stats.mean !== null) statText += '; mean=' + fmt(stats.mean);
    if (stats.sd !== null) statText += '; SD=' + fmt(stats.sd);

    layout.annotations.push({
      text:statText,
      x:(x0+x1)/2,
      y:Math.min(0.995, y1 + 0.045),
      xref:'paper',
      yref:'paper',
      showarrow:false,
      align:'center',
      bgcolor:'rgba(255,255,255,0.78)',
      bordercolor:'rgba(220,198,255,0.65)',
      borderpad:2,
      font:{color:cellColor, size:10}
    });
  });

  const height = Math.max(620, 390 * nRows);
  document.getElementById('plot').style.height = height + 'px';
  Plotly.newPlot('plot', traces, layout, {responsive:true, toImageButtonOptions:{format:'png', filename:'SatelliteQ_cell_panels', scale:3}});
}

function drawTimeLapsePlot(table, info, rows, x, y) {
  // Time-lapse viewer: rows are frame summaries grouped by condition/replicate/frame.
  if (!rows || rows.length === 0) {
    document.getElementById('plot').style.height = '700px';
    Plotly.newPlot('plot', [], {title:{text:'No time-lapse frames match the selected filters', font:{color:SATQ_PURPLE}}, paper_bgcolor:'white', plot_bgcolor:'white'});
    return;
  }
  const plotType = typeSelect.value || 'line';
  const xCol = x || (info.columns.includes('Time_min') ? 'Time_min' : 'Frame');
  const yTitle = label(table, y);
  const xTitle = label(table, xCol);
  const conditionsVisible = uniqueNonEmpty(values(rows, 'Condition').map(v => String(v))).sort((a,b)=>a.localeCompare(b, undefined, {numeric:true}));
  const replicatesVisible = uniqueNonEmpty(values(rows, 'Replicate').map(v => String(v))).sort((a,b)=>a.localeCompare(b, undefined, {numeric:true}));
  let traces = [];
  let allX = [];
  let allY = [];

  conditionsVisible.forEach((cond) => {
    replicatesVisible.forEach((rep) => {
      let rr = rows.filter(r => String(r.Condition || '') === String(cond) && String(r.Replicate || '') === String(rep) && !isNaN(safeNumber(r[xCol])) && !isNaN(safeNumber(r[y])) && isFinite(safeNumber(r[xCol])) && isFinite(safeNumber(r[y])));
      if (rr.length === 0) return;
      rr = rr.slice().sort((a,b)=>safeNumber(a[xCol])-safeNumber(b[xCol]));
      const xx = rr.map(r => safeNumber(r[xCol]));
      const yy = rr.map(r => safeNumber(r[y]));
      allX = allX.concat(xx);
      allY = allY.concat(yy);
      const color = replicateShadeColor(cond, rep, conditionsVisible, replicatesVisible);
      const hover = rr.map((r) => {
        let txt = [];
        txt.push('Condition: ' + (r.Condition || ''));
        txt.push('Replicate: ' + (r.Replicate || ''));
        txt.push('Frame: ' + r.Frame);
        if (r.Time_min !== null && r.Time_min !== undefined && r.Time_min !== '') txt.push('Time min: ' + r.Time_min);
        txt.push(yTitle + ': ' + fmt(safeNumber(r[y])));
        if (r.Frame_file_count !== null && r.Frame_file_count !== undefined && r.Frame_file_count !== '') txt.push('Frame files combined: ' + r.Frame_file_count);
        return txt.join('<br>');
      });
      traces.push({
        type:'scatter',
        mode:plotType === 'line' ? 'lines+markers' : 'markers',
        x:xx,
        y:yy,
        name:String(cond) + ' / ' + String(rep),
        legendgroup:String(cond),
        line:{color:color, width:2},
        marker:{size:8, opacity:0.86, color:color, line:{color:'#333333', width:0.35}},
        text:hover,
        hovertemplate:'%{text}<extra></extra>'
      });
    });
  });

  let yMin = allY.length ? Math.min(...allY) : 0;
  let yMax = allY.length ? Math.max(...allY) : 1;
  let ySpan = Math.max(1e-9, yMax - yMin);
  let xMin = allX.length ? Math.min(...allX) : 0;
  let xMax = allX.length ? Math.max(...allX) : 1;
  let xSpan = Math.max(1e-9, xMax - xMin);
  document.getElementById('plot').style.height = '700px';
  Plotly.newPlot('plot', traces, {
    title:{text:yTitle + ' over time-lapse frames', font:{color:SATQ_PURPLE, size:19}},
    xaxis:{title:{text:xTitle, font:{color:'#333333'}}, range:[xMin - xSpan*0.05, xMax + xSpan*0.05], zeroline:false, linecolor:'#444444', tickfont:{color:'#333333'}, automargin:true},
    yaxis:{title:{text:yTitle, font:{color:'#333333'}}, range:[yMin - ySpan*0.12, yMax + ySpan*0.12], zeroline:false, linecolor:'#444444', tickfont:{color:'#333333'}, automargin:true},
    plot_bgcolor:'white',
    paper_bgcolor:'white',
    hovermode:'closest',
    legend:{title:{text:'Condition / replicate'}, orientation:'v', x:1.01, y:1, xanchor:'left', yanchor:'top', font:{size:11}},
    margin:{l:86, r:230, t:84, b:115}
  }, {responsive:true, toImageButtonOptions:{format:'png', filename:'SatelliteQ_timelapse_graph', scale:3}});
}

function drawPlot() {
  updateCellPanelVisibility();
  const table = tableSelect.value;
  const info = DATA[table];
  const rowsAll = info.rows;
  const type = typeSelect.value;
  let rows = rowsAfterSatelliteFilters(rowsAfterCellSelection(rowsAll));
  const x = xSelect.value;
  const y = ySelect.value;
  const group = '';

  if (isSatelliteLevelTable(table)) {
    drawSatellitePointPlot(table, info, rowsAfterSatelliteFilters(rowsAll), x, y, rowsAll);
    return;
  }

  if (isPerCellTable(table)) {
    drawPerCellPointPlot(table, info, rowsAfterSatelliteFilters(rowsAll), x, y);
    return;
  }

  if (isTimeLapseTable(table)) {
    drawTimeLapsePlot(table, info, rowsAfterSatelliteFilters(rowsAll), x, y);
    return;
  }

  if (type === 'cell_panels') {
    drawCellPanels(table, info, rowsAll, x, y);
    return;
  }

  const xCategorical = isCategoricalAxis(table, x) || type === 'bar' || type === 'box';
  let traces = [];
  const grouped = group ? uniqueNonEmpty(values(rows, group)) : [null];

  grouped.forEach((g, idx) => {
    const rr = group ? rowsForGroup(rows, group, g) : rows;
    const name = group ? String(g) : undefined;
    const c = group ? SATQ_PURPLE_PALETTE[idx % SATQ_PURPLE_PALETTE.length] : SATQ_PURPLE;
    const xsRaw = values(rr, x);
    const xsNum = numericValues(rr, x);
    const xsPlot = xCategorical ? xsRaw : xsNum;
    const ysNum = numericValues(rr, y);

    if (type === 'histogram') {
      traces.push({
        type:'histogram',
        x:numericValues(rr, x),
        name:name,
        opacity:0.80,
        marker:{color:SATQ_LIGHT_PURPLE, line:{color:SATQ_PURPLE, width:1}}
      });
    } else if (type === 'bar') {
      traces.push({
        type:'bar',
        x:xsRaw,
        y:ysNum,
        name:name,
        marker:{color:SATQ_LIGHT_PURPLE, line:{color:SATQ_PURPLE, width:1}}
      });
    } else if (type === 'line') {
      const pts = rr.map(r => [safeNumber(r[x]), r[x], safeNumber(r[y])]).filter(p => !isNaN(p[0]) && !isNaN(p[2])).sort((a,b)=>a[0]-b[0]);
      traces.push({
        type:'scatter',
        mode:'lines+markers',
        x:pts.map(p=>p[1]),
        y:pts.map(p=>p[2]),
        name:name,
        line:{color:c, width:2},
        marker:{color:SATQ_LIGHT_PURPLE, line:{color:c, width:1}, size:8}
      });
    } else if (type === 'box') {
      traces.push({
        type:'box',
        x:xsRaw,
        y:ysNum,
        name:name || label(table,y),
        boxpoints:'outliers',
        marker:{color:SATQ_LIGHT_PURPLE, line:{color:SATQ_PURPLE_DARK, width:1}},
        line:{color:SATQ_PURPLE}
      });
    } else {
      traces.push({
        type:'scatter',
        mode:'markers',
        x:xsPlot,
        y:ysNum,
        name:name,
        marker:{size:8, opacity:0.82, color:SATQ_LIGHT_PURPLE, line:{color:c, width:1}}
      });

      if (!xCategorical && (!isPerCellTable(table) || (perCellShowTrendline && perCellShowTrendline.checked))) {
        const fit = linearRegressionR2(xsNum, ysNum);
        if (fit) {
          const rVal = correlationRFromFit(fit);
          traces.push({
            type:'scatter',
            mode:'lines',
            x:[fit.x1, fit.x2],
            y:[fit.y1, fit.y2],
            name:(name ? name + ' trendline' : 'Trendline') + ' (R=' + (rVal === null ? 'NA' : rVal.toFixed(3)) + ', R²=' + fit.r2.toFixed(3) + ', n=' + fit.n + ')',
            line:{color:c, width:2, dash:'dash'},
            hovertemplate:'R=' + (rVal === null ? 'NA' : rVal.toFixed(3)) + '<br>R²=' + fit.r2.toFixed(3) + '<br>n=' + fit.n + '<extra></extra>'
          });
        }
      }
    }
  });


  if (outlierToggleOn()) {
    const outlierMetric = type === 'histogram' ? x : y;
    const outRows = computeDynamicOutliers(
      rows,
      label(table, outlierMetric),
      r => type === 'histogram' ? safeNumber(r[x]) : safeNumber(r[y]),
      r => dynamicGroupForRow(r, x, table, type)
    );
    let outX = [];
    let outY = [];
    let outText = [];
    outRows.forEach((r) => {
      let xv = xCategorical ? r[x] : safeNumber(r[x]);
      let yv = type === 'histogram' ? 0 : safeNumber(r[y]);
      if (type !== 'histogram' && (isNaN(yv) || !isFinite(yv))) return;
      if (!xCategorical && (isNaN(xv) || !isFinite(xv))) return;
      outX.push(xv);
      outY.push(yv);
      outText.push(
        (r.Plot_Cell ? prettyCellLabel(r.Plot_Cell) + '<br>' : '') +
        label(table, x) + ': ' + (r[x] === undefined ? '' : r[x]) + '<br>' +
        (type === 'histogram' ? label(table, x) + ': ' + fmt(safeNumber(r[x])) + '<br>' : label(table, y) + ': ' + fmt(yv) + '<br>') +
        outlierHoverSuffix(r)
      );
    });
    addOutlierTrace(traces, outX, outY, outText, 'visible graph');
    updateOutlierPanel(outRows, 'Outliers are recalculated from the current graph type, metric, axes, and filters.');
  } else {
    clearDynamicOutliers(rows);
    updateOutlierPanel([], '');
  }

  let yTitle = type === 'histogram' ? 'Number of rows / satellites' : label(table, y);
  let xTitle = label(table, x);
  let title = type.charAt(0).toUpperCase()+type.slice(1) + ': ' + (type === 'histogram' ? xTitle : yTitle + ' vs ' + xTitle);

  let annotations = [];
  let shapes = [];
  const perCellMeanSDOn = isPerCellTable(table) && perCellShowMeanSD && perCellShowMeanSD.checked;
  const perCellRValueOn = !isPerCellTable(table) || (perCellShowRValue && perCellShowRValue.checked);
  if (type === 'point' && !group && !xCategorical && perCellRValueOn) {
    const fitAll = linearRegressionR2(numericValues(rows, x), numericValues(rows, y));
    if (fitAll) {
      const rVal = correlationRFromFit(fitAll);
      annotations.push({
        text:'R = ' + (rVal === null ? 'NA' : rVal.toFixed(3)) + '<br>R² = ' + fitAll.r2.toFixed(3) + '<br>n = ' + fitAll.n,
        xref:'paper',
        yref:'paper',
        x:0.02,
        y:0.98,
        showarrow:false,
        align:'left',
        bgcolor:'rgba(255,255,255,0.9)',
        bordercolor:SATQ_LIGHT_PURPLE,
        font:{color:SATQ_PURPLE_DARK, size:13}
      });
    }
  }
  if (type === 'point' && perCellMeanSDOn) {
    const yStats = meanAndSd(numericValues(rows, y));
    const xValsClean = numericValues(rows, x).filter(v => v !== null && !isNaN(v) && isFinite(v));
    const xMin = xValsClean.length ? Math.min(...xValsClean) : 0;
    const xMax = xValsClean.length ? Math.max(...xValsClean) : 1;
    const xSpan = Math.max(1e-9, xMax - xMin);
    if (yStats.n > 0 && yStats.mean !== null) {
      shapes.push({type:'line', xref:'x', yref:'y', x0:xMin - xSpan*0.03, x1:xMax + xSpan*0.03, y0:yStats.mean, y1:yStats.mean, line:{color:'#111111', width:2}});
      if (yStats.sd !== null) {
        shapes.push({type:'line', xref:'x', yref:'y', x0:xMin - xSpan*0.03, x1:xMax + xSpan*0.03, y0:yStats.mean + yStats.sd, y1:yStats.mean + yStats.sd, line:{color:'#111111', width:1.4, dash:'dot'}});
        shapes.push({type:'line', xref:'x', yref:'y', x0:xMin - xSpan*0.03, x1:xMax + xSpan*0.03, y0:yStats.mean - yStats.sd, y1:yStats.mean - yStats.sd, line:{color:'#111111', width:1.4, dash:'dot'}});
      }
      annotations.push({text:'n = ' + yStats.n + '<br>mean = ' + fmt(yStats.mean) + (yStats.sd !== null ? '<br>SD = ' + fmt(yStats.sd) : ''), xref:'paper', yref:'paper', x:0.98, y:0.98, xanchor:'right', showarrow:false, align:'right', bgcolor:'rgba(255,255,255,0.9)', bordercolor:SATQ_LIGHT_PURPLE, font:{color:SATQ_PURPLE_DARK, size:12}});
    }
  }

  const xAxisConfig = {
    title: {text: xTitle, font: {color:'#333333'}},
    zeroline:false,
    linecolor:'#444444',
    tickfont:{color:'#333333'},
    automargin:true
  };

  if (xCategorical) {
    let categories = uniqueNonEmpty(values(rows, x).map(v => String(v)));
    xAxisConfig.type = 'category';
    xAxisConfig.categoryorder = 'array';
    xAxisConfig.categoryarray = categories;
    xAxisConfig.tickmode = 'array';
    xAxisConfig.tickvals = categories;
    xAxisConfig.ticktext = categories;
    xAxisConfig.tickangle = -45;
  }

  document.getElementById('plot').style.height = '620px';
  Plotly.newPlot('plot', traces, {
    title: {text: title, font: {color: SATQ_PURPLE, size: 20}},
    xaxis: xAxisConfig,
    yaxis: {title: {text: yTitle, font: {color:'#333333'}}, zeroline:false, linecolor:'#444444', tickfont:{color:'#333333'}, automargin:true},
    plot_bgcolor: 'white',
    paper_bgcolor: 'white',
    barmode: 'overlay',
    hovermode: 'closest',
    legend: {orientation:'h'},
    margin: {l:80, r:30, t:70, b:190},
    annotations: annotations,
    shapes: shapes
  }, {responsive:true, toImageButtonOptions:{format:'png', filename:'SatelliteQ_interactive_graph', scale:3}});
}

[tableSelect,typeSelect,xSelect,ySelect,groupSelect,nameDisplayMode,cellPlotType,panelsPerRow,showMean,showSD,satelliteShowMeanSD,satelliteShowOutliers,satelliteNormalizeCondition,perCellShowMeanSD,perCellShowOutliers,perCellShowTrendline,perCellShowRValue].filter(Boolean).forEach(el => el.addEventListener('change', el === tableSelect ? updateColumns : function(){ updateSatelliteNormalizeControl(); drawPlot(); }));
init();
</script>
</body>
</html>"""

    html = html_template.replace("__PLOTLYJS__", plotlyjs).replace("__DATA_JSON__", data_json).replace("__SATELLITEQ_VERSION__", SATELLITEQ_VERSION)

    out = outdir / "Graphs" / "Interactive_Graph_Explorer.html"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(html, encoding="utf-8")
    return out

def experiment_cell_dot_plot(cell_summary: pd.DataFrame, y: str, outdir: Path, title: str) -> Optional[Path]:
    """GraphPad-style condition plot: one point per cell, colored by replicate."""
    needed = {"Condition", "Replicate", y}
    if cell_summary.empty or not needed.issubset(set(cell_summary.columns)):
        return None

    d = cell_summary.copy()
    d[y] = pd.to_numeric(d[y], errors="coerce")
    d = d.dropna(subset=["Condition", "Replicate", y])
    if d.empty or d["Condition"].nunique() < 1:
        return None

    conditions = list(pd.unique(d["Condition"].astype(str)))
    replicates = list(pd.unique(d["Replicate"].astype(str)))
    if not conditions:
        return None

    palette = [
        "#7E3FB5", "#2A9D8F", "#E76F51", "#3A86FF", "#FF9F1C",
        "#06A77D", "#F15BB5", "#6D597A", "#00BBF9", "#B56576"
    ]

    fig_w = max(5.6, min(12.0, 1.1 * len(conditions) + 3.8))
    fig, ax = plt.subplots(figsize=(fig_w, 5.2))

    rng = np.random.default_rng(42)
    rep_offsets = {}
    if len(replicates) == 1:
        rep_offsets[replicates[0]] = 0
    else:
        for i, rep in enumerate(replicates):
            rep_offsets[rep] = np.linspace(-0.18, 0.18, len(replicates))[i]

    for ri, rep in enumerate(replicates):
        gr = d[d["Replicate"].astype(str) == rep]
        xs = []
        ys = []
        for _, row in gr.iterrows():
            cx = conditions.index(str(row["Condition"]))
            jitter = rng.normal(0, 0.025)
            xs.append(cx + rep_offsets[str(rep)] + jitter)
            ys.append(float(row[y]))
        ax.scatter(
            xs, ys,
            s=42,
            alpha=0.86,
            color=palette[ri % len(palette)],
            edgecolors="white",
            linewidths=0.6,
            label=str(rep),
            zorder=3,
        )

    for ci, cond in enumerate(conditions):
        vals = pd.to_numeric(d.loc[d["Condition"].astype(str) == cond, y], errors="coerce").dropna()
        if vals.empty:
            continue
        mean = float(vals.mean())
        sd = float(vals.std(ddof=1)) if len(vals) > 1 else 0.0
        ax.hlines(mean, ci - 0.27, ci + 0.27, color="black", linewidth=2.0, zorder=4)
        if len(vals) > 1:
            ax.vlines(ci, mean - sd, mean + sd, color="black", linewidth=1.4, zorder=4)
            ax.hlines([mean - sd, mean + sd], ci - 0.11, ci + 0.11, color="black", linewidth=1.4, zorder=4)

    ax.set_xticks(np.arange(len(conditions)))
    ax.set_xticklabels(conditions, rotation=45, ha="right")
    ax.set_xlabel("Condition")
    ax.set_ylabel(label_for(y, cell_summary))
    ax.set_title(title)
    graphpad_axes(ax)
    if len(replicates) <= 12:
        ax.legend(title="Replicate", frameon=False, fontsize=8, title_fontsize=9, loc="best")
    return save_fig(fig, outdir / f"{safe_name(title)}.png")

def make_experiment_explorer_html(cell_summary: pd.DataFrame, outdir: Path) -> Optional[Path]:
    """Interactive experiment viewer: cell-level points grouped by condition and colored by replicate."""
    if cell_summary.empty or "Condition" not in cell_summary.columns:
        return None
    try:
        import json
        import plotly.offline as po
    except Exception:
        return None

    _, cell_summary_with_outliers, _, _ = add_satelliteq_outlier_flags(pd.DataFrame(), cell_summary.copy(), None, None)
    df = sanitize_cilium_length_metrics(cell_summary_with_outliers.copy())
    candidate_metrics = [
        "Cell_area_2D", "Cell_volume_3D", "Cell_measure",
        "Satellite_count", "Centrosome_number_per_cell",
        "Mean_area_2D", "Median_area_2D", "Total_area_2D",
        "Mean_volume_3D", "Median_volume_3D", "Total_volume_3D",
        "Mean_size", "Median_size", "Total_size",
        "Mean_intensity", "Total_intensity", "Mean_distance", "Mean_sphericity", "Mean_circularity_2D", "Mean_shape_compactness",
        "Cell_spatial_entropy",
        "Cell_pericentrosomal_clustering_index",
        PCI_CELL_INNER_COL, PCI_CELL_OUTER_COL, PCI_CELL_TOTAL_COL,
        "Mean_distance_to_cell_or_image_center",
        "Cilium_count", "Cilium_length_um", "Mean_cilium_length_um", "Ciliated"
    ]
    metrics = [
        c for c in candidate_metrics
        if c in df.columns and (
            c in PERICENTROSOMAL_CELL_METRICS
            or pd.to_numeric(df[c], errors="coerce").notna().any()
        )
    ]
    if not metrics:
        return None

    outlier_cols = [c for c in ["Outlier_Cell", "Outlier_Score", "Outlier_Reason", "Outlier_Label"] if c in df.columns]
    keep_cols = ["Condition", "Replicate", "Source_Input", "Unique_Cell", "Cell_ID", "Frame", "Mode", "Ciliation_status"] + metrics + outlier_cols
    keep_cols = [c for c in keep_cols if c in df.columns]
    data = df[keep_cols].copy()
    for c in metrics:
        data[c] = pd.to_numeric(data[c], errors="coerce")

    labels = {c: label_for(c, cell_summary) for c in keep_cols}
    has_cilia = (
        "Cilium_length_um" in df.columns and pd.to_numeric(df["Cilium_length_um"], errors="coerce").notna().any()
    ) or (
        "Ciliated" in df.columns and pd.to_numeric(df["Ciliated"], errors="coerce").notna().any()
    )
    payload = {
        "rows": data.where(pd.notnull(data), None).to_dict(orient="records"),
        "metrics": metrics,
        "labels": labels,
        "hasCilia": bool(has_cilia),
    }

    plotlyjs = po.get_plotlyjs()
    data_json = json.dumps(payload, ensure_ascii=False)

    html_template = """<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>SatelliteQ Experiment Explorer</title>
<style>
body { font-family: Arial, sans-serif; margin: 0 auto; padding: 16px 18px; max-width: 1220px; background: #fbf8ff; color: #222; }
h1 { margin: 0 0 10px 0; color: #7E3FB5; font-size: 22px; }
.controls { display: grid; grid-template-columns: minmax(210px, 1.4fr) minmax(150px, 0.8fr) minmax(135px, 0.8fr) minmax(145px, 0.9fr) minmax(140px, 0.8fr) minmax(150px, 0.8fr) minmax(170px, 0.9fr); gap: 9px; background: white; border: 1px solid #DCC6FF; padding: 11px; border-radius: 10px; box-shadow: 0 2px 10px rgba(126,63,181,0.08); margin-top: 8px; max-width: 1220px; }
label { font-weight: bold; font-size: 13px; color: #4B216E; display:block; margin-bottom:4px; }
select { width: 100%; padding: 6px; }
.checklabel { display:flex; align-items:center; gap:7px; height:100%; margin:0; font-weight:bold; }
.checklabel input { width:auto; margin:0; }
.filterbox { background: white; border: 1px solid #DCC6FF; padding: 9px 11px; border-radius: 10px; margin-top: 9px; max-width: 1120px; }
.chips { display:flex; flex-wrap:wrap; gap:7px; margin-top:6px; }
.chip { border: 1px solid #DCC6FF; background:#fff; border-radius:999px; padding: 5px 9px; font-size:12px; }
.chip input { margin-right:4px; }
button { background:#7E3FB5; color:white; border:0; padding:7px 10px; border-radius:7px; cursor:pointer; margin-right:6px; }
.statbar { margin-top:9px; display:flex; gap:14px; flex-wrap:wrap; align-items:center; max-width:1220px; background:white; border:1px solid #DCC6FF; padding:10px 12px; border-radius:10px; box-shadow:0 2px 10px rgba(126,63,181,0.08); color:#4B216E; }
.statbar label { display:flex; align-items:center; gap:6px; margin:0; font-weight:bold; color:#4B216E; }
.statbar select { min-width:210px; padding:6px; }
.graphwrap { display:block; max-width:1220px; margin-top:12px; }
#plot { width:100%; min-width:0; height: 580px; background: white; border-radius: 10px; border: 1px solid #E8D9FF; }
#outlierPanel { margin-top:12px; width:auto; max-width:1220px; background:white; border:1px solid #DCC6FF; border-radius:10px; padding:12px; box-shadow:0 2px 10px rgba(126,63,181,0.08); display:none; }
#outlierPanel h3 { margin:0 0 8px 0; color:#7E3FB5; font-size:15px; }
.outlierlist { max-height:260px; overflow:auto; display:grid; grid-template-columns:repeat(auto-fit, minmax(260px, 1fr)); gap:8px; font-size:12px; line-height:1.35; color:#333; }
.outlieritem { border-left:4px solid #D62728; padding:7px 8px; margin:0; background:#fff6f6; border-radius:7px; }
.outlieritem b { color:#A00000; }
.outlierhint { color:#555; font-size:12px; margin-bottom:8px; grid-column:1 / -1; }
.legendbottom { margin-top:12px; width:auto; max-width:1220px; background:white; border:1px solid #DCC6FF; border-radius:10px; padding:12px; box-shadow:0 2px 10px rgba(126,63,181,0.08); }
.legendbottom h3 { margin:0 0 8px 0; color:#7E3FB5; font-size:15px; }
.legendhint { color:#555; font-size:12px; margin-bottom:8px; }
.legenditems { display:grid; grid-template-columns:repeat(auto-fit, minmax(260px, 1fr)); gap:8px; }
.legenditem { display:flex; align-items:flex-start; gap:8px; padding:7px 8px; background:#fbf8ff; border:1px solid #eadcff; border-radius:7px; font-size:12px; line-height:1.25; min-width:0; }
.legenddot { width:13px; height:13px; border-radius:50%; margin-top:2px; flex:0 0 13px; border:1px solid rgba(0,0,0,0.25); }
.legendtext { overflow-wrap:anywhere; word-break:break-word; }

</style>
<script>__PLOTLYJS__</script>
</head>
<body>
<h1>SatelliteQ Experiment Explorer</h1>
<div style="font-size:12px;color:#555;margin-bottom:6px;">SatelliteQ version __SATELLITEQ_VERSION__</div>
<div class="controls">
  <div><label>Y-axis metric</label><select id="metricSelect"></select></div>
  <div id="normalizeBox" style="display:none;"><label>Normalize to</label><select id="normalizeCondition"></select></div>
  <div><label>Color points by</label><select id="colorBy"><option value="Replicate">Replicate</option><option value="Condition">Condition</option></select></div>
  <div><label>Point display</label><select id="pointMode"><option value="all">All cells</option><option value="mean_only">Mean ± SD only</option></select></div>
  <div><label>Long names</label><select id="nameDisplayMode"><option value="short" selected>Short</option><option value="wrapped">Wrapped</option><option value="full">Full</option></select></div>
  <div><label class="checklabel"><input type="checkbox" id="showMean" checked> Show mean ± SD</label></div>
  <div><label class="checklabel"><input type="checkbox" id="showOutliers"> Show outlier cells/points</label></div>
</div>

<div class="filterbox">
  <b style="color:#4B216E;">Condition filter</b>
  <div id="conditionChips" class="chips"></div>
</div>
<div id="experimentSignificanceStatsControl" class="statbar">
  <b>Statistical comparison</b>
  <label><input type="checkbox" id="expShowSignificance" checked> Show significance on condition graph</label>
  <label>Compare each condition to<select id="expReferenceCondition"></select></label>
  <label>Statistical test<select id="expStatTestSelect"><option value="student_t" selected>Student t-test, unpaired</option><option value="paired_t">Paired t-test by replicate</option><option value="anova">One-way ANOVA, pairwise vs reference</option></select></label>
  <label>Use values<select id="expStatSampleLevel"><option value="replicate_means" selected>Replicate averages</option><option value="each_value">Each value</option></select></label>
  <label>Significance label<select id="expSignificanceLabelMode"><option value="stars_only" selected>Only stars</option><option value="full">Full parameters</option></select></label>
</div>
<div class="filterbox">
  <b style="color:#4B216E;">Replicate filter</b>
  <div id="replicateChips" class="chips"></div>
</div>

<div class="graphwrap">
  <div id="plot"></div>
</div>
<section id="legendPanel" class="legendbottom" style="display:none;">
  <h3>Graph legend</h3>
  <div id="legendSummary"></div>
</section>
<section id="outlierPanel">
  <h3>Outlier cells</h3>
  <div id="outlierSummary" class="outlierlist"></div>
</section>

<script>
const DATA = __DATA_JSON__;
const COLORS = ['#7E3FB5','#2A9D8F','#E76F51','#3A86FF','#FF9F1C','#06A77D','#F15BB5','#6D597A','#00BBF9','#B56576','#70C1B3','#F28482'];
const metricSelect = document.getElementById('metricSelect');
const normalizeBox = document.getElementById('normalizeBox');
const normalizeCondition = document.getElementById('normalizeCondition');
const colorBy = document.getElementById('colorBy');
const pointMode = document.getElementById('pointMode');
const nameDisplayMode = document.getElementById('nameDisplayMode');
const showMean = document.getElementById('showMean');
const showOutliers = document.getElementById('showOutliers');
const outlierPanel = document.getElementById('outlierPanel');
const outlierSummary = document.getElementById('outlierSummary');
const expShowSignificance = document.getElementById('expShowSignificance');
const expReferenceCondition = document.getElementById('expReferenceCondition');
const expStatTestSelect = document.getElementById('expStatTestSelect');
const expStatSampleLevel = document.getElementById('expStatSampleLevel');
const expSignificanceLabelMode = document.getElementById('expSignificanceLabelMode');

function unique(vals) { return [...new Set(vals.filter(v => v !== null && v !== undefined && String(v).trim() !== '').map(String))]; }
function label(col) { return DATA.labels[col] || col; }
function isNormalizedMetric(col) { return String(col || '').startsWith('NORM::'); }
function rawMetric(col) { return isNormalizedMetric(col) ? String(col).replace('NORM::', '') : col; }
function metricLabel(col) {
  const raw = rawMetric(col);
  if (!isNormalizedMetric(col)) return label(raw);
  const ref = normalizeCondition ? normalizeCondition.value : 'reference';
  return 'Normalized ' + label(raw) + ' (reference: ' + ref + ' = 1)';
}
function mean(vals) { return vals.reduce((a,b)=>a+b,0)/vals.length; }
function sd(vals) { if (vals.length < 2) return 0; const m=mean(vals); return Math.sqrt(vals.reduce((a,b)=>a+(b-m)*(b-m),0)/(vals.length-1)); }

function cleanDisplayName(value) {
  let s = String(value === null || value === undefined ? '' : value);
  const suffixes = ['_SatelliteQuantify_Fiji_Output.csv','_timelapse_satellite_rows_combined.csv','_SatelliteQ_Python_Report.xlsx','satellite_3D_results_clean.csv','satellite_3D_results_clean','_SatelliteQuantify_3D_satellites_only','SatelliteQuantify_3D_satellites_only','.csv','.xlsx','.xlsm'];
  suffixes.forEach(suf => { s = s.split(suf).join(''); });
  s = s.replace(/(?:^|[_\\s-])(?:20\\d{2}|19\\d{2})[_\\s-]?\\d{1,2}[_\\s-]?\\d{1,2}[_\\s-]?\\d{1,2}[_\\s-]?\\d{1,2}[_\\s-]?\\d{1,2}(?=$|[_\\s-])/g, '_');
  s = s.replace(/_{2,}/g, '_').replace(/\\s{2,}/g, ' ').trim();
  return s || String(value || '');
}
function shortenDisplayPart(value, maxLen=28) {
  let s = cleanDisplayName(value);
  if (s.length <= maxLen) return s;
  const head = Math.max(8, Math.floor(maxLen * 0.58));
  const tail = Math.max(6, maxLen - head - 1);
  return s.slice(0, head) + '…' + s.slice(s.length - tail);
}
function wrapOneDisplayLine(line, maxLen=30) {
  const s = String(line || '');
  if (s.length <= maxLen) return [s];
  const tokens = s.split(/([_\\s-]+)/);
  const out = [];
  let cur = '';
  tokens.forEach(tok => {
    if (!tok) return;
    if ((cur + tok).length <= maxLen) cur += tok;
    else {
      if (cur.trim()) out.push(cur.trim());
      if (tok.length > maxLen) { for (let i=0; i<tok.length; i+=maxLen) out.push(tok.slice(i, i+maxLen)); cur = ''; }
      else cur = tok;
    }
  });
  if (cur.trim()) out.push(cur.trim());
  return out.length ? out : [s];
}
function wrapDisplayName(value, maxLen=30, maxLines=5) {
  let s = cleanDisplayName(value);
  const roughLines = s.split(/\\s+\\/\\s+|\\s+—\\s+/);
  let lines = [];
  roughLines.forEach(part => { lines = lines.concat(wrapOneDisplayLine(part, maxLen)); });
  if (lines.length > maxLines) { lines = lines.slice(0, maxLines); lines[maxLines-1] = lines[maxLines-1].replace(/…?$/, '') + '…'; }
  return lines.join('<br>');
}
function displayName(value, mode=null) {
  const m = mode || (nameDisplayMode ? nameDisplayMode.value : 'short');
  if (m === 'full') return String(value === null || value === undefined ? '' : value);
  if (m === 'wrapped') return wrapDisplayName(value, 30, 5);
  const s = cleanDisplayName(value);
  const parts = s.split(/\\s+\\/\\s+/);
  if (parts.length > 1) return parts.map(p => shortenDisplayPart(p, 24)).join(' / ');
  return shortenDisplayPart(s, 42);
}
function displayTickLabel(value) {
  const m = nameDisplayMode ? nameDisplayMode.value : 'short';
  if (m === 'full') return String(value === null || value === undefined ? '' : value);
  if (m === 'wrapped') return wrapDisplayName(value, 22, 4);
  return shortenDisplayPart(value, 28);
}
function patchAxisTickLabels(axisObj) {
  if (!axisObj || !Array.isArray(axisObj.ticktext)) return axisObj;
  const copy = Object.assign({}, axisObj);
  copy.ticktext = axisObj.ticktext.map(displayTickLabel);
  return copy;
}
function patchLayoutDisplayNames(layout) {
  if (!layout || !nameDisplayMode || nameDisplayMode.value === 'full') return layout;
  const out = Object.assign({}, layout);
  ['xaxis','yaxis','xaxis2','yaxis2','xaxis3','yaxis3','xaxis4','yaxis4'].forEach(k => { if (out[k]) out[k] = patchAxisTickLabels(out[k]); });
  return out;
}
function patchTraceDisplayNames(traces) {
  if (!Array.isArray(traces) || !nameDisplayMode || nameDisplayMode.value === 'full') return traces;
  return traces.map(tr => {
    if (!tr || tr.name === undefined || tr.name === null) return tr;
    const t = Object.assign({}, tr);
    t.name = displayName(tr.name);
    return t;
  });
}

function legendEscapeHtml(value) {
  return String(value === null || value === undefined ? '' : value)
    .replace(/&/g, '&amp;')
    .replace(/</g, '&lt;')
    .replace(/>/g, '&gt;')
    .replace(/"/g, '&quot;');
}
function legendTraceColor(trace, index) {
  let c = null;
  if (trace && trace.marker && trace.marker.color !== undefined) c = trace.marker.color;
  if ((c === null || c === undefined) && trace && trace.line && trace.line.color !== undefined) c = trace.line.color;
  if (Array.isArray(c)) c = c.find(v => typeof v === 'string') || null;
  if (typeof c !== 'string' || c.length === 0) c = COLORS[index % COLORS.length] || '#7E3FB5';
  return c;
}
function updateCustomLegend(traces, layout) {
  if (!legendPanel || !legendSummary) return;
  const entries = [];
  const seen = new Set();
  (traces || []).forEach((tr, i) => {
    if (!tr || tr.showlegend === false) return;
    if (tr.name === undefined || tr.name === null || String(tr.name).trim() === '') return;
    const name = String(tr.name);
    if (name.toLowerCase() === 'undefined' || name.toLowerCase() === 'null') return;
    const color = legendTraceColor(tr, i);
    const key = name + '|' + color;
    if (seen.has(key)) return;
    seen.add(key);
    entries.push({name:name, color:color});
  });
  if (!entries.length) {
    legendPanel.style.display = 'none';
    legendSummary.innerHTML = '';
    return;
  }
  let title = 'Graph legend';
  try {
    if (layout && layout.legend && layout.legend.title) {
      title = typeof layout.legend.title === 'string' ? layout.legend.title : (layout.legend.title.text || title);
    }
  } catch(e) {}
  legendPanel.style.display = '';
  let html = '<div class="legendhint"><b>' + legendEscapeHtml(title) + '</b> — colors and dot labels for the current graph.</div>';
  html += '<div class="legenditems">';
  entries.forEach(e => {
    html += '<div class="legenditem"><span class="legenddot" style="background:' + legendEscapeHtml(e.color) + '"></span><span class="legendtext">' + legendEscapeHtml(e.name) + '</span></div>';
  });
  html += '</div>';
  legendSummary.innerHTML = html;
}

const SatelliteQOriginalNewPlot = Plotly.newPlot.bind(Plotly);
Plotly.newPlot = function(div, traces, layout, config) {
  const patchedTraces = patchTraceDisplayNames(traces);
  const patchedLayout = patchLayoutDisplayNames(layout);
  updateCustomLegend(patchedTraces, patchedLayout);
  let layoutNoBuiltInLegend = patchedLayout;
  if (legendPanel && patchedLayout) layoutNoBuiltInLegend = Object.assign({}, patchedLayout, {showlegend:false});
  return SatelliteQOriginalNewPlot(div, patchedTraces, layoutNoBuiltInLegend, config);
};
function checkedValues(containerId) {
  return [...document.querySelectorAll('#'+containerId+' input[type=checkbox]')].filter(x=>x.checked).map(x=>x.value);
}
function setAll(containerId, state) {
  document.querySelectorAll('#'+containerId+' input[type=checkbox]').forEach(x => x.checked = state);
  render();
}
function initChips(containerId, values) {
  const el = document.getElementById(containerId);
  el.innerHTML = '';
  values.forEach(v => {
    const lab = document.createElement('label');
    lab.className = 'chip';
    lab.innerHTML = `<input type="checkbox" value="${String(v).replace(/"/g,'&quot;')}" checked>${v}`;
    lab.querySelector('input').addEventListener('change', render);
    el.appendChild(lab);
  });
  const controls = document.createElement('span');
  controls.innerHTML = `<button onclick="setAll('${containerId}', true)">All</button><button onclick="setAll('${containerId}', false)">None</button>`;
  el.appendChild(controls);
}
function updateNormalizeControl() {
  if (!normalizeBox || !normalizeCondition) return;
  if (!isNormalizedMetric(metricSelect.value)) {
    normalizeBox.style.display = 'none';
    return;
  }
  const prev = normalizeCondition.value;
  const conds = unique(DATA.rows.map(r => r.Condition)).sort((a,b)=>a.localeCompare(b, undefined, {numeric:true}));
  normalizeCondition.innerHTML = '';
  conds.forEach(c => {
    const o = document.createElement('option');
    o.value = c;
    o.textContent = c;
    normalizeCondition.appendChild(o);
  });
  if (conds.includes(prev)) normalizeCondition.value = prev;
  else if (conds.includes('WT')) normalizeCondition.value = 'WT';
  else if (conds.includes('Wild Type')) normalizeCondition.value = 'Wild Type';
  else if (conds.length) normalizeCondition.value = conds[0];
  normalizeBox.style.display = '';
}
function normalizationDenominators(rawCol, refCondition) {
  const refRows = DATA.rows.filter(r => String(r.Condition || '') === String(refCondition));
  const reps = unique(refRows.map(r => r.Replicate));
  const byRep = {};
  reps.forEach(rep => {
    const vals = refRows.filter(r => String(r.Replicate || '') === String(rep)).map(r => { const q=r[rawCol]; return (q===null||q===undefined||q===''||String(q).toUpperCase()==='NA') ? NaN : Number(q); }).filter(v => !isNaN(v) && isFinite(v));
    if (vals.length > 0) byRep[String(rep)] = mean(vals);
  });
  const allVals = refRows.map(r => { const q=r[rawCol]; return (q===null||q===undefined||q===''||String(q).toUpperCase()==='NA') ? NaN : Number(q); }).filter(v => !isNaN(v) && isFinite(v));
  return {byRep: byRep, global: allVals.length ? mean(allVals) : null};
}
function metricValue(row, col, denom) {
  const raw = rawMetric(col);
  const rawValue = row[raw];
  if (rawValue === null || rawValue === undefined || rawValue === '' || String(rawValue).toUpperCase() === 'NA') return null;
  const v = Number(rawValue);
  if (isNaN(v) || !isFinite(v)) return null;
  if (!isNormalizedMetric(col)) return v;
  const rep = String(row.Replicate || '');
  const d = denom && denom.byRep && denom.byRep[rep] ? denom.byRep[rep] : (denom ? denom.global : null);
  if (!d || isNaN(d) || !isFinite(d) || d === 0) return null;
  return v / d;
}

function boolish(v) {
  return v === true || v === 1 || String(v).toLowerCase() === 'true' || String(v).toLowerCase() === 'yes';
}
function clearDynamicOutliers(rows) {
  (rows || []).forEach(r => {
    if (!r) return;
    r._Dynamic_Outlier_Cell = false;
    r._Dynamic_Outlier_Score = null;
    r._Dynamic_Outlier_Reason = '';
    r._Dynamic_Outlier_Label = '';
  });
}
function percentileSorted(sorted, p) {
  if (!sorted || sorted.length === 0) return NaN;
  if (sorted.length === 1) return sorted[0];
  const idx = (sorted.length - 1) * p;
  const lo = Math.floor(idx);
  const hi = Math.ceil(idx);
  if (lo === hi) return sorted[lo];
  return sorted[lo] + (sorted[hi] - sorted[lo]) * (idx - lo);
}
function outlierBaseLabel(row) {
  return (row && (row.Unique_Cell || row.Cell_ID)) ? String(row.Unique_Cell || ('Cell ' + row.Cell_ID)) : 'Outlier cell';
}
function computeDynamicOutliers(rows, metricName, valueGetter, groupGetter) {
  clearDynamicOutliers(rows);
  const groups = new Map();
  (rows || []).forEach(r => {
    let v = null;
    try { v = Number(valueGetter(r)); } catch(e) { v = NaN; }
    if (isNaN(v) || !isFinite(v)) return;
    let g = 'visible data';
    try { g = String(groupGetter ? groupGetter(r) : 'visible data'); } catch(e) { g = 'visible data'; }
    if (!groups.has(g)) groups.set(g, []);
    groups.get(g).push({row:r, value:v});
  });
  groups.forEach((items, groupName) => {
    if (items.length < 5) return;
    const vals = items.map(d => d.value).sort((a,b)=>a-b);
    const q1 = percentileSorted(vals, 0.25);
    const med = percentileSorted(vals, 0.50);
    const q3 = percentileSorted(vals, 0.75);
    const iqr = q3 - q1;
    const absDev = vals.map(v => Math.abs(v - med)).sort((a,b)=>a-b);
    const mad = percentileSorted(absDev, 0.50);
    const meanVal = vals.reduce((a,b)=>a+b, 0) / vals.length;
    const sdVal = Math.sqrt(vals.reduce((a,b)=>a + Math.pow(b - meanVal, 2), 0) / Math.max(1, vals.length - 1));
    const low = q1 - 1.5 * iqr;
    const high = q3 + 1.5 * iqr;
    items.forEach(d => {
      const v = d.value;
      let reasons = [];
      let scores = [];
      if (iqr > 0 && (v < low || v > high)) {
        const direction = v > high ? 'high' : 'low';
        const iqrScore = v > high ? (v - high) / iqr : (low - v) / iqr;
        scores.push(Math.abs(iqrScore) + 1.5);
        reasons.push(direction + ' ' + metricName + ' for current graph set');
      }
      let rz = 0;
      if (mad > 0) rz = 0.6745 * (v - med) / mad;
      else if (sdVal > 0) rz = (v - meanVal) / sdVal;
      if (Math.abs(rz) >= 3) {
        scores.push(Math.abs(rz));
        reasons.push('robust z-score ' + rz.toFixed(2) + ' for ' + metricName);
      }
      if (reasons.length === 0) return;
      d.row._Dynamic_Outlier_Cell = true;
      d.row._Dynamic_Outlier_Score = scores.length ? Math.max(...scores) : Math.abs(rz);
      d.row._Dynamic_Outlier_Label = outlierBaseLabel(d.row);
      d.row._Dynamic_Outlier_Reason = reasons.join('; ') + '. Set: ' + groupName + '. Value=' + fmt(v) + ', median=' + fmt(med) + ', IQR=' + fmt(iqr) + ', n=' + vals.length + '.';
    });
  });
  return (rows || []).filter(r => isOutlierRow(r));
}
function experimentOutlierGroup(row) {
  const parts = [];
  if (row && row.Condition) parts.push('condition ' + row.Condition);
  if (row && row.Mode) parts.push('mode ' + row.Mode);
  return parts.join(' / ') || 'current visible experiment';
}
function isOutlierRow(row) {
  return row && boolish(row._Dynamic_Outlier_Cell);
}
function outlierLabel(row) {
  return (row && (row._Dynamic_Outlier_Label || row.Unique_Cell || row.Cell_ID)) ? String(row._Dynamic_Outlier_Label || row.Unique_Cell || ('Cell ' + row.Cell_ID)) : 'Outlier cell';
}
function outlierReason(row) {
  return row && row._Dynamic_Outlier_Reason ? String(row._Dynamic_Outlier_Reason) : 'Flagged dynamically for the current metric and visible condition/set.';
}
function outlierScore(row) {
  const s = row ? Number(row._Dynamic_Outlier_Score) : NaN;
  return (!isNaN(s) && isFinite(s)) ? s.toFixed(2) : 'NA';
}
function escapeHtml(value) {
  return String(value === null || value === undefined ? '' : value)
    .replace(/&/g, '&amp;')
    .replace(/</g, '&lt;')
    .replace(/>/g, '&gt;')
    .replace(/"/g, '&quot;');
}
function outlierHoverSuffix(row) {
  return '<br><b style="color:#D62728">Current-metric outlier</b><br>Reason: ' + outlierReason(row) + '<br>Outlier score: ' + outlierScore(row);
}
function updateOutlierPanel(rows) {
  if (!outlierPanel || !outlierSummary) return;
  if (!showOutliers || !showOutliers.checked) {
    outlierPanel.style.display = 'none';
    outlierSummary.innerHTML = '';
    return;
  }
  const seen = new Set();
  const uniqueRows = [];
  (rows || []).forEach(r => {
    if (!isOutlierRow(r)) return;
    const key = outlierLabel(r) + '|' + outlierReason(r);
    if (seen.has(key)) return;
    seen.add(key);
    uniqueRows.push(r);
  });
  outlierPanel.style.display = 'block';
  if (uniqueRows.length === 0) {
    outlierSummary.innerHTML = '<div class="outlierhint">No outliers are visible for the current metric, filters, and condition/set.</div>';
    return;
  }
  uniqueRows.sort((a,b) => {
    const sa = Number(a._Dynamic_Outlier_Score);
    const sb = Number(b._Dynamic_Outlier_Score);
    return (isNaN(sb) ? -1 : sb) - (isNaN(sa) ? -1 : sa);
  });
  let html = '<div class="outlierhint">Large red dots are recalculated from the current Experiment Explorer metric and visible condition/set.</div>';
  html += '<div class="outlierhint"><b>Visible current-metric outliers:</b> ' + uniqueRows.length + '</div>';
  uniqueRows.slice(0, 80).forEach(r => {
    html += '<div class="outlieritem"><b>' + escapeHtml(outlierLabel(r)) + '</b><br>Score: ' + escapeHtml(outlierScore(r)) + '<br>Reason: ' + escapeHtml(outlierReason(r)) + '</div>';
  });
  if (uniqueRows.length > 80) html += '<div class="outlierhint">Showing first 80 outliers.</div>';
  outlierSummary.innerHTML = html;
}

function init() {
  if (DATA.hasCilia && colorBy && !Array.from(colorBy.options).some(o => o.value === 'Ciliation_status')) {
    const copt = document.createElement('option');
    copt.value = 'Ciliation_status';
    copt.textContent = 'Ciliation status';
    colorBy.appendChild(copt);
  }
  DATA.metrics.forEach(m => {
    const o = document.createElement('option');
    o.value = m;
    o.textContent = label(m);
    metricSelect.appendChild(o);
  });
  const normGroup = document.createElement('optgroup');
  normGroup.label = 'Normalized value';
  DATA.metrics.forEach(m => {
    if (m === 'Ciliated') return;
    const o = document.createElement('option');
    o.value = 'NORM::' + m;
    o.textContent = label(m) + ' — normalized';
    normGroup.appendChild(o);
  });
  metricSelect.appendChild(normGroup);
  if (DATA.metrics.includes('Total_size')) metricSelect.value = 'Total_size';
  else if (DATA.metrics.includes('Satellite_count')) metricSelect.value = 'Satellite_count';

  initChips('conditionChips', unique(DATA.rows.map(r => r.Condition)));
  initChips('replicateChips', unique(DATA.rows.map(r => r.Replicate)));
  [metricSelect, normalizeCondition, colorBy, pointMode, nameDisplayMode, showMean, showOutliers, expShowSignificance, expReferenceCondition, expStatTestSelect, expStatSampleLevel, expSignificanceLabelMode].filter(Boolean).forEach(el => el.addEventListener('change', render));
  updateNormalizeControl();
  render();
}
function jitterFor(index, groupIndex) {
  const s = Math.sin((index + 1) * 12.9898 + groupIndex * 78.233) * 43758.5453;
  return (s - Math.floor(s) - 0.5) * 0.16;
}

function expStarLabel(p) {
  if (p === null || isNaN(p) || !isFinite(p)) return '';
  if (p < 0.0001) return '****';
  if (p < 0.001) return '***';
  if (p < 0.01) return '**';
  if (p < 0.05) return '*';
  return 'ns';
}
function expPLabel(p) {
  if (p === null || isNaN(p) || !isFinite(p)) return 'P = NA';
  if (p < 0.0001) return 'P < 0.0001';
  return 'P = ' + p.toPrecision(3);
}
function expLogGamma(z) {
  const g = 7;
  const p = [0.99999999999980993,676.5203681218851,-1259.1392167224028,771.32342877765313,-176.61502916214059,12.507343278686905,-0.13857109526572012,9.9843695780195716e-6,1.5056327351493116e-7];
  if (z < 0.5) return Math.log(Math.PI) - Math.log(Math.sin(Math.PI * z)) - expLogGamma(1 - z);
  z -= 1;
  let x = p[0];
  for (let i = 1; i < p.length; i++) x += p[i] / (z + i);
  const t = z + g + 0.5;
  return 0.5 * Math.log(2 * Math.PI) + (z + 0.5) * Math.log(t) - t + Math.log(x);
}
function expBetaCf(a,b,x) {
  const MAXIT=100, EPS=3e-7, FPMIN=1e-30;
  let qab=a+b, qap=a+1, qam=a-1, c=1, d=1-qab*x/qap;
  if (Math.abs(d)<FPMIN) d=FPMIN;
  d=1/d;
  let h=d;
  for (let m=1; m<=MAXIT; m++) {
    const m2=2*m;
    let aa=m*(b-m)*x/((qam+m2)*(a+m2));
    d=1+aa*d; if (Math.abs(d)<FPMIN) d=FPMIN;
    c=1+aa/c; if (Math.abs(c)<FPMIN) c=FPMIN;
    d=1/d; h*=d*c;
    aa=-(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
    d=1+aa*d; if (Math.abs(d)<FPMIN) d=FPMIN;
    c=1+aa/c; if (Math.abs(c)<FPMIN) c=FPMIN;
    d=1/d; const del=d*c; h*=del;
    if (Math.abs(del-1.0)<EPS) break;
  }
  return h;
}
function expIBeta(x,a,b) {
  if (x<=0) return 0;
  if (x>=1) return 1;
  const bt=Math.exp(expLogGamma(a+b)-expLogGamma(a)-expLogGamma(b)+a*Math.log(x)+b*Math.log(1-x));
  if (x < (a+1)/(a+b+2)) return bt*expBetaCf(a,b,x)/a;
  return 1 - bt*expBetaCf(b,a,1-x)/b;
}
function expStudentTCdf(t,df) {
  if (df<=0 || !isFinite(t)) return NaN;
  const x=df/(df+t*t);
  const ib=expIBeta(x, df/2, 0.5);
  return t>=0 ? 1 - 0.5*ib : 0.5*ib;
}
function expTPValue(t,df) {
  const cdf=expStudentTCdf(Math.abs(t), df);
  if (isNaN(cdf)) return NaN;
  return Math.max(0, Math.min(1, 2*(1-cdf)));
}
function expFPValue(f,df1,df2) {
  if (!isFinite(f) || f<0 || df1<=0 || df2<=0) return NaN;
  const x=(df1*f)/(df1*f+df2);
  return Math.max(0, Math.min(1, 1-expIBeta(x, df1/2, df2/2)));
}
function expMeanVar(vals) {
  const a=vals.map(Number).filter(v => !isNaN(v) && isFinite(v));
  const n=a.length;
  if (n===0) return {n:0, mean:NaN, variance:NaN};
  const mean=a.reduce((s,v)=>s+v,0)/n;
  const variance=n>1 ? a.reduce((s,v)=>s+Math.pow(v-mean,2),0)/(n-1) : 0;
  return {n:n, mean:mean, variance:variance};
}
function expGroupValues(rows, conditions, yGetter) {
  return conditions.map(c => rows.filter(r => String(r.Condition)===String(c)).map(r => Number(yGetter(r))).filter(v=>!isNaN(v)&&isFinite(v)));
}
function expGroupReplicateMeans(rows, conditions, yGetter) {
  return conditions.map(c => {
    const byRep={};
    rows.filter(r => String(r.Condition)===String(c)).forEach((r, idx) => {
      const v=Number(yGetter(r)); if (isNaN(v)||!isFinite(v)) return;
      const rep = r.Replicate !== undefined && r.Replicate !== null && String(r.Replicate).trim() !== '' ? String(r.Replicate) : ('row_'+idx);
      if (!byRep[rep]) byRep[rep]=[];
      byRep[rep].push(v);
    });
    return Object.keys(byRep).map(k => byRep[k].reduce((a,b)=>a+b,0)/byRep[k].length);
  });
}
function expUnpairedT(a,b,welch=true) {
  const A=expMeanVar(a), B=expMeanVar(b);
  if (A.n<2 || B.n<2) return {p:NaN, test:welch?'Welch t-test':'Student t-test'};
  if (welch) {
    const se2=A.variance/A.n + B.variance/B.n;
    if (se2<=0) return {p:NaN, test:'Welch t-test'};
    const t=(A.mean-B.mean)/Math.sqrt(se2);
    const df=Math.pow(se2,2)/(Math.pow(A.variance/A.n,2)/(A.n-1)+Math.pow(B.variance/B.n,2)/(B.n-1));
    return {p:expTPValue(t, df), test:'Welch t-test'};
  }
  const pooled=((A.n-1)*A.variance+(B.n-1)*B.variance)/(A.n+B.n-2);
  if (pooled<=0) return {p:NaN, test:'Student t-test'};
  const t=(A.mean-B.mean)/Math.sqrt(pooled*(1/A.n+1/B.n));
  return {p:expTPValue(t, A.n+B.n-2), test:'Student t-test'};
}
function expPairedT(rows,conditions,yGetter) {
  if (conditions.length!==2) return {p:NaN, test:'Paired t-test'};
  const maps={}; conditions.forEach(c => maps[String(c)]={});
  rows.forEach((r, idx) => {
    const c=String(r.Condition); if (!(c in maps)) return;
    const y=Number(yGetter(r)); if (isNaN(y)||!isFinite(y)) return;
    const rep = r.Replicate !== undefined && r.Replicate !== null && String(r.Replicate).trim() !== '' ? String(r.Replicate) : ('row_'+idx);
    if (!maps[c][rep]) maps[c][rep]=[];
    maps[c][rep].push(y);
  });
  const reps=Object.keys(maps[String(conditions[0])] || {}).filter(k => maps[String(conditions[1])] && maps[String(conditions[1])][k]);
  const diffs=reps.map(k => {
    const a=maps[String(conditions[0])][k].reduce((s,v)=>s+v,0)/maps[String(conditions[0])][k].length;
    const b=maps[String(conditions[1])][k].reduce((s,v)=>s+v,0)/maps[String(conditions[1])][k].length;
    return a-b;
  });
  const D=expMeanVar(diffs);
  if (D.n<2 || D.variance<=0) return {p:NaN, test:'Paired t-test'};
  return {p:expTPValue(D.mean/Math.sqrt(D.variance/D.n), D.n-1), test:'Paired t-test'};
}
function expAnova(samples) {
  const groups=samples.map(s => s.map(Number).filter(v=>!isNaN(v)&&isFinite(v))).filter(s=>s.length>0);
  const k=groups.length, nTotal=groups.reduce((s,g)=>s+g.length,0);
  if (k<2 || nTotal<=k) return {p:NaN, test:'One-way ANOVA'};
  const all=[].concat(...groups);
  const grand=all.reduce((s,v)=>s+v,0)/all.length;
  let ssb=0, ssw=0;
  groups.forEach(g => {
    const m=g.reduce((s,v)=>s+v,0)/g.length;
    ssb += g.length*Math.pow(m-grand,2);
    g.forEach(v => ssw += Math.pow(v-m,2));
  });
  const df1=k-1, df2=nTotal-k, msb=ssb/df1, msw=ssw/df2;
  if (msw<=0) return {p:NaN, test:'One-way ANOVA'};
  return {p:expFPValue(msb/msw, df1, df2), test:'One-way ANOVA'};
}
function updateExperimentReferenceConditionOptions(conditions) {
  if (!expReferenceCondition) return;
  const previous = expReferenceCondition.value;
  expReferenceCondition.innerHTML = '';
  conditions.forEach(c => {
    const opt = document.createElement('option');
    opt.value = String(c);
    opt.textContent = String(c);
    expReferenceCondition.appendChild(opt);
  });
  if (conditions.includes(previous)) expReferenceCondition.value = previous;
  else if (conditions.length > 0) expReferenceCondition.value = String(conditions[0]);
}

function expSamplesForTwoConditions(rows, conditionA, conditionB, yGetter) {
  const selected = [String(conditionA), String(conditionB)];
  if (expStatSampleLevel && expStatSampleLevel.value === 'each_value') {
    return selected.map(c => rows.filter(r => String(r.Condition) === c).map(r => Number(yGetter(r))).filter(v => !isNaN(v) && isFinite(v)));
  }
  return selected.map(c => {
    const byRep = {};
    rows.filter(r => String(r.Condition) === c).forEach((r, idx) => {
      const v = Number(yGetter(r));
      if (isNaN(v) || !isFinite(v)) return;
      const rep = r.Replicate !== undefined && r.Replicate !== null && String(r.Replicate).trim() !== '' ? String(r.Replicate) : ('row_' + idx);
      if (!byRep[rep]) byRep[rep] = [];
      byRep[rep].push(v);
    });
    return Object.keys(byRep).map(k => byRep[k].reduce((a,b)=>a+b,0) / byRep[k].length);
  });
}

function expPairwiseAnova(a, b) {
  return expAnova([a, b]);
}

function addExperimentSignificance(shapes, annotations, rows, conditions, yGetter, yMin, yMax, ySpan) {
  if (!expShowSignificance || !expShowSignificance.checked || !conditions || conditions.length < 2) return null;
  updateExperimentReferenceConditionOptions(conditions);
  const ref = expReferenceCondition && expReferenceCondition.value ? String(expReferenceCondition.value) : String(conditions[0]);
  if (!conditions.includes(ref)) return null;

  const test = expStatTestSelect ? expStatTestSelect.value : 'student_t';
  const sampleTxt = expStatSampleLevel && expStatSampleLevel.value === 'each_value' ? 'each value' : 'replicate averages';
  const refIndex = conditions.indexOf(ref);
  let topY = yMax;
  let level = 0;

  conditions.forEach((cond, idx) => {
    cond = String(cond);
    if (cond === ref) return;

    let result = null;
    if (test === 'paired_t') {
      result = expPairedT(rows, [ref, cond], yGetter);
    } else {
      const samples = expSamplesForTwoConditions(rows, ref, cond, yGetter);
      if (test === 'anova') result = expPairwiseAnova(samples[0], samples[1]);
      else result = expUnpairedT(samples[0], samples[1], false);
    }
    if (!result) return;

    const x0 = refIndex;
    const x1 = idx;
    const xa = Math.min(x0, x1);
    const xb = Math.max(x0, x1);
    const y = yMax + ySpan * (0.14 + level * 0.13);
    const h = ySpan * 0.04;

    shapes.push({type:'line', xref:'x', yref:'y', x0:xa, x1:xa, y0:y, y1:y+h, line:{color:'#111111', width:1.6}});
    shapes.push({type:'line', xref:'x', yref:'y', x0:xa, x1:xb, y0:y+h, y1:y+h, line:{color:'#111111', width:1.6}});
    shapes.push({type:'line', xref:'x', yref:'y', x0:xb, x1:xb, y0:y+h, y1:y, line:{color:'#111111', width:1.6}});

    const stars = expStarLabel(result.p);
    const labelMode = expSignificanceLabelMode ? expSignificanceLabelMode.value : 'stars_only';
    const labelText = labelMode === 'full'
      ? (stars + '<br><span style="font-size:10px">' + ref + ' vs ' + cond + '<br>' + expPLabel(result.p) + '<br>' + result.test + ', ' + sampleTxt + '</span>')
      : stars;
    annotations.push({
      x:(xa+xb)/2,
      y:y+h*1.7,
      text:labelText,
      showarrow:false,
      align:'center',
      font:{color:'#111111', size: labelMode === 'full' ? 12 : 16},
      bgcolor: labelMode === 'full' ? 'rgba(255,255,255,0.90)' : 'rgba(255,255,255,0)',
      bordercolor: labelMode === 'full' ? '#111111' : 'rgba(255,255,255,0)',
      borderpad: labelMode === 'full' ? 3 : 0
    });

    topY = Math.max(topY, y + h * 2.8);
    level += 1;
  });
  return level > 0 ? topY : null;
}
function render() {
  updateNormalizeControl();
  const metric = metricSelect.value;
  const raw = rawMetric(metric);
  const denom = isNormalizedMetric(metric) ? normalizationDenominators(raw, normalizeCondition.value) : null;
  const conds = checkedValues('conditionChips');
  const reps = checkedValues('replicateChips');
  let rows = DATA.rows.filter(r => conds.includes(String(r.Condition)) && reps.includes(String(r.Replicate)));
  rows = rows.filter(r => metricValue(r, metric, denom) !== null);
  const conditions = conds.filter(c => rows.some(r => String(r.Condition) === String(c)));
  updateExperimentReferenceConditionOptions(conditions);
  const colorField = colorBy.value;
  const groups = unique(rows.map(r => r[colorField]));
  const traces = [];

  if (pointMode.value === 'all') {
    groups.forEach((g, gi) => {
      const rr = rows.filter(r => String(r[colorField]) === String(g));
      const x = rr.map((r, i) => conditions.indexOf(String(r.Condition)) + jitterFor(i, gi));
      const y = rr.map(r => metricValue(r, metric, denom));
      const txt = rr.map(r => `${r.Unique_Cell || ('Cell '+r.Cell_ID)}<br>Condition: ${r.Condition}<br>Replicate: ${r.Replicate}<br>${metricLabel(metric)}: ${metricValue(r, metric, denom).toFixed(3)}`);
      traces.push({
        type:'scatter',
        mode:'markers',
        x:x,
        y:y,
        name:String(g),
        marker:{size:7, color:COLORS[gi % COLORS.length], opacity:0.82, line:{color:'white', width:0.7}},
        text:txt,
        hovertemplate:'%{text}<extra></extra>'
      });
    });
  }


  if (showOutliers && showOutliers.checked) {
    const outRows = computeDynamicOutliers(
      rows,
      metricLabel(metric),
      r => metricValue(r, metric, denom),
      r => experimentOutlierGroup(r)
    ).filter(r => metricValue(r, metric, denom) !== null);
    const outX = outRows.map((r, i) => conditions.indexOf(String(r.Condition)) + jitterFor(i, 777));
    const outY = outRows.map(r => metricValue(r, metric, denom));
    const outText = outRows.map(r =>
      `${r.Unique_Cell || ('Cell '+r.Cell_ID)}<br>Condition: ${r.Condition}<br>Replicate: ${r.Replicate}<br>${metricLabel(metric)}: ${metricValue(r, metric, denom).toFixed(3)}` + outlierHoverSuffix(r)
    );
    if (outRows.length > 0) {
      traces.push({
        type:'scatter',
        mode:'markers',
        x:outX,
        y:outY,
        name:'Current-metric outliers',
        marker:{size:16, color:'#D62728', opacity:0.96, symbol:'circle', line:{color:'#7A0000', width:1.4}},
        text:outText,
        hovertemplate:'%{text}<extra></extra>'
      });
    }
    updateOutlierPanel(outRows);
  } else {
    clearDynamicOutliers(rows);
    updateOutlierPanel([]);
  }

  if (showMean.checked) {
    conditions.forEach((cond, ci) => {
      const vals = rows.filter(r => String(r.Condition) === String(cond)).map(r => metricValue(r, metric, denom)).filter(v => v !== null && !isNaN(v));
      if (vals.length === 0) return;
      const m = mean(vals);
      const s = sd(vals);
      traces.push({type:'scatter', mode:'lines', x:[ci-0.22, ci+0.22], y:[m,m], line:{color:'black', width:3}, showlegend:false, hovertemplate:`${cond}<br>Mean=${m.toFixed(3)}<extra></extra>`});
      if (vals.length > 1) {
        traces.push({type:'scatter', mode:'lines', x:[ci,ci], y:[m-s,m+s], line:{color:'black', width:1.6}, showlegend:false, hovertemplate:`${cond}<br>Mean±SD=${m.toFixed(3)} ± ${s.toFixed(3)}<extra></extra>`});
        traces.push({type:'scatter', mode:'lines', x:[ci-0.08,ci+0.08], y:[m-s,m-s], line:{color:'black', width:1.6}, showlegend:false, hoverinfo:'skip'});
        traces.push({type:'scatter', mode:'lines', x:[ci-0.08,ci+0.08], y:[m+s,m+s], line:{color:'black', width:1.6}, showlegend:false, hoverinfo:'skip'});
      }
    });
  }

  let shapes = [];
  let annotations = [];
  const yVals = rows.map(r => metricValue(r, metric, denom)).filter(v => v !== null && !isNaN(v) && isFinite(v));
  let yMin = yVals.length ? Math.min(...yVals) : 0;
  let yMax = yVals.length ? Math.max(...yVals) : 1;
  let ySpan = Math.max(1e-9, yMax - yMin);
  const sigY = addExperimentSignificance(shapes, annotations, rows, conditions, r => metricValue(r, metric, denom), yMin, yMax, ySpan);
  if (sigY !== null) yMax = Math.max(yMax, sigY);
  ySpan = Math.max(1e-9, yMax - yMin);

  const layout = {
    title:{text: metricLabel(metric) + ' by condition', font:{color:'#7E3FB5', size:19}},
    xaxis:{tickmode:'array', tickvals:conditions.map((_,i)=>i), ticktext:conditions, title:'Condition', zeroline:false, linecolor:'#222'},
    yaxis:{title:metricLabel(metric), range:[yMin-ySpan*0.12, yMax+ySpan*0.14], zeroline:false, linecolor:'#222'},
    paper_bgcolor:'white',
    plot_bgcolor:'white',
    margin:{l:68,r:35,t:92,b:95},
    shapes:shapes,
    annotations:annotations,
    hovermode:'closest',
    legend:{title:{text: colorField, font:{size:12}}}
  };
  Plotly.newPlot('plot', traces, layout, {responsive:true, toImageButtonOptions:{format:'png', filename:'SatelliteQ_experiment_plot', scale:3}});
}
init();
</script>
</body>
</html>
"""
    html = html_template.replace("__PLOTLYJS__", plotlyjs).replace("__DATA_JSON__", data_json).replace("__SATELLITEQ_VERSION__", SATELLITEQ_VERSION)
    out = outdir / "Graphs" / "Experiment_Explorer.html"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(html, encoding="utf-8")
    return out

def make_experiment_graphs(cell_summary: pd.DataFrame, outdir: Path) -> List[Path]:
    paths: List[Path] = []
    if cell_summary.empty or "Condition" not in cell_summary.columns:
        return paths
    exp_dir = outdir / "Graphs" / "Experiment_conditions"
    metrics = [
        ("Total_size", "Total satellite/granule size per cell by condition"),
        ("Satellite_count", "Satellite count per cell by condition"),
        ("Mean_size", "Mean satellite size per cell by condition"),
        ("Total_intensity", "Total satellite integrated intensity per cell by condition"),
        ("Mean_intensity", "Mean satellite integrated intensity per cell by condition"),
        ("Mean_distance", "Mean satellite-to-centrosome distance per cell by condition"),
        ("Cell_measure", "Cell area or volume by condition"),
        ("Cell_spatial_entropy", "Cell spatial entropy by condition"),
        ("Cell_pericentrosomal_clustering_index", "Pericentrosomal clustering index by condition"),
        ("Mean_distance_to_cell_or_image_center", "Mean satellite distance to cell/image center by condition"),
        ("Cilium_length_um", "Cilium length per cell by condition"),
        ("Mean_cilium_length_um", "Mean cilium length per cell by condition"),
        ("Cilium_count", "Cilium count per cell by condition"),
        ("Ciliated", "Ciliation status per cell by condition"),
    ]
    for col, title in metrics:
        p = experiment_cell_dot_plot(cell_summary, col, exp_dir, title)
        if p is not None:
            paths.append(p)
    html = make_experiment_explorer_html(cell_summary, outdir)
    if html is not None:
        paths.append(html)
    return paths

def make_tabbed_graph_dashboard(outdir: Path) -> Optional[Path]:
    """Create one HTML file with two tabs: general viewer and experiment viewer."""
    graphs_dir = outdir / "Graphs"
    general = graphs_dir / "Interactive_Graph_Explorer.html"
    experiment = graphs_dir / "Experiment_Explorer.html"
    if not general.exists() and not experiment.exists():
        return None

    def _read_srcdoc(p: Path) -> str:
        if not p.exists():
            return "<html><body><h2>Viewer not available</h2></body></html>"
        return p.read_text(encoding="utf-8", errors="replace")

    general_doc = html_lib.escape(_read_srcdoc(general), quote=True)
    exp_doc = html_lib.escape(_read_srcdoc(experiment), quote=True)

    html = f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>SatelliteQ Graph Viewers</title>
<style>
body {{ margin:0; font-family: Arial, sans-serif; background:#fbf8ff; color:#222; }}
.header {{ padding:12px 22px 8px 22px; background:white; border-bottom:1px solid #DCC6FF; }}
h1 {{ margin:0; color:#7E3FB5; font-size:21px; }}
.tabs {{ display:flex; gap:8px; padding:9px 22px; background:#F6EFFF; border-bottom:1px solid #DCC6FF; }}
.tabbtn {{ border:1px solid #7E3FB5; background:white; color:#4B216E; padding:9px 14px; border-radius:9px; cursor:pointer; font-weight:bold; }}
.tabbtn.active {{ background:#7E3FB5; color:white; }}
.viewer {{ width:100%; height:calc(100vh - 96px); border:0; display:none; background:white; }}
.viewer.active {{ display:block; }}
</style>
</head>
<body>
<div class="header">
  <h1>SatelliteQ Graph Viewers</h1>
</div>
<div class="tabs">
  <button class="tabbtn active" id="btnGeneral" onclick="showTab('general')">General graph explorer</button>
  <button class="tabbtn" id="btnExperiment" onclick="showTab('experiment')">Experiment viewer</button>
</div>
<iframe id="frameGeneral" class="viewer active" srcdoc="{general_doc}"></iframe>
<iframe id="frameExperiment" class="viewer" srcdoc="{exp_doc}"></iframe>
<script>
function showTab(which) {{
  const gen = document.getElementById('frameGeneral');
  const exp = document.getElementById('frameExperiment');
  const bgen = document.getElementById('btnGeneral');
  const bexp = document.getElementById('btnExperiment');
  if (which === 'general') {{
    gen.classList.add('active'); exp.classList.remove('active');
    bgen.classList.add('active'); bexp.classList.remove('active');
  }} else {{
    exp.classList.add('active'); gen.classList.remove('active');
    bexp.classList.add('active'); bgen.classList.remove('active');
  }}
}}
</script>
</body>
</html>"""
    out = graphs_dir / satelliteq_dated_name(outdir, "Graph_Viewer_Tabs", ext=".html")
    out.write_text(html, encoding="utf-8")
    return out

def make_graphs(df: pd.DataFrame, image_summary: pd.DataFrame, cell_summary: pd.DataFrame, frame_summary: pd.DataFrame, outdir: Path) -> List[Path]:
    paths: List[Path] = []
    gdir = outdir / "Graphs"
    rel = gdir / "Relationships"
    dist = gdir / "Distributions"
    group = gdir / "Per_image_and_cell"
    time = gdir / "Time_lapse"

    for p in [
        hist(df, "Std_Area_2D", dist, "2D satellite area distribution"),
        hist(df, "Std_Volume_3D", dist, "3D satellite volume distribution"),
        hist(df, "Std_Size", dist, "Satellite size distribution"),
        hist(df, "Std_Intensity", dist, "Satellite integrated intensity distribution"),
        hist(df, "Std_Distance", dist, "Satellite-to-centrosome distance distribution"),
        hist(df, "Std_Sphericity", dist, "3D satellite sphericity distribution"),
        hist(df, "Std_Circularity_2D", dist, "2D satellite circularity distribution"),
        hist(df, "Std_Shape_Compactness", dist, "Satellite shape compactness distribution"),
        hist(cell_summary, "Cell_spatial_entropy", dist, "Cell spatial entropy distribution"),
        hist(cell_summary, "Cell_pericentrosomal_clustering_index", dist, "Pericentrosomal clustering index distribution"),
        hist(cell_summary, "Mean_distance_to_cell_or_image_center", dist, "Mean satellite distance to cell or image center distribution"),
        hist(df, "Std_Cell_Measure", dist, "Cell area or volume distribution"),
        bar(image_summary, "Source_Input", "Satellite_count", group, "Number of satellites per image"),
        bar(image_summary, "Source_Input", "Mean_area_2D", group, "Mean 2D satellite area per image"),
        bar(image_summary, "Source_Input", "Mean_volume_3D", group, "Mean 3D satellite volume per image"),
        bar(image_summary, "Source_Input", "Mean_size", group, "Mean satellite size per image"),
        bar(image_summary, "Source_Input", "Mean_intensity", group, "Mean satellite integrated intensity per image"),
        bar(image_summary, "Source_Input", "Mean_sphericity", group, "Mean 3D satellite sphericity per image"),
        bar(image_summary, "Source_Input", "Mean_circularity_2D", group, "Mean 2D satellite circularity per image"),
        bar(image_summary, "Source_Input", "Mean_shape_compactness", group, "Mean satellite shape compactness per image"),
    ]:
        if p is not None:
            paths.append(p)

    if not cell_summary.empty:
        for p in [
            scatter(cell_summary, "Cell_area_2D", "Mean_area_2D", rel, "Mean 2D satellite area vs 2D cell area"),
            scatter(cell_summary, "Cell_volume_3D", "Mean_volume_3D", rel, "Mean 3D satellite volume vs 3D cell volume"),
            scatter(cell_summary, "Cell_measure", "Satellite_count", rel, "Number of satellites vs cell area or volume"),
            scatter(cell_summary, "Cell_measure", "Total_intensity", rel, "Total satellite integrated intensity vs cell area or volume"),
            scatter(cell_summary, "Cell_measure", "Mean_size", rel, "Mean satellite size vs cell area or volume"),
            scatter(cell_summary, "Cell_measure", "Mean_distance", rel, "Mean satellite-to-centrosome distance vs cell area or volume"),
            scatter(cell_summary, "Cell_measure", "Mean_sphericity", rel, "Mean 3D satellite sphericity vs cell area or volume"),
            scatter(cell_summary, "Cell_measure", "Mean_circularity_2D", rel, "Mean 2D satellite circularity vs cell area or volume"),
            scatter(cell_summary, "Cell_measure", "Cell_spatial_entropy", rel, "Cell spatial entropy vs cell area or volume"),
            scatter(cell_summary, "Cell_measure", "Cell_pericentrosomal_clustering_index", rel, "Pericentrosomal clustering index vs cell area or volume"),
            scatter(cell_summary, "Cell_measure", "Mean_distance_to_cell_or_image_center", rel, "Mean satellite distance to cell or image center vs cell area or volume"),
            bar(cell_summary, "Cell_ID", "Satellite_count", group, "Number of satellites per cell"),
            bar(cell_summary, "Cell_ID", "Mean_sphericity", group, "Mean 3D satellite sphericity per cell"),
            bar(cell_summary, "Cell_ID", "Mean_circularity_2D", group, "Mean 2D satellite circularity per cell"),
            bar(cell_summary, "Ciliation_status", "Satellite_count", group, "Satellite number in ciliated and non-ciliated cells"),
            scatter(cell_summary, "Cilium_length_um", "Satellite_count", rel, "Satellite number vs cilium length"),
            scatter(cell_summary, "Cilium_length_um", "Cell_spatial_entropy", rel, "Cell spatial entropy vs cilium length"),
            scatter(cell_summary, "Cilium_length_um", "Cell_pericentrosomal_clustering_index", rel, "Pericentrosomal clustering index vs cilium length"),
        ]:
            if p is not None:
                paths.append(p)

    for p in [
        scatter(df, "Std_Area_2D", "Std_Intensity", rel, "Satellite integrated intensity vs 2D area"),
        scatter(df, "Std_Volume_3D", "Std_Intensity", rel, "Satellite integrated intensity vs 3D volume"),
        scatter(df, "Std_Size", "Std_Intensity", rel, "Satellite integrated intensity vs size"),
        scatter(df, "Std_Cell_Measure", "Std_Size", rel, "Satellite size vs cell area or volume"),
        scatter(df, "Std_Cell_Measure", "Std_Intensity", rel, "Satellite integrated intensity vs cell area or volume"),
        scatter(df, "Std_Distance", "Std_Size", rel, "Satellite size vs distance to centrosome"),
    ]:
        if p is not None:
            paths.append(p)

    if not frame_summary.empty:
        x = "Time_min" if "Time_min" in frame_summary.columns and pd.to_numeric(frame_summary["Time_min"], errors="coerce").notna().any() else "Frame"
        for p in [
            line(frame_summary, x, "Satellite_count", "Source_Input", time, "Number of satellites over time"),
            line(frame_summary, x, "Mean_area_2D", "Source_Input", time, "Mean 2D satellite area over time"),
            line(frame_summary, x, "Mean_volume_3D", "Source_Input", time, "Mean 3D satellite volume over time"),
            line(frame_summary, x, "Mean_size", "Source_Input", time, "Mean satellite size over time"),
            line(frame_summary, x, "Mean_intensity", "Source_Input", time, "Mean satellite integrated intensity over time"),
            line(frame_summary, x, "Mean_sphericity", "Source_Input", time, "Mean 3D satellite sphericity over time"),
            line(frame_summary, x, "Mean_circularity_2D", "Source_Input", time, "Mean 2D satellite circularity over time"),
            line(frame_summary, x, "Mean_shape_compactness", "Source_Input", time, "Mean satellite shape compactness over time"),
            line(frame_summary, x, "Mean_distance", "Source_Input", time, "Mean satellite-to-centrosome distance over time"),
            line(frame_summary, x, "Mean_cell_spatial_entropy", "Source_Input", time, "Mean cell spatial entropy over time"),
            line(frame_summary, x, "Mean_cell_pericentrosomal_clustering_index", "Source_Input", time, "Mean pericentrosomal clustering index over time"),
            line(frame_summary, x, "Mean_distance_to_cell_or_image_center", "Source_Input", time, "Mean satellite distance to cell or image center over time"),
        ]:
            if p is not None:
                paths.append(p)

    experiment_paths = make_experiment_graphs(cell_summary, outdir)
    for p in experiment_paths:
        if p is not None and p.name != "Experiment_Explorer.html":
            paths.append(p)

    html_path = make_interactive_html(df, image_summary, cell_summary, frame_summary, outdir)

    dashboard = make_tabbed_graph_dashboard(outdir)
    if dashboard is not None:
        paths.append(dashboard)

        for tmp in [html_path, outdir / "Graphs" / "Experiment_Explorer.html"]:
            try:
                if tmp is not None and Path(tmp).exists():
                    Path(tmp).unlink()
            except Exception:
                pass
    else:
        if html_path is not None:
            paths.append(html_path)
        for p in experiment_paths:
            if p is not None and p.name == "Experiment_Explorer.html":
                paths.append(p)

    return paths

def _viewer_numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")

def _viewer_mode_for_subset(df_src: pd.DataFrame) -> str:
    vals = df_src.get("Mode", pd.Series(dtype=object)).dropna().astype(str).unique().tolist()
    if vals:
        return vals[0]
    if pd.to_numeric(df_src.get("Std_Z", pd.Series(dtype=float)), errors="coerce").notna().any():
        return "3D"
    return "2D"

def _viewer_time_for_subset(df_src: pd.DataFrame) -> bool:
    fr = pd.to_numeric(df_src.get("Std_Frame", pd.Series(dtype=float)), errors="coerce")
    return fr.notna().sum() > 0 and fr.dropna().nunique() > 1

def make_cilium_viewer(df: pd.DataFrame, cilium_rows: pd.DataFrame, skeleton_points: pd.DataFrame, outdir: Path) -> Optional[Path]:
    """Create a self-contained interactive cilium/satellite viewer."""
    if cilium_rows is None or cilium_rows.empty: return None
    try:
        import json
        import plotly.offline as po
    except Exception:
        return None
    vdir=outdir/"Spatial_Viewers"; vdir.mkdir(parents=True,exist_ok=True)
    sat_cols=[c for c in ["Source_Input","Condition","Replicate","Std_Cell_ID","Std_Frame","Std_X","Std_Y","Std_Z","Std_Satellite_ID"] if c in df.columns]
    sat=df[sat_cols].copy() if sat_cols else pd.DataFrame()
    cilia=cilium_rows.copy()
    # Normalize common cilium coordinate names for the JS viewer.
    ren={}
    for dst,cands in {"Cell_ID":["Cell_ID","Std_Cell_ID"],"Frame":["Frame","Std_Frame"],"X":["Cilium_X_um","X_um","X"],"Y":["Cilium_Y_um","Y_um","Y"],"Z":["Cilium_Z_um","Z_um","Z"],"Length":["Cilium_Length_um","Length_um","Length"],"Cilium_ID":["Cilium_ID","ID"]}.items():
        src=first_existing(list(cilia.columns),cands)
        if src and src!=dst: ren[src]=dst
    cilia=cilia.rename(columns=ren)
    for c in ["Cell_ID","Frame","X","Y","Z","Length","Cilium_ID"]:
        if c not in cilia.columns: cilia[c]=np.nan
    sat_records=sat.where(pd.notna(sat),None).to_dict("records") if not sat.empty else []
    cil_records=cilia.where(pd.notna(cilia),None).to_dict("records")
    sk_records=skeleton_points.where(pd.notna(skeleton_points),None).to_dict("records") if skeleton_points is not None and not skeleton_points.empty else []
    plotly_js=po.get_plotlyjs()
    html=f"""<!doctype html><html><head><meta charset="utf-8"><title>SatelliteQ Cilium Viewer</title><script>{plotly_js}</script>
<style>body{{font-family:Arial,sans-serif;margin:18px}} .controls{{display:flex;gap:12px;flex-wrap:wrap;margin-bottom:10px}} select{{padding:6px}} #plot{{height:760px}}</style></head><body>
<h2>SatelliteQ Cilium Viewer</h2><div class="controls"><label>Cell <select id="cell"></select></label><label>Frame <select id="frame"></select></label></div><div id="plot"></div>
<script>const SAT={json.dumps(sat_records)}; const CIL={json.dumps(cil_records)}; const SK={json.dumps(sk_records)};
const cell=document.getElementById('cell'), frame=document.getElementById('frame');
function norm(v){{if(v===null||v===undefined||String(v)==='nan')return ''; let n=Number(v); return Number.isFinite(n)&&Math.floor(n)===n?String(n):String(v)}}
function vals(rows,k){{return [...new Set(rows.map(r=>norm(r[k])).filter(Boolean))].sort((a,b)=>a.localeCompare(b,undefined,{{numeric:true}}))}}
function fill(sel,arr){{sel.innerHTML='<option value="">All</option>'; arr.forEach(v=>{{let o=document.createElement('option');o.value=v;o.textContent=v;sel.appendChild(o)}})}}
fill(cell, vals(CIL,'Cell_ID')); fill(frame, vals(CIL,'Frame'));
function draw(){{let ck=cell.value,fk=frame.value; const keep=r=>(!ck||norm(r.Cell_ID||r.Std_Cell_ID)===ck)&&(!fk||norm(r.Frame||r.Std_Frame)===fk); let s=SAT.filter(keep), c=CIL.filter(keep), traces=[];
traces.push({{type:'scatter',mode:'markers',name:'Satellites',x:s.map(r=>r.Std_X),y:s.map(r=>r.Std_Y),marker:{{size:6,opacity:.55}},text:s.map(r=>'Satellite '+r.Std_Satellite_ID),hovertemplate:'%{{text}}<br>x=%{{x}}<br>y=%{{y}}<extra></extra>'}});
traces.push({{type:'scatter',mode:'markers',name:'Cilia',x:c.map(r=>r.X),y:c.map(r=>r.Y),marker:{{size:12,symbol:'diamond'}},text:c.map(r=>'Cilium '+r.Cilium_ID+'<br>Length='+r.Length+' µm'),hovertemplate:'%{{text}}<br>x=%{{x}}<br>y=%{{y}}<extra></extra>'}});
// Draw exported skeleton centerlines when Fiji provided path points.
let sk=SK.filter(keep); let ids=[...new Set(sk.map(r=>norm(r.Cilium_ID)).filter(Boolean))]; ids.forEach(id=>{{let q=sk.filter(r=>norm(r.Cilium_ID)===id).sort((a,b)=>Number(a.Point_Index||0)-Number(b.Point_Index||0)); let x=q.map(r=>Number(r.X_um??r.X)), y=q.map(r=>Number(r.Y_um??r.Y)); if(x.length)traces.push({{type:'scatter',mode:'lines',name:'Cilium '+id+' skeleton',x:x,y:y,line:{{width:4}},showlegend:false,hoverinfo:'skip'}})}});
Plotly.newPlot('plot',traces,{{xaxis:{{title:'X (µm)',scaleanchor:'y',scaleratio:1}},yaxis:{{title:'Y (µm)',autorange:'reversed'}},legend:{{orientation:'h'}},margin:{{t:30}}}},{{responsive:true,toImageButtonOptions:{{format:'png',scale:3}}}})}}
cell.onchange=draw; frame.onchange=draw; draw();</script></body></html>"""
    path=vdir/satelliteq_dated_name(outdir,"Cilium_Viewer",ext=".html")
    path.write_text(html,encoding="utf-8")
    return path


def make_spatial_viewers(df: pd.DataFrame, outdir: Path, cilium_rows: Optional[pd.DataFrame] = None, skeleton_points: Optional[pd.DataFrame] = None) -> List[Path]:
    """Create one self-contained Plotly HTML spatial viewer per Source_Input.

    Real Fiji cell ROI boundaries are used when available. These come from
    *_cell_boundaries_for_python.csv files exported by Fiji. If exact boundary
    CSVs are absent, the viewer falls back to estimated circular/spherical
    boundaries. Optional primary-cilium coordinates and skeleton paths are
    integrated as a toggleable layer in the same viewer.
    """
    try:
        import json
        import plotly.offline as po
    except Exception:
        return []

    viewers_dir = outdir
    viewers_dir.mkdir(parents=True, exist_ok=True)
    plotlyjs = po.get_plotlyjs()

    outputs: List[Path] = []
    group_col = "Source_UID" if "Source_UID" in df.columns else ("Source_Input" if "Source_Input" in df.columns else None)
    grouped = df.groupby(group_col, dropna=False) if group_col else [("SatelliteQ", df)]

    def _numeric_series_local(table: pd.DataFrame, col: str) -> pd.Series:
        if col not in table.columns:
            return pd.Series(np.nan, index=table.index, dtype=float)
        return pd.to_numeric(table[col], errors="coerce")

    def _normalize_cilia_for_viewer(table: Optional[pd.DataFrame]) -> pd.DataFrame:
        if table is None or table.empty:
            return pd.DataFrame()
        t = table.copy()
        aliases = {
            "Std_Cell_ID": ["Std_Cell_ID", "Cell_ID"],
            "Std_Frame": ["Std_Frame", "Frame"],
            "Cilium_ID": ["Cilium_ID", "ID"],
            "Cilium_X": ["Cilium_X_um", "X_um", "X"],
            "Cilium_Y": ["Cilium_Y_um", "Y_um", "Y"],
            "Cilium_Z": ["Cilium_Z_um", "Z_um", "Z"],
            "Cilium_Length": ["Cilium_Length_um", "Length_um", "Length"],
        }
        out = pd.DataFrame(index=t.index)
        for dst, candidates in aliases.items():
            src = first_existing(list(t.columns), candidates)
            out[dst] = t[src] if src else np.nan
        for c in ["Std_Frame", "Cilium_X", "Cilium_Y", "Cilium_Z", "Cilium_Length"]:
            out[c] = pd.to_numeric(out[c], errors="coerce")
        return out[out["Cilium_X"].notna() & out["Cilium_Y"].notna()].copy()

    def _normalize_cilium_skeleton_for_viewer(table: Optional[pd.DataFrame]) -> pd.DataFrame:
        if table is None or table.empty:
            return pd.DataFrame()
        t = table.copy()
        aliases = {
            "Std_Cell_ID": ["Std_Cell_ID", "Cell_ID"],
            "Std_Frame": ["Std_Frame", "Frame"],
            "Cilium_ID": ["Cilium_ID", "ID"],
            "Point_Index": ["Point_Index", "Point", "Index"],
            "X": ["X_um", "X"],
            "Y": ["Y_um", "Y"],
            "Z": ["Z_um", "Z"],
        }
        out = pd.DataFrame(index=t.index)
        for dst, candidates in aliases.items():
            src = first_existing(list(t.columns), candidates)
            out[dst] = t[src] if src else np.nan
        for c in ["Std_Frame", "Point_Index", "X", "Y", "Z"]:
            out[c] = pd.to_numeric(out[c], errors="coerce")
        return out[out["X"].notna() & out["Y"].notna()].copy()

    def _viewer_companion_subset(table: Optional[pd.DataFrame], source_key, source_display="") -> pd.DataFrame:
        if table is None or table.empty:
            return pd.DataFrame()
        t = table.copy()
        if "Source_UID" in t.columns:
            mask = t["Source_UID"].astype(str).eq(str(source_key))
            if mask.any():
                return t.loc[mask].copy()
        if "Source_Input" in t.columns and source_display:
            mask = t["Source_Input"].astype(str).eq(str(source_display))
            if mask.any():
                return t.loc[mask].copy()
        return t.copy() if (group_col is None or len(pd.unique(df.get("Source_UID", df.get("Source_Input", pd.Series(["SatelliteQ"]))))) <= 1) else pd.DataFrame()

    def _candidate_boundary_files(gsrc: pd.DataFrame) -> List[Path]:
        files = []
        seen = set()
        if "Source_Path" in gsrc.columns:
            for sp in gsrc["Source_Path"].dropna().astype(str).unique():
                p = Path(sp)
                search_root = p.parent if p.suffix else p
                if search_root.exists():
                    for f in search_root.rglob("*cell_boundaries_for_python.csv"):
                        if f not in seen:
                            files.append(f)
                            seen.add(f)
        return files

    def _load_exact_boundaries(gsrc: pd.DataFrame) -> pd.DataFrame:
        files = _candidate_boundary_files(gsrc)
        if not files:
            return pd.DataFrame()

        frames = []
        for f in files:
            try:
                b = pd.read_csv(f)
                if b.empty:
                    continue
                b["Boundary_Source_File"] = f.name
                b["Boundary_Source_Path"] = str(f)
                frames.append(b)
            except Exception:
                continue
        if not frames:
            return pd.DataFrame()

        bdf = pd.concat(frames, ignore_index=True, sort=False)

        if "Image" in bdf.columns and "Image" in gsrc.columns:
            image_names = set(gsrc["Image"].dropna().astype(str).unique())
            if image_names:
                bdf = bdf[bdf["Image"].astype(str).isin(image_names)].copy()

        if bdf.empty:
            return pd.DataFrame()

        out = pd.DataFrame()
        out["Std_Cell_ID"] = bdf["Cell_ID"] if "Cell_ID" in bdf.columns else np.nan
        out["Std_Frame"] = pd.to_numeric(bdf["Frame"], errors="coerce") if "Frame" in bdf.columns else np.nan
        out["Point_Index"] = pd.to_numeric(bdf["Point_Index"], errors="coerce") if "Point_Index" in bdf.columns else np.arange(1, len(bdf) + 1)
        out["Boundary_X"] = pd.to_numeric(bdf["X_um"], errors="coerce") if "X_um" in bdf.columns else np.nan
        out["Boundary_Y"] = pd.to_numeric(bdf["Y_um"], errors="coerce") if "Y_um" in bdf.columns else np.nan
        out["Boundary_X_pixel"] = pd.to_numeric(bdf["X_pixel"], errors="coerce") if "X_pixel" in bdf.columns else np.nan
        out["Boundary_Y_pixel"] = pd.to_numeric(bdf["Y_pixel"], errors="coerce") if "Y_pixel" in bdf.columns else np.nan
        out["Image"] = bdf["Image"] if "Image" in bdf.columns else ""
        out["Run_ID"] = bdf["Run_ID"] if "Run_ID" in bdf.columns else ""
        out["Boundary_Source_File"] = bdf["Boundary_Source_File"]
        out = out[out["Boundary_X"].notna() & out["Boundary_Y"].notna()].copy()

        return out

    def _build_centrosomes(gsrc: pd.DataFrame, mode: str) -> pd.DataFrame:
        if gsrc.empty:
            return pd.DataFrame()
        t = gsrc.copy()
        cx = _numeric_series_local(t, "Std_Centrosome_X")
        cy = _numeric_series_local(t, "Std_Centrosome_Y")
        cz = _numeric_series_local(t, "Std_Centrosome_Z")
        sid = t.get("Std_Centrosome_ID", pd.Series([None] * len(t), index=t.index))
        cell = t.get("Std_Cell_ID", pd.Series([None] * len(t), index=t.index))
        frame = _numeric_series_local(t, "Std_Frame")
        tmin = _numeric_series_local(t, "Std_Time_min")
        sx = _numeric_series_local(t, "Std_X")
        sy = _numeric_series_local(t, "Std_Y")
        sz = _numeric_series_local(t, "Std_Z")

        tmp = pd.DataFrame({
            "Std_Centrosome_ID": sid,
            "Std_Cell_ID": cell,
            "Std_Frame": frame,
            "Std_Time_min": tmin,
            "Std_Centrosome_X": cx,
            "Std_Centrosome_Y": cy,
            "Std_Centrosome_Z": cz,
            "Std_X": sx,
            "Std_Y": sy,
            "Std_Z": sz,
        })
        tmp["FrameKey"] = tmp["Std_Frame"].fillna(-1)
        tmp["CentKey"] = tmp["Std_Centrosome_ID"].astype(str)
        tmp["CellKey"] = tmp["Std_Cell_ID"].astype(str)

        if tmp[["Std_Centrosome_X", "Std_Centrosome_Y"]].notna().sum().sum() == 0:
            grp = tmp.groupby(["FrameKey", "CellKey", "CentKey"], dropna=False)
            cent = grp.agg({
                "Std_X": "mean",
                "Std_Y": "mean",
                "Std_Z": "mean",
                "Std_Cell_ID": "first",
                "Std_Centrosome_ID": "first",
                "Std_Frame": "first",
                "Std_Time_min": "first",
            }).reset_index(drop=True)
            cent = cent.rename(columns={
                "Std_X": "Std_Centrosome_X",
                "Std_Y": "Std_Centrosome_Y",
                "Std_Z": "Std_Centrosome_Z",
            })
            cent["Coordinate_Source"] = "Derived from assigned satellite positions"
            return cent[[c for c in ["Std_Centrosome_ID", "Std_Cell_ID", "Std_Frame", "Std_Time_min", "Std_Centrosome_X", "Std_Centrosome_Y", "Std_Centrosome_Z", "Coordinate_Source"] if c in cent.columns]]

        explicit = tmp[tmp["Std_Centrosome_X"].notna() & tmp["Std_Centrosome_Y"].notna()].copy()
        if explicit.empty:
            return pd.DataFrame()
        grp = explicit.groupby(["FrameKey", "CellKey", "CentKey"], dropna=False)
        cent = grp.agg({
            "Std_Centrosome_X": "mean",
            "Std_Centrosome_Y": "mean",
            "Std_Centrosome_Z": "mean",
            "Std_Cell_ID": "first",
            "Std_Centrosome_ID": "first",
            "Std_Frame": "first",
            "Std_Time_min": "first",
        }).reset_index(drop=True)
        cent["Coordinate_Source"] = "From Fiji output"
        return cent[[c for c in ["Std_Centrosome_ID", "Std_Cell_ID", "Std_Frame", "Std_Time_min", "Std_Centrosome_X", "Std_Centrosome_Y", "Std_Centrosome_Z", "Coordinate_Source"] if c in cent.columns]]

    def _build_estimated_cells(gsrc: pd.DataFrame, cent_df: pd.DataFrame, mode: str) -> pd.DataFrame:
        if gsrc.empty or "Std_Cell_ID" not in gsrc.columns:
            return pd.DataFrame()
        t = gsrc.copy()
        t = t[t["Std_Cell_ID"].notna() & (t["Std_Cell_ID"].astype(str) != "")].copy()
        if t.empty:
            return pd.DataFrame()

        t["Std_Frame"] = _numeric_series_local(t, "Std_Frame")
        t["Std_X"] = _numeric_series_local(t, "Std_X")
        t["Std_Y"] = _numeric_series_local(t, "Std_Y")
        t["Std_Z"] = _numeric_series_local(t, "Std_Z")
        t["Std_Cell_Measure"] = _numeric_series_local(t, "Std_Cell_Measure")
        t["Std_Time_min"] = _numeric_series_local(t, "Std_Time_min")
        t["FrameKey"] = t["Std_Frame"].fillna(-1)
        t["CellKey"] = t["Std_Cell_ID"].astype(str)

        agg = {
            "Std_X": "mean", "Std_Y": "mean", "Std_Z": "mean",
            "Std_Cell_ID": "first", "Std_Frame": "first", "Std_Time_min": "first",
            "Std_Cell_Measure": "median",
        }
        if "Std_Cell_Measure_Name" in t.columns:
            agg["Std_Cell_Measure_Name"] = "first"

        sat_centers = t.groupby(["FrameKey", "CellKey"], dropna=False).agg(agg).reset_index(drop=True)

        if not cent_df.empty:
            cc = cent_df.copy()
            cc["FrameKey"] = pd.to_numeric(cc.get("Std_Frame", pd.Series(np.nan, index=cc.index)), errors="coerce").fillna(-1)
            cc["CellKey"] = cc.get("Std_Cell_ID", pd.Series([None]*len(cc), index=cc.index)).astype(str)
            cent_centers = cc.groupby(["FrameKey", "CellKey"], dropna=False).agg({
                "Std_Centrosome_X": "mean", "Std_Centrosome_Y": "mean", "Std_Centrosome_Z": "mean"
            }).reset_index()
            sat_centers["FrameKey"] = pd.to_numeric(sat_centers.get("Std_Frame", pd.Series(np.nan, index=sat_centers.index)), errors="coerce").fillna(-1)
            sat_centers["CellKey"] = sat_centers.get("Std_Cell_ID", pd.Series([None]*len(sat_centers), index=sat_centers.index)).astype(str)
            sat_centers = sat_centers.merge(cent_centers, on=["FrameKey", "CellKey"], how="left")
            sat_centers["Center_X"] = sat_centers["Std_Centrosome_X"].fillna(sat_centers["Std_X"])
            sat_centers["Center_Y"] = sat_centers["Std_Centrosome_Y"].fillna(sat_centers["Std_Y"])
            sat_centers["Center_Z"] = sat_centers["Std_Centrosome_Z"].fillna(sat_centers["Std_Z"])
            sat_centers["Center_Source"] = np.where(sat_centers["Std_Centrosome_X"].notna(), "Centrosome-based estimate", "Satellite-centroid estimate")
        else:
            sat_centers["Center_X"] = sat_centers["Std_X"]
            sat_centers["Center_Y"] = sat_centers["Std_Y"]
            sat_centers["Center_Z"] = sat_centers["Std_Z"]
            sat_centers["Center_Source"] = "Satellite-centroid estimate"

        if mode == "3D":
            sat_centers["Boundary_Radius"] = ((3.0 * sat_centers["Std_Cell_Measure"]) / (4.0 * np.pi)).pow(1.0/3.0)
        else:
            sat_centers["Boundary_Radius"] = np.sqrt(sat_centers["Std_Cell_Measure"] / np.pi)

        sat_centers = sat_centers.rename(columns={"Std_Cell_Measure": "Cell_Measure", "Std_Cell_Measure_Name": "Cell_Measure_Name"})
        sat_centers = sat_centers[sat_centers["Center_X"].notna() & sat_centers["Center_Y"].notna()].copy()
        return sat_centers[[c for c in ["Std_Cell_ID", "Std_Frame", "Std_Time_min", "Cell_Measure", "Cell_Measure_Name", "Center_X", "Center_Y", "Center_Z", "Boundary_Radius", "Center_Source"] if c in sat_centers.columns]]

    for src_key, g in grouped:
        g = g.copy()
        src_name = str(g["Source_Input"].iloc[0]) if "Source_Input" in g.columns and not g.empty else str(src_key)
        if g.empty:
            continue
        if "Std_X" not in g.columns or "Std_Y" not in g.columns:
            continue

        x = pd.to_numeric(g["Std_X"], errors="coerce")
        y = pd.to_numeric(g["Std_Y"], errors="coerce")
        valid = x.notna() & y.notna()
        if valid.sum() == 0:
            continue
        g = g.loc[valid].copy()

        mode = _viewer_mode_for_subset(g)
        has_time = _viewer_time_for_subset(g)
        has_z = mode == "3D" and pd.to_numeric(g.get("Std_Z", pd.Series(dtype=float)), errors="coerce").notna().any()

        cent_df = _build_centrosomes(g, mode)
        exact_boundary_df = _load_exact_boundaries(g)
        estimated_cell_df = _build_estimated_cells(g, cent_df, mode)
        cilia_df = _normalize_cilia_for_viewer(_viewer_companion_subset(cilium_rows, src_key, src_name))
        cilium_skeleton_df = _normalize_cilium_skeleton_for_viewer(_viewer_companion_subset(skeleton_points, src_key, src_name))

        keep_cols = [
            "Source_Input", "Source_Path", "Source_Folder", "Mode", "Analysis_Type",
            "Std_Satellite_ID", "Std_Cell_ID", "Std_Centrosome_ID",
            "Std_Frame", "Std_Time_min",
            "Std_X", "Std_Y", "Std_Z",
            "Std_Centrosome_X", "Std_Centrosome_Y", "Std_Centrosome_Z",
            "Std_Area_2D", "Std_Volume_3D", "Std_Size", "Std_Surface", "Std_Sphericity", "Std_Sphericity_Source", "Std_Circularity_2D", "Std_Shape_Compactness",
            "Std_Intensity", "Std_Mean_Intensity", "Std_Distance", "Std_Distance_to_Cell_or_Image_Center_XY",
            "Std_Cell_Spatial_Entropy", "Std_Cell_Pericentrosomal_Clustering_Index",
            PCI_STD_INNER_COL, PCI_STD_OUTER_COL, PCI_STD_TOTAL_COL,
            "Std_Cell_Area_2D", "Std_Cell_Volume_3D", "Std_Cell_Measure", "Std_Cell_Measure_Name"
        ]
        cols = [c for c in keep_cols if c in g.columns]
        g = g[cols].copy()

        g = g.replace({np.nan: None})
        cent_df = cent_df.replace({np.nan: None}) if not cent_df.empty else pd.DataFrame()
        estimated_cell_df = estimated_cell_df.replace({np.nan: None}) if not estimated_cell_df.empty else pd.DataFrame()
        exact_boundary_df = exact_boundary_df.replace({np.nan: None}) if not exact_boundary_df.empty else pd.DataFrame()
        cilia_df = cilia_df.replace({np.nan: None}) if not cilia_df.empty else pd.DataFrame()
        cilium_skeleton_df = cilium_skeleton_df.replace({np.nan: None}) if not cilium_skeleton_df.empty else pd.DataFrame()

        rows = g.to_dict(orient="records")
        cent_rows = cent_df.to_dict(orient="records") if not cent_df.empty else []
        estimated_cell_rows = estimated_cell_df.to_dict(orient="records") if not estimated_cell_df.empty else []
        exact_boundary_rows = exact_boundary_df.to_dict(orient="records") if not exact_boundary_df.empty else []
        cilia_rows_view = cilia_df.to_dict(orient="records") if not cilia_df.empty else []
        cilium_skeleton_rows = cilium_skeleton_df.to_dict(orient="records") if not cilium_skeleton_df.empty else []

        color_candidates = []
        for c in ["Std_Intensity", "Std_Area_2D", "Std_Volume_3D", "Std_Size", "Std_Distance", "Std_Distance_to_Cell_or_Image_Center_XY", "Std_Cell_Spatial_Entropy", "Std_Cell_Pericentrosomal_Clustering_Index", PCI_STD_INNER_COL, PCI_STD_OUTER_COL, PCI_STD_TOTAL_COL, "Std_Shape_Compactness", "Std_Cell_ID", "Std_Centrosome_ID", "Std_Frame", "Source_Folder"]:
            if c in g.columns and pd.Series([r.get(c) for r in rows]).replace("", np.nan).notna().any():
                color_candidates.append(c)
        if not color_candidates:
            color_candidates = ["Std_Size"] if "Std_Size" in g.columns else [cols[0]]

        size_candidates = ["Fixed"]
        for c in ["Std_Area_2D", "Std_Volume_3D", "Std_Size", "Std_Intensity", "Std_Distance", "Std_Shape_Compactness"]:
            if c in g.columns and pd.to_numeric(g[c], errors="coerce").notna().any():
                size_candidates.append(c)

        labels = {
            "Std_Satellite_ID": "Satellite ID",
            "Std_Cell_ID": "Cell ID",
            "Std_Centrosome_ID": "Centrosome ID",
            "Std_Frame": "Frame",
            "Std_Time_min": "Time (min)",
            "Std_X": "X position (µm)",
            "Std_Y": "Y position (µm)",
            "Std_Z": "Z position (µm)",
            "Std_Area_2D": "Satellite area (2D; µm²)",
            "Std_Volume_3D": "Satellite volume (3D; µm³)",
            "Std_Size": "Satellite area (µm²)" if mode == "2D" else "Satellite volume (µm³)",
            "Std_Surface": "Satellite surface area (µm²)",
            "Std_Intensity": "Satellite integrated intensity (a.u.)",
            "Std_Mean_Intensity": "Satellite mean intensity (a.u.)",
            "Std_Distance": "Satellite-to-centrosome distance (µm)",
            "Std_Cell_Spatial_Entropy": "Cell spatial entropy (0–1)",
            "Std_Cell_Pericentrosomal_Clustering_Index": "Pericentrosomal clustering index (0–3 µm / 3–12 µm)",
            PCI_STD_INNER_COL: "Pericentrosomal satellite fluorescence, 0–3 µm (a.u.)",
            PCI_STD_OUTER_COL: "Outer satellite fluorescence, 3–12 µm (a.u.)",
            PCI_STD_TOTAL_COL: "Total satellite fluorescence, 0–12 µm (a.u.)",
            "Std_Distance_to_Cell_or_Image_Center_XY": "Distance to cell/image center (XY)",
            "Std_Cell_Area_2D": "Cell area (2D; µm²)",
            "Std_Cell_Volume_3D": "Cell volume (3D; µm³)",
            "Std_Cell_Measure": "Cell area (µm²)" if mode == "2D" else "Cell volume (µm³)",
            "Fixed": "Fixed marker size",
            "Source_Folder": "Input folder",
        }

        frames_available = sorted(pd.Series([r.get("Std_Frame") for r in rows]).dropna().astype(float).unique().tolist()) if has_time else []

        data_json = json.dumps({
            "rows": rows,
            "centrosomes": cent_rows,
            "estimatedCells": estimated_cell_rows,
            "exactBoundaries": exact_boundary_rows,
            "cilia": cilia_rows_view,
            "ciliumSkeleton": cilium_skeleton_rows,
            "mode": mode,
            "hasZ": bool(has_z),
            "hasTime": bool(has_time),
            "frames": frames_available,
            "colorCandidates": color_candidates,
            "sizeCandidates": size_candidates,
            "labels": labels,
            "sourceName": str(src_name),
        }, ensure_ascii=False)

        html_template = """<!DOCTYPE html>
<html>
<head>
<meta charset=\"utf-8\">
<title>SatelliteQ Spatial Viewer</title>
<style>
body { font-family: Arial, sans-serif; margin: 22px; background: #fbf8ff; color: #222; }
h1 { color: #7E3FB5; margin-bottom: 4px; }
.small { color:#555; font-size:13px; margin-bottom:14px; max-width:1300px; line-height:1.4; }
.controls { display:grid; grid-template-columns: repeat(4, minmax(180px, 1fr)); gap:12px; background:white; border:1px solid #DCC6FF; padding:14px; border-radius:10px; box-shadow:0 2px 10px rgba(126,63,181,0.08); }
.controlbox { background:#fff; }
label { font-weight:bold; font-size:13px; display:block; margin-bottom:4px; color:#4B216E; }
select { width:100%; padding:6px; }
.chkrow { display:flex; gap:12px; align-items:center; flex-wrap:wrap; padding-top:4px; }
.chkrow label { font-weight:normal; color:#333; margin:0; }
#plot { background:white; border:1px solid #DCC6FF; border-radius:10px; margin-top:16px; padding:8px; box-shadow:0 2px 10px rgba(126,63,181,0.08); }
.note { margin-top:10px; color:#555; font-size:12px; }
</style>
<script type=\"text/javascript\">__PLOTLYJS__</script>
</head>
<body>
<h1>SatelliteQ Spatial Viewer</h1>
<div class="small">SatelliteQ version __SATELLITEQ_VERSION__</div>
<div class=\"small\">Interactive spatial viewer generated from SatelliteQ Python output. Satellites, centrosomes, primary cilia, and cell boundaries can be toggled independently. If Fiji exported exact 2D cell ROI coordinates, those real boundaries are shown in red. Otherwise, the viewer falls back to estimated boundaries from cell area/volume.</div>
<div class=\"controls\">
  <div class=\"controlbox\"><label>Color satellites by</label><select id=\"colorSelect\"></select></div>
  <div class=\"controlbox\"><label>Satellite marker size</label><select id=\"sizeSelect\"></select></div>
  <div class=\"controlbox\"><label>Filter by cell</label><select id=\"cellSelect\"></select></div>
  <div class=\"controlbox\" id=\"frameBox\"><label>Frame</label><select id=\"frameSelect\"></select></div>
  <div class=\"controlbox\" style=\"grid-column:1 / span 4;\">
    <label>Visible layers</label>
    <div class=\"chkrow\">
      <label><input type=\"checkbox\" id=\"showSat\" checked> Satellites</label>
      <label><input type=\"checkbox\" id=\"showCent\" checked> Centrosomes</label>
      <label><input type=\"checkbox\" id=\"showCilia\" checked> Cilia</label>
      <label><input type=\"checkbox\" id=\"showCells\" checked> Cell boundaries</label>
    </div>
  </div>
</div>
<div id=\"plot\" style=\"height:800px;\"></div>
<div class=\"note\" id=\"boundaryNote\"></div>
<script>
const DATA = __DATA_JSON__;
const ROWS = DATA.rows;
const CENTS = DATA.centrosomes || [];
const EST_CELLS = DATA.estimatedCells || [];
const BOUNDARIES = DATA.exactBoundaries || [];
const CILIA = DATA.cilia || [];
const CIL_SK = DATA.ciliumSkeleton || [];
const MODE = DATA.mode;
const HAS_Z = DATA.hasZ;
const HAS_TIME = DATA.hasTime;
const FRAMES = DATA.frames || [];
const colorSelect = document.getElementById('colorSelect');
const sizeSelect = document.getElementById('sizeSelect');
const cellSelect = document.getElementById('cellSelect');
const frameSelect = document.getElementById('frameSelect');
const showSat = document.getElementById('showSat');
const showCent = document.getElementById('showCent');
const showCilia = document.getElementById('showCilia');
const showCells = document.getElementById('showCells');
const frameBox = document.getElementById('frameBox');
const boundaryNote = document.getElementById('boundaryNote');
const PURPLE = '#7E3FB5';
const DARK = '#4B216E';
const CYAN = '#16D6FF';
const RED = '#E24A63';
const CILIUM_COLOR = '#F59E0B';
const PALETTE = ['#7E3FB5','#A66BDA','#C39BEA','#6B4FA3','#B57EDC','#8E5CC2','#D6B7FF','#5D3A91','#CB6CE6','#9F86FF','#43AA8B','#4D96FF'];

function lbl(c) { return DATA.labels[c] || c; }
function uniq(arr) { return Array.from(new Set(arr.filter(v => v !== null && v !== undefined && v !== ''))); }
function isNumericValue(v) { return v !== null && v !== '' && !isNaN(Number(v)); }
function num(v) { return isNumericValue(v) ? Number(v) : null; }
function addOptions(sel, values, noneLabel=null) {
  sel.innerHTML = '';
  if (noneLabel !== null) {
    const o = document.createElement('option'); o.value = ''; o.textContent = noneLabel; sel.appendChild(o);
  }
  values.forEach(v => { const o = document.createElement('option'); o.value = v; o.textContent = lbl(v) || String(v); sel.appendChild(o); });
}
function addFilterOptions(sel, values, labelPrefix) {
  sel.innerHTML = '';
  const allOpt = document.createElement('option'); allOpt.value = ''; allOpt.textContent = 'All'; sel.appendChild(allOpt);
  values.forEach(v => { const o = document.createElement('option'); o.value = String(v); o.textContent = labelPrefix + String(v); sel.appendChild(o); });
}
function closeEnough(a,b) {
  if (a === null || a === undefined || b === null || b === undefined || a === '' || b === '') return false;
  return Math.abs(Number(a)-Number(b)) < 1e-9;
}
function initControls() {
  addOptions(colorSelect, DATA.colorCandidates, null);
  addOptions(sizeSelect, DATA.sizeCandidates, null);
  const cells = uniq(ROWS.map(r => r.Std_Cell_ID).concat(CILIA.map(r => r.Std_Cell_ID)));
  addFilterOptions(cellSelect, cells, 'Cell ');
  if (showCilia) {
    const hasCilia = CILIA.length > 0 || CIL_SK.length > 0;
    showCilia.disabled = !hasCilia;
    showCilia.checked = hasCilia;
  }
  if (HAS_TIME && FRAMES.length > 0) {
    frameBox.style.display = 'block';
    frameSelect.innerHTML = '';
    FRAMES.forEach(v => { const o = document.createElement('option'); o.value = String(v); o.textContent = 'Frame ' + String(v); frameSelect.appendChild(o); });
    frameSelect.value = String(FRAMES[0]);
  } else {
    frameBox.style.display = 'none';
  }
  if (DATA.colorCandidates.includes('Std_Size')) colorSelect.value = 'Std_Size';
  if (DATA.sizeCandidates.includes('Std_Size')) sizeSelect.value = 'Std_Size'; else sizeSelect.value = 'Fixed';
  boundaryNote.innerHTML = BOUNDARIES.length > 0
    ? 'Cell boundaries: using <b>exact Fiji ROI coordinates</b> exported from the Fiji macro.'
    : 'Cell boundaries: exact Fiji ROI coordinates were not found, so boundaries are estimated from cell area/volume.';
}
function activeFrameValue() {
  if (!HAS_TIME || !frameSelect.value) return null;
  return Number(frameSelect.value);
}
function inFrame(row) {
  const af = activeFrameValue();
  if (af === null) return true;
  if (row.Std_Frame === -1 || row.Std_Frame === '-1') return false;
  return closeEnough(row.Std_Frame, af);
}
function inCell(row) {
  if (!cellSelect.value) return true;
  return String(row.Std_Cell_ID) === String(cellSelect.value);
}
function filteredRows() { return ROWS.filter(r => inFrame(r) && inCell(r)); }
function filteredCentrosomes() { return CENTS.filter(r => inFrame(r) && inCell(r)); }
function filteredEstimatedCells() { return EST_CELLS.filter(r => inFrame(r) && inCell(r)); }
function filteredBoundaries() { return BOUNDARIES.filter(r => inFrame(r) && inCell(r)); }
function filteredCilia() { return CILIA.filter(r => inFrame(r) && inCell(r)); }
function filteredCiliumSkeleton() { return CIL_SK.filter(r => inFrame(r) && inCell(r)); }
function buildSatelliteHover(r) {
  let out = [];
  if (r.Std_Satellite_ID !== null && r.Std_Satellite_ID !== undefined) out.push('Satellite: ' + r.Std_Satellite_ID);
  if (r.Std_Cell_ID !== null && r.Std_Cell_ID !== undefined && r.Std_Cell_ID !== '') out.push('Cell: ' + r.Std_Cell_ID);
  if (r.Std_Centrosome_ID !== null && r.Std_Centrosome_ID !== undefined && r.Std_Centrosome_ID !== '') out.push('Centrosome: ' + r.Std_Centrosome_ID);
  if (r.Std_Frame !== null && r.Std_Frame !== undefined) out.push('Frame: ' + r.Std_Frame);
  if (r.Std_Time_min !== null && r.Std_Time_min !== undefined) out.push('Time (min): ' + r.Std_Time_min);
  if (r.Std_Size !== null && r.Std_Size !== undefined) out.push(lbl('Std_Size') + ': ' + r.Std_Size);
  if (r.Std_Intensity !== null && r.Std_Intensity !== undefined) out.push(lbl('Std_Intensity') + ': ' + r.Std_Intensity);
  if (r.Std_Distance !== null && r.Std_Distance !== undefined) out.push(lbl('Std_Distance') + ': ' + r.Std_Distance);
  return out.join('<br>');
}
function buildCentrosomeHover(r) {
  let out = ['Centrosome: ' + (r.Std_Centrosome_ID ?? 'NA')];
  if (r.Std_Cell_ID !== null && r.Std_Cell_ID !== undefined && r.Std_Cell_ID !== '') out.push('Cell: ' + r.Std_Cell_ID);
  if (r.Std_Frame !== null && r.Std_Frame !== undefined) out.push('Frame: ' + r.Std_Frame);
  if (r.Std_Time_min !== null && r.Std_Time_min !== undefined) out.push('Time (min): ' + r.Std_Time_min);
  if (r.Std_Centrosome_X !== null && r.Std_Centrosome_X !== undefined) out.push('X (µm): ' + r.Std_Centrosome_X);
  if (r.Std_Centrosome_Y !== null && r.Std_Centrosome_Y !== undefined) out.push('Y (µm): ' + r.Std_Centrosome_Y);
  if (HAS_Z && r.Std_Centrosome_Z !== null && r.Std_Centrosome_Z !== undefined) out.push('Z (µm): ' + r.Std_Centrosome_Z);
  if (r.Coordinate_Source) out.push('Coordinate source: ' + r.Coordinate_Source);
  return out.join('<br>');
}
function buildCiliumHover(r) {
  let out = ['Cilium: ' + (r.Cilium_ID ?? 'NA')];
  if (r.Std_Cell_ID !== null && r.Std_Cell_ID !== undefined && r.Std_Cell_ID !== '') out.push('Cell: ' + r.Std_Cell_ID);
  if (r.Std_Frame !== null && r.Std_Frame !== undefined) out.push('Frame: ' + r.Std_Frame);
  if (r.Cilium_Length !== null && r.Cilium_Length !== undefined) out.push('Length (µm): ' + r.Cilium_Length);
  if (r.Cilium_X !== null && r.Cilium_X !== undefined) out.push('X (µm): ' + r.Cilium_X);
  if (r.Cilium_Y !== null && r.Cilium_Y !== undefined) out.push('Y (µm): ' + r.Cilium_Y);
  if (HAS_Z && r.Cilium_Z !== null && r.Cilium_Z !== undefined) out.push('Z (µm): ' + r.Cilium_Z);
  return out.join('<br>');
}
function buildCellHover(r) {
  let out = ['Cell: ' + (r.Std_Cell_ID ?? 'NA')];
  if (r.Std_Frame !== null && r.Std_Frame !== undefined) out.push('Frame: ' + r.Std_Frame);
  if (r.Std_Time_min !== null && r.Std_Time_min !== undefined) out.push('Time (min): ' + r.Std_Time_min);
  if (r.Cell_Measure !== null && r.Cell_Measure !== undefined) out.push((r.Cell_Measure_Name || (MODE==='3D' ? 'Cell volume (µm³)' : 'Cell area (µm²)')) + ': ' + r.Cell_Measure);
  if (r.Boundary_Radius !== null && r.Boundary_Radius !== undefined) out.push('Estimated boundary radius (µm): ' + r.Boundary_Radius);
  if (r.Center_Source) out.push('Center source: ' + r.Center_Source);
  return out.join('<br>');
}
function markerArrays(rr, colorField, sizeField) {
  let marker = {opacity: 0.84, line:{color: DARK, width:0.4}, sizemode:'diameter'};

  function isDistanceField(field) {
    return String(field || '').toLowerCase().indexOf('distance') >= 0;
  }

  function sizeTransform(field, rawValue) {
    const v = num(rawValue);
    if (v === null || !isFinite(v)) return null;
    const f = String(field || '').toLowerCase();
    if (f.indexOf('volume') >= 0) return Math.cbrt(Math.max(v, 0));
    if (f.indexOf('size') >= 0 || f.indexOf('area') >= 0 || f.indexOf('surface') >= 0) return Math.sqrt(Math.max(v, 0));
    if (f.indexOf('intensity') >= 0 || f.indexOf('intden') >= 0) return Math.sqrt(Math.max(v, 0));
    return v;
  }

  const colorVals = rr.map(r => r[colorField]);
  const numericColor = colorVals.filter(v => isNumericValue(v)).length >= Math.max(3, Math.floor(rr.length * 0.5));
  if (numericColor) {
    marker.color = colorVals.map(v => num(v));
    marker.colorscale = 'Viridis';
    marker.colorbar = {title: {text: lbl(colorField), side: 'right'}};
    marker.showscale = true;
    const finite = marker.color.filter(v => v !== null && isFinite(v));
    if (finite.length > 0) {
      marker.cmin = Math.min(...finite);
      marker.cmax = Math.max(...finite);
    }
  } else {
    const cats = uniq(colorVals.map(v => String(v)));
    const cmap = {};
    cats.forEach((c, i) => cmap[c] = PALETTE[i % PALETTE.length]);
    marker.color = colorVals.map(v => cmap[String(v)] || PURPLE);
    marker.showscale = false;
  }

  const fixed2D = 6.0, fixed3D = 4.4;
  const min2D = 4.0, max2D = 13.5;
  const min3D = 3.2, max3D = 10.0;

  if (sizeField === 'Fixed') {
    marker.size = rr.map(_ => MODE === '3D' ? fixed3D : fixed2D);
  } else {
    const vals = rr
      .map(r => sizeTransform(sizeField, r[sizeField]))
      .filter(v => v !== null && isFinite(v));

    if (vals.length === 0) {
      marker.size = rr.map(_ => MODE === '3D' ? fixed3D : fixed2D);
    } else {
      const vmin = Math.min(...vals), vmax = Math.max(...vals);
      marker.size = rr.map(r => {
        const v = sizeTransform(sizeField, r[sizeField]);
        if (v === null || !isFinite(v) || vmax === vmin) return MODE === '3D' ? fixed3D : fixed2D;
        let norm = Math.max(0, Math.min(1, (v - vmin) / (vmax - vmin)));
        if (isDistanceField(sizeField)) norm = 1 - norm;
        return MODE === '3D'
          ? (min3D + norm * (max3D - min3D))
          : (min2D + norm * (max2D - min2D));
      });
    }
  }

  return marker;
}
function satelliteTrace(rr) {
  const colorField = colorSelect.value;
  const sizeField = sizeSelect.value;
  const marker = markerArrays(rr, colorField, sizeField);
  const hover = rr.map(buildSatelliteHover);
  if (MODE === '3D' && HAS_Z) {
    return {type:'scatter3d', mode:'markers', name:'Satellites', x:rr.map(r=>num(r.Std_X)), y:rr.map(r=>num(r.Std_Y)), z:rr.map(r=>num(r.Std_Z)), marker:marker, text:hover, hovertemplate:'%{text}<extra></extra>'};
  }
  return {type:'scattergl', mode:'markers', name:'Satellites', x:rr.map(r=>num(r.Std_X)), y:rr.map(r=>num(r.Std_Y)), marker:marker, text:hover, hovertemplate:'%{text}<extra></extra>'};
}
function centrosomeTrace(rr) {
  if (MODE === '3D' && HAS_Z) {
    return {type:'scatter3d', mode:'markers', name:'Centrosomes', x:rr.map(r=>num(r.Std_Centrosome_X)), y:rr.map(r=>num(r.Std_Centrosome_Y)), z:rr.map(r=>num(r.Std_Centrosome_Z)), marker:{size:6, color:CYAN, opacity:1, line:{color:'#0C4A6E', width:1}}, text:rr.map(buildCentrosomeHover), hovertemplate:'%{text}<extra></extra>'};
  }
  return {type:'scattergl', mode:'markers', name:'Centrosomes', x:rr.map(r=>num(r.Std_Centrosome_X)), y:rr.map(r=>num(r.Std_Centrosome_Y)), marker:{size:15, color:CYAN, symbol:'x', opacity:1, line:{color:'#0C4A6E', width:1}}, text:rr.map(buildCentrosomeHover), hovertemplate:'%{text}<extra></extra>'};
}
function ciliumTraces() {
  const c = filteredCilia();
  const sk = filteredCiliumSkeleton();
  let traces = [];
  if (c.length > 0) {
    if (MODE === '3D' && HAS_Z) {
      traces.push({type:'scatter3d', mode:'markers', name:'Cilia', x:c.map(r=>num(r.Cilium_X)), y:c.map(r=>num(r.Cilium_Y)), z:c.map(r=>num(r.Cilium_Z)), marker:{size:7, color:CILIUM_COLOR, symbol:'diamond', opacity:1, line:{color:'#92400E', width:1}}, text:c.map(buildCiliumHover), hovertemplate:'%{text}<extra></extra>', legendgroup:'cilia'});
    } else {
      traces.push({type:'scattergl', mode:'markers', name:'Cilia', x:c.map(r=>num(r.Cilium_X)), y:c.map(r=>num(r.Cilium_Y)), marker:{size:13, color:CILIUM_COLOR, symbol:'diamond', opacity:1, line:{color:'#92400E', width:1}}, text:c.map(buildCiliumHover), hovertemplate:'%{text}<extra></extra>', legendgroup:'cilia'});
    }
  }
  let groups = {};
  sk.forEach(r => {
    const key = String(r.Std_Frame ?? '') + '|' + String(r.Std_Cell_ID ?? '') + '|' + String(r.Cilium_ID ?? '');
    if (!groups[key]) groups[key] = [];
    groups[key].push(r);
  });
  Object.keys(groups).forEach(k => {
    const pts = groups[k].sort((a,b) => Number(a.Point_Index || 0) - Number(b.Point_Index || 0));
    const cid = pts[0]?.Cilium_ID ?? 'NA';
    if (MODE === '3D' && HAS_Z) {
      traces.push({type:'scatter3d', mode:'lines', name:'Cilium ' + cid + ' path', x:pts.map(r=>num(r.X)), y:pts.map(r=>num(r.Y)), z:pts.map(r=>num(r.Z)), line:{color:CILIUM_COLOR, width:6}, hoverinfo:'skip', showlegend:false, legendgroup:'cilia'});
    } else {
      traces.push({type:'scattergl', mode:'lines', name:'Cilium ' + cid + ' path', x:pts.map(r=>num(r.X)), y:pts.map(r=>num(r.Y)), line:{color:CILIUM_COLOR, width:4}, hoverinfo:'skip', showlegend:false, legendgroup:'cilia'});
    }
  });
  return traces;
}
function exactBoundaryTraces(rr) {
  let groups = {};
  rr.forEach(r => {
    const key = String(r.Std_Frame ?? '') + '|' + String(r.Std_Cell_ID ?? '') + '|' + String(r.Image ?? '') + '|' + String(r.Boundary_Source_File ?? '');
    if (!groups[key]) groups[key] = [];
    groups[key].push(r);
  });
  let traces = [];
  Object.keys(groups).forEach(k => {
    const pts = groups[k].sort((a,b) => Number(a.Point_Index || 0) - Number(b.Point_Index || 0));
    const cell = pts[0]?.Std_Cell_ID ?? 'NA';
    const frame = pts[0]?.Std_Frame ?? '';
    const hover = 'Cell: ' + cell + '<br>Boundary: exact Fiji ROI coordinates' + (frame !== '' && frame !== -1 ? '<br>Frame: ' + frame : '');
    traces.push({
      type:'scattergl', mode:'lines', name:'Cell ' + cell + ' ROI boundary',
      x:pts.map(r => num(r.Boundary_X)), y:pts.map(r => num(r.Boundary_Y)),
      line:{color:RED, width:2.5}, hovertext:hover, hoverinfo:'text', legendgroup:'cell_boundaries'
    });
  });
  return traces;
}
function estimatedCellTraces(rr) {
  let traces = [];
  rr.forEach((r, idx) => {
    const cx = num(r.Center_X), cy = num(r.Center_Y), cz = num(r.Center_Z), rad = num(r.Boundary_Radius);
    if (cx === null || cy === null || rad === null) return;
    const hover = buildCellHover(r);
    if (MODE === '3D' && HAS_Z) {
      const u = []; const v = [];
      for (let i=0;i<=12;i++) u.push(2*Math.PI*i/12.0);
      for (let j=0;j<=6;j++) v.push(Math.PI*j/6.0);
      const xs=[], ys=[], zs=[];
      u.forEach(uu => v.forEach(vv => { xs.push(cx + rad*Math.cos(uu)*Math.sin(vv)); ys.push(cy + rad*Math.sin(uu)*Math.sin(vv)); zs.push((cz===null?0:cz) + rad*Math.cos(vv)); }));
      traces.push({type:'mesh3d', name:'Cell ' + r.Std_Cell_ID + ' estimated boundary', x:xs, y:ys, z:zs, opacity:0.12, color:RED, hovertext:hover, hoverinfo:'text', showscale:false, legendgroup:'cells'});
    } else {
      const pts = 72; const xs=[], ys=[];
      for (let i=0;i<=pts;i++) { const ang = 2*Math.PI*i/pts; xs.push(cx + rad*Math.cos(ang)); ys.push(cy + rad*Math.sin(ang)); }
      traces.push({type:'scattergl', mode:'lines', name:'Cell ' + r.Std_Cell_ID + ' estimated boundary', x:xs, y:ys, line:{color:RED, width:2}, hovertext:hover, hoverinfo:'text', legendgroup:'cells'});
    }
  });
  return traces;
}
function cellBoundaryTraces() {
  const exact = filteredBoundaries();
  if (exact.length > 0 && MODE === '2D') return exactBoundaryTraces(exact);
  return estimatedCellTraces(filteredEstimatedCells());
}
function render() {
  const sats = filteredRows();
  const cents = filteredCentrosomes();
  let traces = [];
  if (showSat.checked && sats.length > 0) traces.push(satelliteTrace(sats));
  if (showCent.checked && cents.length > 0) traces.push(centrosomeTrace(cents));
  if (showCilia && showCilia.checked) traces = traces.concat(ciliumTraces());
  if (showCells.checked) traces = traces.concat(cellBoundaryTraces());
  if (traces.length === 0) traces = MODE === '3D' && HAS_Z ? [{type:'scatter3d', x:[], y:[], z:[]}] : [{type:'scattergl', x:[], y:[]}];
  let title = DATA.sourceName + ' — Spatial Viewer';
  if (HAS_TIME && frameSelect.value) title += ' (Frame ' + frameSelect.value + ')';
  const layout = MODE === '3D' && HAS_Z ? {
    title:{text:title, font:{color:PURPLE, size:20}},
    scene:{xaxis:{title:'X position (µm)'}, yaxis:{title:'Y position (µm)'}, zaxis:{title:'Z position (µm)'}, aspectmode:'data'},
    paper_bgcolor:'white', plot_bgcolor:'white', margin:{l:0,r:0,t:55,b:0},
    legend:{orientation:'h', y:1.02, x:0}
  } : {
    title:{text:title, font:{color:PURPLE, size:20}},
    xaxis:{title:'X position (µm)', zeroline:false, linecolor:'#444444'},
    yaxis:{title:'Y position (µm)', zeroline:false, linecolor:'#444444', autorange:'reversed', scaleanchor:'x', scaleratio:1},
    paper_bgcolor:'white', plot_bgcolor:'white', margin:{l:70,r:30,t:55,b:70},
    legend:{orientation:'h', y:1.02, x:0}
  };
  Plotly.newPlot('plot', traces, layout, {responsive:true, toImageButtonOptions:{format:'png', filename:'SatelliteQ_spatial_viewer', scale:3}});
}
[colorSelect, sizeSelect, cellSelect, frameSelect, showSat, showCent, showCilia, showCells].forEach(el => { if (el) el.addEventListener('change', render); });
initControls();
render();
</script>
</body>
</html>
"""

        html = html_template.replace("__PLOTLYJS__", plotlyjs).replace("__DATA_JSON__", data_json).replace("__SATELLITEQ_VERSION__", SATELLITEQ_VERSION)
        viewer_stem = safe_name(Path(str(src_name)).stem)
        if group_col == "Source_UID":
            uid_suffix = safe_name(str(src_key))[-10:]
            viewer_stem = viewer_stem + "_" + uid_suffix
        outpath = viewers_dir / f"{viewer_stem}_Spatial_Viewer.html"
        outpath.write_text(html, encoding="utf-8")
        outputs.append(outpath)

    if outputs:
        def _make_spatial_tabs(viewer_paths: List[Path]) -> Path:
            tab_buttons = []
            frames = []
            for i, vp in enumerate(viewer_paths):
                label = vp.stem.replace("_Spatial_Viewer", "").replace("_", " ")
                try:
                    srcdoc = html_lib.escape(vp.read_text(encoding="utf-8"), quote=True)
                except Exception:
                    srcdoc = html_lib.escape("<html><body><h2>Viewer not available</h2></body></html>", quote=True)
                active = " active" if i == 0 else ""
                tab_buttons.append(f'<button class="tabbtn{active}" onclick="showSpatialTab({i})">{html_lib.escape(label)}</button>')
                frames.append(f'<iframe id="spatialFrame{i}" class="spatialFrame{active}" srcdoc="{srcdoc}"></iframe>')
            html_doc = f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>SatelliteQ Spatial Viewers</title>
<style>
body {{ margin:0; font-family:Arial, sans-serif; background:#f8fafc; color:#111827; }}
.header {{ padding:12px 16px; background:white; border-bottom:1px solid #e5e7eb; }}
h1 {{ margin:0 0 8px 0; font-size:20px; }}
.tabs {{ display:flex; gap:6px; flex-wrap:wrap; }}
.tabbtn {{ border:1px solid #d1d5db; background:#f3f4f6; padding:7px 10px; border-radius:8px; cursor:pointer; font-size:13px; }}
.tabbtn.active {{ background:#4b1d7a; color:white; border-color:#4b1d7a; }}
.spatialFrame {{ width:100%; height:calc(100vh - 86px); border:0; display:none; background:white; }}
.spatialFrame.active {{ display:block; }}
</style>
</head>
<body>
<div class="header">
  <h1>SatelliteQ Spatial Viewers</h1>
  <div class="tabs">{''.join(tab_buttons)}</div>
</div>
{''.join(frames)}
<script>
function showSpatialTab(idx) {{
  document.querySelectorAll('.tabbtn').forEach((b,i)=>b.classList.toggle('active', i===idx));
  document.querySelectorAll('.spatialFrame').forEach((f,i)=>f.classList.toggle('active', i===idx));
}}
</script>
</body>
</html>
"""
            out_tabs = viewers_dir / satelliteq_dated_name(outdir, "Spatial_Viewer_Tabs", ext=".html")
            out_tabs.write_text(html_doc, encoding="utf-8")
            for vp in viewer_paths:
                try:
                    vp.unlink()
                except Exception:
                    pass
            return out_tabs

        tabbed_path = _make_spatial_tabs(outputs)
        outputs = [tabbed_path]
    return outputs














def satelliteq_result_column_explanations(columns: List[str], analysis_kind: str = "quantification") -> pd.DataFrame:
    """Explain each column in the output table without including result rows."""
    base = {
        "Run_ID": ("Identifier for the Fiji/Python run.", "Text", "Used to trace rows back to one processed image/run."),
        "Image": ("Original image name exported by Fiji.", "Text", "Useful for checking which microscopy image produced each row."),
        "Mode": ("Analysis mode.", "Text", "Usually 2D, 3D, or time-lapse."),
        "Analysis_Type": ("Standardized analysis category added by Python.", "Text", "Used by interactive viewers and summaries."),
        "Condition": ("Experimental condition inferred from folder structure or metadata.", "Text", "Example: WT, KO, drug-treated, rescue."),
        "Replicate": ("Replicate label inferred from folder structure or metadata.", "Text", "Example: Replicate_1, Replicate_2."),
        "Source_Input": ("Source Fiji output file used by Python.", "Text", "Helps trace combined reports back to individual Fiji outputs."),
        "Source_File": ("Source file name.", "Text", "Used for traceability."),
        "Source_Path": ("Full path to the source file.", "Text", "Useful when troubleshooting combined analyses."),
        "Source_Folder": ("Folder containing the source output.", "Text", "Used for automatic condition/replicate inference."),
        "Cell_ID": ("Cell identifier assigned by Fiji.", "ID", "Cells are numbered within each image/run."),
        "Unique_Cell": ("Python-generated unique cell identifier.", "ID", "Combines source information so cells remain unique after combining files."),
        "Plot_Cell": ("Short display label for cell-level graphs.", "Text", "Used only for graph labels."),
        "Satellite_ID": ("Satellite/granule identifier assigned by Fiji.", "ID", "One row usually corresponds to one detected satellite/granule."),
        "Unique_Satellite": ("Python-generated unique satellite identifier.", "ID", "Keeps satellite identities unique after combining files."),
        "Frame": ("Time-lapse frame number from Fiji.", "Frame", "Only present for time-lapse outputs."),
        "Std_Frame": ("Standardized numeric frame number.", "Frame", "Used by Python time-lapse graphs."),
        "Time_min": ("Time value in minutes if provided/exported.", "Minutes", "Used as time axis when available."),
        "Std_Time_min": ("Standardized time value.", "Minutes", "Used by Python graphs and viewers."),
        "Std_X": ("Satellite X coordinate after calibration.", "µm", "Measured coordinate in the image/spatial viewer."),
        "Std_Y": ("Satellite Y coordinate after calibration.", "µm", "Measured coordinate in the image/spatial viewer."),
        "Std_Z": ("Satellite Z coordinate after calibration.", "µm", "Only available for 3D stacks."),
        "X": ("Raw or Fiji-exported X coordinate.", "Usually pixels or calibrated units", "Standardized to Std_X by Python when possible."),
        "Y": ("Raw or Fiji-exported Y coordinate.", "Usually pixels or calibrated units", "Standardized to Std_Y by Python when possible."),
        "Z": ("Raw or Fiji-exported Z coordinate.", "Usually slice or calibrated units", "Standardized to Std_Z by Python when possible."),
        "Area": ("2D satellite/granule area.", "µm² if calibrated, otherwise pixel²", "Object size measured from the detected 2D mask."),
        "Volume": ("3D satellite/granule volume.", "µm³ if calibrated", "Object size measured from the detected 3D mask."),
        "Std_Size": ("Legacy unified satellite size column.", "µm² for 2D, µm³ for 3D", "Kept for backward compatibility. For mixed 2D/3D experiments, use Std_Area_2D and Std_Volume_3D instead."),
        "Std_Area_2D": ("Standardized 2D satellite area.", "µm²", "Populated only for 2D rows. Remains blank for 3D rows so area and volume are never mixed."),
        "Std_Volume_3D": ("Standardized 3D satellite volume.", "µm³", "Populated only for 3D rows. Remains blank for 2D rows so area and volume are never mixed."),
        "Std_Surface": ("3D satellite surface area when available.", "µm²", "May be produced by 3D Objects Counter outputs."),
        "Sphericity_3D": ("Three-dimensional shape compactness calculated from volume and surface area.", "Dimensionless", "A value of 1 represents a perfect sphere; lower values indicate less spherical or more irregular objects."),
        "Circularity_2D": ("Two-dimensional shape compactness calculated by Fiji.", "Dimensionless", "A value of 1 represents a perfect circle; lower values indicate elongated or irregular objects."),
        "Std_Sphericity": ("Standardized 3D satellite sphericity.", "Dimensionless", "Uses 3D Objects Counter+ Morph_Sphericity when available. A value of 1 is maximally sphere-like. SatelliteQ no longer derives this value from an ambiguous generic Surface column."),
        "Std_Sphericity_Source": ("Source used for standardized 3D sphericity.", "Text", "Normally Morph_Sphericity from 3D Objects Counter+. A formula fallback is used only when an explicitly named calibrated surface-area column exists."),
        "Std_Circularity_2D": ("Standardized 2D satellite circularity.", "Dimensionless", "A value of 1 represents a perfect circle; lower values indicate less circular objects."),
        "Std_Shape_Compactness": ("Unified shape compactness column.", "Dimensionless", "Uses 2D circularity for 2D data and 3D sphericity for 3D data."),
        "Std_Cell_Spatial_Entropy": ("Primary cell-level spatial entropy. It is the mean normalized Shannon entropy across 5x5, 10x10 and 20x20 spatial grids, weighted by each satellite's integrated intensity.", "0–1", "The same value is repeated for every satellite in that cell. A large bright merged aggregate contributes strongly because it carries a large fraction of the cell's total satellite intensity. 0 means concentrated satellite signal; values approaching 1 mean satellite signal is broadly and evenly distributed. The normalization accounts for the observed satellite number."),
        "Std_Cell_Spatial_Entropy_5x5": ("Intensity-weighted normalized cell spatial entropy using a 5x5 grid.", "0–1", "Coarse spatial scale; retained for QC and transparency."),
        "Std_Cell_Spatial_Entropy_10x10": ("Intensity-weighted normalized cell spatial entropy using a 10x10 grid.", "0–1", "Intermediate spatial scale; retained for QC and transparency."),
        "Std_Cell_Spatial_Entropy_20x20": ("Intensity-weighted normalized cell spatial entropy using a 20x20 grid.", "0–1", "Fine spatial scale; retained for QC and transparency."),
        "Std_Cell_Spatial_Entropy_Count_QC": ("Count-weighted multiscale cell spatial entropy.", "0–1", "QC-only comparison in which every detected satellite has equal weight. A large difference from the primary intensity-weighted entropy can indicate that a few large/bright aggregates dominate the cell's satellite material."),
        "Std_Cell_Pericentrosomal_Clustering_Index": ("Fluorescence ratio comparing the centrosome-proximal region with the surrounding shell.", "Dimensionless ratio", "Calculated as total satellite integrated intensity at 0–3 µm divided by total satellite integrated intensity at 3–12 µm from the assigned/nearest centrosome."),
        PCI_STD_INNER_COL: ("Total satellite integrated intensity within 3 µm of the assigned/nearest centrosome.", "A.U.", "Inner component of the pericentrosomal clustering index."),
        PCI_STD_OUTER_COL: ("Total satellite integrated intensity between >3 and <=12 µm from the assigned/nearest centrosome.", "A.U.", "Outer-shell denominator of the pericentrosomal clustering index."),
        PCI_STD_TOTAL_COL: ("Total satellite integrated intensity within 12 µm of the assigned/nearest centrosome.", "A.U.", "Stored for transparency."),
        "Std_Cell_Pericentrosomal_Geometry": ("Geometry used for the pericentrosomal clustering index.", "Text", "2D radial distance or 3D radial distance, according to the analysis mode."),
        "Std_Cell_Pericentrosomal_Clustering_Status": ("Calculation status for the pericentrosomal clustering index.", "Text", "Reports whether the ratio was calculated or why it is undefined, for example when the 3–12 µm denominator is zero."),
        "Std_Cilium_Length_um": ("Longest accepted cilium length assigned to the cell.", "µm", "Fiji derives cilium length from the longest shortest path of the skeletonized cilium mask using a CiliaQ-style workflow."),
        "Std_Ciliated": ("Binary ciliation classification.", "0 or 1", "1 when at least one accepted cilium exceeds the user-defined ciliation length threshold; default threshold is 1 µm."),
        "Std_Distance_to_Cell_or_Image_Center_XY": ("XY distance from each satellite centroid to the geometric center of its assigned cell, or to the image center when cell assignment is not available.", "µm when coordinates are calibrated; otherwise pixels", "For manual/membrane/cytoplasmic cell ROIs, the geometric centroid of the exact Fiji boundary is used. Centrosome-radius/sphere assignments use their defined assignment center. Without cell assignment, the exact image center is used when image dimensions are available."),
        "Std_Center_Reference_X": ("X coordinate of the center used for the center-distance calculation.", "µm or pixels", "Repeated on satellite rows for transparency."),
        "Std_Center_Reference_Y": ("Y coordinate of the center used for the center-distance calculation.", "µm or pixels", "Repeated on satellite rows for transparency."),
        "Std_Center_Reference_Type": ("Description of which center was used.", "Text", "Examples: exact cell-boundary centroid, centrosome-radius/sphere assignment center, exact image center, or estimated image-center fallback."),
        "Std_Center_Distance_Unit": ("Unit of the center-distance measurement.", "Text", "µm for calibrated XY coordinates or pixels for uncalibrated/raw pixel coordinates."),
        "Std_XY_Coordinate_Space": ("Coordinate system used for satellite XY centroids in Python.", "Text", "Usually um for 2D Fiji output and pixel for 3D Objects Counter output. Used internally so spatial entropy uses a cell boundary in the same coordinate system."),
        "Cell_spatial_entropy": ("Cell-level spatial entropy copied once into the per-cell summary.", "0–1", "Use this column for condition comparisons in the Experiment Viewer."),
        "Mean_cell_spatial_entropy": ("Mean cell spatial entropy across cells in one image or frame.", "0–1", "Each cell contributes once, regardless of satellite number."),
        "Cell_pericentrosomal_clustering_index": ("Cell-level pericentrosomal clustering index.", "Dimensionless ratio", "Total satellite integrated intensity within 0–3 µm divided by total satellite integrated intensity in the 3–12 µm shell around the assigned/nearest centrosome(s)."),
        "Mean_cell_pericentrosomal_clustering_index": ("Mean pericentrosomal clustering index across cells in one image or frame.", "Dimensionless ratio", "Each cell contributes once, regardless of satellite number."),
        "Mean_distance_to_cell_or_image_center": ("Mean satellite XY distance to the relevant cell center or image center.", "µm when calibrated; otherwise pixels", "This is the direct average-distance parameter for comparing whether granules lie near or far from the cell/image center."),
        "Mean": ("Mean intensity inside the object.", "A.U.", "Average raw signal per pixel/voxel in the detected object."),
        "Mean_Intensity": ("Mean intensity inside the object.", "A.U.", "Average raw signal per pixel/voxel in the detected object."),
        "Std_Mean_Intensity": ("Standardized mean intensity.", "A.U.", "Used by Python graphs when present."),
        "Integrated_Intensity": ("Sum of pixel/voxel intensities inside the object.", "A.U.", "Also called integrated density/intensity. Larger objects can have larger values because more pixels/voxels are summed."),
        "IntDen": ("Integrated density/intensity exported by Fiji.", "A.U.", "Equivalent to summed signal inside the detected object when available."),
        "Std_Intensity": ("Standardized integrated intensity.", "A.U.", "Used by Python graphs regardless of original Fiji column name."),
        "Distance_to_nearest_centrosome": ("Distance from satellite to nearest centrosome.", "µm if calibrated", "Requires centrosome detection or manual centrosome ROIs."),
        "Distance_to_nearest_centrosome_3D": ("3D distance from satellite to nearest centrosome.", "µm", "Uses X/Y/Z coordinates in 3D mode."),
        "Std_Distance": ("Standardized satellite-to-centrosome distance.", "µm", "Used by Python graphs regardless of original Fiji column name."),
        "Centrosome_ID": ("Centrosome identifier assigned by Fiji.", "ID", "Only present when centrosomes are used."),
        "Std_Centrosome_ID": ("Standardized centrosome identifier.", "ID", "Used by the spatial viewer."),
        "Std_Centrosome_X": ("Centrosome X coordinate.", "µm", "Used by the spatial viewer."),
        "Std_Centrosome_Y": ("Centrosome Y coordinate.", "µm", "Used by the spatial viewer."),
        "Std_Centrosome_Z": ("Centrosome Z coordinate.", "µm", "Used by the spatial viewer for 3D."),
        "Cell_area": ("2D cell area.", "µm² if calibrated", "Measured from cell ROI/mask when cell assignment is used."),
        "Cell_volume": ("3D cell volume or sphere volume.", "µm³", "In 3D sphere mode this can represent the assignment sphere volume."),
        "Cell_measure": ("Cell size measurement exported by Fiji.", "µm² or µm³", "Standardized by Python to Std_Cell_Measure."),
        "Std_Cell_Measure": ("Legacy unified cell area/volume column.", "µm² for 2D, µm³ for 3D", "Kept for backward compatibility."),
        "Std_Cell_Area_2D": ("Standardized 2D cell area.", "µm²", "Populated only for 2D rows."),
        "Std_Cell_Volume_3D": ("Standardized 3D cell volume.", "µm³", "Populated only for 3D rows."),
        "Satellite_count": ("Number of satellites/granules in the group.", "Count", "Used in per-image, per-cell, or time-lapse summaries."),
        "Assigned_satellites": ("Number of satellites assigned to a cell/centrosome.", "Count", "Satellites outside selected ROIs may be unassigned depending on Fiji settings."),
        "Cells": ("Number of cells detected or included.", "Count", "Used in image-level summaries."),
        "Mean_size": ("Legacy mean satellite size in a group.", "µm² or µm³", "Kept for backward compatibility; do not combine 2D area and 3D volume numerically."),
        "Median_size": ("Legacy median satellite size in a group.", "µm² or µm³", "Kept for backward compatibility."),
        "Total_size": ("Legacy sum of satellite sizes in a group.", "µm² or µm³", "Kept for backward compatibility."),
        "Mean_area_2D": ("Mean 2D satellite area in a group.", "µm²", "Calculated only from 2D satellite rows."),
        "Median_area_2D": ("Median 2D satellite area in a group.", "µm²", "Calculated only from 2D satellite rows."),
        "Total_area_2D": ("Total 2D satellite area in a group.", "µm²", "Sum of 2D satellite areas only."),
        "Mean_volume_3D": ("Mean 3D satellite volume in a group.", "µm³", "Calculated only from 3D satellite rows."),
        "Median_volume_3D": ("Median 3D satellite volume in a group.", "µm³", "Calculated only from 3D satellite rows."),
        "Total_volume_3D": ("Total 3D satellite volume in a group.", "µm³", "Sum of 3D satellite volumes only."),
        "Cell_area_2D": ("2D cell area in a per-cell summary.", "µm²", "Blank for 3D cells."),
        "Cell_volume_3D": ("3D cell volume in a per-cell summary.", "µm³", "Blank for 2D cells."),
        "Mean_cell_area_2D": ("Mean 2D cell area in a frame/group.", "µm²", "Uses each cell once."),
        "Total_cell_area_2D": ("Total 2D cell area in a frame/group.", "µm²", "Uses each cell once."),
        "Mean_cell_volume_3D": ("Mean 3D cell volume in a frame/group.", "µm³", "Uses each cell once."),
        "Total_cell_volume_3D": ("Total 3D cell volume in a frame/group.", "µm³", "Uses each cell once."),
        "Mean_intensity": ("Mean integrated intensity in a group.", "A.U.", "Average of satellite integrated intensities."),
        "Total_intensity": ("Sum of integrated intensities in a group.", "A.U.", "Can reflect both satellite number and signal strength."),
        "Mean_distance": ("Mean satellite-to-centrosome distance in a group.", "µm", "Requires centrosome measurements."),
        "Outlier_Cell": ("Whether the cell is marked as an outlier by Python.", "Boolean", "Used only for graph highlighting/filtering."),
        "Outlier_Score": ("Outlier score calculated by Python.", "Numeric", "Higher values indicate more extreme cell-level behavior."),
        "Outlier_Reason": ("Reason a cell was marked as an outlier.", "Text", "For interpretation of graph highlights."),
    }

    rows = []
    for col in columns:
        if col in base:
            meaning, units, notes = base[col]
        elif str(col).startswith("Std_"):
            meaning, units, notes = "Standardized Python column derived from Fiji output.", "Depends on measurement", "Used to make 2D, 3D, time-lapse, and combined files graph consistently."
        elif "Intensity" in str(col) or "intensity" in str(col):
            meaning, units, notes = "Intensity measurement.", "A.U.", "A.U. means arbitrary/raw image intensity units."
        elif "Area" in str(col) or "area" in str(col):
            meaning, units, notes = "Area measurement.", "µm² if calibrated, otherwise pixel²", "Depends on Fiji image calibration or manual pixel size."
        elif "Volume" in str(col) or "volume" in str(col):
            meaning, units, notes = "Volume measurement.", "µm³ if calibrated", "Depends on XY pixel size and Z-step calibration."
        elif "Distance" in str(col) or "distance" in str(col):
            meaning, units, notes = "Distance measurement.", "µm if calibrated", "Usually satellite-to-centrosome distance."
        else:
            meaning, units, notes = "Column exported by Fiji or generated by Python.", "Depends on column", "Kept for traceability, filtering, plotting, or downstream analysis."
        rows.append((col, meaning, units, notes))
    return pd.DataFrame(rows, columns=["Column", "Meaning", "Units_or_range", "Notes"])

def satelliteq_general_explanation_tables(columns: List[str], analysis_kind: str, combined: bool) -> dict:
    readme = pd.DataFrame([
        ("SatelliteQ version", SATELLITEQ_VERSION),
        ("Purpose", "This workbook explains SatelliteQ outputs and metrics. It does not contain quantification result rows."),
        ("Analysis kind", analysis_kind),
        ("Combined analysis", str(combined)),
        ("Use together with", "Use this guide alongside the quantification Excel report, graph viewer, and spatial viewer."),
    ], columns=["Item", "Explanation"])

    workflow = pd.DataFrame([
        ("2D quantification", "Detects satellites/granules from a 2D image or projection; optionally assigns them to cells and centrosomes."),
        ("3D quantification", "Detects satellites/granules in a 3D stack; coordinates, sizes, and distances are interpreted with XY pixel size and Z-step calibration."),
        ("Time-lapse", "Repeats the selected 2D or 3D workflow over frames and creates frame-level summaries when frame information is available."),
        ("Combined reports", "Python combines multiple Fiji outputs and infers Condition and Replicate from the folder structure."),
    ], columns=["Mode", "Explanation"])

    units = pd.DataFrame([
        ("pixel", "2D image element in X/Y.", "Used by Fiji detection and object masks."),
        ("voxel", "3D pixel: one X/Y pixel with Z-slice depth.", "Used internally for 3D object counting."),
        ("µm", "Micrometer.", "Distance and X/Y/Z coordinates when calibrated."),
        ("µm²", "Square micrometer.", "2D satellite/cell area."),
        ("µm³", "Cubic micrometer.", "3D satellite/cell volume."),
        ("A.U.", "Arbitrary/raw image intensity units.", "Intensity values depend on microscope settings, staining, exposure, and preprocessing."),
        ("Integrated intensity", "Sum of all pixel/voxel intensities inside an object.", "Affected by both object size and signal brightness."),
        ("2D circularity", "4π × area / perimeter².", "Dimensionless; 1 is a perfect circle."),
        ("3D sphericity", "π^(1/3) × (6 × volume)^(2/3) / surface area.", "Dimensionless; 1 is a perfect sphere."),
        ("Cell spatial entropy", "Mean of intensity-weighted normalized Shannon entropy calculated at 5x5, 10x10 and 20x20 spatial scales.", "0 = satellite signal concentrated in a small part of the cell; 1 = satellite signal broadly and evenly distributed. Large bright aggregates receive proportionally more weight, and normalization accounts for satellite number."),
        ("Pericentrosomal clustering index", "Total satellite integrated intensity within 0–3 µm divided by total integrated intensity in the 3–12 µm shell around the assigned/nearest centrosome(s).", "Dimensionless."),
        ("Distance to cell/image center", "XY Euclidean distance from each satellite centroid to the geometric cell center, or to image center when no cell is available.", "Per-cell summaries report the average of these satellite distances. Exact Fiji cell boundaries are preferred for cell-center calculation."),
        ("Mean intensity", "Average intensity inside an object.", "Less dependent on object size than integrated intensity."),
    ], columns=["Term", "Meaning", "Notes"])

    folder_structure = pd.DataFrame([
        ("Recommended structure", "Experiment / Condition / Replicate / Fiji_output.csv"),
        ("Condition examples", "WT, KO, Rescue, Drug_A, Control"),
        ("Replicate examples", "Replicate_1, Replicate_2, Replicate_3"),
        ("Condition inference", "Python uses folder names when no metadata CSV is supplied."),
        ("Metadata override", "A metadata CSV can override inferred Condition and Replicate labels when needed."),
    ], columns=["Topic", "Explanation"])

    viewers = pd.DataFrame([
        ("Graph viewer", "Interactive graph tabs for satellite-level, cell-level, condition/replicate, cilium-status/cilium-length, and time-lapse plots."),
        ("Spatial viewer tabs", "Interactive spatial display of satellites, centrosomes, and cell boundaries/spheres using measured coordinates."),
        ("Outlier display", "Outliers are graph annotations for exploration; they do not change Fiji measurements."),
        ("Normalization", "Normalized graphs use a selected reference condition and are intended for visualization/comparison."),
    ], columns=["Viewer_or_feature", "Explanation"])

    return {
        "README": readme,
        "Workflow_modes": workflow,
        "Units_and_terms": units,
        "Column_explanations": satelliteq_result_column_explanations(columns, analysis_kind),
        "Folder_structure": folder_structure,
        "Interactive_viewers": viewers,
    }

def write_satelliteq_explanation_workbook(outdir: Path, columns: List[str], analysis_kind: str, combined: bool, inputs: Optional[List[Path]] = None) -> Optional[Path]:
    try:
        outdir.mkdir(parents=True, exist_ok=True)
        out_xlsx = outdir / satelliteq_dated_name(outdir, "SatelliteQ_Explanation_Guide", inputs=inputs, combined=combined, ext=".xlsx")
        sheets = satelliteq_general_explanation_tables(columns, analysis_kind, combined)
        with pd.ExcelWriter(out_xlsx, engine="openpyxl") as writer:
            for name, table in sheets.items():
                table.to_excel(writer, sheet_name=name[:31], index=False)
        format_excel(out_xlsx)
        print(f"Explanation guide: {out_xlsx}")
        return out_xlsx
    except Exception as e:
        print(f"Explanation guide skipped: {e}")
        return None



def format_excel(path: Path) -> None:
    try:
        from openpyxl import load_workbook
        from openpyxl.styles import Font, PatternFill, Border, Side, Alignment
        from openpyxl.worksheet.table import Table, TableStyleInfo

        wb = load_workbook(path)
        header_fill = PatternFill("solid", fgColor="D9EAF7")
        thin = Side(style="thin", color="DDDDDD")
        border = Border(left=thin, right=thin, top=thin, bottom=thin)

        for ws in wb.worksheets:
            ws.freeze_panes = "A2"
            for cell in ws[1]:
                cell.font = Font(bold=True)
                cell.fill = header_fill
                cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
                cell.border = border
            for row in ws.iter_rows():
                for cell in row:
                    cell.border = border
                    cell.alignment = Alignment(vertical="top")
                    if isinstance(cell.value, float):
                        cell.number_format = "0.000"
            for col_cells in ws.columns:
                col = col_cells[0].column_letter
                max_len = max((len(str(c.value)) for c in col_cells[:100] if c.value is not None), default=10)
                ws.column_dimensions[col].width = min(max(max_len + 2, 10), 45)
            if ws.max_row >= 2 and ws.max_column >= 2:
                try:
                    ref = f"A1:{ws.cell(ws.max_row, ws.max_column).coordinate}"
                    tab = Table(displayName="tbl_" + safe_name(ws.title)[:18], ref=ref)
                    tab.tableStyleInfo = TableStyleInfo(
                        name="TableStyleMedium2",
                        showFirstColumn=False,
                        showLastColumn=False,
                        showRowStripes=True,
                        showColumnStripes=False,
                    )
                    ws.add_table(tab)
                except Exception:
                    pass
        wb.save(path)
    except Exception as e:
        print(f"Excel formatting skipped: {e}")

def write_report(df: pd.DataFrame, outdir: Path, inputs: List[Path], combined: bool) -> Path:
    outdir.mkdir(parents=True, exist_ok=True)
    image_summary = summarize_by_image(df)
    cell_summary = summarize_by_cell(df)
    frame_summary = summarize_by_frame(df)
    cilium_rows, cilium_cells, cilium_skeleton = collect_cilium_companions(df)
    cell_summary = merge_cilium_only_cells_into_summary(cell_summary, cilium_cells)
    # Enforce the biological rule at the final report layer too: no cilium = no length.
    df = sanitize_cilium_length_metrics(df)
    cell_summary = sanitize_cilium_length_metrics(cell_summary)
    spatial_metrics = make_spatial_metrics_table(cell_summary)
    image_summary = update_image_ciliation_from_cell_summary(image_summary, cell_summary)
    graph_paths = make_graphs(df, image_summary, cell_summary, frame_summary, outdir)
    viewer_paths = make_spatial_viewers(df, outdir, cilium_rows, cilium_skeleton)
    all_output_paths = list(graph_paths) + list(viewer_paths)
    overall = summarize_overall(df, input_count=len(inputs), graph_count=len(all_output_paths))

    readme = pd.DataFrame([
        ("SatelliteQ version", SATELLITEQ_VERSION),
        ("What this report is", "SatelliteQ Python report from one or more image outputs."),
        ("Combined report", str(combined)),
        ("Inputs combined", len(inputs)),
        ("Satellite rows", len(df)),
("Graphs and viewers created", len(all_output_paths)),
        ("Units", "Spatial values are reported in µm, µm², or µm³ when Fiji image calibration or manual pixel size was supplied. Intensity is A.U. / raw image intensity units."),
("Cell area / volume", "If present, per-cell summaries and cell area/volume relationship graphs are created."),
        ("Experimental conditions", "Recommended structure: Experiment/Condition/Replicate/Fiji_output.csv. Python infers condition and replicate labels automatically from folders."),
    ], columns=["Item", "Description"])

    input_index = pd.DataFrame({
        "Input_file": [p.name for p in inputs],
        "Input_folder": [str(p.parent) for p in inputs],
        "Full_path": [str(p) for p in inputs],
    })

    graph_index = pd.DataFrame({
        "Output_file": [p.name for p in all_output_paths],
        "Output_folder": [str(p.parent) for p in all_output_paths],
        "Full_path": [str(p) for p in all_output_paths],
    })

    units_legend = pd.DataFrame([
        ("Std_Area_2D", "µm²; populated only for 2D satellite rows"),
        ("Std_Volume_3D", "µm³; populated only for 3D satellite rows"),
        ("Std_Size", "Legacy unified column: µm² for 2D area; µm³ for 3D volume"),
        ("Std_Surface", "Raw 3D surface/surface-area value reported by Fiji; interpretation depends on the source column"),
        ("Std_Sphericity", "Dimensionless; Morph_Sphericity preferred; 1 represents a perfect sphere"),
        ("Std_Sphericity_Source", "Text; identifies the source used for Std_Sphericity"),
        ("Std_Circularity_2D", "Dimensionless; 1 represents a perfect circle"),
        ("Std_Shape_Compactness", "Dimensionless; 2D circularity or 3D sphericity"),
        ("Std_Cell_Spatial_Entropy", "Dimensionless; intensity-weighted multiscale entropy normalized 0–1"),
        ("Std_Cell_Spatial_Entropy_5x5", "Dimensionless; normalized 0–1"),
        ("Std_Cell_Spatial_Entropy_10x10", "Dimensionless; normalized 0–1"),
        ("Std_Cell_Spatial_Entropy_20x20", "Dimensionless; normalized 0–1"),
        ("Std_Cell_Spatial_Entropy_Count_QC", "Dimensionless; normalized 0–1; QC only"),
        ("Std_Cell_Pericentrosomal_Clustering_Index", "Dimensionless; 0–3 µm integrated intensity / 3–12 µm integrated intensity"),
        (PCI_STD_INNER_COL, "A.U.; integrated satellite fluorescence within 0–3 µm"),
        (PCI_STD_OUTER_COL, "A.U.; integrated satellite fluorescence in the 3–12 µm shell"),
        (PCI_STD_TOTAL_COL, "A.U.; integrated satellite fluorescence within 0–12 µm"),
        ("Std_Distance_to_Cell_or_Image_Center_XY", "µm when calibrated; otherwise pixels"),
        ("Std_Cilium_Length_um", "µm; longest accepted cilium length assigned to the cell"),
        ("Std_Cilium_Mean_Length_um", "µm; mean accepted cilium length assigned to the cell"),
        ("Std_Cilium_Count", "Number of accepted cilium objects assigned to the cell"),
        ("Std_Ciliated", "1 = ciliated, 0 = non-ciliated according to the user-defined length threshold"),
        ("Std_Distance", "µm"),
        ("Std_X / Std_Y / Std_Z", "µm coordinates after Fiji calibration/manual pixel-size input"),
        ("Std_Cell_Area_2D", "µm²; populated only for 2D cell rows"),
        ("Std_Cell_Volume_3D", "µm³; populated only for 3D cell rows"),
        ("Std_Cell_Measure", "Legacy unified column: µm² for 2D cell area; µm³ for 3D cell volume"),
        ("Std_Intensity", "A.U. / raw image intensity units"),
        ("Std_Mean_Intensity", "A.U. / raw image intensity units"),
        ("Pixel size source", "Stored in Fiji Setting rows when available: Calibration source, Pixel width um, Pixel height um, Voxel depth um"),
    ], columns=["Column_or_measure", "Unit_or_meaning"])

    out_xlsx = outdir / satelliteq_dated_name(outdir, "SatelliteQ_Quantification_Report", inputs=inputs, combined=combined, ext=".xlsx")

    with pd.ExcelWriter(out_xlsx, engine="openpyxl") as writer:
        readme.to_excel(writer, sheet_name="README", index=False)
        input_index.to_excel(writer, sheet_name="Input_files", index=False)
        overall.to_excel(writer, sheet_name="Summary", index=False)
        df.to_excel(writer, sheet_name="Satellite_rows_clean" if not combined else "Combined_satellite_rows", index=False)
        image_summary.to_excel(writer, sheet_name="Image_summary" if not combined else "Combined_image_summary", index=False)
        if not cell_summary.empty:
            cell_summary.to_excel(writer, sheet_name="Per_cell_summary" if not combined else "Combined_cell_summary", index=False)
        if spatial_metrics is not None and not spatial_metrics.empty:
            spatial_metrics.to_excel(writer, sheet_name="Spatial_metrics" if not combined else "Combined_spatial_metrics", index=False)
        if not frame_summary.empty:
            frame_summary.to_excel(writer, sheet_name="Frame_summary" if not combined else "Combined_frame_summary", index=False)
        if not cilium_rows.empty:
            cilium_rows.to_excel(writer, sheet_name="Cilium_rows", index=False)
        if not cilium_cells.empty:
            cilium_cells.to_excel(writer, sheet_name="Cilium_cell_summary", index=False)
        if not cilium_skeleton.empty:
            cilium_skeleton.to_excel(writer, sheet_name="Cilium_skeleton_points", index=False)
        units_legend.to_excel(writer, sheet_name="Units", index=False)
        graph_index.to_excel(writer, sheet_name="Graph_and_viewer_files", index=False)

    format_excel(out_xlsx)
    explanation_columns = sorted(set(df.columns).union(image_summary.columns).union(cell_summary.columns).union(frame_summary.columns).union(spatial_metrics.columns if spatial_metrics is not None else []))
    write_satelliteq_explanation_workbook(outdir, explanation_columns, "combined quantification" if combined else "single-file quantification", combined=combined, inputs=inputs)
    return out_xlsx

def process_one(input_path: Path) -> Optional[Path]:
    df = read_input(input_path)
    if df is None or df.empty:
        print(f"No usable satellite rows found in: {input_path}")
        return None
    df = attach_cilium_measurements(df, input_path)
    df, meta_path = apply_experiment_metadata(df, [input_path], input_path.parent)
    if meta_path:
        print(f"Applied experiment metadata: {meta_path.name}")
    outdir = input_path.parent / f"{safe_name(input_path.stem)}_Python_Report"
    return write_report(df, outdir, [input_path], combined=False)

def process_combined(inputs: List[Path], root: Path) -> Optional[Path]:
    dfs = []
    good = []
    for f in inputs:
        df = read_input(f)
        if df is not None and not df.empty:
            df = attach_cilium_measurements(df, f)
            dfs.append(df)
            good.append(f)
    if not dfs:
        print("No usable satellite rows found.")
        return None
    combined_df = pd.concat(dfs, ignore_index=True, sort=False)
    combined_df, meta_path = apply_experiment_metadata(combined_df, good, root)
    if meta_path:
        print(f"Applied experiment metadata: {meta_path.name}")
    else:
        print("Using folder names for Condition and Replicate labels.")
        print("Recommended structure: Experiment/Condition/Replicate/Fiji_output.csv.")
    outdir = root / "SatelliteQ_Combined_Report"
    return write_report(combined_df, outdir, good, combined=True)

def main() -> None:
    if len(sys.argv) > 1:
        selected = Path(sys.argv[1].strip().strip('"'))
    else:
        selected = choose_input_gui()

    if selected is None:
        typed = input("Paste CSV/Excel file path or parent folder path: ").strip().strip('"')
        if not typed:
            print("No path selected.")
            return
        selected = Path(typed)

    if not selected.exists():
        print(f"Path not found: {selected}")
        return

    inputs = find_candidate_files(selected)
    if not inputs:
        print("No compatible SatelliteQ CSV/TSV/Excel files found.")
        print("For 3D Fiji results, accepted names include:")
        print("  ImageName_3D_2026_7_8_16_42_10_results.csv")
        print("  ImageName_3D_2026_7_8_16_42_10_results.xlsx")
        print("  ImageName_3D_2026_7_8_16_42_10_results_excel_friendly.tsv")
        return

    quant_inputs = inputs
    outputs = []

    if len(quant_inputs) == 1:
        print(f"Found one segmentation/quantification file: {quant_inputs[0].name}")
        out = process_one(quant_inputs[0])
        if out:
            outputs.append(out)
        print("\nDone.")
        for o in outputs:
            print(f"Report: {o}")
        return

    mode = ask_multiple_mode(len(quant_inputs))
    print(f"Found {len(quant_inputs)} segmentation/quantification files.")
    print(f"Mode: {mode}")

    if mode in ["separate", "both"]:
        for f in quant_inputs:
            print(f"Processing separately: {f.name}")
            out = process_one(f)
            if out:
                outputs.append(out)

    if mode in ["combine", "both"]:
        root = selected if selected.is_dir() else selected.parent
        print("Creating combined segmentation/quantification report...")
        out = process_combined(quant_inputs, root)
        if out:
            outputs.append(out)

    print("\nDone.")
    for o in outputs:
        print(f"Report: {o}")

if __name__ == "__main__":
    main()
