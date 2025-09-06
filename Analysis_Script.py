#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import zipfile
import logging
import subprocess
import json
from math import sqrt
from typing import Dict, List, Tuple, Optional, Set

# Third-party
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
import gemmi
from Bio.PDB import MMCIFParser, PDBIO

# Project
from Classes import Atom

# ---------------- Logging ----------------
logging.basicConfig(
    filename='interaction_analysis.log',
    level=logging.DEBUG,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
log = logging.getLogger("af3-onezip-menu")
log.addHandler(logging.StreamHandler(sys.stdout))

# ---------------- Defaults ----------------
DISTANCE_DEFAULT = 5.0
INTERACTION_COLORS = {
    "vdw_clash": "red",
    "weak_polar": "blue",
    "hydrophobic": "green",
    "polar": "orange",
    "ionic": "purple",
    "hbond": "cyan",
    "vdw": "gray",
    "covalent": "pink"
}
plt.rcParams.update({"figure.autolayout": True})

# ---------------- Utils ----------------
def ensure_dir(path: str):
    os.makedirs(path, exist_ok=True)

def list_local_zips(cwd=".") -> List[str]:
    try:
        return sorted([
            os.path.abspath(f)
            for f in os.listdir(cwd)
            if f.lower().endswith(".zip") and os.path.isfile(f)
        ])
    except Exception:
        return []

def choose_zip_menu() -> Optional[str]:
    """Menu that lists zips; 0 = browse (Tk), t = typed path, q = quit."""
    zips = list_local_zips(".")
    print("\n=== AF3 One-Zip Analyzer ===")
    print("Pick a ZIP from the list, or choose an option:")
    if zips:
        for i, z in enumerate(zips, start=1):
            print(f" {i}) {os.path.basename(z)}")
    else:
        print(" (No .zip files found in this folder.)")
    print(" 0) Browse with file dialog")
    print(" t) Type/paste a path")
    print(" r) Refresh list")
    print(" q) Quit")

    while True:
        choice = input("Choice: ").strip().lower()
        if choice == "q":
            return None
        if choice == "r":
            return choose_zip_menu()
        if choice == "0":
            try:
                import tkinter as tk
                from tkinter import filedialog
                root = tk.Tk(); root.withdraw(); root.update_idletasks()
                path = filedialog.askopenfilename(
                    title="Select AF3 ZIP",
                    filetypes=[("Zip files", "*.zip"), ("All files", "*.*")]
                )
                return path or None
            except Exception as e:
                log.warning(f"File dialog unavailable ({e}). Fallback to typed path.")
                choice = "t"
        if choice == "t":
            path = input("ZIP path: ").strip()
            return path or None
        try:
            idx = int(choice)
            if 1 <= idx <= len(zips):
                return zips[idx-1]
        except Exception:
            pass
        print("Invalid choice. Enter a number from the list, 0, t, r, or q.")

def extract_zip_preserve_tree(zip_path: str) -> str:
    base = os.path.splitext(os.path.basename(zip_path))[0]
    out_dir = os.path.abspath(base)
    ensure_dir(out_dir)
    with zipfile.ZipFile(zip_path, 'r') as zf:
        zf.extractall(out_dir)
    print(f"[extract] {zip_path} → {out_dir}")
    return out_dir

def find_cifs(root: str) -> List[str]:
    cifs = []
    for dp, _, files in os.walk(root):
        for fn in files:
            if fn.lower().endswith(".cif"):
                cifs.append(os.path.join(dp, fn))
    return sorted(cifs)

def find_run_dirs(root: str) -> List[str]:
    """AF3 runs = first-level dirs containing any CIFs; if CIFs live at root, include root."""
    run_dirs = []
    for d in sorted(os.listdir(root)):
        p = os.path.join(root, d)
        if os.path.isdir(p) and find_cifs(p):
            run_dirs.append(p)
    if any(fn.lower().endswith(".cif") for fn in os.listdir(root)):
        run_dirs.insert(0, root)
    return run_dirs or [root]

# ---------------- CIF/PDB helpers ----------------
def convert_cif_to_pdb(cif_file: str, output_pdb_file: str):
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("", cif_file)
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb_file)

def convert_pdb_to_cif(pdb_file: str, output_cif_file: str):
    try:
        st = gemmi.read_structure(pdb_file)
        st.setup_entities()
        block = st.make_mmcif_block()
        with open(output_cif_file, "w") as fw:
            fw.write(block.as_string())
    except Exception as e:
        log.error(f"Error converting {pdb_file} to CIF: {e}")

def parse_pdb_atoms(filename: str, chains: Optional[Set[str]]=None) -> Dict[str, Atom]:
    atoms = {}
    allow = None if chains is None else set(chains)
    with open(filename, "r") as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                a = Atom(pdbfileline=line)
                if allow is None or a.chain_identifier in allow:
                    atoms[a.atom_key] = a
    return atoms

def count_chain_stats_from_pdb(pdb_file: str) -> Dict[str, Dict]:
    stats: Dict[str, Dict] = {}
    with open(pdb_file, "r") as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            a = Atom(pdbfileline=line)
            ch = a.chain_identifier
            resname = a.residue_name
            if ch not in stats:
                stats[ch] = {
                    "res_ids": set(),
                    "atoms": 0,
                    "het_non_hoh": set(),
                    "water_count": 0
                }
            stats[ch]["atoms"] += 1
            if resname == "HOH":
                stats[ch]["water_count"] += 1
            else:
                stats[ch]["res_ids"].add((a.residue_sequence_number, resname))
                if line.startswith("HETATM"):
                    stats[ch]["het_non_hoh"].add((a.residue_sequence_number, resname))
    for ch, d in stats.items():
        d["residues"] = len(d["res_ids"])
        d["het_residues"] = len(d["het_non_hoh"])
        d["is_water_only"] = (d["residues"] == 0 and d["water_count"] > 0)
    return stats

def print_chain_table(stats: Dict[str, Dict]):
    print("\nChains detected (first model in this run):")
    print(f"{'Chain':<6} {'Residues':>8} {'Atoms':>8} {'HET res':>8} {'WaterOnly':>10}")
    for ch, d in sorted(stats.items(), key=lambda kv: kv[0]):
        print(f"{ch:<6} {d['residues']:>8} {d['atoms']:>8} {d['het_residues']:>8} {str(d['is_water_only']):>10}")

# ---------------- Interaction prompts ----------------
def choose_interaction_type() -> str:
    print("\nPick interaction type:")
    print(" 1) Peptide–Protein")
    print(" 2) Protein–Protein")
    print(" 3) Ligand–Receptor")
    while True:
        s = input("Choice [1/2/3]: ").strip()
        if s in {"1","2","3"}:
            return {"1":"pep-pro","2":"prot-prot","3":"lig-rec"}[s]

def prompt_yes_no(msg: str, default_yes: bool = True) -> bool:
    prompt = " [Y/n]: " if default_yes else " [y/N]: "
    s = input(msg + prompt).strip().lower()
    if s == "":
        return default_yes
    return s.startswith("y")


def suggest_chains(kind: str, stats: Dict[str, Dict]) -> Tuple[str, str, str, str]:
    clean = {ch: d for ch, d in stats.items() if not d["is_water_only"]}
    ordered = sorted(clean.items(), key=lambda kv: kv[1]["residues"], reverse=True)
    if not ordered:
        raise RuntimeError("No non-water chains detected.")

    if kind == "pep-pro":
        smallest = min(ordered, key=lambda kv: kv[1]["residues"])[0]
        largest  = max(ordered, key=lambda kv: kv[1]["residues"])[0]
        return smallest, largest, "Peptide", "Protein"

    if kind == "prot-prot":
        if len(ordered) >= 2:
            a = ordered[0][0]; b = ordered[1][0]
        else:
            a = ordered[0][0]; b = ordered[0][0]
        return a, b, "Protein A", "Protein B"

    ligand_candidates = [ch for ch,d in clean.items() if d["residues"] <= 3 and d["atoms"] <= 150]
    if ligand_candidates:
        lig = sorted(ligand_candidates, key=lambda c: (clean[c]["residues"], clean[c]["atoms"]))[0]
    else:
        lig = min(ordered, key=lambda kv: kv[1]["residues"])[0]
    rec = max(ordered, key=lambda kv: kv[1]["residues"])[0]
    if lig == rec and len(ordered) >= 2:
        rec = ordered[1][0]
    return lig, rec, "Ligand", "Receptor"

def prompt_chain_names(defaults: Tuple[str,str,str,str], stats: Dict[str,Dict]) -> Tuple[str,str,str,str]:
    g1_chain, g2_chain, g1_name, g2_name = defaults
    print(f"\nSuggested mapping:")
    print(f"  Group 1: chain [{g1_chain}] → name [{g1_name}]")
    print(f"  Group 2: chain [{g2_chain}] → name [{g2_name}]")
    while True:
        x = input(f"Enter Group 1 chain (blank = {g1_chain}): ").strip() or g1_chain
        y = input(f"Enter Group 2 chain (blank = {g2_chain}): ").strip() or g2_chain
        if x not in stats or y not in stats:
            print(" One or both chains not found in this run. Try again.")
            continue
        if x == y:
            print(" Chains must be different. Try again.")
            continue
        break
    n1 = input(f"Name for Group 1 (blank = {g1_name}): ").strip() or g1_name
    n2 = input(f"Name for Group 2 (blank = {g2_name}): ").strip() or g2_name
    return x, y, n1, n2

def prompt_distance(default_value: float) -> float:
    try:
        s = input(f"Distance threshold Å [default {default_value}]: ").strip()
        return float(s) if s else default_value
    except Exception:
        return default_value

# ---------------- Contacts & plots ----------------
def find_close_atoms(group1_atoms: List[Atom], group2_atoms: List[Atom], distance_threshold: float) -> Set[Atom]:
    close = set()
    for a1 in group1_atoms:
        for a2 in group2_atoms:
            dx = (a1.x - a2.x) ** 2
            dy = (a1.y - a2.y) ** 2
            dz = (a1.z - a2.z) ** 2
            if (dx + dy + dz) ** 0.5 < distance_threshold:
                close.add(a1); close.add(a2)
    return close

def plot_contact_bar(unique_g1: Set[str], unique_g2: Set[str], out_png: str, labels: Tuple[str, str], threshold: float):
    data = [len(unique_g1), len(unique_g2)]
    x = np.arange(2)
    fig, ax = plt.subplots()
    ax.bar(x, data, 0.5)
    ax.set_ylabel(f'Atoms within {threshold} Å')
    ax.set_xticks(x)
    ax.set_xticklabels([labels[0], labels[1]])
    ax.set_title("Close-contact atom counts")
    plt.savefig(out_png, dpi=300)
    plt.close(fig)

# ---------------- Arpeggio & post ----------------
def run_arpeggio(input_cif: str, output_dir: str):
    for exe in ("pdbe-arpeggio", "arpeggio"):
        try:
            subprocess.run([exe, input_cif, "-o", output_dir], check=True)
            return
        except FileNotFoundError:
            continue
        except subprocess.CalledProcessError as e:
            log.error(f"{exe} failed for {input_cif}: {e}")
            return
    log.warning("Arpeggio CLI not found on PATH; skipping interaction typing step.")

def process_json_to_csv(json_dir: str):
    for file_name in os.listdir(json_dir):
        if file_name.endswith("_close_contacts.json"):
            json_file = os.path.join(json_dir, file_name)
            csv_file = os.path.join(json_dir, file_name.replace(".json", ".csv"))
            try:
                with open(json_file, 'r') as f:
                    data = json.load(f)
                rows = []
                for interaction in data:
                    bgn = interaction.get('bgn', {})
                    end = interaction.get('end', {})
                    contacts = interaction.get('contact', [])
                    if not bgn or not end or "proximal" in contacts:
                        continue
                    bgn_chain = bgn.get('auth_asym_id')
                    end_chain = end.get('auth_asym_id')
                    if bgn_chain == end_chain or bgn.get('label_comp_id') == "HOH" or end.get('label_comp_id') == "HOH":
                        continue
                    for contact_type in contacts:
                        rows.append({
                            "Ligand": bgn.get('label_comp_id'),
                            "Ligand_Chain": bgn_chain,
                            "Ligand_Number": bgn.get('auth_seq_id'),
                            "Ligand_Atom": bgn.get('auth_atom_id'),
                            "Contact_Type": contact_type,
                            "Distance": interaction.get('distance', ''),
                            "Residue": end.get('label_comp_id'),
                            "Residue_Chain": end_chain,
                            "Residue_Number": end.get('auth_seq_id'),
                            "Residue_Atom": end.get('auth_atom_id'),
                        })
                if rows:
                    pd.DataFrame(rows).to_csv(csv_file, index=False)
            except Exception as e:
                log.error(f"Error processing JSON file {json_file}: {e}")

def resolve_overlaps(pos, radius=0.2, max_iterations=50):
    nodes = list(pos.keys())
    for _ in range(max_iterations):
        adjusted = False
        for i, n1 in enumerate(nodes):
            for j, n2 in enumerate(nodes):
                if i >= j:
                    continue
                x1, y1 = pos[n1]
                x2, y2 = pos[n2]
                dx, dy = x2 - x1, y2 - y1
                dist = (dx**2 + dy**2) ** 0.5
                if dist < radius and dist > 1e-6:
                    adjusted = True
                    offset = (radius - dist) / 2
                    dx, dy = dx/dist * offset, dy/dist * offset
                    pos[n1] = (x1 - dx, y1 - dy)
                    pos[n2] = (x2 + dx, y2 + dy)
        if not adjusted:
            break

def generate_rin(csv_dir: str):
    for fn in os.listdir(csv_dir):
        if fn.endswith("_close_contacts.csv"):
            try:
                model_number = fn.split("model_")[-1].split("_")[0]
                df = pd.read_csv(os.path.join(csv_dir, fn))
                G = nx.Graph()
                for _, row in df.iterrows():
                    r1 = f"{row['Residue']} {row['Residue_Number']}{row['Residue_Chain']}"
                    r2 = f"{row['Ligand']} {row['Ligand_Number']}{row['Ligand_Chain']}"
                    G.add_node(r1, chain=row['Residue_Chain'])
                    G.add_node(r2, chain=row['Ligand_Chain'])
                    G.add_edge(r1, r2, interaction=row['Contact_Type'], distance=row['Distance'])
                pos = nx.spring_layout(G, k=2.0, iterations=300)
                resolve_overlaps(pos, radius=0.5)
                plt.figure(figsize=(16, 14))
                node_colors = ["#6fa8dc" if G.nodes[n]['chain'] == "A" else "#93c47d" for n in G.nodes]
                node_sizes = [500 + 100 * G.degree[n] for n in G.nodes]
                nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color=node_colors, edgecolors="black")
                for e in G.edges(data=True):
                    edge_color = INTERACTION_COLORS.get(e[2]['interaction'], "black")
                    nx.draw_networkx_edges(G, pos, edgelist=[(e[0], e[1])], edge_color=edge_color, width=2, alpha=0.7)
                nx.draw_networkx_labels(G, pos, font_size=10, font_color="black")
                legend_elements = [plt.Line2D([0], [0], color=color, lw=2, label=interaction)
                                   for interaction, color in INTERACTION_COLORS.items()]
                plt.legend(handles=legend_elements, loc="upper left", title="Interaction Types")
                plt.title(f"Residue Interaction Network for Model {model_number}")
                out_png = os.path.join(csv_dir, f"Residue_Interaction_Network_Model_{model_number}.png")
                plt.axis("off")
                plt.savefig(out_png, bbox_inches="tight", dpi=300)
                plt.close()
            except Exception as e:
                log.error(f"RIN error for {fn}: {e}")

def interaction_summary(csv_dir: str):
    all_frames = []
    for fn in os.listdir(csv_dir):
        if fn.endswith("_close_contacts.csv"):
            try:
                model_number = fn.split("model_")[-1].split("_")[0]
                df = pd.read_csv(os.path.join(csv_dir, fn))
                df['Model'] = f"Model_{model_number}"
                df['Residue_Pair'] = df.apply(
                    lambda row: '-'.join(sorted([
                        f"{row['Residue']}{row['Residue_Number']}{row['Residue_Chain']}",
                        f"{row['Ligand']}{row['Ligand_Number']}{row['Ligand_Chain']}"
                    ])),
                    axis=1
                )
                all_frames.append(df)
            except Exception as e:
                log.error(f"Summary read error for {fn}: {e}")

    if not all_frames:
        return
    combined = pd.concat(all_frames, ignore_index=True)
    grouped = combined.groupby(['Residue_Pair', 'Contact_Type'])['Model'].apply(set).reset_index()
    summary = grouped.groupby('Residue_Pair').agg({'Contact_Type': list, 'Model': list}).reset_index()
    out_csv = os.path.join(csv_dir, "interaction_summary_across_models.csv")
    summary.to_csv(out_csv, index=False)
    print(f"[summary] Saved: {out_csv}")

def generate_graph(unique_atoms_g1: Set[str], unique_atoms_g2: Set[str], pdb_file_path: str,
                   distance_threshold: float, labels: Tuple[str,str]):
    data = [len(unique_atoms_g1), len(unique_atoms_g2)]
    x = np.arange(2)
    fig, ax = plt.subplots()
    ax.bar(x, data, 0.5)
    ax.set_ylabel(f'Number of Atoms within {distance_threshold} Å')
    ax.set_xticks(x)
    ax.set_xticklabels([labels[0], labels[1]])
    plt.savefig(os.path.splitext(pdb_file_path)[0] + "_contact_plot.png", dpi=300)
    plt.close(fig)

# ---------------- Per-CIF processing ----------------
def process_cif(cif_path: str, g1_chain: str, g2_chain: str, g1_name: str, g2_name: str, distance: float):
    dp = os.path.dirname(cif_path)
    base = os.path.splitext(os.path.basename(cif_path))[0]
    pdb_path = os.path.join(dp, base + ".pdb")
    convert_cif_to_pdb(cif_path, pdb_path)

    atoms = parse_pdb_atoms(pdb_path, {g1_chain, g2_chain})
    g1_atoms = [a for a in atoms.values() if a.chain_identifier == g1_chain]
    g2_atoms = [a for a in atoms.values() if a.chain_identifier == g2_chain]

    close_atoms = find_close_atoms(g1_atoms, g2_atoms, distance)

    cc_pdb = os.path.join(dp, base + "_close_contacts.pdb")
    with open(cc_pdb, "w") as fw:
        for a in close_atoms:
            fw.write(a.__repr__())

    cc_cif = os.path.join(dp, base + "_close_contacts.cif")
    convert_pdb_to_cif(cc_pdb, cc_cif)

    unique_g1 = {a.atom_key for a in close_atoms if a.chain_identifier == g1_chain}
    unique_g2 = {a.atom_key for a in close_atoms if a.chain_identifier == g2_chain}
    generate_graph(unique_g1, unique_g2, pdb_path, distance, (g1_name, g2_name))

    run_arpeggio(cc_cif, dp)

# ---------------- Driver ----------------
def main():
    zip_path = choose_zip_menu()
    if not zip_path:
        print("No ZIP selected. Exiting.")
        return
    if not os.path.exists(zip_path):
        print("File not found. Exiting.")
        return

    root_dir = extract_zip_preserve_tree(zip_path)
    run_dirs = find_run_dirs(root_dir)
    if not run_dirs:
        print("No AF3 runs found in ZIP. Exiting.")
        return

    distance_default = prompt_distance(DISTANCE_DEFAULT)
    prev_settings = None  # (g1_chain, g2_chain, g1_name, g2_name, kind, distance)

    for i, run_dir in enumerate(run_dirs, start=1):
        print("\n" + "="*80)
        print(f"RUN {i}/{len(run_dirs)}: {os.path.relpath(run_dir, root_dir)}")
        print("="*80)

        cifs = find_cifs(run_dir)
        if not cifs:
            print("  (No CIFs here; skipping.)")
            continue

        # Chain stats from first CIF
        first_cif = cifs[0]
        tmp_pdb = os.path.splitext(first_cif)[0] + ".pdb"
        convert_cif_to_pdb(first_cif, tmp_pdb)
        stats = count_chain_stats_from_pdb(tmp_pdb)
        print_chain_table(stats)

        # Reuse previous settings?
        if prev_settings:
            reuse = input("Reuse previous run's settings? [Y/n]: ").strip().lower() or "y"
        else:
            reuse = "n"

        if reuse.startswith("y"):
            g1_chain, g2_chain, g1_name, g2_name, kind, distance = prev_settings
            if g1_chain not in stats or g2_chain not in stats or g1_chain == g2_chain:
                print(" Previous chains not valid for this run; re-prompting.")
                reuse = "n"

        if not reuse.startswith("y"):
            # NEW:
            kind = choose_interaction_type()
            defaults = suggest_chains(kind, stats)
            g1_chain, g2_chain, g1_name, g2_name = prompt_chain_names(defaults, stats)

            if prev_settings is None:
                # First run: just use the global value you already entered
                distance = distance_default
                print(f"Using distance {distance} Å for this run.")
            else:
                # Subsequent runs: offer an override, default is to keep current
                if prompt_yes_no(f"Keep distance {distance_default} Å for this run?", default_yes=True):
                    distance = distance_default
                else:
                    distance = prompt_distance(distance_default)
                    # Optionally make the new value the default for later runs
                    if prompt_yes_no("Use this as the default for subsequent runs?", default_yes=True):
                        distance_default = distance


        prev_settings = (g1_chain, g2_chain, g1_name, g2_name, kind, distance)

        print("\nFinal mapping for this run:")
        print(f"  {g1_name}: chain {g1_chain}")
        print(f"  {g2_name}: chain {g2_chain}")
        print(f"  Distance: {distance} Å")
        print("  Processing models...")

        for cif in cifs:
            try:
                process_cif(cif, g1_chain, g2_chain, g1_name, g2_name, distance)
            except Exception as e:
                log.error(f"Top-level error for {cif}: {e}")
                print(f"[error] {os.path.basename(cif)} failed; see interaction_analysis.log")

        # Per-run post
        process_json_to_csv(run_dir)
        generate_rin(run_dir)
        interaction_summary(run_dir)

    print("\nAll done. Outputs are saved next to the CIFs inside the extracted ZIP tree.")
    print(f"Root: {root_dir}")

if __name__ == "__main__":
    main()
