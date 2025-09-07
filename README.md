# AF3 One‑Zip Analyzer

Interactive tool to analyze AlphaFold 3 ZIP exports (one ZIP at a time), compute close contacts between two chains, optionally run PDBe Arpeggio to classify interactions, and generate plots + Residue Interaction Network (RIN) images. It works whether your ZIP contains a single AF3 run or multiple run folders.

---

## What this is

AlphaFold 3 typically gives you a ZIP that contains one or more “run” folders with CIFs like:

- `fold_<name>_model_0.cif` … `fold_<name>_model_4.cif`

This tool lets you:

- Pick a ZIP (quick list in terminal, or a file dialog).
- Inspect chain stats from the first model so you can confidently map **which chain is which**.
- Choose the **interaction type**: *Peptide–Protein*, *Protein–Protein*, or *Ligand–Receptor*.
- Accept sensible defaults (e.g., smaller chain → peptide) or rename groups to match your biology.
- Compute **close contacts**, optionally run **PDBe Arpeggio** to type interactions.
- Produce per-model **plots** + **RIN** images and a per-run **summary CSV**.

All outputs are written next to the CIFs inside the extracted ZIP tree so you always know where files go.

---

## Repo contents

Keep these in your repo root:

- `Analysis_Script.py` — the menu‑driven analyzer
- `Classes.py` — Atom parser used by the analyzer

Optional example ZIPs for demos:

- `fold_mouse_serca314_757_trim59.zip` — single‑run example
- `folds_2025_09_06_20_45.zip` — multi‑run example

Assets referenced by this README (in `./assets/`):

```
assets/
├── closecontacts.png
├── mac_zip_empty.png
├── multi_menu.png
├── multi_run1_chains.png
├── multi_run1_final.png
├── multi_run1_mapping.png
├── multi_run2_mapping.png
├── multi_run2_reuse.png
├── rin.png
├── solo_chains.png
├── solo_final.png
├── solo_mapping.png
└── solo_menu.png
```

---

## One‑minute install (non‑interactive)

Pinned to Python 3.11 because Open Babel wheels often lag the newest Python versions.

```bash
conda install -y -n base -c conda-forge mamba
mamba create -y -n af3-arpeggio -c conda-forge   python=3.11 openbabel biopython gemmi numpy pandas matplotlib networkx tk pip
conda activate af3-arpeggio
pip install -q pdbe-arpeggio
```

Or use the curated env file:

```bash
mamba env create -y -f environment.yml
conda activate af3-arpeggio
```

Windows/macOS/Linux are all supported. The environment includes `tk` for the file dialog. If your Linux image is missing Tk at the system level, install your distro’s Tk package (the conda one usually suffices).

---

## Run an analysis

```bash
python Analysis_Script.py
```

You’ll see a menu like:

```
=== AF3 One-Zip Analyzer ===
 1) fold_mouse_serca314_757_trim59.zip
 2) folds_2025_09_06_20_45.zip
 0) Browse with file dialog
 t) Type/paste a path
 r) Refresh list
 q) Quit
```

Press `1` to analyze the single‑run example, or `2` for the multi‑run ZIP. `0` opens a file dialog; `t` lets you paste a path.

---

## Example A — single‑run ZIP (`fold_mouse_serca314_757_trim59.zip`)

1. Pick the ZIP. It extracts to `./fold_mouse_serca314_757_trim59/`.
2. The first CIF converts to PDB and you see the chain table. This helps assign A/B.
3. Choose interaction type (e.g., Protein–Protein) and confirm/edit the suggested mapping. Defaults:
   - Protein–Protein → the two largest chains
   - Peptide–Protein → smaller chain is the peptide
   - Ligand–Receptor → a small HET‑like chain is suggested as ligand
4. Confirm names and the distance threshold (e.g., **5.0 Å**). The script processes each model.

Per‑model outputs (saved next to each CIF):

- `fold_<name>_model_N.pdb` — CIF → PDB (fast parsing)
- `fold_<name>_model_N_close_contacts.{pdb,cif}` — atoms within cutoff between the two selected chains
- `fold_<name>_model_N_contact_plot.png` — per‑group unique close‑contact atom counts
- `Residue_Interaction_Network_Model_N.png` — RIN colored by interaction type (if typing available)

If Arpeggio is installed:

- `*_close_contacts.json` — interaction annotations
- `*_close_contacts.csv` — CSV converted by the script

Per‑run summary:

- `interaction_summary_across_models.csv` — interaction types per residue pair across models

---

## Example B — multi‑run ZIP (`folds_2025_09_06_20_45.zip`)

Contains multiple runs (e.g., `mouse_serca111_253_nnat`, `mouse_serca314_757_nnat`, `mouse_serca314_757_trim59`). The tool prompts you per run.

1. Choose a global distance threshold (default **5.0 Å**).
2. For Run 1, see chain stats and pick type → map chains.
3. Confirm labels and process models.
4. For later runs, optionally **reuse** the previous run’s settings if chains line up.
5. Ligand–Receptor picks a small chain as ligand; rename either side as needed.

Each run writes outputs in its own subfolder beneath the extracted ZIP.

---

## Where files go

The ZIP is extracted to a folder with the same name (minus `.zip`). Each run folder contains its CIFs and outputs. Example:

```
./folds_2025_09_06_20_45/
├─ mouse_serca111_253_nnat/
│  ├─ fold_mouse_serca111_253_nnat_model_0.cif
│  ├─ fold_mouse_serca111_253_nnat_model_0.pdb
│  ├─ fold_mouse_serca111_253_nnat_model_0_close_contacts.pdb
│  ├─ fold_mouse_serca111_253_nnat_model_0_close_contacts.cif
│  ├─ fold_mouse_serca111_253_nnat_model_0_contact_plot.png
│  ├─ Residue_Interaction_Network_Model_0.png
│  └─ interaction_summary_across_models.csv
├─ mouse_serca314_757_nnat/
└─ mouse_serca314_757_trim59/
```

---

## Requirements recap

Minimum (without Arpeggio typing):

- Python 3.10–3.11
- `biopython` (MMCIFParser, PDBIO)
- `gemmi`
- `numpy`, `pandas`
- `matplotlib`
- `networkx`
- `tk` (file dialog)

Optional (for interaction typing):

- `openbabel` (C++ lib + Python bindings)
- `pdbe-arpeggio` (CLI)

Why 3.11? Open Babel wheels/conda builds often trail the newest Python releases; 3.13 may not have wheels on your OS yet.

---

## Curated conda environment (portable)

Example `environment.yml` to keep in your repo:

```yaml
name: af3-arpeggio
channels: [conda-forge]
dependencies:
  - python=3.11
  - openbabel
  - biopython
  - gemmi
  - numpy
  - pandas
  - matplotlib
  - networkx
  - tk
  - pip
  - pip:
      - pdbe-arpeggio
```

Create & activate:

```bash
mamba env create -y -f environment.yml
conda activate af3-arpeggio
```

---

## Pip‑only option (no conda)

```bash
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
pip install --upgrade pip
pip install biopython gemmi numpy pandas matplotlib networkx
pip install openbabel-wheel      # may not exist for latest Python/OS
pip install pdbe-arpeggio
```

If Open Babel doesn’t install, you can still run the analyzer (contacts/plots/RIN) — Arpeggio typing is simply skipped.

---

## Troubleshooting

- **`ModuleNotFoundError: openbabel`** when Arpeggio runs  
  Install Open Babel via conda-forge: `mamba install -y -c conda-forge openbabel`.  
  Or disable typing by commenting out the `run_arpeggio(...)` call.

- **Arpeggio CLI not found**  
  `pip install -q pdbe-arpeggio`.

- **Performance on very large complexes**  
  The distance loop is O(N×M). For a speed boost, replace it with a KD‑Tree (`scipy.spatial.cKDTree`).

- **macOS “archive is empty”**  
  Ensure you created a standard ZIP. The sample archive in this package extracts cleanly with macOS Archive Utility.

---

## Reproducible environments: which file do I use?

- **Use the curated `environment.yml`** when sharing with colleagues or CI (portable, not over‑pinned).
- To clone your exact workstation (same OS/arch), export your current env:

```bash
conda env export --no-builds > environment.exported.yml
conda env export --from-history > environment.from-history.yml
```

Most teams standardize on the curated `environment.yml`.
