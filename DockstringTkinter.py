import subprocess
import tkinter as tk
from tkinter import ttk
from dockstring import load_target
import rdkit.Chem as Chem

# Install necessary libraries
subprocess.call(["pip", "install", "dockstring"])
subprocess.call(["pip", "install", "gradio"])

# Install necessary libraries
subprocess.call(["conda", "install", "-c", "conda-forge", "openbabel"])

# List of 58 protein targets
protein_targets = [
    "IGF1R", "JAK2", "KIT", "LCK", "MAPK14", "MAPKAPK2", "MET", "PTK2", "PTPN1", "SRC",
    "ABL1", "AKT1", "AKT2", "CDK2", "CSF1R", "EGFR", "KDR", "MAPK1", "FGFR1", "ROCK1",
    "MAP2K1", "PLK1", "HSD11B1", "PARP1", "PDE5A", "PTGS2", "ACHE", "MAOB", "CA2", "GBA",
    "HMGCR", "NOS1", "REN", "DHFR", "ESR1", "ESR2", "NR3C1", "PGR", "PPARA", "PPARD",
    "PPARG", "AR", "THRB", "ADAM17", "F10", "F2", "BACE1", "CASP3", "MMP13", "DPP4",
    "ADRB1", "ADRB2", "DRD2", "DRD3", "ADORA2A", "CYP2C9", "CYP3A4", "HSP90AA1"
]

def perform_docking():
    protein_target = protein_target_var.get()
    smiles = smiles_entry.get()
    target = load_target(protein_target)

    # Split the input smiles by commas
    smiles_list = smiles.split(",")

    docking_results = []
    for smiles in smiles_list:
        # Convert SMILES string to RDKit molecule object
        ligand = Chem.MolFromSmiles(smiles.strip())

        if ligand is None:
            docking_results.append(f"Invalid SMILES string '{smiles.strip()}'. Please enter a valid SMILES.")
        else:
            # Perform docking
            score, _ = target.dock(smiles.strip())
            docking_results.append(f"Docking score for SMILES '{smiles.strip()}': {score:.3f}")

    # Display the docking scores
    result_text.config(state=tk.NORMAL)
    result_text.delete('1.0', tk.END)
    result_text.insert(tk.END, "\n".join(docking_results))
    result_text.config(state=tk.DISABLED)

# Create the main tkinter window
root = tk.Tk()
root.title("Molecular Docking")

# Protein target selection
protein_target_label = ttk.Label(root, text="Select Protein Target:")
protein_target_label.grid(row=0, column=0, padx=10, pady=5, sticky="w")
protein_target_var = tk.StringVar(root)
protein_target_combobox = ttk.Combobox(root, textvariable=protein_target_var, values=protein_targets)
protein_target_combobox.grid(row=0, column=1, padx=10, pady=5)

# SMILES entry
smiles_label = ttk.Label(root, text="Enter SMILES separated by commas:")
smiles_label.grid(row=1, column=0, padx=10, pady=5, sticky="w")
smiles_entry = ttk.Entry(root)
smiles_entry.grid(row=1, column=1, padx=10, pady=5)

# Docking button
dock_button = ttk.Button(root, text="Perform Docking", command=perform_docking)
dock_button.grid(row=2, column=0, columnspan=2, padx=10, pady=10)

# Result display
result_text = tk.Text(root, height=10, width=50, state=tk.DISABLED)
result_text.grid(row=3, column=0, columnspan=2, padx=10, pady=5)

root.mainloop()
