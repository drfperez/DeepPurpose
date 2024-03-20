try:
    import requests
    import pubchempy as pc
    from Bio import SeqIO
    from tkinter import tk
    from tkinter import ttk
    from tkinter import messagebox
    from io import StringIO
    from DeepPurpose import utils
    from DeepPurpose import DTI as models
except ImportError:
    print("Installing required libraries...")
    import sys
    try:
        import subprocess
        subprocess.check_call([sys.executable, "-m", "pip", "install", 'requests', 'pubchempy', 'biopython', 'git+https://github.com/bp-kelley/descriptastorus', 'pandas-flavor'])
        print("Libraries installed successfully!")
    except subprocess.CalledProcessError:
        print("Failed to install required libraries.")
        sys.exit(1)

# Continue with the rest of your code...


# Import required libraries
import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
import requests
from io import StringIO
import pubchempy as pc
from DeepPurpose import utils
from DeepPurpose import DTI as models
from Bio import SeqIO

# Load the pretrained models for different datasets
model_binding = models.model_pretrained(model='MPNN_CNN_BindingDB')
model_kiba = models.model_pretrained(model='MPNN_CNN_KIBA')
model_davis = models.model_pretrained(model='MPNN_CNN_DAVIS')

def fetch_pdb_sequence(pdb_code):
    try:
        # Fetch the PDB file from RCSB PDB website
        pdb_url = f"https://files.rcsb.org/download/{pdb_code.lower()}.pdb"
        response = requests.get(pdb_url)

        # Check if the request was successful
        if response.status_code == 200:
            # Parse the response content as a PDB file
            pdb_content = response.text
            pdb_file = StringIO(pdb_content)

            # Parse the PDB file and extract the target sequence
            target_sequence = ""
            for structure in SeqIO.parse(pdb_file, "pdb-seqres"):
                target_sequence += str(structure.seq)
            return target_sequence
        else:
            return f"Error: PDB code {pdb_code} not found."
    except Exception as e:
        return f"Error fetching sequence for PDB code {pdb_code}: {str(e)}"

def get_smiles_from_name(chemical_name):
    try:
        compound = pc.get_compounds(chemical_name, 'name')
        if compound:
            return compound[0].canonical_smiles
        else:
            return None
    except Exception as e:
        return f"Error fetching SMILES for {chemical_name}: {str(e)}"

def predict_interaction():
    dataset = dataset_var.get()
    drug_type = drug_type_var.get()
    drug_input = drug_input_entry.get()
    target_type = target_type_var.get()
    target_input = target_input_entry.get()

    if target_type == 'PDB Code':
        target = fetch_pdb_sequence(target_input)
        if "Error" in target:
            messagebox.showerror("Error", target)
            return
    else:
        target = target_input

    if drug_type == 'Chemical Names':
        drugs = [get_smiles_from_name(name.strip()) for name in drug_input.split(",")]
        drugs = [smiles for smiles in drugs if smiles]
    else:
        drugs = drug_input.split(",")

    if dataset == 'BindingDB':
        model = model_binding
    elif dataset == 'KIBA':
        model = model_kiba
    elif dataset == 'DAVIS':
        model = model_davis

    predictions = []

    for drug in drugs:
        X_pred = utils.data_process(X_drug=[drug], X_target=[target], y=[0],
                                    drug_encoding='MPNN', target_encoding='CNN',
                                    split_method='no_split')
        y_pred = model.predict(X_pred)
        predictions.append(str(y_pred[0]))

    result_label.config(text=", ".join(predictions))

root = tk.Tk()
root.title("Drug-Target Interaction Prediction")

# Dataset selection
dataset_label = ttk.Label(root, text="Select Dataset:")
dataset_label.grid(row=0, column=0, padx=5, pady=5, sticky="w")
dataset_var = tk.StringVar(value="BindingDB")
dataset_combobox = ttk.Combobox(root, textvariable=dataset_var, values=['BindingDB', 'KIBA', 'DAVIS'])
dataset_combobox.grid(row=0, column=1, padx=5, pady=5, sticky="w")

# Drug input
drug_type_label = ttk.Label(root, text="Select Drug Input Type:")
drug_type_label.grid(row=1, column=0, padx=5, pady=5, sticky="w")
drug_type_var = tk.StringVar(value="SMILES")
drug_type_combobox = ttk.Combobox(root, textvariable=drug_type_var, values=['SMILES', 'Chemical Names'])
drug_type_combobox.grid(row=1, column=1, padx=5, pady=5, sticky="w")

drug_input_label = ttk.Label(root, text="Enter Drug Input:")
drug_input_label.grid(row=2, column=0, padx=5, pady=5, sticky="w")
drug_input_entry = ttk.Entry(root)
drug_input_entry.grid(row=2, column=1, padx=5, pady=5, sticky="w")

# Target input
target_type_label = ttk.Label(root, text="Select Target Input Type:")
target_type_label.grid(row=3, column=0, padx=5, pady=5, sticky="w")
target_type_var = tk.StringVar(value="Sequence")
target_type_combobox = ttk.Combobox(root, textvariable=target_type_var, values=['Sequence', 'PDB Code'])
target_type_combobox.grid(row=3, column=1, padx=5, pady=5, sticky="w")

target_input_label = ttk.Label(root, text="Enter Target Input:")
target_input_label.grid(row=4, column=0, padx=5, pady=5, sticky="w")
target_input_entry = ttk.Entry(root)
target_input_entry.grid(row=4, column=1, padx=5, pady=5, sticky="w")

# Prediction button
predict_button = ttk.Button(root, text="Predict", command=predict_interaction)
predict_button.grid(row=5, column=0, columnspan=2, padx=5, pady=5)

# Result
result_label = ttk.Label(root, text="")
result_label.grid(row=6, column=0, columnspan=2, padx=5, pady=5)

root.mainloop()
