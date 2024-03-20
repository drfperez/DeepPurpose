import tkinter as tk
from tkinter import ttk, messagebox, filedialog
from io import StringIO
from Bio import SeqIO
import requests
import pubchempy as pc
from DeepPurpose import utils
from DeepPurpose import DTI as models

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
            pdb_content = response.content
            pdb_file = StringIO(pdb_content.decode('utf-8'))

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

def parse_pdb_file(file_path):
    try:
        with open(file_path, 'r') as file:
            pdb_content = file.read()
        pdb_file = StringIO(pdb_content)
        target_sequence = ""
        for structure in SeqIO.parse(pdb_file, "pdb-seqres"):
            target_sequence += str(structure.seq)
        return target_sequence
    except Exception as e:
        return f"Error parsing PDB file: {str(e)}"

def DTI_pred(data, drug_type, drug_input, target_type, target_input):
    try:
        if target_type == 'Execute PDB Code':
            target = fetch_pdb_sequence(target_input)
        elif target_type == 'Upload PDB File':
            if not target_input:
                return "Error: No file uploaded."
            pdb_file = parse_pdb_file(target_input)
            if isinstance(pdb_file, str):
                return pdb_file
            target_sequence = ""
            for structure in SeqIO.parse(pdb_file, "pdb-seqres"):
                target_sequence += str(structure.seq)
            target = target_sequence
        elif target_type == 'Sequence':
            target = target_input

        if drug_type == 'SMILES':
            drugs = drug_input.split(",")
        elif drug_type == 'Chemical Names':
            drugs = [get_smiles_from_name(name.strip()) for name in drug_input.split(",")]
            drugs = [smiles for smiles in drugs if smiles]

        if data == 'BindingDB':
            model = model_binding
        elif data == 'KIBA':
            model = model_kiba
        elif data == 'DAVIS':
            model = model_davis

        predictions = []

        for drug in drugs:
            X_pred = utils.data_process(X_drug=[drug], X_target=[target], y=[0],
                                        drug_encoding='MPNN', target_encoding='CNN',
                                        split_method='no_split')

            y_pred = model.predict(X_pred)

            predictions.append(str(y_pred[0]))

        return ", ".join(predictions)
    except Exception as e:
        return f"Error predicting drug-target interaction: {str(e)}"

def upload_pdb():
    file_path = filedialog.askopenfilename(filetypes=[("PDB files", "*.pdb")])
    if file_path:
        pdb_sequence = parse_pdb_file(file_path)
        if isinstance(pdb_sequence, str):
            messagebox.showerror("Error", pdb_sequence)
        else:
            sequence_text.delete(1.0, tk.END)
            sequence_text.insert(tk.END, pdb_sequence)

def predict():
    dataset = dataset_var.get()
    drug_type = drug_type_var.get()
    drug_input = drug_input_entry.get()
    target_type = target_type_var.get()
    target_input = sequence_text.get(1.0, tk.END)

    result = DTI_pred(dataset, drug_type, drug_input, target_type, target_input)
    result_label.config(text=result)

# Create main window
root = tk.Tk()
root.title("Drug-Target Interaction Prediction")

# Dataset selection
dataset_label = ttk.Label(root, text="Training Dataset:")
dataset_label.grid(row=0, column=0, padx=5, pady=5, sticky="w")
dataset_var = tk.StringVar(value="BindingDB")
dataset_combobox = ttk.Combobox(root, textvariable=dataset_var, values=['BindingDB', 'KIBA', 'DAVIS'])
dataset_combobox.grid(row=0, column=1, padx=5, pady=5, sticky="w")

# Drug input
drug_type_label = ttk.Label(root, text="Drug Input Type:")
drug_type_label.grid(row=1, column=0, padx=5, pady=5, sticky="w")
drug_type_var = tk.StringVar(value="SMILES")
drug_type_combobox = ttk.Combobox(root, textvariable=drug_type_var, values=['SMILES', 'Chemical Names'])
drug_type_combobox.grid(row=1, column=1, padx=5, pady=5, sticky="w")

drug_input_label = ttk.Label(root, text="Enter Drug Input:")
drug_input_label.grid(row=2, column=0, padx=5, pady=5, sticky="w")
drug_input_entry = ttk.Entry(root)
drug_input_entry.grid(row=2, column=1, padx=5, pady=5, sticky="w")

# Target input
target_type_label = ttk.Label(root, text="Target Input Type:")
target_type_label.grid(row=3, column=0, padx=5, pady=5, sticky="w")
target_type_var = tk.StringVar(value="Execute PDB Code")
target_type_combobox = ttk.Combobox(root, textvariable=target_type_var, values=['Execute PDB Code', 'Upload PDB File', 'Sequence'])
target_type_combobox.grid(row=3, column=1, padx=5, pady=5, sticky="w")

upload_button = ttk.Button(root, text="Upload PDB", command=upload_pdb)
upload_button.grid(row=4, column=0, columnspan=2, padx=5, pady=5)

sequence_label = ttk.Label(root, text="Target Sequence:")
sequence_label.grid(row=5, column=0, padx=5, pady=5, sticky="w")
sequence_text = tk.Text(root, height=5, width=50)
sequence_text.grid(row=5, column=1, padx=5, pady=5, sticky="w")

# Prediction button
predict_button = ttk.Button(root, text="Predict", command=predict)
predict_button.grid(row=6, column=0, columnspan=2, padx=5, pady=5)

# Result
result_label = ttk.Label(root, text="")
result_label.grid(row=7, column=0, columnspan=2, padx=5, pady=5)

root.mainloop()

