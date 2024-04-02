import tkinter as tk
from tkinter import ttk, messagebox
from Bio import SeqIO
import requests
import pubchempy as pc
from io import BytesIO
from DeepPurpose import utils
from DeepPurpose import DTI as models

# Load the pretrained models for different datasets
model_binding = models.model_pretrained(model='MPNN_CNN_BindingDB')
model_kiba = models.model_pretrained(model='MPNN_CNN_KIBA')
model_davis = models.model_pretrained(model='MPNN_CNN_DAVIS')

def get_smiles_from_name(chemical_name):
    try:
        compound = pc.get_compounds(chemical_name, 'name')
        if compound:
            return compound[0].canonical_smiles
        else:
            return None
    except Exception as e:
        return f"Error fetching SMILES for {chemical_name}: {str(e)}"

def DTI_pred(data, drug_type, drug_input, target_input):
    try:
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
            X_pred = utils.data_process(X_drug=[drug], X_target=[target_input], y=[0],
                                        drug_encoding='MPNN', target_encoding='CNN',
                                        split_method='no_split')

            y_pred = model.predict(X_pred)

            predictions.append(str(y_pred[0]))

        return ", ".join(predictions)
    except Exception as e:
        return f"Error predicting drug-target interaction: {str(e)}"

def predict():
    dataset = dataset_var.get()
    drug_type = drug_type_var.get()
    drug_input = drug_input_entry.get()
    target_input = sequence_text.get(1.0, tk.END)

    result = DTI_pred(dataset, drug_type, drug_input, target_input)
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
sequence_label = ttk.Label(root, text="FASTA sequence:")
sequence_label.grid(row=3, column=0, padx=5, pady=5, sticky="w")
sequence_text = tk.Text(root, height=5, width=50)
sequence_text.grid(row=3, column=1, padx=5, pady=5, sticky="w")

# Prediction button
predict_button = ttk.Button(root, text="Predict", command=predict)
predict_button.grid(row=4, column=0, columnspan=2, padx=5, pady=5)

# Result
result_label = ttk.Label(root, text="")
result_label.grid(row=5, column=0, columnspan=2, padx=5, pady=5)

root.mainloop()
