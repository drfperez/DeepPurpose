!pip install gradio pubchempy rdkit
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import AllChem
import gradio as gr

def iupac_to_mol2(iupac_name):
    # Search for the compound by its IUPAC name
    compounds = pcp.get_compounds(iupac_name, 'name')
    if len(compounds) == 0:
        return "No compound found for the given IUPAC name."

    # Get the compound with the highest CID (most relevant result)
    compound = max(compounds, key=lambda x: x.cid)

    # Get the 2D structure (SMILES) and convert it to an RDKit molecule
    smiles = compound.isomeric_smiles
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "Failed to generate a molecular structure from the given IUPAC name."

    # Add hydrogens, generate 3D coordinates, and optimize the structure
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    AllChem.MMFFOptimizeMolecule(mol)

    # Write the optimized structure to a mol2 file
    output_file = "output.mol2"
    Chem.MolToMolFile(mol, output_file)
    return output_file

iface = gr.Interface(fn=iupac_to_mol2, 
                     inputs="text", 
                     outputs="file", 
                     title="IUPAC Name to Mol2 Converter",
                     description="Enter the IUPAC name of the molecule and download the generated Mol2 file.")
iface.launch(inline=False)