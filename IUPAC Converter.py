!pip install gradio rdkit pubchempy
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import AllChem
import gradio as gr

def convert_iupac(iupac_name, output_format):
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

    # Write the optimized structure to a file based on the selected output format
    if output_format == "PDB":
        output_file = "output.pdb"
        Chem.MolToPDBFile(mol, output_file)
    elif output_format == "Mol2":
        output_file = "output.mol2"
        Chem.MolToMolFile(mol, output_file)
    elif output_format == "SDF":
        output_file = "output.sdf"
        w = Chem.SDWriter(output_file)
        w.write(mol)
        w.close()
    else:
        return "Invalid output format selected."

    return output_file

# Set up the Gradio interface
iface = gr.Interface(
    fn=convert_iupac,
    inputs=[
        gr.Textbox(label="Enter IUPAC name"),
        gr.Radio(choices=["PDB", "Mol2", "SDF"], label="Select Output Format")
    ],
    outputs=gr.File(label="Download", type="filepath"),
    title="IUPAC Name Converter",
    description="Enter the IUPAC name of the molecule and select the output format.",
    examples=[["caffeine", "PDB"], ["aspirin", "Mol2"], ["ginkgolide B", "SDF"]]
)

# Launch the Gradio app
iface.launch(inline=False)