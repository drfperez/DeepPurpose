
import requests
from google.colab import files
import gradio as gr

# Function to download and save PDB or FASTA file
def download_files(pdb_codes, file_type):
    codes = pdb_codes.split(',')
    paths = []
    
    for code in codes:
        code = code.strip().upper()
        if file_type == 'PDB':
            url = f'https://files.rcsb.org/download/{code}.pdb'
        elif file_type == 'FASTA':
            url = f'https://www.rcsb.org/fasta/entry/{code}'
        
        response = requests.get(url)
        
        if response.status_code == 200:
            file_extension = 'pdb' if file_type == 'PDB' else 'fasta'
            filename = f'{code}.{file_extension}'
            with open(filename, 'wb') as file:
                file.write(response.content)
            paths.append(filename)
        else:
            print(f'Error: Could not download the file for {code} in {file_type} format.')
    
    return paths

# Gradio interface
iface = gr.Interface(
    fn=download_files,
    inputs=[
        gr.Textbox(placeholder="Enter PDB codes separated by commas"),
        gr.Dropdown(choices=['PDB', 'FASTA'], label="Select File Type")
    ],
    outputs=gr.File(label="Download Files"),
    title="PDB/FASTA File Downloader",
    description="Enter PDB codes separated by commas and select the file type to download corresponding files."
)

# Run the Gradio app
iface.launch()