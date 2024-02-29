import requests
from Bio import SeqIO
from io import StringIO

def fetch_pdb_sequences(pdb_code):
    try:
        # Fetch the PDB file from RCSB PDB website
        pdb_url = f"https://files.rcsb.org/download/{pdb_code.lower()}.pdb"
        response = requests.get(pdb_url)

        # Check if the request was successful
        if response.status_code == 200:
            # Parse the response content as a PDB file
            pdb_content = response.text
            pdb_file = StringIO(pdb_content)

            # Parse each record in the PDB file and print its sequence
            for structure in SeqIO.parse(pdb_file, "pdb-seqres"):
                print(f"Chain {structure.annotations['chain']}:")
                print(structure.seq)
        else:
            print(f"Error: PDB code {pdb_code} not found.")
    except Exception as e:
        print(f"Error fetching sequences for PDB code {pdb_code}: {str(e)}")

# Prompt the user to input the PDB code
pdb_code = input("Enter PDB Code: ").strip()

# Call the function to fetch and print the sequences
fetch_pdb_sequences(pdb_code)