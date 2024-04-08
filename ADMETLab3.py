import tkinter as tk
from tkinter import messagebox
import requests
import json

baseUrl = 'https://admetlab3.scbdd.com'

def save_as_json(json_response):
    with open('result.json', 'w') as f:
        json.dump(json_response, f)
    print("JSON file saved successfully!")

def process_smiles():
    smiles_input = entry.get("1.0", tk.END).strip()  # Get text from the text widget
    smiles_list = smiles_input.split(',')

    param = {
        'SMILES': smiles_list
    }

    response = requests.post(baseUrl + '/api/admet', json=param)

    if response.status_code == 200:
        json_response = response.json()
        print(json_response)
        save_as_json(json_response)
        messagebox.showinfo("Success", "JSON response saved successfully!")
    else:
        messagebox.showerror("Error", "Failed to retrieve data from the API.")

# Create tkinter window
window = tk.Tk()
window.title("ADMETLab 3.0")
window.geometry('600x400')  # Set width and height

# Explanation text
explanation = "Enter the SMILES strings separated by commas in the input box below.\n"\
              "After entering, click the 'Submit' button to fetch ADMET data from the API.\n"\
              "The JSON response will be saved as 'result.json'."

explanation_label = tk.Label(window, text=explanation, wraplength=550, justify=tk.LEFT)
explanation_label.pack(pady=10)

# Create label and text widget for SMILES input
smiles_label = tk.Label(window, text="SMILES:")
smiles_label.pack()

entry = tk.Text(window, height=10, width=50)  # Adjust the height of the text widget
entry.pack()

# Create button to submit SMILES
button = tk.Button(window, text="Submit", command=process_smiles)
button.pack(pady=10)

# Run the tkinter event loop
window.mainloop()
