import tkinter as tk
from tkinter import messagebox
import requests
import json
import pandas as pd
import csv
import os  # Import os module for file handling


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

    try:
      # Define headers based on JSON structure (modify as needed)
      data = json_response['data']['data']  # Access data list
      headers = list(data[0].keys())  # Extract headers from the first element in the data list

      # Extract data and populate rows
      data_rows = []
      for item in data:
        row = []
        for key in headers:
          row.append(item.get(key))  # Access data using key
        data_rows.append(row)

      # Write data to CSV
      output_filename = "admet_data.csv"

      # Option 1: Clear Input Field (prevents multiple calls)
      # entry.delete(1.0, tk.END)  # Clear input field after successful processing

      # Option 2: Check for Existing File (Prompt User)
      if os.path.exists(output_filename):
        answer = messagebox.askquestion("File Exists", f"A file named '{output_filename}' already exists. Overwrite?")
        if answer.lower() != "yes":
          return  # Exit function if user doesn't want to overwrite

      # Option 3: Generate Unique Name (commented out as unnecessary here)
      # timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
      # output_filename = f"admet_data_{timestamp}.csv"

      with open(output_filename, "w") as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(headers)
        csv_writer.writerows(data_rows)
      print(f"CSV file saved successfully: {output_filename}")
    except Exception as e:
      print(f"Error converting JSON to CSV: {e}")

  else:
    messagebox.showerror("Error", "Failed to retrieve data from the API.")


# Create tkinter window
window = tk.Tk()
window.title("ADMETLab 3.0")
window.geometry('600x400')  # Set width and height

# Explanation text
explanation = "Enter the SMILES strings separated by commas in the input box below.\n"\
                "After entering, click the 'Submit' button to fetch ADMET data from the API.\n"\
                "The JSON response will be saved as 'result.json' and a CSV will be generated with data from nested structures in separate columns."

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
