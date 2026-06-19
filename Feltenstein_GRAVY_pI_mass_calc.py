import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# Load the Excel file (fix path issue with raw strings or forward slashes)
file_path = r'C:\Pycharm_Projects\name_of_your_document.csv'  # Corrected file path, change "name_of_your_document" to your document name
df = pd.read_csv(file_path)

# Function to calculate properties
def calculate_properties(sequence):
    analysed_seq = ProteinAnalysis(sequence)
    gravy = analysed_seq.gravy()
    iso_point = analysed_seq.isoelectric_point()
    aromaticity = analysed_seq.aromaticity()
    instability = analysed_seq.instability_index()
    return gravy, iso_point, aromaticity, instability

# Apply the function to each sequence in the DataFrame
df[['GRAVY', 'Isoelectric Point', 'Aromaticity', 'Instability']] = df['Sequence'].apply(
    lambda seq: pd.Series(calculate_properties(seq))
)

# Save the updated DataFrame to a new Excel file (ensure output file path is correct)
output_file_path = r'C:\Pycharm_Projects\name_of_your_document_finished.csv'  # Corrected file path, adding "finished" to end of document name to distinguish the two
df.to_csv(output_file_path, index=False)

print("Calculations completed and saved to the new Excel file.")

