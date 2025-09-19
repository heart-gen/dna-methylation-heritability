import pandas as pd

# Load the dataset (tab-delimited)
df = pd.read_csv("/projects/p32505/users/elisa/testing/caudate/_m/filtered_file.txt", sep="\t")

# Count occurrences of each diagnosis
diagnosis_counts = df["primarydx"].value_counts()

# Calculate percentages
diagnosis_percentages = diagnosis_counts / diagnosis_counts.sum() * 100

# Print results
for dx, pct in diagnosis_percentages.items():
    print(f"{dx}: {pct:.2f}%")

# Optionally, save results to a file
diagnosis_percentages.to_csv("diagnosis_percentages.txt", sep="\t", header=["Percentage"])
