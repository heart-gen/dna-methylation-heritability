import pandas as pd

# Load the main dataset (tab-delimited)
df = pd.read_csv("/projects/p32505/users/elisa/testing/dlpfc/_m/filtered_file.txt", sep="\t")

# Count values in the 'sex' column
sex_counts = df["sex"].value_counts()

# Print results
print("Number of Females:", sex_counts.get("F", 0))
print("Number of Males:", sex_counts.get("M", 0))

# If you want to save the counts to a file:
sex_counts.to_csv("sex_counts.txt", sep="\t")
