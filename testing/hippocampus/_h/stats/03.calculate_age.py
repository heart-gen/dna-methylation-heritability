import pandas as pd

# Load the dataset (tab-delimited)
df = pd.read_csv("/projects/p32505/users/elisa/testing/hippocampus/_m/filtered_file.txt", sep="\t")

# Calculate mean and standard deviation for 'agedeath'
mean_age = df["agedeath"].mean()
std_age = df["agedeath"].std()

print(f"Mean age at death: {mean_age:.2f}")
print(f"Standard deviation of age at death: {std_age:.2f}")

# Optionally, save the results to a file
results = pd.DataFrame({
    "Statistic": ["Mean", "Standard Deviation"],
    "Age": [mean_age, std_age]
})
results.to_csv("age_stats.txt", sep="\t", index=False)
