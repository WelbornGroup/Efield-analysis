import pandas as pd

# Define file paths
file1_path = "proj_totfield_1_latest.csv"
file2_path = "proj_totfield_2_latest.csv"
file3_path = "proj_totfield_3_latest.csv"

# Define output file path
output_path = "solvent_250.csv"


# Read each file into a separate list of DataFrames
df1  = pd.read_csv(file1_path)
df2  = pd.read_csv(file2_path)
df3  = pd.read_csv(file3_path)

# Concatenate the DataFrames vertically (axis=0)
result_df = pd.concat([df1, df2, df3], axis=1)

# Save the merged DataFrame to a new CSV file
result_df.to_csv(output_path, index=False)

