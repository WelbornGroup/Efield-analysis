
input_file_path = 'big_mol_files/dynamic_input_250.arc'
output_file_path = 'big_mol_files/big_mol_250_formatted.arc'
output_file = open(output_file_path, 'w')

# Read the input file and process each line
with open(input_file_path, 'r') as input_file:
    lines = input_file.readlines()

# Remove leading and trailing whitespaces from each line
lines = [line.strip() for line in lines]

# Populate the DataFrame with the formatted data
for line in lines:
    line_arr = line.split(' ')
    filtered_line = list(filter(None, line_arr))
    to_fill = 10 - len(filtered_line)
    new_line = ','.join(filtered_line) + ','*to_fill
    output_file.write(new_line + '\n')
    # print(new_line)

print(f"Formatted content written to {output_file}")

