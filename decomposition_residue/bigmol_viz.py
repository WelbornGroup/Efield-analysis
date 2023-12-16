import numpy as np
import pandas as pd

df = pd.read_csv('big_proj_250.csv')
file_path = 'big_mol_50_formatted.arc'
log_file_path = 'resid.log'
global global_insert_atom_df

# Converts the CSV file into a hashmap with the thresholded EF values
def filter_df(df):
  
    low, high = -3, 3
    frame_count = len(df.columns) - 1

    # Using vectorized operations to filter values based on thresholds
    filtered_columns = df.iloc[:, 1:].apply(lambda col: col.where((col < low) | (col > high), 0))

    outer_field_map = {}
    highest_keys_df = {}
    highest_keys_count = 0
    pos_count, neg_count = 0,0
    max_pos_cnt, max_neg_cnt = 0, 0

    for k, column_name in enumerate(filtered_columns.columns[0:], start=0):
        # print('Column name and k:', column_name, k)
        field_val = np.array(filtered_columns[column_name])
        field_val = field_val.astype(np.float64)
        pos_count = (np.count_nonzero(field_val[0:1149] > 3))
        neg_count = (np.count_nonzero(field_val[0:1149] < -3))
        max_pos_cnt = max_pos_cnt if pos_count < max_pos_cnt else pos_count
        max_neg_cnt = max_neg_cnt if neg_count < max_neg_cnt else neg_count
        field_map = {i + 1: val for i, val in enumerate(field_val) if val != 0}
        
        keys_number = len(field_map)
        
        if keys_number > highest_keys_count:
            highest_keys_count = keys_number
            highest_keys_df = field_map
        
        outer_field_map[k+1] = field_map

    # print(outer_field_map)
    max_pos_atoms, max_neg_atoms = max_pos_cnt * 18, max_neg_cnt * 18
    max_atom_dict = {
        'max_pos_atoms' : max_pos_atoms,
        'max_neg_atoms' : max_neg_atoms 
    }

    return outer_field_map, frame_count, highest_keys_count, highest_keys_df, max_atom_dict

# Converts the resid.log file into a hashmap
def create_log_map(log_path):
    log_file = open(log_path, 'r')
    lines = log_file.readlines()
    log_map = {}
    residue_map = {}
    for i in range(len(lines)):
        # if i>20:
        #     break
        line_arr = lines[i].strip().split(' ')
        filtered_line = list(filter(None, line_arr))

        if len(filtered_line) == 3:
            log_map[int(filtered_line[1])] = (filtered_line[2])
            residue_map[int(filtered_line[0])] = int(filtered_line[1])
    
    print("create_log_map() done")
    return log_map, residue_map

# Builds the insertion dataframe with dummy values 
def build_insertion_df(max_atom_dict, header_0):

    # print(max_atom_dict)
    global global_insert_atom_df
    new_mol_count = int(header_0) + 1
    local_radium_index_and_df_start = [0,0]
    local_selenium_index_and_df_start = [0,0]
    total_atom_count = 0
    insertion_cnt = 0
    for key in max_atom_dict.keys():
        # After Python 3.7+ dictionaries are iterated in the order in which they are inserted
        if key == 'max_pos_atoms':
            local_radium_index_and_df_start[0] = new_mol_count 
            local_radium_index_and_df_start[1] = insertion_cnt
        else:
            local_selenium_index_and_df_start[0] = new_mol_count
            local_selenium_index_and_df_start[1]  = insertion_cnt
        for i in range(max_atom_dict[key]):
            temp_dict = {
                header_0: [new_mol_count],
                'atom_type': 'Ra' if key == 'max_pos_atoms' else 'Se',
                'x': -60  if key == 'max_pos_atoms' else 60,   # These x,y,z coords are at the box edge
                'y': 70  if key == 'max_pos_atoms' else -70,
                'z': 70  if key == 'max_pos_atoms' else -70,
                'atom_number': '950'  if key == 'max_pos_atoms' else '800',
                'link1': '',
                'link2': '',
                'link3': '',
                'link4': '',
            }
            df = pd.DataFrame(temp_dict)
            global_insert_atom_df = pd.concat([global_insert_atom_df, df], ignore_index=True)
            new_mol_count += 1
            insertion_cnt += 1


    # print("Radium Index start : ", local_radium_index_and_df_start[0])
    # print("Selenium Index start : ", local_selenium_index_and_df_start[0])
    # print("Radium Index DF start : ", local_radium_index_and_df_start[1])
    # print("Selenium Index DF start : ", local_selenium_index_and_df_start[1])
    total_atom_count = local_radium_index_and_df_start[0] + insertion_cnt - 1
    # print("final count: ", total_atom_count)
   
    return local_radium_index_and_df_start, local_selenium_index_and_df_start, total_atom_count
    

# Replaces dummy atom rows in global_insert_atom_df with actual atom values
def insert_atom(row_number, residue_number, field_map, ra_index_and_df_start, se_index_and_df_start):
    
    global global_insert_atom_df
  

    atom_type =  'Ra' if field_map[residue_number] > 3 else 'Se'
    atom_number = '950' if field_map[residue_number] > 3 else '800'
   
    # print(global_insert_atom_df)
    
    insertion_index = ra_index_and_df_start[1] if field_map[residue_number] > 3 else se_index_and_df_start[1]
    mol_index = ra_index_and_df_start[0] if field_map[residue_number] > 3 else se_index_and_df_start[0]
    header_list = [header_0 ,'atom_type', 'x', 'y', 'z', 'atom_number']
    value_list = [mol_index, atom_type, row_number.x, row_number.y, row_number.z, atom_number]
    # print("insertion index : ", insertion_index)
    global_insert_atom_df.loc[insertion_index, header_list] = value_list

# Adds the extra atoms to the original arc file and writes a new file
def write_arc(arc_df, field_map, log_map, output_file, arcfile_number, ra_index_and_df_start, se_index_and_df_start, total_atom_count):
    # print("write_arc() called")

    residue_number = 0
    arc_line_cnt = -1  
    
    for row in arc_df.itertuples():
        arc_line_cnt += 1
        
        if arc_line_cnt in log_map:
            residue_number += 1

        if residue_number < 1146 and residue_number in field_map:
            
            insert_atom(row, residue_number, field_map, ra_index_and_df_start, se_index_and_df_start)
            if field_map[residue_number] > 3:    
                ra_index_and_df_start[0] += 1
                ra_index_and_df_start[1] += 1
            else:
                se_index_and_df_start[0] += 1 
                se_index_and_df_start[1] += 1
            
    print("final count: ", total_atom_count)

    arc_df.rename(columns={header_0: total_atom_count}, inplace=True)
    arc_df = arc_df.astype({total_atom_count: int})

    global global_insert_atom_df
    
    global_insert_atom_df.rename(columns={header_0: total_atom_count}, inplace=True)
    
    arc_df = pd.concat([arc_df, global_insert_atom_df], ignore_index=True)

    # Filtering out global variable header
    arc_df = arc_df.dropna(subset=['atom_type'])

    # Resetting global variable
    global_insert_atom_df = pd.DataFrame(columns=[header_0, 'atom_type', 'x','y','z','atom_number','link1','link2','link3','link4'])
    global_insert_atom_df.reset_index(drop=True)
    
    # Write to arc file
    output_file.write(str(arc_df.columns[0]) + '\n')
    arc_df.to_csv(output_file, sep='\t', index=False, header=False)


# Splits the DataFrame into smaller chunks
def split_dataframe(df, chunk_size): 
    chunks = list()
    num_chunks = (len(df) // chunk_size + 1)
    print("no. of chunks :", num_chunks)
    for i in range(num_chunks):
        chunk = df[i*chunk_size:(i+1)*chunk_size]  
        # Resetting the index
        chunk.reset_index(inplace=True, drop=True)
        chunks.append(chunk)
    return chunks


# ********* CODE ENTRY **********************************

outer_field_map, frame_count, max_keys_count, max_keys_df, max_atom_dict = filter_df(df)

# Loading the text file into a Pandas DataFrame
arc_df_full = pd.read_csv(file_path)

new_column_names = [arc_df_full.columns[0]] + ['atom_type', 'x','y','z','atom_number','link1','link2','link3','link4'] 
arc_df_full.columns = new_column_names
header_0 = new_column_names[0]

global_insert_atom_df = pd.DataFrame(columns=[header_0, 'atom_type', 'x','y','z','atom_number','link1','link2','link3','link4'])

log_map, residue_map = create_log_map(log_file_path)

output_file = open('big_mol_output_250_frames.arc', 'w')
data_chunks = split_dataframe(arc_df_full, int(header_0)+2)

for i in range(50):  
    print(f"{i+1}th iteration : ")
    ra_index_and_df_start, se_index_and_df_start, total_atom_count = build_insertion_df(max_atom_dict, header_0)
    write_arc(data_chunks[i], outer_field_map[i + 1], log_map, output_file, i,ra_index_and_df_start, se_index_and_df_start, total_atom_count)
    
