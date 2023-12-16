# Efield-analysis

The visualization scripts are run on arc files to display the areas in the residues/molecules with significant electrical activity. The scripts make use of the associated electric fields CSV file and the arc file to create a new arc file with added atoms in order to visualize the electric fields. 


Steps for running the visualization scripts:

# Step 1:

Sometimes the CSV files were too big and they were received in multiple files. In such cases, the 'append_csv.py' script was used to join all the parts. In a similar scenario, the 'arc_file_concat.py' script was used to join all the arc file parts.

# Step 2:

The final version of the visualization script uses Pandas Dataframes for the required operations. Since the arc files are essentially text files, a script 'arc_formatter.py' was used to first convert an arc file to a comma separated value file which can be directly read as a Dataframe.

# Step 3:

The command to run the visualization scripts is pretty straightforward, just navigate to the respective script directory and run it as a python command.

-> python bigmol_viz.py







