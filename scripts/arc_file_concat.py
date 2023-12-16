import os

filenames = ["dynamic_water_move_1.arc", "dynamic_water_move_2.arc", "dynamic_water_move_3.arc"]
output_file = "big_solvent_250_frames.arc"

with open(output_file, "w") as outfile:
    for filename in filenames:
        with open(filename, "r") as infile:
            outfile.write(infile.read())

print(f"Combined data written to: {output_file}")