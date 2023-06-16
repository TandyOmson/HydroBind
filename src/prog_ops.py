# Fucntions involved in program control
# Testing and error handling occurs betwen modules
# Restart also occurs betwen each module

from shutil import rmtree
import numpy as np
import os

def delete_and_create_dir(dirname):
    try:
        os.mkdir(f"{dirname}")
    except FileExistsError:
        rmtree(f"{dirname}")
        os.mkdir(f"{dirname}")

def coords_check(coords_array,outfile):
    """ Takes coordinates of a molecule and makes a quick xyz outfile"""

    coords_array = np.array(coords_array)

    with open(outfile,'w') as f:
        f.write(f"{len(coords_array)}\n\n")
        for atom in coords_array:
            f.write(f"X {atom[0]} {atom[1]} {atom[2]}\n")
    
    return "Created xyz file for molecule", outfile