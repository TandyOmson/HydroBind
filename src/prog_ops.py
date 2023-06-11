# Fucntions involved in program control
# Testing and error handling occurs betwen modules
# Restart also occurs betwen each module

from shutil import rmtree
import os

def delete_and_create_dir(dirname):
    try:
        os.mkdir(f"D{dirname}")
    except FileExistsError:
        rmtree(f"D{dirname}")
        os.mkdir(f"D{dirname}")