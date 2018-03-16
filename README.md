# pachinam_obj.py
To be used from within PyMOL, right after cif files have been loaded into objects
Uses the already defined objects to pull pairs of chain identifiers and names from cif files and create separate objects - essentially does the job of split_chains but adds the defined name from the cif file.
Names are of the format:
cifID_name_chainID
Hard-codes destination of cif files to a directory which is called cif and is within the current working directory - this means you need to have the cif files in this directory which is in the working directory (since the script is being ran from within PyMOL this would be the current working directory of PyMOL)
Tries to remove some lengthy parts of names often found in big structure files (e.g. 30S ribosomal protein, ribosomal RNA etc.) for better viewing in PyMOL.
Has been tested to work up to version 2.1.0 in PyMOL
Uses biopython and python3 so have them installed.
