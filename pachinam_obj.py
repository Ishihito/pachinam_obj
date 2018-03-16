#!/usr/bin/env python3
import sys
from sys import argv
import re
import copy
from Bio.PDB import *
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

"""
To be used from within PyMOL script, right after cif files have been loaded into objects 
and alignment has been done
Uses the already defined objects to further define their chains into separate objects
with names pulled from the cif files
Hard-codes the destination of cif files.
"""

parser = MMCIFParser(QUIET=True)

string = "../cif/"
fileext = ".cif"
cifs = list(cmd.get_object_list())
newcifs = [x + fileext for x in cifs]
newercifs = [string + x for x in newcifs]
filenames = newercifs


#Creates the objects with given dictionary of chainids -> names and name of the object parent
def create_objects(parent,element_dict):
	for chain in sorted(element_dict.keys()):
		x = re.sub('30S ribosomal protein ', '',element_dict[chain], flags=re.I)
		x = re.sub('50S ribosomal protein ', '',x, flags=re.I)
		x = re.sub('60S ribosomal protein ', '',x, flags=re.I)
		x = re.sub('40S ribosomal protein ', '',x, flags=re.I)
		x = re.sub('DNA-directed ', '',x, flags=re.I)
		x = re.sub('ribosomal protein ', '',x, flags=re.I)
		x = re.sub('ribosomal RNA', 'rRNA',x, flags=re.I)
		x = re.sub('translation', '',x, flags=re.I)
		x = re.sub(' initiation ', '',x, flags=re.I)
		x = re.sub('RNA POLYMERASES', 'RNApol',x, flags=re.I)
		x = re.sub('RNA POLYMERASE', 'RNApol',x, flags=re.I)
		x = re.sub('SUBUNIT', 'subunit',x, flags=re.I)
		cmd.create(parent + "_" + x + "_" + chain, "chain " + chain + " and " + object)

#Creates a dictionary {chainids -> names} from a given cif file
def assign(filename, dict):
	n=0
	#Grab header
	mmcif_dict = MMCIF2Dict(filename)
	#Grab entity names
	details = mmcif_dict['_entity.pdbx_description']
	#Grab chain ids
	strand_id = mmcif_dict['_entity_poly.pdbx_strand_id'] 
	for x in strand_id:
		dict[x]=details[n]
		n+=1
	return dict
	

for file in filenames:
	temp_chain_names={}
	assign(file, temp_chain_names)
	chain_names=copy.deepcopy(temp_chain_names)
	#Following for loop is for structures which have multiple copies of the structure
	#otherwise the hashes would get overwritten due to the repeating of objects
	for chid in temp_chain_names.keys():
		if re.search(",", chid):
			temp = chid.split(",")
			for x in temp:
				chain_names[x]=chain_names[chid]
			del chain_names[chid]
	cutfile = file.split("/")[-1]
	cutfile = re.sub('\.cif', '', cutfile)
	create_objects(cutfile, chain_names)
