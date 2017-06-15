##################################################################
#                   THE CONTEXT
##################################################################
The goal of this program is to study the sRNP H/ACA complex

This complex catalyzes isomerisation of uridin into pseudo-uridin, which is important for RNA structure, particularly tRNAs.

The project was to study global and local changes using RMSD and frequency of interaction of AAs with another chain.



##################################################################
#                   THE PROGRAMM
##################################################################
This program aims to processes a 3D simulation of a given protein.

More precisely, it calculates the RMSD of the protein and of its chains.
It also calculate interaction frequency between pairs of Amino Acids
In the simulation, coordinates of reference PDB and simulation PDB have been adjusted.

2 programs have been coded : a "global" and a "local". 
They have a modular structure and both call the MDtools.py module, which must be in the same directory



##################################################################
#                      USAGE
##################################################################

MDglobal.py -pdbref -pdb -chaine [-seuil] [-dist] [-RMSD] [-Matrice] [-ofreq] [-oRMSD]


mandatory:
	===========
     
	-pdbref	-> a reference pdb file

	-pdb	-> the simulation pdb file
	
	-chaine	-> reference chain

	optionnal:
	=========
	-seuil   -> threshold under which 2 AAs are reported to interact, in Angstrom, defaults to 4. Indicate desired threshold
     
	-dist    -> the way distance between AAs is calculated :"bar" for AA's mass center, "mini" for minimal distance between atoms of the AAs
	
	-RMSD	 -> the way RMSD is calculated: "all" for each atom, "alpha" for alpha carbons, "backbone" for the AA's backbone (CA, C, O, N).

	-Matrice -> if "yes", returns matrix of RMSDs; defaults to "no"

	-ofreq	-> name of the interaction frequencies file (par defaut freq_interface.txt), WITH extension
	
	-oRMSD	-> name of the  RMSD file (defaults to RMSD), WITHOUT extension. txt format

	"""
