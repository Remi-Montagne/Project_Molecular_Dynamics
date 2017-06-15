#!/usr/bin/python
# -*- coding: utf-8 -*-

###==========================================================###
###			DESCRIPTION
###==========================================================###

# This program processes a 3D simulation of a given protein. It calculates the RMSD 
# of the global protein and of its individual chains. It also calculates
# the interaction frequency between each Amino Acid (AA) and the AAs of a reference chain
# If the user asks (-Matrice yes), it returns the RMSD matrix for heatmap generation
# Input : reference and simulation pdb, proteins must have the same origin coordinates. Name of the reference chain
# Output : a txt file with RMSDs, a txt file for interaction frequencies

###==========================================================###

import sys, MDtools

###==========================================================###
###			USAGE
###==========================================================###

def usage():
	print """
	This program processes a 3D simulation of a given protein. It calculates the RMSD 
	of the global protein and of its individual chains. It also calculates
	the interaction frequency between each Amino Acid (AA) and the AAs of a reference chain
	If the user asks (-Matrice yes), it returns the RMSD matrix for heatmap generation
	
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

###==========================================================###
###			       GET ARGUMENTS
###==========================================================###

#dPDB_ref
try:
	infile_ref = sys.argv[sys.argv.index("-pdbref")+1]
	print "reference pdb", infile_ref
except:    
	usage()
	print "ERROR: Please enter the reference pdb."
	sys.exit()

#dPDB    
try:
	infile = sys.argv[sys.argv.index("-pdb")+1]
	print "simulation pdb", infile
except:    
	usage()
	print "ERROR: Please enter the simulation pdb."
	sys.exit()

#chaine	
try:
	chaine = str(sys.argv[sys.argv.index("-chaine")+1])
except:
	usage()
	print "ERROR: Please indicate the chain to compare."
	sys.exit()

#seuil
try:
	seuil = float(sys.argv[sys.argv.index("-seuil")+1])
except:
	seuil = 4.0
	
#dist how to calculate distance between AAs : mass center or min distance
try:
	dist_mode = sys.argv[sys.argv.index("-dist")+1]
except:
	dist_mode = "bar"

#RMSD how to calculate the RMSD : every atoms (all), alpha C only (alpha) oor backbone's AAs (backbone)
try:
	mode_RMSD = sys.argv[sys.argv.index("-RMSD")+1]
except:
	mode_RMSD = "all"
	
#RMSD Matrix
try:
	matrice = sys.argv[sys.argv.index("-Matrice")+1]
except:
	matrice = "no"

#nom du fichier de sortie des frequences
try:
	outfile_freq = sys.argv[sys.argv.index("-ofreq")+1]
except:
	outfile_freq = "freq_interface.txt"

#début du nom des fichier RMSD sans
try:
	outfile_RMSD = sys.argv[sys.argv.index("-oRMSD")+1]
except:
	outfile_RMSD = "RMSD"


###==========================================================###
###			       MAIN
###==========================================================###

#get coordinates of the atoms in dictionaries
dPDB_ref = MDtools.load_PDB(infile_ref)
dPDB = MDtools.load_PDB(infile)

#calculate RMSDs and RMSD matrix if specified
MDtools.RMSD(dPDB_ref,dPDB, mode_RMSD, outfile_RMSD)
if matrice == "yes" :
	MDtools.Matrice_RMSD(dPDB, mode_RMSD, outfile_RMSD)

#calculate interaction frequency between each residu and the given reference chain
MDtools.freq_interface(dPDB_ref, dPDB, chaine, seuil, dist_mode, outfile_freq)
