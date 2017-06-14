#!/usr/bin/python
# -*- coding: utf-8 -*-

###==========================================================###
###			DESCRIPTION
###==========================================================###

# This program calculates the frequence of interaction between each pair of Amino Acids (AA)
# from 2 proteic chains, in 3D conformation simulation.
# Inputs : a reference conformation and a simulation, at PDB format. 
#		   the name of the 1st and 2nd chains to compare.
#          the distance to consider : distance between AA mass centers ("bar") or minimum distance between atoms of the AAs
#		   the threshold under which 2 AA are said to interact.
# Outputs : frequence of interaction of AAs from 2 chains in a txt file. 1 txt file for each pair of chains
#           graphic representation in a png format


###==========================================================###

import sys, MDtools

###==========================================================###
###			USAGE
###==========================================================###

def usage():
	print """

	mandatory:
	===========
     
	-pdbref	-> reference pdb file

	-pdb	-> simulation pdb file

	-chaine1 -> 1e chain to compare
	
	-chaine2 -> 2e chain to compare

	optional:
	=========
	-seuil   -> threshold under which 2 AAs interac, in Angstrom, defaults to 5
     
	-dist    -> distance to use : "bar" for the distance between residues' mass centers, "mini" for minimum distance between 2 AAs atoms (significantly longer)
	
	-RMSD	 -> the way RMSD is calculated: "all" for each atom, "alpha" for alpha carbons, "backbone" for the AA's backbone (CA, C, O, N).

	-ofreq	-> name of the outfile indicating frequences of interaction (defaults to freq_interface.txt)

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

#chaine1	
try:
	chaine1 = str(sys.argv[sys.argv.index("-chaine1")+1])
except:
	usage()
	print "ERROR: Please enter 1st chain to compare."
	sys.exit()

#chaine2	
try:
	chaine2 = str(sys.argv[sys.argv.index("-chaine2")+1])
except:
	usage()
	print "ERROR: Please enter 2nd chain to compare."
	sys.exit()

#seuil
try:
	seuil = float(sys.argv[sys.argv.index("-seuil")+1])
except:
	seuil = 4.0
	
#dist
try:
	dist_mode = sys.argv[sys.argv.index("-dist")+1]
except:
	dist_mode = "bar"

#RMSD
try:
	mode_RMSD = sys.argv[sys.argv.index("-RMSD")+1]
except:
	mode_RMSD = "all"

#name of the output file
try:
	outfile_freq = sys.argv[sys.argv.index("-ofreq")+1]
except:
	outfile_freq = "freq_resid.txt"


###==========================================================###
###			       MAIN
###==========================================================###

#get coordinates of the atoms in dictionaries
dPDB_ref = MDtools.load_PDB(infile_ref)
dPDB = MDtools.load_PDB(infile)

#calculate the interaction frequency of every residues pair in the specified chains 
#write outpu file freq_interfacech<chaine1>-<chaine2>.txt
#returns the list of interacting AAs
liste_paires = MDtools.interfaces(dPDB_ref, dPDB, chaine1, chaine2, dist_mode, outfile_freq, seuil)

#returns the plots of the distance between interacting residues for each pair of interacting residues. Also returns a graphical representation of RMSD of the AAs.
#renvoie des fichiers .png 
MDtools.distance_resid2resid(dPDB_ref, dPDB, chaine1, chaine2, liste_paires, dist_mode, mode_RMSD)
