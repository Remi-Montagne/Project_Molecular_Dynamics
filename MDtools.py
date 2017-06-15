#!/usr/bin/python
# -*- coding: utf-8 -*-

import re
import matplotlib.pyplot as plt
from math import *

###==========================================================###
###			            PDB PARSER
###==========================================================###

#parse pdb files containing several frames
#returns a dictionnary containing : frame (integer) -> chain (string) -> residue (integer) -> atome (string) -> list of coordinates x, y, z (floats)
def load_PDB (fichier) :
	
	#open and stocks the file
	f = open(fichier, "r")
	lines = f.readlines()
	print("File %s. %d lines" % (fichier, len(lines)))
	

	dico = {}
	dico["liste_frames"] = []

	frame = 0
	#si l on est dans un fichier de reference, sans indication de frame, cette valeur reste a 0

	for line in lines:
		if (line[0:5]=="MODEL"): 
			#si l on est dans un nouveau frame, creation d un nouveau dictionnaire
			frame = int(line[6:14].strip())
			
		elif (line[0:4]=="ATOM"):
			coord = [float(line[31:39].strip()), float(line[39:47].strip()), float(line[47:55].strip()), 0.0]
			if frame not in dico:
				dico[frame] = {}
				dico[frame]["liste_chaines"] = []
				dico["liste_frames"].append(frame)
				#ajout de frame a la liste des frames 
				#et creation du dico des chaines contenues dans ce frame
			
			chaine = line[72:76].strip()
			if chaine not in dico[frame] :
				#creation d un nouveau dictionnaire correspondant a la chaine pour ce frame
				dico[frame][chaine] = {}
				dico[frame]["liste_chaines"].append(chaine)
				dico[frame][chaine]["liste_resids"] = []
			
			resid = int(line[22:27].strip())
			if resid not in dico[frame][chaine]:
				#creation d un nouveau dictionnaire correspondant au residu de cette chaine pour ce frame
				dico[frame][chaine][resid] = {}
				dico[frame][chaine]["liste_resids"].append(resid)
				dico[frame][chaine][resid]["liste_atomes"] = []
				
			nom_atome = line[12:17].strip()
			if nom_atome not in dico[frame][chaine][resid]:
				#stockage des coordonnees de l atome de ce residu, de cette chaine, et pour ce frame
				dico[frame][chaine][resid][nom_atome]=coord
				dico[frame][chaine][resid]["liste_atomes"].append(nom_atome)
				#print frame, chaine, resid, nom_atome, dico[frame][chaine][resid][nom_atome]
	
	f.close()
	
	return dico


###==========================================================###
###			    FUNCTIONS FOR DISTANCE
###==========================================================###

#returns distance between 2 atoms as a float	
def distance(l1, l2):
	return sqrt(((l1[0]-l2[0])**2 + (l1[1] - l2[1])**2 + (l1[2] - l2[2])**2))

#determine the mass center of an AA
#returns a list of 3 coordinates
def barycentre(residu) :
	bar = [0,0,0]
	for clef in residu["liste_atomes"] :
		bar[0] += residu[clef][0]
		bar[1] += residu[clef][1]
		bar[2] += residu[clef][2]
	bar[0] = bar[0] / len(residu["liste_atomes"])
	bar[1] = bar[1] / len(residu["liste_atomes"])
	bar[2] = bar[2] / len(residu["liste_atomes"])
	return bar

#determine the minimal distance between the atomes constituting 2 residues
#retourne un float
def distance_min_resid(resid1, resid2) : 
	compt = 0;
	for clef1 in resid1["liste_atomes"] :
		for clef2 in resid2["liste_atomes"] :
			if compt==0 :
				mini = distance(resid1[clef1],resid2[clef2])
				compt += 1
			else :
				if distance(resid1[clef1],resid2[clef2])< mini :
					mini = distance(resid1[clef1],resid2[clef2])
	return mini


###==========================================================###
###			                RMSD
###==========================================================###

#calculate the RMSD of every chain of the protein and the global RMSD global, for every frame
#returns the file : <outfile_RMSD>.txt
def RMSD(dPDB_ref, dPDB, mode_RMSD, outfile_RMSD) :
	print "Calcul of RMSD  mode:", mode_RMSD
	
	#ouverture du fichier des resultats
	nom = outfile_RMSD + ".txt"
	output = open(nom, "w")
	
	#liste des chaines a traiter pour le RMSD (les A mais pas les B et C)
	chaines_a_traiter = dPDB[10]["liste_chaines"][:]
	chaines_a_traiter.remove("B")
	chaines_a_traiter.remove("C")
	
	#ecriture de la 1e ligne (header)
	output.write("frame (mode:%s)" % mode_RMSD)
	for chaine in chaines_a_traiter :
		output.write(";%s" % chaine)
	output.write(";global\n")
	
	#pour la premiere ligne, le RMSD vaut 0
	output.write("0")
	for chaine in chaines_a_traiter :
		output.write(";0")
	output.write(";0\n")
	
	#pour chaque frame, on calcule le RMSD total et le RMSD de chaque chaine
	for frame in dPDB["liste_frames"] :
		
		output.write("%s" % str(frame))
			
		#initialisation a 0 de la somme globale des distances entre paires d'atomes (RMSD_global) et du nombre d'atomes, N_glob
		RMSD_glob = 0 
		N_glob = 0		
		
		#pour chaque chaine, on calcul le RMSD de la chaine et on implémente les valeurs globales
		for chaine in chaines_a_traiter :
			RMSD_chaine = 0
			N_chaine = 0
			
			for residu in dPDB[frame][chaine]["liste_resids"] :
				for atome in dPDB[frame][chaine][residu]["liste_atomes"] :		
					if mode_RMSD == "all" : #si mode all
						RMSD_chaine += distance(dPDB_ref[0][chaine][residu][atome], dPDB[frame][chaine][residu][atome])
						N_chaine += 1
						RMSD_glob += distance(dPDB_ref[0][chaine][residu][atome], dPDB[frame][chaine][residu][atome])
						N_glob += 1
					
					elif mode_RMSD == "alpha" : #si mode all
						if atome == "CA" :
							RMSD_chaine += distance(dPDB_ref[0][chaine][residu][atome], dPDB[frame][chaine][residu][atome])
							N_chaine += 1
							RMSD_glob += distance(dPDB_ref[0][chaine][residu][atome], dPDB[frame][chaine][residu][atome])
							N_glob += 1
					
					else : #si mode backbone
						if atome in ["CA", "C", "N", "O"] :
							RMSD_chaine += distance(dPDB_ref[0][chaine][residu][atome], dPDB[frame][chaine][residu][atome])
							N_chaine += 1
							RMSD_glob += distance(dPDB_ref[0][chaine][residu][atome], dPDB[frame][chaine][residu][atome])
							N_glob += 1
			
			#calcul et ecriture du RMSD de la chaine
			RMSD_chaine = sqrt(RMSD_chaine / N_chaine)
			output.write(";%s" % str(RMSD_chaine))
		
		#calcul et ecriture du RMSD global
		RMSD_glob = sqrt(RMSD_glob / N_glob)
		output.write(";%s\n" % str(RMSD_glob))
	
	#fermeture du fichier
	output.close()
	
	print("Outfile for RMSDs : %s" % nom)


#calculate the RMSD for a specified residue
#returns a dictionary containing a list of frames and a list of RMSDs
def RMSD_residu(dPDB_ref, dPDB, chaine, residu, mode_RMSD) :
	if chaine not in ["B", "C"] :	
		#creation du dictionnaire
		RMSD_resid = {}
		RMSD_resid["frames"] = []
		RMSD_resid["RMSD"] = []
		
		#pour chaque frame, on calcule le RMSD du residu
		for frame in dPDB["liste_frames"] :
			
			#ajout de la frame a la liste de frames
			RMSD_resid["frames"].append(frame)
			
			#initialisation a 0 de la somme des distances entre paires d'atomes et du nombre d'atomes, N_glob
			RMSD = 0 
			N = 0		
			#calcul du RMSD de la frame
			for atome in dPDB[frame][chaine][residu]["liste_atomes"] :	
				if mode_RMSD == "all" : #mode all
					RMSD += distance(dPDB_ref[0][chaine][residu][atome], dPDB[frame][chaine][residu][atome])
					N += 1
					
				elif mode_RMSD == "alpha" : #mode alpha
					if atome == "CA" :
							RMSD += distance(dPDB_ref[0][chaine][residu][atome], dPDB[frame][chaine][residu][atome])
							N += 1			
							
				else : 	#mode backbone
					if atome in ["CA", "C", "N", "O"] :
							RMSD += distance(dPDB_ref[0][chaine][residu][atome], dPDB[frame][chaine][residu][atome])
							N += 1
			
			#calcul et ajout a la liste RMSD
			RMSD = sqrt(RMSD / N)
			RMSD_resid["RMSD"].append(RMSD)
		
		return RMSD_resid

#calculate the RMSD between two frames for the global protein or a chain.
#retourne le resultat sous forme d'un entier
def RMSD_frame2frame(dPDB, frame_ref, frame_comparee, mode_RMSD, glob_ou_chaine) :
	#initialisation a 0 de la somme des distances entre paires d'atomes (RMSD) et du nombre d'atomes, N
	RMSD = 0
	N = 0
	
	#si le RMSD global est demande
	if glob_ou_chaine == "glob" :
		
		#liste des chaines a traiter pour le RMSD (les A mais pas les B et C)
		chaines_a_traiter = dPDB[10]["liste_chaines"][:]
		chaines_a_traiter.remove("B")
		chaines_a_traiter.remove("C")
			
		for chaine in chaines_a_traiter :
			for residu in dPDB[frame_comparee][chaine]["liste_resids"] :
				for atome in dPDB[frame_comparee][chaine][residu]["liste_atomes"] :
					if mode_RMSD == "all" :	#mode all
						RMSD += distance(dPDB[frame_ref][chaine][residu][atome], dPDB[frame_comparee][chaine][residu][atome])
						N+=1
				
					elif mode_RMSD == "alpha" :	#mode alpha
						if atome == "CA" :
							RMSD += distance(dPDB[frame_ref][chaine][residu][atome], dPDB[frame_comparee][chaine][residu][atome])
							N+=1
							
					else :	#sinon, mode backbone
						if atome in ["CA", "O", "C", "N"] :
							RMSD += distance(dPDB[frame_ref][chaine][residu][atome], dPDB[frame_comparee][chaine][residu][atome])
							N+=1
	
		RMSD = sqrt(RMSD / N)
		
	#si une chaine est indiquee
	else :
		chaine = glob_ou_chaine
		for residu in dPDB[frame_comparee][chaine]["liste_resids"] :
			for atome in dPDB[frame_comparee][chaine][residu]["liste_atomes"] :
				if mode_RMSD == "all" :	#mode all
					RMSD += distance(dPDB[frame_ref][chaine][residu][atome], dPDB[frame_comparee][chaine][residu][atome])
					N+=1
				
				elif mode_RMSD == "alpha" :	#mode alpha
					if atome == "CA" :
						RMSD += distance(dPDB[frame_ref][chaine][residu][atome], dPDB[frame_comparee][chaine][residu][atome])
						N+=1
							
				else :	#sinon, mode backbone
					if atome in ["CA", "O", "C", "N"] :
						RMSD += distance(dPDB[frame_ref][chaine][residu][atome], dPDB[frame_comparee][chaine][residu][atome])
						N+=1
	
		RMSD = sqrt(RMSD / N)
	
	return RMSD


#calculate a RMSD matrix for the global protein and matrices for each chain's RMSD
#taking each frame as a reference
#calls the RMSD_frame2frame function
#returns an outfile <outfile_RMSD>_Matrice_<chAi/globale>.txt
def Matrice_RMSD(dPDB, mode_RMSD, outfile_RMSD) :
	print"Matrices de RMSD"
	
	#Matrice RMSD de la proteine entiere
	
	#ouverture du fichier contenant la matrice
	nom = outfile_RMSD + "_Matrice_globale.txt"
	output = open(nom, "w")
	print "Matrice RMSD de la proteine entiere : ", nom
	
	for i in dPDB["liste_frames"] :			
		for j in dPDB["liste_frames"] :
			RMSD = RMSD_frame2frame(dPDB, i, j, mode_RMSD, "glob")
			output.write("%s;" % str(RMSD))
		output.write("\n")
	output.close()

	
	#Matrice RMSD de chaque chaine
	for chaine in dPDB[10]["liste_chaines"] :
		if chaine not in ["B", "C"] :
			
			#ouverture du fichier contenant la matrice
			nom = outfile_RMSD + "_Matrice_ch" + chaine +".txt"
			print "Matrice RMSD de la chaine " , chaine, " : ", nom
			output = open(nom, "w")
		
			for i in dPDB["liste_frames"] :	
				for j in dPDB["liste_frames"] :
					RMSD = RMSD_frame2frame(dPDB, i, j, mode_RMSD, chaine)
					output.write("%s;" % str(RMSD))
				output.write("\n")
			output.close()


###==========================================================###
###			           INTERFACE
###==========================================================###

#determine, for a list of residues, if they interac with each other (distance smaller than a given threshold in Angstroms)
#returns the degree of proximity, rangin from 0 to 3 with 0 = no interaction and 3 = very close
#to modes available: distance between mass centers or distance between every atoms of the residues
def interface(sequence, compar, seuil, dist_mode) :
	#liste des chaines sans celle comparee
	clef_chaine = sequence["liste_chaines"][:]
	clef_chaine.remove(compar)
	for chaine in clef_chaine :
		for resid1 in sequence[chaine]["liste_resids"] :
			if dist_mode == "bar" :
				bar1 = barycentre(sequence[chaine][resid1])
			for resid2 in sequence[compar]["liste_resids"] :
				if dist_mode == "mini" :
					dist = distance_min_resid(sequence[chaine][resid1], sequence[compar][resid2])
				else :
					dist = distance(bar1, barycentre(sequence[compar][resid2]))
				if dist <= (seuil - 4) :
					for atome in sequence[chaine][resid1]["liste_atomes"]:
						if sequence[chaine][resid1][atome][3] < 3 :
							sequence[chaine][resid1][atome][3] = 3				
				elif dist <= (seuil - 2):
					for atome in sequence[chaine][resid1]["liste_atomes"]:
						if sequence[chaine][resid1][atome][3] < 2 :
							sequence[chaine][resid1][atome][3] = 2
				elif dist <= seuil:
					for atome in sequence[chaine][resid1]["liste_atomes"]:
						if sequence[chaine][resid1][atome][3] < 1 :
							sequence[chaine][resid1][atome][3] = 1

#calculate, for a given number of frames, the frequency where a residue interact with a specified chain and returns the result
def freq_interface(ref_frame, frames, compar, seuil, dist_mode, outfile) :
	output = open(outfile, "w")
	
	#calcul de l appartenance a l interface des residus du pdb de reference
	interface(ref_frame[0], compar, seuil, dist_mode)
	
	for frame in frames["liste_frames"] :
		interface(frames[frame], compar, seuil, dist_mode)
		#calcul de l appartenance a l interface des residus de chaque frame du pdb a etudier
	
	for chaine in ref_frame[0]["liste_chaines"] :
		#s il ne s agit pas de la chaine de reference
		if re.match(compar, chaine) == None:
			for residu in ref_frame[0][chaine]["liste_resids"] :
				for atome in ref_frame[0][chaine][residu]["liste_atomes"]:
					freq = 0.0
					if ref_frame[0][chaine][residu][atome][3] > 0 :
						freq = 1.0
					for frame in frames["liste_frames"] :
						if chaine in frames[frame]["liste_chaines"] and residu in frames[frame][chaine]["liste_resids"] and atome in frames[frame][chaine][residu]["liste_atomes"] :
							if frames[frame][chaine][residu][atome][3] > 0 :
								#pour chaque residu, la frequence est incrementee si la troisieme valeur, qui correspond a l indication d appartenance a l interface retournee par la fonction interface, est non nulle
								freq = freq + 1.0
					
					output.write("%s\t%s\t%s\t%s\n"%(chaine, str(residu), atome, str(freq/float(len(frames)+1))))
	output.close()


#for every pair of AA of two given chains, returns the interaction frequency, between 0 and 1
#results are returned as a txt file giving the list of every interacting pair
def interfaces(dPDB_ref, dPDB, chaine1, chaine2, dist_mode, outfile_freq, seuil) :

	nom = outfile_freq + "ch" + chaine1 + "-" + chaine2 + ".txt"
	l = len(dPDB["liste_frames"]) +1
	output = open(nom, "w")

	print("interaction entre residus des chaines %s et %s, mode de calcul : %s" % (chaine1, chaine2, dist_mode))
	output.write("chaine%s;chaine%s;freq (seuil : %s)\n" % (chaine1, chaine2, seuil))
	
	liste_paires = []
	
	for residu1 in dPDB[10][chaine1]["liste_resids"] :					
		for residu2 in dPDB[10][chaine2]["liste_resids"] :
					
			freq = 0.0;
			
			#est-ce que la distance est inferieure au seuil	dans la reference
			if dist_mode == "bar" : #si dist_mode == bar
				if distance(barycentre(dPDB_ref[0][chaine1][residu1]),barycentre(dPDB_ref[0][chaine2][residu2])) < seuil :
					freq += 1
			else : #si dist_mode == mini
				if distance_min_resid(dPDB_ref[0][chaine1][residu1],dPDB_ref[0][chaine2][residu2]) < seuil :
					freq += 1
			
			#est-ce que la distance est inferieure au seuil	dans les autres frame
			if dist_mode == "bar" : #si dist_mode == bar
				for frame in dPDB["liste_frames"] :
					if distance(barycentre(dPDB[frame][chaine1][residu1]),barycentre(dPDB[frame][chaine2][residu2])) < seuil :
						freq += 1
			else : #si dist_mode == mini
				for frame in dPDB["liste_frames"] :
					if distance_min_resid(dPDB[frame][chaine1][residu1],dPDB[frame][chaine2][residu2]) < seuil :
						freq += 1	
							
			#calcul et ecriture de la frequence d'interaction
			freq = freq / l							
			output.write("%s;%s;%s\n" % (str(residu1), str(residu2), str(freq)))

			#si la frequence de contact n'est pas de 0, on memorise la paire d'atomes dans la liste liste_paire
			if freq != 0 :	
				paire = []
				paire.append(residu1)
				paire.append(residu2)
				liste_paires.append(paire)
	
	output.close()
		
	return liste_paires

#for every pair of residues of a given list, determine the distance between the distance
#and the RMSD of each residue. Returns a png file
def distance_resid2resid(dPDB_ref, dPDB, chaine1, chaine2, liste_paires, dist_mode, mode_RMSD) :
	for paire in liste_paires :	
		nom = "distance_ch" + str(chaine1) + "res" + str(paire[0]) + "-ch" + str(chaine2) + "res" + str(paire[1]) + ".png"
		print nom
		
		distances = distance_resid_all_frames(dPDB_ref, dPDB, chaine1, paire[0], chaine2, paire[1], dist_mode)
		
		plt.figure(1)
		plt.subplot(311)
		plt.plot(distances["frames"], distances["dist"], linestyle = 'solid', color = "red")
		plt.subplot(312)
		RMSD_resid = RMSD_residu(dPDB_ref, dPDB, chaine1, paire[0], mode_RMSD)
		plt.plot(RMSD_resid["frames"], RMSD_resid["RMSD"], linestyle = 'solid')
		if chaine2 not in ["B", "C"] :
			plt.subplot(313)
			RMSD_resid = RMSD_residu(dPDB_ref, dPDB, chaine2, paire[1], mode_RMSD)
			plt.plot(RMSD_resid["frames"], RMSD_resid["RMSD"], linestyle = 'solid')
		plt.savefig(nom)
		plt.close()

#calculate the distance between 2 residus for every frame
#returns the results as a dictionnary containing a list of frames and a list of distances
def distance_resid_all_frames(dPDB_ref, dPDB, chaine1, resid1, chaine2, resid2, dist_mode) : 
	distances = {}
	distances["frames"] = []
	distances["dist"] = []
	
	#distance dans la structure de reference
	distances["frames"].append(0)
	if dist_mode == "bar" : #dist_mode = "bar"
		distances["dist"].append(distance(barycentre(dPDB_ref[0][chaine1][resid1]),barycentre(dPDB_ref[0][chaine2][resid2])))
	else : #dist_mode = "mini"
		distances["dist"].append(distance_min_resid(dPDB_ref[0][chaine1][resid1],dPDB_ref[0][chaine2][resid2]))
	
	#distance pour les autres frames	
	for frame in dPDB["liste_frames"] :		
		distances["frames"].append(frame)
		if dist_mode == "bar" : #dist_mode = "bar"
			distances["dist"].append(distance(barycentre(dPDB[frame][chaine1][resid1]),barycentre(dPDB[frame][chaine2][resid2])))
		else :#dist_mode = "mini"
			distances["dist"].append(distance_min_resid(dPDB[frame][chaine1][resid1],dPDB[frame][chaine2][resid2]))
	
	return distances
