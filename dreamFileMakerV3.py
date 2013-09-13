#/usr/bin/python2.7
#-*- coding: cp1251 -*-
#-*- coding: iso-8859-15 -*-
#-*- coding: utf-8 -*- 

#///////////////////////////////////////////////////////////////////////
__author__ = "BEN HASSINE Najla(Najla.Ben-Hassine@versailles.inra.fr)"#/
__version__ = "2.0"				   		     						  #/
__copyright__ = "Copyright (c) 2013-2014 BHN"                         #/
__license__ = "GROUPE DEV IJPB"			                     		  #/
#///////////////////////////////////////////////////////////////////////


"""
VERSION PYTHON UTILISEE:  Python 2.7.3 (64-bit)
PROGRAMME: TRAITEMENT DU FICHIER .vcf OBTENUS PAR SNPEFF Version (snpEff_development)

MODULE: TRAITEMENT FICHIER """



#------------------------------------------------------------------------------------------------------------------- BIBLIOTHEQUE___DEBUT
#Bibliotheque systeme
import sys
import os
import getopt
import commands
import shutil

#Bibliotheque temps et date
import time
import datetime

#Bibliotheque expression reguliere
import re

#Bibliotheque gestion d'erreurs
import warnings
import logging
import logging.config

#------------------------------------------------------------------------------------------------------------------- BIBLIOTHEQUE___FIN


#------------------------------------------------------------------------------------------------------------------- CONFIGURATION___DEBUT
#DECLARATION DES VARIABLES GLOBALES
global IN_PUT_VCF
global DREAM_FILE
global LOG_FILE
global UNTREATED_CASES_FILE

#CREATION DE REPERTOIRE DE SORTIE  -------------------------------------------------------------
#FONCTION DE RECUPERATION DU NOM DU FICHIER -------------------------------------------------------------
def file_name_extraction (fileName):
	""" ****  FONCTION :  ****  \n\tfile_name_extraction : Extraire le nom du fichier sans l'extension ****  """
	fileNameList=fileName.split("/")
	nameFile= fileNameList[len(fileNameList)-1]
	nameFilewEx=nameFile.split(".")
	nameFileOnly=nameFilewEx[0]
	extension = nameFilewEx[1]
	#print nameFileOnly
	#print extension
	return nameFileOnly  
	
	
#FONCTION CREATION FICHIER LOG ---------------------------------------------------------------------------
def log_report(LOG_FILE_NAME):
 	""" ****  FONCTION :  ****  \n\tlog_report : FORMATAGE DU FICHIER LOG ****  """
	logging.basicConfig(
					level=logging.DEBUG,
                    format='%(asctime)s %(name)-8s %(levelname)-8s %(message)s',
		    		datefmt='%m/%d/%Y %I:%M:%S %p',
                    filename=LOG_FILE_NAME,
                    filemode='w')

#__OBTENIR L HEURE ET LA DATE ACTUELLE
now = datetime.datetime.now()

def header_of_dreamFile(DREAM_FILE):
	""" ****  FONCTION :  ****  \n\theader_of_dreamFile(DREAM_FILE) : Fonction qui ecrit l'entete du fichier de rêve. ****  """
	#ECRITURE DE L ENTETE DANS LE DREAM FILE
	fout_parsVCF=open(DREAM_FILE,"w")
	fout_parsVCF.write("CHROMOSOME\tPOSITION\tREFERENCE\tCHANGE\tCUSTOM_SNP\tCUSTOM_VARIATION\tCHANGE_TYPE\tDP4_REF_FWD\tDP4_REF_REV\tDP4_ALT_FWD\tDP4_ALT_REV\tDP4_TOTAL_REF\tDP4_TOTAL_ALT\tCOVERAGE\tHOM_HET\tQUALITY\tID_GENE\tEFFECT\tNBR_BASE\tOLD_AA/NEW_AA\tOLD_CODON/NEW_CODON\tCODON_NUM(CDS)\tCODON_DEGENERACY\n")
	fout_parsVCF.close()


#FONCTION D AIDE -----------------------------------------------------------------------------------------
def help_dreamFileMaker_script():
	print "**************"
	print  "H E L P  :"
	print "**************"
	print " usage : treat_OLE_snpeff.py  <Input_VCF_OutPutFile_From_Perl_Script>"
	print " -h    : help. "
	print " NB    : OLE, One Line Effect."

#------------------------------------------------------------------------------------------------------------------- CONFIGURATION___FIN


#------------------------------------------------------------------------------------------------------------------- TRAITEMENT___DEBUT
#FONCTION FORMATAGE GENE ID
def geneID_formatClean(gene):	
	ID_GENE = str(gene).replace("set(['","")
	ID_GENE = ID_GENE.replace("'])","")
	ID_GENE = ID_GENE.replace("Gene:","")
	ID_GENE = ID_GENE.replace("Gene_","")
	ID_GENE = ID_GENE.replace("Transcript:","")
	ID_GENE = ID_GENE.replace("exon:","")
	return ID_GENE
	
#FONCTION FORMATAGE GENE ID IN EFF
def geneID_format_EFF_OL(line):
	line=line.replace(".1|","|")
	line = line.replace("-Protein","")

	line = line.replace(".1-Protein","")
	line = line.replace(".2-Protein","")
	line = line.replace(".3-Protein","")
	line = line.replace(".4-Protein","")
	line = line.replace(".5-Protein","")
	line = line.replace(".6-Protein","")
	line = line.replace(".7-Protein","")
	line = line.replace(".8-Protein","")
	line = line.replace(".9-Protein","")

	line = line.replace(".1","")
	line = line.replace(".2","")
	line = line.replace(".3","")
	line = line.replace(".4","")
	line = line.replace(".5","")
	line = line.replace(".6","")
	line = line.replace(".7","")
	line = line.replace(".8","")
	line = line.replace(".9","")
	
	line= line.replace("Gene:","")
	line = line.replace("Gene_","")
	line = line.replace("Transcript:","")
	line = line.replace("exon:","")
	return line

#FONCTION FORMATAGE GENE ID
def geneID_format(info):
	ID_GENE = info[0].replace(")","")
	ID_GENE = ID_GENE.replace("Gene:","")
	ID_GENE = ID_GENE.replace("Gene_","")
	ID_GENE = ID_GENE.replace("Transcript:","")
	ID_GENE = ID_GENE.replace("exon:","")
	return ID_GENE
	


	
#FONCTION DE TRAITEMENT DU FICHIER DE SORTIE DU SCRIPT PERL / EFFET PAR LIGNE ---------------------
def trait_OneLineVcfFile(IN_PUT_VCF):
 	""" ****  FONCTION :  ****  \n\ttrait_OneLineVcfFile : Traitement du fihcier de sortie du script perl OneLineEffet **** \nRole: Eliminer les annotations transcript, exon , gene modele du même gene et garder que le ID_GENE."""

	linesList=[]
	logging.info("DEBUT ___TRAITEMENT DU FICHIER DE SORTIE DU SCRIPT PERL / EFFET PAR LIGNE ")
	print "DEBUT ___ TRAITEMENT DU FICHIER DE SORTIE DU SCRIPT PERL / EFFET PAR LIGNE "

	fin= open(IN_PUT_VCF,"r")
	lines = fin.readlines()

	for line in lines:
		line=geneID_format_EFF_OL(line)
		#print line
		linesList.append(line)
	fin.close()
	logging.info("Lignes de depart " + str(len(sorted(set(lines)))))
	print "Lignes de depart " + str(len(sorted(set(lines))))
	#print sorted(set(lines))
	logging.info("Lignes apres trait " + str(len(sorted(set(linesList)))))
	print "Lignes apres trait " + str(len(sorted(set(linesList))))
	#print sorted(set(linesList))
	logging.info("FIN ___TRAITEMENT DU FICHIER DE SORTIE DU SCRIPT PERL / EFFET PAR LIGNE ")
	print "FIN ___ TRAITEMENT DU FICHIER DE SORTIE DU SCRIPT PERL / EFFET PAR LIGNE "

	return linesList
	


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@------- FONCTION DE FORMATAGE DES COLONNES  -------@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#DEBUT
#FONCTION QUI DETERMINE L HOMOZYGOTIE
def homozygotie_determin(listCln9,UNTREATED_CASES_FILE):
	""" ****  FONCTION :  ****  \n\thomozygotie_determin : Fonction qui determine l'homozygotie. ****  """
	HOMOHET=listCln9[0]
	#print HOMOHET
				
	if HOMOHET =="1/1":
		HOMOZYGOTIE = "Hom"
		#print HOMOZYGOTIE
		return HOMOZYGOTIE
	elif HOMOHET =="0/1":
		HOMOZYGOTIE = "Het"
		#print HOMOZYGOTIE
		return HOMOZYGOTIE
	else:
		HOMOZYGOTIE="Untreated case"
		#print HOMOZYGOTIE
		finrecup = open(UNTREATED_CASES_FILE,"a")
		finrecup.write(str(geneList_details) + "\n")
		logging.info("WARNING : Untreated CASE. Unknown HOMOZYOGOTIE "+ str(listCln9))
		finrecup.close()
		return HOMOZYGOTIE
		
#--- FONCTION QUI CALCUL ET RETOURNE LE DP4 ----------------------------------
def calculDP4(listInfo,k):
	""" ****  FONCTION :  ****  \n\tcalculDP4 : Elle prend en argument une liste et un entier.\n
	-- La liste consiste en listInfo: champs info du fichier .vcf 
	-- L'entier, k et la colonne qui correspond au dp4 """
	listDP4=[]
	listDP4=listInfo[k].split(",")
	#print listDP4
	DP4_REF_FWD=int(listDP4[0].replace("DP4=",""))
	#print DP4_REF_FWD
	DP4_REF_REV=int(listDP4[1])
	#print DP4_REF_REV
	DP4_ALT_FWD=int(listDP4[2])
	#print DP4_ALT_FWD
	DP4_ALT_REV=int(listDP4[3])	
	#printDP4_ALT_REV
					
	DP4_TOTAL_REF=DP4_REF_FWD+DP4_REF_REV
	DP4_TOTAL_ALT=DP4_ALT_FWD+DP4_ALT_REV
	
	recupDP4 = str(DP4_REF_FWD)+"\t"+str(DP4_REF_REV)+"\t"+str(DP4_ALT_FWD)+"\t"+str(DP4_ALT_REV)+"\t"+str(DP4_TOTAL_REF)+"\t"+str(DP4_TOTAL_ALT)+"\t"
	#print recupDP4
	#print DP4_REF_FWD,DP4_REF_REV,DP4_ALT_FWD,DP4_ALT_REV,DP4_TOTAL_REF,DP4_TOTAL_ALT
	return recupDP4
	
#FONCTION QUI DETERMINE LE CHAMPS OLD_AA_NEW_AA DU FICHIER DE REVE
def recupAA(geneList_details):
	"""****  FONCTION :  ****  \n\trecupAA : recupere la liste des acide amines pour le champs OLD_AA_NEW_AA"""
	#---RECUPERATION AA
	#print EFFECT_NAME  # verif non_syno syno
	#print geneList_details
		
	listAA =list(geneList_details[3])
	#print listGene[3]
	AA1=listAA[0]
	#print AA1
	AA2=listAA[len(listAA)-1]
		
	regex = re.compile('[^0-9]')
	recup= regex.findall(listAA[len(listAA)-1]) 
	if  len(recup) == 0:
		#print AA1
		OLD_AA_NEW_AA = AA1 +"/"+ AA1
		return OLD_AA_NEW_AA
	elif  len(recup) == 1:
		#print AA2
		OLD_AA_NEW_AA = AA1 +"/"+ AA2
		return OLD_AA_NEW_AA
	else:
		OLD_AA_NEW_AA = "na"
		return OLD_AA_NEW_AA
		
#FONCTION QUI DETERMINE LE CODON_NUM_CDS DU FICHIER DE REVE
def recupCDS(geneList_details):
	"""****  FONCTION :  ****  \n\trecupCDS : recupere le numero du cds pour le champs CODON_NUM_CDS"""
	listAA =list(geneList_details[3])
	cdsnumVal=geneList_details[3].replace(listAA[0],"")
	cdsnumVal=cdsnumVal.replace(listAA[len(listAA)-1],"")
	CODON_NUM_CDS = cdsnumVal
	return CODON_NUM_CDS
	
#--- FONCTION QUI FORMATTE LE CHAMP EFFECT ----------------------------------
def recupEffect(clnEff,listInfo,UNTREATED_CASES_FILE):
	""" ****  FONCTION :  ****  \n\trecupEffect : Fonction qui recupere les effets par gene. ****  """
	#logging.info("DEBUT ___ ETAPES DE TRAITEMENT DU CHAMPS EFFET .VCF")
	#print ("DEBUT ___ ETAPES DE TRAITEMENT DU CHAMPS EFFET  .VCF")
	nbr_gene=0
	nbr_effet=0
	
	RECUP_Effect =  "Untreated Name Effect. Problem with annotations is suspectd."
	
	listEffdetails = []
	listEffect = []
	len_listEffect = []
	len_listGene = []
	geneList = []
	listLen =[]
	
	
	ID_GENE = "na"
	EFFECT_NAME= "na"
	NBR_BASE = "na"
	OLD_AA_NEW_AA = "na"
	OLD_CODON_NEW_CODON = "na"
	CODON_NUM_CDS = "na"
	CODON_DEGENERACY = "na"
	
	#RECUP_Effect = ID_GENE + "\t" + EFFECT_NAME + "\t"+  NBR_BASE +"\t" + OLD_AA_NEW_AA + "\t" + OLD_CODON_NEW_CODON + "\t" + CODON_NUM_CDS + "\t" + CODON_DEGENERACY

	#TOUT LE CHAMP EFFECT
	listEffect = listInfo[clnEff].replace("EFF=","")
	#print listEffect
	
	list_effect_details_0= ["DOWNSTREAM","UPSTREAM","UTR_3_PRIME","UTR_5_PRIME","UTR_5_DELETED","INTRON","SPLICE_SITE_ACCEPTOR","START_GAINED","SPLICE_SITE_DONOR"]

	list_effect_details_1=["NON_SYNONYMOUS_CODING","SYNONYMOUS_CODING","SYNONYMOUS_START","FRAME_SHIFT","CODON_CHANGE","CODON_INSERTION","CODON_CHANGE_PLUS_CODON_INSERTION","CODON_DELETION","CODON_CHANGE_PLUS_CODON_DELETION","STOP_GAINED","SYNONYMOUS_STOP","START_LOST","STOP_LOST","NON_SYNONYMOUS_START"]

	#list_effect_details_2 = ["SPLICE_SITE_DONOR"]
	
	#RECUPERATION DUN NOM DE L EFFET
	EFFECT_NAME_details = listEffect.split("(")
	EFFECT_NAME = EFFECT_NAME_details[0]
	#print EFFECT_NAME
	
	#RECUPERATION DE LA LISTE DES GENES
	geneList_details = listEffect.split("|")

	#VERIFICATION NBR DE CLN DANS LE CHAMPS EFFET
	#if len(geneList_details) == 11:
		#print geneList_details
		
	if EFFECT_NAME in list_effect_details_1:	
		#---RECUPERATION CODON
		OLD_CODON_NEW_CODON = geneList_details[2]
		#print OLD_CODON_NEW_CODON
		
		#---RECUPERATION CDS
		CODON_NUM_CDS = recupCDS(geneList_details)
		#print CODON_NUM_CDS
	
		#RECUPERATION DU ID GENE 
		geneList = geneList_details[5] +"*"+ geneList_details[8]
		geneList = geneList.split("*")
		geneList_F = set(geneList)
		
		#RECUPERATION AA
		OLD_AA_NEW_AA = recupAA(geneList_details)
		
		#DANS LE CAS OU LA LISTE DE GENE COMPREND QU UN SEUL / ON A VIRER LES AUTRES ANNOTE DU MEME GENE
		if len(geneList_F )< 2:
			#print geneList
			#print geneList_details
			ID_GENE = geneList[0]
			#print ID_GENE
			
			#RECUPERATION DU CODON DEGENERE
			CODON_DEGENERACY = geneList_details[len(geneList_details)-1].replace(")","")
			#print CODON_DEGENERACY
			 
			RECUP_Effect = ID_GENE + "\t" + EFFECT_NAME + "\t" + NBR_BASE +"\t" + OLD_AA_NEW_AA + "\t" + OLD_CODON_NEW_CODON + "\t" + CODON_NUM_CDS + "\t" + CODON_DEGENERACY
			#print RECUP_Effect   #ok
			return RECUP_Effect
			
		elif len(geneList_F) ==2 : 
			#print geneList 
			#RECUPERATION DU CODON DEGENERE
			CODON_DEGENERACY = geneList_details[len(geneList_details)-1].replace(")","")
			
			if	(geneList_details[len(geneList_details)-2]) != CODON_DEGENERACY :
					#print geneList_details
					recupInfoForEdition =  geneList_details[5] + "*" + geneList_details[9] + "--" + geneList_details[8] + "*" + geneList_details[10].replace(")","")
					recupInfoForEditionList =recupInfoForEdition.split("--")
					#print recupInfoForEditionList
					#print recupInfoForEditionList
					for recupinfo in recupInfoForEditionList:
						#print recupinfo
						getGeneCodon = recupinfo.split("*")
						#print getGeneCodon
						ID_GENE = getGeneCodon[0]
						#print ID_GENE
						CODON_DEGENERACY = getGeneCodon[1]
						#print getGeneCodon[1]
						
						RECUP_Effect = ID_GENE + "\t" + EFFECT_NAME + "\t"+  NBR_BASE +"\t" + OLD_AA_NEW_AA + "\t" + OLD_CODON_NEW_CODON + "\t" + CODON_NUM_CDS + "\t" + CODON_DEGENERACY
						#print RECUP_Effect   #ok
						return RECUP_Effect
			else:
					#print geneList_details				
					for ID_GENE in geneList:
						RECUP_Effect = ID_GENE + "\t" + EFFECT_NAME + "\t"+  NBR_BASE +"\t" + OLD_AA_NEW_AA + "\t" + OLD_CODON_NEW_CODON + "\t" + CODON_NUM_CDS + "\t" + CODON_DEGENERACY
						#print RECUP_Effect   #ok
						return RECUP_Effect
					
	elif EFFECT_NAME in list_effect_details_0:
		geneList = geneList_details[5] +"*"+ geneList_details[8]
		geneList = geneList.split("*")
		geneList_F = set(geneList)
		
		#DANS LE CAS OU LA LISTE DE GENE COMPREND QU UN SEUL / ON A VIRER LES AUTRES ANNOTE DU MEME GENE
		if len(geneList_F )< 2:
			#print geneList
			#print geneList_details
			ID_GENE = geneList[0]
			#print ID_GENE
			
			#RECUPERATION DU CODON DEGENERE
			CODON_DEGENERACY = geneList_details[len(geneList_details)-1].replace(")","")
			#print CODON_DEGENERACY
			
			#RECUPERATION DU NOMBRE DE BASE
			NBR_BASE =  geneList_details[2]
			
			if NBR_BASE == "":
				NBR_BASE = "na"
				#print EFFECT_NAME
				RECUP_Effect = ID_GENE + "\t" + EFFECT_NAME + "\t"+  NBR_BASE +"\t" + OLD_AA_NEW_AA + "\t" + OLD_CODON_NEW_CODON + "\t" + CODON_NUM_CDS + "\t" + CODON_DEGENERACY			
				#print RECUP_Effect #ok
				return RECUP_Effect
			else:
				#print EFFECT_NAME
				RECUP_Effect = ID_GENE + "\t" + EFFECT_NAME + "\t"+  NBR_BASE +"\t" + OLD_AA_NEW_AA + "\t" + OLD_CODON_NEW_CODON + "\t" + CODON_NUM_CDS + "\t" + CODON_DEGENERACY
				return RECUP_Effect
		# CAS DE DEUX GENES	
		elif len(geneList_F) ==2 : 
			#print geneList 
			#if "|" not in geneList_details[len(geneList_details)-1]:
			for ID_GENE in geneList:
				#print  ID_GENE
					
				#RECUPERATION DU CODON DEGENERE
				CODON_DEGENERACY = geneList_details[len(geneList_details)-1].replace(")","")
				#print CODON_DEGENERACY
					
				#RECUPERATION DU NOMBRE DE BASE
				NBR_BASE =  geneList_details[2]
				#print NBR_BASE
				if NBR_BASE == "":
					NBR_BASE = "na"
					#print EFFECT_NAME
					RECUP_Effect = ID_GENE + "\t" + EFFECT_NAME + "\t"+  NBR_BASE +"\t" + OLD_AA_NEW_AA + "\t" + OLD_CODON_NEW_CODON + "\t" + CODON_NUM_CDS + "\t" + CODON_DEGENERACY			
					#print RECUP_Effect #ok
					return RECUP_Effect
				else:
					#print EFFECT_NAME
					RECUP_Effect = ID_GENE + "\t" + EFFECT_NAME + "\t"+  NBR_BASE +"\t" + OLD_AA_NEW_AA + "\t" + OLD_CODON_NEW_CODON + "\t" + CODON_NUM_CDS + "\t" + CODON_DEGENERACY			
					#print RECUP_Effect #ok
					return RECUP_Effect
				
	elif EFFECT_NAME ==  "INTERGENIC" :
		#RECUPERATION DU CODON DEGENERE
		CODON_DEGENERACY = geneList_details[len(geneList_details)-1].replace(")","")
		#print geneList_details
		#print EFFECT_NAME
		
		RECUP_Effect = ID_GENE + "\t" + EFFECT_NAME + "\t"+  NBR_BASE +"\t" + OLD_AA_NEW_AA + "\t" + OLD_CODON_NEW_CODON + "\t" + CODON_NUM_CDS + "\t" + CODON_DEGENERACY
		#print RECUP_Effect   #ok	
		return RECUP_Effect
	
	elif EFFECT_NAME == "INTRAGENIC" :
		#print geneList_details
		ID_GENE = geneID_formatClean(geneList_details[5])
		#RECUPERATION DU CODON DEGENERE
		CODON_DEGENERACY = geneList_details[len(geneList_details)-1].replace(")","")
		#print EFFECT_NAME
		
		RECUP_Effect = ID_GENE + "\t" + EFFECT_NAME + "\t"+  NBR_BASE +"\t" + OLD_AA_NEW_AA + "\t" + OLD_CODON_NEW_CODON + "\t" + CODON_NUM_CDS + "\t" + CODON_DEGENERACY
		#print RECUP_Effect   #ok	
		return RECUP_Effect

	else:
		#print EFFECT_NAME
		finrecup = open(UNTREATED_CASES_FILE,"a")
		finrecup.write("Fonction de traitement du champs effet : "+ "\n")
		finrecup.write(str(geneList_details) + "\n")
		logging.info("WARNING : Untreated Name Effect. Problem with annotations is suspected. "+ str(geneList_details))
		finrecup.close()
		
		return RECUP_Effect
		
	#elif EFFECT_NAME == "INTRAGENIC" :
		#print geneList_details
	#logging.info("FIN ___ ETAPES DE TRAITEMENT DU CHAMPS EFFET .VCF")
	#print ("FIN ___ ETAPES DE TRAITEMENT DU CHAMPS EFFET  .VCF")

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@------- FONCTION DE GESTION D'AFFICHAGE  -------@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#DEBUT
def snp_affich(GENOME,POS,REF,ALT,QUALITY,listInfo,listCln9,DREAM_FILE,UNTREATED_CASES_FILE):
	""" ****  FONCTION :  ****  \n\tsnp_affich : Fonction qui recupere les colonnes correspondant à un SNP. ****  """
		
	CHANGE_TYPE ="SNP"
	COVERAGE = listInfo[0].replace("DP=","")
	
	#old version 		
	#CUSTOM_VARIATION = REF + "_" + ALT.replace(REF,"")
	
	#new version 
	CUSTOM_VARIATION = REF + "_" + ALT
	CUSTOM_SNP = GENOME + "_" + POS + "_" + CUSTOM_VARIATION

	RecupFirstClnVCF = GENOME+ "\t" + POS+ "\t" + REF + "\t" + ALT+ "\t" + CUSTOM_SNP + "\t" + CUSTOM_VARIATION + "\t" + CHANGE_TYPE + "\t"
	RecupLastClnVCF = COVERAGE + "\t" + homozygotie_determin(listCln9,UNTREATED_CASES_FILE) + "\t" + QUALITY.replace(".",",") + "\t"
	
	
	#OUVERTURE EN ECRITURE DU FICHIER DE REVE
	parsed_file=open(DREAM_FILE,"a")
	
	if  (len(listInfo)==7) and ("DP4" in listInfo[4]):
		editLine = RecupFirstClnVCF + calculDP4(listInfo,4)  + RecupLastClnVCF + recupEffect(6,listInfo,UNTREATED_CASES_FILE) +"\n"
		#print editLine
		parsed_file.write(editLine)
		
							
	elif  (len(listInfo)==7) and ("DP4" in listInfo[3]):
		editLine = RecupFirstClnVCF + calculDP4(listInfo,3) +  RecupLastClnVCF + recupEffect(6,listInfo,UNTREATED_CASES_FILE) +"\n"
		#print editLine
		parsed_file.write(editLine)
								
	elif  (len(listInfo)==8) and ("DP4" in listInfo[4]):
		editLine = RecupFirstClnVCF + calculDP4(listInfo,4)  +  RecupLastClnVCF + recupEffect(7,listInfo,UNTREATED_CASES_FILE) +"\n"
		#print editLine
		parsed_file.write(editLine)
		
	elif  (len(listInfo)==8) and ("DP4" in listInfo[3]):
		editLine =RecupFirstClnVCF + calculDP4(listInfo,3)  +  RecupLastClnVCF + recupEffect(7,listInfo,UNTREATED_CASES_FILE) +"\n"
		#print editLine
		parsed_file.write(editLine)	
		
								
	elif  (len(listInfo)==9) and ("DP4" in listInfo[4]):
		editLine = RecupFirstClnVCF + calculDP4(listInfo,4)  +  RecupLastClnVCF + recupEffect(8,listInfo,UNTREATED_CASES_FILE) +"\n"
		#print editLine
		parsed_file.write(editLine)

	elif (len(listInfo)==9) and ("DP4" in listInfo[3]):
		editLine = RecupFirstClnVCF + calculDP4(listInfo,3) + RecupLastClnVCF	+ recupEffect(8,listInfo,UNTREATED_CASES_FILE) +"\n"
		#print editLine
		parsed_file.write(editLine)
	else:
		finrecup = open(UNTREATED_CASES_FILE,"a")
		finrecup.write("Fonction affichage SNP : " + "\n")
		finrecup.write(str(RecupFirstClnVCF) + "\n")
		logging.info("WARNING : Untreated cases. "+ str(RecupFirstClnVCF))
		finrecup.close()
	parsed_file.close()
	
#--- FONCTION QUI AFFICHE LES CLN PAR INDEL  ----------------------------------
def affichClnIndel(GENOME,POS,REF,ALT,QUALITY,CHANGE_TYPE,listInfo,listCln9,DREAM_FILE,UNTREATED_CASES_FILE):
	COVERAGE = listInfo[1].replace("DP=","")
			
	#old version 		
	#CUSTOM_VARIATION = REF + "_" + ALT.replace(REF,"")
	
	#new version 
	CUSTOM_VARIATION = REF + "_" + ALT
	CUSTOM_SNP = GENOME + "_" + POS + "_" + CUSTOM_VARIATION

	RecupFirstClnVCF = GENOME+ "\t" + POS+ "\t" + REF + "\t" + ALT+ "\t" + CUSTOM_SNP + "\t" + CUSTOM_VARIATION + "\t" + CHANGE_TYPE + "\t"
	RecupLastClnVCF =  COVERAGE + "\t" +  homozygotie_determin(listCln9,UNTREATED_CASES_FILE) + "\t" + QUALITY.replace(".",",") + "\t"
	
	#OUVERTURE EN ECRITURE DU FICHIER DE REVE
	parsed_file=open(DREAM_FILE,"a")
	
	if  (len(listInfo)==7) and ("DP4" in listInfo[4]):
		editLine = RecupFirstClnVCF + calculDP4(listInfo,4) +  RecupLastClnVCF	+ recupEffect(6,listInfo,UNTREATED_CASES_FILE) +"\n"
		#print editLine
		parsed_file.write(editLine)
		
	elif  (len(listInfo)==8) and ("DP4" in listInfo[4]):
		editLine = RecupFirstClnVCF + calculDP4(listInfo,4)  +  RecupLastClnVCF	+ recupEffect(7,listInfo,UNTREATED_CASES_FILE) +"\n"
		#print editLine
		parsed_file.write(editLine)
								
	elif  (len(listInfo)==9) and ("DP4" in listInfo[5]):
		editLine = RecupFirstClnVCF + calculDP4(listInfo,5)  +  RecupLastClnVCF	+ recupEffect(8,listInfo,UNTREATED_CASES_FILE) +"\n"
		#print editLine
		parsed_file.write(editLine)

	elif  (len(listInfo)==9) and ("DP4" in listInfo[4]):
		editLine = RecupFirstClnVCF + calculDP4(listInfo,4)  +  RecupLastClnVCF	+ recupEffect(8,listInfo,UNTREATED_CASES_FILE) +"\n"
		#print editLine
		parsed_file.write(editLine)
							
	elif  (len(listInfo)==10) and ("DP4" in listInfo[5]):
		editLine = RecupFirstClnVCF + calculDP4(listInfo,5) +  RecupLastClnVCF + recupEffect(9,listInfo,UNTREATED_CASES_FILE) +"\n"
		#print editLine
		parsed_file.write(editLine)
	else:
		finrecup = open(UNTREATED_CASES_FILE,"a")
		finrecup.write("Fonction affichage INDEL : " + "\n")
		finrecup.write(str(RecupFirstClnVCF) + "\n")
		logging.info("WARNING : Untreated cases. "+ str(RecupFirstClnVCF))
		finrecup.close()		
	parsed_file.close()	
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@------- FONCTION DE GESTION D'AFFICHAGE  -------@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#FIN

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@------- FONCTION DE FORMATAGE DES COLONNES  -------@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#DEBUT

#FONCTION DE PARSING DU FICHIER VCF -------------------------------------------------------------
def parsing_vcf(IN_PUT_VCF,DREAM_FILE,UNTREATED_CASES_FILE):
	logging.info("DEBUT ___ ETAPES DE TRAITEMENT DU FICHIER .VCF")
	print ("DEBUT ___ ETAPES DE TRAITEMENT DU FICHIER .VCF")
	
	#FONCTION DE TRAITEMENT DU FICHIER DE SORTIE DU SCRIPT PERL / EFFET PAR LIGNE
	fipvcfWoh=trait_OneLineVcfFile(IN_PUT_VCF)

	#DECLARATION DE VARIABLES
	recupln=[]
	listInfo=[]
	lengthListInfo=[]
	liststartchamp=[]
	
	altList=[]
	listEffectLen =[]
	listEffTrait =[]
	listCln9=[]
	
	nbr_ln_mitochondria=0
	nbr_ln_chloroplast=0
	nbr_ln_vcf_recup=0
	nbr_ln_snp_recup=0

	nbr_ln_snp_1_alt=0
	nbr_ln_snp_2_alt=0	
	
	nbr_ln_indel_recup=0
	
	nbr_ln_ins_recup=0
	nbr_ln_del_recup=0
	nbr_ln_ins_1_alt=0
	nbr_ln_ins_2_alt=0
	nbr_ln_del_1_alt=0
	nbr_ln_del_2_alt=0
	
	nbr_ln_cansNontraite_recup=0
	
		
	#RECUPERATION DES COLONNES
	for line in fipvcfWoh:
		line=line.rstrip("\n\r")
		recupln=line.split("\t")
		
		if recupln[0] == "mitochondria" :
			nbr_ln_mitochondria= nbr_ln_mitochondria+1
				
		elif recupln[0] == "chloroplast":
			nbr_ln_chloroplast= nbr_ln_chloroplast+1
				
		elif recupln[0] != "mitochondria" and recupln[0] != "chloroplast":
			GENOME= recupln[0]
			POS = recupln[1]
			ID = recupln[2]
			REF= recupln[3]
			ALT= recupln[4]
			QUALITY = recupln[5]
			FILTER = recupln[6]
			INFO = recupln[7]
			
			#WHO =  recupln[8]
			#print WHO #:GT:PL:GQ
	
			#Comptage du nombre de lignes pour les chromosomes
			nbr_ln_vcf_recup= nbr_ln_vcf_recup+1
				
			#VERIFICATION DU NOMBRE DE COLONNE DANS LE CHAMPS INFOS
			listInfo=INFO.split(";")
			lengthListInfo.append(len(listInfo))
			
			#RECUPERATION CHAMPS POUR L HOMOZYGOTIE
			HOMOHET_BASE_N1_N2_N3 = recupln[9]
			listCln9=HOMOHET_BASE_N1_N2_N3.split(":")
			
			#MODIF RAPHAEL
			GENOME=GENOME.replace("Chr","")
			REF = REF.replace("*","X")
			
			#CAS UN SEUL ALLEL ALTERNATIF
			#CAS SNP
			if ("INDEL" not in listInfo[0]) and ("," not in ALT):
								
					#RECUPERATION DES COLONNES
					snp_affich(GENOME,POS,REF,ALT,QUALITY,listInfo,listCln9,DREAM_FILE,UNTREATED_CASES_FILE)
					
					#COMPTAGE DE VERIFICATION
					nbr_ln_snp_1_alt = nbr_ln_snp_1_alt + 1
					nbr_ln_snp_recup = nbr_ln_snp_recup + 1
			#CAS DEUX ALLELES ALTERNATIF
			#CAS SNP
			elif ("INDEL" not in listInfo[0]) and (","  in ALT):
										
					#fout_parsVCF
					altList=ALT.split(",")
					for alterallel in altList:
						ALT = alterallel

					#RECUPERATION DES COLONNES
					snp_affich(GENOME,POS,REF,ALT,QUALITY,listInfo,listCln9,DREAM_FILE,UNTREATED_CASES_FILE)
					
					#COMPTAGE DE VERIFICATION
					nbr_ln_snp_2_alt = nbr_ln_snp_2_alt + 1
					nbr_ln_snp_recup = nbr_ln_snp_recup + 1
			
			#CAS UN SEUL ALLEL ALTERNATIF
			#CAS INDEL
			elif ("INDEL" in listInfo[0]) and ("," not in ALT) :
				#CAS INSERTION	
				if len(ALT) > len(REF) :
					CHANGE_TYPE = "INS"	
					#----------
					#AFFICHAGE CLN
					affichClnIndel(GENOME,POS,REF,ALT,QUALITY,CHANGE_TYPE,listInfo,listCln9,DREAM_FILE,UNTREATED_CASES_FILE)
					#----------
					nbr_ln_ins_recup = nbr_ln_ins_recup + 1
					nbr_ln_ins_1_alt = nbr_ln_ins_1_alt +1			
	
				#CAS DELETION
				elif len(ALT) < len(REF) :
					CHANGE_TYPE = "DEL"
					#----------
					#AFFICHAGE CLN
					affichClnIndel(GENOME,POS,REF,ALT,QUALITY,CHANGE_TYPE,listInfo,listCln9,DREAM_FILE,UNTREATED_CASES_FILE)
					#----------												
					nbr_ln_del_1_alt = nbr_ln_del_1_alt +1
					nbr_ln_del_recup = nbr_ln_del_recup + 1						
				nbr_ln_indel_recup = nbr_ln_indel_recup + 1

					
			#CAS DEUX ALLELES ALTERNATIF
			#CAS INDEL	
			elif ("INDEL" in listInfo[0]) and ("," in ALT) :
				#CAS INSERTION
				if len(ALT) > len(REF) :
					CHANGE_TYPE = "INS"
					altList=ALT.split(",")
					for alterallel in altList:
						ALT = alterallel															
						#----------
						#AFFICHAGE CLN
						affichClnIndel(GENOME,POS,REF,ALT,QUALITY,CHANGE_TYPE,listInfo,listCln9,DREAM_FILE,UNTREATED_CASES_FILE)
						#----------
						nbr_ln_ins_2_alt = 	nbr_ln_ins_2_alt +1
						nbr_ln_ins_recup = nbr_ln_ins_recup + 1
				#CAS DELETION				
				elif len(ALT) < len(REF) :	
					CHANGE_TYPE = "DEL"
					altList=ALT.split(",")						
					for alterallel in altList:
						ALT = alterallel																
						#----------
						#AFFICHAGE CLN
						affichClnIndel(GENOME,POS,REF,ALT,QUALITY,CHANGE_TYPE,listInfo,listCln9,DREAM_FILE,UNTREATED_CASES_FILE)
						#----------
						nbr_ln_del_2_alt = nbr_ln_del_2_alt +1
						nbr_ln_del_recup = nbr_ln_del_recup + 1					
				nbr_ln_indel_recup = nbr_ln_indel_recup + 1
				
			else:
				#cas si la colonne 0 du champs info est differente de DP et INDEL
				finrecup = open(UNTREATED_CASES_FILE,"a")
				finrecup.write(str(RecupFirstClnVCF) + "\n")
				logging.info("WARNING : Untreated cases.First Cln In INFO different then DP/INDEL "+ str(RecupFirstClnVCF))
				finrecup.close()	
					
								
		else:
			# cas si le genome est different de (chloro, mito, 1-5)
			finrecup = open(UNTREATED_CASES_FILE,"a")
			finrecup.write(str(RecupFirstClnVCF) + "\n")
			logging.info("WARNING : Untreated cases.Unknown Genome "+ str(RecupFirstClnVCF))
			finrecup.close()				
	#print "CTRL/NOMBRE DE COLONNES DANS LE CHAMPS INFO : " + str(sorted(set(lengthListInfo)))
				
	logging.info("FIN ___ ETAPES DE TRAITEMENT DU FICHIER .VCF")
	print ("FIN ___ ETAPES DE TRAITEMENT DU FICHIER .VCF")	

#FILTRER LE DREAM_FILE
def	dreamFile_Filter(DREAM_FILE,PATH_IN_PUT):
	
	PATH_DREAM_FILE = os.path.dirname(DREAM_FILE).strip("/")
	#*** TRIE DU FICHIER PARSER POUR AVOIR LES UNIQUES
	#PREPARATION DU SORTE UNIQ, TRIER LE FICHIER PAR CHR, POS, ID_GENE
	cmd0="sort -k1n,1n -n -k2n,2n -n -k17n,17n -d "+ PATH_DREAM_FILE + "/" + DREAM_FILE.replace(PATH_IN_PUT,"")
	
	cmd1 = cmd0+"> "+ PATH_DREAM_FILE + "/" + DREAM_FILE.replace(PATH_IN_PUT,"").strip(".txt")+"_sorted.txt"
	os.system(cmd1)
	
	#OBTENTION DES UNIQUES
	cmd2 =" sort -u " + PATH_DREAM_FILE + "/" + DREAM_FILE.replace(PATH_IN_PUT,"").strip(".txt")+"_sorted.txt"
	cmd3 = cmd2 +"> "+ PATH_DREAM_FILE + "/" + DREAM_FILE.replace(PATH_IN_PUT,"").strip(".txt")+"_sorted_Uniq.txt"
	os.system(cmd3)	
	
	#SUPPRIMER FICHIER INTERMEDIAIRE 1
	os.system("rm " + PATH_DREAM_FILE + "/" + DREAM_FILE.replace(PATH_IN_PUT,"").strip(".txt")+"_sorted.txt")
	
	#CREATION DU FICHIER LISIBLE PAR LES BIOLOGISTES
	cmd4 = "sort -k1n,1n -n -k2n,2n -n -k17n,17n -d " + PATH_IN_PUT + "/" + DREAM_FILE.replace(PATH_IN_PUT,"").strip(".txt")+"_sorted_Uniq.txt"
	cmd5 = cmd4 + "> " + PATH_DREAM_FILE + "/" + DREAM_FILE.replace(PATH_IN_PUT,"").strip(".txt") + "_FINAL.txt"
	
	os.system(cmd5)
	
	#SUPPRIMER FICHIER INTERMEDIAIRE 2
	os.system("rm " +  PATH_DREAM_FILE + "/" + DREAM_FILE.replace(PATH_IN_PUT,"").strip(".txt")+"_sorted_Uniq.txt")

	#SUPPRIMER FICHIER DREAM_FILE NON FILTRE
	os.system("rm " +  PATH_DREAM_FILE + "/" + DREAM_FILE.replace(PATH_IN_PUT,""))
#------------------------------------------------------------------------------------------------------------------- TRAITEMENT___FIN

#------------------------------------------------------------------------------------------------------------------- TEST_PRELIM___DEBUT
#FONCTION DE VERIFICATION DES ARGUMENTS ------------------------------------------------------------------
def verif_arg_nbr(argv):
	""" ****  FONCTION :  ****  \n\tverif_arg_nbr : VERIFICATION DES ARGUMENT"""
    	#print sys.argv
    	if  (len(sys.argv) <= 1):	
    		help_dreamFileMaker_script()
    		sys.exit()
   	elif len(sys.argv) == 2:
		if sys.argv[1] == "-h":
			help_dreamFileMaker_script()
			sys.exit()
		elif ".vcf" in sys.argv[1]:
			#print 'Input VCF file is : ', sys.argv[1]
			return  sys.argv[1]
		elif ".vcf" not in sys.argv[1]:
			print 'Input file : '+ sys.argv[1]+' IS NOT A .vcf'
   	elif len(sys.argv) == 3:
		if sys.argv[1] == "-h" and sys.argv[2] == "log_report" :
			print log_report.__doc__
		elif sys.argv[1] == "-h" and sys.argv[2] == "verif_arg_nbr" :
			print verif_arg_nbr.__doc__

#FONCTION VERFICIATION DE L EXISTENCE D UN FICHIER -------------------------------------------------------
def verif_existing_file(fileName):
	""" ****  FONCTION :  ****  \n\tverif_existing_file : VERIFICATION DE L EXISTENCE DU FICHIER .vcf"""
	pathFileName= os.path.exists(fileName)
	checkExistingFileValue=""
	#print pathFileName
	#print fileName
	if str(pathFileName) == 'True':
		checkExistingFileValue="OK"
		#print checkExistingFileValue
		msgExistingFile= "OK : The file : "+fileName+" exist"
		print str(msgExistingFile)
		logging.info(msgExistingFile)
	elif pathFileName != "True":
		msgExistingFile=  "[ERROR] : The file : "+fileName+" does not exist"
		print str(msgExistingFile)
		logging.error(msgExistingFile)
		checkExistingFileValue="ERROR"
		sys.stderr.write(msgExistingFile)
	return checkExistingFileValue

#FONCTION VERFICIATION QUE LE FICHIER N'EST PAS VIDE ----------------------------------------------		
def verif_file_not_empty(fileName):
	""" ****  FONCTION :  ****  \n\tverif_file_not_empty : VERIFICATION QUE LE  FICHIER .vcf N'EST PAS VIDE"""
	sizeNameFile=os.path.getsize(fileName)
	if sizeNameFile == 0:
		msgSizeFile=  "[ERROR] : The file : "+fileName+" is empty"
		print msgSizeFile
		logging.error(msgSizeFile)
		sys.stderr.write(msgExistingFile)
		return sizeNameFile
	else:
		msgSizeFile=  "OK : The file : "+fileName+" is NOT empty"
		print msgSizeFile
		logging.info(msgSizeFile)
		return sizeNameFile

#------------------------------------------------------------------------------------------------------------------- MAIN___DEBUT
# MAIN AVEC UN ARGUMENT A L ENTREE / script.py fichier.vcf
def main(argv):
	#FONCTION DE VERIFICATION DES ARGUMENTS --
	IN_PUT_VCF = verif_arg_nbr(argv)

	PATH_IN_PUT = os.path.dirname(IN_PUT_VCF).strip("/")
	#print PATH_IN_PUT
	if (PATH_IN_PUT == ""):
		PATH_IN_PUT=PATH_IN_PUT.replace("",".")
		
	#print PATH_IN_PUT
	#CREATION DE FICHIER DE SORTIE  ----------------------------------------------------------------
	#__ CREATION DU FICHIER LOG  --
	LOG_FILE=PATH_IN_PUT+"/LogFile_"+file_name_extraction(IN_PUT_VCF)+".log"
	log_report(LOG_FILE)
	#__ CREATION DU FICHIER CONTENANT LES CAS NON TRAITES  --
	UNTREATED_CASES_FILE = PATH_IN_PUT+"/UntreatedCasesFile_"+file_name_extraction(IN_PUT_VCF)+".txt"
	
	#__ CREATION DU FICHIER TXT FORMATE A PARTIR DU VCF  --
	#DREAM_FILE = PATH_IN_PUT+"/"+file_name_extraction(IN_PUT_VCF)+"_"+str(now.year)+str(now.month)+str(now.day)+"_DF.txt"
	DREAM_FILE = PATH_IN_PUT+"/"+file_name_extraction(IN_PUT_VCF)+"_"+"DF.txt"
	
	#__ RECUPERATION DE L ENTETE
	header_of_dreamFile(DREAM_FILE)	
	
	#FONCTION DE VERIFICATION DE LA VALIDITE DU FICHIER .vcf --
	logging.info("DEBUT ___ ETAPES DE VERIFICATION DES ARGUMENTS")
	print "DEBUT ___ ETAPES DE VERIFICATION DES ARGUMENTS"
	checkExistingFileValue = verif_existing_file(IN_PUT_VCF)
	sizeNameFile = verif_file_not_empty(IN_PUT_VCF)
	
	if (checkExistingFileValue  == "OK") and (sizeNameFile != 0):
		logging.info("OK : VCF INPUT IS VALID.")
		print "OK : VCF INPUT IS VALID."
		logging.info("FIN ___ ETAPES DE VERIFICATION DES ARGUMENTS")
		print "FIN ___ ETAPES DE VERIFICATION DES ARGUMENTS"
		
		#FONCTION DE PARSING DU FICHIER VCF
		parsing_vcf(IN_PUT_VCF,DREAM_FILE,UNTREATED_CASES_FILE)

		#FILTRER LE FICHIER FORMATER
		dreamFile_Filter(DREAM_FILE,PATH_IN_PUT)
		
	else:
		logging.info("FIN ___ ETAPES DE VERIFICATION DES ARGUMENTS")
		print "FIN ___ ETAPES DE VERIFICATION DES ARGUMENTS"
		#__ APPEL DE LA FONCTION AIDE 
		help_dreamFileMaker_script()
		sys.exit()	

#------------------------------------------------------------------------------------------------------------------- MAIN___FIN


#------------------------------------------------------------------------------------------------------------------- MAIN_EXECUTION__DEBUT
if __name__ == "__main__":
	argv = sys.argv
	main(sys.argv[1:])
#------------------------------------------------------------------------------------------------------------------- MAIN_EXECUTION___FIN


