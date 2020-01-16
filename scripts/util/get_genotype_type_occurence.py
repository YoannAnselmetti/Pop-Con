#!/usr/bin/env python
"""
Author:              Yoann Anselmetti
Last modification:   2019/12/13

Goal: From a "genotype_profiles_distrib"* file get occurences of the SNP types
		and occurence of the genotype types in function of the SNP type
INPUT:
	+ "genotype_profiles_distrib"* file
	+ output file where GT type count and heterozygosity stats will be written
OUTPUT:
	For all/each individual(s):
		+ Sites categories occurence:
			- All sites
			- Sites with at least 1 individual genotyped
			- Sites with at least 1 alternative allele genotyped
			- Sites where all individuals are genotyped
		+ Genotype types occurence in function of the sites category
		+ Heterozygosity statistics: H(obs) & H(exp)

License: This software is distributed under the CeCILL free software license
			(Version 2.1 dated 2013-06-21)
"""

from __future__ import print_function

### Global python packages import
from sys import argv
from os import path, makedirs
from datetime import datetime
from decimal import Decimal
import errno


# Taken from https://stackoverflow.com/a/600612/119527
def mkdir_p(PATH):
	''' Recusive mkdir (equivalent mkdir -p in bash)
	'''
	try:
		makedirs(PATH)
	except OSError as exc: # Python >2.5
		if exc.errno == errno.EEXIST and path.isdir(PATH):
			pass
		else:
			raise
# Taken from https://stackoverflow.com/a/23794010
def safe_open_w(PATH):
	''' Open "path" for writing, creating any parent directories as needed.
	'''
	try:
		mkdir_p(path.dirname(PATH))
	except:
		pass

	return open(PATH, 'w')



def get_GT(GT_list):
	'''
		GOAL:
			Get occurence of the different genotype type.
			Only consider valid genotype: "0/0", "0/1" and "1/1".
	'''
	### x="0/0" / y="0/1" / z="1/1"
	x,y,z=0,0,0
	for gt in GT_list:
		if gt=="0/0":
			x+=1
		if gt=="0/1":
			y+=1
		if gt=="1/1":
			z+=1
	return Decimal(x),Decimal(y),Decimal(z)



def get_allGT(occurence,GT_list,dict_GT):
	'''
		GOAL:
			Count genotype types in "GT_list" and store it in "dict_GT".
			For all individuals and for each individual 
	'''
	indID=1
	for gt in GT_list:
		######
		### Initialize dict_GT
		######
		### For all individuals
		if not "All" in dict_GT:
			dict_GT["All"]=dict()
		### For individual $(indID)
		if not indID in dict_GT:
			dict_GT[indID]=dict()

		######
		### Add genotype count to dict_GT 
		######
		### For all individuals
		if not gt in dict_GT["All"]:
			dict_GT["All"][gt]=0
		dict_GT["All"][gt]+=occurence
		### For individual $(indID)
		if not gt in dict_GT[indID]:
			dict_GT[indID][gt]=0
		dict_GT[indID][gt]+=occurence

		indID+=1

	return dict_GT



def compute_heterozygosity(occurence,GT_list,dict_hetero):
	'''
		GOAL:
			Compute H(obs) & H(exp) from a list of genotypes.
			!!! WARNING !!!
			With this approach if 1 or more individuals are ungenotyped
			or didn't pass the filtering steps (TAG: lowCov & lowMARF)
			The computation is biased as it considers sites as complete
			(all individuals genotyped) while it is not!!!
			=> Need extension to take into account this bias.

			=> For now, we are only confident with case where
			   all individuals are genotyped and passed the filtering steps!
	'''
	### Get x="0/0", y="0/1" and z="1/1" count in "GT_list"
	x,y,z=get_GT(GT_list)
	if x==y==z==0:
		### Sites where no individual is genotyped 
		### and/or no genotype passed the filtering steps (MRCt & MARFt)
		pass
	else:
		### x="0/0" / y="0/1" / z="1/1"
		Hobs=y/(x+y+z)
		Hexp=(Decimal(2.0)*((Decimal(2.0)*x+y)/(Decimal(2.0)*(x+y+z))) \
				*((Decimal(2.0)*z+y)/(Decimal(2.0)*(x+y+z))))

		if not "All" in dict_hetero:
			dict_hetero["All"]=dict()
		if not "Hobs" in dict_hetero["All"]:
			dict_hetero["All"]["Hobs"]=Decimal(0)
			dict_hetero["All"]["Hexp"]=Decimal(0)
			dict_hetero["All"]["sites"]=0

		dict_hetero["All"]["Hobs"]+=(Hobs*occurence)
		dict_hetero["All"]["Hexp"]+=(Hexp*occurence)
		dict_hetero["All"]["sites"]+=occurence

		######
		### COMPUTE HETEROZYGOSITY PER INDIVIDUAL
		######
		indID=1
		for gt in GT_list:
			list_GT=[gt]
			### Get x="0/0", y="0/1" and z="1/1" count in "gt" -> "list_GT" 
			x,y,z=get_GT(list_GT)
			if x==y==z==0:
				### Individual $(indID) is not genotyped in this profile 
				### or genotype didn't pass the filtering steps (MRCt & MARFt)
				pass
			else:
				### x="0/0" / y="0/1" / z="1/1"
				Hobs=y/(x+y+z)
				Hexp=(Decimal(2.0)*((Decimal(2.0)*x+y)/(Decimal(2.0)*(x+y+z))) \
						*((Decimal(2.0)*z+y)/(Decimal(2.0)*(x+y+z))))

				if not indID in dict_hetero:
					dict_hetero[indID]=dict()
				if not "Hobs" in dict_hetero[indID]:
					dict_hetero[indID]["Hobs"]=Decimal(0)
					dict_hetero[indID]["Hexp"]=Decimal(0)
					dict_hetero[indID]["sites"]=0

				dict_hetero[indID]["Hobs"]+=(Hobs*occurence)
				dict_hetero[indID]["Hexp"]+=(Hexp*occurence)
				dict_hetero[indID]["sites"]+=occurence

			indID+=1

	return dict_hetero



def write_heterozygosity_file(sites_type,dict_GT,dict_hetero,out_file):
	"""
		GOAL:	Write heterozygosity file with heterozygosity statistics
				and genotype type count for all and each individual(s)
		INPUT:
			+ sites_type -> Set of sites used as input for stats computation
			+ dict_GT -> dict() containing stats on genotype type count
			+ dict_hetero -> dict() containing stats on heterozygosity
		OUTPUT:
			+ 4 conditions:
				- with all sites
				- with sites where at least 1 individual is genotyped
				- with sites where at least 1 individual is genotyped
				  with at least 1 alternative alle
				- with sites where all individuals are genotyped
			+ For each condition:
				- #(0/0), #(0/1), #(1/1) for all/each individual(s)
				- H(obs) & H(exp) for all/each individual(s)
	"""
	print("\n\tFor "+sites_type+":")


	for indID in sorted(dict_GT):
		out_file.write(sites_type)
		if isinstance(indID,int):
			print("\t\tIndividual "+str(indID)+":")
			out_file.write("\tIndividual_"+str(indID))
		else:
			print("\t\tAll individuals:")
			out_file.write("\tAll_individuals")

		######
		### PRINT/WRITE GENOTYPE TYPE COUNT
		### FOR ALL INDIVIDUALS AND FOR EACH INDIVIDUALS
		######
		print("\t\t\tGENOTYPE TYPE COUNT:")
		for gt in sorted(dict_GT[indID]):
			occ=dict_GT[indID][gt]
			print("\t\t\t\t"+gt+" -> "+str(occ))
			if gt in ["0/0","0/1","1/1"]:
				out_file.write("\t"+str(occ))

		######
		### COMPUTE AND PRINT/WRITE HETEROZYGOSITY STATISTICS
		### FOR ALL INDIVIDUALS AND FOR EACH INDIVIDUAL
		######
		print("\t\t\tHETEROZYGOSITY STATISTICS:")
		sum_Hobs=dict_hetero[indID]["Hobs"]
		sum_Hexp=dict_hetero[indID]["Hexp"]
		sites=dict_hetero[indID]["sites"]

		Hobs=sum_Hobs/Decimal(sites)
		Hexp=sum_Hexp/Decimal(sites)

		print("\t\t\t\tH(obs): "+'{:.2e}'.format(Hobs))
		out_file.write("\t"+format(Hobs))
		print("\t\t\t\tH(exp): "+'{:.2e}'.format(Hexp))
		out_file.write("\t"+str(Hexp)+"\n")



def parse_genotype_profiles_distrib_file(gt_distrib_file,out_file):
	"""
		GOAL:	Parse genotype profiles distribution file to get statistics
				on heterozygosity/individual (and for all individuals)
		INPUT:
			+ genotype profiles distribution file
		OUTPUT:
			+ 4 conditions:
				- with all sites
				- with sites where at least 1 individual is genotyped
				- with sites where at least 1 individual is genotyped
				  with at least 1 alternative alle
				- with sites where all individuals are genotyped
			+ For each condition:
				- #(0/0), #(0/1), #(1/1) for all/each individual(s)
				- H(obs) & H(exp) for all/each individual(s)
	"""
	totSNP,SNP_1gt,SNP_1gtALT,SNP_allGT=0,0,0,0

	dict_GT_allSites=dict()
	dict_GT_1GT=dict()
	dict_GT_1GTwithALT=dict()
	dict_GT_allGT=dict()
	
	dict_hetero_allSites=dict()
	dict_hetero_1GT=dict()
	dict_hetero_1GTwithALT=dict()
	dict_hetero_allGT=dict()

	### 
	for line in open(gt_distrib_file,"r"):
		allGT=True
		if ("low" in line) or ("./." in line):
			allGT=False
		occurence=int(line.split()[0])
		GT_list=line.split()[1].split(sep)


		### Genotype analysis of the SNP/site to determine the SNP/site type
		for gt in GT_list:
			### For ALL SITES
			totSNP+=occurence
			dict_GT_allSites=get_allGT(occurence,GT_list,dict_GT_allSites)
			dict_hetero_allSites=compute_heterozygosity( \
									occurence, \
									GT_list, \
									dict_hetero_allSites)
			if not "low" in gt and gt!="./.":
				### For sites where at least 1 individuall is genotyped
				SNP_1gt+=occurence
				dict_GT_1GT= \
				get_allGT(occurence,GT_list,dict_GT_1GT)
				
				dict_hetero_1GT=compute_heterozygosity(
										occurence, \
										GT_list, \
										dict_hetero_1GT)
				if gt!="0/0":
					### For sites where at least 1 individual is genotyped
					### and with at least 1 alternative allele
					SNP_1gtALT+=occurence
					dict_GT_1GTwithALT=get_allGT( \
							occurence,GT_list,dict_GT_1GTwithALT)
					
					dict_hetero_1GTwithALT=compute_heterozygosity( \
						occurence,GT_list,dict_hetero_1GTwithALT)
			if allGT:
				### For sites where all individuals are genotyped
				SNP_allGT+=occurence
				dict_GT_allGT=get_allGT(occurence,GT_list,dict_GT_allGT)
				dict_hetero_allGT=compute_heterozygosity( \
									occurence, \
									GT_list, \
									dict_hetero_allGT)


	print("\nSNP TYPE STATISTICS:")
	print("\t",totSNP,"total SNPs")
	print("\t",SNP_1gt,"SNPs with at least 1 individual genotyped")
	print("\t",SNP_1gtALT,"SNPs with at least 1 individual genotyped \
and with at least 1 alternative allele")
	print("\t",SNP_allGT,"SNPs with all individuals genotyped")


	out_file.write("Sites_type\tIndividual_ID\t0/0\t0/1\t1/1\tH(obs)\tH(exp)\n")
	print("\nGENOTYPE TYPE STATISTICS:")
	write_heterozygosity_file("All_sites", \
		dict_GT_allSites,dict_hetero_allSites,out_file)
	write_heterozygosity_file("At_least_1_individual_genotyped", \
		dict_GT_1GT,dict_hetero_1GT,out_file)
	write_heterozygosity_file("At_least_1_alternative_allele", \
		dict_GT_1GTwithALT,dict_hetero_1GTwithALT,out_file)
	write_heterozygosity_file("All_individuals_genotyped", \
		dict_GT_allGT,dict_hetero_allGT,out_file)
	out_file.close()



################
###   MAIN   ###
################
if __name__ == '__main__':

	start_time=datetime.now()

	sep=","

	gt_distrib_file=argv[1]
	output_file=argv[2]


	# dict_allSites,dict_1GT,dict_1GTwithALT,dict_allGT=\
	parse_genotype_profiles_distrib_file(gt_distrib_file,safe_open_w(output_file))


	end_time=datetime.now()
	print('\nDuration: '+str(end_time-start_time))
