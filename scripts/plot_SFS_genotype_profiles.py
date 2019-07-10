#!/usr/bin/env python
###
###   Author:              Yoann Anselmetti
###   Last modification:   2019/07/09
###
###   Goal: Plot SFS with genotype profiles
###
###   License: This software is distributed under the CeCILL free software license (Version 2.1 dated 2013-06-21)
###

from __future__ import print_function

### Global python packages import
from sys import argv
from os import path
from datetime import datetime
from decimal import Decimal
import argparse

filePATH=path.dirname(path.realpath(__file__))
import sys
print(filePATH)
sys.path.append(filePATH+"/..")
### Personnal package import
from popcon import SFS_profiles_plot, util



######################
###   PARAMETERS   ###
######################
parser = argparse.ArgumentParser(prog='SFS plot with genotype profiles', description='Additional Pop-Con module to plot SFS plot with genotype profiles from heterozygosity file produced by Pop-Con with option "--heterozygosity_file".', epilog='''Source code and manual: http://github.com/YoannAnselmetti/Pop-Con\n\n''', formatter_class=argparse.RawDescriptionHelpFormatter)
## REQUIRED
parser.add_argument('-i', "--input_file", dest='heterozygosity_file', type=str, required=True, help='Heterozygosity file produced by Pop-Con.')
parser.add_argument('-var', "--variant_type", dest='variant_type', type=str, required=True, help='Variant type: "indel" or "SNP"')
### output
parser.add_argument('-p', "--prefix", dest='prefix', type=str, default='exp1', help='Experiment name (used as prefix for output files). (Default: "exp1")')
parser.add_argument('-v', "--verbose", dest='verbose', type=int, default=1, help='Verbose intensity. (Default: 1)')
parser.add_argument('-o', "--output", dest='output_dir', type=str, default=".", help='Output directory path. (Default: ./)')
### SFS plot
parser.add_argument('-sep', "--separator", dest='sep', type=str, default=",", help='Separator used in genotype profiles. (Default: ",")')
parser.add_argument('-max', "--max_profiles", dest='max_profiles', type=int, default=10, help='Maximum number of genotype profiles displayed in SFS plot. (Default: 10)')



################
###   MAIN   ###
################
if __name__ == '__main__':
   start_time=datetime.now()

   if len(argv)==1:
      argv.append('-h')

   args=parser.parse_args()

   heterozygosity_file=args.heterozygosity_file
   variant_type=args.variant_type

   prefix=args.prefix
   OUTPUT=args.output_dir

   verbose=args.verbose
   sep=args.sep
   max_profiles=args.max_profiles

   util.mkdir_p(OUTPUT)


   #########
   ### For SNP
   #########
   if variant_type=="SNP":
      if verbose>=1:
         print('\nVariant type is "SNP"')
         print('\t\tA/ Store genotypes infos with filtering TAG')
      ### Set SFS plot in dict()
      dict_SFS_allPos,dict_SFS_allIndGT,individuals_number=SFS_profiles_plot.store_SNP_genotypes_file(heterozygosity_file,sep,verbose)

      ### Case for positions where all individuals are genotyped and passed filtering step
      if verbose>=1:
         print('\t\tB/ Plot SFS with genotypes profiles:')
         print('\t\t\ta/ Case for positions where all individuals are genotyped and passed "read coverage+minor allele read frequency" filtering step')
      dict_SFS_allIndGT_profile_occ,dict_SFS_allIndGT_profile_name=SFS_profiles_plot.get_dict_for_SFS_plot_with_genotypes_profile(individuals_number,dict_SFS_allIndGT,max_profiles,False,verbose)
      output_SFS_plot=OUTPUT+'/'+prefix+'_max'+str(max_profiles)+'.pdf'
      SFS_profiles_plot.SFS_plot_genotypes_profiles(individuals_number,dict_SFS_allIndGT_profile_occ,dict_SFS_allIndGT_profile_name,max_profiles,False,output_SFS_plot)
   
      ### Case for all positions (same positions filtered out)
      if verbose>=1:
         print('\t\t\tb/ Case for all positions (same positions "read coverage+minor allele read frequency" filtered out)')
      dict_SFS_allPos_profile_occ,dict_SFS_allPos_profile_name=SFS_profiles_plot.get_dict_for_SFS_plot_with_genotypes_profile(individuals_number,dict_SFS_allPos,max_profiles,True,verbose)
      output_SFS_plot=OUTPUT+'/'+prefix+'_max'+str(max_profiles)+'_all_positions.pdf'
      SFS_profiles_plot.SFS_plot_genotypes_profiles(individuals_number,dict_SFS_allPos_profile_occ,dict_SFS_allPos_profile_name,max_profiles,True,output_SFS_plot)
   #########
   ### For indel
   #########
   elif variant_type=="indel":
      if verbose>=1:
         print('\nVariant type is "indel"')
         print('\t\tA/ Store genotypes infos with filtering TAG')
      ### Set SFS plot in dict()
      dict_SFS_allPos,dict_SFS_allIndGT=SFS_profiles_plot.store_indel_genotypes_file(heterozygosity_file,sep,verbose)

      ### Case for positions where all individuals are genotyped and passed filtering step
      if verbose>=1:
         print('\t\tB/ Plot SFS with genotypes profiles:')
         print('\t\t\ta/ Case for positions where all individuals are genotyped and passed "read coverage+minor allele read frequency" filtering step')
      for indel_size_interval in dict_SFS_allIndGT:
         SFSplot_DIR=OUTPUT+'/'+indel_size_interval
         util.mkdir_p(SFSplot_DIR)
         if verbose>=1:
            print('\t\t\t\tIndel size interval "'+indel_size_interval+'"')
         dict_SFS_allIndGT_profile_occ,dict_SFS_allIndGT_profile_name=SFS_profiles_plot.get_dict_for_SFS_plot_with_genotypes_profile(individuals_number,dict_SFS_allIndGT,max_profiles,False,verbose)
         output_SFS_plot=SFSplot_DIR+'/'+prefix+'_max'+str(max_profiles)+'.pdf'
         SFS_profiles_plot.SFS_plot_genotypes_profiles(individuals_number,dict_SFS_allIndGT_profile_occ,dict_SFS_allIndGT_profile_name,max_profiles,False,output_SFS_plot)

      ### Case for all positions (same positions filtered out)
      if verbose>=1:
         print('\t\t\tb/ Case for all positions (same positions "read coverage+minor allele read frequency" filtered out)')
      for indel_size_interval in dict_SFS_allPos:
         SFSplot_DIR=OUTPUT+'/'+indel_size_interval
         util.mkdir_p(SFSplot_DIR)
         if verbose>=1:
            print('\t\t\t\tIndel size interval "'+indel_size_interval+'"')
         dict_SFS_allPos_profile_occ,dict_SFS_allPos_profile_name=SFS_profiles_plot.get_dict_for_SFS_plot_with_genotypes_profile(individuals_number,dict_SFS_allPos,max_profiles,True,verbose)
         output_SFS_plot=SFSplot_DIR+'/'+prefix+'_max'+str(max_profiles)+'_all_positions.pdf'
         SFS_profiles_plot.SFS_plot_genotypes_profiles(individuals_number,dict_SFS_allPos_profile_occ,dict_SFS_allPos_profile_name,max_profiles,True,output_SFS_plot)
   else:
      exit('\nERROR, variant type (option "-var/--variant_type") should be equal to "SNP" or "indel" and not "'+variant_type+'"!!!\n')


   end_time=datetime.now()
   print('\nDuration: '+str(end_time-start_time))