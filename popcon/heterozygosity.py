from collections import namedtuple   #New in version 2.6
from os import path
from decimal import Decimal
from util import mkdir_p

GENOTYPE=namedtuple('GENOTYPE',['x','y','z','x_allIndGT','y_allIndGT','z_allIndGT','mono','snp','lowCov','pos_allIndGT'])
FIS=namedtuple('FIS',['nfobs','nfexp','fobs','fexp'])


### Courtesy to Brent Pedersen (https://github.com/brentp/cyvcf2/issues/58)
### Class Genotype for a more ergonomic view of genotypes parsed by cyvcf2
class Genotype(object):
    __slots__ = ('alleles','phased')

    def __init__(self, li):
        self.alleles = li[:-1]
        self.phased = li[-1]

    def __str__(self):
        sep = "/|"[int(self.phased)]
        return sep.join("0123."[a] for a in self.alleles)
    __repr__ = __str__



def get_polymorphism_GATK(variant,cov_threshold,minor_allele_frequency_threshold):
   # Get list of read coverage of the site for allele1 (list_cov1) and allele2 (list_cov2)
   list_cov1,list_cov2=list(),list()
   for ref_depth in variant.gt_ref_depths:
      list_cov1.append(ref_depth)
   for alt_depth in variant.gt_alt_depths:
      list_cov2.append(alt_depth)

   # Get read coverage of the genotypes
   list_cov=list()
   for depth in variant.gt_depths:
      list_cov.append(depth)

   # Get individuals genotypes of the site 
   geno=variant.genotypes
   genotypes=[Genotype(individual) for individual in geno]
   ind=len(geno)

   # x:0/0, y:0/1 and z:1/1
   i,alt,x,y,z=0,0,0,0,0
   all_individuals=True
   listInd=list()
   # For each genotype get allele frequency
   for gt in genotypes:
      # Get read coverage of each genotype
      tot_cov=int(list_cov[i])
      if tot_cov<0:
         x=i-ind
         str_coverage=str(variant).split()[x].split(":")[1]
         if str_coverage==".":
            tot_cov=0
         elif str_coverage=="x":
            tot_cov=float(variant.INFO.get('DP')/ind)
         else:
            tot_cov=int(str_coverage)

      cov1=float(list_cov1[i])
      cov2=float(list_cov2[i])

      # Keep only genotype for which read coverage is upper a given threshold
      allele1=str(gt).split("/")[0]
      allele2=str(gt).split("/")[1]
      if tot_cov>=cov_threshold:
         if allele1==allele2=="0":
            listInd.append(str(gt)+":"+str(tot_cov))
            x+=1
         elif allele1==allele2=="1":
            listInd.append(str(gt)+":"+str(tot_cov))
            alt+=2
            z+=1
         elif (allele1=="0" and allele2=="1") or (allele1=="1" and allele2=="0"):
            if cov1>0.0 and cov2>0.0:
               MARF=min(round(cov1/float(tot_cov),2),round(cov2/float(tot_cov),2))
               if MARF>=minor_allele_frequency_threshold:
                  listInd.append("0/1:"+str(tot_cov))
                  alt+=1
                  y+=1
               else:
                  listInd.append("lowMARF:0/1:"+str(MARF))
                  all_individuals=False
            else:
               listInd.append("./.:"+str(tot_cov))
               all_individuals=False
         # If allele is null
         elif allele1==allele2==".":
            listInd.append("./.:"+str(tot_cov))
            all_individuals=False
         else:
            exit("ERROR, allele should be equal to \"0\" or \"1\" and not "+str(allele1)+" or "+str(allele2)+"!!!")
      else:
         all_individuals=False
         if allele1==allele2==".":
            listInd.append("./.:"+str(tot_cov))
         else:
            listInd.append("lowCov:"+str(gt)+":"+str(tot_cov))
      i+=1

   if not all_individuals:
      alt=-alt
   
   return alt,x,y,z,listInd



def get_polymorphism_read2snp(variant):
   # Get individuals genotypes of the site 
   geno=variant.genotypes
   genotypes=[Genotype(ind) for ind in geno]

   # x:0/0, y:0/1 and z:1/1
   alt,x,y,z=0,0,0,0
   all_individuals=True
   listInd=list()
   # For each genotype get allele frequency
   for gt in genotypes:
      allele1=str(gt).split("|")[0]
      allele2=str(gt).split("|")[1]

      if allele1==allele2=="0":
         listInd.append("0/0")
         x+=1
      elif allele1==allele2=="1":
         listInd.append("1/1")
         alt+=2
         z+=1
      elif (allele1=="0" and allele2=="1") or (allele1=="1" and allele2=="0"):
         listInd.append("0/1")
         alt+=1
         y+=1
      # If allele is null
      elif allele1==allele2==".":
         listInd.append("./.")
         all_individuals=False
      else:
         exit("ERROR, allele should be equal to \"0\", \"1\" or \".\" and not "+str(allele1)+" or "+str(allele2)+"!!!")

   if not all_individuals:
      alt=-alt

   return alt,x,y,z,listInd



def store_SFS(variant,dict_SFS_allPos,dict_SFS_allIndGT,read_coverage_threshold,MARFt,bool_INDEL,variant_caller,verbose):
   ##########
   ### Get INDELs polymorphism infos and store them in dict(): dict_SFS_allIndGT & dict_SFS_allPos 
   ##########
   variant_size,alt,x,y,z,listInd=1,0,0,0,0,list()
   if bool_INDEL:
      ##############
      ### Get size of the indel
      ##############
      shortest,longest=0,0
      # if variant_caller=='read2snp':
      ### Get shortest sequence of current indel
      if not variant.ALT:
         shortest=len(variant.REF)
      elif variant.ALT[0]=='*' or variant.REF=='*':
         shortest=0
      else:
         shortest=min(len(variant.REF),len(variant.ALT[0]))
      ### Get longest sequence of current indel
      if not variant.ALT:
         longest=len(variant.REF)
      else:
         longest=max(len(variant.REF),len(variant.ALT[0]))
      # elif variant_caller=='GATK':
      #    shortest=min(len(variant.REF),len(variant.ALT))
      #    longest=max(len(variant.REF),len(variant.ALT))

      ### Set current indel size
      variant_size=longest-shortest

      if verbose>2:
         print('\t\t\t\t-> Indel size:',variant_size,'| variant.ALT:',variant.ALT,'| variant.REF:',variant.REF)

      if not variant_size in dict_SFS_allIndGT[read_coverage_threshold][MARFt]:
         dict_SFS_allIndGT[read_coverage_threshold][MARFt][variant_size]=dict()
         dict_SFS_allPos[read_coverage_threshold][MARFt][variant_size]=dict()

      ###############
      ### Get polymorphism infos and store them in dict()
      ###############
      if variant_caller=='read2snp':
         alt,x,y,z,listInd=get_polymorphism_read2snp(variant)
      elif variant_caller=='GATK':
         alt,x,y,z,listInd=get_polymorphism_GATK(variant,read_coverage_threshold,MARFt)
      if alt>0:# All individuals are genotyped
         if not alt in dict_SFS_allIndGT[read_coverage_threshold][MARFt][variant_size]:
            dict_SFS_allIndGT[read_coverage_threshold][MARFt][variant_size][alt]=0
         dict_SFS_allIndGT[read_coverage_threshold][MARFt][variant_size][alt]+=1
         if not alt in dict_SFS_allPos[read_coverage_threshold][MARFt][variant_size]:
            dict_SFS_allPos[read_coverage_threshold][MARFt][variant_size][alt]=0
         dict_SFS_allPos[read_coverage_threshold][MARFt][variant_size][alt]+=1
      elif alt<0:# NOT all individuals are genotyped
         if not -alt in dict_SFS_allPos[read_coverage_threshold][MARFt][variant_size]:
            dict_SFS_allPos[read_coverage_threshold][MARFt][variant_size][-alt]=0
         dict_SFS_allPos[read_coverage_threshold][MARFt][variant_size][-alt]+=1

   ##########
   ### Get SNPs polymorphism infos and store them in dict(): dict_SFS_allIndGT & dict_SFS_allPos 
   ##########
   else:
      if variant_caller=='read2snp':
         alt,x,y,z,listInd=get_polymorphism_read2snp(variant)
      elif variant_caller=='GATK':
         alt,x,y,z,listInd=get_polymorphism_GATK(variant,read_coverage_threshold,MARFt)
      # print '\t\t\t\t',alt,x,y,z,listInd
      if alt>0:# All individuals are genotyped
         if not alt in dict_SFS_allIndGT[read_coverage_threshold][MARFt]:
            dict_SFS_allIndGT[read_coverage_threshold][MARFt][alt]=0
         dict_SFS_allIndGT[read_coverage_threshold][MARFt][alt]+=1
         if not alt in dict_SFS_allPos[read_coverage_threshold][MARFt]:
            dict_SFS_allPos[read_coverage_threshold][MARFt][alt]=0
         dict_SFS_allPos[read_coverage_threshold][MARFt][alt]+=1
      elif alt<0:# NOT all individuals are genotyped
         if not -alt in dict_SFS_allPos[read_coverage_threshold][MARFt]:
            dict_SFS_allPos[read_coverage_threshold][MARFt][-alt]=0
         dict_SFS_allPos[read_coverage_threshold][MARFt][-alt]+=1

   return variant_size,x,y,z,listInd,dict_SFS_allPos,dict_SFS_allIndGT



def get_morphism(x,y,z,list_Individuals,verbose):
   x_allIndGT,y_allIndGT,z_allIndGT,mono,snp,lowCov,pos_allIndGT=0,0,0,0,0,0,0
   if x or y or z:
      if int(y+z+x)==len(list_Individuals):
         pos_allIndGT+=1
         x_allIndGT=x
         y_allIndGT=y
         z_allIndGT=z
      if y==0 and z==0:
         mono+=1
         if verbose>3:
            if read_coverage_threshold==0:
               print('WARNING MONO!!!')
      # Polymorphism site (SNP or indel)
      else:
         snp+=1 
   else:
      lowCov+=1

   return x_allIndGT,y_allIndGT,z_allIndGT,mono,snp,lowCov,pos_allIndGT



def write_heterozygosity_SNP(output_file,variant,x,y,z,listInd,list_Individuals):
   if not path.exists(output_file):
      output=open(output_file,'w')
      output.write('#chromosome/scaffold\tposition\t0/0\t0/1\t1/1')
      for ind in list_Individuals:
         output.write('\t'+str(ind))
      output.write('\n')
      output.write(variant.CHROM+'\t'+str(variant.start+1)+'\t'+str(x)+'\t'+str(y)+'\t'+str(z)+'\t'+'\t'.join(listInd)+'\n')
      output.close()
   else:
      output=open(output_file,'a')
      output.write(variant.CHROM+'\t'+str(variant.start+1)+'\t'+str(x)+'\t'+str(y)+'\t'+str(z)+'\t'+'\t'.join(listInd)+'\n')
      output.close()



def write_heterozygosity_indel(output_file,indel_size,variant,x,y,z,listInd,list_Individuals):
   if not path.exists(output_file):
      output=open(output_file,'w')
      output.write('#chromosome/scaffold\tposition\tindel_size\t0/0\t0/1\t1/1')
      for ind in list_Individuals:
         output.write('\t'+str(ind))
      output.write('\n')
      output.write(variant.CHROM+'\t'+str(variant.start+1)+'\t'+str(indel_size)+'\t'+str(x)+'\t'+str(y)+'\t'+str(z)+'\t'+'\t'.join(listInd)+'\n')
      output.close()
   else:
      output=open(output_file,'a')
      output.write(variant.CHROM+'\t'+str(variant.start+1)+'\t'+str(indel_size)+'\t'+str(x)+'\t'+str(y)+'\t'+str(z)+'\t'+'\t'.join(listInd)+'\n')
      output.close()


def get_indel_size_interval(indel_size):
   indel_size_interval=''
   if indel_size<10:
      indel_size_interval='[1,10['
   ### Set indel SFS dict() for indels with size>=25bp
   elif indel_size>=25:
      indel_size_interval='[25,infty['
   ### Set indel SFS dict() for 10bp<indel size<=25bp
   else:
      indel_size_interval='[10,25['

   return indel_size_interval


def get_gt_profile(listInd,sep):
   gt_profile=''
   for gt in listInd:
      if not gt_profile:
         if "low" in gt:
            gt_profile+=(gt.split(":")[0]+":"+gt.split(":")[1])
         elif ":" in gt:
            gt_profile+=gt.split(":")[0]
         else:
            gt_profile+=gt
      else:
         if "low" in gt:
            gt_profile+=(sep+gt.split(":")[0]+":"+gt.split(":")[1])
         elif ":" in gt:
            gt_profile+=(sep+gt.split(":")[0])
         else:
            gt_profile+=(sep+gt)

   return gt_profile



def fill_dict_genotypes(dict_genotypes,read_coverage_threshold,MARFt,bool_INDEL,variant_size,listInd,sep):
   gt_profile=get_gt_profile(listInd,sep)

   if not read_coverage_threshold in dict_genotypes:
      dict_genotypes[read_coverage_threshold]=dict()
   if not MARFt in dict_genotypes[read_coverage_threshold]:
      dict_genotypes[read_coverage_threshold][MARFt]=dict()

   if bool_INDEL:
      indel_size_interval=get_indel_size_interval(variant_size)
      if not indel_size_interval in dict_genotypes[read_coverage_threshold][MARFt]:
         dict_genotypes[read_coverage_threshold][MARFt][indel_size_interval]=dict()
      if not gt_profile in dict_genotypes[read_coverage_threshold][MARFt][indel_size_interval]:
         dict_genotypes[read_coverage_threshold][MARFt][indel_size_interval][gt_profile]=0
      dict_genotypes[read_coverage_threshold][MARFt][indel_size_interval][gt_profile]+=1
   else:
      if not gt_profile in dict_genotypes[read_coverage_threshold][MARFt]:
         dict_genotypes[read_coverage_threshold][MARFt][gt_profile]=0
      dict_genotypes[read_coverage_threshold][MARFt][gt_profile]+=1

   return dict_genotypes



def store_polymorphism(variant,list_Individuals,dict_genotypes,dict_SFS_allPos,dict_SFS_allIndGT,dict_hetero,dict_FIS,OUTPUT,variantType,prefix,list_read_coverage_threshold,list_MARFt,bool_INDEL,variant_caller,write_heterozygosity_file,sep,verbose):
   for read_coverage_threshold in list_read_coverage_threshold:
      ########
      ### 
      ########
      if not read_coverage_threshold in dict_SFS_allIndGT:
         dict_SFS_allIndGT[read_coverage_threshold]=dict()
         dict_SFS_allPos[read_coverage_threshold]=dict()

      #############
      ### Set SFS & heterozygosity for the different MARFt (Minor Allele Frequency threshold) values 
      #############
      for MARFt in list_MARFt:
         ### Create MARFt key for dict(): dict_SFS_allIndGT & dict_SFS_allPos 
         if not MARFt in dict_SFS_allIndGT[read_coverage_threshold]:
            dict_SFS_allIndGT[read_coverage_threshold][MARFt]=dict()
            dict_SFS_allPos[read_coverage_threshold][MARFt]=dict()

         ###########
         ### Get SFS infos for current variant and add it to dict(): 'dict_SFS_allIndGT' & 'dict_SFS_allPos'
         ###########
         variant_size,x,y,z,listInd,dict_SFS_allPos,dict_SFS_allIndGT=store_SFS(variant,dict_SFS_allPos,dict_SFS_allIndGT,read_coverage_threshold,MARFt,bool_INDEL,variant_caller,verbose)

         ###########
         ### Get type of current variant site
         ###########
         x_allIndGT,y_allIndGT,z_allIndGT,mono,snp,lowCov,pos_allIndGT=get_morphism(x,y,z,list_Individuals,verbose)

         dict_genotypes=fill_dict_genotypes(dict_genotypes,read_coverage_threshold,MARFt,bool_INDEL,variant_size,listInd,sep)

         ###########
         ### Write genotypes of each position
         ###########
         if write_heterozygosity_file:
            # Write current genotype in 'output_file'
            output_dir=OUTPUT+'/'+variantType+'/MARFt'+str(MARFt)+'/heterozygosity/read'+str(read_coverage_threshold)
            mkdir_p(output_dir)
            output_file=output_dir+'/heterozygosity_read'+str(read_coverage_threshold)+'_'+prefix+'_'+variantType+'_MARFt'+str(MARFt)+'.tab'
            if bool_INDEL:
               write_heterozygosity_indel(output_file,variant_size,variant,x,y,z,listInd,list_Individuals)
            else:
               write_heterozygosity_SNP(output_file,variant,x,y,z,listInd,list_Individuals)

         ###########
         ### Set F-stats values in dict() dict_FIS
         ###########
         x_allIndGT,x,y_allIndGT,y,z_allIndGT,z=Decimal(x_allIndGT),Decimal(x),Decimal(y_allIndGT),Decimal(y),Decimal(z_allIndGT),Decimal(z)
         hetero_obs,hetero_exp,hetero_obs_allIndGT,hetero_exp_allIndGT=Decimal(0.0),Decimal(0.0),Decimal(0.0),Decimal(0.0)
         # Set the value of dict_FIS to compute F-stat:
         if read_coverage_threshold in dict_FIS:
            if MARFt in dict_FIS[read_coverage_threshold]:
               hetero_obs=dict_FIS[read_coverage_threshold][MARFt].nfobs
               hetero_exp=dict_FIS[read_coverage_threshold][MARFt].nfexp
               hetero_obs_allIndGT=dict_FIS[read_coverage_threshold][MARFt].fobs
               hetero_exp_allIndGT=dict_FIS[read_coverage_threshold][MARFt].fexp
         else:
            dict_FIS[read_coverage_threshold]=dict()

         if x+y+z>0:
            if x_allIndGT+y_allIndGT+z_allIndGT>0:
               dict_FIS[read_coverage_threshold][MARFt]=FIS(hetero_obs+(y/(x+y+z)),hetero_exp+(Decimal(2.0)*((Decimal(2.0)*x+y)/(Decimal(2.0)*(x+y+z)))*((Decimal(2.0)*z+y)/(Decimal(2.0)*(x+y+z)))),hetero_obs_allIndGT+(y_allIndGT/(x_allIndGT+y_allIndGT+z_allIndGT)),hetero_exp_allIndGT+(Decimal(2.0)*((Decimal(2.0)*x_allIndGT+y_allIndGT)/(Decimal(2.0)*(x_allIndGT+y_allIndGT+z_allIndGT)))*((Decimal(2.0)*z_allIndGT+y_allIndGT)/(Decimal(2.0)*(x_allIndGT+y_allIndGT+z_allIndGT)))))
            else:
               dict_FIS[read_coverage_threshold][MARFt]=FIS(hetero_obs+(y/(x+y+z)),hetero_exp+(Decimal(2.0)*((Decimal(2.0)*x+y)/(Decimal(2.0)*(x+y+z)))*((Decimal(2.0)*z+y)/(Decimal(2.0)*(x+y+z)))),hetero_obs_allIndGT,hetero_exp_allIndGT)
         else:
            dict_FIS[read_coverage_threshold][MARFt]=FIS(hetero_obs,hetero_exp,hetero_obs_allIndGT,hetero_exp_allIndGT)

         ###########
         ### Set heterozygosity values in dict() dict_hetero
         ###########
         if read_coverage_threshold in dict_hetero:
            if MARFt in dict_hetero[read_coverage_threshold]:
               x+=dict_hetero[read_coverage_threshold][MARFt].x
               y+=dict_hetero[read_coverage_threshold][MARFt].y
               z+=dict_hetero[read_coverage_threshold][MARFt].z
               x_allIndGT+=dict_hetero[read_coverage_threshold][MARFt].x_allIndGT
               y_allIndGT+=dict_hetero[read_coverage_threshold][MARFt].y_allIndGT
               z_allIndGT+=dict_hetero[read_coverage_threshold][MARFt].z_allIndGT
               mono+=dict_hetero[read_coverage_threshold][MARFt].mono
               snp+=dict_hetero[read_coverage_threshold][MARFt].snp
               lowCov+=dict_hetero[read_coverage_threshold][MARFt].lowCov
               pos_allIndGT+=dict_hetero[read_coverage_threshold][MARFt].pos_allIndGT
         else:
            dict_hetero[read_coverage_threshold]=dict()

         dict_hetero[read_coverage_threshold][MARFt]=GENOTYPE(x,y,z,x_allIndGT,y_allIndGT,z_allIndGT,mono,snp,lowCov,pos_allIndGT)

   return dict_genotypes,dict_SFS_allPos,dict_SFS_allIndGT,dict_hetero,dict_FIS
