from math import pow,ceil
from collections import OrderedDict,namedtuple   #New in version 2.6

### Matplotlib packages import
import numpy as np 
# import matplotlib
# matplotlib.use('Agg')
import matplotlib.cm as cm
from matplotlib import pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.gridspec import GridSpec
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

### home-made package
from heterozygosity import get_indel_size_interval



############
### DEFINITIONS for genotype profiles
###########
GT=namedtuple("GT",["hetero","homoALT"])
CAT=namedtuple("CAT",["ratio","occ"])
PROFILE=namedtuple('PROFILE',['profile','gt_color','hwe_color'])
OCC=namedtuple('OCC',['obs','hwe','gt_color','hwe_color'])

def Factoriel(N):
   if N>1:
      return N*Factoriel(N-1)
   else:
      return 1


def HW_expectations(ind):
   altNb,dict_HWE_altNb_cat=1,dict()
   while altNb<=(ind*2):
      ### Get "0/1" maximal number for a category with $altNb alternative allele sites and $ind individuals
      x,hetero=0,altNb
      while (hetero+x)>ind:
         hetero-=2
         x+=1

      ### Get occurence of the different category (category: combination of "0/1" and "1/1" for $altNb alternative alleles number and $ind individuals)
      size_tot,dict_GT=0,dict()
      while hetero>=0:
         homoALT=(altNb-hetero)/2
         homoREF=ind-(hetero+homoALT)
         occ_cat,size=0,pow(2,hetero)
         if homoREF:
            occ_cat=Factoriel(ind)/(Factoriel(homoREF)*Factoriel(homoALT)*Factoriel(hetero))
         else:
            occ_cat=Factoriel(ind)/(Factoriel(homoALT)*Factoriel(hetero))
         size_tot+=(occ_cat*size)
         gtALT=GT(hetero,homoALT)
         dict_GT[gtALT]=occ_cat
         hetero-=2

      ### Get the 
      for gt in dict_GT:
         occ_cat=dict_GT[gt]
         hetero=gt.hetero
         size=pow(2,hetero)
         dict_GT[gt]=CAT(float(size)/float(size_tot),occ_cat)
      dict_HWE_altNb_cat[altNb]=dict_GT
      altNb+=1

   return dict_HWE_altNb_cat


def get_alt_nb(genotypes,sep,bool_allPOS):
   altNb=0
   if bool_allPOS or (not bool_allPOS and (not "low" in genotypes) and (not "./." in genotypes)):
      for gt in genotypes.split(sep):
         if not "low" in gt:
            allele1=gt.split("/")[0]
            allele2=gt.split("/")[1]
            if allele1=="1":
               altNb+=1
            if allele2=="1":
               altNb+=1
   return altNb


def store_indel_genotypes_file(GENOTYPES,PROFILES,sep,verbose):
   step,y,x=100000,0,0
   dict_genotypes=dict()
   genotype_file=open(GENOTYPES,"r")
   for line in genotype_file:
      if not "#chromosome/scaffold" in line:
         columns=len(line.split())
         indel_size=int(line.split()[2])
         z,genotypes=0,""
         while 6+z<columns:
            col=line.split()[6+z]
            if not genotypes:
               if "low" in col:
                  genotypes+=(col.split(":")[0]+":"+col.split(":")[1])
               elif ":" in col:
                  genotypes+=col.split(":")[0]
               else:
                  genotypes+=col
            else:
               if "low" in col:
                  genotypes+=(sep+col.split(":")[0]+":"+col.split(":")[1])
               elif ":" in col:
                  genotypes+=(sep+col.split(":")[0])
               else:
                  genotypes+=(sep+col)
            z+=1


         if not 'all_size' in dict_genotypes:
            dict_genotypes['all_size']=dict()
         if not genotypes in dict_genotypes['all_size']:
            dict_genotypes['all_size'][genotypes]=0
         dict_genotypes['all_size'][genotypes]+=1

         indel_size_interval=get_indel_size_interval(indel_size)

         if not indel_size_interval in dict_genotypes:
            dict_genotypes[indel_size_interval]=dict()
         if not genotypes in dict_genotypes[indel_size_interval]:
            dict_genotypes[indel_size_interval][genotypes]=0
         dict_genotypes[indel_size_interval][genotypes]+=1

   output_full_distrib_file=open(PROFILES,"w")
   dict_SFS_allIndGT,dict_SFS_allPos=dict(),dict()
   for indel_size_interval in dict_genotypes:
      dict_SFS_allIndGT[indel_size_interval],dict_SFS_allPos[indel_size_interval]=dict(),dict()
      for gt in OrderedDict(sorted(dict_genotypes[indel_size_interval].items(), key=lambda t: t[1], reverse=True)):
         # Get list of genotype profiles with occurences / SFS value 
         occ=dict_genotypes[indel_size_interval][gt]

         ### Write genotypes profiles distribution (with lowCov, lowMARFt and "./." genotypes)
         output_full_distrib_file.write(str(occ)+"\t"+gt+"\n")

         ### Set SFS dictionary for all positions with at least 1 individual is genotyped
         altNb=get_alt_nb(gt,sep,True)
         if altNb>0:
            if not altNb in dict_SFS_allPos[indel_size_interval]:
               dict_SFS_allPos[indel_size_interval][altNb]=dict()
            if not occ in dict_SFS_allPos[indel_size_interval][altNb]:
               dict_SFS_allPos[indel_size_interval][altNb][occ]=list()
            dict_SFS_allPos[indel_size_interval][altNb][occ].append(gt)

         ### Set SFS dictionary for positions where all individuals are genotyped
         altNb=get_alt_nb(gt,sep,False)
         if altNb>0:
            if not altNb in dict_SFS_allIndGT[indel_size_interval]:
               dict_SFS_allIndGT[indel_size_interval][altNb]=dict()
            if not occ in dict_SFS_allIndGT[indel_size_interval][altNb]:
               dict_SFS_allIndGT[indel_size_interval][altNb][occ]=list()
            dict_SFS_allIndGT[indel_size_interval][altNb][occ].append(gt)

   output_full_distrib_file.close()

   return dict_SFS_allPos,dict_SFS_allIndGT








def store_SNP_genotypes_file(GENOTYPES,PROFILES,sep,verbose):
   step,y,x=100000,0,0
   dict_genotypes=dict()
   genotype_file=open(GENOTYPES,"r")
   for line in genotype_file:
      if not "#chromosome/scaffold" in line:
         columns=len(line.split())
         z,genotypes=0,""
         while 5+z<columns:
            col=line.split()[5+z]
            if not genotypes:
               if "low" in col:
                  genotypes+=(col.split(":")[0]+":"+col.split(":")[1])
               elif ":" in col:
                  genotypes+=col.split(":")[0]
               else:
                  genotypes+=col
            else:
               if "low" in col:
                  genotypes+=(sep+col.split(":")[0]+":"+col.split(":")[1])
               elif ":" in col:
                  genotypes+=(sep+col.split(":")[0])
               else:
                  genotypes+=(sep+col)
            z+=1
         if not genotypes in dict_genotypes:
            dict_genotypes[genotypes]=0
         dict_genotypes[genotypes]+=1

   output_full_distrib_file=open(PROFILES,"w")
   dict_SFS_allIndGT,dict_SFS_allPos=dict(),dict()
   for gt in OrderedDict(sorted(dict_genotypes.items(), key=lambda t: t[1], reverse=True)):
      # Get list of genotype profiles with occurences / SFS value 
      occ=dict_genotypes[gt]

      ### Write genotypes profiles distribution (with lowCov, lowMARFt and "./." genotypes)
      output_full_distrib_file.write(str(occ)+"\t"+gt+"\n")

      ### Set SFS dictionary for all positions with at least 1 individual is genotyped
      altNb=get_alt_nb(gt,sep,True)
      if altNb>0:
         if not altNb in dict_SFS_allPos:
            dict_SFS_allPos[altNb]=dict()
         if not occ in dict_SFS_allPos[altNb]:
            dict_SFS_allPos[altNb][occ]=list()
         dict_SFS_allPos[altNb][occ].append(gt)

      ### Set SFS dictionary for positions where all individuals are genotyped
      altNb=get_alt_nb(gt,sep,False)
      if altNb>0:
         if not altNb in dict_SFS_allIndGT:
            dict_SFS_allIndGT[altNb]=dict()
         if not occ in dict_SFS_allIndGT[altNb]:
            dict_SFS_allIndGT[altNb][occ]=list()
         dict_SFS_allIndGT[altNb][occ].append(gt)

   output_full_distrib_file.close()

   return dict_SFS_allPos,dict_SFS_allIndGT









def HWE_others_occ_color(low,dict_HWE_altNb_cat,altNb,occ_others,occ_tot):
   HWE_others_occ,HWE_others_color=0,"grey"
   if not low:
      ##########
      ### HWE expectations computations
      ##########
      ratio_obs_on_HWE_others_occ=0
      if occ_others!=0:
         HWE_others_occ=0
         for gt in dict_HWE_altNb_cat[altNb]:
            ratio=dict_HWE_altNb_cat[altNb][gt].ratio
            occCAT=dict_HWE_altNb_cat[altNb][gt].occ
            HWE_others_occ+=ratio*float(occ_tot)*occCAT
            # print(ratio,occCAT,HWE_others_occ)

         ratio_obs_on_HWE_others_occ=float(occ_others)/HWE_others_occ
         if ratio_obs_on_HWE_others_occ<=0.5:
            HWE_others_color="blue"
         elif ratio_obs_on_HWE_others_occ>=2.0:
            HWE_others_color="red"
         else:
            HWE_others_color="grey"
   return HWE_others_occ,HWE_others_color



def get_dict_for_SFS_plot_with_genotypes_profile(individuals,dict_SFS,max_profiles,low,verbose):
   colors=cm.Spectral(np.linspace(0,1,max_profiles+2))
   dict_HWE_altNb_cat=HW_expectations(individuals)
   list_=list()
   dict_SFS_profile_occ,dict_SFS_profile_name=dict(),dict()
   for altNb in range(1,(individuals*2)+1): 
      list_profiles_names=list()
      # if verbose>1:
      #    print("\n"+str(altNb)+":")
      if altNb in dict_SFS:
         occ_tot=0
         ### Get total occurence of $altNb pic
         for occ in dict_SFS[altNb]:
            for profile in sorted(dict_SFS[altNb][occ]):
               occ_tot+=occ

         i,stacked_occ=1,0
         for occ in sorted(dict_SFS[altNb], key=int, reverse=True):
            # Several profiles can have the same #occurences
            for profile in sorted(dict_SFS[altNb][occ]):
               if i<=max_profiles:
                  # if verbose>1:
                  #    print("\t"+str(occ)+":\t"+profile)

                  HWE_occ,hwe_color=0,"grey"
                  if not low:
                     ##########
                     ### HWE expectations computations
                     ##########
                     hetero,homoALT=0,0
                     for gt in profile.split(","):
                        if gt=="0/1":
                           hetero+=1
                        elif gt=="1/1":
                           homoALT+=1
                     ### Get ratio of current profile in Hardy-Weinberg Expectations (HWE)
                     ratio=dict_HWE_altNb_cat[altNb][GT(hetero,homoALT)].ratio
                     ### Remove 1 occurence to the category (hetero,homoALT)
                     ### -> Remaining catzegories will be used to compute HWE count for "others" category in the SFS plot with genotype profiles
                     occCAT=dict_HWE_altNb_cat[altNb][GT(hetero,homoALT)].occ
                     dict_HWE_altNb_cat[altNb][GT(hetero,homoALT)]=CAT(ratio,(occCAT-1))

                     HWE_occ=ratio*float(occ_tot)
                     ratio_obs_on_HWE_occ=float(occ)/HWE_occ
                     if ratio_obs_on_HWE_occ<=0.5:
                        hwe_color="blue"
                     elif ratio_obs_on_HWE_occ>=2.0:
                        hwe_color="red"
                     else:
                        hwe_color="grey"

                  ### Fill list $i that will be stacked with the others list
                  if not i in dict_SFS_profile_occ:
                     dict_SFS_profile_occ[i]=list()
                  dict_SFS_profile_occ[i].append(OCC(occ,HWE_occ,colors[i],hwe_color))

                  ### Get list of profiles present in $altNb stack
                  list_profiles_names.append(PROFILE(profile,colors[i],hwe_color))

                  stacked_occ+=occ
                  i+=1
         
         while i<=max_profiles:
            ### Fill list $i that will be stacked with the others list
            if not i in dict_SFS_profile_occ:
               dict_SFS_profile_occ[i]=list()
            dict_SFS_profile_occ[i].append(OCC(0,0,colors[i],"grey"))
            i+=1
         ##############
         ### SET "OTHERS" CATEGROY
         ##############
         occ_others=occ_tot-stacked_occ
         HWE_others_occ,hwe_others_color=HWE_others_occ_color(low,dict_HWE_altNb_cat,altNb,occ_others,occ_tot)
         if not i in dict_SFS_profile_occ:
            dict_SFS_profile_occ[i]=list()
         dict_SFS_profile_occ[i].append(OCC(occ_others,HWE_others_occ,colors[i],hwe_others_color))
         if occ_others:
            ### Get list of profiles present in $altNb stack
            list_profiles_names.append(PROFILE('others',colors[i],hwe_others_color))
      else:
         i=1
         while i<=max_profiles:
            ### Fill list $i that will be stacked with the others list
            if not i in dict_SFS_profile_occ:
               dict_SFS_profile_occ[i]=list()
            dict_SFS_profile_occ[i].append(OCC(0,0,colors[i],"grey"))
            i+=1
         ##############
         ### SET "OTHERS" CATEGORY
         ##############
         if not i in dict_SFS_profile_occ:
            dict_SFS_profile_occ[i]=list()
         dict_SFS_profile_occ[i].append(OCC(0,0,colors[i],"grey"))

      dict_SFS_profile_name[altNb]=list_profiles_names

   return dict_SFS_profile_occ,dict_SFS_profile_name



def get_legend_col_nb(dict_SFS_profile_name,max_lines_per_column,individuals_number):
   colNb,nb_lines,altNb=1.0,0,1
   while altNb<=(individuals_number*2):
      profiles_nb_in_altNb=0
      if altNb in dict_SFS_profile_name:
         profiles_nb_in_altNb=len(dict_SFS_profile_name[altNb])
      if (profiles_nb_in_altNb+nb_lines+1)>max_lines_per_column:
         colNb=colNb+1.0
         nb_lines=0
      nb_lines=nb_lines+profiles_nb_in_altNb+1
      altNb+=1
   return colNb



def print_SFS_gt_profiles_legend(Wrect,dict_SFS_profile_name,ax_legend,yinit,legendfontsize,legendlineheight,max_lines_per_column,charwidth,colwidth,HWE):
   ### Legend positions in columns and raws 
   col,raws=0.0,0
   ### Y axis parameter
   ypad,ydec=legendlineheight,yinit
   ### X axis parameter
   xpad,xdec=charwidth,0

   ### Rectangle dimensions
   Hrect=ypad*2.0/3.0

   for altNb in sorted(dict_SFS_profile_name):
      profiles_nb_in_altNb=len(dict_SFS_profile_name[altNb])
      if (raws+profiles_nb_in_altNb+1)>max_lines_per_column:
         col=col+1.0
         raws=0
         ydec=yinit
         xdec=col*colwidth

      ### Set position on Y axis
      raws=raws+1
      ydec=ydec-(ypad*1.5)

      ### Write #(alternative allele) category
      ax_legend.text(x=xdec,y=ydec,s=str(altNb),fontsize=legendfontsize*1.5)

      for profile_colors in dict_SFS_profile_name[altNb]:
         gt_profile=profile_colors.profile
         color=""
         if HWE:
            color=profile_colors.hwe_color
         else:
            color=profile_colors.gt_color

         ### Set position on Y axis
         raws=raws+1
         ydec=ydec-ypad

         ############
         ### LEGEND
         ############
         # Add a rectangle
         rect=mpatches.Rectangle([xdec,ydec-(.05*ypad)],Wrect,Hrect,color=color)
         ax_legend.add_patch(rect)
         ### Add text for each square
         ax_legend.text(x=xdec+(1.5*Wrect),y=ydec,s=gt_profile,fontsize=legendfontsize)



def SFS_plot_genotypes_profiles(individuals_number,dict_SFS_profile_occ,dict_SFS_profile_name,max_profiles,read_coverage_threshold,MARFt,prefix,low,OUTPUT):
   ### Set legend parameters
   max_lines_per_column=(max_profiles+2)*4
   colNb=get_legend_col_nb(dict_SFS_profile_name,max_lines_per_column,individuals_number)

   legendfontsize=12.0
   ratioPointvsInch=0.013837
   charwidth=legendfontsize*ratioPointvsInch
   indWidth=2
   if low:
      indWidth=4.8
   colwidth=(individuals_number*indWidth*charwidth)+3*charwidth
   legendWidth=(colNb*colwidth)   
   legendlineheight=legendfontsize*ratioPointvsInch*1.25
   
   ### Set plot/figure Height
   plotHeight=int(ceil(legendlineheight*(max_lines_per_column+2)))
   figHeight=2*plotHeight
   if low:
      figHeight=plotHeight

   ### Parameter for legend plot
   colwidth=1.0/float(colNb)
   Wrect=colwidth/(indWidth*(individuals_number)+3)*10.0*legendWidth

   ### Set legend Width
   plotWidth=individuals_number
   legendWidth=int(ceil(legendWidth))
   figWidth=plotWidth+legendWidth

   ### Plot parameters
   width=.75

   ### Intitialize the figure plot
   fig=plt.figure(figsize=(figWidth,figHeight))
   gs=GridSpec(figHeight,figWidth)

   ##########
   ### Axe plot 1
   ##########
   # print("Axe plot 1")
   ax_plot1=fig.add_subplot(gs[:plotHeight-1,:plotWidth-1])
   ax_plot1.clear()
   # print("Axe legend 1")
   ax_legend1=fig.add_subplot(gs[:plotHeight-1,plotWidth-1:])
   ax_legend1.set_axis_off()
   legendWidth=10.0*legendWidth
   ax_legend1.axis([0,legendWidth,0,plotHeight])
   colwidth=legendWidth/float(colNb)
   print_SFS_gt_profiles_legend(Wrect,dict_SFS_profile_name,ax_legend1,plotHeight,legendfontsize,legendlineheight,max_lines_per_column,charwidth,colwidth,False)

   ##########
   ### Axe plot 2
   ##########
   ax_plot2=fig.add_subplot()
   if not low:
      # print("Axe plot 2")
      ax_plot2=fig.add_subplot(gs[plotHeight:,:plotWidth-1])
      ax_plot2.clear()
      # print("Axe legend 2")
      ax_legend2=fig.add_subplot(gs[plotHeight:,plotWidth-1:])
      ax_legend2.set_axis_off()
      ax_legend2.axis([0,legendWidth,0,plotHeight])
      print_SFS_gt_profiles_legend(Wrect,dict_SFS_profile_name,ax_legend2,plotHeight,legendfontsize,legendlineheight,max_lines_per_column,charwidth,colwidth,True)


   ##########
   ### SFS plot 1 & 2
   ##########
   N=individuals_number*2
   index=np.arange(N)
   # print(dict_SFS_profile_occ)
   bottom_ax_plot1,bottom_ax_plot2=list(),list()
   for stack in sorted(dict_SFS_profile_occ,reverse=True):
      stack_list_gt,stack_list_gt_color,stack_list_hwe,stack_list_hwe_color=list(),list(),list(),list()
      for occ in dict_SFS_profile_occ[stack]:
         stack_list_gt.append(occ.obs)
         stack_list_gt_color.append(occ.gt_color)
         stack_list_hwe.append(occ.hwe)
         stack_list_hwe_color.append(occ.hwe_color)

      ### Plot current stack
      sidecolor,Lwidth='black',.1
      if bottom_ax_plot1:
         ax_plot1.bar(index,stack_list_gt,width=width,bottom=bottom_ax_plot1,color=stack_list_gt_color,edgecolor=sidecolor,linewidth=Lwidth)
         new_bottom_ax_plot1=list()
         for x in zip(bottom_ax_plot1,stack_list_gt):
            # print(sum(x))
            new_bottom_ax_plot1.append(sum(x))
         bottom_ax_plot1=new_bottom_ax_plot1
      else:
         bottom_ax_plot1=stack_list_gt
         ax_plot1.bar(index,stack_list_gt,width=width,color=stack_list_gt_color,edgecolor=sidecolor,linewidth=Lwidth)

      if not low:
         if bottom_ax_plot2:
            ax_plot2.bar(index,stack_list_hwe,width=width,bottom=bottom_ax_plot2,color=stack_list_hwe_color,edgecolor=sidecolor,linewidth=Lwidth)
            new_bottom_ax_plot2=list()
            for x in zip(bottom_ax_plot2,stack_list_hwe):
               # print(sum(x))
               new_bottom_ax_plot2.append(sum(x))
            bottom_ax_plot2=new_bottom_ax_plot2
         else:
            ax_plot2.bar(index,stack_list_hwe,width=width,color=stack_list_hwe_color,edgecolor=sidecolor,linewidth=Lwidth)
            bottom_ax_plot2=stack_list_hwe

   ### ax_plot1 title
   ax_plot1.set_title('Observed',fontstyle='italic',fontweight='bold')
   ### X axis
   ax_plot1.set_xlabel("#alternative allele")
   ax_plot1.set_xticks(index)
   ax_plot1.set_xticklabels([str(altNb+1) for altNb in index])
   ### Y axis
   ax_plot1.set_ylabel("Frequency")
   # ax_plot1.set_yticks(fontsize=24)
   ax_plot1.set_ylim(ymin=0)

   if not low:
      ### ax_plot2 title
      ax_plot2.set_title('Expected\n(Hardy-Weinberg Equilibrium)',fontstyle='italic',fontweight='bold')
      ### X axis
      ax_plot2.set_xlabel("#alternative allele")
      ax_plot2.set_xticks(index)
      ax_plot2.set_xticklabels([str(altNb+1) for altNb in index])
      ### Y axis
      ax_plot2.set_ylabel("Frequency")
      # ax_plot1.set_yticks(fontsize=24)
      ax_plot2.set_ylim(ymin=0)
      ### Hardy-Weinberg Equilibrium comparison legend 
      rect1=Line2D([],[],marker="s",markersize=10,linewidth=0,color="red")
      rect2=Line2D([],[],marker="s",markersize=10,linewidth=0,color="grey")
      rect3=Line2D([],[],marker="s",markersize=10,linewidth=0,color="blue")
      ax_plot2.legend([rect1,rect2,rect3],['Excess','Equilibrium','Deficit'])

   ### Write output plot
   plt.tight_layout()
   plt.savefig(OUTPUT,format='pdf')
   plt.close()








#########
### def to write SFS genotypes profiles file readable with Rscript (scripts/SFs)
#########
def write_profile_ID(i):
   if i<10:
      return "0"+str(i)
   else:
      if i<100:
         return str(i)
      else:
         exit("ERROR: maximum number of genotype profiles displayed should be < 100")


def write_R_SFS_profiles(individuals,dict_SFS,OUTPUT,max_profiles,verbose):
   dict_HWE_altNb_cat=HW_expectations(individuals)
   output_distrib_per_altNb_file=open(OUTPUT,"w")
   output_distrib_per_altNb_file.write("alt\tocc\tID\tprofile\n")
   for altNb in range(1,(individuals*2)+1): 
      if altNb in dict_SFS:
         occ_tot=0
         ### Get total occurence of $altNb pic
         for occ in dict_SFS[altNb]:
            for profile in sorted(dict_SFS[altNb][occ]):
               occ_tot+=occ
         i,stacked_occ=1,0
         for occ in sorted(dict_SFS[altNb], key=int, reverse=True):
            # Several profiles can have the same #occurences
            for profile in sorted(dict_SFS[altNb][occ]):
               if i<=max_profiles:
                  output_distrib_per_altNb_file.write(str(altNb)+"\t"+str(occ)+"\t"+write_profile_ID(i)+"\t"+profile+"\n")
                  stacked_occ+=occ
                  i+=1
         while i<=max_profiles:
            output_distrib_per_altNb_file.write(str(altNb)+"\t0\t"+write_profile_ID(i)+"\tN/A\n")
            i+=1
         output_distrib_per_altNb_file.write(str(altNb)+"\t"+str(occ_tot-stacked_occ)+"\tothers\tothers\n")

      else:
         i=1
         while i<=max_profiles:
            output_distrib_per_altNb_file.write(str(altNb)+"\t0\t"+write_profile_ID(i)+"\tN/A\n")
            i+=1
         output_distrib_per_altNb_file.write(str(altNb)+"\t0\tothers\tothers\n")
   output_distrib_per_altNb_file.close()