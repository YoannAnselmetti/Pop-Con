#!/usr/bin/env Rscript
###
###   GOAL: Draw SFS plot with genotype profiles
###      
###   INPUT:
###      1- INPUT file: Genotype profiles file
###         (../example/results/Lineus_longissimus/SNP/MARFt0.0/heterozygosity/read0/genotype_profiles_per_altNb_read0_Lineus_longissimus_SNP_MARFt0.0_max20_with_positions_with_all_individuals_genotyped.tab)
###      2- Maximum number of genotypes profiles displayed in SFS plot
###         (20)
###      3- #individuals in SFS plot
###         (6)
###      4- Title for th SFS plot with genotypes profiles
###         (Lineus_longissimus)
###      5- Path of the SFS plot with the genotypes profiles
###         (../example/results/Lineus_longissimus/SNP/MARFt0.0/SFS_profiles/read0/SFSplot_genotypes_profiles_read0_Lineus_longissimus_SNP_MARFt0.0_max20.pdf)
###      6- Boolean to get the correct -> 1 if lowCov/lowMARF TAG in genotypes profiles, 0 else
###         (0)
###
###   OUTPUT:
###      - SFS plot with genotypes profiles
###
###   Name: SFS_genotypes_profiles_plot.R   Author: Thibault Leroy + Yoann Anselmetti
###   Creation date: 2018/11/13             Last modification: 2019/07/04
###

library(ggplot2)
library(jcolors)
require(grid)


### Get parameters
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=6){
  stop("!!! 6 arguments must be supplied !!!", call.=FALSE)
} else{
  inputFile=args[1]             ### INPUT file path
  max_profiles=strtoi(args[2])  ### Number of the most represented genotype profiles displayed on SFS plot for each pic 
  indNb=strtoi(args[3])         ### Number of individuals in the VCF file
  title=args[4]                 ### Title of the plot
  outputFile=args[5]            ### OUTPUT file path
  low=strtoi(args[6])           ### 0 IF no genotype with lowCov/lowMARF TAG, ELSE 1
}



### Parameters for graph and legend
allTOT=indNb*2
max_lines_per_column=(max_profiles+2)*4


# jcolors palette
colfunc<-jcolors_contin(palette="pal12") # pal12 or rainbow
jcols<-colfunc(max_profiles+1)


### Store genotype profiles in varaible "SFS"
SFS=read.table(inputFile,h=T)


### Get columns and lines number in legend
colNb=1.0
nb_lines=0
for (i in 1:nrow(SFS)){
  if (SFS[i,3]=="01"){ # for each group of combinaisons add separator and text
    SFSI=subset(SFS,alt==SFS[i,1])
    if (length(SFSI$occ[SFSI$occ>0])+1+nb_lines>max_lines_per_column){
      colNb=colNb+1.0
      nb_lines=0
    }
    nb_lines=nb_lines+1
  }
  if (SFS[i,2]>0){
    nb_lines=nb_lines+1
  }
}


### Set SFS plot size + Create SFS plot output file
ratioPointvsInch=0.013837
legendfontsize=12.0
charwidth=ratioPointvsInch*legendfontsize
nbCharperInd=2.0
if (low){
  nbCharperInd=5.0  
}
colwidth=(indNb*nbCharperInd*charwidth)+(3.0*charwidth)
legendWidth=1.0+(colNb*colwidth)
plotWidth=indNb+legendWidth
legendlineheight=12.0*ratioPointvsInch*1.25
plotHeight=legendlineheight*max_lines_per_column+2.0*legendlineheight
pdfSFSplot<-c(outputFile)
pdf(pdfSFSplot, width=plotWidth, height=plotHeight)


### Set SFS barplot (ggplot)
ggplot(data=SFS, aes(x=alt, y=occ, fill=ID)) +
  ggtitle(title) +
  geom_bar(stat="identity", show.legend=FALSE) +
  xlab("#alternative alleles") + ylab("Frequency") +
  scale_x_discrete(name="#alternative alleles", limits=c(1:allTOT)) + theme_bw() +
  theme(axis.line=element_line(colour = 'black', size=1.25), axis.ticks=element_line(colour='black',size=1),
        plot.margin=unit(c(0,legendWidth,0,0),"inches"),
        plot.title=element_text(colour="black", size=16, hjust=.5, face="bold.italic"),
        axis.text.x=element_text(colour="black", size=12, angle=0, hjust=.5, vjust=-.5, face="plain"),
        axis.text.y=element_text(colour="black", size=12, angle=0, hjust=.5, vjust=.5, face="plain"),
        axis.title.x=element_text(colour="black", size=14, angle=0, hjust=.5, vjust=.2, face="italic"),
        axis.title.y=element_text(colour="black", size=14, angle=90, hjust=.5, vjust=.5, face="italic")) +
  scale_fill_manual(name=allTOT/2,labels=SFS$profile[SFS$alt==allTOT/2],values=jcols)


### Print personalized legend
yinit=plotHeight
ypad=legendlineheight
ydec=yinit
xpad=charwidth
xdec=plotWidth-(colNb*colwidth)
column=0.0
nb_lines=0
matrixcolor=rep(jcols,allTOT) # vector of colors

# grid.text("LEGEND:", x=unit(xdec,"inches"), y=unit(yinit,"inches"),just="left",gp=gpar(fontsize=14))

for (i in 1:nrow(SFS)){
  currentcolor=matrixcolor[i]
  if (SFS[i,3]=="01"){ # for each group of combinaisons add separator and text
    SFSI=subset(SFS,alt==SFS[i,1])
    if (length(SFSI$occ[SFSI$occ>0])+1+nb_lines>max_lines_per_column){
      column=column+1.0
      nb_lines=0
      ydec=yinit
    }
    nb_lines=nb_lines+1
    ydec=ydec-ypad
    grid.text(SFS[i,1], x=unit(xdec+(column*colwidth),"inches"), y=unit(ydec,"inches"),just="left",gp=gpar(fontsize=14))
  }
  if (SFS[i,2]>0){
    nb_lines=nb_lines+1
    ydec=ydec-ypad
    # add square
    grid.polygon(x=c(xdec+(column*colwidth), xdec+(column*colwidth), xdec+(column*colwidth)+2.0*charwidth, xdec+(column*colwidth)+2.0*charwidth), y=c(ydec-(legendlineheight/3.0),ydec+(legendlineheight/3.0),ydec+(legendlineheight/3.0),ydec-(legendlineheight/3.0)), default.units="inches", gp=gpar(fill=currentcolor))
    # add text for each square
    grid.text(SFS[i,4], x=unit(xdec+(column*colwidth)+3.0*charwidth,"inches"), y=unit(ydec,"inches"), just="left", gp=gpar(col="darkgrey", fontsize=legendfontsize))
  }
}

dev.off()
