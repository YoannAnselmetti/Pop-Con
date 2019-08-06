Pop-Con: a tool to visualize genotype profiles on Site Frequency Spectrum
=========================================================================

INTRODUCTION
------------
Pop-Con is a tool designed to visualize genotype profiles of a Site Frequency Spectrum (SFS) from a Variant Calling Format (VCF) file [1] for SNP and indel variant.
We define a *genotype profile* as the list of genotypes for all individuals at a given site.
For example, **"0/1,1/1,0/1,0/1"** is the genotype profile corresponding to the VCF line below:
```
Dst_b1v03_scaf00008 54329 . G A 1970.68 PASS AC=6;AF=0.500;AN=12;BaseQRankSum=-6.140e-01;ClippingRankSum=0.00;DP=144;ExcessHet=4.8280;FS=6.723;MLEAC=5;MLEAF=0.417;MQ=60.00;MQRankSum=0.00;QD=15.28;ReadPosRankSum=-8.460e-01;SOR=0.440 GT:AD:DP:GQ:PL 0/1:16,22:38:99:558,0,409 1/1:1,11:12:3:309,3,0 0/1:26,28:54:99:750,0,617 0/1:8,11:19:99:286,0,203
```
Where individuals 1, 3 and 4 are heterozygote ("0/1") and individual 3 is homozygote for the alternative allele ("1/1") for this site.

Pop-Con takes as input a VCF file and uses the **cyvcf2** python package [2] to read and parse the VCF file
 To be readable and parsable by **cyvcf2** the VCf file has to be compressed and indexed with one of the 2 following command lines:
```bash
bgzip input_vcf_file.vcf
tabix -p input_vcf_file.vcf.gz
```
OR
```bash
bcftools view input_vcf_file.vcf -Oz -o input_vcf_file.vcf.gz
bcftools index input_vcf_file.vcf.gz
```

Pop-Con has been developped on VCF files produced with **GATK** [3] and **read2snp** [4] variant callers. Variant caller used to produced VCF file has to be given with option `--tool` (**GATK** is the default variant caller).

`--read` and `--marft` filtering options can be used to test different set of genotype filtering:
* `--read` option take as input a list of read coverage threshold. Genotype with read coverage lower than threshold will be tagged *lowCov*
* `--marft` option take as input a list of minor allele read frequency threshold. Genotype with minor allele read frequency lower than threshold will be tagged *lowMARF*

For each filtering combination ("read coverage" + "minor allele read frequency"), Pop-Con plots 2 figures (for **SNP** and **indel**):
1. with sites where all individuals are genotyped and passed read coverage and minor allele read frequency filtering.
2. with all sites where at least 1 individual is genotyped and passed read coverage and minor allele read frequency filtering.

For figure 1, two plots are produced:
* an upper plot representing the observed proportion of the **X** most represented genotype profiles for each SFS pic (**X** is set with the `-max/--max_profiles` option)
* a lower plot representing for each pic the expected proportion of genotype profile under the Hardy-Weinberg equilibrium. Genotype profiles are colored in red if observed genotype profile is **Y** times upper Hardy-Weinberg expectation, blue if **Y** times lower and white if at equilibrum (between **Y** times upper and **Y** times lower). Where **Y** is set with the option `-fold/--hwe_fold_change`

For figure 2, only SFS plot with observed genotype profiles proportion is produced. This plot can be used to detect abnormal porpotion of low read coverage, low minor allele read frequency or ungenotyped individuals but display is not satisfying for now...

Pop-Con is supported on Linux with python2 (version==2.7) and python3 (versions>=3.4) and surely supported on Macintosh and Windows with the pypi python packages manager but not tested. 


INSTALLATION
------------
## Source from GitHub
This assumes that you have the python modules cyvcf2, numpy and matplotlib installed (cf *requirements.txt*).
```bash
git clone https://github.com/YoannAnselmetti/Pop-Con.git
cd Pop-Con
python ./Pop-Con
```

## Using pip (recommended)
```bash
sudo pip install Pop-Con
```

## Manual installation
In the **Pop-Con** folder where "setup.py" is located, run:
```bash
sudo python setup.py install
```

If you encountered problem during the installation of cyvcf2 or get this following error message when executing Pop-Con:
```
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "/usr/local/lib/python3.6/site-packages/cyvcf2/__init__.py", line 1, in <module>
    from .cyvcf2 import (VCF, Variant, Writer, r_ as r_unphased, par_relatedness,
ImportError: /usr/local/lib/python3.6/site-packages/cyvcf2/cyvcf2.cpython-36m-x86_64-linux-gnu.so: undefined symbol: EVP_sha1
```
You have to install it manually from the source:
https://github.com/brentp/cyvcf2


USAGE
-----
```
usage: Pop-Con [-h] -i VCF_FILE [-t VARIANT_CALLER]
               [-r READ_COVERAGE_THRESHOLD [READ_COVERAGE_THRESHOLD ...]]
               [-m MARFT [MARFT ...]] [-f VARIANT_CALLER_FILTERING]
               [-fmono MONOMORPH_FILTERING] [-p PREFIX] [-v VERBOSE]
               [-o OUTPUT_DIR] [-hf WRITE_HETEROZYGOSITY_FILE] [-sep SEP]
               [-max MAX_PROFILES] [-fold HWE_FOLD_CHANGE]

Pop-Con - A tool for Population genomic Conflicts detection.

optional arguments:
  -h, --help            show this help message and exit
  -i VCF_FILE, --input_file VCF_FILE
                        Variant Calling Format file containing variant calling
                        data.
  -t VARIANT_CALLER, --tool VARIANT_CALLER
                        Variant calling tool used to call variant. Values:
                        "read2snp" or "GATK". (Default: "GATK")
  -r READ_COVERAGE_THRESHOLD [READ_COVERAGE_THRESHOLD ...], --read READ_COVERAGE_THRESHOLD [READ_COVERAGE_THRESHOLD ...]
                        List of read coverage threshold filtering for each
                        genotype. Values range: [0,infinity[. (Default: 0)
  -m MARFT [MARFT ...], --marft MARFT [MARFT ...]
                        List of Minor Allele Read Frequency (MARF) threshold
                        for heterozygous genotypes filtering. Values range:
                        [0.0,0.5]. (Default: 0.0)
  -f VARIANT_CALLER_FILTERING, --variant_caller_filtering VARIANT_CALLER_FILTERING
                        Boolean to set if variant caller filtering is consider
                        or not. (if "True", sites with TAG column are not
                        empty are filtered out from analysis). (Default: True)
  -fmono MONOMORPH_FILTERING, --monomorph_filtering MONOMORPH_FILTERING
                        Boolean to filter out monomorph sites. ("True":
                        monomorph sites are filtered out from analysis. Reduce
                        drastically the execution time!!!). (Default: True)
  -p PREFIX, --prefix PREFIX
                        Experiment name (used as prefix for output files).
                        (Default: "exp1")
  -v VERBOSE, --verbose VERBOSE
                        Verbose intensity. (Default: 1)
  -o OUTPUT_DIR, --output OUTPUT_DIR
                        Output directory path. (Default: ./)
  -hf WRITE_HETEROZYGOSITY_FILE, --heterozygosity_file WRITE_HETEROZYGOSITY_FILE
                        Boolean to set if the heterozygosity files
                        (summarizing VCF file for each combination of read
                        coverage and MARF filtering) have to be written or
                        not. File used by script
                        "scripts/plot_SFS_geneotype_profiles.py" to plot SFS
                        without re-parsing VCF file. (Default: True)
  -sep SEP, --separator SEP
                        Separator used in genotype profiles. (Default: ",")
  -max MAX_PROFILES, --max_profiles MAX_PROFILES
                        Maximum number of genotype profiles displayed in SFS
                        plot. (Default: 10)
  -fold HWE_FOLD_CHANGE, --hwe_fold_change HWE_FOLD_CHANGE
                        Fold change value to define when an observed genotype
                        profile proportion is in excess/deficit compare to the
                        expected value under Hardy-Weinberg Equilibrium.
                        (Default: 2.0)

Source code and manual: http://github.com/YoannAnselmetti/Pop-Con
```

EXAMPLE
-------
Below, the minimal command line to run Pop-Con:
```bash
python Pop-Con -i input_vcf_file.vcf.gz
```
If the VCF file contains SNP and indel variants, Pop-Con will create two directories *"SNP/"* and *"indel/"* in current directory.

Below, an example on VCF file present in directory *"example/data/Lineus_longissimus/"*:
```bash
python Pop-Con -i example/data/Lineus_longissimus/Lineus_longissimus_read10_par0.vcf.gz -p Lineus_longissimus -o example/results/Lineus_longissimus/ -t read2snp -fmono True -max 20
```

Below the architecture of the output directory *"example/results/Lineus_longissimus/"* produced by Pop-Con with command line above:
```
example/results/Lineus_longissimus/
└── SNP
    ├── heterozygosity
    │   └── heterozygosity_allFilter_Lineus_longissimus.tab
    └── MARFt0.0
        ├── heterozygosity
        │	  └── read0
        │		    ├── genotype_profiles_distrib_read0_Lineus_longissimus_SNP_MARFt0.0.tab
        │	      ├── genotype_profiles_per_altNb_read0_Lineus_longissimus_SNP_MARFt0.0_max20_with_all_positions.tab
        │		    ├── genotype_profiles_per_altNb_read0_Lineus_longissimus_SNP_MARFt0.0_max20_with_positions_with_all_individuals_genotyped.tab
        │		    └── heterozygosity_read0_Lineus_longissimus_SNP_MARFt0.0.tab
        └── SFS_profiles
            └── read0
                ├── SFSplot_genotypes_profiles_read0_Lineus_longissimus_SNP_MARFt0.0_max20_all_positions.pdf
                └── SFSplot_genotypes_profiles_read0_Lineus_longissimus_SNP_MARFt0.0_max20.pdf
```

The 2 SFS plot with genotype profiles produced:
### For sites where all individuals are genotyped and passed read coverage + minor allele read frequency filtering
File: "SFSplot_genotypes_profiles_read0_Lineus_longissimus_SNP_MARFt0.0_max20.pdf"
![Lineus longissimus SFS plot with genotype profiles (only individuals genotyped)](docs/img/SFSplot_genotypes_profiles_read0_Lineus_longissimus_SNP_MARFt0.0_max20.pdf)

### For sites where at least 1 individual is genotyped and passed read coverage + minor allele read frequency filtering
File: "SFSplot_genotypes_profiles_read0_Lineus_longissimus_SNP_MARFt0.0_max20_all_positions.pdf"
![Lineus longissimus SFS plot with genotype profiles (all positions)](docs/img/SFSplot_genotypes_profiles_read0_Lineus_longissimus_SNP_MARFt0.0_max20_all_positions.pdf)


BIBLIOGRAPHY
-------------
[1] Danecek, P. et al. The variant call format and VCFtools. Bioinformatics 27, 2156–2158 (2011).
[2] Pedersen, B. S. & Quinlan, A. R. cyvcf2: fast, flexible variant analysis with Python. Bioinformatics 33, 1867–1869 (2017).
[3] McKenna, A. et al. The Genome Analysis Toolkit: A MapReduce framework for analyzing next-generation DNA sequencing data. Genome Res. 20, 1297–1303 (2010).
[4] Gayral, P. et al. Reference-Free Population Genomics from Next-Generation Transcriptome Data and the Vertebrate–Invertebrate Gap. PLOS Genetics 9, e1003457 (2013).


TROUBLESHOOTING
---------------
You get the following error message, when you forgot to set the option `--tool read2snp`. 
```bash
Traceback (most recent call last):
  File "./Pop-Con", line 404, in <module>
    individuals_number,dict_SNP_genotypes,dict_indel_genotypes=read_VCF_file(vcf_file,prefix,OUTPUT,list_read_coverage_threshold,list_MARFt,monomorph_filtering,sep,verbose)
  File "./Pop-Con", line 127, in read_VCF_file
    dict_SNP_genotypes,dict_SNP_SFS_allPos,dict_SNP_SFS_allIndGT,dict_SNP_hetero,dict_FIS=heterozygosity.store_polymorphism(variant,list_Individuals,dict_SNP_genotypes,dict_SNP_SFS_allPos,dict_SNP_SFS_allIndGT,dict_SNP_hetero,dict_FIS,OUTPUT,'SNP',prefix,list_read_coverage_threshold,list_MARFt,False,variant_caller,write_heterozygosity_file,sep,verbose)
  File "/home/yanselmetti/Bureau/Pop-Con/popcon/heterozygosity.py", line 365, in store_polymorphism
    variant_size,x,y,z,listInd,dict_SFS_allPos,dict_SFS_allIndGT=store_SFS(variant,dict_SFS_allPos,dict_SFS_allIndGT,read_coverage_threshold,MARFt,bool_INDEL,variant_caller,verbose)
  File "/home/yanselmetti/Bureau/Pop-Con/popcon/heterozygosity.py", line 212, in store_SFS
    alt,x,y,z,listInd=get_polymorphism_GATK(variant,read_coverage_threshold,MARFt)
  File "/home/yanselmetti/Bureau/Pop-Con/popcon/heterozygosity.py", line 54, in get_polymorphism_GATK
    str_coverage=str(variant).split()[x].split(":")[1]
IndexError: list index out of range
```
This error can also arise if you used an other variant caller tool than *GATK* or *read2snp* as Pop-Con was only tested on VCF files produced by those. If you encountered such problem can you please write an issue, I'll try to shortly produce a new module to parse the VCF produced by the variant caller you used.


MISCELLANEOUS
-------------
For R lovers, there is an extra script (*"scripts/SFS_genotypes_profiles_plot.R"*) to plot SFS with genotype profiles from file.
Example:
```bash
Rscript ./scripts/SFS_genotypes_profiles_plot.R example/results/Lineus_longissimus/SNP/MARFt0.0/heterozygosity/read0/genotype_profiles_per_altNb_read0_Lineus_longissimus_SNP_MARFt0.0_max20_with_positions_with_all_individuals_genotyped.tab 20 6 Lineus_longissimus example/results/Lineus_longissimus/SNP/MARFt0.0/SFS_profiles/read0/SFSplot_genotypes_profiles_read0_Lineus_longissimus_SNP_MARFt0.0_max20_with_Rscript.pdf 0
```
This script only plot SFS with observed porportion of genotype profiles and the one with comparison to the expected values under the Hardy-Weinberg equilibrium:
![Lineus longissimus SFS plot with genotype profiles (only individuals genotyped) with Rscript](docs/img/SFSplot_genotypes_profiles_read0_Lineus_longissimus_SNP_MARFt0.0_max20_with_Rscript.pdf)


#### **Contact**
Yoann Anselmetti : *yoann.anselmetti@umontpellier.fr*