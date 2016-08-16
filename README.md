# cisASE
A likelihood-based method for detecting putative cis-regulated allele-specific expression in RNA sequencing data

Version 1.0.2

Updated date: 2016.08.16

#####
## Author and License

Author: Zhi Liu

Email: liuzhi_2015@outlook.com, liuzhi2011@sibs.ac.cn

Licensed under the GNU Affero General Public License version 3 or later

#####

##DEPENDENCIES

perl 5.10

R 3.2.2

perl modules  Statistics::R, Statistics::Descriptive

Currently Unix/Linux is supported 
#####
##USAGE

To get the usage of cisASE, run the script

######perl cisASE.1.0.2.pl

To run examples, use the test script in the examples directory:

######cd examples
######sh run.test.sh

###Input prepare
#####I. SNV list for analysis
Data should be a table separated txt consisting of four columns, which are chromosome, position, reference base and alternative base. 
 
#####II. pileup files for RNA-seq and DNA-seq for analysis
######samtools mpileup -C50 –B –f ref.fa –L SNVlist.txt sample.bam > cell.chr1.pileup
pileup file of DNA-seq is optional but recommended. 

#####III. Phasing file
Phasing file is optional for exon and gene level ASE detection.
File should be a table separated txt consisting of three columns, which are SNP id, base of haplotype1 ad base of haplotype2. 

#####VI. Annotation file
Annotation file should be contained in exon and gene level detection. Human and mouse annotation file can be downloaded from this page, or derived from Ensemble biomart. Or users can generate their own annotation file containing the following columns,

Column 1, chromosome(without prefix "chr")

Column 2, gene id

Column 3, exon id

Column 4, constitutive status, that is whether the exon is a constitutive exon (1) of the gene or not (0)

Column 5, exon start

Column 6, exon end

###Command
######perl cisASE.pl  [SNV|EXON|GENE]  - X RNA.pileup -L SNV.list -f * -o outfile [options]*

#####Required arguments
-X	RNA pileup file

-L	List of SNVs for ASE dection 

-f	RNA fastq encoding format. 64=Solexa,Illumina 1.3+,Illumina 1.5+; 33=Sanger,Illumina 1.8+

-o	output file 

#####Optional arguments
-Y	DNA pileup file

-F	DNA fastq encoding format. 64=Solexa,Illumina 1.3+,Illumina 1.5+; 33=Sanger,Illumina 1.8+(required, if -Y is specified)

-B	DNA bias

-M	Method. 1 for chi-square test; 2 for cisASE; 3 for both. Default:2

-D	minimum DNA read depth. Default:20

-R	minimum DNA Reference Allele counts. Default:5

-A	minimum DNA Alternative Allele counts Default:5

-d	minimum RNA read depth Default:20 

-r	minimum RNA Reference Allele counts Default:0

-a	minimum RNA Alternative Allele counts Default:0

-S	times for simulation. Default: 2000

-p	phasing file. Only be available for EXON and GENE level. 

-N	Annotation file. Only be available for EXON and GENE level. 

-C	Consider exon statue (1) or not (0); Default: 0. Specific for GENE level. 


#####-X/-Y/-B
DNA pileup file is the same as that of RNA as described before. Though DNA pileup file is strongly recommended, cisASE can also deal with situation when it is not available, in this case and no specified parameters for option –B, cisASE calculates the mean allele ratio from RNA data, and takes it as DNA ratio. As we presented in the paper, the accuracy will be much higher when adding DNA information in ASE detection. And when no DNA data is available, setting DNA bias –B as 0.5 is a better choice than leave it as default (substituted by RNA ratio).

When some SNVs in SNV list are not included DNA pileup files, cisASE automatically assigns them a DNA ratio as the mean DNA ratio of all the other SNVs having DNA information. As the results presented in the paper, setting DNA bias as mean DNA ratio performed slightly better than setting it as 0.5. For the SNV level detection, if users don’t want to take these SNVs into further consideration, they can be removed from the output file. But for the exon and gene level ASE detection, we don’t removed them by default, if you are sure using mean ratio substitution for these SNVs are questionable, please filter them out from the SNV list before processing them to cisASE. 

#####-f/-F
For the base quality encoding format, if you run FastQC on your file, they will detect the encoding automatically. Then you can assign a corresponding number to –f and –F according to the reported encoding format. We also provide a script GetFastqQualityEncoding.pl for this purpose. 

######perl GetFastqQualityEncoding.pl input.fastq

#####-S
Times of simulation to decide a decision threshold for likelihood test. In all of our analysis, we set the simulation times as 2000. When parameter –S is not specified, no simulation will be performed, but we have given a recommended LOD at significant level 0.05 and 0.01 based on the application on our real dataset. However, please notice that since the LOD threshold is generated based on real dataset, variance exists for dataset of different quality and coverage, if computational ability allowed, please always set –S, or at least performing simulation on several samples out of your whole panel of dataset.

###Quick Examples
#####SNV level

Parallel DNA pileup file is strongly recommended for each level.

######perl cisASE.pl SNV -L snv.list -X RNA.pileup -Y DNA.pileup -F 64 -f 33 -o outfile

If DNA data is not available, users can predefine DNA bias by setting -B or leave it as default.

######perl cisASE.pl SNV -L snv.list -X RNA.pileup -B 0.50 -f 33 -o outfile

Both likelihood-based method and chi-square test are provided in cisASE. Set -M to select.

#####EXON level
For EXON-level and GENE-level ASE detection, option for annotation file -N should be specified.

######perl cisASE.pl EXON -L snv.list -X RNA.pileup -Y DNA.pileup -F 64 -f 33 -N annotation file -o outfile 

Phasing file is optional for EXON-level and GENE-level ASE detection. If no phasing file is specified, cisASE will perform a pseudo-phasing as described in our paper.

######perl cisASE.pl EXON -L snv.list -X RNA.pileup -Y DNA.pileup -P phased.txt -F 64 -f 33 -N annotation file -o outfile

#####GENE level
GENE-level ASE detection can accept another option -C, to specify whether to consider the constitutive statue of an exon when combining exons to a gene. cisASE will take account  constitutive exons only, and calculate the difference between each alternative exons and the gene when setting -C 1, and take all exons into account when the default setting(-C 0) is applied.

######perl cisASE.pl GENE -L snv.list -X RNA.pileup -Y DNA.pileup -P phased.txt -F 64 -f 33 - N annotation file -C 1 -o outfile

##RELEASE NOTES

V1.0.2, 2016.08.12, fixed bugs in parsing reference skip symbol in pileup file

v1.0.1, 2015.10, origal version
