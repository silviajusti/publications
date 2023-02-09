## This pipelines assumes there is a text file with the sample names within the same foler as the scripts. Sample names should be one per line. 
### Sample name is just the identifier for the specimen without any file extensions or designations e.g. An_gambiae_1234.fastq.gz - The samplename is An_gambiae_1234). The text file with thebsmaple names will be referred here as sample_list.txt. This can contain only one or several samples. 

## if this pipeline was helpful, please cite **https://doi.org/10.1371/journal.pone.0247068** and the corresponding citations for the listed software. 
### This is just a guide on how I do these analyses. I seriously encourage you to read through the manuals and original papers of the software listed below, in dorder to be able to adapt this to your own samples /  research questions. :) 
### Note that each server is going to require distinct headers for memory requirements, etc. Also, the server I use, uses qsub to submit the jobs, make sure you comply with your servers requirements for job submission. 
### I am also listing a few alternative analytical tools for each step, but the pipeline is based on the ones I use more frequently
### I will only list the commands for the analyses and not the memory requirements, etc.
###  I always assume the commands will be given from the cwd
### This is an updated/easier way to run "https://github.com/silviajusti/publications/blob/master/SCO_pipeline.md"

### Sample name (referred to throughout as "samplename" is just the identifier for the specimen without any file extensions or designations e.g. An_gambiae_1234.fastq.gz - The samplename is An_gambiae_1234)

## Illumina reads QC
Input must be fastq format!\
First we run [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).
Create a run_fastqc.sh file containing the following commands 

```bash
fastqc $1_R1.fastq.gz
fastqc $1_R2.fastq.gz
```

# To submit
```bash
for filename in $(cat sample_list.txt )
do  
qsub -o $filename-fastqc-log run_fastqc.sh $filename
done
```

## Now trim and filter sequences according to the observed fastqc results
Will use [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic).
Create a run_trimmomatic.sh file containing the following commands 

Alternatively, use [trim_galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)

## The MAXINFO flag will depend on your samples and overall quality and size distribution. Because I wrote this initially for collection samples, the small fragments are quite abundant and important, so I kept the minimum length at 30. 

```bash
runtrimmomatic PE -phred33 $1_R1.fastq.gz $1_R2.fastq.gz $1_R2_paired.fq.gz \
$1_R1_unpaired.fq.gz $1_R2_paired.fq.gz $1_R2_unpaired.fq.gz MAXINFO:30:0.8
```

## To submit
```bash
for filename in $(cat sample_list.txtt )
do  
qsub -o $filename-trimmomati-log run_trimmomatic.sh $filename
done
```

* trimmomatic will output paired and unpaired reads. Only the paired reads will be used for further analyses.


## Genome assembly
After trimming and filtering the reads assembly is performed with [gatb-minia pipeline](https://github.com/GATB/gatb-minia-pipeline)

Alternatively use [AbySS](https://github.com/bcgsc/abyss) or [MaSuRCA](https://github.com/alekseyzimin/masurca). MaSuRCA is especially good for a combo of long and short reads, but also give good results with either. 
 
Create a run_gatbminia.sh file containing the following commands 
*remember to set the number of cores, as the default takes up all cores available!*

```bash
gatb -1 $1_R1_paired.fq.gz -2 $1_R2_paired.fq.gz --nb-cores 12
```

## To submit
```bash
for filename in $(cat sample_list.txtt )
do  
qsub -o $filename-gatb-log run_gatb.sh $filename
done
```
Note that the resulting assembly.fasta is an alias, and it is a good idea to copy it to a new file (cp assembly.fasta samplename_assembly.fasta) Do not use mv, as this will only change the name of the alias. Also, the pipeline creates a lot of intermediate files that I don't use, so after results inspections, I delete those to save storage space.


## Genome annotation/gene prediction
Now, run [Augustus](http://bioinf.uni-greifswald.de/augustus/) to annotate the assembly.
Create a run_augustus.sh file containing the following commands 
Make sure to set the species to the closest available for your sample, or use [BUSCO](https://busco.ezlab.org/) to generate training gene sets

```bash
augustus --species=aedes $1_assembly.fasta > $1.augustus.out
```

## To submit
```bash
for filename in $(cat sample_list.txtt )
do  
qsub -o $filename-augustus-log run_augustus.sh $filename
done
```

## Now convert the the sequences from the annotation gff files into fasta using [getAnnoFasta.pl]( https://raw.githubusercontent.com/nextgenusfs/augustus/master/scripts/getAnnoFasta.pl)

```bash
perl getAnnoFasta.pl --seqfile=$1_assembly.fasta $1.augustus.out
```

## To submit
```bash
for filename in $(cat sample_list.txtt )
do  
qsub -o $filename-getAnnoFasta-log run_getAnnoFasta.sh $filename
done
```

## To run orthograph (assuming you already have the whole folder set up - see bellow). We will use the \*.codingseq output from getAnnoFasta.pl

```bash
cp samplename.codinseq /Path_to_orthograph_folder/input/samplename.codinseq.fasta
```



## Single-copy orthologs (SCO) search
Assuming the ortograph database was already created, see [orthoset_construction](https://github.com/jsoghigian/orthoset_construction)

### To run orthograph 
* The directory should contain everything needed to run Orthograph on an input transcriptome, gene set, or sequence capture set.
* It expects the following:
1) You have Orthograph and SQLite installed.
2) Input files are .fasta and in the input directory.
3) You have properly configured the .conf file for your particular paths to programs. 
It should not be necessary to change relative paths to files, because that is handled in the script below.
4) You will use the launch loop script below to identify orthologs across a set of fasta input - you can modify it as needed, 
or take individual bits and pieces as you need.
5) You change the paths to the programs as appropriate for your system in the launch loop script below.
6) You have a ref_tax.txt list of taxa used to create the database
## You can find an example folder containing everything needed to run orthograph, as described, [here](https://doi.org/10.25573/data.22060787.v1)
* The script below does the following:
1) It creates a new config file for an input fasta file.
2) It runs Orthograph Analyzer, which identifies putative orthologs in the input file. Output is written to output/taxname
3) It runs Orthograph Reporter, which chooses the best hit based on HMMER3 and BLAST. Output is written to output/taxname
4) It summarizes the Orthograph results, which produces .nt.summarized and .aa.summarized files for each ortholog... and removes the 
reference species from the fasta file outputs.  These files will be located within done/taxname/nt_summarized and done/taxname/aa_summarized .  
5) It then compresses the output directory, because it contains a very large number of files, and a large .sqlite file that takes up a lot of space  
and then removes the uncompressed directory.
The database folder will be the relative path to the database sqlite-database on orthograph.config

```bash
for taxa in $(ls input/*.fasta);
do
taxa=$(basename ${taxa} .fasta)
sed "s/basetax/${taxa}/g" orthograph.conf > configs/orthograph_${taxa}.conf
orthograph-analyzer --configfile configs/orthograph_${taxa}.conf
orthograph-reporter --configfile configs/orthograph_${taxa}.conf
mkdir done/${taxa}
perl summarize_orthograph_results.pl -i output -o done/${taxa} -c -t -s -u -d ref_tax.txt
cd output
tar -czvf ${taxa}.tar.gz ${taxa}/
rm -r ${taxa}/
cd ..
done
```
