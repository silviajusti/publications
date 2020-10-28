## The pipeline was written to accomodate several samples at once (i.e. for batch analyses). 

# Illumina reads QC
Needs a file called list.txt, with one column containing the sample names without spaced.\
Input must be fastq format!\
Assumes commands will always be given within the same folder.\
First we run [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).


```bash
for filename in $(cat list.txt)
do 
fastqc ${filename}_R1_001.fastq.gz
fastqc ${filename}_R2_001.fastq.gz
done
```


Now trim and filter sequences according to the observed fastqc results
Will use [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
Create a run_trimmomatic_sample.sh file containing the following commands 

```bash
runtrimmomatic PE -phred33 sample_R1.fastq sample_R1.fastq sample_R2_paired.fq.gz \
sample_R1_unpaired.fq.gz sample_R2_paired.fq.gz sample_R2_unpaired.fq.gz MAXINFO:50:0.8
```


Now, replace "sample" with the sample names using the list.txt file.

```bash
for filename in $(cat list.txt)
        do
sed "s/sample/$filename/g" run_trimmomatic_sample.sh >run_trimmomatic_${filename}.sh
mkdir ${filename}
mv run_trimmomatic_${filename}.sh *${filename}/
done
```


And then submmit all jobs

```bash
for filename in $(cat list.txt)
        do
        cd ${filename}/
  qsub run_trimmomatic_${filename}.sh   
  cd ..  
  done
```
* trimmomatic will output paired and unpaired reads. Only the paired reads will be used.

# Genome assembly
After trimming and filtering the reads assembly is performed with [gatb-minia pipeline](https://github.com/GATB/gatb-minia-pipeline) \
Create a run_gatbminia_sample.sh file containing the following commands 

```bash
gatb -1 sample_R1_paired.fq.gz -2 sample_R2_paired.fq.gz --nb-cores 12
```


Now, replace "sample" with the sample names using the list.txt file.

```bash
for filename in $(cat list.txt)
        do
sed "s/sample/$filename/g" run_gatbminia_sample.sh >run_gatbminia_${filename}.sh
mv run_gatbminia_${filename}.sh ${filename}/
done
```

And then submmit all jobs

```bash
for filename in $(cat list.txt)
        do
        cd ${filename}/
  qsub run_gatbminia_${filename}.sh  
  cd ../  
  done
```


When gatbminia finishes running, the file name for the assembly is consistently assembly.fasta, now, we change it to the name of the sample. 

```bash
for filename in $(cat list.txt)
        do
        cp ${filename}/assembly.fasta ${filename}/${filename}_assembly.fasta
    done

```

# Genome annotation/gene prediction
Now, run [Augustus](http://bioinf.uni-greifswald.de/augustus/) to annotate the assembly. \
Create a run_augustus_sample.sh file containing the following commands 

```bash
augustus --species=aedes sample_assembly.fasta > sample.augustus.out
```

Now, replace "sample" wwith the sample names using the list.txt file.

```bash
for filename in $(cat list.txt)
        do
sed "s/sample/$filename/g" run_augustus_sample.sh > run_augustus_${filename}.sh
mv run_augustus_${filename}.sh ${filename}/
done
```

And then submmit all jobs

```bash
for filename in $(cat list.txt)
        do
        cd ${filename}/
  qsub run_augustus_${filename}.sh 
  cd ../  
  done
```


Now convert the the sequences from the annotation gff files into fasta using [getAnnoFasta.pl]( https://raw.githubusercontent.com/nextgenusfs/augustus/master/scripts/getAnnoFasta.pl)

```bash
for filename in $(cat list.txt)
do
perl getAnnoFasta.pl --seqfile=${filename}/${filename}_assembly.fasta ${filename}/${filename}.augustus.out
 done
```

# Single-copy orthologs (SCO) search
Assuming the ortograph database was already created, see [orthoset_construction](github.com/jsoghigian/orthoset_construction)
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
mkdir configs done input output database

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

# Results summary
Run blastn and blastp for all the orthologs identified by orthograph
Create a run_blast_sample.sh file containing the following commands 

```bash
for dir in `ls done/sample/nt_summarized`
do
for filename in `ls $dir`
do   
blastn -query ${dir}/${filename} -db nt -out ${dir}/blastn_${filename}  -outfmt "6 std staxid " -max_target_seqs 1
done
done

for dir in `ls done/sample/aa_summarized`
do
for filename in `ls $dir`
do   
blastp -query ${dir}/${filename} -db nr -out ${dir}/blastp_${filename}  -outfmt "6 std staxid " -max_target_seqs 1
done
done
```

Now, replace "sample" with the sample names using the list.txt file.

```bash
for filename in $(cat list.txt)
        do
sed "s/sample/$filename/g" run_blast_sample.sh > run_blast_${filename}.sh
mv run_blast_${filename}.sh ${filename}/
done
```
Concatenate the blast results by specimen

```bash
mkdir blast_results

for filename in `ls path/to/ortograph/done/`
do
cat path/to/ortograph/done/${filename}/nt_summarized/blastn* >> ${filename}_blastn.txt
cat path/to/ortograph/done/${filename}/aa_summarized/blastp* >> ${filename}_blastp.txt
done
```

To batch check the blast results, download the taxonomy database from NCBI and extract only the target taxon (in this case, Arthropoda)

```bash
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
#extract only the file with the full lineage names and make a file containing only the arthropod lineages (or whatever group you are working on)
tar --extract --file=new_taxdump.tar.gz fullnamelineage.dmp
grep Arthropoda fullnamelineage.dmp | awk '{print $1}'  > tax_list.txt
```


For each of the files, find the orthologs that are, actually, from the target taxon. 
For this, list all taxid on the blast results that correspond to the target taxa

```bash
for filename in `ls *blastn*` ; do awk '{print $10}' ${filename} | grep -w  -F -f - tax_list.txt; done >> temp_ID.txt
for filename in `ls *blastp*` ; do awk '{print $13}' ${filename} | grep -w  -F -f - tax_list.txt; done >> temp_ID.txt
sort temp_ID.txt | uniq| sed -z 's/\n/\\|/g'|sed 's/..$//' >> sample_taxid.txt
```


Create a bash file to output the counts and the SCO of interest

```bash
nano give_target_taxa.sh
```


And include the following script

```bash
#!/usr/bin/env bash
filename=${1?Error: input blast results file}
 
grep $(cat sample_taxid.txt) $filename >> $filename.out
wc -l $filename.out>> $filename.count
awk -F"|" '{print $1}' $filename.out>> $filename.target_orthologs
```

Make it executable

```bash
chmod +x give_target_taxa.sh>> 
```

And then run for all samples


```bash
for filename in $(cat sample_list.txt )
do 
mkdir ${filename}
mv *${filename}* ${filename}/
./give_target_taxa.sh ${filename}/${filename}.txt
done
```

