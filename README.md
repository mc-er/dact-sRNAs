# Differentially targeting analysis of small RNA between _Dactylorhiza_ species

Code and scripts used in the analysis for the publication: <br>
[info to be inserted] // TODO: add citation and doi


## Sources

**BamIndexDecoder** <br>
version: 1.03 <br>
git: https://github.com/gq1/illumina2bam 

**CLC Genomics Workbench** <br> 
version: 8.0 <br>
source: [QIAGEN](https://digitalinsights.qiagen.com/products-overview/discovery-insights-portfolio/analysis-and-visualization/qiagen-clc-genomics-workbench/)
 

**STAR** <br>
version: 2.7.3a <br>
git: https://github.com/alexdobin/STAR <br>
publication: https://doi.org/10.1093/bioinformatics/bts635


**Samtools** <br>
version: 1.10 <br>
git: https://github.com/samtools/samtools <br>
publication: https://doi.org/10.1093/bioinformatics/btr330


**deepTools2** <br>
version: 3.5.1 <br>
git: https://github.com/deeptools/deepTools <br>
publication: https://doi.org/10.1093/nar/gkw257


## Notes


## Extract and preprocess smallRNA reads
// TODO: Add command used to demultiplex the samples

// TODO: Outline how 20-24 nt reads were extracted using CLC GW



## Mapping smallRNA reads to genome

```
STAR \
    --runThreadN {threads} \
    --genomeDir {/path/to/genomeDir} \
    --readFilesIn {sampleID.fastq} \  // TODO: check, was data paired end?
    --outFileNamePrefix {sampleID} \
    --alignIntronMax 1 \
    --outFilterMismatchNoverLmax 0.05 \
    --outFilterMulitmapNmax 100

samtools sort \
    --threads {threads} \
    -o {sampleID_sorted.bam} \
    {sampleID.bam}

samtools index \
    --threads {threads} \
    {sampleID_sorted.bam}
```

## Windowd read counts

```
bamCoverage \
    -b {sampleID_sorted.bam} \
    -o {sampleID_windowed-read-counts-100bp.txt} \
    -of bedgraph \
    --normalizeUsing None \
    -bs 100
```

Outputs a file that looks like this: // TODO: add example
```

```

## Extract smallRNA windows of interest
A window of interest is defined as a 100 bp window that contains at least 10 read in any sample. These windows were extracted using the following steps:
- merging all sample `bedgraph` files
- filtering for windows with at least 10 reads in any sample (column 4)
- extracting the first three columns (chromosome, start, end)
- sorting the output by chromosome and start position
- removing duplicates
- as a last step, neighboring windows were merged into a larger continuous window
- save as a `bed` file

```
cat {*_windowed-read-counts-100bp.txt} | \
awk -F'\t' '$4>=10' | \
cut -f1,2,3 | \
sort -k1,1 -k2,2 -V | \
uniq | \
bedtools merge -i - > {genome-wide_100bp-windows-of-interest.bed}
```

## Convert to gff3
// TODO: add script in script folder

Path to script: // TODO: add path
// TODO: add description of the script
// TODO: double check command and if its actually gff3
```
python3 bedcov2gff.py \ 
    -i {genome-wide_100bp-windows-of-interest.bed} \
    -o {genome-wide_100bp-windows-of-interest.gff3}
```



