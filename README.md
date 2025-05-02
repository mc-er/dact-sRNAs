# Differentially targeting analysis of small RNA between _Dactylorhiza_ species

Code and scripts used in the analysis for the publication: <br>
[info to be inserted] // TODO: add citation and doi


## Sources

**BamIndexDecoder** <br>
version: 1.03 <br>
github: https://github.com/gq1/illumina2bam

**CLC Genomics Workbench** <br> 
version: 8.0 <br>
source: [QIAGEN](https://digitalinsights.qiagen.com/products-overview/discovery-insights-portfolio/analysis-and-visualization/qiagen-clc-genomics-workbench/)

**deepTools2** <br>
function: bamCoverage <br>
version: 3.5.1 <br>
github: https://github.com/deeptools/deepTools <br>
publication: https://doi.org/10.1093/nar/gkw257

**EdgeR** <br>
version: 3.34.0 <br>
publication: https://doi.org/10.1093/bioinformatics/btp616; https://doi.org/10.1093/nar/gks042

**R** <br>
version: 4.1.1

**Samtools** <br>
version: 1.10 <br>
github: https://github.com/samtools/samtools <br>
publication: https://doi.org/10.1093/bioinformatics/btr330

**STAR** <br>
version: 2.7.3a <br>
github: https://github.com/alexdobin/STAR <br>
publication: https://doi.org/10.1093/bioinformatics/bts635

**Subread** <br>
function: featureCounts <br>
version: 2.0.0 <br>
doc: https://subread.sourceforge.net/ <br>
publication: https://doi.org/10.1093/bioinformatics/btt656

**topGO** <br>
version: 2.44.0 <br>
doi: https://doi.org/doi:10.18129/B9.bioc.topGO


## Notes


## Extract and preprocess smallRNA reads
// TODO: Add command used to demultiplex the samples

// TODO: Outline how 20-24 nt reads were extracted using CLC GW



## Mapping smallRNA reads to genome
Mapping of the 20-24nt long reads to the _Dactylorhiza incarnata_ reference genome v.1.0 was done using STAR. With the following parameters:
- runThreadN: threads to use
- genomeDir: directory with the reference genome
- readFilesIn: fastq file/-s for the sample
- outFileNamePrefix: prefix to use for the output files
- alignIntronMax: maximum intron length to allow for splicing, set to one since with small RNA we do not expect any splicing
- outFilterMismatchNoverLmax: maximum number of mismatches allowed in the read, set to 0.05
- outFilterMismatchNmax: maximum number of multimappers allowed, set to 100

```
STAR \
    --runThreadN {threads} \
    --genomeDir {/path/to/genomeDir} \
    --readFilesIn {sampleID.fastq} \  
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

## Generate a file with regions of interest
### Windowd read counts
We want to tackle the smallRNA quantification in an annotation free approach, to also get intergenic regions targeted by small RNAs. To do this, we first count the number of reads in a 100 bp window across the genome for each sample. This is done using the `bamCoverage` function from the `deepTools2` package. We use `--normalizeUsing None` to not normalize the read counts, as we want to keep the raw counts for later filtering. 
The output is a `bedgraph` file with the following columns:
- chromosome
- start position
- end position
- number of reads in the window

```
bamCoverage \
    --bam {sampleID_sorted.bam} \
    --outFileName {sampleID_windowed-read-counts-100bp.txt} \
    --outFileFormat bedgraph \
    --normalizeUsing None \
    --binSize 100
```

Outputs a file that looks like this: // TODO: add example
```

```

### Extract smallRNA windows of interest
We then define a window of interest as a 100 bp window that contains at least 10 read in any sample. These windows were extracted using the following steps:
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

### Convert to gff3
As a last step we convert the bed file to a gff3 file. This is done using a custom script `bedcov2gff.py` that was written for this purpose. The script takes the bed file as input and outputs a gff3 file with the following columns:
- scaffold
- source
- feature
- start
- end
- score
- strand
- frame
- attributes (the ID of the feature)

// TODO: add script in script folder

Path to script: // TODO: add path <br>
// TODO: add description of the script <br>
// FIXME: change rtracklayer to bed2gff3 in the script <br>
// FIXME: change script name to bed2gff3

```
python3 bedcov2gff.py \ 
    -i {genome-wide_100bp-windows-of-interest.bed} \
    -o {genome-wide_100bp-windows-of-interest.gff3}
```

## Quantifying smallRNA reads in windows of interest (i.e. features of interest)
// TODO: add description
```
featureCounts \
    -t peak \
    -g ID \
    -T {threads} \
    -M \
    --fraction \
    -a {genome-wide_100bp-windows-of-interest.gff3} \
    -o {featureCounts_MM-fraction.txt} \
    {*_sorted.bam}
```

## Differential targeting (DT) by small RNAs
// TODO: make notebook for DT <br>
**Jupyter notebook:** `.ipynb` 


## Genomic interactions
// TODO: make notebook for genomic interactions <br>
**Jupyter notebook:** `.ipynb` 


## Functional interpretation of differential targeted (DT) genes
// TODO: make notebook for functional interpretation <br>
**Jupyter notebook:** `.ipynb` 

## Intersection with RNA-seq data ([Wolfe et al. 2023](https://doi.org/10.1111/mec.17070))
// TODO: make notebook for RNA seq intersection <br>
**Jupyter notebook:** `.ipynb` 