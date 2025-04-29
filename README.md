Dactylorhiza small RNAs

Code and scripts used in the analysis for the publication:
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
    --readFilesIn {sampleID.fastq} \  // FIXME: do we have paired end reads?
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
    -o {sampleID_windowedReadCounts100bp.txt} \
    -of bedgraph \
    --normalizeUsing None \
    -bs 100
```

Outputs a file that looks like this: // TODO: add example
```

```

## Extract smallRNA windows of interest
A window of interest is defined as a 100 bp window that contains at least 10 read. The following code extracts windows of interest from the bedgraph file and outputs a new bedgraph file containing only windows containing >= 10 reads.

```