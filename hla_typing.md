# HLA Typing exercise using Polysolver and Optitype

The samples can be found at `/etc/ace-data/CancerGenomicsWG/VariantCalling/samples/Validation`. For demonstration, I will call my sample as `sample`, and its corresponding reads as sample_R1.fastq and sample_R2.fastq for forward and reverse reads respectively. Modify this according to your assigned sample.

## PolySolver

### Steps 1: Alignment with BWA:

```
module load bwa-0.7.17-gcc-8.5.0-h3gmzcf samtools-1.13-gcc-8.5.0-hx66cfb

mkdir -p Results/Alignment_hg19

bwa mem -t 16 ~/References/hg19.fa  sample_R1.fastq sample_R2.fastq  -o Results/Alignment_hg19/sample.sam
samtools sort -@ 12 Results/Alignment_hg19/sample.sam -o Results/Alignment_hg19/sample.bam
samtools index Results/Alignment_hg19/sample.bam

mkdir -p Results/polysolver
```
