# HLA Typing exercise using Polysolver and Optitype

The samples can be found at `/etc/ace-data/CancerGenomicsWG/VariantCalling/samples/Validation`. For demonstration, I will call my sample as `sample`, and its corresponding reads as sample_R1.fastq and sample_R2.fastq for forward and reverse reads respectively. Modify this according to your assigned sample.

## PolySolver

### Steps 1: Alignment with BWA:
We first align our samples to the human reference genome, version 19 (hg19). A copy of this reference and its indices have been added to our shared folder under the location: `/etc/ace-data/CancerGenomicsWG/Tools/References`. However, feel free to use your own copy of the reference in case you have one.

```bash
# Load the bwa and samtools modules
module load bwa-0.7.17-gcc-8.5.0-h3gmzcf samtools-1.13-gcc-8.5.0-hx66cfb

mkdir -p Results/Alignment_hg19

bwa mem -t 16 /etc/ace-data/CancerGenomicsWG/Tools/References/hg19.fa  sample_R1.fastq sample_R2.fastq  -o Results/Alignment_hg19/sample.sam
samtools sort -@ 12 Results/Alignment_hg19/sample.sam -o Results/Alignment_hg19/sample.bam
samtools index Results/Alignment_hg19/sample.bam
```

### Step 2: Create a Polysolver conda environment
We will need to first create a conda environment for this. 

```bash
conda env create --name polysolver  --file /etc/ace-data/CancerGenomicsWG/Tools/polysolver.yml --yes
```

### Step 3: Run the Polysolver script

When the environment is done being created, go ahead and activate it and run the polysolver hla typing script. 

```bash
#Activate the environmnet
conda activate polysolver

mkdir -p Results/polysolver

#Run the polysolver script
shell_call_hla_type Results/Alignment_hg19/sample.bam Black 1 hg19 STDFQ 0 Results/polysolver/sample
```

When the script is done running, check for the file called `winners.hla.txt` in the output directory specified (ie Results/polysolver/sample) in our command above.
Its contents will look something like this:

```
HLA-A   hla_a_01_04n    hla_a_01_04n
HLA-B   hla_b_07_03     hla_b_07_03
HLA-C   hla_c_01_03     hla_c_01_03
```

## Optitype
Read more about optitype at their [Github repo](https://github.com/FRED-2/OptiType)

### Step 1: Creating an Optitype Conda environment

```bash
conda env create --name optitype  --file /etc/ace-data/CancerGenomicsWG/Tools/optitype.yml --yes
```

## Step 2: Running the Optitype tool

```bash
conda activate optitype

python /etc/ace-data/CancerGenomicsWG/Tools/OptiType/OptiTypePipeline.py \
              -i sample_R1.fastq  sample_R2.fastq \
              --dna -v -o ./Results/Optitype -c /etc/ace-data/CancerGenomicsWG/Tools/OptiType/OptiType/config.ini --prefix sample
```

A few things to note:
- Make a copy of the `config.ini` file and include your own version the razers3 program installed in your space, that is `razers3=/path/to/razers3`
- Modifify the input files and sample ids to match your sample assigned.
