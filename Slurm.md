# ACE HPC Cluster: SLURM Usage Guide
Welcome to this SLURM Usage Guide! This document is meant to introduce new HPC users on how to effectively and efficiently submit their analysis jobs using SLURM, as well as broaden knowledge for seasoned users on the different options available to better streamline their job submission to the ACE HPC cluster. 

This document will be divided into two major sections to cater for different levels of experience with SLURM:
- SLURM for Absolute Beginners
- SLURM for Advanced Users

## Section 1: SLURM for Absolute Beginners
This section is designed to introduce new users to SLURM, with the aim of helping beginners understand the basic concepts and commands necessary to start using SLURM for managing and running their computational tasks on the ACE HPC cluster. We however assume the basic familiarity of the user with the LINUX systems and basic shell commands. 

### What is SLURM:
SLURM (Simple Linux Utility for Resource Management) is a powerful and flexible job scheduler used on many HPC (High-Performance Computing) systems to efficiently manage and allocate resources (like CPU, memory, GPUs, and more) among users.

### Basic SLURM Concepts
Before we go any further, let us define some key terms that will be used for the rest of this guide tutorial. 
- **Node**: A single computational unit, often equivalent to a single physical server in the cluster. Think of the HPC a cluster a collection of many different computers (nodes) stacked together, where each of these nodes is capable of performing a computational analysis. 
For the case of the ACE, we have types of nodes: the login node (master node) and the compute node. The login node is where everyone is placed when they log into the HPC and the compute nodes are where your analyses (jobs) will be executed when submitted.
- **Job**: A unit of work submitted to SLURM, which can be a single task such as alignment or a complex workflow.
- **Job Script (Slurm Script)**: A script file containing both the SLURM directives and commands to execute the job. It specifies the resources required and the commands to run.

### Creating our first SLURM script
We will create our first slurm script that performs a basic analysis, in this case, performs a quality check on the fastq sample using the fastqc program.

Let us start by downloading a test file to use, `sample.fastq.gz` which can be found [here](https://github.com/acecancer/acecancer.github.io/raw/main/sample.fastq.gz). Create a directory called `slurm_tutorial` and work from that. 

```bash
mkdir slurm_tutorial
cd slurm_tutorial

wget https://github.com/acecancer/acecancer.github.io/raw/main/sample.fastq.gz
```

We can now proceed to create a very basic slurm script for our intended analysis. As mentioned previously, this will contain a **SLURM directives** section which specifies the resources to use for the task and commands to run. Slurm directives start with **#SBATCH** which is how the system is able to distinguish them from your normal commands to perform a task. 

Let us create a script called `01slurmJob.sh`

```bash
#!/bin/bash

#SBATCH --job-name=FastQC
#SBATCH --nodes=1

fastqc sample.fastq.gz
```

**Before we submit our script above, let us first go over what we have**
-  `#!/bin/bash`: This line specifies which shell interpreter we are to use for our script, which, in this case, is bash.
-  `#SBATCH --job-name=FastQC`: This is our first slurm directive line (it starts with #SBATCH), and it defines the name of the job we are to run in this case **FastQC**, using the **--job-name** option. This name definition can be anything you choose to define your job, however, shorter and more precise names are preferred.
-  `#SBATCH --nodes=1`: This is the second slurm directive (also starts with #SBATCH), which specifies how many compute nodes we need to run our particular job, using the **--nodes** option. In this case, we request for only **1 node** which is sufficient for our task. Every ACE node has 32GB of RAM and 16 CPUs.
-  `fastqc sample.fastq.gz`: This specifies the actual job we are to run. In this case, we are running a FastQC quality control check on our fastq file.

At this moment, we should see these two files in our working directory when we list its contents.

```
(base) [fkakembo@kla-ac-hpc-01 slurm_tutorial]$ ls 
01slurmJob.sh  sample.fastq.gz
```

### Submitting a SLURM Job
To submit our analysis job using slurm, we use the **`sbatch`** command followed by our script name. 

```bash
(base) [fkakembo@kla-ac-hpc-01 slurm_tutorial]$ sbatch 01slurmJob.sh
Submitted batch job 6000
```
When a job is submitted successfully, it is assigned an analysis job identifier (JOBID) number. In our case, our JOBID is `6000`. Please note that this will be different for every job submitted, therefore you will have a different JOBID when you submit your job. 

At the same time slurm will create a new file in our working directory named as `slurm-<JOBID>.out`, which is a slurm file to store the log outputs from your analysis. In our case this file will be `slurm-6000.out` since our JOBID was 6000. Again keep in mind that this will be different for you as well. 

```bash
(base) [fkakembo@kla-ac-hpc-01 slurm_tutorial]$ ls
01slurmJob.sh  sample.fastq.gz  slurm-6000.out
```


### Why use SLURM?
- Efficiently manages computational resources.
- Supports a wide range of job types (batch, interactive, parallel).
- Provides tools for job monitoring and management.
- Ensures fair sharing of resources among users.



