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
- **Node**: A single computational unit, often equivalent to a single physical server in the cluster. Think of the HPC a cluster a collection of many different computers (nodes) stacked together, where each of these nodes are capable of performing a computational analysis. 
For the case of the ACE, we have types of nodes that is the login node (master node) and the compute node. The login node is where everyone is placed when they log into the HPC and the compute nodes are where your analyses (jobs) will be executed when submitted.
- **Job**: A unit of work submitted to SLURM, which can be a single task such as alignment or a complex workflow.
- **Job Script (Slurm Script)**: A script file containing both the SLURM directives and commands to execute the job. It specifies the resources required and the commands to run.


