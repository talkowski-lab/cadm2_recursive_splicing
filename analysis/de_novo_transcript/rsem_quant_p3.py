#!/usr/bin/env python3

import os
import subprocess

def kallisto(a, b, rsf):
    custom_output_root = "/data/talkowski/Samples/cadm2/data/de_novo/cadm2_cap_rna/rsem_expression"
    job_dir = os.path.join(custom_output_root, "jobs")
    os.makedirs(job_dir, exist_ok=True)

    with open(a, "r") as infile:
        for line in infile:
            line = line.strip()
            t = line.split("/")
            folder0 = "/".join(t[:-1])
            sampleName = t[-1]
            refinedName = sampleName

            refindex = os.path.join(b, "cadm2_iN_all_transcripts.fa")
            fq1 = os.path.join(folder0, sampleName, sampleName + ".R1.fq")
            fq2 = os.path.join(folder0, sampleName, sampleName + ".R2.fq")

            rsemOutFolder = os.path.join(custom_output_root, sampleName, rsf)
            os.makedirs(rsemOutFolder, exist_ok=True)

            # Output and error logs
            stde = os.path.join(job_dir, f"{rsf}_{refinedName}_RSEM.err")
            stdo = os.path.join(job_dir, f"{rsf}_{refinedName}_RSEM.out")
            slurm_script_path = os.path.join(job_dir, f"{rsf}_{refinedName}_RSEM.slurm")

            with open(slurm_script_path, "w") as slurm_script:
                slurm_script.write(f"""#!/bin/bash
#SBATCH --job-name={refinedName}
#SBATCH --output={stdo}
#SBATCH --error={stde}
#SBATCH --partition=short
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00

module load Perl/5.28.0-GCCcore-7.3.0
module load rsem/1.3.3
module load bowtie2/2.4.2

rsem-calculate-expression \\
    --forward-prob=0 \\
    --bowtie2 \\
    --bowtie2-mismatch-rate 0.05 \\
    --estimate-rspd \\
    --keep-intermediate-files \\
    --paired-end \\
    {fq1} \\
    {fq2} \\
    {refindex} \\
    {rsemOutFolder}/RSEM
""")

            # Submit SLURM job
            subprocess.call(["sbatch", slurm_script_path])

    return 1

# Example usage
kallisto(
    "/data/talkowski/Samples/cadm2/data/de_novo/cadm2_cap_rna/rsem_index/aligned_bam_files.txt",
    "/data/talkowski/Samples/cadm2/data/de_novo/cadm2_cap_rna/rsem_index",
    "cap_rna_trinity_rsem00"
)
