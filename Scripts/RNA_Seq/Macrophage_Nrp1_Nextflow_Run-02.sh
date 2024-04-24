#!/bin/bash
#!Nextflow bash script 
#!Two batches Macrophages with iWD and cWD, NRP1+/+ vs NRP1 -/- under cWD. 
#!source .bashrc (export nextflow, and singularity)
#! Xiaohui Zhao (xz289@cam.ac.uk)

## GRCm39 .fa and .gtf file need to be compressed.
Basedir=/home/CAD_Projects/CAD_mt709_0001
Ensemdir=/home/Genomes/Mus_musculus/GRCm39
Outdir1=/home/CAD-Projects/CAD_mt709_0001/Macrophage_Batch1_NextFlow_GRCm39
Outdir2=/home/CAD-Projects/CAD_mt709_0001/Macrophage_Batch2_NextFlow_GRCm39
Outdir3=/home/CAD-Projects/CAD_mt709_0001/SLX-22500_NextFlow_GRCm39
export PATH=$PATH:$HOME/packages/NextFlow
NXF_OPTS='-Xms1g -Xmx4g'
export NXF_SINGULARITY_CACHEDIR=$HOME/.singularity

nextflow run nf-core/rnaseq -bg -resume -profile singularity -r 3.2 --skipBiotypeQC \
--input ${Basedir}/Data/First_Batch/Macrophage_FirstBatch_Nextflow_SampleTable.csv \
--aligner star_salmon --fasta ${Ensemdir}/GRCm39.fa.gz --gtf ${Ensemdir}/GRCm39.gtf.gz \
--gtf_extra_attributes 'gene_name' --outdir  ${Outdir1} --multiqc_title Macrophage_FirstBatch_NF_Ensembl_GRCm39 \
--email xz289@cam.ac.uk -with-report ${Outdir1}/report.html &> ${Outdir1}/nextflow_command.log &

nextflow run nf-core/rnaseq -bg -resume -profile singularity -r 3.2 --skipBiotypeQC \
--input ${Basedir}/Data/Second_Batch/Macrophage_SecondBatch_Nextflow_SampleTable.csv \
--aligner star_salmon --fasta ${Ensemdir}/GRCm39.fa.gz --gtf ${Ensemdir}/GRCm39.gtf.gz \
--gtf_extra_attributes 'gene_name' --outdir  ${Outdir2} --multiqc_title NRP1_NF_Ensembl_GRCm39 \
--email xz289@cam.ac.uk -with-report ${Outdir2}/report.html &> ${Outdir2}/nextflow_command.log &

nextflow run nf-core/rnaseq -bg -resume -profile singularity -r 3.2 --skipBiotypeQC \
--input ${Basedir}/Data/SLX-22500/NRP1_Nextflow_SampleTable_cWD.csv \
--aligner star_salmon --fasta ${Ensemdir}/GRCm39.fa.gz --gtf ${Ensemdir}/GRCm39.gtf.gz \
--gtf_extra_attributes 'gene_name' --outdir  ${Outdir3} --multiqc_title NRP1_NF_Ensembl_GRCm39 \
--email xz289@cam.ac.uk -with-report ${Outdir3}/report.html &> ${Outdir3}/nextflow_command.log &