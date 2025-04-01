#!/bin/bash
#SBATCH --job-name=bcftools-filter
#SBATCH -p sonmi
#SBATCH -n 32
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --nodelist=compute-0-1
##SBATCH -mem=200G
##SBATCH -t 72:00:00

ulimit -s unlimited
ulimit -l unlimited

cd $SLURM_SUBMIT_DIR


for i in BSA1-GKO-PRM BSA2-GKO-PRM BSA2-GKO-noPRM BSA2-Neo-PRM;
do
/share/home/hty/software/bcftools-1.10/bcftools view -S ${i}.txt all_dip_extract_addid.vcf -Ov > ${i}_dip_extract_addid.vcf
bgzip ${i}_dip_extract_addid.vcf
tabix -p vcf ${i}_dip_extract_addid.vcf.gz
/share/home/hty/software/bcftools-1.10/bcftools view -R Marker_1065.txt ${i}_dip_extract_addid.vcf.gz > ${i}_dip_Marker_1065.vcf
done
exit
