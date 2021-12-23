#qsub -pe make 10 -l h_vmem=25G freebayes_all.sh

source /home/soft/miniconda3/bin/activate /home/soft/miniconda3/envs/freebayes/
freebayes -f /piec1/dades/genomes/ccaretta/map_feb2021_70/GENOME.populations.loci.34.fa -L /piec1/c.pegueroles/pau/mapping/BAM/listOfBAM.txt > /piec1/c.pegueroles/pau/mapping/SNP_calling/pau_freebayes_all.vcf
deactivate
