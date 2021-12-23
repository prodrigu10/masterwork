DIR="/piec1/c.pegueroles/pau/sequences/Trimmed"
GENOME_DIR="/piec1/dades/genomes/ccaretta/map_feb2021_70"
MAPPING="/piec1/c.pegueroles/pau/mapping"

N_CORES=5

rm -f batch_commands_mapping.sh

for file in ${DIR}/*fasta #single end
do
	#echo $file
	withpath="${file}"
	filename=${withpath##*/}
	#echo $filename
	base="${filename%.fasta}"
	#echo "${base}"
	echo "qsub -pe make ${N_CORES} -l h_vmem=15G ${MAPPING}/${base}_mapping.sh" >> batch_commands_mapping.sh

echo "/home/soft/hisat2-2.2.1/hisat2 -p ${N_CORES} -x ${GENOME_DIR}/GENOME.populations.loci.34 -f ${DIR}/$filename -S ${MAPPING}/${base}_mapping.sam --rg-id ${base} --rg SM:${base} &> ${MAPPING}/${base}_mapping.sam.info" > ${base}_mapping.sh

done

#to launch all jobs, connect to cluster and 
#bash batch_commands_mapping.sh
