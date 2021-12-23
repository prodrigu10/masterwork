#qsub -pe make 1 -l h_vmem=15G /piec1/c.pegueroles/pau/mapping/sam2sortedbam.sh 

MAPPING="/piec1/c.pegueroles/pau/mapping"

for file in ${MAPPING}/SAM/*sam #single end
do
	#echo $file
	withpath="${file}"
	filename=${withpath##*/}
	echo $filename
	base="${filename%.sam}"
	echo "${base}"
	/home/soft/samtools-0.1.19/samtools view -S -b ${MAPPING}/SAM/$filename > ${MAPPING}/BAM/${base}.bam
	/home/soft/samtools-0.1.19/samtools sort ${MAPPING}/BAM/${base}.bam ${base}.sorted
	mv /home/usuaris/c.pegueroles/${base}.sorted.bam ${MAPPING}/BAM
done

#rm *sam
