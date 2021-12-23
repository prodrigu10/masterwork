#!/bin/bash

# Cargar el pipeline: TRIMMING de secuencias

# qsub -l h_vmem=50G trimSeq.sh 
bash /piec1/dades/pipelines/2b_deNovo/trimming_denovo/deNovoPipe_AlfI_uncompressed.sh /piec1/c.pegueroles/pau/sequences/
