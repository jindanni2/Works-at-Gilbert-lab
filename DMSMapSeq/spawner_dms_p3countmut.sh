#/bin/bash
D=/home/dj448/palmer_scratch/DMS_v2_020524
for i in ${D}/*__barcode_Bar2_1.fastq.gz; do
  name=${i}
  name2=${i:0:-10}2.fastq.gz
  echo ${name}
  echo ${name2}
  sbatch wrapper_dms_p3countmut.slurm ${name} ${name2}
done
