## gene analyse

for var in *_R2.fq.gz; do mv "$var" "${var%_R2.fq.gz}.R2.fastq.gz"; done