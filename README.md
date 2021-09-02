## gene analyse

for var in *_R2.fq.gz; do mv "$var" "${var%_R2.fq.gz}.R2.fastq.gz"; done


snp indel注释
java -Xmx10G -jar snpEff.jar eff -c snpEff.config AT_10  input/merge.Spring_P50.indel.filter.vcf > positive.snp.eff.vcf -csvStats positive.csv -stats positive.htmlls