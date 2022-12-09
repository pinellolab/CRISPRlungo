echo '>chr11_sm' > genome.fa
samtools faidx /srv/local/kendell/genomes/Bowtie2Index/genome.fa chr11:5225364-5225863 | tail -n +2 >> genome.fa
echo '>chr21_sm' >> genome.fa
samtools faidx /srv/local/kendell/genomes/Bowtie2Index/genome.fa chr21:17004905-17005444 | tail -n +2 >> genome.fa
samtools faidx genome.fa
bowtie2-build genome.fa genome
