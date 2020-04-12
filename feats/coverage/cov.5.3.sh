bedtools genomecov -bga -strand + -5 -ibam ../data/SRR307903.star.sort.bam  > genome.plus.5.cov
bedtools genomecov -bga -strand - -5 -ibam ../data/SRR307903.star.sort.bam  > genome.minus.5.cov
bedtools genomecov -bga -strand + -3 -ibam ../data/SRR307903.star.sort.bam  > genome.plus.3.cov
bedtools genomecov -bga -strand - -3 -ibam ../data/SRR307903.star.sort.bam  > genome.minus.3.cov
