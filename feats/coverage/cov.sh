bedtools genomecov -bga -strand + -split -ibam ../data/SRR307903.star.sort.bam  > genome.plus.split.cov
bedtools genomecov -bga -strand - -split -ibam ../data/SRR307903.star.sort.bam  > genome.minus.split.cov
bedtools genomecov -bga -strand + -ibam ../data/SRR307903.star.sort.bam  > genome.plus.cov
bedtools genomecov -bga -strand - -ibam ../data/SRR307903.star.sort.bam  > genome.minus.cov
