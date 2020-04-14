'''
    A script for TSS/TES identification with SVM
'''
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import pickle


def label_extraction():
    genome_gtf_file = '../data/GRCh38.97.gtf'
    # try load results
    try:
        with open('../processed_data/tss.tes.plu.minus.labels.pkl', 'rb') as labels_file:
            tss_locs_plus_strd, tes_locs_plus_strd, tss_locs_minus_strd, tes_locs_minus_strd = pickle.load(labels_file)
            return tss_locs_plus_strd, tes_locs_plus_strd, tss_locs_minus_strd, tes_locs_minus_strd
    except FileNotFoundError:
        pass
    tss_locs_plus_strd = []
    tes_locs_plus_strd = []
    tss_locs_minus_strd = []
    tes_locs_minus_strd = []
    extracted_table = []
    with open(genome_gtf_file, 'r') as gtf:
        for line in gtf.readlines():
            # for each non-comment line
            if not line.startswith('#'):
                fields = line.strip().split('\t')
                # if it's a coding transcript from curated havana database, and on autosomes
                if fields[2] == 'transcript' and ('havana' in fields[1]) and ('transcript_biotype \"protein_coding\"' in fields[-1]) and fields[0].isnumeric():
                    # TODO: save extracted lines in a table
                    # depends on the strand, TSS/TES are reversed
                    if fields[6] == '+':
                        tss_locs_plus_strd.append((fields[0], int(fields[3])))
                        tes_locs_plus_strd.append((fields[0], int(fields[4])))
                    elif fields[6] == '-':
                        tss_locs_minus_strd.append((fields[0], int(fields[4])))
                        tes_locs_minus_strd.append((fields[0], int(fields[3])))
                    else:
                        print("Warning, strand information wrong", line)
    with open('../processed_data/tss.tes.plu.minus.labels.pkl', 'wb') as labels_file:
        pickle.dump([tss_locs_plus_strd, tes_locs_plus_strd, tss_locs_minus_strd, tes_locs_minus_strd], labels_file)
    return tss_locs_plus_strd, tes_locs_plus_strd, tss_locs_minus_strd, tes_locs_minus_strd

def cov_extraction(sites):
    cov_file = '../feat/'


def sites_read():
    pass
# class Site:
#     def __init__(self, chr_num, loc, cov_file):
#         self.cov_file = cov_file
#         self._chr_num = chr_num
#         self._loc = loc
#     def feat_cov(self):


if __name__ == "__main__":
    tss_plus, tes_plus, tss_minus, tes_minus = label_extraction()
    sites = sites_read()




