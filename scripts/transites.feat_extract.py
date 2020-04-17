'''
    A script for TSS/TES identification with SVM
'''
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import pickle
from sklearn.svm import LinearSVC
from sklearn.kernel_approximation import Nystroem


# extract true labels from GRCh38.gtf
def label_extraction():
    genome_gtf_file = '../data/GRCh38.97.gtf'
    # try load results
    try:
        with open('../processed_data/tss.tes.plu.minus.labels.pkl', 'rb') as labels_file:
            labels = pickle.load(labels_file)
            return labels
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
                if fields[2] == 'transcript' and ('havana' in fields[1]) and (
                        'transcript_biotype \"protein_coding\"' in fields[-1]) and fields[0].isnumeric():
                    # TODO: save extracted lines in a table
                    # depends on the strand, TSS/TES are reversed
                    if fields[6] == '+':
                        tss_locs_plus_strd.append((int(fields[0]), int(fields[3])))
                        tes_locs_plus_strd.append((int(fields[0]), int(fields[4])))
                    elif fields[6] == '-':
                        tss_locs_minus_strd.append((int(fields[0]), int(fields[4])))
                        tes_locs_minus_strd.append((int(fields[0]), int(fields[3])))
                    else:
                        print("Warning, strand information wrong", line)
    labels = pd.DataFrame(
        np.concatenate((tss_locs_plus_strd, tes_locs_plus_strd, tss_locs_minus_strd, tes_locs_minus_strd)))
    labels.columns = ['chr', 'loc']
    labels['strand'] = ['plus'] * len(tss_locs_plus_strd) * 2 + ['minus'] * len(tss_locs_minus_strd) * 2
    labels['TSS/TES'] = ['tss'] * len(tss_locs_plus_strd) + ['tes'] * len(tes_locs_plus_strd) + \
                        ['tss'] * len(tss_locs_minus_strd) + ['tes'] * len(tes_locs_minus_strd)
    labels = labels.sort_values(by=['chr', 'loc'])
    labels = labels.reset_index(drop=True)
    with open('../processed_data/tss.tes.plu.minus.labels.pkl', 'wb') as labels_file:
        pickle.dump(labels, labels_file)
    return labels


# extract features - coverage, chip-seq
def feat_extraction(sites, chr_col=0, loc_col=1, bin_size=20):
    def coverage():
        # read coverage in a 2d-array
        try:
            with open('../processed_data/genome.bga.coverage.pkl', 'rb') as f:
                cov = pickle.load(f)
        except FileNotFoundError:
            cov = pd.read_csv('../feats/coverage/genome.bga.coverage', sep='\t',
                              header=None, names=['chr', 'start', 'end', 'cov'])
            # keep autosomes only
            cov = cov[cov.chr.apply(lambda x: type(x) == int)]
            # sort accoding to chr:loc, and reset index
            cov = cov.sort_values(by=['chr', 'start'])
            cov = cov.reset_index(drop=True)
            cov = cov.values
            with open('../processed_data/genome.bga.coverage.pkl', 'wb') as f:
                pickle.dump(cov, f)
        start_col = 1
        end_col = 2
        cov_col = 3
        # iterate through rows
        i = 0
        j = 0
        q = 1
        cov_feat_bin = np.zeros((sites_vals.shape[0], bin_size*2))
        while i < cov.shape[0] and j < sites_vals.shape[0]:
            site_val = sites_vals[j]
            row = cov[i]
            if row[chr_col] == site_val[chr_col]:
                if row[start_col] < site_val[loc_col]:
                    if row[end_col] >= site_val[loc_col]:
                        cov_feat_bin[j, q] = row[cov_col]
                        j = j + 1
                    else:
                        i = i + 1
                        continue
                else:
                    print("Should not arrive here, coverage, loc error\n", site_val, row)
                    cov_feat_bin[j, q] = np.nan
                    j = j + 1
            elif row[chr_col] < site_val[chr_col]:
                i = i + 1
                continue
            else:
                print('Should not arrive here, coverage, chr num error\n', site_val, row)
                cov_feat_bin[j, q] = np.nan
                j = j + 1
        return cov_feat_bin

    def chip_seq():
        files = ['GSM1869138_BJ_PolII_MeDiChi-seq_peaks.bed', 'GSM1869137_BJ_H3K27me3_SICER.bed',
                 'GSM1869134_BJ_H3K4me3_MeDiChi-seq_peaks.bed', 'GSM1869136_BJ_H3K27ac_MeDiChi-seq_peaks.bed',
                 'GSM1869135_BJ_H3K9ac_MeDiChi-seq_peaks.bed']
        for file in files:
            signal = pd.read_csv('../feats/chip-seq/'+file, sep='\t', header=0)
            signal = signal.rename(columns={'end ': 'end', 'chrom': 'chromosome', 'ChIP_island_read_count':'intensity'})
            signal.chromosome = signal.chromosome.apply(lambda x: int(x.split('chr')[-1]) if x.split('chr')[-1].isnumeric() else x )
            signal = signal[signal.chromosome.apply(lambda x: type(x) == int)]
            # sort according to chr:loc, and reset index
            signal = signal.sort_values(by=['chromosome', 'position'])
            signal = signal.reset_index(drop=True)
            col_name = file.split('_')[2]
            i = 0
            j = 0
            while i < signal.shape[0] and j < sites.shape[0]:
                site = sites.loc[j, :]
                print(i, j)
                row = signal.loc[i, :]
                if row['chromosome'] == site['chr']:
                    if row['start'] < site['loc']:
                        if row['end'] >= site['loc']:
                            sites.at[j, col_name] = row['intensity']
                            j = j + 1
                            # site = sites.loc[j, :]
                        else:
                            i = i + 1
                            continue
                    else:
                        sites.at[j, col_name] = 0
                        j = j + 1
                        # site = sites.loc[j, :]
                elif row['chromosome'] < site['chr']:
                    i = i + 1
                    continue
                else:
                    sites.at[j, col_name] = 0
                    j = j + 1
                    site = sites.loc[j, :]

    sites_vals = sites.values
    coverage()


# get candidate sites from bam
def sites_read():
    # try load results
    try:
        with open('../processed_data/sites.pkl', 'rb') as sites_file:
            sites = pickle.load(sites_file)
            return sites
    except FileNotFoundError:
        pass
    cov_fs = ['../feats/coverage/genome.plus.5.cov', '../feats/coverage/genome.plus.3.cov',
              '../feats/coverage/genome.minus.5.cov', '../feats/coverage/genome.minus.3.cov']
    suspected_labels = ['tss', 'tes', 'tss', 'tes']
    strds = ['plus', 'plus', 'minus', 'minus']
    sites = []
    for i in range(4):
        with open(cov_fs[i], 'r') as cov_f:
            for line in cov_f.readlines():
                fields = line.split()
                if fields[0].isnumeric():
                    sites.append([int(fields[0]), int(fields[2]), strds[i], suspected_labels[i]])
    sites = pd.DataFrame(sites)
    sites.columns = ['chr', 'loc', 'strand', 'TSS/TES']
    # sort accoding to chr:loc, and reset index
    sites = sites.sort_values(by=['chr', 'loc'])
    sites = sites.reset_index(drop=True)
    with open('../processed_data/sites.pkl', 'wb') as sites_file:
        pickle.dump(sites, sites_file)
    return sites


# get candidate sites from scallop
def sites_read_scallop():
    # try load results
    try:
        with open('../processed_data/sites_scallop.pkl', 'rb') as sites_file:
            scallop_labels = pickle.load(sites_file)
            return scallop_labels
    except FileNotFoundError:
        pass
    tss_locs_plus_strd = []
    tes_locs_plus_strd = []
    tss_locs_minus_strd = []
    tes_locs_minus_strd = []
    with open('../data/SRR307903.scallop.gtf', 'r') as scallop_f:
        for line in scallop_f.readlines():
            fields = line.split()
            if fields[2] == 'transcript' and fields[0].isnumeric():
                # depends on the strand, TSS/TES are reversed
                if fields[6] == '+':
                    tss_locs_plus_strd.append((int(fields[0]), int(fields[3])))
                    tes_locs_plus_strd.append((int(fields[0]), int(fields[4])))
                elif fields[6] == '-':
                    tss_locs_minus_strd.append((int(fields[0]), int(fields[4])))
                    tes_locs_minus_strd.append((int(fields[0]), int(fields[3])))
                else:
                    print("Warning, strand information wrong", line)
    scallop_labels = pd.DataFrame(
        np.concatenate((tss_locs_plus_strd, tes_locs_plus_strd, tss_locs_minus_strd, tes_locs_minus_strd)))
    scallop_labels.columns = ['chr', 'loc']
    scallop_labels['strand'] = ['plus'] * len(tss_locs_plus_strd) * 2 + ['minus'] * len(tss_locs_minus_strd) * 2
    scallop_labels['TSS/TES'] = ['tss'] * len(tss_locs_plus_strd) + ['tes'] * len(tes_locs_plus_strd) + \
                                ['tss'] * len(tss_locs_minus_strd) + ['tes'] * len(tes_locs_minus_strd)
    scallop_labels = scallop_labels.sort_values(by=['chr', 'loc'])
    scallop_labels = scallop_labels.reset_index(drop=True)
    with open('../processed_data/sites_scallop.pkl', 'wb') as sites_file:
        pickle.dump(scallop_labels, sites_file)
    return scallop_labels


if __name__ == "__main__":
    labels = label_extraction()
    scallop_labels = sites_read_scallop()
    for i in range(scallop_labels.shape[0]):
        row = scallop_labels.loc[i, :]
        try:
            index = \
            labels[labels['chr'] == row['chr']][labels['strand'] == row['strand']][labels['TSS/TES'] == row['TSS/TES']]. \
                index[labels[labels['chr'] == row['chr']][labels['strand'] == row['strand']][
                labels['TSS/TES'] == row['TSS/TES']]['loc'].searchsorted(row['loc'])]
            if index > 0:
                left_val = labels.loc[index - 1, 'loc']
            else:
                left_val = labels.loc[index - 1, 'loc']
            right_val = labels.loc[index, 'loc']
            if abs(right_val - scallop_labels.loc[i, 'loc']) <= 50 or abs(
                    left_val - scallop_labels.loc[i, 'loc']) <= 50:
                scallop_labels.at[i, 'Real?'] = True
            else:
                scallop_labels.at[i, 'Real?'] = False

        except IndexError:
            scallop_labels.at[i, 'Real?'] = False
    with open('../processed_data/sites_scallop_labeled.pkl', 'wb') as sites_file:
        pickle.dump(scallop_labels, sites_file)

    # remove duplicates
    rmdup = scallop_labels.drop_duplicates(subset=['chr', 'loc'], keep='first', inplace=False).reset_index(drop=True)
    feat_extraction(rmdup)

    labels_nodup = labels.drop_duplicates(subset=['chr', 'loc'], keep='first', inplace=False).reset_index(drop=True)
    feat_extraction(labels_nodup)

    clf = LinearSVC(class_weight = 'balanced')
    sl = sites.loc[:, ['chr', 'loc', 'cov', 'PolII', 'H3K27me3', 'H3K4me3', 'H3K27ac', 'H3K9ac']]
    clf.fit(sl, sites['Real?'])

    kernel = Nystroem('sigmoid', n_components=300)
    kernel = kernel.fit_transform(sl)
    clf.fit(kernel,  sites['Real?'])

