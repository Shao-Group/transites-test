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
                if fields[2] == 'transcript' and ('havana' in fields[1]) and ('transcript_biotype \"protein_coding\"' in fields[-1]) and fields[0].isnumeric():
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
    labels = pd.DataFrame(np.concatenate((tss_locs_plus_strd, tes_locs_plus_strd, tss_locs_minus_strd, tes_locs_minus_strd)))
    labels.columns = ['chr', 'loc']
    labels['strand'] = ['plus'] * len(tss_locs_plus_strd) * 2 + ['minus'] *  len(tss_locs_minus_strd) *2
    labels['TSS/TES'] = ['tss'] * len(tss_locs_plus_strd) + ['tes'] * len(tes_locs_plus_strd) + \
                        ['tss'] * len(tss_locs_minus_strd) + ['tes'] * len(tes_locs_minus_strd)
    labels = labels.sort_values(by=['chr', 'loc'])
    labels = labels.reset_index(drop=True)
    with open('../processed_data/tss.tes.plu.minus.labels.pkl', 'wb') as labels_file:
        pickle.dump(labels, labels_file)
    return labels


def feat_extraction(sites, bin_size=50):
    def coverage():
        cov = pd.read_csv('../feats/coverage/genome.bga.coverage', sep='\t',
                          header=None, names=['chr', 'start', 'end', 'cov'])
        # keep autosomes only
        cov = cov[cov.chr.apply(lambda x: type(x) == int)]
        # sort accoding to chr:loc, and reset index
        cov = cov.sort_values(by=['chr', 'start'])
        cov = cov.reset_index(drop=True)
        # iterate through rows
        i = 0
        j = 0
        site = sites.loc[j, :]
        while i < cov.shape[0]:
            row = cov.loc[i, :]
            if row['chr'] == site['chr']:
                if row['start'] < site['loc']:
                    if row['end'] >= site['loc']:
                        sites.at[j, 'cov'] = row['cov']
                        j = j + 1
                        site = sites.loc[j, :]
                    else:
                        i = i + 1
                        continue
                else:
                    print("Should not arrive here, coverage, loc error\n", site, row)
                    assert 0
                    sites.at[j, 'cov'] = np.nan
            elif row['chr'] < site['chr']:
                i = i + 1
                continue
            else:
                print('Should not arrive here, coverage, chr num error\n', site, row)
                assert 0
                sites.at[j, 'cov'] = np.nan
        return sites['cov']
    feats = [coverage()]
    return feats


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
            index = labels[labels['chr'] == row['chr']][labels['strand'] == row['strand']][labels['TSS/TES'] == row['TSS/TES']].\
                        index[labels[labels['chr'] == row['chr']][labels['strand'] == row['strand']][labels['TSS/TES'] == row['TSS/TES']]['loc'].searchsorted(row['loc'])]
            if index > 0:
                left_val = labels.loc[index - 1, 'loc']
            else:
                left_val = labels.loc[index - 1, 'loc']
            right_val = labels.loc[index, 'loc']
            if abs(right_val - scallop_labels.loc[i, 'loc']) <= 50 or abs(left_val - scallop_labels.loc[i, 'loc']) <= 50:
                scallop_labels.at[i, 'Real?'] = True
            else:
                scallop_labels.at[i, 'Real?'] = False

        except IndexError:
            scallop_labels.at[i, 'Real?'] = False
    with open('../processed_data/sites_scallop_labeled.pkl', 'wb') as sites_file:
        pickle.dump(scallop_labels, sites_file)

    feat_extraction(scallop_labels)
    # remove duplicates
    rmdup = scallop_labels.drop_duplicates(subset=['chr', 'loc'], keep='first', inplace=False).reset_index(drop=True)
    feat_extraction(rmdup)

    labels_nodup = labels.drop_duplicates(subset=['chr', 'loc'], keep='first', inplace=False).reset_index(drop=True)
    feat_extraction(labels_nodup)

    clf = LinearSVC
    clf_logis = SVC





