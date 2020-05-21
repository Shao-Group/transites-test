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
from sklearn.preprocessing import normalize
from sklearn.model_selection import cross_validate
from sklearn.metrics import recall_score
from sklearn.decomposition import PCA
from sklearn.metrics import roc_auc_score
from sklearn.metrics import roc_curve
from sklearn.metrics import auc


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
def feat_extraction(sites, chr_col=0, loc_col=1, bin_size=100): # checkbin size
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
        cov_feat_bin = np.zeros((sites_vals.shape[0], bin_size*2+1))
        cov_feat_bin[:] = np.NaN
        while i < cov.shape[0] and j < sites_vals.shape[0]:
            site_val = sites_vals[j]
            q = 0
            p = i
            # check if need to increment i
            if cov[i][end_col] < site_val[loc_col] - bin_size and cov[i][chr_col] <= site_val[chr_col]:
                i = i + 1
                continue
            while q < bin_size * 2 +1 and p < cov.shape[0]:
                loc = site_val[loc_col] - bin_size + q
                if loc < 0:
                    cov_feat_bin[j, q] = np.nan
                    q = q + 1
                else:
                    row = cov[p]
                    if row[chr_col] == site_val[chr_col]:
                        if row[start_col] < loc:
                            if row[end_col] >= loc:
                                cov_feat_bin[j, q] = row[cov_col]
                                q = q + 1
                            else:
                                p = p + 1
                        else:
                            print("Should not arrive here, coverage, loc error\n", site_val, row)
                            cov_feat_bin[j, q] = np.nan
                            p = p + 1
                    elif row[chr_col] < site_val[chr_col]:
                        i = i + 1
                        p = i
                    else:
                        print('Should not arrive here, coverage, chr num error\n', site_val, row)
                        cov_feat_bin[j, q] = np.nan
                        # need to increment j, break
                        # j = j + 1
                        break
            j = j + 1
        for i in range(len(sites_vals)):
            if sites_vals[i][2] == 'minus':
                cov_feat_bin[i] = np.flip(cov_feat_bin[i], axis=0)
        return cov_feat_bin

    def chip_seq():
        files = ['GSM1869138_BJ_PolII_MeDiChi-seq_peaks.bed', 'GSM1869137_BJ_H3K27me3_SICER.bed',
                 'GSM1869134_BJ_H3K4me3_MeDiChi-seq_peaks.bed', 'GSM1869136_BJ_H3K27ac_MeDiChi-seq_peaks.bed',
                 'GSM1869135_BJ_H3K9ac_MeDiChi-seq_peaks.bed']
        chip_list = []
        for file in files:
            signal = pd.read_csv('../feats/chip-seq/' + file, sep='\t', header=0)
            signal = signal.rename(
                columns={'end ': 'end', 'chrom': 'chromosome', 'ChIP_island_read_count': 'intensity'})
            signal.chromosome = signal.chromosome.apply(
                lambda x: int(x.split('chr')[-1]) if x.split('chr')[-1].isnumeric() else x)
            signal = signal[signal.chromosome.apply(lambda x: type(x) == int)]
            # sort according to chr:loc, and reset index
            signal = signal.sort_values(by=['chromosome', 'start'])
            signal = signal.reset_index(drop=True)
            col_name = file.split('_')[2]
            cov = signal.values
            start_col = 1
            end_col = 2
            if file == '../feats/chip-seq/GSM1869137_BJ_H3K27me3_SICER.bed':
                cov_col = 3
            else:
                cov_col = 4
            # iterate through rows
            i = 0
            j = 0
            cov_feat_bin = np.zeros((sites_vals.shape[0], bin_size * 2 + 1))
            cov_feat_bin[:] = np.NaN
            while i < cov.shape[0] and j < sites_vals.shape[0]:
                site_val = sites_vals[j]
                q = 0  # q-th neighbor
                p = i  # index of cov
                # check if need to increment i
                if cov[i][end_col] < site_val[loc_col] - bin_size and cov[i][chr_col] <= site_val[chr_col]:
                    i = i + 1
                    continue
                while q < bin_size * 2 + 1 and p < cov.shape[0]:
                    loc = site_val[loc_col] - bin_size + q
                    if loc < 0:
                        cov_feat_bin[j, q] = np.nan
                        q = q + 1
                    else:
                        row = cov[p]
                        if row[chr_col] == site_val[chr_col]:
                            if row[start_col] < loc:
                                if row[end_col] >= loc:
                                    cov_feat_bin[j, q] = row[cov_col]
                                    q = q + 1
                                else:
                                    p = p + 1
                            else:
                                # print("Should not arrive here, coverage, loc error\n", site_val, row)
                                cov_feat_bin[j, q] = np.nan
                                q = q + 1
                                p = p + 1
                        elif row[chr_col] < site_val[chr_col]:
                            i = i + 1
                            p = i
                        else:
                            # print('Should not arrive here, coverage, chr num error\n', site_val, row)
                            cov_feat_bin[j, q] = np.nan
                            # need to increment j, break
                            # j = j + 1
                            break
                j = j + 1
            # print(i,j)
            for i in range(len(sites_vals)):
                if sites_vals[i][2] == 'minus':
                    cov_feat_bin[i] = np.flip(cov_feat_bin[i], axis=0)
            chip_list.append(cov_feat_bin)
        return chip_list

    sites_vals = sites.values
    cov_bin = coverage()
    chip_list = chip_seq()
    feats = [cov_bin, chip_list]
    return feats


# get candidate sites from bam
# this is used if extends model to broader candidates
# def sites_read():
#     # try load results
#     try:
#         with open('../processed_data/sites.pkl', 'rb') as sites_file:
#             sites = pickle.load(sites_file)
#             return sites
#     except FileNotFoundError:
#         pass
#     cov_fs = ['../feats/coverage/genome.plus.5.cov', '../feats/coverage/genome.plus.3.cov',
#               '../feats/coverage/genome.minus.5.cov', '../feats/coverage/genome.minus.3.cov']
#     suspected_labels = ['tss', 'tes', 'tss', 'tes']
#     strds = ['plus', 'plus', 'minus', 'minus']
#     sites = []
#     for i in range(4):
#         with open(cov_fs[i], 'r') as cov_f:
#             for line in cov_f.readlines():
#                 fields = line.split()
#                 if fields[0].isnumeric():
#                     sites.append([int(fields[0]), int(fields[2]), strds[i], suspected_labels[i]])
#     sites = pd.DataFrame(sites)
#     sites.columns = ['chr', 'loc', 'strand', 'TSS/TES']
#     # sort accoding to chr:loc, and reset index
#     sites = sites.sort_values(by=['chr', 'loc'])
#     sites = sites.reset_index(drop=True)
#     with open('../processed_data/sites.pkl', 'wb') as sites_file:
#         pickle.dump(sites, sites_file)
#     return sites


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


def add_labels_for_scallop(labels, scallop_labels):
    chr_col = 0
    loc_col = 1
    strand_col = 2
    ts = 3
    label_val = labels.values
    scallop_labels_val = scallop_labels.values
    to_label = []
    j = 0
    for i in range(scallop_labels_val.shape[0]):
        j = max(0, j - 100)
        # print(i,j)
        while j < label_val.shape[0]:
            row = scallop_labels_val[i]
            if label_val[j][chr_col] == row[chr_col]:
                diff = label_val[j][loc_col] - row[loc_col]
                if diff < -50:
                    j = j + 1
                    continue
                elif diff < 50:
                    if label_val[j][strand_col] == row[strand_col] and label_val[j][ts] == row[ts]:
                        to_label.append(True)
                        break
                    else:
                        j = j + 1
                        continue
                else:
                    # diff > 50
                    to_label.append(False)
                    break
            elif label_val[j][chr_col] < row[chr_col]:
                j = j + 1
                continue
            else:
                # label_val[j][chr_col] > row[chr_col]
                to_label.append(False)
                break
        if j >=label_val.shape[0]:
            to_label.append(False)
    return to_label


if __name__ == "__main__":
    # get candidate sites and assign labels
    labels = label_extraction()
    scallop_labels = sites_read_scallop()
    scallop_labels_list = add_labels_for_scallop(labels, scallop_labels)
    scallop_labels['Real?'] = scallop_labels_list
    # remove duplicated lines
    scallop_labels = scallop_labels.drop_duplicates(subset=['chr', 'loc'], keep='first', inplace=False).reset_index(drop=True)
    # extract features
    feats = feat_extraction(scallop_labels)
    feats_cat = np.concatenate([feats[0], np.concatenate(feats[1], axis=1)], axis=1)
    feats_cat = np.nan_to_num (feats_cat)
    # this is the column numbers of corresponding identifiers
    chr_col = 0
    loc_col =1
    strand_col = 2
    ts = 3

    # use RNA reads of 200 neighbors and chip-seq of the locations of interest
    mask = np.concatenate([np.arange(201), [301, 502, 703, 904, 1105]])
    y = list(scallop_labels['Real?'])
    # train model and get scores
    for i in [100]:
        clf = LinearSVC(class_weight='balanced', C=2, dual=True)
        sl2 = feats_cat[:, mask]
        # normalize data and use chi-square transformation
        sl2 = normalize(sl2)
        kernel = Nystroem('chi2', n_components=201)
        sl2 = kernel.fit_transform(sl2)
        # select top 100 features with PCA
        pca = PCA(n_components=i)
        pca.fit(sl2)
        sl2 = pca.transform(sl2)
        # perform 5 cross-validation
        scoring = ['precision_macro', 'recall_macro', 'f1', 'accuracy', 'roc_auc']
        scores = cross_validate(clf, sl2, np.array(y), scoring=scoring, cv=5, return_train_score=True)
        for x in scores:
            print(x, scores[x], sum(scores[x]) / 5, np.std(scores[x]))


    # Visualize with PCA
    sl2 = feats_cat[:, mask]
    sl2 = normalize(sl2)
    kernel = Nystroem('chi2', n_components=201)
    sl2 = kernel.fit_transform(sl2)
    pca = PCA(n_components=2)
    pca.fit(sl2)
    sl2 = pca.transform(sl2)
    print(pca.explained_variance_ratio_)   # 0.36897021, 0.08642462
    sns.scatterplot(sl2[:,0], sl2[:,1], y)
    plt.show()

    # make roc plot
    # reference https://scikit-learn.org/stable/auto_examples/model_selection/plot_roc.html#sphx-glr-auto-examples-model-selection-plot-roc-py
    x_train = sl2[:-len(sl2) // 5]
    y_train = y[:-len(sl2) // 5]
    x_test = sl2[-len(sl2) // 5:]
    y_test = y[-len(sl2) // 5:]
    clf.fit(x_train, y_train);
    y_score = clf.decision_function(x_test)
    falposit, prec, _ = roc_curve(y_test, y_score)
    auc_roc = auc(falposit, prec)  # 0.8617685034357059
    plt.plot(falposit, prec, label='ROC curve')
    plt.show()


