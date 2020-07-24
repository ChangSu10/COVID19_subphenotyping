"""

Clustering analysis


1. Select features of interest, 6 SOFA sub-scores;

2. Imputing the sub-scores;

3. Generating the total score;

4. Visualization

5. DTW

6. Clustering

"""

import pandas as pd
import numpy as np
from scipy import stats
from tslearn.metrics import cdist_dtw
import scipy.spatial as sp, scipy.cluster.hierarchy as hc
from sklearn import cluster

import json
import os
import pickle as pkl
import matplotlib.pyplot as plt
import seaborn as sns
import datetime
import time
import copy




VERSION = "COVID-0512"




def mkdir(path):
    folder = os.path.exists(path)

    if not folder:
        os.makedirs(path)

    else:
        pass

def load_feature_unit(version=VERSION):
    unit_dict = {}
    try:
        df = pd.read_csv(version + "/feature_unit.csv", sep=",", header=0).dropna()
        for idx, row in df.iterrows():
            feature, unit = row["FeatureName"], row["Unit"]
            unit_dict[feature] = unit
        return unit_dict
    except:
        print("Did not find feature_unit.csv")
        return {}


def load_death_info(version=VERSION, window=1):
    intubation_info = {}
    event_df = pd.read_csv(version + '/' + version + '_event.csv', sep=',', header=0, dtype={'clientguid':str}).fillna('None')
    for idx, row in event_df.iterrows():
        p, intubation_start = row['clientguid'], row['intubation_start']
        intubation_info[p] = intubation_start

    death_info = {}
    death_df = pd.read_csv(version + '/CEDAR/patient_death.csv', sep=',', header=0, dtype={'CLIENTGUID':str})
    for idx, row in death_df.iterrows():
        p, death_time = row['CLIENTGUID'], row['DEATH_DATE']
        if p not in intubation_info:
            continue
        else:
            intubation_time = intubation_info[p]
            if intubation_time == 'None':
                continue
            intubation_time = time.mktime(time.strptime(intubation_time, "%Y-%m-%d %H:%M:%S"))
            death_time = time.mktime(time.strptime(death_time, "%Y-%m-%d %H:%M:%S"))

            intubation_death_window = (death_time - intubation_time) / 60 / 60

            if intubation_death_window <= window * 24:
                death_info[p] = "Death within %sh of intubation" % (window * 24)
            else:
                death_info[p] = "Death after %sh of intubation"  % (window * 24)

    # format as DataFrame
    death_info_df = pd.DataFrame(columns=['clientguid', 'group'])
    idx = 0
    for p in death_info:
        death_info_df.loc[idx, ['clientguid', 'group']] = [p, death_info[p]]
        idx += 1
    return death_info, death_info_df


def length_distribution(patient_matrices_folder, version=VERSION):
    # load patient data
    infile = open(version + "/" + patient_matrices_folder + "/_all.pkl", 'rb')
    patient_matrices = pkl.load(infile)
    infile.close()

    data_length = []
    for p in patient_matrices:
        patient_mtx = patient_matrices[p].dropna(subset=['Visit'])
        patient_mtx = patient_matrices[p].dropna(how='all', subset=patient_mtx.columns.tolist()[1:])
        max_visit = max(patient_mtx['Visit'].values.tolist()) + 1
        data_length.append(max_visit)

    plt.hist(data_length)
    plt.xlabel("Visit length")
    plt.savefig(version + "/" + patient_matrices_folder + "/_visit_length.pdf")
    plt.close()







def is_all_nan(seq):
    """
    :param seq: 1-D numpy array
    """
    if np.count_nonzero(seq == seq) == 0:
        return True
    else:
        return False

def LOCF_FOCB(seq):
    """
    :param seq: 1-D sequence (pandas series)
    :return: imputed sequence
    """
    N = len(seq)
    seq = seq.fillna(-1).values  # fill NaN as -1

    # Do LOCF
    mask_idx = np.where(seq == -1)[0]
    mask_array = np.ones(N, dtype='int')
    mask_array[mask_idx] = 0
    missing_num = len(mask_idx)

    for i in range(missing_num):
        idx = mask_idx[i]
        if idx != 0:   # is the first entry
            if int(seq[idx - 1]) != -1:
                seq[idx] = seq[idx - 1]

    # Do FOCB
    mask_idx = np.where(seq == -1)[0]
    mask_array = np.ones(N, dtype='int')
    mask_array[mask_idx] = 0
    missing_num = len(mask_idx)
    for i in range(missing_num):
        idx = mask_idx[missing_num - i - 1]
        if idx != N - 1:
            if int(seq[idx + 1]) != -1:
                seq[idx] = seq[idx + 1]
    return seq


def filtering(patient_matrices, features, min_visit_length=4, counting_point=0):
    """
    To filter out patients whose number of visits (after counting point) < min_visit_length.
    """
    new_patient_matrices = {}
    for p in patient_matrices:

        v_count = 0
        patient_df = patient_matrices[p][features]

        for v in range(counting_point, len(patient_df)):
            visit_data = patient_df.iloc[v].values[1:]
            if np.count_nonzero(visit_data == visit_data) != 0:
                v_count += 1

        if v_count >= min_visit_length:
            new_patient_matrices[p] = patient_df
    print("Filtering out patients has less than %s visit...\n" % min_visit_length)
    print("Total patients: %s\n" % len(patient_matrices))
    print("Filtered out: %s\n" % (len(patient_matrices) - len(new_patient_matrices)))
    print("Patients left: %s\n" % len(new_patient_matrices))
    return new_patient_matrices

def imputing(patient_matrices_folder, version=VERSION, save=True):
    # load patient data
    infile = open(VERSION + "/" + patient_matrices_folder + "/_all.pkl", 'rb')
    patient_matrices = pkl.load(infile)
    infile.close()

    for p in patient_matrices:
        patient_matrices[p] = patient_matrices[p].dropna(subset=['Visit'])

    # load unit info
    unit_dict = load_feature_unit(version=version)



    features = ['Coagulation', 'Liver', 'Central_nervous_system', 'Renal', 'Cardiovascular', 'Respiration']


    # filtering   --  filter out patients whose SOFA records <= 4 visits
    patient_matrices = filtering(patient_matrices, features=features, min_visit_length=4, counting_point=0)


    # calculate median of the features:
    feature_median = {}
    for p in patient_matrices:

        patient_df = patient_matrices[p]
        for f in features:
            values = patient_df[f].values.tolist()
            if f not in feature_median:
                feature_median[f] = values
            else:
                feature_median[f] += values
    for f in feature_median:
        feature_median[f] = np.nanmedian(feature_median[f])


    # imputation
    new_patient_matrices = {}
    for p in patient_matrices:

        print("Imputing patient: %s..." % p)
        patient_df = patient_matrices[p][features]

        for f in features:
            feature_column = patient_df[f]
            if is_all_nan(feature_column.values):  # Missing all values of the feature
                if f == "Respiration":
                    patient_df[f] = 0.0
                else:
                    patient_df[f] = feature_median[f]

            else:   # Missing partial values of the feature

                imputed_seq = LOCF_FOCB(feature_column)

                patient_df[f] = imputed_seq
                if p != "9000400429900200":
                    print(f)
                    print(feature_column.values)
                    print(imputed_seq)

        new_patient_matrices[p] = patient_df

    print("Total %s patients." % len(new_patient_matrices))


    # clculate sofa total score
    print("Generating sofa total score...\n")
    for p in new_patient_matrices:
        patient_df = new_patient_matrices[p]
        patient_df['Sofa_total_score'] = patient_df['Cardiovascular'] + patient_df['Respiration'] + patient_df['Coagulation'] + \
                                         patient_df['Liver'] + patient_df['Central_nervous_system'] + patient_df['Renal']
        new_patient_matrices[p] = patient_df

    features.append('Sofa_total_score')

    # Saving
    save_dir = version + '/' + patient_matrices_folder + "/_SOFA_imputed"
    mkdir(save_dir)
    if save == True:
        for p in new_patient_matrices:
            new_patient_matrices[p].to_csv(save_dir + '/' + p + ".csv", sep=",", index=False, header=True)
        with open(save_dir + '/' + '_all.pkl', 'wb') as wf:
            pkl.dump(new_patient_matrices, wf)
        with open(save_dir + '/' + '_features.txt', 'w') as wf:
            for feature in features:
                wf.write("%s\n" % feature)



def baseline_strata(patient_matrices_folder, baseline_visit, version=VERSION):
    # load feature list
    feature_list = []
    with open(VERSION + "/" + patient_matrices_folder + "/_SOFA_imputed/_features.txt") as f:
        all_lines = f.readlines()
        for line in all_lines:
            feature_list.append(line.strip())

    # load patient data
    infile = open(VERSION + "/" + patient_matrices_folder + "/_SOFA_imputed/_all.pkl", 'rb')
    patient_matrices = pkl.load(infile)
    infile.close()
    print("Total %s patients.\n" % len(patient_matrices))

    data = pd.DataFrame(columns=['clientguid'] + feature_list)
    idx = 0
    patients = []
    for p in patient_matrices:
        patient_data = patient_matrices[p]
        baseline = patient_data.loc[baseline_visit].values.tolist()

        data.loc[idx] = [p] + baseline

        patients.append(p)
        idx += 1

    save_dir = "results/" + version + "/baseline_cluster/"
    mkdir(save_dir)

    label_df = pd.DataFrame(columns=['clientguid', 'group'])
    idx = 0
    for p in patient_matrices:
        baseline_sofa = patient_matrices[p].loc[baseline_visit, 'Sofa_total_score']
        if baseline_sofa <= 10:
            label = 0
        elif baseline_sofa <= 12:
            label = 1
        else:
            label = 2

        label_df.loc[idx] = [p, label]
        idx += 1

    label_df.to_csv(save_dir + 'partitioning.csv', sep=',', index=False, header=True)



def subtyping(patient_matrices_folder, baseline_group, cohort, cluster_num, version=VERSION, data_type='both', decentring=False):
    """
    :param data_type: univariate, multivariate, both
    """
    # load feature list
    feature_list = []
    with open(VERSION + "/" + patient_matrices_folder + "/_SOFA_imputed/_features.txt") as f:
        all_lines = f.readlines()
        for line in all_lines:
            feature_list.append(line.strip())

    # load patient data
    infile = open(VERSION + "/" + patient_matrices_folder + "/_SOFA_imputed/_all.pkl", 'rb')
    patient_matrices = pkl.load(infile)
    infile.close()
    print("Total %s patients.\n" % len(patient_matrices))

    # death info
    death_info, _ = load_death_info(version=version)

    # save dir
    save_path = "results/" + version + "/" + patient_matrices_folder + "/_clustering_combine/"
    mkdir(save_path)


    # ---------------- univariate version
    if data_type == 'both' or data_type == 'univariate':
        for feature in feature_list:
            if feature not in ['Sofa_total_score']:
                continue
            print("Clustering with univariate: %s...\n" % feature)

            # formatting
            data = []
            patients = []
            for p in patient_matrices:
                if p not in cohort:
                    continue

                

                patient_data = patient_matrices[p][feature].values
                if decentring == True:  # decentering
                    patient_data_mean = np.nanmean(patient_data, axis=0)
                    patient_data_max = np.nanmax(patient_data, axis=0)
                    patient_data_min = np.nanmin(patient_data, axis=0)


                    if patient_data_max == patient_data_min:
                        patient_data = np.zeros(len(patient_data))
                    else:
                        patient_data = (patient_data - patient_data_min) / (patient_data_max - patient_data_min)


                    # patient_data = patient_data - patient_data_mean

                data.append(patient_data)
                patients.append(p)
            print("%s patients remaining for analysis...\n" % len(data))

            # calculating similarity by DTW
            dist_mtx = cdist_dtw(data)

            dist_mtx_df = pd.DataFrame(dist_mtx)
            dist_mtx_df.to_csv(save_path+'distance_mtx_baselinegroup_%s.csv' % baseline_group, index=True, header=True)


            # draw clustergram
            for method in ['average']:   # 'average', 'complete', 'ward', 'single'
                linkage = hc.linkage(sp.distance.squareform(dist_mtx, checks=False), method=method)
                ns_plot = sns.clustermap(dist_mtx, row_linkage=linkage, col_linkage=linkage, standard_scale=1)
                plt.savefig(save_path+'univariate_%s_decentring_%s_method_%s_baselinegroup_%s.png' % (feature, decentring, method, baseline_group), dpi=300)
                plt.close()

                if feature in ['Sofa_total_score']:  # ['Respiration' , 'Sofa_total_score']
                    label_df = pd.DataFrame(columns=['clientguid', 'group'])
                    label_df['clientguid'] = patients
                    model = cluster.AgglomerativeClustering(n_clusters=cluster_num, affinity='precomputed', linkage=method)
                    model.fit(dist_mtx)
                    labels = model.labels_

                    label_df['group'] = labels
                    label_df.to_csv(save_path+'univariate_%s_decentring_%s_method_%s_baselinegroup_%s.csv' % (feature, decentring, method, baseline_group), sep=',', index=False, header=True)

                print("%s ----- %s ------- finished...\n" % (feature, method))








def main():
    aligned_event = 'intubation_start'  # intubation_start, AdmitDtm
    slide_window = 0
    interval = 24
    method = 'worse' 
    max_visit_before_event = 0
    baseline_visit = 0 
    length_day = 7
    max_visit_after_event = 24 * length_day
    min_visit_length = 4


    plot_length_before_event = 0  # <= max_visit_before_event
    plot_length_after_event = int((24 / interval) * length_day)  # <= max_visit_after_event

    patient_matrices_folder = "patient_longitudinal_tables_bucked_alignedEvent_%s_slide_%s_interval_%s_method_%s_maxBef_%s_maxAft_%s_minLen_%s" % (aligned_event, slide_window, interval, method, max_visit_before_event, max_visit_after_event, min_visit_length)

    length_distribution(patient_matrices_folder=patient_matrices_folder, version=VERSION)

    imputing(patient_matrices_folder=patient_matrices_folder, version=VERSION)

    baseline_strata(patient_matrices_folder=patient_matrices_folder, baseline_visit=baseline_visit, version=VERSION)

    baseline_clus_df = pd.read_csv("results/" + VERSION + "/" + "baseline_cluster" + "/" + "partitioning" + ".csv", sep=',', header=0).astype(str)
    baseline_groups = baseline_clus_df.drop_duplicates(subset=['group'])['group'].values.tolist()
    cluster_num_dic = {"0": 2, "1": 2, "2": 2}

    for bg in baseline_groups:
        sub_cohort = baseline_clus_df[baseline_clus_df['group'] == bg]['clientguid'].values.tolist()
        k = cluster_num_dic[bg]
        subtyping(patient_matrices_folder=patient_matrices_folder, baseline_group=bg, cohort=sub_cohort, cluster_num=k, version=VERSION, data_type='both', decentring=True)


main()