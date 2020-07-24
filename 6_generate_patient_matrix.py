"""

Generate a DataFrame-style matrix for each patient, visit X feature.

"""

import pandas as pd
import numpy as np
from scipy import stats
import json
import os
import pickle as pkl
import matplotlib.pyplot as plt


VERSION = 'COVID-0519'


def mkdir(path):
    folder = os.path.exists(path)

    if not folder:
        os.makedirs(path)

    else:
        pass


def load_intubation_status(version=VERSION, save=True):
    df = pd.DataFrame(columns=['clientguid', 'is_intubated'])

    event_df = pd.read_csv(version + '/' + version+'_event.csv', sep=',', header=0, dtype={'clientguid':str}).fillna(value='None')

    for idx, row in event_df.iterrows():
        clientguid, intubation_start = row['clientguid'], row['intubation_start']
        if intubation_start != 'None':
            intubation_status = 1
        else:
            intubation_status = 0
        df = df.append(pd.DataFrame({'clientguid':clientguid, 'is_intubated':intubation_status}, index=["0"]), ignore_index=True)

    if save == True:
        df.to_csv(version+'/'+version+'_intubation_status.csv', index=False, header=True)
    return df


def load_intubation_status_dict(version=VERSION):
    intubation_status_dic = {}
    event_df = pd.read_csv(version + '/' + version+'_event.csv', sep=',', header=0, dtype={'clientguid':str}).fillna(value='None')

    for idx, row in event_df.iterrows():
        clientguid, intubation_start = row['clientguid'], row['intubation_start']
        if intubation_start != 'None':
            intubation_status = 1
        else:
            intubation_status = 0
        intubation_status_dic[clientguid] = intubation_status
    return intubation_status_dic


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



def filtering(patient_matrices, min_visit_length, counting_point=0):
    """
    To filter out patients whose number of visits (after counting point) < min_visit_length.
    """
    new_patient_matrices = {}
    print(len(patient_matrices))
    for p in patient_matrices:

        v_count = 0
        patient_df_whole = patient_matrices[p]
        patient_df = patient_df_whole[['Coagulation', 'Liver', 'Central_nervous_system', 'Renal', 'Cardiovascular', 'Respiration']]

        for v in range(counting_point, len(patient_df)):
            visit_data = patient_df.iloc[v].values[1:]
            if np.count_nonzero(visit_data == visit_data) != 0:
                v_count += 1

        if v_count >= min_visit_length:
            new_patient_matrices[p] = patient_df_whole
    print("Filtering out patients has less than %s visit...\n" % min_visit_length)
    print("Total patients: %s\n" % len(patient_matrices))
    print("Filtered out: %s\n" % (len(patient_matrices) - len(new_patient_matrices)))
    print("Patients left: %s\n" % len(new_patient_matrices))
    return new_patient_matrices


def append_additional_seq(bucked_data, add_feature, seq_data_file, interval, version):
    if "AdmitDtm" in seq_data_file:
        event = "AdmitDtm"
    else:
        event = "intubation_start"

    jsonf = open(version + '/additional_seq/' + version + "_%s_%s_interval_%s.json" % (add_feature, event, interval), 'r')
    new_data = json.load(jsonf)
    jsonf.close()

    bucked_data[add_feature] = {"Observation_name":add_feature, "Unit": "", "data":new_data}
    return bucked_data




def generate_patient_matrix(seq_data_file, version=VERSION, slide_window=30, interval=24, max_visit_before_event=30, max_visit_after_event=42, min_visit_length=10, save=True):
    save_dir = version + "/patient_longitudinal_tables_" + seq_data_file[:-5] + "_maxBef_%s_maxAft_%s_minLen_%s/" % (max_visit_before_event, max_visit_after_event, min_visit_length)

    mkdir(save_dir)

    #--------- load patient list, containing patients who have admission time
    patient_list = []
    event_df = pd.read_csv(version + '/' + version+'_event.csv', sep=',', header=0, dtype={'clientguid':str}).fillna(value='None')
    for idx, row in event_df.iterrows():
        clientguid = row['clientguid']
        patient_list.append(clientguid)

    #--------- load bucked sequence data
    jsonf = open(version + '/' + seq_data_file, 'r')
    bucked_data = json.load(jsonf)
    jsonf.close()

    #--------- load additional sequence data
    bucked_data = append_additional_seq(bucked_data=bucked_data, add_feature="Cardiovascular", seq_data_file=seq_data_file, interval=interval, version=version)
    bucked_data = append_additional_seq(bucked_data=bucked_data, add_feature="Respiration", seq_data_file=seq_data_file, interval=interval, version=version)

    # get feature list
    feature_list = list(bucked_data.keys())


    #--------- generate patient matrix
    # Initialization
    patient_matrices = {}  # key: clientguid, value: visit X feature matrix
    for p in patient_list:
        patient_df = pd.DataFrame(columns=['Visit']+feature_list)
        # print([i for i in range(-max_visit_before_event, int(max_visit_after_event/interval))])
        patient_df["Visit"] = [i for i in range(-int(max_visit_before_event/interval), int(max_visit_after_event/interval))]
        patient_matrices[p] = patient_df

    # Updating
    for feature in bucked_data:
        print("Processing feature: %s...\n" % feature)
        feature_data = bucked_data[feature]['data']
        for p in feature_data:

            patient_feature_data = feature_data[p]


            patient_feature_data = patient_feature_data[int((slide_window - max_visit_before_event)/interval):]


            for v in range(min(len(patient_feature_data), int((max_visit_after_event+max_visit_before_event)/interval))):
                value = float(patient_feature_data[v])
                if p not in patient_matrices:
                    # print(p)
                    continue
                patient_matrices[p].loc[v, feature] = value


    # clculate sofa total score
    print("Generating sofa total score...\n")
    for p in patient_matrices:
        patient_df = patient_matrices[p]
        patient_df['Sofa_total_score'] = patient_df['Cardiovascular'] + patient_df['Respiration'] + patient_df['Coagulation'] + \
                                         patient_df['Liver'] + patient_df['Central_nervous_system'] + patient_df['Renal']
        patient_matrices[p] = patient_df
    feature_list.append('Sofa_total_score')

    # Filtering
    patient_matrices = filtering(patient_matrices, min_visit_length, counting_point=int(max_visit_before_event/interval))

    # Saving
    if save == True:

            patient_matrices[p].to_csv(save_dir + '/' + p + ".csv", sep=",", index=False, header=True)
        with open(save_dir + '/' + '_all.pkl', 'wb') as wf:
            pkl.dump(patient_matrices, wf)
        with open(save_dir + '/' + '_features.txt', 'w') as wf:
            for feature in feature_list:
                wf.write("%s\n" % feature)






def main():
    aligned_event = 'intubation_start'  # intubation_start, AdmitDtm
    slide_window = 24 * 15  # hour
    interval = 24
    method = 'worse'  # 'np.mean'
    max_visit_before_event = 24 * 15 # <= slide_window
    max_visit_after_event = 24 * 60  # hour
    min_visit_length = 1

    plot_length_before_event = 24  # <= max_visit_before_event
    plot_length_after_event = 24 * 30  # <= max_visit_after_event

    seq_data_file = "bucked_alignedEvent_%s_slide_%s_interval_%s_method_%s.json" % (aligned_event, slide_window, interval, method)

    generate_patient_matrix(seq_data_file=seq_data_file, version=VERSION, slide_window=slide_window, interval=interval,
                            max_visit_before_event=max_visit_before_event, max_visit_after_event=max_visit_after_event,
                            min_visit_length=min_visit_length)

    

main()