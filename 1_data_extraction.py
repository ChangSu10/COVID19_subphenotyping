"""
Extract feature from raw data

    Entire data: list of length M
        Every feature is a list
        You need to have an attribute encoding the name of the feature

    For each feature, construct a patient list
        The length of the list could be smaller than N because the feature could be missing for specific patients
        Every patient also need to have an attribute encoding the MRN/ClientGuid
            Investigate the issue of one MRN associated with multiple ClientGuid

    For each patient, construct a feature value list
        Each entry is a two-tuple
            Feature value
            Timestamp
        Sort the entries according to timestamps

"""

import pandas as pd
import numpy as np
import json
import time
import math

VERSION = 'COVID-0519'
FILE_NAME = VERSION + '/CEDAR/observations.csv'
FILE_NAME_INO = VERSION + '/CEDAR/observations_ino.csv'


# to determine if an object is digit (int or float)
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        pass

    try:
        import unicodedata
        unicodedata.numeric(s)
        return True
    except (TypeError, ValueError):
        pass

    return False


def load_feature_of_interest(file_name = 'feature_of_interest.txt'):
    features = []
    with open(file_name) as f:
        all_lines = f.readlines()
        for line in all_lines:
            line = line.strip()
            features.append(line)
    return features




def extract_feature(file_name, file_name_ino, feature_list=None):
    df = pd.read_csv(file_name, sep=',', header=0, dtype={'clientguid': str, 'valuetext':str, 'unitofmeasure':str})
    # df = df[df["Observation_name"].isin(["PCO2 ARTERIAL {NYP}", "PCO2 (ARTERIAL) - EPOC {NYP}", "PO2 ARTERIAL {NYP}", "PO2 (ARTERIAL) - EPOC {NYP}"])]
    N0 = df.shape[0]
    df = df.dropna(subset=['clientguid'])
    N = df.shape[0]
    print('Dropped %s rows with empty clientguid, and %s rows remaining...' % (N0-N, N))

    if feature_list == None:
        feature_list = list(set(df['Observation_name'].values))
    print("Generating data for %s features...\n" % len(feature_list))

    data = {}
    for var in feature_list:
        data[var] = {}

    feature_unit = {}

    for idx, row in df.iterrows():
        clientguid, dtm, doc, obs_name, value, unit = row['clientguid'], row['AuthoredDtm'], row['documentname'], \
                                                      row['Observation_name'], row['valuetext'], row['unitofmeasure']
        dtm_timestamp = time.mktime(time.strptime(dtm, "%Y-%m-%d %H:%M:%S"))


        # only focus on features of interest
        if obs_name not in feature_list:
            continue

        # count feature unit overlap
        if obs_name in feature_unit:
            if unit not in feature_unit[obs_name]:
                feature_unit[obs_name].append(unit)
        else:
            feature_unit[obs_name] = [unit]

        # convert string to int, if possible
        if is_number(value) == True:
            value = float(value)

        if is_number(clientguid) == False:
            # print(idx, clientguid, obs_name, value)
            continue

        clientguid = int(clientguid)
        if clientguid not in data[obs_name]:
            data[obs_name][clientguid] = [(dtm, dtm_timestamp, value, unit)]
        else:
            data[obs_name][clientguid].append((dtm, dtm_timestamp, value, unit))

    print("Regular feature extracted...\n")

    # extract Urine data
    if 'io_output_indwell_urine_cath' not in feature_list:
        feature_list.append('io_output_indwell_urine_cath')
        data['io_output_indwell_urine_cath'] = {}
    df_ino = pd.read_csv(file_name_ino, sep=',', header=0, dtype={'clientguid': str, 'OutValue': str, 'BagVolumeUnit': str}).dropna(subset=['OutValue', 'BagVolumeUnit'])  # .head(300)
    df_ino = df_ino[df_ino['Observation_name'] == 'io_output_indwell_urine_cath'][['clientguid', 'AuthoredDtm', 'Observation_name', 'OutValue', 'BagVolumeUnit']]
    for idx, row in df_ino.iterrows():

        clientguid, dtm, observation_name, value, unit = row['clientguid'], row['AuthoredDtm'], row['Observation_name'], row['OutValue'], row['BagVolumeUnit']

        dtm_timestamp = time.mktime(time.strptime(dtm, "%Y-%m-%d %H:%M:%S"))

        # only focus on features of interest
        if observation_name not in feature_list:
            continue

        # count feature unit overlap
        if observation_name in feature_unit:
            if unit not in feature_unit[observation_name]:
                feature_unit[observation_name].append(unit)
        else:
            feature_unit[observation_name] = [unit]

        # convert string to int, if possible
        if is_number(value) == True:
            value = float(value)

        if is_number(clientguid) == False:
            continue

        clientguid = int(clientguid)
        if clientguid not in data['io_output_indwell_urine_cath']:
            data['io_output_indwell_urine_cath'][clientguid] = [(dtm, dtm_timestamp, value, unit)]
        else:
            data['io_output_indwell_urine_cath'][clientguid].append((dtm, dtm_timestamp, value, unit))

    print("Urine feature extracted...\n")

    print('Finished initial generating...\n')

    # order the sequence of each patient by time
    for feature in data:
        for p in data[feature]:
            temp_time = []
            temp_data = []
            for i in range(len(data[feature][p])):
                temp_time.append(data[feature][p][i][1])
                temp_data.append(data[feature][p][i])

            temp_time = np.array(temp_time)
            temp_data = np.array(temp_data)
            order = np.argsort(temp_time)
            temp_time = temp_time[order].tolist()
            temp_data = temp_data[order].tolist()
            data[feature][p] = temp_data

            # print(feature, p, temp_time, temp_data)
            # print(temp_time)
            # print(temp_data)
            # temp_time, temp_data = (list(t) for t in zip(*sorted(zip(temp_time, temp_data), reverse=False)))
            # temp_data = [x for _, x in sorted(zip(temp_time, temp_data), reverse=False)]
            # data[feature][p] = temp_data

    print('Finished ordering by time...\n')



    # Combine "PCO2 ARTERIAL {NYP}" and "PCO2 (ARTERIAL) - EPOC {NYP}" to construct PCO2
    if 'PCO2_combined' not in feature_list:
        feature_list.append('PCO2_combined')
        data['PCO2_combined'] = {}
    PCO2_art = data["PCO2 ARTERIAL {NYP}"]
    PCO2_art_EPOC = data["PCO2 (ARTERIAL) - EPOC {NYP}"]
    patient_union = set(PCO2_art).union(set(PCO2_art_EPOC))
    for p in patient_union:
        if (p in PCO2_art) and (p not in PCO2_art_EPOC):
            patient_PCO2 = PCO2_art[p]
        if (p not in PCO2_art) and (p in PCO2_art_EPOC):
            patient_PCO2 = PCO2_art_EPOC[p]
        if (p in PCO2_art) and (p in PCO2_art_EPOC):
            A = PCO2_art[p]
            B = PCO2_art_EPOC[p]

            patient_PCO2 = []
            i, j = 0, 0
            while i < len(A) and j < len(B):
                if A[i][1] < B[j][1]:
                    patient_PCO2.append(A[i])
                    i += 1
                else:
                    patient_PCO2.append(B[j])
                    j += 1
            patient_PCO2 += A[i:]
            patient_PCO2 += B[j:]
        data['PCO2_combined'][p] = patient_PCO2

    print("Finished combining PCO2...\n")

    # Combine "PO2 ARTERIAL {NYP}" and "PO2 (ARTERIAL) - EPOC {NYP}" to construct PO2
    if 'PO2_combined' not in feature_list:
        feature_list.append('PO2_combined')
        data['PO2_combined'] = {}
    PO2_art = data["PO2 ARTERIAL {NYP}"]
    PO2_art_EPOC = data["PO2 (ARTERIAL) - EPOC {NYP}"]
    patient_union = set(PO2_art).union(set(PO2_art_EPOC))
    for p in patient_union:
        if (p in PO2_art) and (p not in PO2_art_EPOC):
            patient_PO2 = PO2_art[p]
        if (p not in PO2_art) and (p in PO2_art_EPOC):
            patient_PO2 = PO2_art_EPOC[p]
        if (p in PO2_art) and (p in PO2_art_EPOC):
            A = PO2_art[p]
            B = PO2_art_EPOC[p]

            patient_PO2 = []
            i, j = 0, 0
            while i < len(A) and j < len(B):
                if A[i][1] < B[j][1]:
                    patient_PO2.append(A[i])
                    i += 1
                else:
                    patient_PO2.append(B[j])
                    j += 1
            patient_PO2 += A[i:]
            patient_PO2 += B[j:]

        data['PO2_combined'][p] = patient_PO2

    print("Finished combining PO2...\n")




    # save json file
    file_name = file_name.split('/')
    with open(file_name[0] + "/" +file_name[0] + "_data.json", "w") as jsonf:
        json.dump(data, jsonf, indent=4)


    # save feature dictionary
    with open(file_name[0] + "/" + file_name[0] + "_data_feature_dictionary.txt", 'w') as wf:
        for feature in data.keys():
            wf.write(feature)
            wf.write('\n')


    # save features that have multiple unit
    multi_unit_feature = {}
    for var in feature_unit:
        if len(feature_unit[var]) > 1:
            multi_unit_feature[var] = feature_unit[var]
    with open(file_name[0] + "/" + "multi_unit_features.json", "w") as jsonf:
        json.dump(multi_unit_feature, jsonf, indent=4)

    print('Finished saving.')

def main():
    # feature_list = load_feature_of_interest()
    extract_feature(file_name=FILE_NAME, file_name_ino=FILE_NAME_INO)

main()