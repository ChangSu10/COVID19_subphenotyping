"""

Construct time sequence-like data of each feature, for each patient.

    1. All sequences are alined at the time of intubation;

    2. For each feature, we pick a value every 24-hour bucket, by selecting the worst value within this bucket.

"""

import json
import pandas as pd
import numpy as np

import datetime
import time

VERSION = 'COVID-0519'


def load_event_time(event='AdmitDtm', version=VERSION):
    file_name = version + '/' + version + '_event.csv'

    event_time_dict = {}

    df = pd.read_csv(file_name, sep=',', header=0, dtype={'clientguid':str}).fillna('None')
    for idx, row in df.iterrows():
        clientguid, event_time = row['clientguid'], row[event]
        if event_time == 'None':
            continue
        event_time_dict[clientguid] = event_time

    return event_time_dict


def buckting(aligned_event='AdmitDtm', slide_window=0, max_len=None, interval=24, method='np.mean', version=VERSION):
    """
    :param aligned_event: time point used to align the longitudinal data;
    :param slide_window: allow to capture data within a window before event; default 0 hour
    :param max_len: maximal length (days) of each sequence. Default None, allowing any length
    :param interval: time interval;
    :param method: method used to aggregate data within time window. Options: min, max, np.mean, worse;
    :param version: data version.

    :return: bucked data of patients
    """
    bucked_data = {}

    #--------- load patient sequence data
    load_f = open(VERSION + '/' + version + '_data.json', 'r')
    data = json.load(load_f)
    load_f.close()

    #--------- load unit map info
    load_f = open(version + '/' + '/feaure_units_convert_list.json', 'r')
    unit_map = json.load(load_f)
    load_f.close()

    #--------- load feature of interest
    feature_info_df = pd.read_csv(version + '/' + 'feature_info.csv', sep=',', header=0)

    #--------- load event time dict
    event_time_dict = load_event_time(event=aligned_event, version=version)

    #--------- load feature aggregation method
    if method == 'worse':
        feature_worse_dic = load_feature_worse_method()

    #--------- generate bucked data for the single features
    for idx, row in feature_info_df.iterrows():
        featureName, feature_item, is_continuous, unit = row['FeatureName'], row['Feature_item'], row['Is_continuous'], row['Unit']
        # if feature_item != 'WHITE BLOOD CELL {NYP}':   # WHITE BLOOD CELL {NYP}, D-DIMER {NYP}
        #     continue
        print("Profcessing feature: %s -- %s..." % (featureName, feature_item))
        if is_continuous == 1:
            # initialization
            bucked_data[featureName] = {'Observation_name': feature_item, 'Unit':unit, 'data': {}}
            feature_data = data[feature_item]

            if featureName != 'Urine':
                for p in feature_data:
                    patient_feature_data = feature_data[p]


                    if p in event_time_dict:  # drop the patient if he/she has no event time information

                        patient_event_time = time.mktime(time.strptime(event_time_dict[p], "%Y-%m-%d %H:%M:%S"))


                        # convert sequence info: time and unit
                        patient_feature_converted = []
                        for entry in patient_feature_data:
                            entry_time, entry_value, entry_unit = float(entry[1]), float(entry[2]), entry[3]

                            # unit transformation
                            if feature_item in unit_map and entry_unit in unit_map[feature_item]['outlier_units']:
                                entry_value = entry_value * unit_map[feature_item]['outlier_units'][entry_unit]
                                entry_unit = unit_map[feature_item]['ref_unit'] + '_transformed'

                            # time transformation
                            aligned_time = ((entry_time - patient_event_time) / 60 / 60 + slide_window) / interval

                            patient_feature_converted.append([aligned_time, entry_value, entry_unit])

                        LEN = int(patient_feature_converted[-1][0]) + 1


                        if LEN > 0:
                            if max_len != None:
                                LEN = min(LEN, int(max_len * (24 / interval)))

                            bucked_seq = []
                            for i in range(LEN):
                                bucked_seq.append([])

                            for entry in patient_feature_converted:
                                entry_time, entry_value = entry[0], entry[1]
                                if entry_time < 0:  # ignore records before event
                                    continue
                                else:
                                    buck = int(entry_time)
                                    if max_len != None:
                                        if buck < LEN:
                                            bucked_seq[buck].append(entry_value)
                                    else:
                                        bucked_seq[buck].append(entry_value)


                            # aggregation with specific method
                            if bucked_seq == [[]]:
                                continue

                            # print(bucked_seq)
                            for i in range(len(bucked_seq)):
                                if len(bucked_seq[i]) == 0:
                                    bucked_seq[i] = np.nan
                                else:
                                    if featureName == 'Urine':
                                        bucked_seq[i] = min(bucked_seq[i]) * 2
                                    else:
                                        if method == 'worse':
                                            if featureName in feature_worse_dic:
                                                worse_method = feature_worse_dic[featureName]
                                            else:
                                                worse_method = 'np.mean'
                                            bucked_seq[i] = eval(worse_method)(bucked_seq[i])
                                        else:
                                            bucked_seq[i] = eval(method)(bucked_seq[i])

                            bucked_data[featureName]['data'][p] = bucked_seq

            else:  # addressing Urine
                for p in feature_data:

                    patient_feature_data = feature_data[p]

                    if p in event_time_dict:  # drop the patient if he/she has no event time information
                        patient_event_time = time.mktime(time.strptime(event_time_dict[p], "%Y-%m-%d %H:%M:%S"))

                        # convert sequence info: time and unit
                        patient_feature_converted = []
                        for entry in patient_feature_data:
                            entry_time, entry_value, entry_unit = float(entry[1]), float(entry[2]), entry[3]

                            # unit transformation
                            if feature_item in unit_map and entry_unit in unit_map[feature_item]['outlier_units']:
                                entry_value = entry_value * unit_map[feature_item]['outlier_units'][entry_unit]
                                entry_unit = unit_map[feature_item]['ref_unit'] + '_transformed'

                            # time transformation

                            aligned_time = ((entry_time - patient_event_time) / 60 / 60 + slide_window)

                            patient_feature_converted.append([aligned_time, entry_value, entry_unit])

                        LEN = int(patient_feature_converted[-1][0] / 24) + 1
                        bucked_seq_by_day = []
                        for i in range(LEN):
                            bucked_seq_by_day.append([])

                        if max_len != None:
                            LEN = min(LEN, int(max_len * (24 / interval)))

                        for entry in patient_feature_converted:
                            entry_time, entry_value = entry[0], entry[1]
                            if entry_time < 0:  # ignore records before event
                                continue
                            else:
                                buck = int(entry_time / 24)
                                bucked_seq_by_day[buck].append(entry_value)

                        for i in range(len(bucked_seq_by_day)):
                            if len(bucked_seq_by_day[i]) == 0:
                                bucked_seq_by_day[i] = np.nan
                            else:
                                bucked_seq_by_day[i] = sum(bucked_seq_by_day[i])

                        bucked_seq_by_day[-1] = np.nan


                        interval_buck = int(24 / interval)
                        bucked_seq = []
                        for i in range(len(bucked_seq_by_day)):
                            bucked_seq += [bucked_seq_by_day[i]] * interval_buck

                        bucked_seq = bucked_seq[:LEN*interval_buck]
                        bucked_data[featureName]['data'][p] = bucked_seq




    #--------- generate bucked data for the combined features
    bucked_data = generate_Lymphocyte_count(bucked_data)
    bucked_data = generate_Neutrophil_count(bucked_data)
    bucked_data = generate_Coagulation(bucked_data)
    bucked_data = generate_Liver(bucked_data)
    bucked_data = generate_Central_nervous_system(bucked_data)
    bucked_data = generate_Renal_v2(bucked_data)

    clean_plateau_pressure(bucked_data)
    clean_PEEP(bucked_data)
    clean_resp_peak(bucked_data)
    clean_Tidal_volume(bucked_data)
    clean_Minute_ventilation(bucked_data)

    bucked_data = generate_driving_pressure(bucked_data)
    bucked_data = generate_Static_compliance(bucked_data)
    bucked_data = generate_Tidal_PBW_ratio(bucked_data)
    bucked_data = generate_Ventilator_ratio(bucked_data)

    with open(version + "/bucked_alignedEvent_%s_slide_%s_interval_%s_method_%s.json" % (aligned_event, slide_window, interval, method), 'w') as jsonf:
        json.dump(bucked_data, jsonf, indent=4)


def generate_Lymphocyte_count(bucked_data):
    print("Profcessing feature: Lymphocyte_count...\n")
    if "Lymphocyte" not in bucked_data:
        print("Combination failed. Lymphocyte not included.")
        return bucked_data
    if "WBC" not in bucked_data:
        print("Combination failed. WBC not included.")
        return bucked_data

    Lymphocyte_count_data = {}
    Lymphocyte_data, WBC_data = bucked_data['Lymphocyte']['data'], bucked_data['WBC']['data']
    for p in Lymphocyte_data:
        patient_Lymphocyte = Lymphocyte_data[p]
        if p in WBC_data:
            patient_WBC = WBC_data[p]
            len_L, len_W = len(patient_Lymphocyte), len(patient_WBC)
            len_combined = min(len_L, len_W)
            patient_Lymphocyte, patient_WBC = patient_Lymphocyte[:len_combined], patient_WBC[:len_combined]
            patient_Lymphocyte_count = (np.multiply(patient_Lymphocyte, patient_WBC) / 100).tolist()
            Lymphocyte_count_data[p] = patient_Lymphocyte_count
    bucked_data["Lymphocyte_count"] = {'Observation_name': "Lymphocyte_count", 'Unit':"x10(9)/L", 'data': Lymphocyte_count_data}
    return bucked_data


def generate_Neutrophil_count(bucked_data):
    print("Profcessing feature: Neutrophil_count...\n")
    if "Neutrophil" not in bucked_data:
        print("Combination failed. Neutrophil not included.")
        return bucked_data
    if "WBC" not in bucked_data:
        print("Combination failed. WBC not included.")
        return bucked_data

    Neutrophil_count_data = {}
    Neutrophil_data, WBC_data = bucked_data['Neutrophil']['data'], bucked_data['WBC']['data']
    for p in Neutrophil_data:
        patient_Neutrophil = Neutrophil_data[p]
        if p in WBC_data:
            patient_WBC = WBC_data[p]
            len_N, len_W = len(patient_Neutrophil), len(patient_WBC)
            len_combined = min(len_N, len_W)
            patient_Neutrophil, patient_WBC = patient_Neutrophil[:len_combined], patient_WBC[:len_combined]
            patient_Neutrophil_count = (np.multiply(patient_Neutrophil, patient_WBC) / 100).tolist()
            Neutrophil_count_data[p] = patient_Neutrophil_count
    bucked_data["Neutrophil_count"] = {'Observation_name': "Neutrophil_count", 'Unit': "x10(9)/L", 'data': Neutrophil_count_data}
    return bucked_data



def generate_Coagulation(bucked_data):
    print("Profcessing feature: Coagulation...\n")
    if "Platelet" not in bucked_data:
        print("Combination failed. Platelet not included.")
        return bucked_data

    Coagulation_data = {}
    Platelet_data = bucked_data['Platelet']['data']
    for p in Platelet_data:
        patient_Platelet = Platelet_data[p]
        patient_Coagulation = []
        for i in range(len(patient_Platelet)):
            if np.isnan(patient_Platelet[i]):
                patient_Coagulation.append(np.nan)
            else:
                if patient_Platelet[i] >= 150:
                    patient_Coagulation.append(0.0)
                elif patient_Platelet[i] >= 100:
                    patient_Coagulation.append(1.0)
                elif patient_Platelet[i] >= 50:
                    patient_Coagulation.append(2.0)
                elif patient_Platelet[i] >= 20:
                    patient_Coagulation.append(3.0)
                else:
                    patient_Coagulation.append(4.0)
        Coagulation_data[p] = patient_Coagulation
    bucked_data["Coagulation"] = {'Observation_name': "Coagulation", 'Unit': "", 'data': Coagulation_data}
    return bucked_data

def generate_Liver(bucked_data):
    print("Profcessing feature: Liver...\n")
    if "Bilirubin" not in bucked_data:
        print("Combination failed. Bilirubin not included.")
        return bucked_data

    Liver_data = {}
    Bilirubin_data = bucked_data['Bilirubin']['data']
    for p in Bilirubin_data:
        patient_Bilirubin = Bilirubin_data[p]
        patient_Liver = []
        for i in range(len(patient_Bilirubin)):
            if np.isnan(patient_Bilirubin[i]):
                patient_Liver.append(np.nan)
            else:
                if patient_Bilirubin[i] < 1.2:
                    patient_Liver.append(0.0)
                elif patient_Bilirubin[i] < 2.0:
                    patient_Liver.append(1.0)
                elif patient_Bilirubin[i] < 6.0:
                    patient_Liver.append(2.0)
                elif patient_Bilirubin[i] < 12.0:
                    patient_Liver.append(3.0)
                else:
                    patient_Liver.append(4.0)
        Liver_data[p] = patient_Liver
    bucked_data["Liver"] = {'Observation_name': "Liver", 'Unit': "", 'data': Liver_data}
    return bucked_data

def generate_Central_nervous_system(bucked_data):
    print("Profcessing feature: Central_nervous_system...\n")
    if "GCS" not in bucked_data:
        print("Combination failed. GCS not included.")
        return bucked_data

    Central_nervous_system_data = {}
    GCS_data = bucked_data['GCS']['data']
    for p in GCS_data:

        patient_GCS = GCS_data[p]
        patient_Central_nervous_system = []
        for i in range(len(patient_GCS)):
            if np.isnan(patient_GCS[i]):
                patient_Central_nervous_system.append(np.nan)
            else:
                if patient_GCS[i] >= 15:
                    patient_Central_nervous_system.append(0.0)
                elif patient_GCS[i] >= 13:
                    patient_Central_nervous_system.append(1.0)
                elif patient_GCS[i] >= 10:
                    patient_Central_nervous_system.append(2.0)
                elif patient_GCS[i] >= 6:
                    patient_Central_nervous_system.append(3.0)
                else:
                    patient_Central_nervous_system.append(4.0)
        if p == "9000366926600200":
            print(GCS_data[p])
            print(patient_Central_nervous_system)
        Central_nervous_system_data[p] = patient_Central_nervous_system
    bucked_data["Central_nervous_system"] = {'Observation_name': "Central_nervous_system", 'Unit': "", 'data': Central_nervous_system_data}
    return bucked_data

def generate_Renal(bucked_data):
    print("Profcessing feature: Renal...\n")
    if "Creatinine" not in bucked_data:
        print("Combination failed. Creatinine not included.")
        return bucked_data

    Renal_data = {}
    Creatinine_data = bucked_data['Creatinine']['data']

    for p in Creatinine_data:
        patient_Creatinine = Creatinine_data[p]
        patient_Renal = []
        for i in range(len(patient_Creatinine)):
            if np.isnan(patient_Creatinine[i]):
                patient_Renal.append(np.nan)
            else:
                if patient_Creatinine[i] < 1.2:
                    patient_Renal.append(0.0)
                elif patient_Creatinine[i] < 1.9:
                    patient_Renal.append(1.0)
                elif patient_Creatinine[i] < 3.4:
                    patient_Renal.append(2.0)
                elif patient_Creatinine[i] < 4.9:
                    patient_Renal.append(3.0)
                else:
                    patient_Renal.append(4.0)
        Renal_data[p] = patient_Renal
    bucked_data["Renal"] = {'Observation_name': "Renal", 'Unit': "", 'data': Renal_data}
    return bucked_data

def generate_Renal_v2(bucked_data):
    print("Profcessing feature: Renal...\n")
    if "Creatinine" not in bucked_data:
        print("Combination failed. Creatinine not included.")
        return bucked_data

    Renal_data = {}
    Creatinine_data = bucked_data['Creatinine']['data']
    Urine_data = bucked_data['Urine']['data']


    for p in Creatinine_data:
        patient_Creatinine = Creatinine_data[p]
        patient_Renal = []
        for i in range(len(patient_Creatinine)):
            if np.isnan(patient_Creatinine[i]):
                patient_Renal.append(np.nan)
            else:
                if patient_Creatinine[i] < 1.2:
                    patient_Renal.append(0.0)
                elif patient_Creatinine[i] <= 1.9:
                    patient_Renal.append(1.0)
                elif patient_Creatinine[i] <= 3.4:
                    patient_Renal.append(2.0)
                elif patient_Creatinine[i] <= 4.9:
                    patient_Renal.append(3.0)
                else:
                    patient_Renal.append(4.0)
        Renal_data[p] = {'Creatinine':patient_Renal}

    for p in Urine_data:
        patient_Urine = Urine_data[p]
        patient_Renal = []
        for i in range(len(patient_Urine)):
            if np.isnan(patient_Urine[i]):
                patient_Renal.append(np.nan)
            else:
                if patient_Urine[i] >= 500:
                    patient_Renal.append(np.nan)
                elif patient_Urine[i] >= 200:
                    patient_Renal.append(3.0)
                else:
                    patient_Renal.append(4.0)
        if p in Renal_data:
            Renal_data[p]['Urine'] = patient_Renal
        else:
            Renal_data[p] = {'Urine': patient_Renal}

    for p in Renal_data:
        patient_Renal = Renal_data[p]
        if len(patient_Renal) == 1:
            if 'Creatinine' in patient_Renal:
                patient_Renal = patient_Renal['Creatinine']
            if 'Urine' in patient_Renal:
                patient_Renal = patient_Renal['Urine']
            Renal_data[p] = patient_Renal
        else:
            creatinine = patient_Renal['Creatinine']
            urine = patient_Renal['Urine']
            L = max(len(creatinine), len(urine))
            if len(creatinine) < L:
                creatinine += [np.nan] * (L - len(creatinine))
            if len(urine) < L:
                urine += [np.nan] * (L - len(urine))
            combination = []

            for i in range(L):
                combination.append(np.nanmax([creatinine[i], urine[i]]))

            Renal_data[p] = combination

    bucked_data["Renal"] = {'Observation_name': "Renal", 'Unit': "", 'data': Renal_data}
    return bucked_data


def clean_plateau_pressure(bucked_data):
    print("Processing feature: Plateau_pressure by removing abnormal values...\n")
    if "Plateau_pressure" not in bucked_data:
        print("Processing failed. Plateau_pressure not included.")
        return bucked_data

    plateau_pressure_data = bucked_data['Plateau_pressure']['data']
    for p in plateau_pressure_data:

        for i in range(len(plateau_pressure_data[p])):

            if plateau_pressure_data[p][i] < 10:
                plateau_pressure_data[p][i] = np.nan
            if plateau_pressure_data[p][i] > 50:
                plateau_pressure_data[p][i] = np.nan
    bucked_data['Plateau_pressure']['data'] = plateau_pressure_data
    return bucked_data

def clean_PEEP(bucked_data):
    print("Processing feature: PEEP by removing abnormal values...\n")
    if "PEEP" not in bucked_data:
        print("Processing failed. PEEP not included.")
        return bucked_data

    PEEP_data = bucked_data['PEEP']['data']
    for p in PEEP_data:

        for i in range(len(PEEP_data[p])):
            if PEEP_data[p][i] > 30:
                PEEP_data[p][i] = np.nan
    bucked_data['PEEP']['data'] = PEEP_data
    return bucked_data

def clean_resp_peak(bucked_data):
    print("Processing feature: resp_peak_insp_pres by removing abnormal values...\n")
    if "resp_peak_insp_pres" not in bucked_data:
        print("Processing failed. resp_peak_insp_pres not included.")
        return bucked_data

    resp_peak_data = bucked_data['resp_peak_insp_pres']['data']
    for p in resp_peak_data:

        for i in range(len(resp_peak_data[p])):
            if resp_peak_data[p][i] < 10:
                resp_peak_data[p][i] = np.nan
    bucked_data['resp_peak_insp_pres']['data'] = resp_peak_data
    return bucked_data

def clean_Tidal_volume(bucked_data):
    print("Processing feature: Tidal_volume by removing abnormal values...\n")
    if "Tidal_volume" not in bucked_data:
        print("Processing failed. Tidal_volume not included.")
        return bucked_data

    Tidal_volume_data = bucked_data['Tidal_volume']['data']
    for p in Tidal_volume_data:

        for i in range(len(Tidal_volume_data[p])):
            if Tidal_volume_data[p][i] > 800:
                Tidal_volume_data[p][i] = np.nan
    bucked_data['Tidal_volume']['data'] = Tidal_volume_data
    return bucked_data

def clean_Minute_ventilation(bucked_data):
    print("Processing feature: Minute_ventilation by removing abnormal values...\n")
    if "Minute_ventilation" not in bucked_data:
        print("Processing failed. Minute_ventilation not included.")
        return bucked_data

    Minute_ventilation_data = bucked_data['Minute_ventilation']['data']
    for p in Minute_ventilation_data:

        for i in range(len(Minute_ventilation_data[p])):
            if Minute_ventilation_data[p][i] > 25:
                Minute_ventilation_data[p][i] = np.nan
    bucked_data['Minute_ventilation']['data'] = Minute_ventilation_data
    return bucked_data

def generate_driving_pressure(bucked_data):
    print("Processing feature: Driving_pressure...\n")
    if "Plateau_pressure" not in bucked_data:
        print("Combination failed. Plateau_pressure not included.")
        return bucked_data
    if "PEEP" not in bucked_data:
        print("Combination failed. PEEP not included.")
        return bucked_data

    driving_pressure_data = {}
    Plateau_pressure_data = bucked_data['Plateau_pressure']['data']
    PEEP_data = bucked_data['PEEP']['data']

    for p in Plateau_pressure_data:
        if p in PEEP_data:
            patient_Plateau_pressure = Plateau_pressure_data[p]
            patient_PEEP = PEEP_data[p]
            patient_driving_pressure = []
            for i in range(min(len(patient_Plateau_pressure), len(patient_PEEP))):
                if patient_Plateau_pressure[i] >= patient_PEEP[i]:
                    patient_driving_pressure.append(patient_Plateau_pressure[i] - patient_PEEP[i])
                else:
                    patient_driving_pressure.append(np.nan)
            driving_pressure_data[p] = patient_driving_pressure
    bucked_data["Driving_pressure"] = {'Observation_name': "Driving_pressure", 'Unit': "",
                                             'data': driving_pressure_data}
    return bucked_data

def generate_Static_compliance(bucked_data):
    print("Processing feature: Static_compliance...\n")
    if "Plateau_pressure" not in bucked_data:
        print("Combination failed. Plateau_pressure not included.")
        return bucked_data
    if "PEEP" not in bucked_data:
        print("Combination failed. PEEP not included.")
        return bucked_data
    if "Tidal_volume" not in bucked_data:
        print("Combination failed. Tidal_volume not included.")
        return bucked_data
    if "resp_peak_insp_pres" not in bucked_data:
        print("Combination failed. resp_peak_insp_pres not included.")
        return bucked_data

    Static_compliance_data = {}
    Plateau_pressure_data = bucked_data['Plateau_pressure']['data']
    PEEP_data = bucked_data['PEEP']['data']
    Tidal_volume_data = bucked_data['Tidal_volume']['data']
    resp_peak_insp_pres_data = bucked_data['resp_peak_insp_pres']['data']

    for p in Tidal_volume_data:
        patient_Tidal_volume = Tidal_volume_data[p]
        L = len(patient_Tidal_volume)
        if p in PEEP_data:
            patient_PEEP = PEEP_data[p]
            patient_PEEP += [np.nan] * (L - len(patient_PEEP))
            if p in Plateau_pressure_data:
                patient_Plateau_pressure = Plateau_pressure_data[p]
                patient_Plateau_pressure += [np.nan] * (L - len(patient_Plateau_pressure))
            else:
                patient_Plateau_pressure = [np.nan] * L

            if p in resp_peak_insp_pres_data:
                pateint_resp_peak = resp_peak_insp_pres_data[p]
                pateint_resp_peak += [np.nan] * (L - len(pateint_resp_peak))
            else:
                pateint_resp_peak = [np.nan] * L
            

            patient_Cstat = []
            for i in range(L):
                Cstat = np.nan
                if patient_Plateau_pressure[i] > patient_PEEP[i]:
                    Cstat = patient_Tidal_volume[i] / (patient_Plateau_pressure[i] - patient_PEEP[i])
                elif pateint_resp_peak[i] > patient_PEEP[i]:
                    Cstat = patient_Tidal_volume[i] / (pateint_resp_peak[i] - patient_PEEP[i])

                if Cstat > 100:
                    Cstat = np.nan
                patient_Cstat.append(Cstat)


            Static_compliance_data[p] = patient_Cstat

    bucked_data["Static_compliance"] = {'Observation_name': "Static_compliance", 'Unit': "",
                                       'data': Static_compliance_data}
    return bucked_data

def generate_Tidal_PBW_ratio(bucked_data):
    print("Processing feature: Tidal_PBW_ratio...\n")
    if "Tidal_volume" not in bucked_data:
        print("Combination failed. Tidal_volume not included.")
        return bucked_data

    PBW_data = load_height()

    Tidal_volume_data = bucked_data['Tidal_volume']['data']

    Tidal_PBW_ratio_data = {}
    for p in Tidal_volume_data:
        if p not in PBW_data:
            continue

        patient_Tidal_volume = Tidal_volume_data[p]

        patient_Tidal_PBW_ratio = []
        for i in range(len(patient_Tidal_volume)):
            Tidal_PBW_ratio = patient_Tidal_volume[i] / PBW_data[p]
            patient_Tidal_PBW_ratio.append(Tidal_PBW_ratio)


        Tidal_PBW_ratio_data[p] = patient_Tidal_PBW_ratio

    bucked_data["Tidal_PBW_ratio"] = {'Observation_name': "Tidal_PBW_ratio", 'Unit': "ml/KG",
                                        'data': Tidal_PBW_ratio_data}
    return bucked_data

def load_height(version=VERSION):
    df = pd.read_csv(version + "/CEDAR/icu_visits.csv", sep=',', header=0, dtype={"ClientGUID": str}).dropna(subset=["ClientGUID", "HEIGHT"])
    demographic_df = pd.read_csv(version + "/CEDAR/patients.csv", sep=',', header=0, dtype={"clientguid": str}).dropna(subset=["clientguid", "sex"])

    df = pd.merge(df, demographic_df, how='left', left_on='ClientGUID', right_on='clientguid')

    PBW_data = {}
    for idx, row in df.iterrows():
        p, height, sex, q = row['ClientGUID'], row['HEIGHT'], row['sex'], row['clientguid']
        if height <= 100 or height > 230:
            continue
        if p not in PBW_data:
            if sex == 'M':
                PBW_data[p] = 50 + 0.91 * (height - 152.4)
            if sex == 'F':
                PBW_data[p] = 45.5 + 0.91 * (height - 152.4)
    return PBW_data

def generate_Ventilator_ratio(bucked_data):
    print("Processing feature: Ventilator_ratio by removing abnormal values...\n")
    if "Minute_ventilation" not in bucked_data:
        print("Processing failed. Minute_ventilation not included.")
        return bucked_data

    if "PCO2" not in bucked_data:
        print("Processing failed. PCO2 not included.")
        return bucked_data

    Minute_ventilation_data = bucked_data['Minute_ventilation']['data']
    PCO2_data = bucked_data['PCO2']['data']

    PBW_data = load_height()


    Ventilator_ratio_data = {}
    for p in Minute_ventilation_data:
        if (p in PCO2_data) and (p in PBW_data):
            patient_Ventilator_ratio = []

            patient_Minute_ventilation = Minute_ventilation_data[p]
            patient_PCO2 = PCO2_data[p]

            for i in range(min(len(patient_Minute_ventilation), len(patient_PCO2))):

                VR = patient_Minute_ventilation[i] * 1000 * patient_PCO2[i] / (PBW_data[p] * 100 * 37.5)
                patient_Ventilator_ratio.append(VR)

            Ventilator_ratio_data[p] = patient_Ventilator_ratio

    bucked_data["Ventilator_ratio"] = {'Observation_name': "Ventilator_ratio", 'Unit': "",
                                      'data': Ventilator_ratio_data}
    return bucked_data





def load_feature_worse_method():
    feature_worse_method = {}
    df = pd.read_csv("sofa_feature_worse.csv", sep=',', header=0)
    for idx, row in df.iterrows():
        feature, method = row['FeatureName'], row['method']
        feature_worse_method[feature] = method
    return feature_worse_method

def main():
    buckting(aligned_event='intubation_start', slide_window=0, max_len=None, interval=24, method='worse', version=VERSION)


main()



