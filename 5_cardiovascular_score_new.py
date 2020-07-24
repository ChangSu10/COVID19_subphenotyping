# compute cardiovascular SOFA score

import json
import time
import numpy as np
import pandas as pd
import os


VERSION = "COVID-0519"


FILE_NAME_map = VERSION + '/' + VERSION + '_data.json'

FILE_NAME_med = VERSION + '/' + VERSION + '_ino_data.json'

EVENT_FILE = VERSION + '/' + VERSION + '_event.csv'








def judge_map(id_j,start_date_j,end_date_j):
    # map
    start_date = start_date_j
    end_date = end_date_j

    S_value = []
    try:
        map_v = data_map['vs_bp_noninvasive (m)'][id_j]
        for j in range(len(map_v)):
            temp_time = float(map_v[j][1])
            if (temp_time >= float(start_date)) and (temp_time <= float(end_date)):
                S_value.append(float(map_v[j][2]))
    except:
        print('there is no map value.')

    return S_value

def MAP_score(id,start_date,end_date):
    map_values = judge_map(id,start_date,end_date)
    print('map_values', map_values)
    if len(map_values)>0:
        min_map = np.nanmin(map_values)
        if min_map >=70:
            map_score = 0
        else:
            map_score =1
    else:
        map_score = 'None'

    return map_score


def judge_dopamine(id_j, start_date_j, end_date_j):
    # extract dopamine from covid-0414_ino_data, data.med
    start_date = start_date_j
    end_date = end_date_j
    S_InValue = []
    S_ValueDose = []
    try:
        dop_v = data_med['Dopamine'][id_j]
        for j in range(len(dop_v)):
            unit = dop_v[j]['InitialDoseUnit']
            if unit == "UNIT/kg/min":
                continue

            temp_time = float(dop_v[j]['Time_float'])
            BagVolumeUnit = dop_v[j]['BagVolumeUnit']
            if (temp_time >= float(start_date)) and (temp_time <= float(end_date)) and (BagVolumeUnit=='ml'):

                S_InValue.append(float(dop_v[j]['InValue']))
                S_ValueDose.append(float(dop_v[j]['ValueDose']))
    except:
        print('there is no dop value.')

    s_value = S_InValue+S_ValueDose

    return s_value

def Dopamine_score(id_j,start_date_j,end_date_j):
# compute score

    dop_values = judge_dopamine(id_j, start_date_j,end_date_j)
    print('dop_values', dop_values)

    dop_score = 'None'
    if len(dop_values) > 0:
        max_dop = np.nanmax(dop_values)
        if max_dop <=5:
            if max_dop > 0:
                dop_score = 2
        elif (max_dop > 5) and (max_dop <= 15):
            dop_score = 3
        else:
            dop_score = 4
    else:
        dop_score = 'None'

    dopamine_score = dop_score

    return dopamine_score

def judge_dobutamine(id_j,start_date_j,end_date_j):
    # extract dobutamine from covid-0414_ino_data, data.med
    start_date = start_date_j
    end_date = end_date_j
    S_InValue = []
    S_ValueDose = []
    try:
        dop_v = data_med['dobutamine'][id_j]
        for j in range(len(dop_v)):
            temp_time = float(dop_v[j]['Time_float'])
            BagVolumeUnit = dop_v[j]['BagVolumeUnit']
            if (temp_time >= float(start_date)) and (temp_time <= float(end_date)) and (BagVolumeUnit=='ml'):

                S_InValue.append(float(dop_v[j]['InValue']))
                S_ValueDose.append(float(dop_v[j]['ValueDose']))
    except:
        print('there is no dob value.')

    s_value = S_InValue+S_ValueDose

    return s_value

def Dobutamine_score(id_j, start_date_j, end_date_j):
    dob_values = judge_dobutamine(id_j, start_date_j, end_date_j)

    print('dob_values',dob_values)
    if len(dob_values) == 0:
        return 'None'

    if np.nanmin(dob_values) > 0:
        any_dob_score = 2
    else:
        any_dob_score = 'None'

    dobutamine_score = any_dob_score

    return dobutamine_score

def judge_epinephrine(id_j,start_date_j,end_date_j):
    # extract epinephrine from covid-0414_ino_data, data.med
    start_date = start_date_j
    end_date = end_date_j
    S_InValue = []
    S_ValueDose = []
    try:
        dop_v = data_med['Epinephrine'][id_j]
        for j in range(len(dop_v)):
            temp_time = float(dop_v[j]['Time_float'])
            BagVolumeUnit = dop_v[j]['BagVolumeUnit']
            if (temp_time >= float(start_date)) and (temp_time <= float(end_date)) and (BagVolumeUnit=='ml'):
                # print dop_v[j]
                # S_time.append(temp_time)
                S_InValue.append(float(dop_v[j]['InValue']))
                S_ValueDose.append(float(dop_v[j]['ValueDose']))
    except:
        print('there is no epinephrine value.')

    s_value = S_InValue+S_ValueDose

    return s_value

def Epinephrine_score(id, start_date, end_date):
# need weight to compute
    try:
        weight = float(data_weight[id])
    except:
        weight = 0

    epi_values = judge_epinephrine(id, start_date, end_date)

    print('epi_values,without weight',epi_values)

    epi_score = 'None'

    if len(epi_values) > 0 and weight>0:
        max_epi = np.nanmax(epi_values)
        max_epi = max_epi/float(weight)
        if max_epi <= 0.1:
            if max_epi > 0:
                epi_score = 3
        else:
            epi_score = 4
    else:
        epi_score ='None'

    epinephrine_score = epi_score

    return epinephrine_score

def judge_norepinephrine(id_j,start_date_j,end_date_j):
    # extract norepinephrine from covid-0414_ino_data, data.med
    start_date = start_date_j
    end_date = end_date_j
    S_InValue = []
    S_ValueDose = []
    try:
        dop_v = data_med['Norepinephrine'][id_j]
        for j in range(len(dop_v)):
            temp_time = float(dop_v[j]['Time_float'])
            BagVolumeUnit = dop_v[j]['BagVolumeUnit']
            if (temp_time >= float(start_date)) and (temp_time <= float(end_date)) and (BagVolumeUnit=='ml'):
                S_InValue.append(float(dop_v[j]['InValue']))
                S_ValueDose.append(float(dop_v[j]['ValueDose']))
    except:
        print('there is no Norepinephrine value.')

    s_value = S_InValue + S_ValueDose
    print('Norepinephrine', s_value)

    return s_value

def Norepinephrine_score(id, start_date, end_date):
    # need weight to compute
    try:
        weight = float(data_weight[id])
    except:
        weight = 0

    norepi_values = judge_norepinephrine(id, start_date, end_date)

    norepi_score = 'None'
    if len(norepi_values) > 0 and weight > 0:
        max_epi = np.nanmax(norepi_values)
        max_epi = max_epi / float(weight)

        print(norepi_values, np.nanmax(norepi_values), weight, max_epi)

        if max_epi <= 0.1:
            if max_epi > 0:
                norepi_score = 3
        else:
            norepi_score = 4
    else:
        norepi_score = 'None'

    norepinephrine_score = norepi_score

    return norepinephrine_score

def judge_vasopressin(id_j,start_date_j,end_date_j):
    # extract norepinephrine from covid-0414_ino_data, data.med
    start_date = start_date_j
    end_date = end_date_j
    S_InValue = []
    S_ValueDose = []
    try:
        dop_v = data_med['Vasopressin'][id_j]
        # print('dop_v',dop_v)

        for j in range(len(dop_v)):
            unit = dop_v[j]['InitialDoseUnit']
            if unit !='UNIT/hr':
                continue
            temp_time = float(dop_v[j]['Time_float'])
            BagVolumeUnit = dop_v[j]['BagVolumeUnit']
            if (temp_time >= float(start_date)) and (temp_time <= float(end_date)) and (BagVolumeUnit=='ml'):

                S_InValue.append(float(dop_v[j]['InValue']))
                S_ValueDose.append(float(dop_v[j]['ValueDose']))
    except:
        print('there is no Vasopressin value.')

    s_value = S_InValue+S_ValueDose

    print('Vasopressin',s_value)

    return s_value

def Vasopressin_score(id, start_date, end_date):
    # compute Vasopressin
    get_vaso_value = 0
    vaso_score = 'None'
    Vaso_values = judge_vasopressin(id, start_date, end_date)
    if len(Vaso_values) > 0:
        max_vaso = np.nanmax(Vaso_values)
        max_vaso = max_vaso / 60.0 # unit, hour->minute
        if max_vaso > 0:
            vaso_score = 3

        get_vaso_value = max_vaso
    else:
        vaso_score = 'None'

    vasopressin_score = vaso_score

    return vasopressin_score, get_vaso_value

def judge_Phenyephrine(id_j,start_date_j,end_date_j):
    # extract phenyephrine from covid-0414_ino_data, data.med
    start_date = start_date_j
    end_date = end_date_j
    S_InValue = []
    S_ValueDose = []
    try:
        dop_v = data_med['Phenylephrine'][id_j]
        for j in range(len(dop_v)):
            temp_time = float(dop_v[j]['Time_float'])
            BagVolumeUnit = dop_v[j]['BagVolumeUnit']
            if (temp_time >= float(start_date)) and (temp_time <= float(end_date)) and (BagVolumeUnit=='ml'):

                S_InValue.append(float(dop_v[j]['InValue']))
                S_ValueDose.append(float(dop_v[j]['ValueDose']))
    except:
        print('there is no Phenylephrine value.')

    s_value = S_InValue + S_ValueDose

    print('Phenyephrine', s_value)

    return s_value

def Phenyephrine_score(id, start_date, end_date):
    # compute phenyephrine_score
    phen_score = 'None'
    phen_values = judge_Phenyephrine(id, start_date, end_date)
    if len(phen_values) > 0:
        max_phen = np.nanmax(phen_values)
        if max_phen <= 200:
            if max_phen > 0:
                phen_score = 2
        else:
            phen_score = 3
    else:
        phen_score = 'None'

    phenyephrine_score = phen_score

    return phenyephrine_score

def compute_cardiovascular_score(data_map, data_med, data_weight, id, start_date, end_date):
    cardi_scores=[]
# map
    print('map.............')
    car_map_score = MAP_score(id, start_date, end_date)
    print(car_map_score)
    if car_map_score !='None':
        cardi_scores.append(car_map_score)
# dopamine
    print('car_dop_score.............')
    car_dop_score = Dopamine_score(id, start_date, end_date)
    print(car_dop_score)
    if car_dop_score!='None':
        cardi_scores.append(car_dop_score)
# Dobutamine
    print('car_dobu_score.............')
    car_dobu_score = Dobutamine_score(id, start_date, end_date)
    print(car_dobu_score)
    if car_dobu_score !='None':
        cardi_scores.append(car_dobu_score)

# Epinephrine_score
    print('car_epi_score.............')
    car_epi_score = Epinephrine_score(id, start_date, end_date)
    print(car_epi_score)
    if car_epi_score!='None':
        cardi_scores.append(car_epi_score)

# Norepinephrine_score
    print('car_nore_score.............')
    car_nore_score = Norepinephrine_score(id, start_date, end_date)
    print(car_nore_score)
    if car_nore_score !='None':
        cardi_scores.append(car_nore_score)

# Vasopressin_score
    print('Vasopressin_score.............')
    car_vaso_score, vaso_value = Vasopressin_score(id, start_date, end_date) # return two values
    print(car_vaso_score)
    if car_vaso_score!='None':
        cardi_scores.append(car_vaso_score)

# Phenyephrine_score
    print('car_phe_score.............')
    car_phe_score = Phenyephrine_score(id, start_date, end_date)
    print(car_phe_score)
    if car_phe_score !='None':
        cardi_scores.append(car_phe_score)

    print('---', cardi_scores)
    if len(cardi_scores)>0:
        max_car_score = float(np.nanmax(np.array(cardi_scores)))
    else:
        max_car_score = np.nan

    cardiovascular_score = max_car_score
    print('cardiovascular_score.................',cardiovascular_score)

# additional socre, note that car_vaso_score==3.5 , vaso>0.4 units/min
    if (car_nore_score==3 or car_epi_score==3 or car_dop_score==3) and (vaso_value > 0.4 or car_phe_score == 3):

        cardiovascular_score = 4.0

    print('final cardiovascular_score....................',cardiovascular_score)

    return cardiovascular_score




def load_weight(version=VERSION):
    df = pd.read_csv(version + '/CEDAR/icu_visits.csv', sep=',', header=0, dtype={'ClientGUID': str})
    patient_weight = {}
    for idx, row in df.iterrows():
        p, weight = row['ClientGUID'], row['WEIGHT']
        if p in patient_weight:
            patient_weight[p].append(weight)
        else:
            patient_weight[p] = [weight]

    for p in patient_weight:
        patient_weight[p] = np.nanmean(patient_weight[p])

    print(patient_weight)
    return patient_weight


def get_admission_info(event_file, end_point="2020-05-19 00:00:00"):
    patient_admission_info = {}
    df = pd.read_csv(event_file, sep=',', header=0, dtype={'clientguid':str}).fillna('None')
    for idx, row in df.iterrows():
        p, admitDtm, dischargeDtm = row['clientguid'], row['AdmitDtm'], row['DischargeDtm']
        if dischargeDtm == 'None':
            dischargeDtm = end_point
        patient_admission_info[p] = [admitDtm, dischargeDtm]
    return patient_admission_info

def get_intubation_info(event_file, end_point="2020-05-19 00:00:00"):
    patient_intubation_info = {}
    df = pd.read_csv(event_file, sep=',', header=0, dtype={'clientguid':str}).fillna('None')
    for idx, row in df.iterrows():
        p, intubation_start, intubation_end = row['clientguid'], row['intubation_start'], row['intubation_end']
        if intubation_start != 'None':
            if intubation_end == 'None':
                intubation_end = end_point
            patient_intubation_info[p] = [intubation_start, intubation_end]
    return patient_intubation_info






aligned_event = "intubation_start"  # AdmitDtm, intubation_start
interval = 24  
len_before_event = 0
end_point = 60  # None or days after event


with open(FILE_NAME_map, 'r') as load_f:
    data_map = json.load(load_f)
with open(FILE_NAME_med, 'r') as load_f2:
    data_med = json.load(load_f2)


data_weight = load_weight(version=VERSION)
print(data_weight)


if aligned_event == "AdmitDtm":
    data_patient = get_admission_info(event_file=EVENT_FILE)
elif aligned_event == "intubation_start":
    data_patient_intu = get_intubation_info(event_file=EVENT_FILE)
    data_patient = data_patient_intu
    print(data_patient)



pat_list = data_patient.keys()
pat_car_score = {}

for i in pat_list:

    print('i_patient:',i)
    tem_score = []
    pat = i
    AdmitDtm = data_patient[pat][0]
    DischargeDtm = data_patient[pat][1]
    AdmitDtm = time.mktime(time.strptime(AdmitDtm, "%Y-%m-%d %H:%M:%S"))  # time to seconds; start_date_o: "2020-04-09 01:40:00",
    if end_point == None:
        DischargeDtm = time.mktime(time.strptime(DischargeDtm, "%Y-%m-%d %H:%M:%S"))
    else:
        DischargeDtm = AdmitDtm + end_point * 24 * 60 * 60

# admit time
    start_date = AdmitDtm - len_before_event * 24 * 60 * 60  # time before event
# intubation
    while start_date <= DischargeDtm:

        end_date = start_date + interval * 60 * 60
        print(start_date, end_date, '-------')
        score = compute_cardiovascular_score(data_map,data_med,data_weight,pat,start_date,end_date) # 86400 seconds=1 day
        tem_score.append(score)

        start_date = end_date



    pat_car_score[i] = tem_score

sample = pat_car_score

with open(VERSION+'/' + VERSION + '_Cardiovascular_' + aligned_event + "_interval_%s.json" % interval, 'w') as fp:
    json.dump(sample, fp, indent=4)




