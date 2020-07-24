"""

Extract event (admission, intubation, extubation, death) time of each patient.

"""

import json
import pandas as pd
import time
import numpy as np


VERSION = 'COVID-0519'


def load_admi_time_from_icu_vists(version=VERSION):
    file_name = version + '/CEDAR/icu_visits.csv'


    event_data = {}

    df = pd.read_csv(file_name, sep=',', header=0, dtype={'ClientGUID':str}).fillna('None')

    df = df.dropna(subset=['ClientGUID'])

    for idx, row in df.iterrows():
        clientguid, admi_time, discharge_time, icu_overall_start_date, icu_overall_end_date = row['ClientGUID'], row['AdmitDtm'], row['DischargeDtm'],\
                                                                                              row['icu_overall_start_date'], row['icu_overall_end_date']

        if clientguid not in event_data:
            event_data[clientguid] = {'AdmitDtm':[], 'DischargeDtm':[], 'icu_overall_start_date':[], 'icu_overall_end_date':[]}
        if admi_time != 'None':
            event_data[clientguid]['AdmitDtm'].append(time.strptime(admi_time, "%Y-%m-%d %H:%M:%S"))
        if discharge_time != 'None':
            event_data[clientguid]['DischargeDtm'].append(time.strptime(discharge_time, "%Y-%m-%d %H:%M:%S"))
        if icu_overall_start_date != 'None':
            event_data[clientguid]['icu_overall_start_date'].append(time.strptime(icu_overall_start_date, "%Y-%m-%d %H:%M:%S"))
        if icu_overall_end_date != 'None':
            event_data[clientguid]['icu_overall_end_date'].append(time.strptime(icu_overall_end_date, "%Y-%m-%d %H:%M:%S"))

    with open(version + "/raw_event.json", 'w') as jsonf:
        json.dump(event_data, jsonf, indent=4)

    for p in event_data:
        if len(event_data[p]['AdmitDtm']) != 0:
            # if min(event_data[p]['AdmitDtm']) != max(event_data[p]['AdmitDtm']):
            #     print(p)
            event_data[p]['AdmitDtm'] = time.strftime("%Y-%m-%d %H:%M:%S", min(event_data[p]['AdmitDtm']))
        else:
            event_data[p]['AdmitDtm'] = ''

        if len(event_data[p]['DischargeDtm']) != 0:
            event_data[p]['DischargeDtm'] = time.strftime("%Y-%m-%d %H:%M:%S", max(event_data[p]['DischargeDtm']))
        else:
            event_data[p]['DischargeDtm'] = ''

        if len(event_data[p]['icu_overall_start_date']) != 0:
            event_data[p]['icu_overall_start_date'] = time.strftime("%Y-%m-%d %H:%M:%S", min(event_data[p]['icu_overall_start_date']))
        else:
            event_data[p]['icu_overall_start_date'] = ''

        if len(event_data[p]['icu_overall_end_date']) != 0:
            event_data[p]['icu_overall_end_date'] = time.strftime("%Y-%m-%d %H:%M:%S", max(event_data[p]['icu_overall_end_date']))
        else:
            event_data[p]['icu_overall_end_date'] = ''

    # print(event_data['9000001147000200'])
    return event_data


def load_admi_time_from_icu_and_ed(version=VERSION):
    event_data = {}


    # icu admission info
    file_name = version + '/CEDAR/icu_visits.csv'
    df = pd.read_csv(file_name, sep=',', header=0, dtype={'ClientGUID': str}).fillna('None')
    df = df.dropna(subset=['ClientGUID'])

    for idx, row in df.iterrows():
        clientguid, admi_time, discharge_time, icu_overall_start_date, icu_overall_end_date = row['ClientGUID'], row[
                            'AdmitDtm'], row['DischargeDtm'],  row['icu_overall_start_date'], row['icu_overall_end_date']

        if clientguid not in event_data:
            event_data[clientguid] = {'AdmitDtm': [], 'DischargeDtm': [], 'icu_overall_start_date': [],
                                      'icu_overall_end_date': []}
        if admi_time != 'None':
            event_data[clientguid]['AdmitDtm'].append(time.strptime(admi_time, "%Y-%m-%d %H:%M:%S"))
        if discharge_time != 'None':
            event_data[clientguid]['DischargeDtm'].append(time.strptime(discharge_time, "%Y-%m-%d %H:%M:%S"))
        if icu_overall_start_date != 'None':
            event_data[clientguid]['icu_overall_start_date'].append(
                time.strptime(icu_overall_start_date, "%Y-%m-%d %H:%M:%S"))
        if icu_overall_end_date != 'None':
            event_data[clientguid]['icu_overall_end_date'].append(
                time.strptime(icu_overall_end_date, "%Y-%m-%d %H:%M:%S"))

    # ed admission info
    file_name = version + '/MBOSS/ed_visits.csv'
    df = pd.read_csv(file_name, sep=',', header=0, dtype={'ClientGUID': str}).fillna('None')
    df = df.dropna(subset=['ClientGUID'])

    for idx, row in df.iterrows():
        clientguid, admi_time, discharge_time = row['ClientGUID'], row['AdmitDtm'], row['DischargeDtm']
        if clientguid not in event_data:
            event_data[clientguid] = {'AdmitDtm': [], 'DischargeDtm': [], 'icu_overall_start_date': [],
                                      'icu_overall_end_date': []}
        if admi_time != 'None':
            event_data[clientguid]['AdmitDtm'].append(time.strptime(admi_time, "%Y-%m-%d %H:%M:%S"))
        if discharge_time != 'None':
            # print(idx, clientguid, discharge_time)
            event_data[clientguid]['DischargeDtm'].append(time.strptime(discharge_time, "%Y-%m-%d %H:%M:%S"))

    with open(version + "/raw_event.json", 'w') as jsonf:
        json.dump(event_data, jsonf, indent=4)

    for p in event_data:
        if len(event_data[p]['AdmitDtm']) != 0:
            # if min(event_data[p]['AdmitDtm']) != max(event_data[p]['AdmitDtm']):
            #     print(p)
            event_data[p]['AdmitDtm'] = time.strftime("%Y-%m-%d %H:%M:%S", min(event_data[p]['AdmitDtm']))
        else:
            event_data[p]['AdmitDtm'] = ''

        if len(event_data[p]['DischargeDtm']) != 0:
            event_data[p]['DischargeDtm'] = time.strftime("%Y-%m-%d %H:%M:%S", max(event_data[p]['DischargeDtm']))
        else:
            event_data[p]['DischargeDtm'] = ''

        if len(event_data[p]['icu_overall_start_date']) != 0:
            event_data[p]['icu_overall_start_date'] = time.strftime("%Y-%m-%d %H:%M:%S", min(event_data[p]['icu_overall_start_date']))
        else:
            event_data[p]['icu_overall_start_date'] = ''

        if len(event_data[p]['icu_overall_end_date']) != 0:
            event_data[p]['icu_overall_end_date'] = time.strftime("%Y-%m-%d %H:%M:%S", max(event_data[p]['icu_overall_end_date']))
        else:
            event_data[p]['icu_overall_end_date'] = ''

    # print(event_data['9000001147000200'])
    return event_data

def load_intubation_time_from_sequence(version=VERSION):
    file_name = version + '/' + version + '_data.json'

    intubation_data = {}

    with open(file_name, 'r') as load_f:
        json_data = json.load(load_f)

        vent_device_data = json_data['resp_non vent device']

        for p in vent_device_data:
            # if p == '9000401012800200':
            #     print(vent_device_data[p])
            # forward search
            vent_start = 'None'
            for i in range(len(vent_device_data[p])):
                if vent_device_data[p][i][2] == 'Ventilator':
                    vent_start = vent_device_data[p][i][0]
                    break

            # backward search
            vent_end = 'None'
            last_timestamp = vent_device_data[p][-1][1]
            for i in range(1, len(vent_device_data[p])+1):
                if vent_device_data[p][-i][1] == last_timestamp:  # at last timestamp
                    if vent_device_data[p][-i][2] == 'Ventilator' or vent_device_data[p][-i][2] == 'Intubated via ETT':  # Ventilator present in the last timestamp
                        break
                    else:
                        continue
                else:  # not at last timestamp
                    if vent_device_data[p][-i][2] == 'Ventilator' or vent_device_data[p][-i][2] == 'Intubated via ETT':
                        vent_end = vent_device_data[p][-i][0]
                        break

            if vent_start == vent_end and vent_start != 'None':
                print(p, vent_start, vent_end)
                vent_end = 'None'
            intubation_data[p] = {'intubation_start':vent_start, 'intubation_end':vent_end}

    return intubation_data


def load_intubation_time_from_sequence_mechvent_invas(version=VERSION):
    file_name = version + '/' + version + '_data.json'

    intubation_data = {}

    with open(file_name, 'r') as load_f:
        json_data = json.load(load_f)

        resp_mechvent_data = json_data['resp_mechvent_invas/noninvas']

        for p in resp_mechvent_data:
            #forward search
            vent_start = 'None'
            for i in range(len(resp_mechvent_data[p])):
                if resp_mechvent_data[p][i][2] == 'Invasive':
                    vent_start = resp_mechvent_data[p][i][0]
                    break
            # backward search
            vent_end = 'None'
            last_timestamp = resp_mechvent_data[p][-1][1]
            for i in range(1, len(resp_mechvent_data[p]) + 1):
                if resp_mechvent_data[p][-i][1] == last_timestamp:  # at last timestamp
                    if resp_mechvent_data[p][-i][2] == 'Invasive':  # Invasive present in the last timestamp
                        break
                    else:
                        continue
                else:  # not at last timestamp
                    if resp_mechvent_data[p][-i][2] == 'Invasive':
                        vent_end = resp_mechvent_data[p][-i][0]
                        break
            intubation_data[p] = {'intubation_start': vent_start, 'intubation_end': vent_end}
    return intubation_data

def load_COVID_status(version=VERSION):
    df = pd.read_csv(version + '/CEDAR/patients.csv', sep=',', header=0, dtype={'clientguid':str})
    df = df[['clientguid', 'ObservationValue']]
    return df

def load_location(version=VERSION):
    df = pd.read_csv(version + '/CEDAR/icu_visits.csv', sep=',', header=0, dtype={'ClientGUID':str})
    patients = list(set(df['ClientGUID'].values))
    location_dic = {}
    for p in patients:
        location_dic[p] = False
    for idx, row in df.iterrows():
        p, name = row['ClientGUID'], row['NAME']
        if ('LMH' in name) or ('lmh' in name) or ('Lmh' in name):
            location_dic[p] = True

    location_df = pd.DataFrame(columns=['clientguid', 'is_LMH'])
    for p in location_dic:
        location_df = location_df.append({'clientguid':p, 'is_LMH':location_dic[p]}, ignore_index=True)

    return location_df


def main():
    # event_data = load_admi_time_from_icu_vists()
    event_data = load_admi_time_from_icu_and_ed()

    intubation_data = load_intubation_time_from_sequence()

    # if '9000401728000200' in intubation_data:
    #     print('yes')


    intubation_mechvent_data = load_intubation_time_from_sequence_mechvent_invas()

    patients = list(event_data.keys() | intubation_data.keys())  # union of patients

    # print(len(event_data), len(intubation_data), len(patients))

    time_colums = ['AdmitDtm', 'DischargeDtm', 'icu_overall_start_date', 'icu_overall_end_date',
                    'intubation_start', 'intubation_end', 'intubation_start_mechvent', 'intubation_end_mechvent']
    event_df = pd.DataFrame({'clientguid':patients})
    event_df = pd.concat([event_df, pd.DataFrame(columns=time_colums)], sort=True)
    event_df = event_df.reindex(columns=['clientguid']+time_colums)


    for idx, row in event_df.iterrows():
        p = row['clientguid']

        if p in event_data:
            event_df.loc[idx, 'AdmitDtm'] = event_data[p]['AdmitDtm']
            event_df.loc[idx, 'DischargeDtm'] = event_data[p]['DischargeDtm']
            event_df.loc[idx, 'icu_overall_start_date'] = event_data[p]['icu_overall_start_date']
            event_df.loc[idx, 'icu_overall_end_date'] = event_data[p]['icu_overall_end_date']
        if p in intubation_data:
            event_df.loc[idx, 'intubation_start'] = intubation_data[p]['intubation_start']
            event_df.loc[idx, 'intubation_end'] = intubation_data[p]['intubation_end']
        if p in intubation_mechvent_data:
            event_df.loc[idx, 'intubation_start_mechvent'] = intubation_mechvent_data[p]['intubation_start']
            event_df.loc[idx, 'intubation_end_mechvent'] = intubation_mechvent_data[p]['intubation_end']



    # filter out patients of non-detected or located at LMH
    COVID_status = load_COVID_status(version=VERSION)

    # COVID_status_2 = load_COVID_status(version='COVID-0501')
    #
    # print(COVID_status)
    # print(COVID_status_2)
    #
    # COVID_status = pd.concat([COVID_status, COVID_status_2])
    # COVID_status = COVID_status.drop_duplicates(subset=['clientguid'], keep='first')
    #
    # print(COVID_status)





    event_df = pd.merge(event_df, COVID_status, how='left', on='clientguid')

    location_df = load_location(version=VERSION)
    event_df = pd.merge(event_df, location_df, how='left', on='clientguid')

    event_df = event_df[event_df['ObservationValue'].isin(['Detected','DETECTED'])]

    LMH_event_df = event_df[event_df['is_LMH'] == True]
    # LMH_event_df.to_csv(VERSION + '/' + VERSION + '_event_LMH.csv', index=False, header=True)


main()
