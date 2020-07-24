# compute respiration SOFA score

import json
import time
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import datetime
from math import isnan

VERSION = "COVID-0519"


FILE_NAME = VERSION + '/' + VERSION + '_data.json'

EVENT_FILE = VERSION + '/' + VERSION + '_event.csv'


def obtain_ven_status(id_o, start_date_o, end_date_o):
    mech_invas = 0
    vent_dev = 0
    non_rebreather = 0
    OptiFlow = 0
    Nasal_mask = 0
    combine_NRB_Opti = 0
    combine_NRB_NC = 0

    try:
        non_vent_divice = data['resp_non vent device'][id_o]

        for i in range(len(non_vent_divice)):
            temp_time = float(non_vent_divice[i][1])
            if (temp_time >= float(start_date_o)) and (temp_time <= float(end_date_o)):
                con = non_vent_divice[i][2]

                if ('Non-Rebreather Mask' in con) or ('Non- Rebreather Mask'in con) or ('NRB' in con) or ('Non rebreather'in con) or ('Nonrebreather'in con) or ('Non-rebreather'in con):
                    non_rebreather = 1

                if ('High Flow Nasal Cannula'in con) or ('High flow NC' in con) or ('HFNC' in con) or('HF NC' in con):
                    OptiFlow =1

                if ('Nasal Cannula' in con) or ('nasal cannula' in con) or ('NC' in con): # 'Nasal Cannula, and Venti-Mask' use same flag, Nasal_mask
                    Nasal_mask =1

                if 'Venti-Mask' in con:
                    Nasal_mask = 1

                if (('Non-Rebreather Mask' in con) or ('Non- Rebreather Mask' in con) or ('NRB' in con) or ('Non rebreather' in con) or ('Nonrebreather' in con) or ('Non-rebreather' in con)) and (('High Flow Nasal Cannula' in con) or ('High flow NC' in con) or ('HFNC' in con) or ('HF NC' in con)):
                    combine_NRB_Opti =1

                if (('Non-Rebreather Mask' in con) or ('Non- Rebreather Mask' in con) or ('NRB' in con) or ('Non rebreather' in con) or ('Nonrebreather' in con) or ('Non-rebreather' in con)) and (('Nasal Cannula' in con) or ('nasal cannula' in con) or ('NC'in con)):
                    combine_NRB_NC =1

        if (non_rebreather + Nasal_mask + OptiFlow) >0:
            vent_dev = 1

    except:

        vent_dev = 0

# resp_mechvent_invas/noninvas
    try:
        in_noninvas = data['resp_mechvent_invas/noninvas'][id_o]
        for j in range(len(in_noninvas)):
            temp_time_2 = float(in_noninvas[j][1])
            if (temp_time_2 >= float(start_date_o)) and (temp_time_2 <= float(end_date_o)):
                mech_invas = 1
        

    except:

        mech_invas = 0

    return mech_invas,vent_dev,non_rebreather,OptiFlow, Nasal_mask,combine_NRB_Opti,combine_NRB_NC

def judge_ABG(id_j,start_date_j,end_date_j):
    # PO2
    start_date =start_date_j
    end_date = end_date_j
    p_time = []
    p_value = []

    try:
        PO2 = data['PO2 ARTERIAL {NYP}'][id_j]

        for j in range(len(PO2)):
            temp_time = float(PO2[j][1])
            if (temp_time >= float(start_date)) and (temp_time <= float(end_date)):
                p_time.append(temp_time)
                p_value.append(float(PO2[j][2]))
    except:
        print('there is no PO2')
    return p_time, p_value

def judge_Spo2(id_j,start_date_j,end_date_j):
    # SpO2
    start_date =start_date_j
    end_date = end_date_j
    S_time = []
    S_value = []
    try:
        SPO2 = data['xp_resp_spo2'][id_j]
        for j in range(len(SPO2)):
            temp_time = float(SPO2[j][1])
            if (temp_time >= float(start_date)) and (temp_time <= float(end_date)):
                S_time.append(temp_time)
                S_value.append(float(SPO2[j][2]))
    except:
        print('there is no SPO2')

    return S_time, S_value

def judge_Fio2(id_j,start_date_j,end_date_j):
    # FiO2
    start_date_o = start_date_j
    end_date_o = end_date_j
    F_time = []
    F_value = []
    try:
        FIO2 = data['resp_fio2'][id_j]
        for j in range(len(FIO2)):
            temp_time = float(FIO2[j][1])
            if (temp_time >= float(start_date_o)) and (temp_time <= float(end_date_o)):
                F_time.append(temp_time)
                F_value.append(float(FIO2[j][2]))
    except:
        print('there is no FIO2')

    return F_time, F_value

def judge_Pco2(id_j,start_date_j,end_date_j):
    # PCO2
    start_date =start_date_j
    end_date = end_date_j
    p_time = []
    p_value = []
    try:
        PCO2 = data['PCO2 ARTERIAL {NYP}'][id_j]
        for j in range(len(PCO2)):
            temp_time = float(PCO2[j][1])
            if (temp_time >= float(start_date)) and (temp_time <= float(end_date)):
                p_time.append(temp_time)
                p_value.append(float(PCO2[j][2]))
    except:
        print('there is no PCO2')
    return p_time, p_value

def judge_oxygen(id_j,start_date_j,end_date_j):
    # oxygen
    start_date_o = start_date_j
    end_date_o = end_date_j
    oxygen_time = []
    oxygen_value = []

    try:
        oxygen = data['resp_insp_gas_label_oxygen'][id_j]
        for j in range(len(oxygen)):
            temp_time = float(oxygen[j][1])
            if (temp_time >= float(start_date_o)) and (temp_time <= float(end_date_o)):
                oxygen_time.append(temp_time)
                oxygen_value.append(float(oxygen[j][2]))
    except:
        print('there is no oxygen')

    return oxygen_time, oxygen_value

def select_near_fio2(Fi_time,Fi_value,min_time):
    dif_time_max = 24*60*60  # seconds, time difference betwee FIO2 and spo2 or pao2 no more than 1 day
    nearest_time_index = 0

    for i in range(len(Fi_time)):
        dif_time = abs(Fi_time[i]-min_time)
        if dif_time<dif_time_max:
            dif_time_max = dif_time
            nearest_time_index = i
    nearest_value = Fi_value[nearest_time_index]

    return nearest_value

def multi_select_ratio(ABG_Spo2_value,ABG_Spo2_time,Fio2_time,Fio2_value):
    temp_ratio =[]
    ABG_value = ABG_Spo2_value
    ABG_time = ABG_Spo2_time

    pao2_invas = min(np.array(ABG_value))
    min_pao2_time_invas = ABG_time[np.argmin(ABG_value)]

    fio2_invas = select_near_fio2(Fio2_time, Fio2_value, min_pao2_time_invas)

    temp_ratio.append(pao2_invas / fio2_invas)

    start_index = np.argmin(ABG_value)
    for i in range(len(ABG_value)):
        if start_index < len(ABG_value) - 1:
            min_pao2_time_invas_2 = ABG_time[start_index + 1]
            pao2_invas_2 = ABG_value[start_index + 1]
            fio2_invas_2 = select_near_fio2(Fio2_time, Fio2_value, min_pao2_time_invas_2)

            temp_ratio.append(pao2_invas_2 / fio2_invas_2)
            start_index = start_index+1

    return min(np.array(temp_ratio))

def apache_ABG_worse_select(ABG_value,ABG_time, Fio2_value, Fio2_time, Pco2_value, Pco2_time):
    weighe_Aado2 = 0
    weighe_Pao2 = 0
    # if intubation_flag ==0: # 0, non intubation, 1 intubation
    #     pao2 = min(np.array(ABG_value))
    #     min_pao2_time = ABG_time[np.argmin(ABG_value)]
    #     final_pao2 = pao2
    # else:
    if len(Fio2_value)>0 and max(np.array(Fio2_value))>=0.5: #
        if len(ABG_value)>0 and len(Pco2_value)>0:
            Aado2 = max(np.array(Fio2_value))*713-min(np.array(Pco2_value))-min(np.array(ABG_value))
        else:
            Aado2 = 0 # no ABG or Pco2

        if Aado2<100:
            weighe_Aado2 = 0
        elif (Aado2>=100) and (Aado2<250):
            weighe_Aado2 =7
        elif (Aado2>=250) and (Aado2<350):
            weighe_Aado2 =9
        elif (Aado2 >= 350) and (Aado2 < 500):
            weighe_Aado2 = 11
        else:
            weighe_Aado2 =14
    else:
        pao2 = min(np.array(ABG_value))
        if pao2<50:
            weighe_Pao2 = 15
        elif (pao2>=50) and (pao2<70):
            weighe_Pao2 =5
        elif (pao2>=70) and (pao2<80):
            weighe_Pao2 =2
        else:
            weighe_Pao2 =0

    if weighe_Aado2 > weighe_Pao2:
        final_pao2 = Aado2
    else:
        pao2 = min(np.array(ABG_value))
        final_pao2 = pao2

    pao2 = final_pao2
    min_pao2_time = ABG_time[np.argmin(ABG_value)]

    return pao2,min_pao2_time

def multi_select_ratio_Pao2(ABG_Spo2_value,ABG_Spo2_time,Fio2_value,Fio2_time, Pco2_value,Pco2_time):
    temp_ratio =[]
    ABG_value = ABG_Spo2_value
    ABG_time = ABG_Spo2_time

    pao2_invas, min_pao2_time_invas = apache_ABG_worse_select(ABG_Spo2_value, ABG_Spo2_time, Fio2_value, Fio2_time, Pco2_value, Pco2_time)

    fio2_invas = select_near_fio2(Fio2_time, Fio2_value, min_pao2_time_invas)

    temp_ratio.append(pao2_invas / fio2_invas)

    start_index = np.argmin(ABG_value)
    for i in range(len(ABG_value)):
        if start_index < len(ABG_value) - 1:
            min_pao2_time_invas_2 = ABG_time[start_index + 1]
            pao2_invas_2 = ABG_value[start_index + 1]
            fio2_invas_2 = select_near_fio2(Fio2_time, Fio2_value, min_pao2_time_invas_2)

            temp_ratio.append(pao2_invas_2 / fio2_invas_2)
            start_index = start_index+1

    return min(np.array(temp_ratio))

def select_near_ox(ox_time,ox_value,near_Po2_spo2_time):
    min_time = near_Po2_spo2_time
    dif_time_max = 24*60*60  # seconds, time difference betwee FIO2 and spo2 or pao2 no more than 1 day
    nearest_time_index = 0

    for i in range(len(ox_time)):
        dif_time = abs(ox_time[i]-min_time)
        if dif_time<dif_time_max:
            dif_time_max = dif_time
            nearest_time_index = i

    nearest_time = ox_time[nearest_time_index-1]
    nearest_value = ox_value[nearest_time_index-1]

    return nearest_time, nearest_value

def judge_room_air(id_o,start_date,end_date):
    start_date_o = start_date
    end_date_o = end_date
    tem_room_time =[]
    try:
        non_vent_divice = data['resp_non vent device'][id_o]
        for i in range(len(non_vent_divice)):
            temp_time = float(non_vent_divice[i][1])
            if (temp_time>= float(start_date_o)) and (temp_time<=float(end_date_o)):
                con = non_vent_divice[i][2]
                if 'Room Air' in con:
                    tem_room_time.append(temp_time)
    except:
        print('there is no room air:')


    return tem_room_time

def apache_estimate_fio2(oxygen_time,oxygen_value):

    
    if len(oxygen_time)>0:
        nearest_ox_time = oxygen_time[0]
        nearest_ox_value = np.nanmean(np.array(oxygen_value))
        ox_volumn = nearest_ox_value #((end_date - nearest_ox_time) / 60) * nearest_ox_value
    else:
        ox_volumn = 0 # room air

    if ox_volumn ==0:
        ox_fio2 = 0.21
    elif ox_volumn >0 and ox_volumn<=1:
        ox_fio2=0.23
    elif ox_volumn>1 and ox_volumn<=2:
        ox_fio2 = 0.25
    elif ox_volumn>2 and ox_volumn<=3:
        ox_fio2 = 0.27
    elif ox_volumn>3 and ox_volumn<=4:
        ox_fio2 = 0.30
    elif ox_volumn>4 and ox_volumn<=5:
        ox_fio2 = 0.35
    elif ox_volumn>5 and ox_volumn<=7:
        ox_fio2 = 0.40
    elif ox_volumn>7and ox_volumn<=10:
        ox_fio2 = 0.49
    else:
        ox_fio2 = 1.0

    estimated_fio2 = ox_fio2
    return estimated_fio2

def compute_respiration_score(data,id,start_date,end_date):

    mech_invas_status,vent_dev_status, Non_rebreather, OptiFlow, Nasal_Can_Venti,combine_NRB_Opti,combine_NRB_NC = obtain_ven_status(id, start_date, end_date)

    if mech_invas_status == 0 and vent_dev_status == 0 and Non_rebreather == 0 and OptiFlow == 0 and Nasal_Can_Venti == 0 and combine_NRB_Opti == 0 and combine_NRB_NC == 0:
        print("#############################")
        print(id)

    ABG_value = Spo2_value = Fio2_value = oxygen_value = []
    ABG_time = Spo2_time = Fio2_time = oxygen_time = room_time = []

    ABG_time, ABG_value = judge_ABG(id, start_date, end_date)
    Spo2_time, Spo2_value = judge_Spo2(id, start_date, end_date)
    Fio2_time, Fio2_value = judge_Fio2(id, start_date, end_date)
    Pco2_time, Pco2_value = judge_Pco2(id, start_date, end_date)

    oxygen_time, oxygen_value = judge_oxygen(id,start_date,end_date)
    room_time = judge_room_air(id, start_date,end_date)

    print('ABG', ABG_time, '\n', ABG_value)
    print('Spo2', Spo2_time, '\n', Spo2_value)
    print('Fio2', Fio2_time, '\n', Fio2_value)
    print('Pco2', Pco2_time, '\n', Pco2_value)

    print(mech_invas_status,vent_dev_status, Non_rebreather, OptiFlow, Nasal_Can_Venti,combine_NRB_Opti,combine_NRB_NC)


    if len(Fio2_value)>0: # Fio2 unit %
        Fio2_value = np.true_divide(np.array(Fio2_value), 100)
        Fio2_value = list(Fio2_value)

    pao2 = 'None'
    spo2 = 'None'
    fio2 = 'None'

    ratios = []

    if vent_dev_status == 1: # 1 supplenment O2;
        if len(ABG_value)>0:
        # non apache method
            pao2= min(np.array(ABG_value))
            min_pao2_time = ABG_time[np.argmin(ABG_value)]
        # apache method for choosing min pao2 ,considering inbutation
        #     pao2,min_pao2_time = apache_ABG_worse_select(intubation_flag,ABG_time, ABG_value,Fio2_time, Fio2_value,Pco2_time, Pco2_value)
        else:
            if len(Spo2_value)>0:
                spo2 = min(np.array(Spo2_value))
                min_spo2_time = Spo2_time[np.argmin(Spo2_value)]

        if OptiFlow ==1:
            fio2_OF = 'None'
            if len(ABG_value)>0:
                if len(Fio2_value)>0:
                    fio2_OF = select_near_fio2(Fio2_time, Fio2_value, min_pao2_time)
            else:
                if len(Spo2_value)>0 and len(Fio2_value)>0:
                    fio2_OF = select_near_fio2(Fio2_time, Fio2_value, min_spo2_time)


        if Nasal_Can_Venti ==1:
            fio2_NC = apache_estimate_fio2(oxygen_time,oxygen_value)

        if len(ABG_value) > 0 and Non_rebreather == 1: #if Non_rebreather ==1: fio2_rebreather = 0.65 # 65%
            ratios.append(pao2/0.65)
        if len(ABG_value) > 0 and combine_NRB_Opti == 1: # if NRB and OptiFlow , fio2 = 1
            ratios.append(pao2 / 1.0)
        if len(ABG_value) > 0 and combine_NRB_NC == 1:  #if NRB and NC , fio2 = 0.80
            ratios.append(pao2 / 0.80)

        if len(ABG_value)>0 and OptiFlow==1 and fio2_OF!='None':
            ratios.append(pao2/fio2_OF)
        if len(ABG_value)>0 and Nasal_Can_Venti==1:
            ratios.append(pao2/fio2_NC)

        if len(ABG_value)==0 and len(Spo2_value)>0 and Non_rebreather == 1:
            ratios.append(((spo2/0.65)*1.19)-76.19)
        if len(ABG_value)==0 and len(Spo2_value)>0 and combine_NRB_Opti == 1:
            ratios.append(((spo2/1.0)*1.19)-76.19)
        if len(ABG_value)==0 and len(Spo2_value)>0 and combine_NRB_NC == 1:
            ratios.append(((spo2/0.80)*1.19)-76.19)

        if len(ABG_value)==0 and len(Spo2_value)>0 and OptiFlow==1 and fio2_OF != 'None':
            ratios.append(((spo2/fio2_OF)*1.19)-76.19)
        if len(ABG_value)==0 and len(Spo2_value)>0 and Nasal_Can_Venti ==1:
            ratios.append(((spo2/fio2_NC) * 1.19) - 76.19)

    if mech_invas_status==1:  # invasive or non-invasive ventilation
        if len(ABG_value)>0 and len(Fio2_value)>0:
            ratio_ABG = multi_select_ratio(ABG_value,ABG_time,Fio2_value,Fio2_time)

            ratios.append(ratio_ABG)

        elif len(Spo2_value)>0 and len(Fio2_value)>0:
            ratio_Spo2 = multi_select_ratio(Spo2_value, Spo2_time, Fio2_time, Fio2_value)
            ratios.append((ratio_Spo2 * 1.19) - 76.19)

    if (mech_invas_status+vent_dev_status)==0: # without mech_invas_status,vent_dev_status,
        if len(ABG_value)>0 and len(Fio2_value)>0:
            pao2_final = min(np.array(ABG_value))
            min_pao2_time_final = ABG_time[np.argmin(ABG_value)]
            fio2_final = select_near_fio2(Fio2_time, Fio2_value, min_pao2_time_final)

            ratios.append(pao2_final/fio2_final)

        elif len(Spo2_value)>0 and len(Fio2_value)>0:
            spo2_final = min(np.array(Spo2_value))
            min_spo2_time_final = Spo2_time[np.argmin(Spo2_value)]
            fio2_final = select_near_fio2(Fio2_time, Fio2_value, min_spo2_time_final)

            ratios.append(spo2_final / fio2_final)

    if len(ratios)>0:
        print('--------------------------------------------------')
        min_ratio = min(np.array(ratios))
        print('ratios', ratios, min_ratio)

        if min_ratio >=400:
            resp_score = 0
        if min_ratio <400:
            resp_score = 1
        if min_ratio <300:
            resp_score = 2
        if min_ratio <200 and mech_invas_status==1:
            resp_score = 3
        if min_ratio<100 and mech_invas_status==1:
            resp_score = 4
        if isnan(min_ratio):
            resp_score='None'
    else:
        resp_score = 'None'  # 'None' and None are different in these codes

    return resp_score

def combile_Pao2(combile_data):
    # combile PO2 ARTERIAL {NYP},PO2 (ARTERIAL) - EPOC {NYP}
    Pao2_art = combile_data['PO2 ARTERIAL {NYP}']
    Pao2_Epoc = combile_data['PO2 (ARTERIAL) - EPOC {NYP}']
    for k in Pao2_Epoc.keys():
        if k in Pao2_art.keys():
            art_data = Pao2_art[k]
            Epoc_data = Pao2_Epoc[k]
            A = art_data
            B = Epoc_data

            result = []
            i, j = 0, 0
            while i < len(A) and j < len(B):
                if A[i][1] < B[j][1]:
                    result.append(A[i])
                    i += 1
                else:
                    result.append(B[j])
                    j += 1
            result += A[i:]
            result += B[j:]
            d1 ={k:result}
            Pao2_art.update(d1)

        else:
            Pao2_art[k] = Pao2_Epoc[k] # add new item

    updated_result = {'PO2 ARTERIAL {NYP}':Pao2_art}
    combile_data.update(updated_result)

    return combile_data

def combile_Pco2(combile_data):
    # combile PO2 ARTERIAL {NYP},PO2 (ARTERIAL) - EPOC {NYP}
    Pao2_art = combile_data['PCO2 ARTERIAL {NYP}']
    Pao2_Epoc = combile_data['PCO2 (ARTERIAL) - EPOC {NYP}']
    for k in Pao2_Epoc.keys():
        if k in Pao2_art.keys():
            art_data = Pao2_art[k]
            Epoc_data = Pao2_Epoc[k]
            A = art_data
            B = Epoc_data

            result = []
            i, j = 0, 0
            while i < len(A) and j < len(B):
                if A[i][1] < B[j][1]:
                    result.append(A[i])
                    i += 1
                else:
                    result.append(B[j])
                    j += 1
            result += A[i:]
            result += B[j:]
            d1 ={k:result}
            Pao2_art.update(d1)

        else:
            Pao2_art[k] = Pao2_Epoc[k] # add new item

    updated_result = {'PCO2 ARTERIAL {NYP}':Pao2_art}
    combile_data.update(updated_result)

    return combile_data

def obtain_intubation_flag(pat,start_date,end_date):

    intu_patient_list = data_patient_intu.keys()
    if pat in intu_patient_list:
        start_date_intu = data_patient_intu[pat][0]
        end_date_intu = data_patient_intu[pat][1]
        start_date_intu = time.mktime(time.strptime(start_date_intu, "%Y-%m-%d %H:%M:%S"))  # time to seconds; start_date_o: "2020-04-09 01:40:00",
        end_date_intu = time.mktime(time.strptime(end_date_intu, "%Y-%m-%d %H:%M:%S"))

        if end_date_intu<AdmitDtm or start_date_intu> DischargeDtm:
            intu_flag = 0
        else:
            intu_flag =1
    else:
        intu_flag = 0
    return intu_flag


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
end_point = 60  



with open(FILE_NAME, 'r') as load_f:
    data = json.load(load_f)


if aligned_event == "AdmitDtm":
    data_patient = get_admission_info(event_file=EVENT_FILE)
elif aligned_event == "intubation_start":
    data_patient_intu = get_intubation_info(event_file=EVENT_FILE)
    data_patient = data_patient_intu


data = combile_Pao2(data)
data = combile_Pco2(data)
pat_list = data_patient.keys()
pat_resp_score = {}


for pat in pat_list:
    tem_score = []



    AdmitDtm = data_patient[pat][0]
    DischargeDtm = data_patient[pat][1]
    print('AdmitDtm and DischargeDtm', AdmitDtm, DischargeDtm)

    AdmitDtm = time.mktime(time.strptime(AdmitDtm, "%Y-%m-%d %H:%M:%S"))  # time to seconds; start_date_o: "2020-04-09 01:40:00",
    if end_point == None:
        DischargeDtm = time.mktime(time.strptime(DischargeDtm, "%Y-%m-%d %H:%M:%S"))
    else:
        DischargeDtm = AdmitDtm + end_point * 24 * 60 * 60.0


    start_date = AdmitDtm - len_before_event * 24 * 60 * 60  # time before event
    while start_date <= DischargeDtm:
        end_date = start_date + interval * 60 * 60

        score = compute_respiration_score(data, pat, start_date, end_date) # 86400 seconds=1 day
        start_date = end_date
        if score !='None':
            tem_score.append(float(score))
        else:
            tem_score.append(np.nan)
    pat_resp_score[pat] = tem_score

sample = pat_resp_score


with open(VERSION+'/' + VERSION + '_Respiration_' + aligned_event + "_interval_%s.json" % interval, 'w') as fp:
    json.dump(sample, fp, indent=4)


