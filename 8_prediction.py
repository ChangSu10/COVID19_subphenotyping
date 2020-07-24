"""

prediction modeling

"""
seed = 42
random_state= seed

import numpy as np
np.random.seed()
import random
random.seed()

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import json
import os
import sys
import pickle as pkl
import pandas as pd
import scipy
import shap
import math
import numpy as np, scipy.stats as st
from xgboost import plot_importance
from sklearn.ensemble import RandomForestClassifier

from sklearn import preprocessing
from sklearn.preprocessing import OneHotEncoder
from sklearn.compose import ColumnTransformer
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.metrics import auc
import matplotlib.pyplot as plt
from sklearn.model_selection import StratifiedKFold
from scipy import interp

from sklearn.model_selection import KFold
from sklearn.tree import DecisionTreeClassifier
import datetime

import pickle




VERSION = "COVID-0501"



def mean_confidence_interval(data, confidence=0.95):
    # https://www.mathsisfun.com/data/confidence-interval-calculator.html
    if confidence==0.95:
        z = 1.960

    a = 1.0 * np.array(data)
    n = len(a)
    a_mean = np.mean(a)
    a_std = np.std(a)

    CI_h = z * a_std / math.sqrt(n)
    CI_low = a_mean - CI_h
    CI_high = a_mean + CI_h



    return CI_low, CI_high,CI_h


def normalize(all_conti_x_fitted, all_categ_x):

    conti_columns = all_conti_x_fitted.columns.values
    # scaler = preprocessing.StandardScaler()
    scaler = preprocessing.MinMaxScaler()
    all_conti_x_fitted = scaler.fit_transform(all_conti_x_fitted)
    all_conti_x_fitted = pd.DataFrame(all_conti_x_fitted,columns=conti_columns)
    all_conti_x_fitted.to_csv(file_path + 'xx_conti_0.csv')

    onehot_enc = OneHotEncoder(sparse=False)  # dense format
    all_categ_x_fitted = onehot_enc.fit_transform(all_categ_x)
    columns = list(onehot_enc.categories_)
    column_list = []
    for term in columns:
        term = term.tolist()
        column_list+=term
    all_categ_x_fitted=pd.DataFrame(all_categ_x_fitted,columns=column_list)

    all_categ_x_fitted.to_csv(file_path+'xx_0.csv')

    return all_conti_x_fitted, all_categ_x_fitted

def normalize_rf(all_conti_x_fitted, all_categ_x):

    all_conti_x_fitted = all_conti_x_fitted

    onehot_enc = OneHotEncoder(sparse=False)  # dense format
    all_categ_x_fitted = onehot_enc.fit_transform(all_categ_x)
    columns = list(onehot_enc.categories_)
    column_list = []
    for term in columns:
        term = term.tolist()
        column_list+=term

    all_categ_x_fitted=pd.DataFrame(all_categ_x_fitted,columns=column_list)

    all_categ_x_fitted.to_csv(file_path+'xx_cat_rf_0.csv')

    return all_conti_x_fitted, all_categ_x_fitted


def load_ethnicity(file_path):
    print("Processing demographics...\n")

    # load patient.csv table
    df = pd.read_csv(file_path+'COVID-0512_demographics.csv', sep=',', header=0) #, dtype={'clientguid':str}
    merged_df = df

    merged_df['Ethnicity'] = np.nan  # add a new column, 'Race' # the first letter of "Race", not "race"


    for idx, row in merged_df.iterrows():
        # race
        race = row['ethnicity']

        race_split = race
        try:
            if 'ASIAN / PACIFIC ISLANDER' in race_split:
                final_race = 'ASIAN / PACIFIC ISLANDER'
            elif 'AFRICAN AMERICAN' in race_split:
                final_race = 'AFRICAN AMERICAN'
            elif 'CAUCASIAN' in race_split:
                final_race = 'CAUCASIAN'
            elif 'MULTI-RACIAL' in race_split:
                final_race = 'MULTI-RACIAL'
            else:
                final_race = 'OTHER'
        except:
            final_race = 'OTHER'

        merged_df.loc[idx, 'Ethnicity'] = final_race

    merged_df = merged_df[['clientguid', 'Age', 'BMI', 'Sex_female','Ethnicity']]

    merged_df.to_csv(file_path+'XX_before_BMI.csv')

    merged_df['BMI'] = df['BMI'].fillna(value=df['BMI'].median()) # fillna   BMI


    return merged_df

def generate_data_df(group,lab_ven_data,demogrphic_data,sofa_data,day_flag,lab_vent_feature,sofa_feature,demographic_feature):
    # group is df including clientguid,grouplabel(0/1)
    if day_flag=='1st_day':
        lab_value_loc = 4
        sofa_value_loc = 0
    if day_flag=='2nd_day':
        lab_value_loc = 5
        sofa_value_loc =1
    if day_flag=='3rd_day':
        lab_value_loc = 6
        sofa_value_loc = 2
    if day_flag=='4th_day':
        lab_value_loc = 7
        sofa_value_loc =3
    if day_flag=='5th_day':
        lab_value_loc = 8
        sofa_value_loc =4

    tem_lab_ven_feature = lab_vent_feature
    tem_sofa_feature = sofa_feature
    tem_demographic_feature = demographic_feature

    obj_column_name = ['clientguid','group']+tem_lab_ven_feature + tem_sofa_feature + tem_demographic_feature

    tem_df = pd.DataFrame(columns=obj_column_name)

    for i in range(len(group)):
        p_id = group.loc[i,'clientguid']
        # print('')
        p_label = group.loc[i,'group']
# extract lab, ven values
        tem_lab_ven_df = lab_ven_data[str(p_id)]
        # add PF_ratio column
        tem_lab_ven_df['PF_ratio'] = tem_lab_ven_df['Pao2']/(0.01*tem_lab_ven_df['Fio2'])

        tem_lab_ven_values = list(tem_lab_ven_df.loc[lab_value_loc,tem_lab_ven_feature])
# extract sofa values
        tem_sofa_df = sofa_data[str(p_id)]
        # print('sofa',tem_sofa_df)
        tem_sofa_df_values = list(tem_sofa_df.loc[sofa_value_loc,tem_sofa_feature])
# extract demographics values

        tem_dem_df = demogrphic_data.loc[demogrphic_data['clientguid']==p_id,tem_demographic_feature]
        tem_dem_values = list(tem_dem_df.iloc[0])

        all_fea_list = [p_id,p_label] + tem_lab_ven_values + tem_sofa_df_values + tem_dem_values

# add records of a patients as a row
        tem_df.loc[len(tem_df)] = all_fea_list
        # print('s')
    df = tem_df
    return df






def main():
    file_path = 'covid/'

    ethnicity_flag = 1
    k_time = 1

    if ethnicity_flag == 1:
        demographic_feature = ['Age','Gender','Ethnicity','BMI'] # demographic_feature = ['Age','Gender', 'BMI','Ethnicity'] other classifier
    else:
        demographic_feature = ['Age','Gender','BMI']


    # load lab, vent,
    with open(file_path + '_patient_matrices_interval_24h_linear_imputed.pkl',
              'rb') as f:  # _patient_matrices_interval_24h.pkl  _patient_matrices_interval_24h_linear_imputed.pkl
        lab_ven_data = pickle.load(f)
    for p in lab_ven_data:
        tem_df = lab_ven_data[p]
        tem_df['PF_ratio']=tem_df['Pao2'].values/tem_df['Fio2'].values
        lab_ven_data[p] = tem_df


    all_lab_feature = ['WBC','Neutrophil_count' ,'Pao2', 'Fio2','PF_ratio', 'Lymphocyte', 'Platelet', 'Hemoglobin', 'Albumin', 'Triglycerides', 'Procalcitonin',
                       'Ferritin', 'Tni: Troponin', 'D-dimer', 'C-Reactive Protein (CRP)', 'LDH', 'Creatinine', 'Sodium',
                       'Potassium', 'Temperature', 'AST', 'ALT', 'Spo2', 'GCS', 'Bilirubin', 'MAP', 'Neutrophil', 'ESR', #'HbA1c',
                       'Glucose', 'GLOBULIN', 'CK', 'Carbon_dioxide', 'Lactic_acid_level', 'Urine', 'PEEP', 'Plateau_pressure',
                       'Minute_ventilation', 'PH', 'PCO2', 'resp_peak_insp_pres', 'Tidal_volume', 'Lymphocyte_count',
                       'Driving_pressure', 'Static_compliance', 'Tidal_PBW_ratio',
                       'Ventilator_ratio']

    # load demographic,
    demogrphic_data = load_ethnicity(file_path)
    demogrphic_data = demogrphic_data.rename(columns={'Sex_female': 'Gender'})

    # load medication
    all_medication_feature = ['tocilizumab', 'hydroxychloroquine', 'predniSONE', 'methylPREDNISolone',
                              'dexamethasone', 'hydrocortisone', 'cefTRIAXone', 'azithromycin', 'piperacillin-tazobactam',
                              'meropenem', 'vancomycin', 'doxycycline', 'enoxaparin_category', 'heparin_category'] # 'enoxaparin_category', 'heparin_category'

    all_medication_feature_binary = ['tocilizumab', 'hydroxychloroquine', 'predniSONE', 'methylPREDNISolone',
                              'dexamethasone', 'hydrocortisone', 'cefTRIAXone', 'azithromycin', 'piperacillin-tazobactam',
                              'meropenem', 'vancomycin', 'doxycycline'] #,'enoxaparin', 'heparin'

    # load comorbidity
    comorbidity = pd.read_csv(file_path+'all_patients_comorbidity_without_index.csv',index_col=False)

    all_comorbidity_feature = ['cad', 'hf',  'dm', 'htn', 'pulm',
                                'active_cancer',  'ibd',  'immunosuppressed'] 


    # load sofa
    with open(file_path + '_all_sofa.pkl', 'rb') as f:
        sofa_data = pickle.load(f)
    all_sofa_feature = ['Coagulation', 'Liver', 'Central_nervous_system', 'Renal', 'Cardiovascular', 'Respiration']


    path = file_path+'classification_results-significant_lab_features_and_others_with_Neutrophil_count.txt'


    time_stamp = datetime.datetime.now()
    k_time= time_stamp.strftime('%Y.%m.%d-%H:%M:%S')

    print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
    print('this is the experiment time:',k_time)
    print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')

    print('demographic:',demographic_feature)
    print('all_comorbidity_feature:',all_comorbidity_feature)

    for group_flag in [0, 1, 2]: # testing 0 ,1 and 2

        auc_df = pd.DataFrame(columns=['time_points', 'mean', 'ci_lower', 'ci_upper', 'ci_h', 'std', 'sem'])

        group_flag = group_flag
        # load label
        group = pd.read_csv(os.path.join(file_path,'univariate_Sofa_total_score_decentring_True_method_average_baselinegroup_' + str(group_flag) + '.csv'), index_col=False)
        
        # change label: 0->1; 1->0, worsening -> recovering
        group['group'] = group['group'].replace({0: 2})
        group['group'] = group['group'].replace({1: 0})
        group['group'] = group['group'].replace({2: 1})

        for day_flag in ['1st_day','2nd_day','3rd_day','4th_day','5th_day']: #'first_third_day'

            print('-----------group_%d in %s----------' % (group_flag, day_flag))

            lab_vent_feature = all_lab_feature
            sofa_feature = all_sofa_feature

            data_df = generate_data_df(group,lab_ven_data,demogrphic_data,sofa_data,day_flag,lab_vent_feature,sofa_feature,demographic_feature)

            # merge medication
            medication = pd.read_csv(file_path + 'COVID-0512_window_-3_%s.csv'%day_flag, index_col=False)
            data_df = pd.merge(data_df, medication, how='left', on='clientguid')
            # merge comorbidity
            data_df = pd.merge(data_df, comorbidity, how='left', on='clientguid')

            data_df.to_csv(file_path+str(group_flag)+'_'+str(day_flag)+'.csv')

            # data_df = data_df.fillna(data_df.median())
            data_df = data_df.fillna('No')
            data_df['enoxaparin_category'] = data_df['enoxaparin_category'].replace({'No': 'enoxaparin_category_No'})
            data_df['heparin_category'] = data_df['heparin_category'].replace({'No': 'heparin_category_No'})

            data_df.to_csv(file_path + str(group_flag) + '_' + str(day_flag) + 'yy.csv')

            feature_cols = list(data_df.drop(['clientguid', 'group'], axis=1).columns)

            # ------- training  --------
            all_conti_x = data_df[all_lab_feature + all_sofa_feature + ['Age','BMI']] # ,'BMI'

            if ethnicity_flag == 1:
                all_categ_x = data_df[['enoxaparin_category','heparin_category','Ethnicity']] # 'Ethnicity',
                # all_categ_x = data_df[['Ethnicity']]  # 'Ethnicity',
            else:
                all_categ_x = data_df[['enoxaparin_category', 'heparin_category']]  # 'Ethnicity',

            all_binary_x = data_df[['Gender'] + all_medication_feature_binary + all_comorbidity_feature]
            all_binary_x = all_binary_x.replace({'Yes':1,'No':0})

            # print(all_conti_x)
            all_conti_x.to_csv(file_path + 'xxxx_conti.csv', index=False, header=True)
            all_categ_x.to_csv(file_path + 'xxxx_categ.csv', index=False, header=True)
            all_binary_x.to_csv(file_path + 'xxxx_binary_1.csv', index=False, header=True)


            all_conti_x_fitted, all_categ_x_fitted = normalize_rf(all_conti_x, all_categ_x)

            combine_df = pd.concat([all_conti_x_fitted, all_categ_x_fitted,all_binary_x],axis=1)
            combine_df.to_csv(file_path+'group_flag_'+str(group_flag)+'_day_flag_'+day_flag+'combine_df.csv')



            combine_df.to_csv(file_path + 'group_flag_' + str(group_flag) + '_day_flag_' + day_flag + 'combine_df_drop.csv')

            X = combine_df
            feature_cols = X.columns.values.tolist()

            X = X.values
            y = data_df['group']

            params = {'n_estimators': [20, 40, 60, 80, 100, 120, 140, 160, 180, 200],
                      'max_depth': [2, 4, 6, 8, 10, 12, 14, 16, 18, 20],
                      'criterion': ['entropy', 'gini']}
            grid = GridSearchCV(classifier, param_grid=params, scoring='roc_auc', cv=5, n_jobs=4)
            grid.fit(X, y)
            scores = pd.DataFrame(grid.cv_results_)

            grid_aucs = scores['mean_test_score'].values.tolist()
            best_auc = grid.best_score_
            best_auc_std = scores.loc[grid.best_index_, 'std_test_score']

            classifier = grid.best_estimator_

            print('Best by searching: %s, Std: %s' % (best_auc, best_auc_std))

            sem = best_auc_std/math.sqrt(5)
            z = 1.960 # (confidence 0.95)
            CI_h = z*sem
            CI_low = best_auc - CI_h
            CI_high = best_auc + CI_h

            result = {'time_points':day_flag,
                      'mean':best_auc,
                      'ci_lower':CI_low,
                      'ci_upper':CI_high,
                      'ci_h':CI_h,
                      'std': best_auc_std,
                      'sem':sem}
        auc_df = auc_df.append(result,ignore_index=True,sort=False)

main()

