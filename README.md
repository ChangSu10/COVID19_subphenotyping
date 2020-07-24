# COVID19_subphenotyping
Implementation of article "Identifying organ dysfunction trajectory-based subphenotypes in critically ill patients with COVID-19"


1_data_extraction.py ---- Extract data from raw data files.

2_extract_patient_event.py ---- Extract event time of patient, including admission, intubation, extubation, death, etc.

3_construct_sequence_data.py ---- Construct time sequence-like data for each feature of each patient. All sequences are alined at the time of intubation. For each feature, we pick a value within every 24-hour bucket, by selecting the worst value within this bucket.

4_respiration_score.py ---- Derive SOFA respiration score

5_cardiovascular_score_new.py ---- Derive SOFA cardiovascular score

6_generate_patient_matrix.py ---- Generate feature matrix for each patient, time X feature

7_subphenotyping_model.py ---- Subphenotyping model

8_prediction.py ---- Prediction model
