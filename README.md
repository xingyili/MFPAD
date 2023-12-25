#### MFPAD

A sequence-based machine learning model for predicting antigenic distance for H3N2 influenza virus

#### DATA

The data includes HA and HI data from 1968 to 2022.

#### Code

feature.py : based on the HA file, perform feature engineering.

MC_main.m: Perform multi-task low-rank matrix completion on the HI data. 

model.py : Train the model based on the established features and antigen distances.
