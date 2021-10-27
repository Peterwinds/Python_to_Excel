#======================Block 0: Adding Libraries=============================
import sys
import pandas as pd
import numpy as np
from datetime import datetime as dt, timedelta, date
import matplotlib.pyplot as plt
import pymysql
pymysql.install_as_MySQLdb()
from sqlalchemy import create_engine
#import WINDSfunctionsandclasses13.py

#exec(open("WINDSfunctionsandclasses_July25.py").read())
pathprefix='/home/ecoslacker/Documents/WINDS_Data/'
# sys.path.append(pathprefix)

# db=create_engine('mysql://UofABEWINDS:WINDSAWSPort2020@windsdatabase-1.cdzagwevzppe.us-west-1.rds.amazonaws.com:3306/winds_test')

# #======================Block 0: Adding Libraries=============================

# import WINDSfunctionsandclasses  as wmf

# Planting_Array=pd.read_sql('SELECT * from plantings',con=db)  #Reads all data from mysql db
# Field_Array=pd.read_sql('SELECT * from fields',con=db)
# Status_Array=pd.read_sql('SELECT * from Status',con=db)
# Soil_Array=pd.read_sql('SELECT * from field_soil_layers',con=db)
# Irrigation_Array=pd.read_sql('SELECT * from irrigation_activity',con=db)
# ET_Daily_Array=pd.read_sql('SELECT * from et_daily',con=db)
# RS_Daily_Array=pd.read_sql('SELECT * from RS_daily',con=db)
# Output_Array=pd.read_sql('SELECT * from Output',con=db)
# Output_Layer_Array=pd.read_sql('SELECT * from LayerOutput',con=db)
# WettingArray=pd.read_sql('SELECT * from Wetted_fractions', con=db)
# ET_frac_Array = pd.read_sql('SELECT * from ET_fractions', con=db)
# Status_Array = pd.read_sql('SELECT * from Status', con=db)

# db=create_engine('mysql://UofABEWINDS:WINDSAWSPort2020@windsdatabase-1.cdzagwevzppe.us-west-1.rds.amazonaws.com:3306/winds_test')

#======================Block 0: Adding Libraries=============================


excel_path = 'WINDS_database_format3.xlsm' 

Planting_Array=pd.read_excel(excel_path, sheet_name = 'Plantings')  #Reads all data from mysql db
Field_Array=pd.read_excel(excel_path, sheet_name = 'Fields')
Status_Array=pd.read_excel(excel_path, sheet_name = 'Status')
Soil_Array=pd.read_excel(excel_path, sheet_name = 'Field soil layers')
Irrigation_Array=pd.read_excel(excel_path, sheet_name = 'Irrigation')
ET_Daily_Array=pd.read_excel(excel_path, sheet_name = 'ET_daily')
RS_Daily_Array=pd.read_excel(excel_path, sheet_name = 'RS_daily')
RS_Daily_Array=pd.read_excel(excel_path, sheet_name = 'RS_daily')
WettingArray=pd.read_excel(excel_path, sheet_name = 'Wetting_fractions_db')
ET_frac_Array = pd.read_excel(excel_path, sheet_name = 'ET_fractions_db')
Status_Array = pd.read_excel(excel_path, sheet_name = 'Status')

#excel_path = pathprefix + 'WINDS guayule_guarNM.xlsx' 
#excel_output_path = pathprefix + 'output_test3.xlsx' 

import WINDSfunctionsandclasses_Excel2  as wmf
P = wmf.plantings(Planting_Array) #This sets up a planting object called P with all of the plantings and their associated information
Num_plantings = len(P.PlantingDate) #This is the total number of plantings in the database based on the length of the planting fields.
for i in range(0, Num_plantings): #this loop cycles through all of the plantings
    if (P.RunPlanting[i] == 1): #checks to see whether a planting is active and should be analyzed
        Planting_Array_in = Planting_Array.loc[i]
        WeatherArray = pd.read_excel(excel_path, sheet_name = P.WeatherSheetName[i])
        Weather_Array_in = WeatherArray.loc[(WeatherArray['Date'] >= dt.strptime(P.StartDate[i], '%Y-%m-%d')-timedelta(days=1))]
        Status_in = Status_Array.loc[Status_Array['Planting num']==P.PlantingNum[i]]
        Field_Array_in = Field_Array.loc[(Field_Array['Field num']==P.FieldNum[i])]
        Soil_Array_in = Soil_Array.loc[(Soil_Array['Field num']==P.FieldNum[i])]
        Wetting_Array_in = WettingArray.loc[(WettingArray['Wetting_fractions_number']==P.WettingFractionsNum[i])]
        ETfrac_Array_in = ET_frac_Array.loc[(ET_frac_Array['ET_fractions_number']==P.ETFractionsNum[i])]
        Irrigation_Array_in = Irrigation_Array.loc[(Irrigation_Array['Irr_num']==P.IrrNum[i])]
        ET_Daily_Array_in = ET_Daily_Array.loc[(ET_Daily_Array['Planting num']==P.PlantingNum[i])]                                                                                                          
        RS_Daily_Array_in = RS_Daily_Array.loc[(RS_Daily_Array['Planting num']==P.PlantingNum[i])]       
        
        M = wmf.model(Planting_Array_in,
              Weather_Array_in, #reads in the weather values for the season for the station associated with the planting from Weather table
              Field_Array_in,
              Status_in,
              Soil_Array_in, #reads in the soil layer parameters from the layer_data table for the loc ation
              ETfrac_Array_in,
              Wetting_Array_in,
              Irrigation_Array_in,
              ET_Daily_Array_in,
              RS_Daily_Array_in)
        
        M.run_model()


        #plt.figure(2)
        #plt.plot(M.
        
        # plt.figure(2)
        # plt.plot(M.rainfall)
        # plt.plot(M.Rain_Infilt)
        # plt.ylabel('mm')
        # plt.xlabel('DOY')
        # plt.show()

        
       
        # temp_layer_output.to_sql('TempLayerOutput', db, if_exists='append', index = 0)
        # temp_output.to_sql('TempOutput',db,if_exists='append', index = 0)

