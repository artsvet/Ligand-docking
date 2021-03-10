import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import re
import os
from statistics import mean, stdev

targets = pd.read_csv(str(os.getcwd() + '\\targetsTableNew.csv')).fillna('')
targets['id'].apply(lambda x: x.lower())
sbir_reports = os.walk(
    'C:\\Users\\starch\\PycharmProjects\\molDocking\\finished docks\\').__next__()
csv_final = pd.concat([pd.read_csv(sbir_reports[0]+'\\'+file)
             for file in sbir_reports[2]
             if file[-10:] == '_final.csv'])
csv_out = pd.concat([pd.read_csv(sbir_reports[0]+'\\'+file)
             for file in sbir_reports[2]
             if file[-8:] == '_out.csv'])
