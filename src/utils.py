import os
import pandas as pd


def read_soil_type_from_table(sand: float, silt: float, clay: float):
    current_dir = os.path.dirname(os.path.realpath(__file__))
    data_file_path = os.path.join(current_dir, '../data/DIN4220_Soil_Types_vanGenuchten.csv')
    df = pd.read_csv(data_file_path, sep='\t')
    min_diff = 999
    determined_soil_type = ""
    determined_soil_name = ""
    for _, row in df.iterrows():
        name, soil_type, fClay, fSilt, fSand, theta_r, theta_s, alpha, n, l, K0 = row.values
        diff = abs(fClay / 100 - clay) + abs(fSilt / 100 - silt) + abs(fSand / 100 - sand)
        if diff < min_diff:
            min_diff = diff
            determined_soil_type = soil_type
            determined_soil_name = name

    print("Bestimmte Bodenart: '{0:} ({1:})'".format(determined_soil_name, determined_soil_type))
    return soil_type
