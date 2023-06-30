import os
import numpy as np
import pandas as pd
import utils


class vanGenuchten:

    theta_r: float  # m3/m3
    theta_s: float  # m3/m3
    alpha: float    # 1/cm
    n: float        # -
    l: float        # -
    K0: float       # cm/d

    def read_parameters_from_table(self, soil_type_short_name: str):
        current_dir = os.path.dirname(os.path.realpath(__file__))
        data_file_path = os.path.join(current_dir, '../data/DIN4220_Soil_Types_vanGenuchten.csv')
        df = pd.read_csv(data_file_path, sep='\t')
        # read table with vG Params
        for _, row in df.iterrows():
            name, soil_type, f_clay, f_silt, f_sand, theta_r, theta_s, alpha, n, l, K0 = row.values
            if soil_type == soil_type_short_name:
                print(f'{name} ({soil_type})')
                self.theta_r = theta_r
                self.theta_s = theta_s
                self.alpha = alpha
                self.n = n
                self.l = l
                self.K0 = K0

                return name, f_sand, f_silt, f_clay

        else:
            raise Exception("Soil type '{}' not found!".format(soil_type_short_name))


    def calculate_parameters_from_equations(self, sand: float, silt: float, clay: float, density: float):
        # we determine the soil type and then read all paramerts from the table
        # we actually only need l and K0, the other ones are calculated from the equations below
        soil_type = utils.read_soil_type_from_table(sand, silt, clay)
        self.read_parameters_from_table(soil_type)
        # the equations use % and g/cm³
        density /= 1000
        sand *= 100
        silt *= 100
        clay *= 100
        if sand < 66.5:
            self.theta_r = 0
            self.theta_s = 0.788 + 0.001*clay - 0.263*density
            self.alpha = np.exp(-0.648 + 0.023*sand + 0.044*clay - 3.168*density)
            self.n = 1.392 - 0.418*sand**(-0.024) + 1.212*clay**(-0.704)
        else:
            self.theta_r = 0
            self.theta_s = 0.89 - 0.001*clay - 0.322*density
            self.alpha = np.exp(-4.197 + 0.013*sand + 0.076*clay - 0.276*density)
            self.n = -2.562 + 7e-9*sand**4.004 + 3.75*clay**(-0.016)

