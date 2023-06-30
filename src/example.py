from srcSoilGenerator import SoilGenerator
from vanGenuchten import vanGenuchten


# Parameter in SI units
density = 1500 # kg/m³
f_sand = 0.3 # m3/m3
f_silt = 0.5 # m3/m3
f_clay = 0.2 # m3/m3

m6_file_path = 'test_soil.m6'


"""
Option 1: Use custom van-Genuchten paramters
Caution: van-Genuchten parameters always in common non-SI units (cm/d, 1/cm, ...)!
"""
vg = vanGenuchten()
vg.theta_r = 0.0 # m3/m3
vg.theta_s = 0.4214 # m3/m3
vg.alpha = 0.18023 # 1/cm
vg.n = 1.1323
vg.l = -3.42
vg.K0 = 305.8 # cm/d
# All properties should be SI units: no '%' but '-'
sg = SoilGenerator(f_sand, f_silt, f_clay, density, vg)
sg.create_m6_file(m6_file_path)


"""
Option 2: Read van-Genuchten Parameters from Table of DIN 4220 by specifying soil type name ('Sl2', 'Lt3', ...)
Density is constant 1500 kg/m3 for moisture related functions (MRC, Kl)
For thermal conductivity, the given density is considered 
All properties should be SI units: no '%' but '-'
"""
vg = vanGenuchten()
vg.read_parameters_from_table('Sl2')
sg = SoilGenerator(f_sand, f_silt, f_clay, density, vg)
sg.create_m6_file(m6_file_path)


"""
Option 3: van-Genuchten Parameter are determined from fitting equations 
Density is now considered for moisture related functions (MRC, Kl)
All properties should be SI units: no '%' but '-'
"""
vg = vanGenuchten()
vg.calculate_parameters_from_equations(f_sand, f_silt, f_clay, density)
sg = SoilGenerator(f_sand, f_silt, f_clay, density, vg)
sg.create_m6_file(m6_file_path)

