# DelphinSoilGenerator
A Python Module for the generation of DELPHIN *.m6 material files for soils.

## Data Structure

The `src` directory contains these classes:

`DelPhynMat.py`: Module for reading / writing DELPHIN m6 files 

`SoilGenerator.py`: Module containing all functions for creating the material properties

`vanGnuchten.py`: Holds the van Genuchten parameters and offers possibility to read them from a table or calculate them from empirical equations 

`examples.py`: Demonstrates the usage of the SoilGenerator with different options


## Dependencies
- pandas (1.3.5)
- numpy (1.21.5)
- matplotlib (3.5.1)

## Usage

You may generate m6 file for the soil type 'Sl2' with density of 1350 kg/m³

````
from SoilGenerator import SoilGenerator
from vanGenuchten import vanGenuchten

density = 1350 # kg/m³
soil_type = 'Sl2'

vg = vanGenuchten()
name, f_sand, f_silt, f_clay = vg.read_parameters_from_table(soil_type)
sg = SoilGenerator(f_sand, f_silt, f_clay, density, vg)
sg.create_m6_file(m6_file_path)

````

Or you want to specify the grain size distribution and density

```
from SoilGenerator import SoilGenerator
from vanGenuchten import vanGenuchten

density = 1600 # kg/m³
f_sand = 0.3 # m3/m3
f_silt = 0.5 # m3/m3
f_clay = 0.2 # m3/m3

vg = vanGenuchten()
vg.calculate_parameters_from_equations(f_sand, f_silt, f_clay, density)
sg = SoilGenerator(f_sand, f_silt, f_clay, density, vg)
sg.create_m6_file(m6_file_path)
```