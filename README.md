# ThermStack
 A thermal simulation and data fit that allows multilayered stacks
useage:
python Thermstack ConfigFile.yaml
where the ConfigFile.yaml has a number of main keywords defines the 
*Material, Heater, Amplifier, Measured and Layers.  Each main keyword defines:
*Material: one or more materials than can be used to defin the Heater material and the Layers
*Heater: which define the characteristics of the heater layer
*Amplifer: which defines the frequency response files (or a constant) 
*Measured: defines the measured data which the simulation is compared against
*Layers: the stack of layers defined using the Heater and materials defined in the Materials keyword
Within each keyword sub-keywords are defined, namely:
Material:
  mat:
    url: string http website 
    thermal_conductivity: float
    porosity: float (optional)
    thickness: float
    density: float
    temp: float or [list]
    cp: float or [list]
    pressure: float or [list]
    boundary: one of  isothermal ,semi-infinite, adiabatic
optimize: thermal_conductivity and/or thermal_diffusivity and/or  thermal_boundary_resistance 

Heater:
  material: must be in Material list
  R: float in Ohms
  l: float in metres
  width: float in metres
  thickness: float in metres
  Prms: float in Watts
  #thermal_boundary_resistance: float in ...
  capacitance_flag: boolean
 
Amplifier:
  gain_nominal: float (constant)
  use_calibration: boolean
  Calibration_in: string (file in Calibration folder)
  Calibration_out: string (file in Calibration folder)

Measured:
  data: string  (csv file of measurements), can be a separated comma list
  downsample: integer  # use downsamples >1 to select every nth value from the frequecies measured, to reduce computation time if required
  output: string= v =write data to the output file as voltages, otherwise temperatures are written 
  range: float (0-1) # typially 1 but a value <1 for example 0.2 will plot only a sub-range from 0 to fmax*range of the measured data
  ambientTinC: float in Celcius
  model: string of type 'Borca' 'Kim or 'Cahill' 
  plot: booleans
  integration:
    min: integration is min*2/width
    max: integration is max*2/width

Layers:
  - List1
  - List2
  - List3

installation
the code requires
conda install matplotlib
conda install scipy
conda install conda-forge::quadpy

