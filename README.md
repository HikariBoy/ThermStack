# ThermStack
 A thermal simulation and data fit program that allows multi-layered stacks to be analysed using the 3-omega method.

 ## Installation
 the code requires<br>
 <code>conda install matplotlib</code>
 <code>conda install scipy</code>
 <code>conda install conda-forge::quadpy</code>

 ## Useage:<br>
To run the program from the correct conda environment type:
<code> python Thermstack ConfigFile.yaml</code><br>
where the <b>ConfigFile.yaml</b> consists of a number top-level and sub-level keywords followed by defined parameters.  The top-level keywords define the <code>Material, Heater, Amplifier, Measured and Layers</code>.  Each main keyword defines:<br>
* <code>Material</code>: one or more materials than can be used to defin the Heater material and the Layers
* <code>Heater</code>: which define the characteristics of the heater layer
* <code>Amplifer</code>: which defines the frequency response files (or a constant)
* <code>Measured</code>: defines the measured data which the simulation is compared against
* <code>Layers</code>: the stack of layers defined using the Heater and materials defined in the Materials keyword

Under each top-level keyword there sub-keywords are defined where the parameters associated with each sub-keyword are defined below (string, integer, float, boolean or list) and optional arguements are shown with {}. The format follows:<br>
<code>Material</code>:
  &ensp;<code>mat</code>:
    &emsp;<code>url</code>: string http website
    &emsp;<code>thermal_conductivity</code>: float
    &emsp;<code>thermal_diffusivity</code>: float
  &emsp;{<code>porosity</code>: float }
    &emsp;<code>thickness</code>: float
    &emsp;<code>density</code>: float
    &emsp;<code>temp</code>: float or [list]
    &emsp;<code>cp</code>: float or [list]
    &emsp;{<code>pressure</code>: float or [list]}
    &emsp;<code>boundary</code>: one of  isothermal ,semi-infinite, adiabatic
&emsp;{<code>optimize</code>: <b>thermal_conductivity, thermal_diffusivity,  thermal_boundary_resistance</b>}

<code>Heater</code>:
  &emsp;<code>material</code>: must be in Material list
  &emsp;<code>R</code>: float in Ohms
  &emsp;<code>l</code>: float in metres
  &emsp;<code>width</code>: float in metres
  &emsp;<code>thickness</code>: float in metres
  &emsp;<code>Prms</code>: float in Watts
  &emsp;{<code>thermal_boundary_resistance</code>: float in }
  &emsp;<code>capacitance_flag</code>: boolean

<code>Amplifier</code>:
  &emsp;<code>gain_nominal</code>: float (constant)
  &emsp;<code>use_calibration</code>: boolean
  &emsp;<code>Calibration_in</code>: string (file in Calibration folder)
  &emsp;<code>Calibration_out</code>: string (file in Calibration folder)

<code>Measured</code>:
  &emsp;<code>data</code>: string  (csv file of measurements), can be a separated comma list
  &emsp;<code>downsample</code>: integer  # use downsamples >1 to select every nth value from the frequecies measured, to reduce computation time if required
  &emsp;<code>output</code>: string= v =write data to the output file as voltages, otherwise temperatures are written
  &emsp;<code>range</code>: float (0-1) # typially 1 but a value <1 for example 0.2 will plot only a sub-range from 0 to fmax*range of the measured data
  &emsp;<code>ambientTinC</code>: float in Celcius
  &emsp;<code>model</code>: string of type 'Borca' 'Kim or 'Cahill'
  &emsp;<code>plot</code>: boolean
  &emsp;<code>integration</code>:
    &emsp;&emsp;<code>min</code>: integration is min*2/width
    &emsp;&emsp;<code>max</code>: integration is max*2/width

<code>Layers</code>:
  &emsp;<code>- List1</code>
  &emsp;<code>- List2</code>
  &emsp;<code>- List3</code>

--------------------------
Note that the preference is for <code>thermal_diffusivity</code> to be defined for material paremeters.  However, alteratively     <code>density, temp,cp and pressure</code> can be defined from which <code>thermal_diffusivity</code> will be internally calculated.  In this case, keyword <code>pressure</code> is only required if there is a variation if the variation of temperature and pressure is known, see for exampel the definition of Silicon (Si1) in the example <code>thin.yaml</code> file

## Operation
if the optional sub-keyword <code>optimize</code> in included under the top-level <code>Material</code> keyword, then the code attempts to fit the file defined under <code>Measured>data</code> and provide the optimized value of  <b>thermal_conductivity, thermal_diffusivity,  thermal_boundary_resistance</b>, depending on which of these were specified.  Specifying more parameters requires more computation time and may result in an undefined solution.
if the optional sub-keyword <code>optimize</code> in not included:
* the parameters given are simply used to plot the expected 3-omega response
* more than List of layers can be specified under the  <code>Layers</code>  top-level keyword, allowing the impact of different layers and parameters to be observed.  
* more than one .csv data file can be specified under <code>Measured>data</code> if the files are seperated by a comma
If the <code>optimize</code> keyword in enabled, only the first .csv dataset specified under the <code>Measured>data</code> keywords and first stack sepecified under the <code>Layers</code>


## Examples


## Outstanding documentation
Requires link to 3-omega papers and others
Should reference UWA
Image of screen output under Examples
Complete Examples text
