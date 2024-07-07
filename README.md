# ThermStack
 A thermal simulation and data fit program that allows multi-layered stacks to be analysed using the 3-omega method.

 ## Installation
 the code has been tested on python 3.9 under minconda.  If you don't have python 3.9 installed, I suggest creating a new conda environment during testing using:<br>
 <code>conda create -n thermal python=3.9</code><br>
 and then enter the new environment using <code>conda activate thermal</code>.  In the new <b>thermal</b> environment the following installtions are required:<br>
 <code>conda install matplotlib</code><br>
 <code>conda install scipy</code><br>
 <code>conda install conda-forge::quadpy</code><br>
 <code>conda install -c anaconda pyyaml</code><br>
 <code>conda remove pillow</code>  you need to remove pillow if installed as version 10 caused issues in python 3.9 <br> 
 <code>conda install pillow=9.2.0</code><br>


 ## Useage:<br>
To run the program from the correct conda environment type:<br>
<code> python Thermstack ConfigFile.yaml</code>  (or just <code> python Thermstack </code> for the default .yaml file)<br> 
where the <b>ConfigFile.yaml</b> consists of a number top-level and sub-level keywords followed by defined parameters.  The top-level keywords define the <code>Material, Heater, Amplifier, Measured and Layers</code>.  Each main keyword defines:<br>
* <code>Material</code>: one or more materials than can be used to defin the Heater material and the Layers
* <code>Heater</code>: which define the characteristics of the heater layer
* <code>Amplifer</code>: which defines the frequency response files (or a constant)
* <code>Measured</code>: defines the measured data which the simulation is compared against
* <code>Layers</code>: the stack of layers defined using the Heater and materials defined in the Materials keyword

Under each top-level keyword there sub-keywords are defined where the parameters associated with each sub-keyword are defined below (string, integer, float, boolean or list) and optional arguements are shown with {}. The format follows:<br>
<code>Material</code>:<br>
  &ensp;<code>mat</code>:<br>
    &emsp;<code>url</code>: string http website<br>
    &emsp;<code>thermal_conductivity</code>: float<br>
    &emsp;<code>thermal_diffusivity</code>: float<br>
  &emsp;{<code>porosity</code>: float }<br>
    &emsp;<code>thickness</code>: float<br>
    &emsp;<code>density</code>: float<br>
    &emsp;<code>temp</code>: float or [list]<br>
    &emsp;<code>cp</code>: float or [list]<br>
    &emsp;{<code>pressure</code>: float or [list]}<br>
    &emsp;<code>boundary</code>: one of  isothermal ,semi-infinite, adiabatic<br>
&emsp;{<code>optimize</code>: <b>thermal_conductivity, thermal_diffusivity,  thermal_boundary_resistance</b>}<br>

<code>Heater</code>:<br>
  &emsp;<code>material</code>: must be in Material list<br>
  &emsp;<code>R</code>: float in Ohms<br>
  &emsp;<code>l</code>: float in metres<br>
  &emsp;<code>width</code>: float in metres<br>
  &emsp;<code>thickness</code>: float in metres<br>
  &emsp;<code>Prms</code>: float in Watts<br>
  &emsp;{<code>thermal_boundary_resistance</code>: float in }<br>
  &emsp;<code>capacitance_flag</code>: boolean<br>

<code>Amplifier</code>:<br>
  &emsp;<code>gain_nominal</code>: float (constant)<br>
  &emsp;<code>use_calibration</code>: boolean<br>
  &emsp;<code>Calibration_in</code>: string (file in Calibration folder)<br>
  &emsp;<code>Calibration_out</code>: string (file in Calibration folder)<br>

<code>Measured</code>:<br>
  &emsp;<code>data</code>: string  (csv file of measurements), can be a separated comma list<br>
  &emsp;<code>downsample</code>: integer  # use downsamples >1 to select every nth value from the frequecies measured, to reduce computation time if required<br>
  &emsp;<code>output</code>: string= v =write data to the output file as voltages, otherwise temperatures are written<br>
  &emsp;<code>range</code>: float (0-1) # typially 1 but a value <1 for example 0.2 will plot only a sub-range from 0 to fmax*range of the measured data<br>
  &emsp;<code>ambientTinC</code>: float in Celcius<br>
  &emsp;<code>model</code>: string of type 'Borca' 'Kim or 'Cahill'<br>
  &emsp;<code>plot</code>: boolean<br>
  &emsp;<code>integration</code>:<br>
    &emsp;&emsp;<code>min</code>: integration is min*2/width<br>
    &emsp;&emsp;<code>max</code>: integration is max*2/width<br>

<code>Layers</code>:<br>
  &emsp;<code>- List1</code><br>
  &emsp;<code>- List2</code><br>
  &emsp;<code>- List3</code><br>

--------------------------
Note that the preference is for <code>thermal_diffusivity</code> to be defined for material paremeters.  However, alteratively     <code>density, temp,cp and pressure</code> can be defined from which <code>thermal_diffusivity</code> will be internally calculated.  In this case, keyword <code>pressure</code> is only required if there is a variation if the variation of temperature and pressure is known, see for exampel the definition of Silicon (Si1) in the example <code>thin.yaml</code> file

## Operation
if the optional sub-keyword <code>optimize</code> in included under the top-level <code>Material</code> keyword, then the code attempts to fit the file defined under <code>Measured>data</code> and provide the optimized value of  <b>thermal_conductivity, thermal_diffusivity,  thermal_boundary_resistance</b>, depending on which of these were specified.  Specifying more parameters requires more computation time and may result in an undefined solution.<br>
if the optional sub-keyword <code>optimize</code> in not included:
* the parameters given are simply used to plot the expected 3-omega response
* more than List of layers can be specified under the  <code>Layers</code>  top-level keyword, allowing the impact of different layers and parameters to be observed.  
* more than one .csv data file can be specified under <code>Measured>data</code> if the files are seperated by a comma.  If the <code>optimize</code> keyword in enabled, only the first .csv dataset specified under the <code>Measured>data</code> keywords and first stack sepecified under the <code>Layers</code> is used to find the optimized thermal parameters (from the data fit).


## Examples


## Outstanding documentation
* Requires link to 3-omega papers and others
* Should reference UWA
* Image of screen output under Examples
* Complete Examples text
* Units need to be defined for all parameetrs
