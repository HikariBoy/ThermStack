# ThermStack
<img src=".\UWA-Full-Ver-CMYK3.png" alt="UWA logo"  align="right" width="150"/><br>
 A thermal simulation and data fit program that allows multi-layered stacks to be analysed using the 3-omega method.  This code should be used in conjuntion with the 3-omega electrical circuit detailed by [Dr. Sobhan Erfantalab](https://github.com/Sobhan10100101/3omega-method-signal-conditioning-circuit-PCB ) which provides a wide frequency band, low third-harmonic noise circuit which enables high quality 3-omega measurements of thin films and bulk substrates.  

 ## Installation
 the code has been tested on python 3.9 under minconda.  If you don't have python 3.9 installed, I suggest creating a new conda environment during testing using:<br>
 <code>conda create -n thermal python=3.9</code><br>
 and then enter the new environment using <code>conda activate thermal</code>.  In the new <b>thermal</b> environment the following installations are required:<br>
 <code>conda install matplotlib</code><br>
 <code>conda install scipy</code><br>
 <code>conda install conda-forge::quadpy</code><br>
 <code>conda install -c anaconda pyyaml</code><br>
 <code>conda remove pillow</code>  you need to remove pillow if installed as version 10 caused issues in python 3.9 <br>
 <code>conda install pillow=9.2.0</code><br>


 ## Useage:<br>
To run the program from the correct conda environment (<b>thermal</b> if following the above example) type:<br>
<code> python ThermStack ConfigFile.yaml</code>  (or just <code> python ThermStack </code> for the default .yaml file)<br>
where the <b>ConfigFile.yaml</b> consists of a number top-level and sub-level keywords followed by defined parameters.  The top-level keywords define the <code>Material, Heater, Amplifier, Measured and Layers</code>.  Each main keyword defines:<br>
* <code>Material</code>: one or more materials than can be used to define the Heater material and the Layers
* <code>Heater</code>: which define the characteristics of the heater layer used in the 3-omega method
* <code>Amplifer</code>: which defines the frequency response files (or a constant)
* <code>Measured</code>: defines the measured data which the simulation is compared against
* <code>Layers</code>: the stack of layers defined using the Heater and materials defined in the Materials keyword

Under each top-level keyword sub-keywords are defined where the parameters associated with each sub-keyword are defined below (string, integer, float, boolean or list) and optional arguements are shown with {}. The format follows:<br>
<code>Material</code>:<br>
  &ensp;<code>mat</code>:<br>
    &emsp;{}<code>url</code>: string http website<br>}
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
  &emsp;<code>material</code>: must be a <code>mat</code> defined under the <code>Material</code> keyword<br>
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
  &emsp;<code>downsample</code>: integer  # use downsamples >1 to select every nth value from the frequencies measured, to reduce computation time if required<br>
  &emsp;<code>output</code>: string= v =write data to the output file as voltages, otherwise temperatures are written<br>
  &emsp;<code>range</code>: float (0-1) # typially 1 but a value <1 will use only a sub-range of frequencies, for example 0.2 will use only a sub-range from 0 to (fmax*range) of the measured data<br>
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
* more than one List of layers can be specified under the  <code>Layers</code>  top-level keyword, allowing the impact of different layers and parameters to be observed.  
* more than one .csv data file can be specified under <code>Measured>data</code> if the files are separated by a comma.  If the <code>optimize</code> keyword in enabled, only the first .csv dataset specified under the <code>Measured>data</code> keywords and first stack specified under the <code>Layers</code> is used to find the optimized thermal parameters (from the data fit).


## Examples
### Bulk Glass
Running the code:<br>

<code> python ThermStack bulk.yaml</code>  <br>
loads measured 3-omega data for a [bulk glass substrate](Glass.csv) obtained by [Dr. Sobhan Erfantalab as part of his Ph.D.](https://doi.org/10.26182/qtxb-2f91) and attempts to fit the model proposed by [D. Cahil - Review of Scientific Instruments 61, 802 (1990)](https://doi.org/10.1063/1.1141498) to the measured data.  <img src=".\ThermStackGlass.png" alt="Bulk Material Fit"  align="right" width="250"/>The resulting measured data and fit are shown along with an illustration of the various layers in the model being fit.  From the data fit shown the thermal parameters extracted from the [bulk glass substrate](Glass.csv) indicated a thermal conductivity of 1.23 W/mK+/-1.27% and thermal diffusivity of 0.43 $mm^2/s$+/-4.64%  for this glass. The original approach to extract thermal parameters by [D. Cahil](https://doi.org/10.1063/1.1141498) was to assume a linear variation of temperature over the low frequency region and use the linearity over this region to allow the thermal conductivity to be directly extracted.  Inspection of the data shown right indicate why this works so for bulk materials at low frequencies.  However the real and imaginary spectra over the entire frequency range is rich with information and [our findings previously reported](https://doi.org/10.1016/j.applthermaleng.2022.119965) showed that fitting the complex data over the entire frequency range allows both thermal conductivity and thermal diffusivity to be extracted simultaneously. <br><br>
Running the code:<br>

<code> python ThermStack thin.yaml</code>  <br>
<img src=".\ThermStackThin.png" alt="ThinFilmFit"  align="right" width="250"/><br>
loads measured 3-omega data for a [thin film of porous silicon](PS77_data.csv), nominally 3-microns thick, which was fabricated by [Dr. Sobhan Erfantalab as part of his Ph.D.](https://doi.org/10.26182/qtxb-2f91) and attempts to fit the thin film stack model proposed by  [T. Borca-Tasciuc in  Review of Scientific Instruments 72, 2139 (2001)](https://doi.org/10.1063/1.1353189) to the data. The resulting measured data and fit (stack#1) is shown along with a 2nd model (stack#2) for a different set of thermal parameters, allowing the effect of different values for thermal conductivity and thermal diffusivity to be compared with optimized values (k(opt) and d(opt) values indicated for stack#1).  From the data fit the thermal parameters extracted from the [thin film of porous silicon](PS77_data.csv)  indicated a thermal conductivity of 0.81 W/mK+/-1.01% and thermal diffusivity of 1.23 $mm^2/s$+/-8.06%.  Observing the complex thermal spectra shown right, it can be seen the looking at only the low frequency region typical of [Cahil's3-omega analysis](https://doi.org/10.1063/1.1141498) would ignore important information throughout the spectra.  This rich spectral data allows [both thermal conductivity and thermal diffusivity to be extracted simultaneously](https://doi.org/10.1016/j.applthermaleng.2022.119965) for thin films<br>

## Outstanding documentation
* Complete Examples text
* Units need to be defined for all parameetrs
