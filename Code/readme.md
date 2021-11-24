## Description

### Simulink models

*  All the Simulink models include the actuator model, sensors, compensation method, DDMRT, and necessary communication blocks.
	* `OpenSeesPy_Raspberry.slx`: Simulink model used for external simulations. In this case it worked with a Raspberry Pi 4 Model B.
	* `OpenSeesPy_Simulink.slx`: Simulink model used for local simulations.
	* `Substructuring.slx`: Simulink model used for substructuring. It is used in local mode, but extendable to external mode.

### Matlab Files
* `initializeSimulation.m`: This function is called by the Simulink model and from this file it loads the necessary variables for correct execution. In this, the variables of the actuator, compensation method, sensors, integration times, among others, are entered.
* `Freq_Resp_Tong.m`: This function is inside the 'initializeSimulation' script. This function based on an input signal, output signal and the sampling frequency, gives the amplitude error between the input/output, phase error, equivalent frequency and delay between the signals.
* `WienerFilter.m`: This function performs the calculation of the W parameters from an input signal.

### Numerical Model (OpenSeesPy scripts)
* `3Stories2Bays_NumericalModel`: 
* `elcentro_0.02.txt`: Text file containing El Centro's 1940 seismic record. This file is read by the numeric substructure model in OpenSees (SubEstNum). In case another seismic record wants to be used, you just have to create the text file and modify the name in the model in OpenSees by the name of the new file. You must know the time spacing of the record and the amount of data to perform enough integration steps to include the entire record.

### Post Processing
* `Output.m`: This script takes the reference structure .mat files, and simulink results and traces them. Within this script, the Freq_Resp_Tong function is called to calculate the delay between input / output signals. In addition, performance indicators J2 and J4 are calculated.
