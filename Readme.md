# Energy pump simulation

Simulation of energy pump. To launch simulation, first set sweep list with SetSweep.py to generate list of parameters to simulate with. Then run LaunchSweep $Nprocesses to start simulation over the parameter sets specified. Runs with $Nprocesses in parallel. 
If running on server, first use SetSweep, and then simply run ./StartSweep and follow the instructions to set the number of parallel processes, and assign server. 

##### Dependencies for simulation

* `PumpFunctions`	 		    : 	Core module (see below)
* `BasicFunctions`			    : 	Helpful functions for logging and handling data
* `TightBindingModel`		    : 	Package for generating tight-binding models
* `DataAnalysisFunctions`		: 	Helpful functions for data analysis

##### Scripts for parallelizing
* `SetSweep` 					: 	Generate file `SweepLists/SweepList.npz`, which contains the parameter sets to run parallelized sweep with 
* `LaunchSweep $Nprocesses`	:  	Sweeps over `SweepLists/SweepList.npz` (by calling `PumpFunctions.EnergyAbsorptionSweep`) by parallelizing to `$Nprocess processes`. Data saves to `../Data` folder.  

##### Scripts for data handling
* `DataCleanUp`				: 	Reads through raw data and sorts it to easily-readable data in `../DataClean` Folder
* `DataAnalysisPcolor`		: 	Plots energy absorption and correlation length vs. disorder strength and frequency
* `DataAnalysisFiniteSize`	: 	1d plots of finite-size scaling 

##### Scripts for running on server
* `UploadScript`			: 	Upload current state of code (and current sweeplist) to server
* `DownloadScript`			: 	Download current state of raw data on server to ../Data folder
* `StartSweep`				: 	Launch parallel sweep on a server, with current sweeplist. 


## PumpFunctions
Core module, where the simulation is done. Contains two functions: `FindEnergyAbsorption`, and `EnergyAbsorptionSweep`. 

### FindEnergyAbsorption
`
[E,B,C,DisorderVec] = FindEnergyAbsorption([L,Nperiods,PBC,omega2,W],OutputProgress=True)
`

The core function. Constructs the model and simulates the model over `Nperiods` periods. Extracts energy absorption and correlation length.  
  
##### Input:
* `L`             : number of unit cells
* `Nperiods`      : number of periods of mode 1 (step cycles) to run over
* `PBC`           : periodic boundary conditions. If 1, periodic boundary condtitions, if 0, open boundary conditions
* `omega2`        : angular frequency of ramping cycle
* `W`             : Disorder strength
* `OutputProgress`: If 1, display progress on console, if 0, do not display progress on console
   
##### Returns: 
* `E`            : average energy absorption over the cycle, when initially filling the left half of the system (`L` leftmost sites)
* `B`            : standard deviation of energy absorption over all sites (measure of bulk absoption)
* `[Cmean,Cmax]` : mean and max correlation length of U(Nperiods *T), where U is the time-evolution operator (i.e. the correlation length of the collumns of U)
* `DisorderVec`  : Disorder potential used (can be used to reproduce results)


### EnergyAbsorptionSweep: 

` 
EnergyAbsorptionSweep(ParameterList,FileString)
`

Core sweep function: Evaluates `FindEnergyAbsoprtion` over a list of parameter-sets, and saves all data into a single file.

##### Input:
* `ParameterList`	: List of parameter sets to be run over. Each element in `ParameterList` (parameter set) a list of the format `[L,Nperiods,PBC,omega2,W]`.
* `FileString` 		: String to append to datafile name. 

Executes EnergyAbsorptionSweep for each parameter set in ParameterList serially (through a foor loop). 
* Saves arrays `
Elist, Blist, Clist, ParameterList
` (see below) to data file `Data/FileString_$Run_ID.npz`, where `$Run_ID=YYMMDD_HHMM-SS.MS` gives the time of the execution of `EnergyAbsorptionSweep` in millisecond precsion. Here `ParameterList` was given as input, and 

`[Elist[n],Blist[n],Clist[n]]=FindEnergyAbsorption(ParameterList[n])`.



 

 