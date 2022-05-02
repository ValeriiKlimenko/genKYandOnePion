# genKYandOnePion

detailed description (CLAS12-Note 2021-003) is available at 
https://misportal.jlab.org/mis/physics/clas12/index.cfm?note_year=2021

This version has wider Q2 coverage up to Q2 = 30 GeV2 

EG is based on:
1) at Q2 < 5GeV2 interpolation of the available CLAS data 
2) at Q2 > 5GeV2 extrapolation of the available CLAS data that is based on Operator Product Expansion 

Requirments: Cern root.

Need to define the path for the data files using env variable DataKYandOnePion to define the path to the "data" directory

Using (if you are using executable file skip following steps: 2-4):

1) install root (https://root.cern.ch/building-root) or type command(it is an example): use root/6.10.02
2) git clone of the EG: git clone https://github.com/ValeriiKlimenko/genKYandOnePion/ 
3) go to the directory: cd genKYandOnePion/
4) Then type command: chmod +x COMPILE_gen
5) To compile: ./COMPILE_gen
6) Set enviroment varible(In tcsh): setenv DataKYandOnePion /WAY/TO/THE/DATA/FOLDER

An example: setenv DataKYandOnePion data

6) run exe file: ./genKYandOnePion --channel=(here should be channel name)

An example: ./genKYandOnePion --channel=KLambda

Options:

--channel is channel name "KLambda" or "KSigma" or "Pi0P" or "Pi0P_2g" or "PiN" ("KLambda" is default value)

Pi0P - EG will generate e,p,pi0

Pi0P_2g is Pi0P with pi0 decay, EG will generate e,p,g,g

--lambda1520=yes It works only for "KLambda" channel. It sets Lambda mass = 1518.8 MeV (default masss is 1115.683 MeV)

--ebeam is Ebeam that should be >0 and less than 30GeV (Energy of the beam) (default value is 10.6 GeV)

--q2min is Q2 min in GeV2 is more than 0.01 and less than 30 GeV2 (default value is 0.05 GeV2)

--q2max is Q2 min in GeV2 >0.01 and less then 30 GeV2 (default value is 5 GeV2)
This is the range in Q2, where the events will be generated.

--w_min is W min value in GeV. wmin should be >0,  if Wmin less than treshold, treshold value will be used. (default value is threshold GeV) 

--w_max is W min value in GeV less than 4 GeV. (default value is 4 GeV)

--v_z_min is vertex z min value in cm (default value is 0 cm)

--v_z_max is vertex z max value in cm (default value is 0 cm)

--targDiameter is target diameter in XY plane in cm (default value is 0 cm, it means that x and y will be 0 for all events)

--trig is number of events to generate >0  (default value is 100000)

--seed it can be unsigned value, the default is current time in microseconds

--outname is the name of the output file. (default value is genKYandOnePion.dat)

--docker is an option that sets all the parametes to default values 

An example of all options use:

./genKYandOnePion --channel=KLambda --ebeam=10.6 --q2min=0.5 --q2max=11 --w_min=0.7 --w_max=3.7 --v_z_min=-9 --v_z_max=5 --targDiameter=3 --trig=5000 --outname=result.dat

	
The output will be a lund file with the name that you set in outname option (genKYandOnePion.dat is default).


Contact: valerii@jlab.org
