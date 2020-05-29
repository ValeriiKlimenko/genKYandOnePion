# genKYandOnePion
EG is based on:
1) at Q2 < 5GeV2 interpolation of the available CLAS data 
2) at Q2 > 5GeV2 extrapolation of the available CLAS data that is based on Operator Product Expansion 

Requirments: Cern root.

Need to define the path for the data files using env variable DataKYandOnePion to define the path to the "data" directory

Using (if you are using executable file skip following steps: 2-4):

1) install root (https://root.cern.ch/building-root) or type command: use root/6.10.02
2) git clone of the EG: git clone https://github.com/ValeriiKlimenko/genKYandOnePion/ 
3) Then type command: chmod +x COMPILE_gen
4) To compile: ./COMPILE_gen
5) setenv DataKYandOnePion /WAY/TO/THE/DATA/FOLDER

6) run exe file: ./genKYandOnePion --channel=(here should be channel name)

Options:

--channel is chanel name "KLambda" or "KSigma" or "Pi0P" or "PiN" (No default value!)

--ebeam is Ebeam that should be >0 and less than 12GeV (Energy of the beam) (default value is 10.6 GeV)

--q2min is Q2 min in GeV2 is more than 0.01 and less than 12 GeV2 (default value is 0.05 GeV2)

--q2max is Q2 min in GeV2 >0.01 and less then 12 GeV2 (default value is 11 GeV2)
This is the range in Q2, where the events will be generated.

--w_min is W min value in GeV. wmin should be >0,  if Wmin less than treshold, treshold value will be used. (default value is threshold GeV) 

--w_max is W min value in GeV less than 4 GeV. (default value is 4 GeV)

--trig is number of events to generate >0  (default value is 100000)

--outname is the name of the output file. (default value is genKYandOnePion.dat)

An example of all options use:

./genKYandOnePion --channel=KLambda --ebeam=10.6 --q2min=0.5 --q2max=11 --w_min=0.7 --w_max=3.7 --trig=5000 --outname=result.dat

	
The output will be a lund file with the name that you set in outname option (genKYandOnePion.dat is default).

## Online submitions to OSG

For online submissions to the OSG the output file name and number of events are controlled by the online submission form. Because of this the online form should include:


?????


Contact: valerii@jlab.org
