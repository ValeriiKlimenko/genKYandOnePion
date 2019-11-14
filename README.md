# genKYandOnePion
EG is based on:
1) at Q2 < 5GeV2 interpolation of the available CLAS data 
2) at Q2 > 5GeV2 extrapolation of the available CLAS data that is based on Operator Product Expansion 

Requirments: Cern root.

Using (if you are using executable file skip following steps: 2-4):

1) install root (https://root.cern.ch/building-root) or type command: use root/6.10.02
2) git clone of the EG: git clone https://github.com/ValeriiKlimenko/genKYandOnePion/ 
3) Then type command: chmod +x COMPILE_eg_ky
4) To compile: ./COMPILE_eg_ky
5) setenv DataKYandOnePion /WAY/TO/THE/DATA/FOLDER

Here is an example for my local PC, do not copy this command: 
setenv DataKYandOnePion /home/CLAS12/EG/data

6) to run exe file: ./eg_ky arg1 arg2 arg3 arg4 arg5 arg6 arg7 arg8

To set the parameters use arg1-arg8
 Its format is as follows.
 - arg1 is the string "KLambda" or "KSigma" or "Pi0P" or "PiN"
 - arg2 is Ebeam that should be >0 and less than 12GeV (Energy of the beam)
 - arg3 and arg4 contain two numbers Q2min and Q2max.
   This is the range in Q2, where the events will be generated.
 - arg5 and arg6 are two numbers Wmin and Wmax.
   This is the range in W, where the events will be generated.
   Wmin more then zero (GeV); (if Wmin less than treshold, treshold value will be used)
   Wmax less then 4 GeV;
 - arg7 is the number of events to be generated.
 - arg8 is the name of the output file.
 
 An example: ./eg_ky KSigma 11. 2. 11.999 1.5 4.0 5000 lund_KS.lund
	
The output will be a lund file with the name that you set in outputFileName variable.

## Online submitions to OSG

For online submissions to the OSG the output file name and number of events are controlled by the online submission form. Because of this the online form should include:

```sh
type beam_energy Q2min Q2max Wmin Wmax
```

For example:

```sh
KSigma 11. 2. 11.999 1.5 4.0
```


Contact: valerii@jlab.org
