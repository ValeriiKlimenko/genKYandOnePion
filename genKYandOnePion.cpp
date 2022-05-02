//SYS LIBRARIES
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <vector>

//ROOT LIBRARIES
#include <TLorentzVector.h>

//#include "sigmaKY.h"
#include "sigmaValera.h"
#include "constants.h" 
#include "kinematics.h"
#include "evGenerator.h"

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>

#include <TRandom3.h>
#include <TRandomGen.h>
#include <sys/time.h>


#include <chrono>
#include <ctime> 

using namespace std;


int main(int argc, char *argv[]) {

  int channel;
	string channelName="KLambda", outputFileName="genKYandOnePion.dat",dataPath;
	double Ebeam=24., Q2min=2., Q2max=20., Wmin=1., Wmax=4.,V_z_min=0.,V_z_max=0.;
	double target_diameter = 0.;
	bool isLam1520 = false;
	int nEventMax = 100000;
	
	auto time = std::chrono::system_clock::now();
	std::chrono::microseconds ms =
					std::chrono::duration_cast<std::chrono::microseconds>(time.time_since_epoch()); 
	long long unsigned seed_value = ms.count(); 
	
    double jr, mr, gr, a12, a32, s12, onlyres;
    
  
    char* short_options = (char*)"a:b:c:d:e:f:g:h:i:j:k:l:p:r:s::";
    const struct option long_options[] = {
        {"channel",required_argument,NULL,'a'},
        {"ebeam",required_argument,NULL,'b'},
        {"q2min",required_argument,NULL,'c'},
        {"q2max",required_argument,NULL,'d'},
        {"w_min",required_argument,NULL,'e'},
        {"w_max",required_argument,NULL,'f'},
        {"trig",required_argument,NULL,'g'},
        {"v_z_min",required_argument,NULL,'i'},
        {"v_z_max",required_argument,NULL,'j'},
        {"targDiameter",required_argument,NULL,'k'},
        {"outname",required_argument,NULL,'l'},
        {"docker",optional_argument,NULL,'p'},
        {"seed",optional_argument,NULL,'r'},
        {"lambda1520",optional_argument,NULL,'s'},
        {NULL,0,NULL,0}
    };

		int rez;
		int option_index;
		
		cout<<'\n';

   while ((rez=getopt_long(argc,argv,short_options,
		long_options,&option_index))!=-1){

		switch(rez){
			case 'h': {
				printf("This is demo help. Try -h or --help.\n");
				break;
			};
			case 'a': {
				if (optarg!=NULL){
					cout<<"channel is set to "<<optarg<<endl;
					channelName=optarg;
				}
				else{
					printf("found channel without value\n");
					cout<<"default value (KLambda) will be used"<<endl;
				}
				break;
			};
	
			case 'b': {
				if (optarg!=NULL){
					cout<<"ebeam is set to "<<optarg<<" GeV"<<endl;
					Ebeam=atof(optarg);
				}
				else{
					printf("found ebeam without value\n");
					cout<<"default value will be used"<<endl;
					}
				break;
			};
			
			case 'c': {
				if (optarg!=NULL){
					cout<<"q2min is set to "<<optarg<<" GeV2"<<endl;
					Q2min=atof(optarg);
				}
				else{
					printf("found q2min without value\n");
					cout<<"default value will be used"<<endl;
					}
				break;
			};
			
			case 'd': {
				if (optarg!=NULL){
					cout<<"q2max is set to "<<optarg<<" GeV2"<<endl;
					Q2max=atof(optarg);
				}
				else{
					printf("found q2max without value\n");
					cout<<"default value will be used"<<endl;
					}
				break;
			};
			
			case 'e': {
				if (optarg!=NULL){
					cout<<"w_min is set to "<<optarg<<" GeV"<<endl;
					Wmin=atof(optarg);
				}
				else{
					printf("found w_min without value\n");
					cout<<"default value will be used"<<endl;
					}
				break;
			};
			
			case 'f': {
				if (optarg!=NULL){
					cout<<"w_max is set to "<<optarg<<" GeV"<<endl;
					Wmax=atof(optarg);
				}
				else{
					printf("found w_max without value\n");
					cout<<"default value will be used"<<endl;
					}
				break;
			};
			
			case 'g': {
				if (optarg!=NULL){
					cout<<"trig (number of events) is set to "<<optarg<<" events"<<endl;
					nEventMax=atoi(optarg);
				}
				else{
					printf("trig (number of events) without value\n");
					cout<<"default value will be used"<<endl;
					}
				break;
			};
			
			case 'i': {
				if (optarg!=NULL){
					cout<<"min Vertex-z is set to "<<optarg<<" cm"<<endl;
					V_z_min=atof(optarg);
				}
				else{
					printf("min Vertex-z without value\n");
					cout<<"default value (zero) will be used"<<endl;
					}
				break;
			};
			
			
			case 'j': {
				if (optarg!=NULL){
					cout<<"max Vertex-z is set to "<<optarg<<" cm"<<endl;
					V_z_max=atof(optarg);
				}
				else{
					printf("max Vertex-z without value\n");
					cout<<"default value (zero) will be used"<<endl;
					}
				break;
			};
			
			case 'k': {
				if (optarg!=NULL){
					cout<<"target diameter in XY in cm: "<<optarg<<endl;
					target_diameter=atof(optarg);
				}
				else{
					printf("target diameter without value\n");
					cout<<"default value will be used"<<endl;
					}
				break;
			};
			
			
			case 'l': {
				if (optarg!=NULL){
					cout<<"outname (name of output file) is set to "<<optarg<<endl;
					outputFileName=(string)optarg;
				}
				else{
					printf("outname (name of output file) without value\n");
					cout<<"default value will be used"<<endl;
					}
				break;
			};
			
			case 'p': {
				if (optarg!=NULL){
					cout<<"All other parameters are set to default value "<<endl;
				}
				else{
					cout<<"All other parameters are set to default value "<<endl;
					}
				break;
			};
			
			case 'r': {
				if (optarg!=NULL){
					seed_value=atoi(optarg);
					cout<<"seed is set to "<<optarg<<endl;
				}
				else{
					cout<<"random seed (sys time in microseconds) is used "<<endl;
					}
				break;
			};
			
			case 's': {
				if (optarg!=NULL){
					string strLam1520 = (string)optarg;
					if (strLam1520 == "yes"){
						isLam1520 = true;
						cout<<"lambda mass: "<<massLambda1520<<" MeV"<<endl;
					}
					
				}
				else{
					}
				break;
			};
			
			case '?': default: {
				printf("found unknown option\n");
				break;
			};
		};
	};
    
  if (!(channelName=="KLambda" || channelName=="KSigma" || channelName=="Pi0P" || channelName=="Pi0P_2g" 
  			|| channelName=="PiN")) {
  cout<<"\nERROR: Chanel name is wrong. Please use option --channel=KLambda or (KSigma,Pi0P and PiN)\n ";
  return 1;
  }
  cout<<'\n'<<"End of option reading\n";
    
    
    
    
    
    
    

////////////////////////////////////////////
/*    if(argc == 10) {
        cout << "Running with docker options" << endl;
    	channelName=argv[4];  
	    Ebeam=atof(argv[5]);
	    Q2min=atof(argv[6]);
	    Q2max=atof(argv[7]);
	    Wmin=atof(argv[8]);
	    Wmax=atof(argv[9]);
	    nEventMax=atoi(argv[2]);
	    outputFileName="genKYandOnePion.dat";
    } else if(argc == 9){
    	channelName=argv[1];  
	    Ebeam=atof(argv[2]);
	    Q2min=atof(argv[3]);
	    Q2max=atof(argv[4]);
	    Wmin=atof(argv[5]);
	    Wmax=atof(argv[6]);
	    nEventMax=atoi(argv[7]);
	    outputFileName=argv[8];  
    }*/

    if(getenv("DataKYandOnePion") != NULL)
        	dataPath=getenv("DataKYandOnePion");
    else {
        cerr << "ERROR! Set DataKYandOnePion environment variable" << endl;
        return 1;
    }
//cout<<"dataPath: "<<dataPath<<endl;

	cout << "\nEvent generator started. " <<  endl;


	

/*
	// read config. file
	ifstream input;
	string fName = argv[1];
	input.open(fName.c_str());
	if(!input.good()) {
	  cerr << " eg_ky.cpp: can not open file: " << fName.c_str() << endl;
	  cerr << " STOP!" << endl;
	  return 1;
	}
	std::getline(input,channelName); 
	cout << " Configuration file is " << fName << endl;
	input >> Ebeam >> Q2min >> Q2max >> Wmin >> Wmax >> nEventMax;
*/

	if (V_z_min>V_z_max) swap(V_z_min,V_z_max);

	cout<<" Input parameters: "<<endl;
	cout << " Channel is " << channelName << endl;
	
	
	struct timeval tv;
	gettimeofday(&tv,NULL);
	unsigned long long rand_start = tv.tv_usec + 1000000ull*tv.tv_sec;
	gRandom = new TRandomMT64();
	gRandom->SetSeed(rand_start);
	cout << " Random seed: " << rand_start << endl;
	cout << " Random Engine = " << gRandom->ClassName() << endl;
	/*
	long long unsigned max_llint = 0;
	max_llint = ~max_llint;
	cout << "max value: " <<  max_llint << endl;
	*/
	
	cout << " Ebeam is " << Ebeam << " Gev"<<endl;
	cout << " Q2min is " << Q2min << " Gev2"<<endl;
	cout << " Q2max is " << Q2max << " Gev2"<<endl;
	cout << " Wmin is " << Wmin << " Gev"<<endl;
	cout << " Wmax is " << Wmax << " Gev"<<endl;
	
	cout << " Vertex-z min: " << V_z_min << " cm"<<endl;
	cout << " Vertex-z max: " << V_z_max << " cm"<<endl;
	cout << " target diameter (XY plane): " << target_diameter << " cm"<<endl;
	
	cout << " nEvents is " << nEventMax << endl;
//	input.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
//	std::getline(input,outputFileName); outputFileName.erase(remove(outputFileName.begin(),outputFileName.end(),' '),outputFileName.end());
	cout << " outputFileName is " << outputFileName << endl;
//	std::getline(input,dataPath); dataPath.erase(remove(dataPath.begin(),dataPath.end(),' '),dataPath.end());
	cout << " dataPath is " << dataPath << endl;
	cout<<endl;
//	input.close();


	if(Wmax>4){
		Wmax=4; 
		cout<<"the maxium W value, that can be used in this EG is 4 GeV. So the Wmax is set to be 4 GeV"<<endl;
		cout << "!!! new Wmax is " << Wmax << " Gev"<<endl;
		cout<<endl;
	}
	
	bool pion_decay = (channelName == "Pi0P_2g") ? true : false;
	
	string channelName_for_test = (channelName == "Pi0P_2g") ? "Pi0P" : channelName;
	
	bool check_in_data = check_input_data(dataPath,channelName_for_test,Ebeam,Q2min,Q2max,Wmin,Wmax,nEventMax);
	if (check_in_data==0) {return 0;}

	// initilize event generator
	evGenerator eg(dataPath,channelName_for_test, Ebeam,  Q2min, Q2max, Wmin, Wmax, rand_start, isLam1520);
	channel=num_chanel(channelName_for_test);

        cout << endl
	     << " eg_ky initialized for channel " << channelName << endl  
	     << " Beam energy: " << Ebeam << endl
	     << " Q2 range: " << Q2min <<"  "<< Q2max << ";  "  
	     << " W range: "  << Wmin  <<"  "<< Wmax << ";  " << endl
	     << " Vertex-z range: "  << V_z_min  <<"  "<< V_z_max << ";  " << endl
	     << "   Number of events " << nEventMax << endl; 
	     if (channelName_for_test == "KLambda"){
	     	if(isLam1520) cout << " Lam is Lam(1520), mass: "<< massLambda1520 << endl;
	     	else cout << " Lam mass: "<< massLambda << endl;
	     
	     }
	     
 
        // 4-momenta of final electron, Kaon and Lambda, proton, pi- in LAB frame
	TLorentzVector Peini( TVector3(0.,0.,Ebeam), sqrt(Ebeam*Ebeam+0.005*0.005) );
	TLorentzVector Ppini( TVector3(0.,0.,0.),    0.9383);
	TLorentzVector Pefin, PK, PL, Ppfin, Ppim, Pgam;
	double Q2,W;
	// output  
	ofstream output(outputFileName.c_str()); 
	// Loop through events
        for (int i=0; i<nEventMax; i++) {	
		
	  // get event. 4-moeg_ky/menta of final state particle. 
	  // Values of Q2 and W are also returned.  
          eg.getEvent(Q2, W, Pefin, PK, PL);//, Ppfin, Ppim, Pgam);


          // output in lund format
	  // header
	  //output << "4 1 1 0 0 0 0 "
	  

	  
	  
	 // V-z calculating:
	 double vz_for_event=V_z_max;
	 if ((V_z_max-V_z_min)>0.01){
	 	vz_for_event=V_z_min + (V_z_max-V_z_min) * 	gRandom->Uniform(0.,1.);
	 					//(static_cast<double> (rand())/static_cast<double>(RAND_MAX));
	 }
	 double vx_event = 0, vy_event = 0;
	 if (target_diameter > 0.){
	 	double rad = target_diameter * gRandom->Uniform(0.,1.); 
	 	//(static_cast<double> (rand())/static_cast<double>(RAND_MAX));
	 	double angle = constantPi2 * gRandom->Uniform(0.,1.);
	 	//(static_cast<double> (rand())/static_cast<double>(RAND_MAX));
	 	vx_event = rad*cos(angle);
	 	vy_event = rad*sin(angle);
	 }
	  
	  
	output << ( (pion_decay) ? "4" : "3") << " 1 1 0 0 0 0 "
	   <<" "<< W <<" "<< Q2 <<" "<< getomega(Q2, W) 
		 << endl;
	  // electron
	  output 
	    << "1 -1 1 " << lundIdElectron << " 0 0 "
	    << Pefin.Px() <<" "<< Pefin.Py() <<" "<< Pefin.Pz() 
	    <<" "<< Pefin.E() <<" 0.0005"
	   // <<" 0 0 0 "
	    << " "<<vx_event<<" "<<vy_event<<" "<<vz_for_event<<" "
	    << endl;

	if (channel==1){
	  // Kaon+
	  output 
	    << "2 1 1 " << lundIdKaonPlus << " 0 0 "
	    << PK.Px() <<" "<< PK.Py() <<" "<< PK.Pz() 
	    <<" "<< PK.E() <<" 0.4936"
	   // <<" 0 0 0 "
	    << " "<<vx_event<<" "<<vy_event<<" "<<vz_for_event<<" "
	    << endl;
	    
	    double tmpMass = (isLam1520) ? massLambda1520 : massLambda;
	    string strMass = to_string(tmpMass).substr(0, 6);

	  output 
	   // << "3 0 1 " << lundIdSigmaZero << " 0 0 "
		<< "3 0 1 " << lundIdLambda << " 0 0 "
	    << PL.Px() <<" "<< PL.Py() <<" "<< PL.Pz() 
	    <<" "<< PL.E() <<" " << strMass
	   // <<" 0 0 0 "
	    << " "<<vx_event<<" "<<vy_event<<" "<<vz_for_event<<" "
	    << endl;
	}

	if (channel==2){

	  output 
	    << "2 1 1 " << lundIdKaonPlus << " 0 0 "
	    << PK.Px() <<" "<< PK.Py() <<" "<< PK.Pz() 
	    <<" "<< PK.E() <<" 0.4936"
	   // <<" 0 0 0 "
	    << " "<<vx_event<<" "<<vy_event<<" "<<vz_for_event<<" "
	    << endl;

	  output 
	    << "3 0 1 " << lundIdSigmaZero << " 0 0 "
	    << PL.Px() <<" "<< PL.Py() <<" "<< PL.Pz() 
	    <<" "<< PL.E() <<" 1.192"
	   // <<" 0 0 0 "
	    << " "<<vx_event<<" "<<vy_event<<" "<<vz_for_event<<" "
	    << endl;
	}

	if (channel==3){
		if (pion_decay){
		
			TLorentzVector gamma1, gamma2;
            getPi0decayProd(PK, gamma1, gamma2, gRandom);
		
		output 
	    << "2 0 1 " << 22 << " 0 0 "
	    << gamma1.Px() <<" "<< gamma1.Py() <<" "<< gamma1.Pz() 
	    <<" "<< gamma1.E() <<" 0.0"
			  << " "<<vx_event<<" "<<vy_event<<" "<<vz_for_event<<" "
	    << endl;
	    
	  output 
	    << "3 0 1 " << 22 << " 0 0 "
	    << gamma2.Px() <<" "<< gamma2.Py() <<" "<< gamma2.Pz() 
	    <<" "<< gamma2.E() <<" 0.0"
			  << " "<<vx_event<<" "<<vy_event<<" "<<vz_for_event<<" "
	    << endl;
	    
	  output 
	    << "4 1 1 " << lundIdProton << " 0 0 "
	    << PL.Px() <<" "<< PL.Py() <<" "<< PL.Pz() 
	    <<" "<< PL.E() <<" 0.9382"
			  << " "<<vx_event<<" "<<vy_event<<" "<<vz_for_event<<" "
	    << endl;
		
		}else{
			output 
			  << "2 0 1 " << lundIdPiZero << " 0 0 "
			  << PK.Px() <<" "<< PK.Py() <<" "<< PK.Pz() 
			  <<" "<< PK.E() <<" 0.134"
			 // <<" 0 0 0 "
			  << " "<<vx_event<<" "<<vy_event<<" "<<vz_for_event<<" "
			  << endl;
			  
			 output 
			  << "3 1 1 " << lundIdProton << " 0 0 "
			  << PL.Px() <<" "<< PL.Py() <<" "<< PL.Pz() 
			  <<" "<< PL.E() <<" 0.9382"
			 // <<" 0 0 0 "
			  << " "<<vx_event<<" "<<vy_event<<" "<<vz_for_event<<" "
			  << endl;
		}
	}

	if (channel==4){

	  output 
	    << "2 1 1 " << lundIdPiPlus << " 0 0 "
	    << PK.Px() <<" "<< PK.Py() <<" "<< PK.Pz() 
	    <<" "<< PK.E() <<" 0.1395"
	   // <<" 0 0 0 "
	    << " "<<vx_event<<" "<<vy_event<<" "<<vz_for_event<<" "
	    << endl;

	  output 
	    << "3 0 1 " << lundIdNeutron << " 0 0 "
	    << PL.Px() <<" "<< PL.Py() <<" "<< PL.Pz() 
	    <<" "<< PL.E() <<" 0.939"
	   // <<" 0 0 0 "
	    << " "<<vx_event<<" "<<vy_event<<" "<<vz_for_event<<" "
	    << endl;
	}

  
	  
	  if( i % 10000 == 0) {
	    cout << " Event # " << i << endl; 
	  }
	}
	output.close();
	

	return 0;
	
	
}
