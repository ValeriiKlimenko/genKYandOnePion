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


using namespace std;


int main(int argc, char *argv[])
{

    

    int channel;
	string channelName, outputFileName,dataPath;
	double Ebeam, Q2min, Q2max, Wmin, Wmax;
	int nEventMax;
    double jr, mr, gr, a12, a32, s12, onlyres;

    if(argc == 10) {
        cout << "Running with docker options" << endl;
    	channelName=argv[4];  
	    Ebeam=atof(argv[5]);
	    Q2min=atof(argv[6]);
	    Q2max=atof(argv[7]);
	    Wmin=atof(argv[8]);
	    Wmax=atof(argv[9]);
	    nEventMax=atoi(argv[2]);
	    outputFileName="genKYandOnePion.dat";
    } else if (argc < 9 || argc > 9)  { 
        cerr<<endl;
	    cerr << " Enter 8 arguments only eg.\"./genKYandOnePion arg1 arg2 arg3 arg4 arg5 arg6 arg7 arg8\"" << endl;
        cerr<<endl;
	    cerr << " arg1 is \"KLambda\" or \"KSigma\" or \"Pi0P\" or \"PiN\"" << endl;
	    cerr << " arg2 is Ebeam, arg3 is Q2min, arg4 is Q2max, arg5 is Wmin, arg6 is Wmax " << endl;
	    cerr << " arg7 is nEvents, arg8 is outputFileName" << endl;
        cerr<<endl;
	    cerr << " Example: .\"./genKYandOnePion KSigma 11. 2. 11.999 1.5 4.0 5000 lund_KS.lund\"" << endl;
        cerr<<endl;
	  cerr << " STOP!" << endl;
	  return 1;
    } else {
    	channelName=argv[1];  
	    Ebeam=atof(argv[2]);
	    Q2min=atof(argv[3]);
	    Q2max=atof(argv[4]);
	    Wmin=atof(argv[5]);
	    Wmax=atof(argv[6]);
	    nEventMax=atoi(argv[7]);
	    outputFileName=argv[8];  
    }

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
	cout << " Channel is " << channelName << endl;
	cout << " Ebeam is " << Ebeam << " Gev"<<endl;
	cout << " Q2min is " << Q2min << " Gev2"<<endl;
	cout << " Q2max is " << Q2max << " Gev2"<<endl;
	cout << " Wmin is " << Wmin << " Gev"<<endl;
	cout << " Wmax is " << Wmax << " Gev"<<endl;
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
	bool check_in_data=check_input_data(dataPath,channelName,Ebeam,Q2min,Q2max,Wmin,Wmax,nEventMax);
	if (check_in_data==0) {return 0;}

	// initilize event generator
	evGenerator eg(dataPath,channelName, Ebeam,  Q2min, Q2max, Wmin, Wmax);
	channel=num_chanel(channelName);

        cout << endl
	     << " eg_ky initialized for channel " << channelName << endl  
	     << " Beam energy: " << Ebeam << endl
	     << " Q2 range: " << Q2min <<"  "<< Q2max << ";  "  
	     << " W range: "  << Wmin  <<"  "<< Wmax << ";  " << endl
	     << "   Number of events " << nEventMax << endl; 
	     
 
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
	output << "3 1 1 0 0 0 0 "
	         <<" "<< W <<" "<< Q2 <<" "<< getomega(Q2, W) 
		 << endl;
	  // electron
	  output 
	    << "1 -1 1 " << lundIdElectron << " 0 0 "
	    << Pefin.Px() <<" "<< Pefin.Py() <<" "<< Pefin.Pz() 
	    <<" "<< Pefin.E() <<" 0.0005"
	    <<" 0 0 0 "
	    << endl;

	if (channel==1){
	  // Kaon+
	  output 
	    << "2 1 1 " << lundIdKaonPlus << " 0 0 "
	    << PK.Px() <<" "<< PK.Py() <<" "<< PK.Pz() 
	    <<" "<< PK.E() <<" 0.4936"
	    <<" 0 0 0 "
	    << endl;

	  output 
	   // << "3 0 1 " << lundIdSigmaZero << " 0 0 "
	<< "3 0 1 " << lundIdLambda << " 0 0 "
	    << PL.Px() <<" "<< PL.Py() <<" "<< PL.Pz() 
	    <<" "<< PL.E() <<" 1.115"
	    <<" 0 0 0 "
	    << endl;
	}

	if (channel==2){

	  output 
	    << "2 1 1 " << lundIdKaonPlus << " 0 0 "
	    << PK.Px() <<" "<< PK.Py() <<" "<< PK.Pz() 
	    <<" "<< PK.E() <<" 0.4936"
	    <<" 0 0 0 "
	    << endl;

	  output 
	    << "3 0 1 " << lundIdSigmaZero << " 0 0 "
	    << PL.Px() <<" "<< PL.Py() <<" "<< PL.Pz() 
	    <<" "<< PL.E() <<" 1.192"
	    <<" 0 0 0 "
	    << endl;
	}

	if (channel==3){

	  output 
	    << "2 0 1 " << lundIdPiZero << " 0 0 "
	    << PK.Px() <<" "<< PK.Py() <<" "<< PK.Pz() 
	    <<" "<< PK.E() <<" 0.134"
	    <<" 0 0 0 "
	    << endl;

	  output 
	    << "3 1 1 " << lundIdProton << " 0 0 "
	    << PL.Px() <<" "<< PL.Py() <<" "<< PL.Pz() 
	    <<" "<< PL.E() <<" 0.9382"
	    <<" 0 0 0 "
	    << endl;
	}

	if (channel==4){

	  output 
	    << "2 1 1 " << lundIdPiPlus << " 0 0 "
	    << PK.Px() <<" "<< PK.Py() <<" "<< PK.Pz() 
	    <<" "<< PK.E() <<" 0.1395"
	    <<" 0 0 0 "
	    << endl;

	  output 
	    << "3 0 1 " << lundIdNeutron << " 0 0 "
	    << PL.Px() <<" "<< PL.Py() <<" "<< PL.Pz() 
	    <<" "<< PL.E() <<" 0.939"
	    <<" 0 0 0 "
	    << endl;
	}



/*
	  output 
	    << "2 1 1 " << lundIdKaonPlus << " 0 0 "
	    << PK.Px() <<" "<< PK.Py() <<" "<< PK.Pz() 
	    <<" "<< PK.E() <<" 0.4936"
	    <<" 0 0 0 "
	    << endl;
	  if(channel==1 || channel==2 || channel==4){
	
	  // Lambda
	      //output 
	      //  << "3 0 1 " << lundIdLambda << " 0 0 "
	      //  << PL.Px() <<" "<< PL.Py() <<" "<< PL.Pz() 
	      //  <<" "<< PL.E() <<" 1.115"
	      //  <<" 0 0 0 "
	      //  << endl;
	  output 
	    << "3 1 1 " << lundIdProton << " 0 0 "
	    << Ppfin.Px() <<" "<< Ppfin.Py() <<" "<< Ppfin.Pz() 
	    <<" "<< Ppfin.E() <<" 0.9383"
	    <<" 0 0 0 "
	    << endl;

	
	  output 
	    << "4 -1 1 " << lundIdPiMinus << " 0 0 "
	    << Ppim.Px() <<" "<< Ppim.Py() <<" "<< Ppim.Pz() 
	    <<" "<< Ppim.E() <<" 0.1396"
	    <<" 0 0 0 "
	    << endl;
	


	  } else if(channel==3) {
	  // Sigma zero
	  output 
	    << "2 0 1 " << lundIdSigmaZero << " 0 0 "
	    << PL.Px() <<" "<< PL.Py() <<" "<< PL.Pz() 
	    <<" "<< PL.E() <<" 1.192"
	    <<" 0 0 0 "
	    << endl;	  
	  }
	*/  
	  // test 
	  //TLorentzVector PfinSum = Pefin + PK + Ppfin + Ppim;
	  //cout << " Pfin " << PfinSum.E() <<" "<< PfinSum.Px() <<" "<< PfinSum.Py() <<" "<< PfinSum.Pz()
	  //     << endl; 
	  //cout << " miss 0 "  << (Peini + Ppini - Pefin - PK - Ppfin - Ppim).M() << endl;
	  //cout << " miss e "  << (Peini + Ppini - PK - Ppfin - Ppim).M() << endl;
	  //cout << " miss p "  << (Peini + Ppini - Pefin - PK - Ppim).M() << endl;
	  //cout << " miss K "  << (Peini + Ppini - Pefin - Ppfin - Ppim).M() << endl;
	  //cout << " miss pi " << (Peini + Ppini - Pefin - PK - Ppfin).M() << endl;
	  //cout << " Mass L "  << (Ppfin + Ppim).M() << endl;
	  
	  
	  if( int((i+1)/1000)*1000 == (i+1) | (i+1)==1) {
	    cout << " Event # " << i+1 << endl; 
	  }
	}
	output.close();
	

	return 0;
	
	
}
