
#include "string.h"
#include <time.h>

#include <TLorentzVector.h>
#include <TRandom.h>

#include "utils.h"
#include "constants.h"
#include "kinematics.h"
//#include "sigmaKY.h"
#include "sigmaValera.h"


class evGenerator {

  string type;
  double m1;
  double m2;  
  double Ebeam;
  double Q2min;
  double Q2max;
  double Wmin;
  double Wmax;
  double d5sigmaMax;
  
  TRandom* gRandom = new TRandomMT64();
  
  // get random number in the interval [min, max].
  double inline randomIntv(double min, double max){

  return gRandom->Uniform(min, max);
  //return min + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(max-min)));
};
	
  Sigma  *model;

  //double cos_min;//only for test
  //double cos_max;//only for test
  
  int nEvent;

public:

evGenerator(string dataPath, string t, double E, 
            double q2min, double q2max,
	    double wmin,  double wmax, unsigned long long rand_start, bool isL1520,
		//double cosmin,  double cosmax,//only for test
	    double jr=-1, double mr=0, double gr=0, 
            double a12=0., double a32=0., double s12=0., 
	    int onlyres=0)
  {
  
  gRandom->SetSeed(rand_start);
  type = t;
  Ebeam = E;
  Q2min = q2min;      
  Q2max = q2max;
  Wmin = wmin;
  Wmax = wmax;
  nEvent = 0;

  //cos_min=cosmin;//only for test
 // cos_max=cosmax;//only for test

    int channel=0;//1-KL 2-KS 3-PiOP 4-PiN

	if(type == "KLambda"){ 
		m1 = massKaon;
		if (isL1520) m2 = massLambda1520;
		else m2 = massLambda;
		channel=1;
	}else{
		if(type == "KSigma"){ 
			m1 = massKaon; 
    			m2 = massSigma0;
			channel=2;
		}else{
			if(type == "Pi0P"){ 
			m1 = massPion0;
    			m2 = massProton;
			channel=3;
			cout<<"chanel:"<<channel<<endl;
			}else{
				if(type == "PiN"){ 
				m1 = massPion; 
    				m2 = massNeutron;
				channel=4;
				}else{
					cerr << " Error! evGenerator::evGenerator Wrong reaction type!";
					cerr << " It must be KLambda or KSigma or Pi0P or PiN.";  
					exit(1);
				}
			}
		}
	}

  // Initilize: read data file
  //model = new sigmaKY(type);

  model = new Sigma(dataPath,channel);


  // Find maximum of the cross section
  if(type == "KLambda" || type == "KSigma" || type == "Pi0P" || type == "PiN") {
    //cout<<"Ebeam=: "<< Ebeam<<" Q2min: "<<Q2min<<" Q2max: "<<Q2max<<" Wmin: "<<Wmin<<" Wmax: "<<Wmax<<endl;
    d5sigmaMax = model->d5sigma_max(Ebeam, Q2min, Q2max, Wmin, Wmax);
    cout<<" d5sigmaMax: "<<d5sigmaMax<<endl;
	if(d5sigmaMax==0){
		cerr << " incorrect kinematic region, pls check input Q2 and Energy";  
		exit(1);
	}
  }
	//cout<<"2"<<endl;
  // initialize random seed: 
  srand (time(NULL));
 
};

~evGenerator(){
    delete model;
};




		 
		 
// all parameters are output parameters
// 4-momenta of particles in the final state in LAB frame:
//    Pefin, PK, PL -> electron, Kaon, Lambda/Sigma
// 4 momenta of 

void getEvent(double &Q2, double &W, 
              TLorentzVector &Pefin, TLorentzVector &PK, TLorentzVector &PY)//,
	   //   TLorentzVector &Ppfin, TLorentzVector &Ppim, TLorentzVector &Pgam) //only for test
										{

  int nTry=0;
  nEvent=0;
  while(true) {
 
    Q2 = randomIntv(Q2min, Q2max);
    W  = randomIntv(Wmin, Wmax);

    //double cosThetaK=randomIntv(cos_min, cos_max);//only for test
    double cosThetaK = randomIntv(-0.999999,0.999999);//real
//cout<<"cosThetaK="<<cosThetaK<<endl;
    double thetaK = acos(cosThetaK);
    double phiK = randomIntv(0.0, 2*constantPi);
//cout<<"cosThetaK="<<cosThetaK<<endl;
    double d5sigma;
    if(type == "KLambda" || type == "KSigma"|| type == "Pi0P"|| type == "PiN") {
//cout<<"Q2: "<<Q2<<" W: "<<W<<" thetaK: "<< thetaK<<" phiK: "<< phiK<<endl;
      d5sigma = model->d5sigma(Ebeam, Q2, W, thetaK, phiK);
    } else {
      cout <<" ::getEvent " <<  "Wrong reaction type " << type << endl;
    }
        
    nTry++;

	//if (cos(acos(arg))<0.9) cout<<" ="<<" d5/max= "<<" d5sigma="<<d5sigma<<" d5sigmaMax="<<d5sigmaMax<<" theta: "<<acos(arg)<<" cos: "<<cos(acos(arg))<<endl;

    if(randomIntv(0.,1.) < d5sigma/d5sigmaMax) {

       nEvent++;

       double W2 = W*W;
       double omega = getomega(Q2,W);
       if(omega <=0. ) continue;
       double Ee = Ebeam - omega;
       if(Ee<massElectron) continue;
       double arg = 1. - Q2/(2.*Ebeam*Ee);
       if( fabs(arg) > 1. ) continue;
       double theta = acos(arg);
       double phi = randomIntv(0.0, 2*constantPi);
       



       // 4-momentum of final electron in LAB frame
       Pefin.SetXYZT(1.,1.,1.,Ee);
       Pefin.SetRho(sqrt(Ee*Ee - massElectron2)); 
       Pefin.SetTheta(theta);
       Pefin.SetPhi(phi); 
       
       // 4-momenta of K and Lambda/Sigma //проверить м1 и м2
       cms2lab(W, Q2, phi, Ebeam, thetaK, phiK, m1, m2, PK, PY);     

       
       //cout << " PY " << PY.E() <<" "<< PY.Px() <<" "<< PY.Py() <<" "<< PY.Pz() << endl; 
       //cout << " PS " << (Ppfin+Ppim).E() <<" "<< (Ppfin+Ppim).Px() <<" "<< (Ppfin+Ppim).Py() <<" "<< (Ppfin+Ppim).Pz() << endl; 
       //cout << " PK " << (PK).E() <<" "<< (PK).Px() <<" "<< (PK).Py() <<" "<< (PK).Pz() << endl; 
       //cout << " Pp " << (Ppfin).E() <<" "<< (Ppfin).Px() <<" "<< (Ppfin).Py() <<" "<< (Ppfin).Pz() << endl; 
       //cout << " Pi " << (Ppim).E()  <<" "<< (Ppim).Px()  <<" "<< (Ppim).Py()  <<" "<< (Ppim).Pz()  << endl; 
       
       return;
    
    }

  
  }

}; //end getEvent(...)

	
};
