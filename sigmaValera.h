#ifndef _CLASS_EV_GEN_H
#define _CLASS_EV_GEN_H

#include <vector>
#include "utils.h"
#include "kinematics.h"
#include <stdexcept>

enum channel{
	KLambda = 0,
	KSigma = 1,
	Pi0P = 2,
	piN = 3
};

struct CS_at_max_Q2{
	double Q2;
	double BeamEner; // Has to be checked, it should be taken from an. notes
	double cosTh, p0,p1,p2;
	double cr1, cr2; // two correcion factors
};

struct f1_params{
	double W, p0,p1,p2;
};

double naive_lerp(const double a,const  double b,const  double t)
{
    return a + t * (b - a);
}

bool check_kin(double Q,double W, double Ebeam)
{
 double MP=massProton;
 double OMEGA=(W*W+Q-MP*MP)/(2*MP);
 double sinus = sqrt(Q/(4*Ebeam*(Ebeam-OMEGA)));
 double TETE=2*asin(sqrt(Q/(4*Ebeam*(Ebeam-OMEGA))));
 double EE2=Ebeam*(Ebeam-getomega(Q,W));
 if  ( ((Q/(4*EE2))<=1)&& ((2-Q/(2*EE2)) >0) && ((-Q/(2*EE2))<=0)&&
  (OMEGA>=0)&&(Ebeam-OMEGA>=0)&&(abs(sinus)<=1)&&(abs(TETE)<=3.1415927/2) )
 {return 1;}
 else return 0;
}

bool check_input_data(string dataPath, string __channelName,double __Ebeam, double __Q2min, double __Q2max, 
						double __Wmin, double __Wmax, int __nEventMax){
 if ((__channelName!="KLambda")&&(__channelName!="KSigma")&&(__channelName!="PiN")&&(__channelName!="Pi0P")){
  cout<<"incorrect name of chanel, pls try: KLambda or KSigma or Pi0P or PiN"<< endl;
  return 0;}
 if ((__Ebeam<0)||(__Ebeam>30)){
  cout<<"incorrect Energy, right value is Energy>0 and Energy <30 "<< endl;
  return 0;}
 if (__nEventMax<0){
  cout<<"incorrect nEventMax, right value >0 "<< endl;
  return 0;}
 if (__Q2max<__Q2min){
  cout<<"incorrect Q2_max or Q2_min"<< endl;
  return 0;}
 if ((__Q2max>30)||(__Q2max<0)){
  cout<<"incorrect Q2_max, Q2_max<30"<< endl;
  return 0;}
 if (__Wmax>5.){cout<<"incorrect W_max, try less than 5 GeV ;"<< endl; return 0;}


 ifstream interp_right;
 interp_right.open(dataPath+"/KL_interp.txt");
//cout<<" wau:"<<dataPath+"/KL_interp.txt"<<endl;
 if (!interp_right.is_open()) { cout<<"ERROR: wrong way to the data"<<endl;return 0;}

 return 1;

 
}

class Sigma{
	private:
/////////////////////////////data 1:///////////////////
 vector<double> Q_Qmax,W_Qmax,costeta_Qmax, p0_Qmax,p1_Qmax,p2_Qmax,param;
 vector<double> Q_F1, W_F1,F1_F1;
 vector<double> CS_inter,W_inter;
 vector<double> costeta_Qmax_int;

 double max_W_F1=0,Q_max_channel=-1,max_W_Qmax=0,max_W_inter=0;
 double min_W_inter=100000, min_W_F1=10000,min_W_Qmax=10000,ph_fac=1;
 double max_cos,min_cos;
 int num_str4=0,n_str_Ev=0,num_str3=0,num_str2=0,num_str=0,type_chanel=0,num_costeta=-1,m=0;

 int range_fi=-1;//1: -180 +180 2: 0 360
 int range_cos=-1;//1: == 2: n !=
////////////////////////Q2 extrapolation 5-12 GeV2 functions:////////////////////////

/// updated 04/21/2022: ////

		vector<pair<double, vector<CS_at_max_Q2>>> vCSmaxQ2;
		vector<f1_params> vF1;
		
		void read_maxQ2(const string& dataPath, const channel& channelName){
			const vector<string> files_fitted_cs_at_maxQ2 = {"/KLambda_CS_Qmax_Fit_gladk.txt",
														"/KSigma_CS_Qmax_Fit_gladk.txt",
														"/Pi0P_CS_Qmax_Fit_gladk.txt",
														"/PiN_CS_Qmax_Fit_gladk.txt"};
			ifstream f_Q2max;
			const string filePath = dataPath + files_fitted_cs_at_maxQ2.at(static_cast<int>(channelName));
			f_Q2max.open(filePath);
			
			if (!f_Q2max.is_open())
				throw invalid_argument("file is not open: " + filePath);
			
			//Get rid of header
			string tmpStr;
			getline(f_Q2max, tmpStr);
    		getline(f_Q2max, tmpStr);
			
			
			double previousW = 0;
			bool firstEnter = true;
			vector<CS_at_max_Q2> tmpVbuff;
			while(f_Q2max.good()){
				double tmpW = 0, tmpQ2 = 0, tmpBE = 0,  tmpCosTh = 0;
				double tmpP0 = 0, tmpP1 = 0, tmpP2 = 0;
				double tmpCr1 = 0, tmpCr2 = 0;
				f_Q2max >> tmpW;
				f_Q2max >> tmpQ2;
				f_Q2max >> tmpBE;
				f_Q2max >> tmpCosTh;
				
				//cout<<tmpW << ' ' << tmpQ2 << ' ' <<tmpBE << ' ' << tmpCosTh << endl;;
				//check that data is reasonable:
				if (tmpW < 1. || tmpW > 5. || tmpQ2 < 0. || tmpQ2 > 30. || tmpCosTh < -1. || tmpCosTh > 1.)
					continue;
				
				f_Q2max >> tmpP0;
				f_Q2max >> tmpP1;
				f_Q2max >> tmpP2;
				
				f_Q2max >> tmpCr1;
				f_Q2max >> tmpCr2;
				double fullCorr = tmpCr1 * tmpCr2;
				//check fullCorr
				
				//cout << tmpW << ' ';
				if (abs(tmpW - previousW )< 0.00001 || firstEnter){
								//cout << '1' << endl;
					tmpVbuff.push_back({tmpQ2, tmpBE, tmpCosTh, tmpP0 * fullCorr, tmpP1 * fullCorr, tmpP2 * fullCorr, tmpCr1, tmpCr1});
				}else{
								//cout << '2' << endl;
					vCSmaxQ2.push_back({previousW, tmpVbuff});
					tmpVbuff.clear();
					tmpVbuff.push_back({tmpQ2, tmpBE, tmpCosTh, tmpP0 * fullCorr, tmpP1 * fullCorr, tmpP2 * fullCorr, tmpCr1, tmpCr1});
				}
				previousW = tmpW;
				firstEnter = false;
			}
			
			
			if ((previousW - vCSmaxQ2.back().first) > 0.001)
				vCSmaxQ2.push_back({previousW, tmpVbuff});
				
			if (vCSmaxQ2.size() < 1)
					throw invalid_argument("something wrong (1) with file format: " + filePath);
			
		};
		
		void read_StrFunF1(const string& dataPath, const channel& channelName){
			const vector<string> files_f1 = {"/KL_Fit_F1_Q2_30.dat",
											"/KS_Fit_F1_Q2_30.dat",
											"/Pi0P_Fit_F1_Q2_30.dat",
											"/PiN_Fit_F1_Q2_30.dat"};
			ifstream f_f1;
			const string filePath = dataPath + "/highQ2" + files_f1.at(static_cast<int>(channelName));
			f_f1.open(filePath);
			
			if (!f_f1.is_open())
				throw invalid_argument("file is not open: " + filePath);
			
			//Get rid of header
			string tmpStr;
			getline(f_f1, tmpStr);
				
			while(f_f1.good()){
				double tmpW = 0,tmpP0 = 0, tmpP1 = 0, tmpP2 = 0;
				f_f1 >> tmpW;
				
				//cout<<tmpW << ' ' << tmpQ2 << ' ' <<tmpBE << ' ' << tmpCosTh << endl;;
				//check that data is reasonable:
				if (tmpW < 1. || tmpW > 5.)
					continue;
				
				f_f1 >> tmpP0;
				f_f1 >> tmpP1;
				f_f1 >> tmpP2;
				
				vF1.push_back({tmpW, tmpP0, tmpP1, tmpP2});
			}
			if (vF1.size() < 1)
					throw invalid_argument("something wrong (1) with file format: " + filePath);
		};

		double interpol(const double W, const double cosTh, const double phi){
		
			//cout<<W<<' '<< cosTh << ' '<< phi << " type:";
			// naive temp.func re-do with binary search:
			for (size_t iW = 0; iW < vCSmaxQ2.size() - 1; iW++){
				if (W == vCSmaxQ2[iW].first || W == vCSmaxQ2[iW + 1].first){
					if (W == vCSmaxQ2[iW + 1].first) 
						iW++;
					//cout << '1' <<endl;
					return getCSfromVbyWposAndCosThValue(iW, cosTh, phi);
				}
				if ( W > vCSmaxQ2[iW].first && W < vCSmaxQ2[iW + 1].first){
					//cout << '2' <<endl;
					// add both W and cos interp
					double cs_minW = getCSfromVbyWposAndCosThValue(iW, cosTh, phi);
					double cs_maxW = getCSfromVbyWposAndCosThValue(iW+1, cosTh, phi);
					double l = (W - vCSmaxQ2[iW].first) / ( vCSmaxQ2[iW + 1].first - vCSmaxQ2[iW].first);
					return naive_lerp(cs_minW, cs_maxW, l);
				}	
			}
			
			//throw error
			return 0;
		};
		
		double getF1(const double W, const double Q2){
			//Get F1 at smaller W:
			int binN_low = (W - vF1.front().W) * (vF1.size() - 1) / (vF1.back().W  - vF1.front().W);
			double minF1 = vF1.at(binN_low).p0 + vF1.at(binN_low).p1 / Q2 + vF1.at(binN_low).p2 / (Q2 *Q2);
			//cout << (W - vF1.front().W) << ' ' << (vF1.back().W  - vF1.front().W) << ' ' << (vF1.back().W  - vF1.front().W) / (vF1.size() - 1) << " M:" << binN_low << endl;
			if (binN_low == vF1.size() - 1)
				return minF1;
			else{
				int binN = binN_low + 1;
				double maxF1 = vF1.at(binN).p0 + vF1.at(binN).p1 / Q2 + vF1.at(binN).p2 / (Q2 * Q2);
				double l = (W - vF1.at(binN_low).W) * (vF1.size() - 1) / (vF1.back().W  - vF1.front().W);
				//cout<<" W:" << W << ' ' <<  minF1 << ' ' << maxF1 << ' ' <<  l << ' ' << binN_low << endl;
				return naive_lerp(minF1, maxF1, l);
			}
		}
		
		//legacy name from sigmaValera.h
		double getCS_fit(const double Ebeam,const  double Q2, const double Q2_max,const double W){
			//F1 and CS at max data Q2:
			double sigma_t_max = getF1(W, Q2_max) * 4 * constantPi2 *
									constantAlpha / (getK(Q2_max,W) * massProton);
			double sigma_l_max = sigma_t_max * 0.2;
			double sigma_max = sigma_t_max + getEpsilon(Ebeam,Q2_max,W) * sigma_l_max;
			//F1 and CS at current Q2:
			double sigma_t = getF1(W, Q2) * 4 * constantPi2 *
									constantAlpha / (getK(Q2,W) * massProton);
			double sigma_l = sigma_t * 0.2;
			double sigma = sigma_t + getEpsilon(Ebeam,Q2,W) * sigma_l;
			return sigma / sigma_max;
		};
		
		
		double getCSfromVbyPos(const size_t iW, const size_t iCosTh, const double phi){
			return vCSmaxQ2[iW].second[iCosTh].p0 +
					vCSmaxQ2[iW].second[iCosTh].p1 * cos(2*phi/57.29578049) +
					vCSmaxQ2[iW].second[iCosTh].p2 * cos(phi/57.29578049);
		}
		
		double getCSfromVbyWposAndCosThValue(const size_t iW, const double cosTh, const double phi){
			//cout<<"ch W:" << vCSmaxQ2[iW].first<<endl;
			for (size_t iCosTh = 0; iCosTh < vCSmaxQ2[iW].second.size() - 1; iCosTh++){
			
				if (cosTh == vCSmaxQ2[iW].second[iCosTh].cosTh || cosTh == vCSmaxQ2[iW].second[iCosTh +1].cosTh){
					if (cosTh == vCSmaxQ2[iW].second[iCosTh + 1].cosTh) iCosTh++;
						return getCSfromVbyPos(iW, iCosTh, phi);
				}
			
				// cos intrep
				if (cosTh > vCSmaxQ2[iW].second[iCosTh].cosTh && cosTh < vCSmaxQ2[iW].second[iCosTh + 1].cosTh){
					double l = (cosTh - vCSmaxQ2[iW].second[iCosTh].cosTh) / 
								(vCSmaxQ2[iW].second[iCosTh + 1].cosTh - vCSmaxQ2[iW].second[iCosTh].cosTh);
					return naive_lerp(getCSfromVbyPos(iW, iCosTh, phi), getCSfromVbyPos(iW, iCosTh + 1, phi), l);
				}
			}
			
			//throw error
			return 0;
		}

/////////// old 2018: ///////////////
 double getK(double Q, double W);
 double fun_points(double Q, double W, const vector<double>& Q_F1, const vector<double>& W_F1, const vector<double>& F1_F1,int num_str);
 double anti_Fit(double fi, double p0,double p1,double p2)
		{return p0+p1*cos(2*fi/57.29578049)+p2*cos(fi/57.29578049);}

 double lineal_interp_1(double x_r, const vector<double>& x,const vector<double>& y,int num_str);

 //double interpol(double Q, double W, double fi, const vector<double>& Q_F1, const vector<double>& W_F1,
 //   const vector<double>& p0_Qmax, const vector<double>& p1_Qmax, const vector<double>& p2_Qmax,double num_str);

 //double getCS_fit(double Ebeam, double Q, double Qmax,double W, 
 //  const vector<double>& Q_F1, const vector<double>& W_F1, const vector<double>&F1_F1,double num_str);
			//,vector<double> Q_F2, vector<double> W_F2, vector<double> F2_F2,double num_str2);
 double check_input_param(double Q,double W, double Ebeam);
 double check_cos(double costeta, double W);
 double get_CS_int_fi(double Q, double W, double costeta, double Ebeam);
 double get_CS_f_i(double Q,double W, double costeta, double fi, double Ebeam);
///////////////////////////////interpolation Q2 0-5 GeV2 functions://///////////////////////////////////////////////////////////////////
//data 2:
 vector<double> _Q2,_W,_cos,_p0,_p1,_p2,param_vec,_Q_int,_W_int,_CS_int;
 vector<double> _CS_ph,_W_ph,W_vec_ph,CS_vec_ph,costeta_vec_ph;
 vector<double> _CS_ph_Ev, _W_ph_Ev;
 vector<double> _W_max,_W_min,_Q_for_ext_point;
 int is_there_glad=0,num_ext_p=0,n_str_CS=0,n_str_ph=0,n_str_ph_int=0,n_str_CS_int=0,interp_num1=0,interp_num2=0;
 double Qmin=10000,Qmax=0,max_W=0,min_W=10000, max_W_ph=0,min_W_ph=10000,_W_min_all=1000,_W_max_all=0,W_ext_min=1000,W_ext_max=0;
//photo func main:
 double get_CS_ph_int(double W);// photo CS as func W
 double get_CS_ph(double W, double cos);// photo CS as func W and cos
//photo func medium
 double ph_int(double W, double cos);
 double ph_int_int_Ev(double W);
 double ph_int_int(double W);
 double get_CS_ph_int_Ev(double W);
 double cos_in_ph(int sp, double cos);
//prepration to work
 void search_exterm_points();
//internal func:
 double lin_interp(double x,double point1, double point2, double value_point1, double value_point2);
 double anti_Fit(double fi, double p0,double p1,double p2, int val);//anti_fit and integr po fi 1-anti_fit, 2-untegr po fi
 double intrep_CS(double Q,double W,double cos,double fi, int type_CS);
 double W_in(int sp,double W, double cos,double fi, int type);
 double cos_in(int sp, double cos,double fi, int type);
 int check_possibil_inter_W(double Q, double W);
 int check_possibil_inter_Q2(double Q);
 double change_Q2(double Q);
 int range_ph(double W);
 double intrep_CS_part(double Q,double W,double cos,double fi, int type_CS);
 double W_in_part(int sp,double W, double cos,double fi, int type);
//////////////////////dont use, it should be checked:////////////////////////////////////////////////////////
 double get_d4CS(double Q,double W, double cos, double Ebeam);
 double int_cos(double Q, double W);
 double W_in_cs(int sp, double W);
 double getCS_d3CS(double Q,double W, double Ebeam);
 double get_CS_int_fi_cos(double Q, double W, double Ebeam);//3-dimens Cross Section
//////////////////can be used, if it is nessas../////
 double get_CS(double Q,double W, double costeta, double fi, double Ebeam);//CS in Q2 range from max Q2 achivable experement data(depends on channel) to 12
 double get_d5CS(double Q,double W, double cos, double fi, double Ebeam);//CS in Q2 range from 0 to max Q2 achivable experement data(depends on channel)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	public:
///Cross Sections:
 double d5sigma2(double Ebeam, double Q2, double W, double costhetaK, double phiK);//5-dimens Cross Section, one of the vars is COS(theta)
 double d5sigma(double beam_energy, double _Q2, double _W, 
  		double thetaK, double phiK);//5-dimens Cross Section, one of the vars is theta
 double dsigma_dcos(double _beam_energy, double _Q2, double _W, double teta);//4-dimens Cross Section

 double int_get_d5CS(double Q,double W, double Ebeam);//3-dimens CS

 double d5sigma_max(double Ebeam, double Q2min, double Q2max,
          double Wmin,  double Wmax  );//max CS in that region
////other:
 double porog_ch(int num_chanel);//threshold of the reaction
 Sigma(string dataPath, int chanel);//constr
	protected:
};
///////////////realization:///////////////////////////////////////////////////////////
double Sigma::d5sigma(double Ebeam, double Q2, double W, double thetaK, double phiK){
 bool ch=check_kin(Q2,W,Ebeam);
 if (ch==0) { //cout<<"uncorrect input Q and W"<<endl;
  return 0;}
 if ((W<porog_ch(type_chanel))||(Q2<0.0001)||(Q2>30)||(W>5)||(phiK<0)||(phiK>6.284)) return 0;
 double mp=massProton;
 double mp2 = mp*mp;
 double pi=constantPi;
 double alpha=constantAlpha;
 double W2 = W*W;
 double omega = (W2 + Q2 -mp2)/(2.*mp);
 if(omega<0. || omega>Ebeam) return 0.;
 double eps2 = Ebeam - omega;
 double k2 = Q2 + omega*omega;
 if(k2<0) return 0.;
 double arg = 1. - Q2/(2.*Ebeam*eps2);
 if(arg>1. || arg<-1.) return 0.;
 double theta = 2.*acos(arg);
 double epsilon = 1./( 1 + (2.*k2/Q2)*tan(theta/2.) );
 double gamma = getGamma(Ebeam, Q2, W);
 double d5sig=0;

 if (Q2<Q_max_channel){
    d5sig = gamma * 
	this->get_d5CS(Q2, W, cos(thetaK), phiK*57.2957, Ebeam);
 }else  d5sig = gamma *
           this->get_CS(Q2, W, cos(thetaK), phiK*57.2957, Ebeam);
  return d5sig;
}

double Sigma::d5sigma2(double Ebeam, double Q2, double W, double costhetaK, double phiK){
 bool ch=check_kin(Q2,W,Ebeam);
 if (ch==0) { //cout<<"uncorrect input Q and W"<<endl;
  return 0;}
 if ((W<porog_ch(type_chanel))||(Q2<0.0001)||(W>5)||(phiK<0)||(phiK>6.284)||(costhetaK>1)||(costhetaK<-1)) return 0;
  double mp=massProton;
  double mp2 = mp*mp;
  double pi=constantPi;
  double alpha=constantAlpha;
  double W2 = W*W;
  double omega = (W2 + Q2 -mp2)/(2.*mp);
  if(omega<0. || omega>Ebeam) {
   return 0.;}
  double eps2 = Ebeam - omega;
  double k2 = Q2 + omega*omega;
  if(k2<0) {
	return 0.;}
  double arg = 1. - Q2/(2.*Ebeam*eps2);
  if(arg>1. || arg<-1.) { 
		return 0.;}
  double theta = 2.*acos(arg);
  double epsilon = 1./( 1 + (2.*k2/Q2)*tan(theta/2.) );
  double gamma = getGamma(Ebeam, Q2, W);
  double d5sig=0;

 if (Q2<Q_max_channel){ 
   d5sig = gamma * 
			this->get_d5CS(Q2, W, costhetaK, phiK*57.2957, Ebeam);
 }else{
  d5sig = gamma *
		this->get_CS(Q2, W, costhetaK, phiK*57.2957, Ebeam);
}    
  return d5sig;
}

double Sigma::d5sigma_max(double Ebeam, double Q2min, double Q2max,
                                       double Wmin,  double Wmax  )
{
        double pi=constantPi;

	// Find maximum of the cross section

        int nQ2=12;
        int nW=40;
        int nCosThetaK=60;
        int nPhiK=60;
	double d5sigmaMax=0.;

Q2max=Q2max-(Q2max-Q2min)/2;
        for(int iQ2=0; iQ2<nQ2; iQ2++) {
	double Q2 = Q2min + (Q2max-Q2min)*iQ2/(nQ2-1);
	cout<<"Find maximum of the cross section: "<<100*iQ2/nQ2<<"%"<<endl;
        for(int iW=0;   iW<nW;  iW++) {
	double W =  Wmin +  (Wmax-Wmin)*iW/(nW-1);
	for(int iCosThK=0; iCosThK<nCosThetaK; iCosThK++) {
	//double cosThetaK=0;
	double cosThetaK = -0.9999 + (0.9999-(-0.9999))*iCosThK/(nCosThetaK-1);
	//cout<<"costeta="<<cosThetaK<<endl;
	double thetaK = acos(cosThetaK);
        for(int iPhiK=0; iPhiK<nPhiK;  iPhiK++) {
	double phiK = 0. + (2.*pi-0.)*iPhiK/(nPhiK-1);
	//cout<<" Ebeam="<<Ebeam<<" Q2="<<Q2<<" W="<<W<<" thetaK="<<thetaK<<" cosThetaK: "<<cosThetaK<<" phiK="<<phiK<<endl; 
          double d5sig = d5sigma(Ebeam, Q2, W, thetaK, phiK);
		//cout<<"d5sig: "<<d5sig<<endl;
          if(d5sig>d5sigmaMax) d5sigmaMax=d5sig;
    
	}
        }
	}
        }

if (d5sigmaMax==0){cout<<"incorrect kinematic region, pls check input Q2 and Energy"<<endl;}
	cout<<"Find maximum of the cross section: "<<"complete"<<endl;
        return d5sigmaMax;

}
	double Sigma::get_CS_f_i(double Q,double W, double costeta, double fi, double Ebeam){

			double real_W=W;
			W=check_input_param(Q,W,Ebeam);
			if (W==0) {cout<<"pls check input Q and W"<<endl;  return 0;}
			bool ch=check_kin(Q,W,Ebeam);
			if (ch==0) {
				cout<<"uncorrect input Q and W"<<endl;
				return 0;
			}


			double new_costeta=check_cos(costeta,W);
			if (ch==1){
			//cout<<"W="<<W<<endl;
			

			
			double CS=interpol(W,new_costeta,fi);
			//double CS=interpol(W,new_costeta,fi,W_Qmax,costeta_Qmax,p0_Qmax,p1_Qmax,p2_Qmax,num_str3);
			getCS_fit(Ebeam,Q,Q_Qmax[2],W);
			//getCS_fit(Ebeam,Q,Q_Qmax[2],W,Q_F1,W_F1,F1_F1,num_str);//,Q_F2,W_F2,F2_F2,num_str2);

			if (real_W<W){
				double an_dw=1/(W-porog_ch(type_chanel));
				double CS2=CS*an_dw*real_W-CS*porog_ch(type_chanel)*an_dw;
				CS=CS2;
			}

			if (real_W>W){
				if (real_W>=max_W_inter){

					double real_W=max_W_inter;

					if (W<min_W_inter) {W=min_W_inter;}
					double CS3=CS*
					lineal_interp_1(real_W,W_inter,CS_inter,num_str4)/
								lineal_interp_1(W,W_inter,CS_inter,num_str4);
					CS=CS3;


				}else{

					if (real_W<=min_W_inter){
					}
	 				else{
						if (W<min_W_inter) {W=min_W_inter;}
						double CS3=CS*
						lineal_interp_1(real_W,W_inter,CS_inter,num_str4)/
									lineal_interp_1(W,W_inter,CS_inter,num_str4);
						CS=CS3;
					}				
				}

			}

			return CS;
			}
			return 0.;
		}
		
		double Sigma::get_CS( double Q, double W, double costeta, double fi, double Ebeam)
		{
			double real_W=W;

			W=check_input_param(Q,W,Ebeam);
			if (W==0) {//cout<<"pls check input Q and W"<<endl;  
				return 0;}
			/*bool ch=check_kin(Q,W,Ebeam);
			if (ch==0) {
				//cout<<"uncorrect input Q and W"<<endl;
				return 0;
			}*/
			bool ch=1;

			if (range_fi==1){
				fi=fi-180;
				if (abs(fi+180)<0.001) {fi=-179.9999;}
				if (abs(fi-180)<0.001) {fi=179.9999;}
			}
			double new_costeta=check_cos(costeta,W);
			if (ch==1){
			

			
			 double part1=interpol(W,new_costeta,fi);
			 //double part1=interpol(W,new_costeta,fi,W_Qmax,costeta_Qmax,p0_Qmax,p1_Qmax,p2_Qmax,num_str3);
			 

			 double part2=getCS_fit(Ebeam,Q,Q_Qmax[2],W);
			 //double part2=getCS_fit(Ebeam,Q,Q_Qmax[2],W,Q_F1,W_F1,F1_F1,num_str);//,Q_F2,W_F2,F2_F2,num_str2);
			 double CS=part1*part2;
			if (real_W<W){
				double an_dw=1/(W-porog_ch(type_chanel));
				double CS2=CS*an_dw*real_W-CS*porog_ch(type_chanel)*an_dw;
				CS=CS2;
			}
			if (real_W>W){CS=CS*get_CS_ph_int_Ev(real_W)/get_CS_ph_int_Ev(W);}
			if ((type_chanel==1)&&(real_W<1.635)) return 1.8*CS;
			if ((type_chanel==1)&&(real_W>1.635)&&(real_W<1.66)) return 1.2*CS;
			if ((type_chanel==1)&&(real_W>2.35)&&(real_W<2.46)) return 1.08*CS;
			if ((type_chanel==1)&&(real_W>2.35)) return 1.2*CS;
			
			if ((type_chanel==2)&&(real_W>2.20)&&(real_W<2.3)) return 1.07*CS;
			if ((type_chanel==2)&&(real_W>2.3)) return 1.14*CS;
			if ((type_chanel==2)&&(real_W<1.748)) return 1.12*CS;
			if ((type_chanel==2)&&(real_W>1.748)&&(real_W<1.752)) return 1.08*CS;
			if ((type_chanel==2)&&(real_W>1.752)&&(real_W<1.775)) return 1.02*CS;
			return CS;
			}
			return 0;
		}

		double Sigma::get_CS_int_fi(double Q, double W, double costeta, double Ebeam)
		{
			bool ch=check_kin(Q,W,Ebeam);
			if (ch==0) 
			{
				cout<<"uncorrect input (kinematic)"<<endl;
				return 0;
			}

			double l_p=0;
			double r_p=0;
			double s_p=0;
			double step=10;

			if (range_fi==1)
			{
			l_p=-180;
			r_p=180;
			s_p=l_p+step/2;
			}
			if (range_fi==2)
			{
			l_p=0;
			r_p=360;
			s_p=l_p+step/2;
			}

			double bin=0;
			double CS_int_fi=0; 

			for (double i=s_p;i<r_p;i+=step)
			{
				if (abs(i-s_p)<0.1) {bin=abs(abs(l_p)-abs(i))+(step/2);}
				else 
				{
					if ((i+step)>r_p) {bin=(step/2)+abs(r_p-i);}
					else bin=step;
				}
				CS_int_fi+=bin*get_CS_f_i(Q,W,costeta,i,Ebeam)/57.29578049;
			}
			return CS_int_fi;
		}

		double Sigma::get_CS_int_fi_cos(double Q, double W, double Ebeam)
		{
		
			if (Q<4.16) return getCS_d3CS(Q,W,Ebeam);
			bool ch=check_kin(Q,W,Ebeam);
			if (ch==0){
				cout<<"uncorrect input"<<endl;
				return 0;
			}

			//cout<<" 111"<<endl;
			double CS_int_fi=0;	
			double CS_int_fi_costeta=0;
			int edge2=0;
			double bin2=0;
			double cache2=0;

			if ((range_fi!=2)||(range_cos==1))
			{
			//cout<<" num_costeta="<<num_costeta<<endl;
				for(int k=0;k<num_costeta;k++)
				{
					edge2=num_costeta-1;
					CS_int_fi=get_CS_int_fi(Q,W,costeta_Qmax_int[k],Ebeam);	

					if ((range_fi==1)||(range_cos==1))
					{
					if (k==0){bin2=((1+costeta_Qmax_int[k])-(costeta_Qmax_int[k]-costeta_Qmax_int[k+1])/2);}
					if (k==edge2){bin2=((1-costeta_Qmax_int[k])+(costeta_Qmax_int[k]-costeta_Qmax_int[k-1])/2);}
					if ((k!=0)&&(k!=edge2))
					{bin2=(abs(costeta_Qmax_int[k]-costeta_Qmax_int[k-1])/2
							+abs(costeta_Qmax_int[k]-costeta_Qmax_int[k+1])/2);}	
					}

					cache2=CS_int_fi*bin2;
					CS_int_fi_costeta=CS_int_fi_costeta+cache2;

					cache2=0;
					bin2=0;
					CS_int_fi=0;
				}
			}
			//cout<<" 112"<<endl;
			if ((range_fi==2)&&(range_cos==2))
			{
				int start_cicle=0;
				if (W<1.47){start_cicle=1;}
				else {start_cicle=0;}

				double cos_for_int[10]={-0.9,-0.7,-0.5,-0.3,-0.1,0.1,0.3,0.5,0.7,0.9};
				num_costeta=10;

				for(int k=start_cicle;k<num_costeta;k++)
				{
					edge2=num_costeta-1;
					CS_int_fi=get_CS_int_fi(Q,W,cos_for_int[k],Ebeam);	

					if (k==start_cicle){bin2=((1+cos_for_int[k])-(cos_for_int[k]-cos_for_int[k+1])/2);}
					if (k==edge2){bin2=((1-cos_for_int[k])+(cos_for_int[k]-cos_for_int[k-1])/2);}
					if ((k!=start_cicle)&&(k!=edge2))
					{bin2=(abs(cos_for_int[k]-cos_for_int[k-1])/2
							+abs(cos_for_int[k]-cos_for_int[k+1])/2);}	

					cache2=CS_int_fi*bin2;
					CS_int_fi_costeta=CS_int_fi_costeta+cache2;

					cache2=0;
					bin2=0;
					CS_int_fi=0;
				}
			}
			//cout<<" CS_int_fi_costeta: "<<CS_int_fi_costeta<<endl;
			return CS_int_fi_costeta;
		}

		Sigma::Sigma(string dataPath, int chanel)
		{
			ifstream interp_right;
			ifstream file_F1;
			ifstream file_F2;
			ifstream file_Qmax;
			ifstream CS_data;
			ifstream CS_data_int;
			ifstream file_low_photo_data;
			ifstream file_ph_int;
			ifstream from_Evgen;
			ifstream data_error;

			type_chanel=chanel;
//cout<<"000000000000000000000";
			switch (chanel)
			{
				case 1:{
					Q_max_channel=3.44999;
					max_cos=0.949999;
					min_cos=-0.7749999;
					range_fi=1;range_cos=1; 

					interp_right.open(dataPath+"/KL_interp.txt");
					if (!interp_right.is_open()) cout<<"ERROR: wrong way to the data, pls check eg_config_test.txt";
					file_F1.open(dataPath+"/KLambda_Fit_F1.txt");
					file_F2.open(dataPath+"/KLambda_Fit_F2.txt");
					file_Qmax.open(dataPath+"/KLambda_CS_Qmax_Fit_gladk.txt");

					CS_data.open(dataPath+"/KL_s4.txt");
					CS_data_int.open(dataPath+"/CS_PiN_int.txt");//ISPRAVIT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
					file_low_photo_data.open(dataPath+"/KL_ph_out.txt");
					file_ph_int.open(dataPath+"/KL_ph_int.txt");
					from_Evgen.open(dataPath+"/KL_ph_int.txt");

				}break;
				case 2:{ 
					Q_max_channel=3.4499;
					max_cos=0.949999;
					min_cos=-0.7749999;
					range_fi=1;range_cos=1; 

					interp_right.open(dataPath+"/KS_interp.txt");
					if (!interp_right.is_open()) cout<<"ERROR: wrong way to the data, pls check eg_config_test.txt";
					file_F1.open(dataPath+"/KSigma_Fit_F1.txt");
					file_F2.open(dataPath+"/KSigma_Fit_F2.txt");
					file_Qmax.open(dataPath+"/KSigma_CS_Qmax_Fit_gladk.txt");

					CS_data.open(dataPath+"/KS_s4.txt");
					CS_data_int.open(dataPath+"/CS_PiN_int.txt");//ISPRAVIT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
					file_low_photo_data.open(dataPath+"/KS_ph_out_new.txt");
					file_ph_int.open(dataPath+"/KS_ph_int_new.txt");
					from_Evgen.open(dataPath+"/KS_ph_int_new.txt");
					data_error.open(dataPath+"/CSKS_Theory.txt");
				}break;
				case 3:{ 
					Q_max_channel=6;
   					range_fi=2;range_cos=1; 
					max_cos=0.89999;
					min_cos=-0.89999;
					interp_right.open(dataPath+"/Pi0P_interp.txt");
					if (!interp_right.is_open()) cout<<"ERROR: wrong way to the data, pls check eg_config_test.txt";
					file_F1.open(dataPath+"/Pi0P_Fit_F1.txt");
					file_F2.open(dataPath+"/Pi0P_Fit_F2.txt");
					file_Qmax.open(dataPath+"/Pi0P_CS_Qmax_Fit_gladk.txt");

					CS_data.open(dataPath+"/CS_pi0p.txt");
					CS_data_int.open(dataPath+"/CS_PiN_int.txt");//ISPRAVIT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
					file_low_photo_data.open(dataPath+"/pi0p_ph_out.txt");
					file_ph_int.open(dataPath+"/pi0p_ph_int.txt");
					from_Evgen.open(dataPath+"/pi0p.data");
				}break;
				case 4:{ 
					Q_max_channel=4.16;
					range_fi=2;range_cos=2; 
					max_cos=0.899999;
					min_cos=-0.699999;
					interp_right.open(dataPath+"/PiN_interp.txt");
					if (!interp_right.is_open()) cout<<"ERROR: wrong way to the data, pls check eg_config_test.txt";
					file_F1.open(dataPath+"/PiN_Fit_F1.txt");
					file_F2.open(dataPath+"/PiN_Fit_F2.txt");
					file_Qmax.open(dataPath+"/PiN_CS_Qmax_Fit_gladk.txt");

					CS_data.open(dataPath+"/CS_PiN.txt");
					CS_data_int.open(dataPath+"/CS_PiN_int.txt");
					file_low_photo_data.open(dataPath+"/piN_ph_out.txt");
					file_ph_int.open(dataPath+"/piN_ph_int_out.txt");
					from_Evgen.open(dataPath+"/pipn.data");

				}break;
				defult:{ 
					cout<<" error in consturctor Sigma"<<endl;
				}
			}

			char name_F[8];//skip first column
			double chislo;//buffer

			double Current_Q;
			double Current_W;
			double Current_F;

			string str01;
			getline(file_F1, str01);
			while(!file_F1.eof())
			{
				num_str++;
				file_F1>>name_F;//F1

				file_F1>>chislo;//Q
				if(chislo!=0)
				{
					Current_Q=chislo;
					Q_F1.push_back(chislo);
				}

				file_F1>>chislo;//W
				if(chislo!=0)
				{
					Current_W=chislo;
					if (chislo>max_W_F1) {max_W_F1=chislo;}
					if (chislo<min_W_F1) {min_W_F1=chislo;}
					W_F1.push_back(chislo);
				}

				file_F1>>chislo;//F1
				if(chislo!=0){
					Current_F=chislo;
					F1_F1.push_back(chislo);
				}
				file_F1>>chislo;//error of F1

			}
			num_str--;

			/*
			Current_Q=0;
			Current_W=0;
			Current_F=0;

			chislo=0;

			getline(file_F2, str01);
			while(!file_F2.eof())
			{
				num_str2++;
				file_F2>>name_F;//F2
				//cout<<name_F<<endl;

				file_F2>>chislo;//Q
				if(chislo!=0)
				{
					Current_Q=chislo;
					//cout<<"Current_Q_F2="<<Current_Q<<endl;
					Q_F2.push_back(chislo);
				}

				file_F2>>chislo;//W
				if(chislo!=0)
				{
					Current_W=chislo;
					if (chislo>max_W_F2) {max_W_F2=chislo;}
					if (chislo<min_W_F2) {min_W_F2=chislo;}
					//cout<<"Current_W_F2="<<Current_W<<endl;
					W_F2.push_back(chislo);
				}

				file_F2>>chislo;//F2
				if(chislo!=0)
				{
					Current_F=chislo;
					//cout<<"Current_F2_F2="<<Current_F<<endl;
					F2_F2.push_back(chislo);
				}

				file_F2>>chislo;//error of F2

			}

			num_str2--;

			*/
			chislo=-1;

			Current_Q=0;
			Current_W=0;
			Current_F=0;
			double Current_costeta_Qmax=0;
			double Current_p0=0;
			double Current_p1=0;
			double Current_p2=0;
			double Current_param=0;
			double Current_glad=0,Current_costeta=0;


			getline(file_Qmax, str01);
			getline(file_Qmax, str01);

			vector<double> Ebeam_Qmax;
			double Current_Ebeam_Qmax=0;

			while(!file_Qmax.eof()){
				num_str3++;
				file_Qmax>>chislo;//W
				if(chislo!=-1)
				{
					Current_W=chislo;
					if (chislo>max_W_Qmax) {max_W_Qmax=chislo;}
					if (chislo<min_W_Qmax) {min_W_Qmax=chislo;}
					W_Qmax.push_back(chislo);
				}

				file_Qmax>>chislo;//Q
				if(chislo!=-1)
				{
					Current_Q=chislo;
					Q_Qmax.push_back(chislo);
				}
		
				file_Qmax>>chislo;//Ebeam
				if(chislo!=-1)
				{
					Current_Ebeam_Qmax=chislo;
					Ebeam_Qmax.push_back(chislo);
				}

				file_Qmax>>chislo;//cos(teta)
				if(chislo!=-1)
				{
					Current_costeta_Qmax=chislo;
					costeta_Qmax.push_back(chislo);
				}

				file_Qmax>>chislo;//p0
				if(chislo!=-1){Current_p0=chislo;}

				file_Qmax>>chislo;//p1
				if(chislo!=-1){Current_p1=chislo;}

				file_Qmax>>chislo;//p2
				if(chislo!=-1){Current_p2=chislo;}

				file_Qmax>>chislo;//param
				if(chislo!=-1)
				{
					//cout<<"pram: "<<chislo<<endl;
					param.push_back(chislo);
					Current_param=chislo;
				}


				file_Qmax>>chislo;//gladk
				if(chislo!=-1)
				{
					Current_glad=chislo;
					//cout<<"Current_p2="<<Current_p2<<endl;

					p2_Qmax.push_back(Current_p2*Current_param*Current_glad);
					p1_Qmax.push_back(Current_p1*Current_param*Current_glad);
					p0_Qmax.push_back(Current_p0*Current_param*Current_glad);
				}

				chislo=-1;
		
			}

			num_str3--;

			chislo=-1;
			getline(interp_right, str01);
			getline(interp_right, str01);

			while(!interp_right.eof())
			{
				num_str4++;
				bool tmp_W=0;
				interp_right>>chislo;//W
				if(chislo!=-1)
				{
					Current_W=chislo;
					if (chislo>max_W_inter) {max_W_inter=chislo;}
					if (chislo<min_W_inter) {min_W_inter=chislo;}
					//cout<<"Current_W(inter right)="<<Current_W<<endl;
					if (((chislo>2.39)&&(chanel==2))||((chislo>2.47)&&(chanel==1))||((chislo>1.479)&&(chanel==3))) tmp_W=1;
					else tmp_W=0;
					W_inter.push_back(chislo);
				}
				interp_right>>chislo;//CS interp right
				if(chislo!=-1)
				{
					//cout<<"CS_inter="<<chislo<<endl;
					if (tmp_W==1){
						if (chanel==2) CS_inter.push_back(chislo*1.08);
						if (chanel==1) CS_inter.push_back(chislo*1.13);
						if ((chanel==3)&&(Current_W>2.131)) CS_inter.push_back(chislo*2.2);
						if ((chanel==3)&&(Current_W<1.51)) CS_inter.push_back(chislo*1.15);
						if ((chanel==3)&&(Current_W>1.51)&&(Current_W<1.60)) CS_inter.push_back(chislo*1.3);
						if ((chanel==3)&&(Current_W>1.65)&&(Current_W<1.8)) CS_inter.push_back(chislo*1.6);
						if ((chanel==3)&&(Current_W>1.85)&&(Current_W<2.31)) CS_inter.push_back(chislo*2);
					}else CS_inter.push_back(chislo);
				}

				chislo=-1;
			}
			num_str4--;


			if (chanel==1||chanel==2||chanel==3) {

				for(int i=0;(costeta_Qmax[i+1]>=costeta_Qmax[i]);i++)
				{
					if (costeta_Qmax[i+1]>costeta_Qmax[i])
					{costeta_Qmax_int.push_back(costeta_Qmax[i]);}
	
					if (!(costeta_Qmax[i+1+1]>=costeta_Qmax[i+1]))
					{
						costeta_Qmax_int.push_back(costeta_Qmax[i+1]);
					}
				}

				num_costeta=costeta_Qmax_int.size();
			}



			while(!CS_data.eof())
			{
				n_str_CS++;

				CS_data>>chislo;//W
				if(chislo!=-3)
				{
					Current_W=chislo;
					//cout<<" W: "<<Current_W;
					if (chislo>max_W) {max_W=chislo;}
					if (chislo<min_W) {min_W=chislo;}
					_W.push_back(chislo);
				}

				CS_data>>chislo;//Q
				//cout<<" "<<chislo;
				if (chislo==0) break;
				if(chislo!=-3)
				{
					if (Qmin>chislo) {Qmin=chislo;}
					if (Qmax<chislo) {Qmax=chislo;}
					Current_Q=chislo;
					//cout<<" Q2: "<<Current_Q;
					_Q2.push_back(chislo);
				}
				if(chanel==3) CS_data>>chislo;//E

				CS_data>>chislo;//cos(teta)
				//cout<<" cos: "<<chislo;
				if(chislo!=-3){
					Current_costeta=chislo;
					//cout<<" cos: "<<Current_costeta;
					_cos.push_back(Current_costeta);
				}

				CS_data>>chislo;//p0
				//cout<<" p0: "<<chislo;
				if(chislo!=-3)
				{
					Current_p0=chislo;
					if(chislo<0) {cout<<"ERROR <0, CS: "<<chislo<<endl;}
				}

				CS_data>>chislo;//p1
				//cout<<" p1: "<<chislo;
				if(chislo!=-3){Current_p1=chislo;}

				CS_data>>chislo;//p2
				//cout<<" p2: "<<chislo;
				if(chislo!=-3){Current_p2=chislo;}

				CS_data>>chislo;//param
				//cout<<" par: "<<chislo;
				if(chislo!=-3)
					{param_vec.push_back(chislo);
					Current_param=chislo;
				}

				if (is_there_glad==1){
					CS_data>>chislo;//gladk
					if(chislo!=-3)
					{
						Current_glad=chislo;
						//cout<<"Current_p2="<<Current_p2<<endl;
						_p2.push_back(Current_p2*Current_param*Current_glad);
						_p1.push_back(Current_p1*Current_param*Current_glad);
						_p0.push_back(Current_p0*Current_param*Current_glad);
					}
				}else{
						_p2.push_back(Current_p2*Current_param);
						_p1.push_back(Current_p1*Current_param);
						_p0.push_back(Current_p0*Current_param);
				}

				chislo=-3;
		
			}

			n_str_CS--;

			chislo=-1;

			search_exterm_points();
			num_ext_p=_Q_for_ext_point.size();

			//getline(file_low_photo_data, str01);
			//getline(file_low_photo_data, str01);
			chislo=-1;

			while(!file_low_photo_data.eof())
			{

				n_str_ph++;

				file_low_photo_data>>chislo;//W
				if(chislo!=-1)
				{
					Current_W=chislo;
					if (chislo>max_W_ph) {max_W_ph=chislo;}
					if (chislo<min_W_ph) {min_W_ph=chislo;}
					W_vec_ph.push_back(chislo);
				}

				file_low_photo_data>>chislo;//cos(teta)
				if(chislo!=-2)
				{
					Current_costeta=chislo;
					costeta_vec_ph.push_back(chislo);
				}

				file_low_photo_data>>chislo;//p0
				if(chislo!=-1)
				{
					Current_p0=chislo;
					CS_vec_ph.push_back(chislo);
				}
				//file_low_photo_data>>chislo;//error

				chislo=-1;
		
			}
			n_str_ph--;

			chislo=-1;

			while(!CS_data_int.eof()){
				n_str_CS_int++;

				CS_data_int>>chislo;//Q
				if(chislo!=-1){
					Current_W=chislo;
					_Q_int.push_back(chislo);
				}

				CS_data_int>>chislo;//W
				if(chislo!=-2){
					Current_costeta=chislo;
					_W_int.push_back(chislo);
				}

				CS_data_int>>chislo;//CS
				if(chislo!=-1)
				{
					Current_p0=chislo;
					_CS_int.push_back(chislo);

				}
				//file_low_photo_data>>chislo;//error
		
				chislo=-1;
		
			}
			n_str_CS_int--;

			//getline(file_ph_int, str01);
			//getline(file_ph_int, str01);
			chislo=-1;

			while(!file_ph_int.eof())
			{
				n_str_ph_int++;
				file_ph_int>>chislo;//W
				if(chislo!=-1){
					Current_W=chislo;
					if (chislo>max_W_ph) {max_W_ph=chislo;}
					if (chislo<min_W_ph) {min_W_ph=chislo;}
					_W_ph.push_back(chislo);
				}
				file_ph_int>>chislo;//p0
				if(chislo!=-1){
					Current_p0=chislo;
					_CS_ph.push_back(chislo);
				}
				chislo=-1;
			}
			n_str_ph_int--;

			chislo=-1;

			while(!from_Evgen.eof())
			{
				n_str_Ev++;
				from_Evgen>>chislo;//W
				if(chislo!=-1){
					Current_W=chislo;
					_W_ph_Ev.push_back(chislo);
				}
				from_Evgen>>chislo;//p0
				if(chislo!=-1){
					Current_p0=chislo;
					_CS_ph_Ev.push_back(chislo);
				}
				chislo=-1;
			}
			n_str_Ev--;
			
			//added 04/21/2022
			
			channel channelName = static_cast<channel>(chanel - 1);
			read_maxQ2(dataPath, channelName);
			read_StrFunF1(dataPath, channelName);
			////
		
		}
		double Sigma::check_cos(double costeta,double W){
			if ((costeta>max_cos)&&(costeta<=1.00001)){return max_cos;}
			if ((costeta<min_cos)&&(costeta>=-1.00001)){return min_cos;}
			return costeta;
		}

double Sigma::porog_ch(int num_chanel){//1-KL 2-KS 3-Pi0P 4-PiN
	switch (num_chanel)
	{
		case 1: return massLambda+massKaon; break;
		case 2: return massSigma0+massKaon; break;
		case 3: return massPion0+massProton; break;
		case 4: return massPion+massNeutron; break;
		defult: cout<<" uncorrect param for calc porog energy"<<endl;
	}
	return 0;
}

double  Sigma::getK(double Q2, double W){
	return (2*getomega(Q2,W)*massProton-Q2)/(2*massProton);
}

double Sigma::fun_points(double Q, double W, const vector<double>& Q_F1, const vector<double>& W_F1, const vector<double>& F1_F1,int num_str){
	double start_point_Q=-1;
	int num_interp=2,only_one=0;

	for (int i=1;i<num_str;i++)
	{
		if ((Q_F1[i]>Q)&&(Q_F1[i-1]<Q)) {start_point_Q=i;}
		if ((Q_F1[i]==Q)&&(only_one==0))
		{
			if (i==1) {start_point_Q=0;}
			if (i!=1) {start_point_Q=i;}
			num_interp=1; 
			only_one=1;
		}
	}
	if (start_point_Q==-1){cout<<"Error: we cant interpolate with such input parameters1"<<endl;}

	int num_W_l=0;
	for(int i=start_point_Q-1;Q_F1[i]==Q_F1[start_point_Q-1];i--)
	{
		num_W_l++;
	}

	int num_W_h=0;
	for(int i=start_point_Q;Q_F1[i]==Q_F1[start_point_Q];i++)
	{
		num_W_h++;
	}

	double setka_W_h=0;
	double setka_W_l=0;
	double F1_W_h_Q_h=0;//znachenie F1 v tochke W,Q high
	double F1_W_l_Q_h=0;
	double F1_W_h_Q_l=0;
	double F1_W_l_Q_l=0;
	
	if (num_interp!=1)
	{
		for(int i=start_point_Q-1-num_W_l;i<start_point_Q;i++)
		{
			if ((W_F1[i]>W)&&(W_F1[i-1]<W)){F1_W_h_Q_l=F1_F1[i]; F1_W_l_Q_l=F1_F1[i-1];}
			if (W_F1[i]==W)
			{
				num_interp=0;
				F1_W_h_Q_l=F1_F1[i];
				F1_W_l_Q_l=F1_F1[i];
			}
		}
	}
	for (int i=start_point_Q;i<num_W_h+start_point_Q;i++){
		if ((W_F1[i]>W)&&(W_F1[i-1]<W)) {setka_W_h=W_F1[i]; setka_W_l=W_F1[i-1]; F1_W_h_Q_h=F1_F1[i]; F1_W_l_Q_h=F1_F1[i-1]; break;}
		if (W_F1[i]==W)
		{
			setka_W_h=W_F1[i]; 
			setka_W_l=W_F1[i]; 
			F1_W_h_Q_h=F1_F1[i]; 
			F1_W_l_Q_h=F1_F1[i];
			if (num_interp==1) {return F1_F1[i];}
		}
	}


	double points[10];

	if (start_point_Q-1<0){start_point_Q=1;}
	points[0]=(Q_F1[start_point_Q-1]);
	points[1]=(Q);
	points[2]=(Q_F1[start_point_Q]);
	points[3]=(setka_W_l);//3
	points[4]=(W);
	points[5]=(setka_W_h);//5
	points[6]=(F1_W_l_Q_l);
	points[7]=(F1_W_h_Q_l);//7
	points[8]=(F1_W_l_Q_h);
	points[9]=(F1_W_h_Q_h);//9

	double fR1=0;
	double fR2=0;
	//y-W x-Q
	if (num_interp==0){return points[6]+(points[8]-points[6])*(points[1]-points[0])/(points[2]-points[0]);}
	if (num_interp==1){return points[8]+(points[9]-points[8])*(points[4]-points[3])/(points[5]-points[3]);}
	if (num_interp==2)
	{
		fR1=points[6]*(points[2]-points[1])/(points[2]-points[0])+points[8]*(points[1]-points[0])/(points[2]-points[0]);
		fR2=points[7]*(points[2]-points[1])/(points[2]-points[0])+points[9]*(points[1]-points[0])/(points[2]-points[0]);
		return ( fR1*(points[5]-points[4])/(points[5]-points[3])+fR2*(points[4]-points[3])/(points[5]-points[3]) );
	}
	return 0;

}
double Sigma::lineal_interp_1(double x_r, const vector<double>& x,const vector<double>& y,int num_str){
	int start_point_x=-1;

	for (int i=0;i<num_str;i++)
	{
		if ((x[i]<x_r)&&(x[i+1]>x_r)&&(i+1,num_str))
		{
			start_point_x=i;
			return (y[i]+(y[i+1]-y[i])*(x_r-x[i])/(x[i+1]-x[i]));
		}
		if (abs(x[i]-x_r)<0.0001) {start_point_x=i; return y[i];}
	}
	
	return 0;
}
/*
double Sigma::interpol(double Q, double W, double fi, const vector<double>& Q_F1, const vector<double>& W_F1,
			 const vector<double>& p0_Qmax, const vector<double>& p1_Qmax, const vector<double>& p2_Qmax,double num_str)
	{
		double start_point_Q=-1;

		int num_interp=2;
		int only_one=0;

		for (int i=1;i<num_str;i++)
		{
			if ((Q_F1[i]>Q)&&(Q_F1[i-1]<Q)) {start_point_Q=i;}
			if ((Q_F1[i]==Q)&&(only_one==0))
			{
				if (i==1) {start_point_Q=0;}
				if (i!=1) {start_point_Q=i;}
				num_interp=1; 
				only_one=1;
			}
		}

		if (start_point_Q==-1){cout<<"Error: we cant interpolate with such input parameters2"<<endl;}

		int num_W_l=0;
		for(int i=start_point_Q-1;Q_F1[i]==Q_F1[start_point_Q-1];i--)
		{
			num_W_l++;
		}


		int num_W_h=0;
		for(int i=start_point_Q;Q_F1[i]==Q_F1[start_point_Q];i++)
		{
			num_W_h++;
		}

		double setka_W_h=0;
		double setka_W_l=0;
		double F1_W_h_Q_h=0;//znachenie F1 v tochke W,Q high
		double F1_W_l_Q_h=0;
		double F1_W_h_Q_l=0;
		double F1_W_l_Q_l=0;
			if((abs(Q-4.3)<0.01)&&(W<1.5)) cout<<" 2! "<<endl;
		if (num_interp!=1)
		{
			for(int i=start_point_Q-1-num_W_l;i<start_point_Q;i++)
			{
				if ((W_F1[i]>W)&&(W_F1[i-1]<W))
				{
					F1_W_h_Q_l=anti_Fit(fi,p0_Qmax[i],p1_Qmax[i],p2_Qmax[i]); 
					F1_W_l_Q_l=anti_Fit(fi,p0_Qmax[i-1],p1_Qmax[i-1],p2_Qmax[i-1]);
				}
				if (W_F1[i]==W)
				{
					num_interp=0;
					F1_W_h_Q_l=anti_Fit(fi,p0_Qmax[i],p1_Qmax[i],p2_Qmax[i]); 
					F1_W_l_Q_l=anti_Fit(fi,p0_Qmax[i],p1_Qmax[i],p2_Qmax[i]); 
				}
			}
		}
	
		if((abs(Q-4.3)<0.01)&&(W<1.5)) cout<<" 3! "<<endl;
		for (int i=start_point_Q;i<num_W_h+start_point_Q;i++)
		{
			if ((W_F1[i]>W)&&(W_F1[i-1]<W)) 
			{
				setka_W_h=W_F1[i]; setka_W_l=W_F1[i-1]; 
				F1_W_h_Q_h=anti_Fit(fi,p0_Qmax[i],p1_Qmax[i],p2_Qmax[i]);  
				F1_W_l_Q_h=anti_Fit(fi,p0_Qmax[i-1],p1_Qmax[i-1],p2_Qmax[i-1]);
			}
	
			if (W_F1[i]==W)
			{
				setka_W_h=W_F1[i]; 
				setka_W_l=W_F1[i]; 
				F1_W_h_Q_h=anti_Fit(fi,p0_Qmax[i],p1_Qmax[i],p2_Qmax[i]);  
				F1_W_l_Q_h=anti_Fit(fi,p0_Qmax[i],p1_Qmax[i],p2_Qmax[i]); 
				if (num_interp==1) {return anti_Fit(fi,p0_Qmax[i],p1_Qmax[i],p2_Qmax[i]); ;}
			}

		}
		if((abs(Q-4.3)<0.01)&&(W<1.5)) cout<<" 4! "<<endl;
		double points[10];

		if (start_point_Q-1<0){start_point_Q=1;}
		points[0]=(Q_F1[start_point_Q-1]);//0
		points[1]=(Q);
		points[2]=(Q_F1[start_point_Q]);
		points[3]=(setka_W_l);//3
		points[4]=(W);
		points[5]=(setka_W_h);//5
		points[6]=(F1_W_l_Q_l);
		points[7]=(F1_W_h_Q_l);//7
		points[8]=(F1_W_l_Q_h);
		points[9]=(F1_W_h_Q_h);//9

		double fR1=0;
		double fR2=0;
		//y-W x-Q
		if (num_interp==0){return points[6]+(points[8]-points[6])*(points[1]-points[0])/(points[2]-points[0]);}
		if (num_interp==1){return points[8]+(points[9]-points[8])*(points[4]-points[3])/(points[5]-points[3]);}
		if (num_interp==2)
		{
			fR1=points[6]*(points[2]-points[1])/(points[2]-points[0])+
						points[8]*(points[1]-points[0])/(points[2]-points[0]);
			fR2=points[7]*(points[2]-points[1])/(points[2]-points[0])+
						points[9]*(points[1]-points[0])/(points[2]-points[0]);
			return ( fR1*(points[5]-points[4])/(points[5]-points[3])+fR2*(points[4]-points[3])/(points[5]-points[3]) );
		}
		cout<<" I dont know why"<<endl;
		return 0;

	}
*/
/*
double Sigma::getCS_fit(double Ebeam, double Q, double Qmax,double W, 
			const vector<double>& Q_F1, const vector<double>& W_F1, const vector<double>& F1_F1,double num_str)
			//,vector<double> Q_F2, vector<double> W_F2, vector<double> F2_F2,double num_str2)
	{
		double Qmax2=Qmax;//*Qmax;
		//cout<<"Q="<<Q<<" Qmax: "<<Qmax<<" W: "<<W<<endl;
		double interpol_F1_max=fun_points(Qmax,W,Q_F1,W_F1,F1_F1,num_str);//calc of interpol of F1;

		//double interpol_F2_max=fun_points(Qmax,W,Q_F2,W_F2,F2_F2,num_str2);//calc of interpol of F2;

		double Gt_max=interpol_F1_max*4*constantPi2*constantAlpha/(getK(Qmax,W)*massProton);

		//double Gl_max=(4*constantPi2*constantAlpha*interpol_F2_max*2*massProton*(-Qmax-getomega(Qmax2,W)*getomega(Qmax2,W))/
		//					(getomega(Qmax2,W)*(-Qmax)*(2*getomega(Qmax2,W)*massProton-Qmax)))-Gt_max;
	
	
		double Gl_max=Gt_max*0.2;

		double G_max=Gt_max+getEpsilon(Ebeam,Qmax2,W)*Gl_max;//?


		double Q2=Q;//*Q;
		double interpol_F1=fun_points(Q,W,Q_F1,W_F1,F1_F1,num_str);//calc of interpol of F1;

		//double interpol_F2=fun_points(Q,W,Q_F2,W_F2,F2_F2,num_str2);//calc of interpol of F2;

		double Gt=interpol_F1*4*constantPi2*constantAlpha/(getK(Q,W)*massProton);
		//double Gl=(4*constantPi2*constantAlpha*interpol_F2*2*massProton*(-Q-getomega(Q2,W)*getomega	(Q2,W))/
		//				(getomega(Q2,W)*(-Q)*(2*getomega(Q2,W)*massProton-Q)))-Gt;
		

		double Gl=Gt*0.2;	

		double G=Gt+getEpsilon(Ebeam,Q2,W)*Gl;

		return G/G_max;
	}
*/
		double Sigma::check_input_param(double Q,double W, double Ebeam){
			double a=W;
			if ((Q<2)||(W<porog_ch(type_chanel))||(W>5.000000001)){return 0;}
			if ((W>max_W_Qmax)||(W>max_W_F1)){
				if (max_W_Qmax>=max_W_F1){a=max_W_F1;}
				else {a=max_W_Qmax;}
			}
			if ((W<min_W_Qmax)||(W<min_W_F1)){
				if (min_W_Qmax>=min_W_F1){a=min_W_Qmax;}
				else {a=min_W_F1;}
			}
			return a;
		}

double Sigma::dsigma_dcos(double _beam_energy, double _Q2, double _W, double costeta){
 double cache1=0;
 for (double fi=0.05;fi<=6.28;fi+=0.1) {cache1+=0.1*d5sigma2(_beam_energy,_Q2,_W,costeta,fi)/getGamma(_beam_energy, _Q2, _W);}
 //cout<<"Q: "<<_Q2<<" W: "<<_W<<" cos: "<<costeta<<"D_cos: "<<cache1<<endl;
 return cache1;
}

void Sigma::search_exterm_points(){
	double W_st1=1000,W_st2=0,Q2_st=_Q2[0];
	for (int i=0;i<n_str_CS;i++){
		if (_Q2[i]==Q2_st){
			if (_W[i]<W_st1) {W_st1=_W[i];}
			if (_W[i]>W_st2) {W_st2=_W[i];}
		}else{
			_W_max.push_back(W_st2);
			_W_min.push_back(W_st1);
			_Q_for_ext_point.push_back(_Q2[i-1]);
			W_st1=1000;
			W_st2=0;
		}
		
		Q2_st=_Q2[i];	
	}
	_W_max.push_back(W_st2);
	_W_min.push_back(W_st1);
	_Q_for_ext_point.push_back(Q2_st);

	_W_min_all=_W_min[0];
	for (int i=0;i<_W_min.size();i++){
		if (_W_min[i]>_W_min_all) {_W_min_all=_W_min[i];}
		if (_W_min[i]<W_ext_min) {W_ext_min=_W_min[i];}
		if (_W_max[i]>W_ext_max) {W_ext_max=_W_max[i];}
		//if (_W_max[i]<_W_max[0])
	}
	if(type_chanel==4){_W_max_all=1.51;}
	if(type_chanel==1){_W_max_all=2.05;}
	if(type_chanel==2){_W_max_all=2.05;}
	if(type_chanel==3){_W_min_all=1.13;}

}

int Sigma::check_possibil_inter_Q2(double Q){
	if (Q<Qmin) return 1;
	if (Q>Qmax) return 2;
	return 0;
}
double Sigma::change_Q2(double Q){
	if (Q<Qmin) return Qmin;
	if (Q>Qmax) cout<<" ERROR it large Q2"<<endl;
	return Q;
}

int Sigma::check_possibil_inter_W(double Q, double W){

	for (int i=0;i<num_ext_p;i++){
		if(Q==_Q_for_ext_point[i]){
			interp_num1=i;
			interp_num2=i;
			if (W>W_ext_max) return 8;// all data have alredy used
			if (W<W_ext_min) return 7;// all data have alredy used
			if ((W<=W_ext_max)&&(W>_W_max_all)) return 6;// one point exterp
			if ((W>=W_ext_min)&&(W<_W_min_all)) return 5;// one point exterp
			//if ((W>_W_max[i])||(W>_W_max_all)) return 4;
			//if ((W<_W_min[i])||(W<_W_min_all)) return 3;
			if (W>=_W_min_all) return 0;
			if (W<=_W_max_all) return 0;
		}
		if((Q>_Q_for_ext_point[i])&&(Q<_Q_for_ext_point[i+1])){
			interp_num1=i;
			interp_num2=i+1;
			if (W>W_ext_max) return 8;// all data have alredy used
			if (W<W_ext_min) return 7;// all data have alredy used
			if ((W<=W_ext_max)&&(W>_W_max_all)) return 6;// one point exterp
			if ((W>=W_ext_min)&&(W<_W_min_all)) return 5;// one point exterp
			//if ((W>_W_max[i])||(W>_W_max[i+1])||(W>_W_max_all)) return 4;
			//if ((W<_W_min[i])||(W<_W_min[i+1])||(W<_W_min_all)) return 3;
			if ((W<=_W_max[i])&&(W<=_W_max[i+1])&&(W<=_W_max_all)) return 0;
			if ((W>=_W_min[i])&&(W>=_W_min[i+1])&&(W>=_W_min_all)) return 0;
		}
	}
cout<<" ERROR check_possibil_inter_W"<<endl;
return -1;
}

int Sigma::range_ph(double W){
	//cout<<" min_W_ph: "<<min_W_ph<<" max_W_ph: "<<max_W_ph<<endl;
	if (W<min_W_ph) return 1;
	if (W>max_W_ph) return 2;
	return 0;
}
double Sigma::int_get_d5CS(double Q,double W, double Ebeam){
 double cache2=0;
 for (double cos=-1;cos<=1;cos+=0.02){
  //cout<<" cos: "<<cos<<endl;
  cache2+=0.02*dsigma_dcos(Ebeam,Q,W,cos);
 } 
 //cout<<cache2<<endl;
	return cache2;
}
double Sigma::get_d5CS(double Q,double W, double cos, double fi, double Ebeam){
 if (range_fi==1){
  fi=fi-180;
  if (abs(fi+180)<0.001) {fi=-179.9999;}
  if (abs(fi-180)<0.001) {fi=179.9999;}
 }
 int test_Q2=check_possibil_inter_Q2(Q);
 int test_W=check_possibil_inter_W(change_Q2(Q),W);
 double changeW=W;
 int trash1=interp_num1;
 int trash2=interp_num2;
 double tmp_val1,tmp_val2,tmp_val3,tmp_val0,tmp_val;
 if ((trash1==trash2)&&(test_W==7)){changeW=W_ext_min;}
 if ((trash1==trash2)&&(test_W==8)){changeW=W_ext_max;}
 if ((trash1!=trash2)&&(test_W==7)){changeW=W_ext_min;}
 if ((trash1!=trash2)&&(test_W==8)){changeW=W_ext_max;}
//cout<<" Q: "<<Q<<" W:"<<W<<" cos: "<<cos<<" fi: "<<fi<<" Eb: "<<Ebeam<<endl;

 if ((test_Q2==0)&&((test_W==5)||(test_W==6))){return intrep_CS_part(Q,W,cos,fi,1);}

 if ((test_Q2==1)&&((test_W==5)||(test_W==6))) { 
  tmp_val1=intrep_CS_part(change_Q2(Q),W,cos,fi,1);
  tmp_val2=get_CS_ph(W,cos)/6.283;
  return lin_interp(Q,0,Qmin,tmp_val2,tmp_val1);
 }

 if ((test_Q2==0)&&(test_W==0)) { 
  double tmp_res=intrep_CS(Q,W,cos,fi,1);
  if(tmp_res>=0) return tmp_res;
  else {/*
			double t_fi=fi;
			if (fi>40) {
				t_fi=fi-40;
				double res1=intrep_CS(Q,W,cos,t_fi,1);
				if (res1>=0) return res1;
				else{
					if (t_fi<320){
						t_fi=fi+40;
						double res2=intrep_CS(Q,W,cos,t_fi,1);
						if (res2>=0) return res2;
					}
				}
			}
			else t_fi=fi+40;
			//cout<<" change_fi ";
			double res=intrep_CS(Q,W,cos,t_fi,1);
			if (res>=0) return res;
			if (abs(cos+0.1)<=1) res=intrep_CS(Q,W,cos+0.1,fi,1);
			if (res>=0) return res;
			if (abs(cos-0.1)<=1) res=intrep_CS(Q,W,cos-0.1,fi,1);
			if (res>=0) return res;
			if (abs(W+0.1)<=5) res=intrep_CS(Q,W+0.1,cos,fi,1);
			if (res>=0) return res;
			if (abs(Q+0.1)<=12) res=intrep_CS(Q+0.1,W,cos,fi,1);
			if (res>=0) return res;
			if (abs(W+0.15)<=5) res=intrep_CS(Q,W+0.15,cos,fi,1);
			if (res>=0) return res;*/
			//cout<<"ERROR:( CS<0 value: "<<tmp_res<<" Q: "<<Q<<" W: "<<W<<" cos: "<<cos<<" fi: "<<fi<<endl; return 0;
			return abs(tmp_res);
		}
	
 }
 if ((test_Q2==1)&&(test_W==7)){
//cout<<"tp1"<<endl;
  tmp_val1=intrep_CS_part(change_Q2(Q),changeW,cos,fi,1);
  tmp_val2=get_CS_ph(W,cos)/6.283;
  tmp_val1=tmp_val1/(changeW-porog_ch(type_chanel));
  tmp_val3=tmp_val1*(W-porog_ch(type_chanel));
  return lin_interp(Q,0,Qmin,tmp_val2,tmp_val3);
  }//photo and W, 2-dimensial lineal interpolation

  if ((test_Q2==1)&&(test_W==0)) {
   tmp_val1=intrep_CS_part(change_Q2(Q),changeW,cos,fi,1);
   tmp_val2=get_CS_ph(W,cos)/6.283;
   return lin_interp(Q,0,Qmin,tmp_val2,tmp_val1);
  }//photo interp with form like Qmin

 if ((test_W==7)&&(test_Q2==0)) {
//cout<<"tp2"<<endl;
  tmp_val=intrep_CS_part(change_Q2(Q),changeW,cos,fi,1);
  tmp_val=tmp_val/(changeW-porog_ch(type_chanel));
  //cout<<" W: "<<W<<" changeW: "<<changeW<<" porog: "<<porog_ch(type_chanel)<<endl;
  return tmp_val*(W-porog_ch(type_chanel));
 }//lin do poroga only W

 if ((test_W==8)&&(test_Q2==0)) {
  tmp_val1=intrep_CS_part(change_Q2(Q),changeW,cos,fi,1);
  tmp_val2=get_CS_ph(W,cos)/6.283;
  tmp_val3=get_CS_ph(changeW,cos)/6.283;
  return tmp_val2*tmp_val1/tmp_val3;
 }//max and photo then lin

 if ((test_W==8)&&(test_Q2==1)) {
  //cout<<" test_W==4 && test_Q2==1 "<<" changeW: "<<changeW<<endl;
  tmp_val1=intrep_CS_part(change_Q2(Q),changeW,cos,fi,1);
  tmp_val2=get_CS_ph(W,cos)/6.283;
  tmp_val0=get_CS_ph(W,cos)/6.283;
  tmp_val3=get_CS_ph(changeW,cos)/6.283;
  tmp_val2=tmp_val2*tmp_val1/tmp_val3;
  return lin_interp(Q,0,Qmin,tmp_val0,tmp_val2);
 }//max and photo then lin plus Q2 lineal

 if (test_Q2==2) {cout<<"there is an error in connection CS in the two part of CS"<<endl; return 0; }
 cout<<" error at low Q2, i dont know why"<<endl;
 return 0;
}

double Sigma::get_d4CS(double Q,double W, double cos, double Ebeam){
	//cout<<" cos2: "<<_cos[3380]<<endl;

	int test_Q2=check_possibil_inter_Q2(Q);
	int test_W=check_possibil_inter_W(change_Q2(Q),W);
	double changeW=W;
	int trash1=interp_num1;
	int trash2=interp_num2;
	double tmp_val1,tmp_val2,tmp_val3,tmp_val0,tmp_val;
	if ((trash1==trash2)&&(test_W==3)){changeW=_W_min_all;}
	if ((trash1==trash2)&&(test_W==4)){changeW=_W_max[trash1];}
	if ((trash1!=trash2)&&(test_W==3)){changeW=_W_min_all;}
	if ((trash1!=trash2)&&(test_W==4)){changeW=_W_max_all;}

	if ((test_Q2==0)&&(test_W==0)) { 
		double tmp_res=intrep_CS(Q,W,cos,0,2);
		if(tmp_res>=0) return tmp_res;
		else {cout<<" STRANGE error tmp_res: "<<tmp_res<<" Q: "<<Q<<" W: "<<W<<"cos: "<<cos<<endl; return 0;}
	}
	if ((test_Q2==1)&&(test_W==3)){
		tmp_val1=intrep_CS(change_Q2(Q),changeW,cos,1,2);
		tmp_val2=get_CS_ph(W,cos);
		tmp_val1=tmp_val1/(changeW-porog_ch(type_chanel));
		tmp_val3=tmp_val1*(W-porog_ch(type_chanel));
		return lin_interp(Q,0,Qmin,tmp_val2,tmp_val3);

	}//photo and W, 2-dimensial lineal interpolation
	if ((test_Q2==1)&&(test_W==0)) {
		tmp_val1=intrep_CS(change_Q2(Q),changeW,cos,1,2);
		tmp_val2=get_CS_ph(W,cos);
		return lin_interp(Q,0,Qmin,tmp_val2,tmp_val1);
	}//photo interp with form like Qmin
if ((test_W==3)&&(test_Q2==0)) {
		tmp_val=intrep_CS(change_Q2(Q),changeW,cos,1,2);
		tmp_val=tmp_val/(changeW-porog_ch(type_chanel));
		//cout<<" W: "<<W<<" changeW: "<<changeW<<" porog: "<<porog_ch(type_chanel)<<endl;
		return tmp_val*(W-porog_ch(type_chanel));
	}//lin do poroga only W
	if ((test_W==4)&&(test_Q2==0)) {
		tmp_val1=intrep_CS(change_Q2(Q),changeW,cos,1,2);
		tmp_val2=get_CS_ph(W,cos);
		tmp_val3=get_CS_ph(changeW,cos);
		return tmp_val2*tmp_val1/tmp_val3;
	}//max and photo then lin
	if ((test_W==4)&&(test_Q2==1)) {
		//cout<<" test_W==4 && test_Q2==1 "<<" changeW: "<<changeW<<endl;
		tmp_val1=intrep_CS(change_Q2(Q),changeW,cos,1,2);
		tmp_val2=get_CS_ph(W,cos);
		tmp_val0=get_CS_ph(W,cos);
		tmp_val3=get_CS_ph(changeW,cos);
		tmp_val2=tmp_val2*tmp_val1/tmp_val3;
		return lin_interp(Q,0,Qmin,tmp_val0,tmp_val2);
	}//max and photo then lin plus Q2 lineal
	if (test_Q2==2) {cout<<"there is an error in connection CS in the two part of CS"<<endl; return 0; }
	cout<<" error at low Q2, i dont know why"<<endl;
	return 0;
}

double Sigma::getCS_d3CS(double Q,double W, double Ebeam){
//cout<<" Q="<<Q<<" W="<<W<<endl;
	bool check=check_kin(Q,W,Ebeam);
	if (check==0) {cout<<" incorrect input_data"<<endl; return 0;}
	if (porog_ch(type_chanel)>W){cout<<" too low W. threshold: "<<porog_ch(type_chanel)<<endl; return 0;}
	double tmp_cache=0,tmp_cache2=0;
	//cout<<" cos3: "<<_cos[3380]<<endl;
	int test_Q2=check_possibil_inter_Q2(Q);
	int test_W=check_possibil_inter_W(change_Q2(Q),W);
	double changeW=W;
	int trash1=interp_num1;
	int trash2=interp_num2;
	double tmp_val1,tmp_val2,tmp_val3,tmp_val0,tmp_val;
	if ((trash1==trash2)&&(test_W==3)){changeW=_W_min_all;}
	if ((trash1==trash2)&&(test_W==4)){changeW=_W_max[trash1];}
	if ((trash1!=trash2)&&(test_W==3)){changeW=_W_min_all;}
	if ((trash1!=trash2)&&(test_W==4)){changeW=_W_max_all;}

	tmp_cache=0;
	tmp_cache+=int_cos(change_Q2(Q),changeW);
	//cout<<"tmp_cache: "<<tmp_cache<<endl;
	if ((test_Q2==0)&&(test_W==0)) { 
		//cout<<" )-"<<endl;
		return tmp_cache;
	}
	if ((test_Q2==1)&&(test_W==3)){
		//cout<<"1"<<endl;
		tmp_val1=tmp_cache;
		tmp_val2=get_CS_ph_int(W);
		tmp_val1=tmp_val1/(changeW-porog_ch(type_chanel));
		tmp_val3=tmp_val1*(W-porog_ch(type_chanel));
		return lin_interp(Q,0,Qmin,tmp_val2,tmp_val3);

	}//photo and W, 2-dimensial lineal interpolation
	if ((test_Q2==1)&&(test_W==0)) {
		//cout<<"2"<<endl;
		tmp_val1=tmp_cache;
		tmp_val2=get_CS_ph_int(W);
		return lin_interp(Q,0,Qmin,tmp_val2,tmp_val1);
	}//photo interp with form like Qmin
if ((test_W==3)&&(test_Q2==0)) {
		//cout<<"3"<<endl;
		tmp_val=tmp_cache;
		tmp_val=tmp_val/(changeW-porog_ch(type_chanel));
		//tmp_val=tmp_val/delta_W(W);
		//cout<<" W: "<<W<<" changeW: "<<changeW<<" porog: "<<porog_ch(type_chanel)<<endl;
		//cout<<" val: "<<tmp_cache<<" changeW: "<<changeW<<" w-por: "<<W-porog_ch(type_chanel)<<" w2-w1: "<<changeW-porog_ch(type_chanel)<<endl;
		return tmp_val*(W-porog_ch(type_chanel));
	}//lin do poroga only W
	if ((test_W==4)&&(test_Q2==0)) {
		//cout<<"4"<<endl;
		tmp_val1=tmp_cache;
		tmp_val2=get_CS_ph_int(W);
		tmp_val3=get_CS_ph_int(changeW);
		return tmp_val2*tmp_val1/tmp_val3;
	}//max and photo then lin
	if ((test_W==4)&&(test_Q2==1)) {
		//cout<<"5"<<endl;
		//cout<<" test_W==4 && test_Q2==1 "<<" changeW: "<<changeW<<endl;
		tmp_val1=tmp_cache;
		tmp_val2=get_CS_ph_int(W);
		tmp_val0=get_CS_ph_int(W);
		tmp_val3=get_CS_ph_int(changeW);
		tmp_val2=tmp_val2*tmp_val1/tmp_val3;
		return lin_interp(Q,0,Qmin,tmp_val0,tmp_val2);
	}//max and photo then lin plus Q2 lineal
	return 0;
}


double Sigma::anti_Fit(double fi, double p0,double p1,double p2, int val){//anti_fit and integr po fi 1-anti_fit, 2-untegr po fi
	if(val==1) return p0+p1*cos(2*fi/57.29578049)+p2*cos(fi/57.29578049);
	if(val==2) return 2*3.1415926*p0;
	cout<<"incorrect param anti_fit"<<endl;
	return 0;
}

//ERROR int 3-30 apperas if I run  W_test and Q2 test
double Sigma::cos_in(int sp, double cos,double fi, int type){
	for (int i=sp;_W[i]==_W[sp];i++){
		if (abs(cos-_cos[i])<0.01) {// cout<<" ct=1 "; 
				double tmp= anti_Fit(fi,_p0[i],_p1[i],_p2[i],type);
				//cout<<" cos: "<<cos<<" fi: "<<fi<<" ty: "<<type<<" cos_val: "<<tmp<<" _p0[i]: "<<_p0[i]<<endl;
				return tmp;
		}
		if ((cos>_cos[i])&&(cos<_cos[i+1])) {
			//cout<<" cos type=2 "<<" cos: "<<cos<<" fi: "<<fi<<" ty: "<<type<<" cos_val(interp): "<<lin_interp(cos,_cos[i],_cos[i+1],anti_Fit(fi,_p0[i],_p1[i],_p2[i],type),anti_Fit(fi,_p0[i+1],_p1[i+1],_p2[i+1],type))<<" _p0[i]: "<<_p0[i]<<endl;
			return lin_interp(cos,_cos[i],_cos[i+1],anti_Fit(fi,_p0[i],_p1[i],_p2[i],type),anti_Fit(fi,_p0[i+1],_p1[i+1],_p2[i+1],type));
		}
	}
	/*
	cout<<"ERROR int 3-30"<<" cos: "<<cos<<" fi: "<<fi<<" int_sp: "<<sp<<endl; 
	for (int i=sp;_W[i]==_W[sp];i++){
			cout<<"cos: "<<cos<<" _cos["<<i<<"]: "<<_cos[i]<<" Q[i]: "<<_Q2[i]<<" W[i]: "<<_W[i]<<" _p0[i] "<<_p0[i]<<endl;
	}
	*/
	return 0;
}
double Sigma::W_in(int sp,double W, double cos,double fi, int type){
	for (int i=sp;_Q2[i]==_Q2[sp];i++){
		if (W==_W[i]) {//cout<<" W type=1 W: "<<W<<endl;
			return cos_in(i,cos,fi,type);}
		if ((W>_W[i])&&(W<_W[i+1])) {//cout<<" W type=2 _W[i]: "<<_W[i]<<" _W[i+1]"<<_W[i+1]<<endl; 
			int tmp=0;
			for (int j=0;_W[i-j]==_W[i];j++){tmp=j;}
			double tmp2=cos_in(i+1,cos,fi,type),tmp1=cos_in(i-tmp,cos,fi,type);
			
			return lin_interp(W,_W[i],_W[i+1],tmp1,tmp2);
		}
	}
	cout<<"ERROR int 3-4"<<endl; return 0;
}
double Sigma::intrep_CS_part(double Q,double W,double cos,double fi, int type_CS){
	int sp_Q=-1,sp_Q1=-1,sp_Q2=-1,tp_int0=-1;
	double W1=-1,W2=-1;
	for (int i=0;i<n_str_CS;i++){
		if (Q==_Q2[i]) {sp_Q=i; tp_int0=1; break;}
		if ((Q>_Q2[i])&&(Q<_Q2[i+1])) {sp_Q2=i+1; tp_int0=2; 
			int tmp=0;
			for (int j=0;((_Q2[i-j]==_Q2[i])&&((i-j)>=0));j++){tmp=j;}
			sp_Q1=i-tmp;
			//cout<<" SP: "<<sp_Q1<<" i: "<<i<<" tmp: "<<tmp<<" _Q2[i]: "<<_Q2[i]<<" _Q2[sp_Q1]: "<<_Q2[sp_Q1]<<" _Q2[sp_Q1-1]: "<<_Q2[sp_Q1-1]<<endl;
		break;}
	}
	if (((sp_Q==-1)&&(sp_Q1==-1))||(tp_int0==-1)) {cout<<" ERROR int3-1 "<<endl; return 0;}
	if (tp_int0==1) {
		return W1=W_in_part(sp_Q,W,cos,fi,type_CS);
	}
	if (tp_int0==2) {
		//cout<<"!!!!!!!!!!!!!!! ERROR !!!!!!!!!!!!!!!!!!!!!!"<<endl;
		W1=W_in_part(sp_Q1,W,cos,fi,type_CS);
		W2=W_in_part(sp_Q2,W,cos,fi,type_CS);
		return lin_interp(Q,_Q2[sp_Q1],_Q2[sp_Q2],W1,W2);
	}
	cout<<"ERROR int final Q: "<<Q<<" W: "<<W<<" cos: "<<cos<<" fi: "<<fi<<endl; 
	cout<<" sp_Q: "<<sp_Q<<" sp_Q1: "<<sp_Q1<<" tp_int0: "<<tp_int0<<endl; 
	return 0;
}
double Sigma::W_in_part(int sp,double W, double cos,double fi, int type){
	if(W<_W[sp]){

		double tmp_val1=cos_in(sp,cos,fi,type);
		double tmp_val2=get_CS_ph(W,cos)/6.283;
		double tmp_val3=get_CS_ph(_W[sp],cos)/6.283;
		return tmp_val2=tmp_val2*tmp_val1/tmp_val3;
		//tmp_val1=tmp_val1/(_W[sp]-porog_ch(type_chanel));
		//return tmp_val1*(W-porog_ch(type_chanel));

	}
	int chet=-1,chet2=0,pr=sp;
	for (int i=sp;_Q2[i]==_Q2[sp];i++){
		//cout<<" tp2=1"<<endl;
		chet++;
		//cout<<" W[i]: "<<_W[i]<<" _Q2[i]: "<<_Q2[i]<<" cos[i]: "<<_cos[i]<<endl;
		if (W==_W[i]) {//cout<<" W type=1 W: "<<W<<endl;
			double tmp_ch=cos_in(i,cos,fi,type);
			//cout<<" tp2=1 final"<<endl;
			return tmp_ch;}
		if ((W>_W[i])&&(W<_W[i+1])) {//cout<<" W type=2 _W[i]: "<<_W[i]<<" _W[i+1]"<<_W[i+1]<<endl; 
			int tmp=0;
			for (int j=0;_W[i-j]==_W[i];j++){tmp=j;}
			double tmp2=cos_in(i+1,cos,fi,type),tmp1=cos_in(i-tmp,cos,fi,type);

			return lin_interp(W,_W[i],_W[i+1],tmp1,tmp2);
		}
		if (_W[pr]==_W[i]){chet2++;}
		else{chet2=1;}
		//cout<<" W[i]: "<<_W[i]<<" _W[pr]: "<<_W[pr]<<" chet2: "<<chet2<<endl;
		pr=i;

	}
	chet2--;
	//cout<<"cos: "<<cos<<" _cos["<<sp+chet-chet2<<"]: "<<_cos[sp+chet-chet2]<<" Q[i]: "<<_Q2[sp+chet-chet2]<<" W[i]: "<<_W[sp+chet-chet2]<<" _p0[i] "<<_p0[sp+chet-chet2]<<" ch2: "<<chet2<<endl;
	if ((_Q2[sp+chet-chet2]==_Q2[sp])&&(_Q2[sp+chet+1]!=_Q2[sp])&&(_W[sp+chet-chet2-1]!=_W[sp+chet-chet2])){

		if(W>_W[sp+chet-chet2]){
			double tmp_val1=cos_in(sp+chet-chet2,cos,fi,type);
			double tmp_val2=get_CS_ph(W,cos)/6.283;
			double tmp_val3=get_CS_ph(_W[sp+chet-chet2],cos)/6.283;
			//if (W<W_ext_max) ph_fac=1.3;
			//else ph_fac=1;
			return tmp_val2=tmp_val2*tmp_val1/tmp_val3;
		}
	}else {cout<<"error in caclucl chet or chet2"<<endl;}

	cout<<"ERROR int 3-4"<<endl; return 0;
}
double Sigma::intrep_CS(double Q,double W,double cos,double fi, int type_CS){
	int sp_Q=-1,sp_Q1=-1,sp_Q2=-1,tp_int0=-1;
	double W1=-1,W2=-1;

	for (int i=0;i<n_str_CS;i++){
		if (Q==_Q2[i]) {sp_Q=i; tp_int0=1; break;}
		if ((Q>_Q2[i])&&(Q<_Q2[i+1])) {sp_Q2=i+1; tp_int0=2; 
			int tmp=0;
			for (int j=0;_Q2[i-j]==_Q2[i];j++){tmp=j;}
			sp_Q1=i-tmp;
		break;}
	}

	//cout<<" Q: "<<Q<<" W: "<<W<<" cos: "<<cos<<" fi: "<<fi<<endl; 
	//cout<<" sp_Q: "<<sp_Q<<" sp_Q1: "<<sp_Q1<<" tp_int0: "<<tp_int0<<endl; 

	if (((sp_Q==-1)&&(sp_Q1==-1))||(tp_int0==-1)) {cout<<" ERROR int3-1 "<<endl; return 0;}

	if (tp_int0==1) {
		return W1=W_in(sp_Q,W,cos,fi,type_CS);
	}
	if (tp_int0==2) {
		//cout<<"!!!!!!!!!!!!!!! ERROR !!!!!!!!!!!!!!!!!!!!!!"<<endl;
		W1=W_in(sp_Q1,W,cos,fi,type_CS);
		W2=W_in(sp_Q2,W,cos,fi,type_CS);
		return lin_interp(Q,_Q2[sp_Q1],_Q2[sp_Q2],W1,W2);
	}
	cout<<"ERROR int final Q: "<<Q<<" W: "<<W<<" cos: "<<cos<<" fi: "<<fi<<endl; 
	cout<<" sp_Q: "<<sp_Q<<" sp_Q1: "<<sp_Q1<<" tp_int0: "<<tp_int0<<endl; 
	
	return 0;
}

double min2(double a, double b){
	if (a<=b) return a; 
	else return b;
return 0;
}
double max2(double a, double b){
	if (a>=b) return a; 
	else return b;
return 0;
}
double min3(double a, double b, double c)
{
	if (a<=b)
	{
		if (a<=c) {return a;}
		else return c;
	}
	else 
	{
		if (b<=c) {return b;}
		else return c;
	}
}
double Sigma::lin_interp(double x,double point1, double point2, double value_point1, double value_point2){

	if (point1==point2){cout<<"ERROR 1"<<endl; return 0;}
	if ((point2>=x)&&(x>=point1)) {double tmp_r1=value_point1+(x-point1)*(value_point2-value_point1)/(point2-point1);
		//cout<<" x: "<<x<<" point1: "<<point1<<" point2: "<<point2<<" res: "<<tmp_r1<<endl;
		return tmp_r1;
	}
	cout<<" incorrect input data in Lin_interp"<<endl;
	cout<<" x: "<<x<<" point1: "<<point1<<" point2: "<<point2<<endl;
	return 0;
}
double Sigma::get_CS_ph(double W, double cos){
	int tmp0=range_ph(W);
	if (tmp0==0) return ph_int(W,cos)*get_CS_ph_int_Ev(W)/get_CS_ph_int(W);
	if (tmp0==1){
		double tmp_val=ph_int(min_W_ph,cos);
		tmp_val=tmp_val/(min_W_ph-porog_ch(type_chanel));
		//cout<<" W: "<<W<<" changeW: "<<changeW<<" porog: "<<porog_ch(type_chanel)<<endl;
		double tmp2= tmp_val*(W-porog_ch(type_chanel));
		//cout<<" TMP2: "<<tmp2<<endl; return tmp2;
		return tmp2*get_CS_ph_int_Ev(W)/get_CS_ph_int(W);
	}
	cout<<"ERROR in ph inter num 0"<<endl;
	return 0;
}
double Sigma::get_CS_ph_int(double W){
	int tmp0=range_ph(W);
	if (tmp0==0) return ph_int_int(W);
	if (tmp0==1){
		double tmp_val=ph_int_int(min_W_ph);
		tmp_val=tmp_val/(min_W_ph-porog_ch(type_chanel));
		//cout<<" W: "<<W<<" changeW: "<<changeW<<" porog: "<<porog_ch(type_chanel)<<endl;
		double tmp2= tmp_val*(W-porog_ch(type_chanel));
		//cout<<" TMP2: "<<tmp2<<endl; return tmp2;
		return tmp2;
	}
	cout<<"ERROR in ph inter int num 0"<<endl;
	return 0;
}
double Sigma::get_CS_ph_int_Ev(double W){
	int tmp0=range_ph(W);
	if (tmp0==0) return ph_int_int_Ev(W);
	if (tmp0==1){
		double tmp_val=ph_int_int_Ev(min_W_ph);
		tmp_val=tmp_val/(min_W_ph-porog_ch(type_chanel));
		//cout<<" W: "<<W<<" changeW: "<<changeW<<" porog: "<<porog_ch(type_chanel)<<endl;
		double tmp2= tmp_val*(W-porog_ch(type_chanel));
		//cout<<" TMP2: "<<tmp2<<endl; return tmp2;
		return tmp2;
	}
	cout<<"ERROR in ph inter int num 0"<<endl;
	return 0;
}
double Sigma::cos_in_ph(int sp, double cos){
	for (int i=sp;W_vec_ph[i]==W_vec_ph[sp];i++){
		if (abs(cos-costeta_vec_ph[i])<0.01) { //cout<<" cos type=1 "; 
				double tmp = CS_vec_ph[i];
				//cout<<" cos: "<<cos<<" fi: "<<fi<<" ty: "<<type<<" cos_val: "<<tmp<<endl;
				return tmp;
		}
		if ((cos>costeta_vec_ph[i])&&(cos<costeta_vec_ph[i+1])) {
//cout<<" cos type=2 ";
			return lin_interp(cos,costeta_vec_ph[i],costeta_vec_ph[i+1],CS_vec_ph[i],CS_vec_ph[i+1]);
		}
	}
	cout<<"ERROR int ph 3-3"<<" cos: "<<cos<<endl; 
	for (int i=sp;_W[i]==_W[sp];i++){
			//cout<<"cos: "<<cos<<" _cos["<<i<<"]: "<<_cos[i]<<" Q[i]: "<<_Q2[i]<<" W[i]: "<<_W[i]<<" _p0[i] "<<_p0[i]<<endl;
	}
	return 0;
}
double Sigma::ph_int(double W, double cos){

	int tmp0=range_ph(W);
	if (tmp0!=0){
		//if (tmp0==1) 
		{cout<<"error in ph int, too short W range in the photo data"<<endl; return 0;}
	}

	int sp_W=-1,sp_W1=-1,sp_W2=-1,tp_int0=-1;
	double W1=-1,W2=-1;
	for (int i=0;i<n_str_ph;i++){
		if (W==W_vec_ph[i]) {sp_W=i; tp_int0=1; break;}
		if ((W>W_vec_ph[i])&&(W<W_vec_ph[i+1])) {sp_W2=i+1; tp_int0=2; 
			int tmp=0;
			for (int j=0;W_vec_ph[i-j]==W_vec_ph[i];j++){tmp=j;}
			sp_W1=i-tmp;
		break;}
	}

	if (tp_int0==1) {
		return W1=cos_in_ph(sp_W,cos);
	}
	if (tp_int0==2) {
		//cout<<"!!!!!!!!!!!!!!! ERROR !!!!!!!!!!!!!!!!!!!!!!"<<endl;
		W1=cos_in_ph(sp_W1,cos);
		W2=cos_in_ph(sp_W2,cos);;
		return lin_interp(W,W_vec_ph[sp_W1],W_vec_ph[sp_W2],W1,W2);
	}
	cout<<"ERROR int final PH: "<<" W: "<<W<<" cos: "<<cos<<endl; 
	cout<<" sp_W: "<<sp_W<<" sp_W1: "<<sp_W1<<" tp_int0: "<<tp_int0<<endl; 
	
	return 0;
}
double Sigma::W_in_cs(int sp, double W){
	for (int i=sp;_Q_int[i]==_Q_int[sp];i++){
		if (abs(W-_W_int[i])<0.01) { //cout<<" cos type=1 "; 
				double tmp = _CS_int[i];
				//cout<<" cos: "<<cos<<" fi: "<<fi<<" ty: "<<type<<" cos_val: "<<tmp<<endl;
				return tmp;
		}
		if ((W>_W_int[i])&&(W<_W_int[i+1])) {
//cout<<" cos type=2 ";
			return lin_interp(W,_W_int[i],_W_int[i+1],_CS_int[i],_CS_int[i+1]);
		}
	}
	cout<<"ERROR int ph 3-3"<<" W: "<<W<<endl; 
	//for (int i=sp;_W[i]==_W[sp];i++){
			//cout<<"cos: "<<cos<<" _cos["<<i<<"]: "<<_cos[i]<<" Q[i]: "<<_Q2[i]<<" W[i]: "<<_W[i]<<" _p0[i] "<<_p0[i]<<endl;
	//}
	return 0;
}
double Sigma::int_cos(double Q, double W){

	int sp_Q=-1,sp_Q1=-1,sp_Q2=-1,tp_int0=-1;
	double Q1=-1,Q2=-1;
	for (int i=0;i<n_str_CS_int;i++){
		if (Q==_Q_int[i]) {sp_Q=i; tp_int0=1; break;}
		if ((Q>_Q_int[i])&&(Q<_Q_int[i+1])) {sp_Q2=i+1; tp_int0=2; 
			int tmp=0;
			for (int j=0;_Q_int[i-j]==_Q_int[i];j++){tmp=j;}
			sp_Q1=i-tmp;
		break;}
	}

	if (tp_int0==1) {
		return Q1=W_in_cs(sp_Q,W);
	}
	if (tp_int0==2) {
		//cout<<"!!!!!!!!!!!!!!! ERROR !!!!!!!!!!!!!!!!!!!!!!"<<endl;
		Q1=W_in_cs(sp_Q1,W);
		Q2=W_in_cs(sp_Q2,W);;
		return lin_interp(Q,_Q_int[sp_Q1],_Q_int[sp_Q2],Q1,Q2);
	}
	cout<<"ERROR int cos cs: "<<" Q: "<<Q<<" W: "<<W<<endl; 
	cout<<" sp_Q: "<<sp_Q<<" sp_Q1: "<<sp_Q1<<" tp_int0: "<<tp_int0<<endl; 
	
	return 0;
}
double Sigma::ph_int_int(double W){
	for(int i=0;i<n_str_ph_int;i++){
		if (W==_W_ph[i]) {return _CS_ph[i];}
		if ((W>_W_ph[i])&&(W<_W_ph[i+1])) return lin_interp(W,_W_ph[i],_W_ph[i+1],_CS_ph[i],_CS_ph[i+1]);
	}
	return 0;
}
double Sigma::ph_int_int_Ev(double W){
	for(int i=0;i<n_str_Ev;i++){
		if (W==_W_ph_Ev[i]) {return _CS_ph_Ev[i];}
		if ((W>_W_ph_Ev[i])&&(W<_W_ph_Ev[i+1])) return lin_interp(W,_W_ph_Ev[i],_W_ph_Ev[i+1],_CS_ph_Ev[i],_CS_ph_Ev[i+1]);
	}
	return 0;
}
#endif 
