#ifndef _DSS_PACQUISITION_H_
#define _DSS_PACQUISITION_H_

#include <algorithm>
#include <math.h>
#include "stdhdr.h"
#include <boost/math/special_functions/gamma.hpp>

using namespace std;

inline void pacquisitionbyroute(parameters param, int I, int N, double E, double &endogenous, double &exogenous, double &env, double &pinfection, double error);
inline void foi1(parameters param,int I,int N,double E,int addingenvironment,double &foi,int bedwise,int previous,double delta_t);
inline void foi0(parameters param,int I,int N,double E,int addingenvironment,double &foi,int bedwise,int previous,double delta_t);

double pacquisition(parameters param,int I,int N,double E, int addingenvironment,int bedwise,int previous,double delta_t);
double pnoacquisition(parameters param,int I,int N,double E, int addingenvironment,int bedwise,int previous,double delta_t);
double logpacquisition(parameters param,int I,int N,double E, int addingenvironment,int bedwise,int previous,double delta_t);
double logpnoacquisition(parameters param,int I,int N,double E, int addingenvironment,int bedwise,int previous,double delta_t);

double pacquisitionward(parameters param,int I,int N,int unit,double E, int addingenvironment);
double pnoacquisitionward(parameters param,int I,int N,int unit,double E, int addingenvironment);
double logpacquisitionward(parameters param,int I,int N,int unit,double E, int addingenvironment);
double logpnoacquisitionward(parameters param,int I,int N,int unit,double E, int addingenvironment);

double pacquisition(parameters param,int I,int N);
double pnoacquisition(parameters param,int I,int N);
double logpacquisition(parameters param,int I,int N);
double logpnoacquisition(parameters param,int I,int N);

double pacquisitionward(parameters param,int I,int N,int unit);
double pnoacquisitionward(parameters param,int I,int N,int unit);
double logpacquisitionward(parameters param,int I,int N,int unit);
double logpnoacquisitionward(parameters param,int I,int N,int unit);

inline void foi1(parameters param,int I,int N,double E,int addingenvironment,double &foi,int bedwise,int previous, double delta_t){
  if(param.mu==0.0 && param.a==0.0 && param.b==0.0 && param.c==0.0){
	foi = 0.0;
  } else{
    double A, B;
    double diff = ((param.nu/param.mu)*(double)I/(double)N-E);

    A = param.a + (param.b + (double)addingenvironment*param.c*param.nu/param.mu)*(double)I/(double)N + (double)bedwise*param.p*exp(-delta_t)*(double)previous;
    B = (param.c/param.mu)*diff;
    foi = A+(double)addingenvironment*B*(exp(-param.mu)-1.0);
    if(isnan(A) || isnan(B) || isnan(foi)){
	    cout << "E="<<E<<"; A="<<A<<"; B="<<B<<"; foi="<<foi<<"; param.a="<<param.a<<"; N="<<N<<endl;
	    cout << "param.b="<<param.b<<" ;param.c="<<param.c<<"; param.nu="<<param.nu<<"; param.mu="<<param.mu<<"; I="<<I<<endl;
	    exit(0);
    }
  }
}

inline void foi0(parameters param,int I,int N,double E,int addingenvironment,double &foi,int bedwise,int previous, double delta_t){
	double A, B;
	double diff = ((param.nu/param.mu)*(double)I-E);

	A = param.a + (param.b + (double)addingenvironment*param.c*param.nu/param.mu)*(double)I+ (double)bedwise*param.p*exp(-delta_t)*(double)previous;
	B = (param.c/param.mu)*diff;
  foi = A+(double)addingenvironment*B*(exp(-param.mu)-1.0);
}

#endif //_DSS_PACQUISITION_H_
