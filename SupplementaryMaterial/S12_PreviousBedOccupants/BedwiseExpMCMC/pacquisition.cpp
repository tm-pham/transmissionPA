#include <algorithm>
#include <math.h>
#include "stdhdr.h"
#include <boost/math/special_functions/gamma.hpp>

using namespace std;
// ADDED environmental contamination E
// ADDED exact computation of probabilities

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

double pacquisition(parameters param,int I,int N,double E,int addingenvironment,int bedwise,int previous,double delta_t)
{
	double foi;
    if(param.frequencydependence==1)  foi1(param,I,N,E,addingenvironment,foi,bedwise,previous,delta_t);
    if(param.frequencydependence==0)  foi0(param,I,N,E,addingenvironment,foi,bedwise,previous,delta_t);
    //cout << "Result in pacquisition = " << 1.0 -exp(-foi) << endl;
    return 1.0 -exp(-foi);
}

double pnoacquisition(parameters param,int I,int N,double E,int addingenvironment,int bedwise,int previous,double delta_t)
{
	double foi;
    if(param.frequencydependence==1)  foi1(param,I,N,E,addingenvironment,foi,bedwise,previous,delta_t);
    if(param.frequencydependence==0)  foi0(param,I,N,E,addingenvironment,foi,bedwise,previous,delta_t);
    //cout << "Result in pnoacquisition = " << exp(-foi) << endl;
    return exp(-foi);
}

double logpacquisition(parameters param,int I,int N,double E,int addingenvironment,int bedwise,int previous,double delta_t)
{
	double foi;
    if(param.frequencydependence==1)  foi1(param,I,N,E,addingenvironment,foi,bedwise,previous,delta_t);
    if(param.frequencydependence==0)  foi0(param,I,N,E,addingenvironment,foi,bedwise,previous,delta_t);
    //cout << "FOI in logpacquisition = " << foi << endl;
    //cout << "Result in logpacquisition = " << log(1.0 -exp(-foi)) << endl;
	if(foi == 0.0){
	    return -1000000.0;
	}
	double logl = log(1.0 -exp(-foi));
    if(isnan(foi) || isinf(foi) || isnan(logl) || isinf(-foi) || isinf(-logl)){
		cout <<" foi="<<foi<<"; logl = "<< logl<< endl;
	}
    return log(1.0 -exp(-foi));
}

double logpnoacquisition(parameters param,int I,int N,double E,int addingenvironment,int bedwise,int previous,double delta_t)
{
	double foi;
    if(param.frequencydependence==1)  foi1(param,I,N,E,addingenvironment,foi,bedwise,previous,delta_t);
    if(param.frequencydependence==0)  foi0(param,I,N,E,addingenvironment,foi,bedwise,previous,delta_t);
    //cout << "FOI in logpnoacquisition = " << foi << endl;
    //cout << "Result in logpnoacquisition = " << -foi << endl;
    if(isnan(foi) || isinf(foi)){
		cout <<" foi="<<foi<<endl;
	}
    return -foi;
}


/*
double pacquisition(parameters param,int I,int N,double E,int addingenvironment)
{
	double endogenous=0.0, exogenous=0.0, env=0.0, pinfection=0.0;
	pacquisitionbyroute(param, I, N, E, endogenous, exogenous, env, pinfection);
    return pinfection;
}

double pnoacquisition(parameters param,int I,int N,double E,int addingenvironment)
{
	double endogenous=0.0, exogenous=0.0, env=0.0, pinfection=0.0;
	pacquisitionbyroute(param, I, N, E, endogenous, exogenous, env, pinfection);
    return 1.0- pinfection;
}

double logpacquisition(parameters param,int I,int N,double E,int addingenvironment)
{
	double endogenous=0.0, exogenous=0.0, env=0.0, pinfection=0.0;
	pacquisitionbyroute(param, I, N, E, endogenous, exogenous, env, pinfection);
    return log(pinfection);
}

double logpnoacquisition(parameters param,int I,int N,double E,int addingenvironment)
{
	double endogenous=0.0, exogenous=0.0, env=0.0, pinfection=0.0;
	pacquisitionbyroute(param, I, N, E, endogenous, exogenous, env, pinfection);
    return log(1.0- pinfection);
}
* */

// Adjustments for pacquisitionward still missing!

// Simple computation
/*
double pacquisition(parameters param,int I,int N,double E,int addingenvironment)
{
    if(param.frequencydependence==1)  return 1.0-exp(-param.a-param.b*(double)I/(double)N-(double)addingenvironment*param.c*E);
    if(param.frequencydependence==0)  return 1.0-exp(-param.a-param.b*(double)I-(double)addingenvironment*param.c*E);
}

double pnoacquisition(parameters param,int I,int N,double E,int addingenvironment)
{
    if(param.frequencydependence==1) return exp(-param.a-param.b*(double)I/(double)N-(double)addingenvironment*param.c*E);
    if(param.frequencydependence==0) return exp(-param.a-param.b*(double)I-(double)addingenvironment*param.c*E);
}

double logpacquisition(parameters param,int I,int N,double E,int addingenvironment)
{
    if(param.frequencydependence==1)return log(1.0-exp(-param.a-param.b*(double)I/(double)N-(double)addingenvironment*param.c*E));
    if(param.frequencydependence==0)return log(1.0-exp(-param.a-param.b*(double)I-(double)addingenvironment*param.c*E));
}

double logpnoacquisition(parameters param,int I,int N,double E,int addingenvironment)
{
    if(param.frequencydependence==1)return -param.a-param.b*(double)I/(double)N-(double)addingenvironment*param.c*E;
    if(param.frequencydependence==0)return -param.a-param.b*(double)I-(double)addingenvironment*param.c*E;
}
*/

double pacquisitionward(parameters param,int I,int N,int unit,double E,int addingenvironment)
{
    if(param.frequencydependence==1)return 1.0-exp(-param.alphavector.at(unit)-param.betavector.at(unit)*(double)I/(double)N-(double)addingenvironment*param.c*E);
    if(param.frequencydependence==0)return 1.0-exp(-param.alphavector.at(unit)-param.betavector.at(unit)*(double)I-(double)addingenvironment*param.c*E);
}

double pnoacquisitionward(parameters param,int I,int N,int unit,double E,int addingenvironment)
{
    if(param.frequencydependence==1)return exp(-param.alphavector.at(unit)-param.betavector.at(unit)*(double)I/(double)N-(double)addingenvironment*param.c*E);
    if(param.frequencydependence==0)return exp(-param.alphavector.at(unit)-param.betavector.at(unit)*(double)I-(double)addingenvironment*param.c*E);
}

double logpacquisitionward(parameters param,int I,int N,int unit,double E,int addingenvironment)
{
    if(param.frequencydependence==1)return log(1.0-exp(-param.alphavector.at(unit)-param.betavector.at(unit)*(double)I/(double)N-(double)addingenvironment*param.c*E));
    if(param.frequencydependence==0)return log(1.0-exp(-param.alphavector.at(unit)-param.betavector.at(unit)*(double)I-(double)addingenvironment*param.c*E));
}

double logpnoacquisitionward(parameters param,int I,int N,int unit,double E,int addingenvironment)
{
    if(param.frequencydependence==1)return -param.alphavector.at(unit)-param.betavector.at(unit)*(double)I/(double)N-(double)addingenvironment*param.c*E;
    if(param.frequencydependence==0)return -param.alphavector.at(unit)-param.betavector.at(unit)*(double)I-(double)addingenvironment*param.c*E;
}

// OLD: DO NOT DELETE
double pacquisition(parameters param,int I,int N)
{
    if(param.frequencydependence==1)  return 1.0-exp(-param.a-param.b*(double)I/(double)N);
    if(param.frequencydependence==0)  return 1.0-exp(-param.a-param.b*(double)I);
}

double pnoacquisition(parameters param,int I,int N)
{
    if(param.frequencydependence==1) return exp(-param.a-param.b*(double)I/(double)N);
    if(param.frequencydependence==0) return exp(-param.a-param.b*(double)I);
}

double logpacquisition(parameters param,int I,int N)
{
    if(param.frequencydependence==1)return log(1.0-exp(-param.a-param.b*(double)I/(double)N));
    if(param.frequencydependence==0)return log(1.0-exp(-param.a-param.b*(double)I));
}

double logpnoacquisition(parameters param,int I,int N)
{
    if(param.frequencydependence==1)return -param.a-param.b*(double)I/(double)N;
    if(param.frequencydependence==0)return -param.a-param.b*(double)I;
}


double pacquisitionward(parameters param,int I,int N,int unit)
{
    if(param.frequencydependence==1)return 1.0-exp(-param.alphavector.at(unit)-param.betavector.at(unit)*(double)I/(double)N);
    if(param.frequencydependence==0)return 1.0-exp(-param.alphavector.at(unit)-param.betavector.at(unit)*(double)I);
}

double pnoacquisitionward(parameters param,int I,int N,int unit)
{
    if(param.frequencydependence==1)return exp(-param.alphavector.at(unit)-param.betavector.at(unit)*(double)I/(double)N);
    if(param.frequencydependence==0)return exp(-param.alphavector.at(unit)-param.betavector.at(unit)*(double)I);
}

double logpacquisitionward(parameters param,int I,int N,int unit)
{
    if(param.frequencydependence==1)return log(1.0-exp(-param.alphavector.at(unit)-param.betavector.at(unit)*(double)I/(double)N));
    if(param.frequencydependence==0)return log(1.0-exp(-param.alphavector.at(unit)-param.betavector.at(unit)*(double)I));
}

double logpnoacquisitionward(parameters param,int I,int N,int unit)
{
    if(param.frequencydependence==1)return -param.alphavector.at(unit)-param.betavector.at(unit)*(double)I/(double)N;
    if(param.frequencydependence==0)return -param.alphavector.at(unit)-param.betavector.at(unit)*(double)I;
}
