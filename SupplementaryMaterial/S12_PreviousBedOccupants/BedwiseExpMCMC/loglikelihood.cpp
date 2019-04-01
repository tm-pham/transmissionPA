#include <iostream>
#include <vector>
#include <math.h>
#include <algorithm>
#include "stdhdr.h"
#include "pacquisition.h"
using namespace std;
// ADDED: Vector for Environmental contamination E
double loglikelihood(vector <int> & acquisitionday,vector <int> & acquisition,vector< vector<int> >& culturearray,vector< vector<int> >& N,vector< vector<int> >& I,vector< vector<int> >& numberofacquisitions, vector <int> & admissionstate,parameters param,int startperiod, int endperiod,
                     vector <int> & admissionday,int multipleparameters,vector<vector<double> >& E, int addingenvironment, int bedwise, vector<int>& previousCol, vector<vector<vector<int> > >& present)
  {
    long double logl=0;
    double loglacq=0,loglnoacq=0,logladm=0,loglsens=0;
    double foipr=0, loglacq1=0, loglacq0=0, logl1=0;
    double delta_t;
    int patid,day,j,sumPat=0,pat=0,numnoacq=0;
    unsigned int unit,patient,culture;
    double ntotal=0;
    double itotal=0;
    if(multipleparameters==0){
      for(day=startperiod;day<endperiod;day++){
        for(unit=0;unit<N.size();unit++){
          itotal+=I.at(unit).at(day);
          ntotal+=N.at(unit).at(day);
          if(N.at(unit).at(day)-I.at(unit).at(day)>0){
            sumPat = present.at(unit).at(day).size();
            for(j=0;j<sumPat;j++){
              pat = present.at(unit).at(day).at(j);
              if(admissionstate.at(pat)==0){
                delta_t = day-admissionday.at(pat);
                foi1(param,I.at(unit).at(day),N.at(unit).at(day),E.at(unit).at(day),addingenvironment,foipr,bedwise,previousCol.at(pat),delta_t);
                if(acquisition.at(pat)==1){
                  if(acquisitionday.at(pat)+admissionday.at(pat)==day){
                    loglacq += log(1.0-exp(-foipr)); 
                  }else {
                    if(acquisitionday.at(pat)+admissionday.at(pat)>day) loglnoacq -= foipr;
                  }
                }else{
                  loglnoacq -= foipr;
                }
              }
            }
            //logl+= loglnoacq*(N.at(unit).at(day)-I.at(unit).at(day)-numberofacquisitions.at(unit).at(day));
            //loglnoacq=logpnoacquisition(param,I.at(unit).at(day),N.at(unit).at(day));
            //cout << "No acquisition logl="<<loglnoacq<<", ";
            //logl+= logpacquisition(param,I.at(unit).at(day),N.at(unit).at(day),E.at(unit).at(day),addingenvironment)*numberofacquisitions.at(unit).at(day);
            
            //loglacq=logpacquisition(param,I.at(unit).at(day),N.at(unit).at(day));
            //cout << "Acquisition logl="<<loglacq<<endl;
            if(isnan(logl)){
              cout << "Day="<<day << endl;
              cout << "I.at(unit).at(day)=" << I.at(unit).at(day)<< "; N.at(unit).at(day)=" << N.at(unit).at(day)<< "; E.at(unit).at(day)="<< E.at(unit).at(day)<< endl;
              cout << "param.a="<<param.a<<"; param.b="<<param.b<<"; param.c="<<param.c<<"; param.nu="<<param.nu<<"; param.mu="<<param.mu << endl;
              //cout << "No acquisition logl = " << logpnoacquisition(param,I.at(unit).at(day),N.at(unit).at(day),E.at(unit).at(day),addingenvironment)*(N.at(unit).at(day)-I.at(unit).at(day)-numberofacquisitions.at(unit).at(day))<<endl;
              //cout << "Acquisition logl = "<<logpacquisition(param,I.at(unit).at(day),N.at(unit).at(day),E.at(unit).at(day),addingenvironment) << endl;
              cout << "Number of acquisitions = " << numberofacquisitions.at(unit).at(day) << endl;
              exit(0);
            }
          }
        }
      }
      logl += loglnoacq + loglacq;
      //logl1 = logl;
    }else{
      //cout << "multiple"<<endl;
      for(day=startperiod;day<endperiod;day++)
      {
        //cout<<"day="<<day<<endl;
        for(unit=0;unit<N.size();unit++)
        {
          //cout<<"unit="<<unit<<endl;
          if(N.at(unit).at(day)-I.at(unit).at(day)>0)
          {
            //cout<<"a"<<endl;
            /** day added */  
            logl+= logpnoacquisitionward(param,I.at(unit).at(day),N.at(unit).at(day),unit,E.at(unit).at(day),addingenvironment)*(N.at(unit).at(day)-I.at(unit).at(day)-numberofacquisitions.at(unit).at(day));
            //cout<<"b"<<endl;
            /** day added */  
            logl+= logpacquisitionward(param,I.at(unit).at(day),N.at(unit).at(day),unit,E.at(unit).at(day),addingenvironment)*numberofacquisitions.at(unit).at(day);
            if(isnan(logl)){
              cout << "N.at(unit).at(day)=" << N.at(unit).at(day)<< endl;
            }
          }
        }
      }
    }
    /** Calclulate contribution being colonized on admission to loglikelihood */
    for(patient=0;patient<admissionstate.size();patient++)
    {
      if(admissionstate.at(patient)==1)
      {
        logl+=log(param.f);logladm+=log(param.f);
      }
      else
      {
        logl+=log(1.-param.f);logladm+=log(1.-param.f);
      }
    }
    for(culture=0;culture<culturearray.size();culture++)
    {
      if(culturearray.at(culture).at(1)>=startperiod&&culturearray.at(culture).at(1)<=endperiod)
      {
        if(culturearray.at(culture).at(2)==0)//culture negative
        {
          patid=culturearray.at(culture).at(0);
          if(admissionstate.at(patid)==1||(acquisition.at(patid)==1&&culturearray.at(culture).at(1)>admissionday.at(patid)+acquisitionday.at(patid)))
          {
            logl+=log(1.0-param.phi);//false negative culture
            loglsens+=log(1.0-param.phi);
            //cout << "False negative loglsens=" << loglsens << endl;
          }
          else
          {
            //culture negative and patient negative
          }
        }
        else //culture positive
        {
          patid=culturearray.at(culture).at(0);
          
          if(admissionstate.at(patid)==1||(acquisition.at(patid)==1&&culturearray.at(culture).at(1)>admissionday.at(patid)+acquisitionday.at(patid)))
          {
            logl+=log(param.phi);//positive culture for colonized patient
            loglsens+=log(param.phi);
            //cout << "True positive loglsens=" << loglsens << endl;
          }
          else
          {
            logl-=1000000.0; //culture positive and patient uncolonized
            cout << "case logl-=1000000.0: logl=" << logl << endl;
            //loglsens-=1000000.0;
          }
        }
      }
      
    }
    if(logl>0){
      cout <<  "logl=" << logl << endl;
      cout << "param.a, param.b, param.p, param.f, param.phi: " << endl;
      cout << param.a << ", " << param.b << ", " << param.p << ", " << param.f << ", " << param.phi << endl;
      cout << "I per day:" << endl;
      for(int d=startperiod;d<endperiod;d++){
        cout << I.at(0).at(d)<< ",";
      }
      cout << endl;
      cout << "N per day:" << endl;
      for(int d=startperiod;d<endperiod;d++){
        cout << N.at(0).at(d)<< ",";
      }
      cout<<endl;
      cout << "Numberofacquisitions per day:" << endl;
      for(int d=startperiod;d<endperiod;d++){
        cout << numberofacquisitions.at(0).at(d)<< ",";
      }
      cout<<endl;
      int numberofpatients = acquisitionday.size();
      cout << "Previous colonized bed opccupant?" << endl;
      for(pat=0;pat<numberofpatients;pat++){
        cout << previousCol.at(pat) << ",";
      }
      cout << endl;
      
      cout << "Importation?" << endl;
      for(pat=0;pat<numberofpatients;pat++){
        cout << admissionstate.at(pat) << ",";
      }
      cout << endl;
      
      cout << "Acquisition time:" << endl;
      for(pat=0;pat<numberofpatients;pat++){
        if(acquisition.at(pat)==1) cout << admissionday.at(pat) + acquisitionday.at(pat) <<",";
        else cout << "-1,";
      }
      cout << endl;
      cout<<"logl1="<<logl1<<", logladm="<<logladm<<", loglsens="<<loglsens<<endl;
      
    }
    //cout <<"mean prev="<< itotal/ntotal<<endl;
    //cout<<"logl="<<logl<<endl;//", loglacq="<<loglacq<<", loglnoacq="<<loglnoacq<<"logladm="<<logladm<<", loglsens="<<loglsens<<endl;
    return logl;
  }