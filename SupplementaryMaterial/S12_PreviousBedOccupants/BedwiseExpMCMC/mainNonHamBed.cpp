#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <fstream>
#include <time.h>
#include "stdhdr.h"
#include "readcsvfileto2Darrayofint.h"
#include "relabel.h"
#include "pacquisition.h"
#include "loglikelihood.h"
#include "changeloglandI.h"
#include "createusefuldata.h"
#include "truepositiveandfalsenegative.h"
#include <boost/math/distributions/exponential.hpp>
#include <boost/math/distributions/uniform.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include "samplefromnormal.h"
#include "environment.h"
#include <time.h>

// MCMC SIMULATION FOR NEW ENVIRONMENTAL CONTRIBUTION

using namespace std;
/***********************************************************************
* Sorting function
* sort based on second element
***********************************************************************/
bool sortFunc1( const vector<int>& p1,const vector<int>& p2 )
{
  return p1.at(1) < p2.at(1);
}

/***********************************************************************
* Main function
***********************************************************************/

int main()
{
  /**********************************************************************
  * DECLARE VARIABLES
  **********************************************************************/
  time_t timer = time(NULL);
  // extern int parallelism_enabled=0; // Multithread if = 1
  
  int inew;
  int numberimportations;
  double forceofinfection;
  int nchange=0;
  double error = 1e-11;
  
  /// First defined in PROCESSING INPUT DATA
  int numberofpatients,numberofwards,maxdate,unit;
  
  /// Declare variables for data augmentation
  int factoraugmented,changeaugmented,newfactoraugmented;
  
  double oldlogl;
  
  int newadmissionstate,newacquisitionday,newacquisition,oldadmissionstate,oldacquisitionday,oldacquisition;
  
  int dummy,truepositive,falsenegative,posatadmission,negatadmission,event,change;
  double newa,newb,newp,newlogl,p,r,x,y,q;//logp,
  double newc, newnu, newmu; /// ADDED to Martins code
  
  int thesize,donothing,index,i,ward;
  vector<double> newalpha(numberofwards),newbeta(numberofwards), newgamma(numberofwards), newnuvector(numberofwards), newmuvector(numberofwards);
  
  int acqday; /// For loop over number of patients
  double contralpha, contrbeta, contrE, contrCrossT, contrEnv, contrprev;
  
  
  // For startline output
  int out = 0; 
  /**********************************************************************
  * DEFINE PARAMETERS FOR MCMC
  **********************************************************************/
  int thinning=10,thinningscreen=10000;
  int thinningfoi=5000;
  int onlyaugmented=0; ///if 1 only augmented data are updated
  /**frequencydependence=0 means density dependence,
  *  frequencydependence=1 means frequency dependence **/
  int frequencydependece =1;
  long int numberofupdates=530000;
  
  /**********************************************************************
  // Determine prior distributions
  * ---------------------------------------------------------------------
  * Sensitivity test has as prior distribution a beta distribution with
  * parameters aa and bb.
  * Parameters endogenous route and cross transmisison (alpha and beta)
  * have as prior an exponential ditribution with parameter lambda
  **********************************************************************/
  double aa=1.,bb=1.,lambda=0.001;
  
  /* parameter w defines the probability that a new colonization
  * time is an admission */
  double w=0.2;
  
  /**********************************************************************
  * Initial standarddeviation proposal distribution for
  * alpha, beta
  * gamma(c), mu, nu, E0 (added to Martins code)
  **********************************************************************/
  double sigmaa=0.2; double sigmab=0.12;
  double sigmap=0.15;
  // Do Not know which value here to take:
  double sigmac=0.2, sigmanu = 0.2, sigmamu=0.2, sigmaE0=0.2;
  // For adaptation of acceptance ratio for gamma in Hamiltonian MCMC
  double newsigmac = 0.2, newmuc = 0.0, propratio = 1.0, newmub=0.0, newsigmab=0.2;
  /** sigmaa and sigmab are modified such that acceptance rate more or
  *  less between lowerboundacceptance and upperboundacceptance **/
  double lowerboundacceptance=0.05;
  double upperboundacceptance=0.5;
  /** updatetochangesigma determines number of updates to check whether
  * acceptancerate for a and b are reasonable **/
  int updatestochangesigma=100;
  
  /// Define seed for the random number generator
  int seed=40;
  /*******************************
  change to
  ENG eng(std::time(0));
  to have random starting of chain
  *******************************/
  /**********************************************************************
  * Initial values of MCMC chain
  * --------------------------------------------------------------------
  * alpha = a, beta = b, importation rate = f, sensitivity = phi
  * gamma = c, nu, mu, E0 (added to Martins code)
  **********************************************************************/
  parameters param;
  param.a=0.01;
  param.b=0.0;
  param.f=0.1;
  param.phi=0.8;
  param.frequencydependence=frequencydependece;
  
  param.p=0.01;
  param.c=0.01;
  param.nu=0.1;
  param.mu=0.01;
  double meanPrev = 0.15;
  param.E0 = (param.nu/param.mu)*meanPrev;
  
  
  /**********************************************************************
  * ADDED: Parameters for environmental contamination
  * --------------------------------------------------------------------
  * addingenvironment = 1 if environmental contamination should be added
  * to force of infection, 0 otherwise
  **********************************************************************/
  int addingenvironment = 0;
  int endogenous = 1, exogenous = 1;
  int bedwise = 1, importation=1;
  double newE0; double meanE;
  int day;
  
  /// Define startperiod and endperiod
  int startperiod=0; /// defined with first date set to day zero
  int endperiod=200; /// defined with first date set to day zero
  
  /// Added caccepted =0, nuaccepted, muaccepted
  /// counts accepted updates for alpha and beta
  long int aaccepted=0,baccepted=0,caccepted=0, nuaccepted=0, muaccepted=0, paccepted=0;
  
  /**********************************************************************
  * PRIOR DISTRIBUTION for
  * --------------------------------------------------------------------
  * alpha, beta
  * gamma, nu, mu, E0 (added to Martins code)
  **********************************************************************/
  boost::math::exponential_distribution<>priorbeta(lambda);
  boost::math::exponential_distribution<>prioralpha(lambda);
  
  boost::math::exponential_distribution<>priorp(lambda);
  //boost::math::uniform_distribution<>priorp(0,1.0);
  boost::math::exponential_distribution<>priorgamma(lambda);
  boost::math::exponential_distribution<>priornu(lambda);
  boost::math::exponential_distribution<>priormu(lambda); string priormustring="Exponential(0.001)";
  //boost::math::uniform_distribution<>priormu(0,1); string priormustring="U(0,1)";
  boost::math::exponential_distribution<>priorE0(lambda);
  
  /// Does each ward has its own parameters for alpha, beta or not?
  int multipleparameters=0;
  
  /**********************************************************************
  * INPUT DATA
  /*********************************************************************/
  
  /*
  string env = "Env"; string icu = "D7121"; string hosp = "Besancon";
  string resultspath="../../"+hosp+"Hospital/Results/7121/UnitD/Test/";
  string savefile = resultspath + "mytestresults" + hosp + icu + env + ".txt";
  string envfile = resultspath + "myenvresults" + hosp + icu + env + ".txt";
  string envstayfile=resultspath + "myenvstayresults" + hosp + icu + env + ".txt";
  string envdisfile=resultspath + "myenvdisresults" + hosp + icu + env + ".txt";
  string prevfile = resultspath + "myprevresults" + hosp + icu + env + ".txt";
  string Nfile = resultspath + "myNresults" + hosp + icu + env + ".txt";
  string acqfile = resultspath + "myacqresults" + hosp + icu + env + ".txt";
  string startline = resultspath + "mystartline" + hosp + icu + env + ".txt";
  string patientsfile = resultspath + "mypatients" + hosp + icu + env + ".txt";
  string colstatusfile = resultspath + "colstatus" + icu + env + ".txt";
  // Input File
  string datapath = "../../" + hosp + "Hospital/Data/";
  string admissiondata = datapath + "admDates" + hosp + icu + ".txt";
  string culturedata = datapath + "cultureresults" + hosp + icu + ".txt";
  */
  
  
  /**
  string savefile="../../BesanconHospital/Results/mytestresults7441.txt";
  //string foifile="../../BesanconHospital/Results/myfoiresults7441.txt";
  string envfile="../../BesanconHospital/Results/myenvresults7441.txt";
  string prevfile="../../BesanconHospital/Results/myprevresults7441.txt";
  string Nfile="../../BesanconHospital/Results/myNresults7441.txt";
  string startline="../../BesanconHospital/Results/mystartline7441.txt";
  string acqfile="../../BesanconHospital/Results/myacqresults7441.txt";
  // Input File
  string admissiondata="../../BesanconHospital/Data/admDates7441.txt";
  string culturedata="../../BesanconHospital/Data/cultureresults7441.txt";
  **/
  
  /**********************************************************************
  * FOR RUNNING ON THE HPC CLUSTER
  /*********************************************************************/
  
  /**
  string resultspath="/home/julius_id/tpham/UtrechtHospital/Results/HPC07122017/";
  string env = "Env"; string icu = "1";
  string savefile = resultspath + "mytestresultsUtrecht" + icu + env + ".txt";
  string envfile = resultspath + "myenvresultsUtrecht" + icu + env + ".txt";
  string prevfile = resultspath + "myprevresultsUtrecht" + icu + env + ".txt";
  string Nfile = resultspath + "myNresultsUtrecht" + icu + env + ".txt";
  string acqfile = resultspath + "myacqresultsUtrecht" + icu + env + ".txt";
  string startline = resultspath + "mystartlineUtrecht" + icu + env + ".txt";
  string patientsfile = resultspath + "mypatientsUtrecht" + icu + env + ".txt";
  // Input File
  string datapath = "../../UtrechtHospital/Data/";
  string admissiondata = datapath + "admDatesUtrecht" + icu + ".txt";
  string culturedata = datapath + "cultureresultsUtrecht" + icu + ".txt";
  **/
  
  
  string env = ""; string icu = "7121BeforeAfter"; string hosp = "Besancon";
  string resultspath="/home/julius_id/tpham/"+hosp+"Hospital/Results/7121/BedwiseExp/";
  string savefile = resultspath + "mytestresults" + hosp + icu + env + ".txt";
  string meanParam = resultspath + "meanParam" + hosp + icu + env + ".txt";
  string foifile = resultspath + "myfoiresults" + hosp + icu + env + ".txt";
  string envfile = resultspath + "myenvresults" + hosp + icu + env + ".txt";
  string envstayfile=resultspath + "myenvstayresults" + hosp + icu + env + ".txt";
  string envdisfile=resultspath + "myenvdisresults" + hosp + icu + env + ".txt";
  string prevfile = resultspath + "myprevresults" + hosp + icu + env + ".txt";
  string impfile = resultspath + "myimpresults" + hosp + icu + env + ".txt";
  string Nfile = resultspath + "myNresults" + hosp + icu + env + ".txt";
  string acqfile = resultspath + "myacqresults" + hosp + icu + env + ".txt";
  string startline = resultspath + "mystartline" + hosp + icu + env + ".txt";
  string patientsfile = resultspath + "mypatients" + hosp + icu + env + ".txt";
  string colstatusfile = resultspath + "colstatus" + hosp + icu + env + ".txt";
  // Input File
  string datapath = "/home/julius_id/tpham/" + hosp + "Hospital/Data/";
  string admissiondata = datapath + "admDates" + hosp + icu + ".txt";
  string culturedata = datapath + "cultureresults" + hosp + icu + ".txt";
  string previouspatdata = datapath + "previousPat" + hosp + icu + ".txt";
  
  
  /*
  string env = ""; string icu = "Mixed_2"; string hosp = "";
  string resultspath="/home/julius_id/tpham/ArtifSim/Results/Mixed/All3/PriorExp0.001/";
  string savefile = resultspath + "mytestresults" + icu + env + ".txt";
  string meanParam = resultspath + "meanParam" + icu + env + ".txt";
  string foifile = resultspath + "myfoiresults" + icu + env + ".txt";
  string envfile = resultspath + "myenvresults" + icu + env + ".txt";
  string envstayfile=resultspath + "myenvstayresults" + icu + env + ".txt";
  string envdisfile=resultspath + "myenvdisresults" + icu + env + ".txt";
  string prevfile = resultspath + "myprevresults" + icu + env + ".txt";
  string Nfile = resultspath + "myNresults" + icu + env + ".txt";
  string acqfile = resultspath + "myacqresults" + icu + env + ".txt";
  string startline = resultspath + "mystartline" + icu + env + ".txt";
  string patientsfile = resultspath + "mypatients" + icu + env + ".txt";
  string colstatusfile = resultspath + "colstatus" + icu + env + ".txt";
  // Input File
  string datapath = "/home/julius_id/tpham/ArtifSim/Data/Mixed/";
  string admissiondata = datapath + "admDates" + hosp + icu + ".txt";
  string culturedata = datapath + "cultureresults" + hosp + icu + ".txt";
  */
  
  /*
  string env = ""; string icu = "17"; string hosp = "Bedwise_";
  string resultspath="/home/julius_id/tpham/ArtifSim/Results/Bedwise/" + icu + "/BedwiseExp/";
  string savefile = resultspath + "mytestresults" + hosp + icu + env + ".txt";
  string meanParam = resultspath + "meanParam" + hosp + icu + env + ".txt";
  string foifile = resultspath + "myfoiresults" + hosp + icu + env + ".txt";
  string envfile = resultspath + "myenvresults" + hosp + icu + env + ".txt";
  string envstayfile=resultspath + "myenvstayresults" + hosp + icu + env + ".txt";
  string envdisfile=resultspath + "myenvdisresults" + hosp + icu + env + ".txt";
  string prevfile = resultspath + "myprevresults" + hosp + icu + env + ".txt";
  string impfile = resultspath + "myimpresults" + hosp + icu + env + ".txt";
  string Nfile = resultspath + "myNresults" + hosp + icu + env + ".txt";
  string acqfile = resultspath + "myacqresults" + hosp + icu + env + ".txt";
  string startline = resultspath + "mystartline" + hosp + icu + env + ".txt";
  string patientsfile = resultspath + "mypatients" + hosp + icu + env + ".txt";
  string colstatusfile = resultspath + "colstatus" + hosp + icu + env + ".txt";
  // Input File
  string datapath = "/home/julius_id/tpham/ArtifSim/Data/Bedwise/";
  string admissiondata = datapath + "admDates" + hosp + icu + ".txt";
  string culturedata = datapath + "cultureresults" + hosp + icu + ".txt";
  string previouspatdata = datapath + "previousPat" + hosp + icu + ".txt";
  */
  
  
  /**********************************************************************
  * PROCESSING INPUT DATA
  /*********************************************************************/
  
  /**********************************************************************
  * culturedata is of form: patid, date, result
  * admissiondata is of form: patid, date, roomnumber
  * Assume rooms are labeled with integers
  * Dates are in integers
  * Patient ID are integers
  * day with the label i runs from day i: 12.00 -> day i+1, 11.59
  **********************************************************************/
  cout << "Admission File: " << admissiondata << endl; cout << "Culture File: " << culturedata << endl;
  cout << "Number of updates: " << numberofupdates << endl;
  cout << "Mean prevalence = " << meanPrev << endl;
  
  vector< vector<int> > admissionarray,culturearray,bedarray;
  READCSVFILETO2DARRAYOFINT(admissiondata,admissionarray);
  READCSVFILETO2DARRAYOFINT(culturedata,culturearray);
  
  // Contains the NEXT patient
  vector<int> previouspatarray(numberofpatients);
  if(bedwise==1){// Assume that the vector was already sorted by patient number!
    ifstream previousfile(previouspatdata.c_str());
    int tmp;  
    while(previousfile >> tmp){
      if(tmp!=-1) tmp-=1;
      previouspatarray.push_back(tmp);
    }
  }
  
  //printarray(admissionarray, admissiondata); cout <<"BB"<<endl;
  //printarray(culturearray,culturedata); cout <<"AA"<<endl;
  
  cout << "**************************************"<<endl;
  cout << "Processing data for MCMC simulation" << endl;
  cout << "**************************************"<<endl;
  /// OUTPUT: Print which routes are considered for simulation
  cout << "Processing data or MCMC simulation. " << endl;
  if(addingenvironment){
    cout << "Estimating parameters for background, cross-transmission AND environmental contamination." << endl;
    cout << "Prior mu: " << priormustring << endl;
  } else{
    if(bedwise){
      cout << "Estimating parameters for background, cross-transmission and impact of prior bed occupants." << endl;
    }else cout << "Estimating parameters for background and cross-transmission." << endl;
  }
  
  
  cout << "Relabel wards" << endl;
  /// relabelwards: Wards are relabeled with integers: 0, 1,2,3,...,N-1
  relabelwards(admissionarray,&numberofwards);
  cout << "numberofpatients="<<numberofpatients<<", numberofwards="<<numberofwards<<", maxdate="<<maxdate<<"\n";
  /// printarray(admissionarray, admissiondata);
  /// relabelpatientID: Patient ID are relabeled: 0, 1, 2, 3, ...., M-1
  relabelpatientID(admissionarray,culturearray,&numberofpatients);
  cout << "numberofpatients="<<numberofpatients<<", numberofwards="<<numberofwards<<", maxdate="<<maxdate<<"\n";
  /// relabeldates: Dates are relabeled 0, 1, 2, ...., T-1
  relabeldates(admissionarray,culturearray,&maxdate);
  
  cout << "numberofpatients="<<numberofpatients<<", numberofwards="<<numberofwards<<", maxdate="<<maxdate<<"\n";
  ///printarray(admissionarray, admissiondata);printarray(culturearray,culturedata);
  cout <<"length culturearray="<<culturearray.size()<<endl;
  cout <<"length admissionarray="<<admissionarray.size()<<endl;
  
  /** factoraugmented: Number of augmentation steps
  * changed later on if numberofiteration >changeaugmented **/
  factoraugmented=numberofpatients*10; // factoraugmented = 5;
  //factoraugmented = 20;
  changeaugmented=2;
  newfactoraugmented=20; 
  // newfactoraugmented=1+numberofpatients/4;
  cout <<"newfactoraugmented="<<newfactoraugmented<<endl;
  
  /// Initialize parameter vectors now since unit is defined
  for(unit=0;unit<numberofwards;unit++)
  {
    (param.alphavector).push_back(param.a);
    (param.betavector).push_back(param.b);
    
    (param.gammavector).push_back(param.c);
    (param.nuvector).push_back(param.nu);
    (param.muvector).push_back(param.mu);
  }
  ///cout <<"lengte alpha="<<param.alphavector.size()<<endl;
  ///cout <<"lengte beta="<<param.betavector.size()<<endl;
  
  
  vector<vector<int> > N(numberofwards, vector<int>(maxdate,0));
  /** whereabouts:
  * For each patient it includes a vector with days on which the patient
  * is present. Note whereabouts of patient i startperiod at
  * admissionday of i **/
  vector<vector<int> > whereabouts(numberofpatients);
  /** colstatus:
  * note colstatus of patient i start at admissionday of i
  * -1, already colonized,
  * 0,  no acquisition during that day,
  * 1,  acquisition during that day. **/
  vector<vector<int> > colstatus(numberofpatients);
  vector<int> admissionday(numberofpatients,maxdate+10),dischargeday(numberofpatients,-1);
  std::sort (admissionarray.begin(), admissionarray.end(), sortFunc1);
  numberofpatientsperward(N,admissionarray,whereabouts,admissionday,dischargeday);
  /** present:
  * per ward, per day a list of which patients are present in the unit.**/
  vector<vector<vector<int> > > present (numberofwards,vector<vector<int> >(maxdate,vector <int>()));
  whopresentwhen(present, admissionarray);
  
  vector<int> listneverpos,listdefpos,listtemppos;
  vector<vector<int> > cultureperpatient(numberofpatients);
  vector<int> admissionstate(numberofpatients, 1);/* 0 uncolonized upon admission, 1 colonized upon admission */
  vector<int> acquisition(numberofpatients,0);/* 0=no acquisition, 1 acquisition */
  vector<int> acquisitionday(numberofpatients,-1);/* acquisition is during that day (at end of that specific day (11.59 next day)) */
  vector<vector<int> > I(numberofwards, vector<int>(maxdate,0));
  vector<vector<int> > numberofacquisitions(numberofwards, vector<int>(maxdate,0));
  vector<int> daysatrisk(numberofwards);
  vector<int> totalacquisitions(numberofwards);
  
  
  cultureresultsperpatient(culturearray,cultureperpatient);
  ///cout << "Culture per patient:" << endl;
  ///printarray(cultureperpatient,"");
  poscultureposatadmission(culturearray,cultureperpatient,admissionstate);
  ///printarray(colstatus,"");
  numberofinfectiouspatientsperward(I,whereabouts,colstatus,admissionday,acquisitionday,admissionstate,acquisition,numberofacquisitions);
  ///printarray(colstatus,"");
  int patient,indexday;
  defineneverposdefpostemppos(listneverpos,listdefpos,listtemppos,admissionstate,acquisition,culturearray,cultureperpatient);
  cout<<"listneverpos,listdefpos,listtemppos="<<listneverpos.size()<<", "<<listdefpos.size()<<", "<<listtemppos.size()<<endl;
  
  // cout << "listneverpos:" <<endl;
  // for(int k=0;k<listneverpos.size();k++){
  //   cout << listneverpos.at(k) << ",";
  // }
  // cout << endl;
  // 
  // cout << "listdefpos:" <<endl;
  // for(int k=0;k<listdefpos.size();k++){
  //   cout << listdefpos.at(k) << ",";
  // }
  // cout << endl;
  // 
  // cout << "listtemppos:" <<endl;
  // for(int k=0;k<listtemppos.size();k++){
  //   cout << listtemppos.at(k) << ",";
  // }
  // cout << endl;
  
  //Generate Random number for exponential prior distribution for environmental contamination
  //boost::mt19937 gen;
  
  // Initialize array with whether previous bed occupant was colonized 
  // cout<<"Print previouspatarray: " << endl;
  // for(int tt=0; tt<previouspatarray.size(); tt++){
  //   cout << previouspatarray.at(tt) << ", " << endl;
  // }
  vector<int> previousCol(numberofpatients,0);
  if(bedwise==1) previouspatientstatus(previouspatarray,previousCol,acquisition,admissionstate);
  
  
  /**********************************************************************
  * ADDED: Initilaization of environmental contamination
  * --------------------------------------------------------------------
  * E: Vector for environmental contamination for each ward and each day
  * E_stay: Vector for env.cont. for each ward and each day contributed
  * by present patients
  * E_dis: Vector for env.cont. for each ward and each day contributed
  * by discharged patients
  * *******************************************************************/
  vector<vector<double> > E(numberofwards, vector<double>(maxdate+1,0));
  vector<vector<double> > E_stay(numberofwards, vector<double>(maxdate+1,0));
  vector<vector<double> > E_dis(numberofwards, vector<double>(maxdate+1,0));
  if(addingenvironment == 1){
    cout << "E0 = " << param.E0 << "\n";
    ///updateenvironment(param,I,N,E,numberofwards,maxdate); // Calculate environmental load with initial number of infectious patients, inital parameters and inital bacterial load E0
    updateenvironment2(param,I,N,admissionstate,acquisition,acquisitionday,admissionday,dischargeday,E,E_stay,E_dis,numberofwards,maxdate);
    cout << "Initial update of environment was successful." << endl;
  }
  
  /// Update endperiod (added to Martins code)
  endperiod = maxdate;
  cout <<"voor logl: maxdate="<<maxdate<<" startperiod="<<startperiod<<" endperiod="<<endperiod<<endl;
  
  /// First logllikelihood computation
  double logl=loglikelihood(acquisitionday,acquisition,culturearray,N,I,numberofacquisitions,admissionstate,param,startperiod,endperiod,admissionday,multipleparameters,E,addingenvironment,bedwise,previousCol,present);
  cout << "logl = " << logl << "\n";
  cout <<"na logl"<<endl;
  
  // Why again?
  numberofinfectiouspatients2(I,whereabouts,colstatus,admissionday,acquisitionday,admissionstate,acquisition,numberofacquisitions,maxdate,numberofwards);
  
  /** possibleacquisitiondays:
  * determine possible days of acquisition per patient */
  vector<int> possibleacquistiondays(numberofpatients);
  determinepossibleacquisitiondays(possibleacquistiondays,culturearray,whereabouts,admissionday);
  
  typedef boost::mt19937  ENG;    /* Mersenne Twister: NOT SEEDED YET */
  ENG eng(seed);
  
  parameters paramnew;
  for(unit=0;unit<numberofwards;unit++)
  {
    (paramnew.alphavector).push_back(param.a);
    (paramnew.betavector).push_back(param.b);
    
    (paramnew.gammavector).push_back(param.c);
    (paramnew.nuvector).push_back(param.nu);
    (paramnew.muvector).push_back(param.mu);
  }
  
  /**********************************************************************
  * Definition of file streams for saving output
  *********************************************************************/
  ofstream myfile;
  myfile.open(savefile.c_str());
  ofstream mymeanfile;
  mymeanfile.open(meanParam.c_str());
  
  ofstream mystartfile;
  mystartfile.open(startline.c_str());
  
  ofstream myfoifile;
  myfoifile.open(foifile.c_str());
  ofstream myprevfile;
  myprevfile.open(prevfile.c_str());
  ofstream myNfile;
  myNfile.open(Nfile.c_str());
  ofstream myacqfile;
  myacqfile.open(acqfile.c_str());
  ofstream myimpfile;
  myimpfile.open(impfile.c_str());
  
  ofstream myenvfile;
  myenvfile.open(envfile.c_str());
  ofstream myenvstayfile;
  myenvstayfile.open(envstayfile.c_str());
  ofstream myenvdisfile;
  myenvdisfile.open(envdisfile.c_str());
  
  
  ofstream mypatients;
  mypatients.open(patientsfile.c_str());
  ofstream colstat;
  colstat.open(colstatusfile.c_str());
  
  myfile <<"i,f,phi,logl,nImp,meanPrev,alpha,beta,p,eps,nu,mu,E0,TP,FN,aacceptrate,bacceptrate,pacceptrate,cacceptrate,muacceptrate,sigmaa,sigmab,sigmap,sigmac,sigmamu,acq,daysatr,contralpha,contrbeta,contrprev,contrE,contrCrossT,contrEnv"<<"\n";
  /// Computing true positives and false negatives
  truepositiveandfalsenegative(&truepositive,&falsenegative,culturearray,admissionday,colstatus,startperiod,endperiod);
  cout << "True positives="<<truepositive<<endl;
  cout << "False negatives="<<falsenegative<<endl;
  if(addingenvironment == 1){
    cout << "E0 = " << param.E0 << "\n";
    ///updateenvironment(param,I,N,E,numberofwards,maxdate); // Calculate environmental load with initial number of infectious patients, inital parameters and inital bacterial load E0
    updateenvironment2(param,I,N,admissionstate,acquisition,acquisitionday,admissionday,dischargeday,E,E_stay,E_dis,numberofwards,maxdate);
    cout << "Initial update of environment was successful." << endl;
  }
  
  // Print the patients that are present at day i
  int pat;
  // for(int d=startperiod;d<endperiod;d++){
  //   cout << "Day " << d << ": ";
  //   for(int pt=0;pt<present.at(0).at(d).size();pt++){
  //     pat = present.at(0).at(d).at(pt);
  //     cout << pat << ",";
  //   }
  //   cout << "\n";
  // }
  
  // Print the colonization status of the prior bed occupants for each patient:
  mypatients << "I per day:" << endl;
  for(int d=startperiod;d<endperiod-1;d++){
    mypatients << I.at(0).at(d)<< ",";
  }
  mypatients << I.at(0).at(endperiod-1)<< "\n";
  
  mypatients << "N per day:" << endl;
  for(int d=startperiod;d<endperiod-1;d++){
    mypatients << N.at(0).at(d)<< ",";
  }
  mypatients << N.at(0).at(endperiod-1)<< "\n";
  
  //int pat;
  mypatients << "Previous colonized bed opccupant?" << endl;
  for(pat=0;pat<numberofpatients-1;pat++){
    mypatients << previousCol.at(pat) << ",";
  }
  mypatients << previousCol.at(numberofpatients-1) << "\n";
  
  mypatients << "Importation?" << endl;
  for(pat=0;pat<numberofpatients-1;pat++){
    mypatients << admissionstate.at(pat) << ",";
  }
  mypatients << admissionstate.at(numberofpatients-1)<<"\n";
  
  mypatients << "Acquisition time:" << endl;
  for(pat=0;pat<numberofpatients-1;pat++){
    if(acquisition.at(pat)==1) mypatients << admissionday.at(pat) + acquisitionday.at(pat) <<",";
    else mypatients << "-1,";
  }
  if(acquisition.at(numberofpatients-1)==1) mypatients<<admissionday.at(numberofpatients-1)+acquisitionday.at(numberofpatients-1)<<"\n";
  else mypatients << "-1\n";
  
  mypatients << "Acquisition time End." << endl;
  
  // Print force of infection per day
  // cout << "Force of infection per day:" << endl;
  // double foip1, foip2;
  // for(int d=startperiod;d<endperiod;d++){
  //   foip1 = 0.0; foip2=0.0;
  //   cout << "Day " << d << ": "<<endl;
  //   foi1(param,I.at(0).at(d),N.at(0).at(d),E.at(0).at(d),addingenvironment,foip1,bedwise,1);
  //   foi1(param,I.at(0).at(d),N.at(0).at(d),E.at(0).at(d),addingenvironment,foip2,bedwise,0);
  //   cout << "FOI, previous=1: " << foip1 << ", FOI, previous=0: " << foip2 << endl;
  // }
  // 
  //   cout << "I per day:" << endl;
  //   for(int d=startperiod;d<endperiod;d++){
  //     cout << I.at(0).at(d)<< ",";
  //   }
  //   cout << endl;
  //   cout << "N per day:" << endl;
  //   for(int d=startperiod;d<endperiod;d++){
  //     cout << N.at(0).at(d)<< ",";
  //   }
  //   cout << endl;
  
  
  cout << "**************************"<<endl;
  cout <<"start MCMC\n";
  cout << "**************************"<<endl;
  
  
  /**********************************************************************
  *  start the MCMC algorithm.
  *********************************************************************/
  i=0;
  /// this part ensures that the acceptance rate of alpha and beta is approximately between 0.05 and 0.5
  while(i<numberofupdates)
  {
    if(onlyaugmented!=1)
    {
      change=0;
      
      if(i==updatestochangesigma)
      {
        if(endogenous==1){
          if((double)aaccepted/(double)updatestochangesigma<lowerboundacceptance)
          {
            sigmaa=sigmaa/samplefromrealuniform(1.5, 2.5, &eng);change=1;
            ///cout<<"Change a. Lower bound. \n";
          }
          if((double)aaccepted/(double)updatestochangesigma>upperboundacceptance)
          {
            sigmaa=sigmaa*samplefromrealuniform(1.5, 2.5, &eng);change=1;
            ///cout<<"Change a. Upper bound. \n";
          }
        }
        
        if(exogenous==1){
          if((double)baccepted/(double)updatestochangesigma<lowerboundacceptance)
          {
            sigmab=sigmab/samplefromrealuniform(1.5, 2.5, &eng);change=1;
            ///cout<<"Change b. Lower bound.\n";
          }
          if((double)baccepted/(double)updatestochangesigma>upperboundacceptance)
          {
            sigmab=sigmab*samplefromrealuniform(1.5, 2.5, &eng);change=1;
            ///cout<<"Change b. Upper bound. \n";
          }
        }
        
        if(bedwise==1){
          if((double)paccepted/(double)updatestochangesigma<lowerboundacceptance)
          {
            sigmap=sigmap/samplefromrealuniform(1.5, 2.5, &eng);change=1;
            ///cout<<"Change p. Lower bound.\n";
          }
          if((double)paccepted/(double)updatestochangesigma>upperboundacceptance)
          {
            sigmap=sigmap*samplefromrealuniform(1.5, 2.5, &eng);change=1;
            ///cout<<"Change p. Upper bound. \n";
          }
        }
        
        if(addingenvironment == 1){
          /// ADDED
          if((double)caccepted/(double)updatestochangesigma<lowerboundacceptance)
          {
            sigmac=sigmac/samplefromrealuniform(1.5, 2.5, &eng);change=1;
            ///cout<<"Change c. Lower bound. \n";
          }
          if((double)caccepted/(double)updatestochangesigma>upperboundacceptance)
          {
            sigmac=sigmac*samplefromrealuniform(1.5, 2.5, &eng);change=1;
            ///cout<<"Change c. Upper bound. \n";
          }
          if((double)muaccepted/(double)updatestochangesigma<lowerboundacceptance)
          {
            sigmamu=sigmamu/samplefromrealuniform(1.5, 2.5, &eng);change=1;
            ///cout<<"Change mu. Lower bound. \n";
          }
          if((double)muaccepted/(double)updatestochangesigma>upperboundacceptance)
          {
            sigmamu=sigmamu*samplefromrealuniform(1.5, 2.5, &eng);change=1;
            ///cout<<"Change mu. Upper bound.\n";
          }
          
        }
        /// END ADDED
        if(change==1){i=0;aaccepted=0;baccepted=0;caccepted=0;nuaccepted=0;muaccepted=0;paccepted=0;nchange+=1;}
      }
      if(i==0)
      {
        inew=updatestochangesigma;
      }else
      {
        inew=i;
      }
      
      if(i==changeaugmented)factoraugmented=newfactoraugmented;
      
      /**update param.a;*****/
      if(i>0)
      {
        ///cout << "i>0" <<endl;
        if(multipleparameters==0)
        {
          if(addingenvironment == 0){
            
            /**update param.a ***/
            
            paramnew.a=param.a;
            paramnew.b=param.b;
            /// ADDED
            paramnew.p=param.p;
            paramnew.c=param.c;
            paramnew.nu=param.nu;
            paramnew.mu=param.mu;
            paramnew.E0=param.E0;
            
            paramnew.frequencydependence= param.frequencydependence;
            paramnew.phi=param.phi;
            paramnew.f=param.f;
            
            if(bedwise){
              newp=abs(param.p+samplefromnormal(0,sigmap,&eng));
              paramnew.p=newp;
              // ADDED Environmental contamination E
              newlogl=loglikelihood(acquisitionday,acquisition,culturearray,N,I,numberofacquisitions,admissionstate,paramnew,startperiod,endperiod,admissionday,multipleparameters,E,addingenvironment,bedwise,previousCol,present);
              //cout <<"param.p="<<param.p<<", logl="<<logl<<", newparam.p="<<paramnew.p<<", newlogl="<<newlogl<<endl;
              //logp=min(0.,newlogl+log(pdf(priorbeta,newb))-logl-log(pdf(priorbeta,param.b)));
              //p=exp(logp);
              
              p=min(1.,exp(newlogl-logl)*pdf(priorp,newp)/pdf(priorp,param.p));
              r=samplefromrealuniform(0, 1, &eng);
              if(r<p){param.p=newp;logl=newlogl;paccepted+=1;}
              
              // mypatients << endl;
              // mypatients << "i=" << i << endl;
              // mypatients <<"param.p="<<param.p<<", logl="<<logl<<", newparam.p="<<paramnew.p<<", newlogl="<<newlogl<<endl;
              // mypatients << "r=" << r << ", p="<< p << endl;
              // truepositiveandfalsenegative(&truepositive,&falsenegative,culturearray,admissionday,colstatus,startperiod,endperiod);
              // mypatients << "TP="<< truepositive << "; FN=" << falsenegative << endl;
              // mypatients << "alpha<-"<<param.a << "; beta<-" << param.b << "; p<-" << param.p << "; f<-" << param.f << "; phi<-" << param.phi << endl;
              // mypatients << "I per day:" << endl;
              // mypatients << "inf<-c(";
              // for(int d=startperiod;d<endperiod-1;d++){
              //   mypatients << I.at(0).at(d)<< ",";
              // }
              // mypatients << I.at(0).at(endperiod-1)<< ")\n";
              // mypatients << "N per day:" << endl;
              // mypatients << "npat<-c(";
              // for(int d=startperiod;d<endperiod;d++){
              //   mypatients << N.at(0).at(d)<< ",";
              // }
              // mypatients<<N.at(0).at(endperiod-1)<< ")\n";
              // 
              // mypatients << "Previous colonized bed opccupant?" << endl;
              // mypatients << "previousCol<-c(";
              // for(pat=0;pat<numberofpatients-1;pat++){
              //   mypatients << previousCol.at(pat) << ",";
              // }
              // mypatients << previousCol.at(numberofpatients-1) << ")\n";
              // 
              // mypatients << "Importation?" << endl;
              // mypatients << "imp<-c(";
              // for(pat=0;pat<numberofpatients-1;pat++){
              //   mypatients << admissionstate.at(pat) << ",";
              // }
              // mypatients << admissionstate.at(numberofpatients-1) << ")\n";
              // 
              // mypatients << "Acquisition time:" << endl;
              // mypatients << "colTime<-c(";
              // for(pat=0;pat<numberofpatients-1;pat++){
              //   if(acquisition.at(pat)==1) mypatients << admissionday.at(pat) + acquisitionday.at(pat) <<",";
              //   else mypatients << "-1,";
              // }
              // if(acquisition.at(numberofpatients-1)==1) mypatients<<admissionday.at(numberofpatients-1)+acquisitionday.at(numberofpatients-1)<<")\n";
              // else mypatients << "-1)\n";
              // 
              // mypatients << "Numberofacquisitions per day:" << endl;
              // mypatients << "numAcq<-c(";
              // for(int d=startperiod;d<endperiod-1;d++){
              //   mypatients << numberofacquisitions.at(0).at(d)<< ",";
              // }
              // mypatients << numberofacquisitions.at(0).at(endperiod-1)<< ")\n";
            }
            paramnew.p=param.p;
            
            if(endogenous){
              newa=abs(param.a+samplefromnormal(0,sigmaa,&eng));
              paramnew.a=newa;
              
              /// ADDED Environmental contamination E
              newlogl=loglikelihood(acquisitionday,acquisition,culturearray,N,I,numberofacquisitions,admissionstate,paramnew,startperiod,endperiod,admissionday,multipleparameters,E,addingenvironment,bedwise,previousCol,present);
              //cout <<"param.a="<<param.a<<", logl="<<logl<<", newparam.a="<<paramnew.a<<", newlogl="<<newlogl<<endl;
              
              //double logp=min(0.,newlogl+log(pdf(prioralpha,newa))-logl-log(pdf(prioralpha,param.a)));
              //p=exp(logp);
              p=min(1.,exp(newlogl-logl)*pdf(prioralpha,newa)/pdf(prioralpha,param.a));
              r=samplefromrealuniform(0, 1, &eng);
              if(r<p){param.a=newa;logl=newlogl;aaccepted+=1;}
            }
            /** End update a **/
            
            
            /**update param.b ***/
            
            paramnew.a=param.a;
            if(exogenous){
              newb=abs(param.b+samplefromnormal(0,sigmab,&eng));
              paramnew.b=newb;
              // ADDED Environmental contamination E
              newlogl=loglikelihood(acquisitionday,acquisition,culturearray,N,I,numberofacquisitions,admissionstate,paramnew,startperiod,endperiod,admissionday,multipleparameters,E,addingenvironment,bedwise,previousCol,present);
              //cout <<"param.b="<<param.b<<", logl="<<logl<<", newparam.b="<<paramnew.b<<", newlogl="<<newlogl<<endl;
              //logp=min(0.,newlogl+log(pdf(priorbeta,newb))-logl-log(pdf(priorbeta,param.b)));
              //p=exp(logp);
              
              p=min(1.,exp(newlogl-logl)*pdf(priorbeta,newb)/pdf(priorbeta,param.b));
              r=samplefromrealuniform(0, 1, &eng);
              if(r<p){param.b=newb;logl=newlogl;baccepted+=1;}
            }
            /** End update b **/
            
            /**Update param.p ***/
            paramnew.b=param.b;
            
            /** End update p **/
          }
          
          if(addingenvironment == 1){
            
            /**update param.a ***/
            paramnew.b=param.b;
            /// ADDED
            paramnew.p=param.p;
            paramnew.c=param.c;
            paramnew.nu=param.nu;
            paramnew.mu=param.mu;
            paramnew.E0=param.E0;
            
            paramnew.frequencydependence= param.frequencydependence;
            paramnew.phi=param.phi;
            paramnew.f=param.f;
            
            if(endogenous){
              newa=abs(param.a+samplefromnormal(0,sigmaa,&eng));
              paramnew.a=newa;
              
              /// ADDED Environmental contamination E
              newlogl=loglikelihood(acquisitionday,acquisition,culturearray,N,I,numberofacquisitions,admissionstate,paramnew,startperiod,endperiod,admissionday,multipleparameters,E,addingenvironment,bedwise,previousCol,present);
              
              //double logp=min(0.,newlogl+log(pdf(prioralpha,newa))-logl-log(pdf(prioralpha,param.a)));
              //p=exp(logp);
              p=min(1.,exp(newlogl-logl)*pdf(prioralpha,newa)/pdf(prioralpha,param.a));
              r=samplefromrealuniform(0, 1, &eng);
              if(r<p){param.a=newa;logl=newlogl;aaccepted+=1;}
            }
            /** End update a **/
            
            
            /**update param.b ***/
            paramnew.a=param.a;
            if(exogenous){
              newb=abs(param.b + samplefromnormal(0,sigmab,&eng));
              paramnew.b=newb;
              // ADDED Environmental contamination E
              newlogl=loglikelihood(acquisitionday,acquisition,culturearray,N,I,numberofacquisitions,admissionstate,paramnew,startperiod,endperiod,admissionday,multipleparameters,E,addingenvironment,bedwise,previousCol,present);
              //cout <<"param.b="<<param.b<<", logl="<<logl<<", newparam.b="<<paramnew.b<<", newlogl="<<newlogl<<endl;
              //logp=min(0.,newlogl+log(pdf(priorbeta,newb))-logl-log(pdf(priorbeta,param.b)));
              //p=exp(logp);
              
              p=min(1.,exp(newlogl-logl)*pdf(priorbeta,newb)/pdf(priorbeta,param.b));
              r=samplefromrealuniform(0, 1, &eng);
              if(r<p){param.b=newb;logl=newlogl;baccepted+=1;}
            }
            /** End update b **/
            
            /**Update param.p ***/
            paramnew.b=param.b;
            if(bedwise){
              newp=abs(param.p+samplefromnormal(0,sigmap,&eng));
              paramnew.p=newp;
              // ADDED Environmental contamination E
              newlogl=loglikelihood(acquisitionday,acquisition,culturearray,N,I,numberofacquisitions,admissionstate,paramnew,startperiod,endperiod,admissionday,multipleparameters,E,addingenvironment,bedwise,previousCol,present);
              //cout <<"param.b="<<param.b<<", logl="<<logl<<", newparam.b="<<paramnew.b<<", newlogl="<<newlogl<<endl;
              //logp=min(0.,newlogl+log(pdf(priorbeta,newb))-logl-log(pdf(priorbeta,param.b)));
              //p=exp(logp);
              
              p=min(1.,exp(newlogl-logl)*pdf(priorp,newp)/pdf(priorp,param.p));
              r=samplefromrealuniform(0, 1, &eng);
              if(r<p){param.p=newp;logl=newlogl;paccepted+=1;}
            }
            /** End update p **/
            
            /**update param.nu **/
            /// param.nu is fixed throughout the simulation
            /** End update nu **/
            
            /**update param.mu ***/
            paramnew.p = param.p;
            paramnew.nu = param.nu;
            
            newmu=abs(param.mu + samplefromnormal(0,sigmamu,&eng));
            paramnew.mu=newmu;
            
            // ADDED
            //if(addingenvironment ==1) updateenvironment(paramnew,I,N,E,numberofwards,maxdate);
            if(addingenvironment ==1) updateenvironment2(paramnew,I,N,admissionstate,acquisition,acquisitionday,admissionday,dischargeday,E,E_stay,E_dis,numberofwards,maxdate);
            
            // ADDED Environmental contamination E
            newlogl=loglikelihood(acquisitionday,acquisition,culturearray,N,I,numberofacquisitions,admissionstate,paramnew,startperiod,endperiod,admissionday,multipleparameters,E,addingenvironment,bedwise,previousCol,present);
            //cout <<"param.mu="<<param.mu<<", logl="<<logl<<", newparam.mu="<<paramnew.mu<<", newlogl="<<newlogl<<endl;
            p=min(1.,exp(newlogl-logl)*pdf(priormu,newmu)/pdf(priormu,param.mu));
            r=samplefromrealuniform(0, 1, &eng);
            if(r<p){param.mu=newmu;logl=newlogl;muaccepted+=1;}
            
            /** End update mu **/
            
            //ADDED
            /**update param.c ***/
            paramnew.mu = param.mu;
            newc=abs(param.c+samplefromnormal(0,sigmac,&eng));
            paramnew.c=newc;
            // ADDED Environmental contamination E
            newlogl=loglikelihood(acquisitionday,acquisition,culturearray,N,I,numberofacquisitions,admissionstate,paramnew,startperiod,endperiod,admissionday,multipleparameters,E,addingenvironment,bedwise,previousCol,present);
            //cout <<"param.c="<<param.c<<", logl="<<logl<<", newparam.c="<<paramnew.c<<", newlogl="<<newlogl<<endl;
            p=min(1.,exp(newlogl-logl)*pdf(priorgamma,newc)/pdf(priorgamma,param.c));
            r=samplefromrealuniform(0, 1, &eng);
            if(r<p){param.c=newc;logl=newlogl;caccepted+=1;}
            
            /** End update c **/
            
            /**update E0 **/
            param.E0=(paramnew.nu/paramnew.mu)*meanPrev;
            
            /** End update E0 **/
          }
          // END ADDED
          
        }//end single parameter
        else//multiple parameters
        {
          /** start update a **/
          
          for(unit=0;unit<numberofwards;unit++)
          {
            newalpha.at(unit)=abs(param.alphavector.at(unit)+samplefromnormal(0,sigmaa,&eng));
            
            for(ward=0;ward<numberofwards;ward++)
            {
              paramnew.alphavector.at(ward)=param.alphavector.at(ward);
              paramnew.betavector.at(ward)=param.betavector.at(ward);
              paramnew.gammavector.at(ward)=param.gammavector.at(ward);
              paramnew.nuvector.at(ward)=param.nuvector.at(ward);
              paramnew.muvector.at(ward)=param.muvector.at(ward);
            }
            paramnew.alphavector.at(unit)=newalpha.at(unit);
            //cout <<"unita="<<unit<<endl;
            //paramnew.b=param.b;
            paramnew.phi=param.phi;
            paramnew.frequencydependence= param.frequencydependence;
            paramnew.f=param.f;
            // ADDED Environmental contamination E
            newlogl=loglikelihood(acquisitionday,acquisition,culturearray,N,I,numberofacquisitions,admissionstate,paramnew,startperiod,endperiod,admissionday,multipleparameters,E,addingenvironment,bedwise,previousCol,present);
            //cout <<"param.a="<<param.a<<", logl="<<logl<<", newparam.a="<<paramnew.a<<", newlogl="<<newlogl<<endl;
            //logp=min(0.,newlogl+log(pdf(prioralpha,newa))-logl-log(pdf(prioralpha,param.a)));
            //p=exp(logp);
            p=min(1.,exp(newlogl-logl)*pdf(prioralpha,newalpha.at(unit))/pdf(prioralpha,param.alphavector.at(unit)));
            r=samplefromrealuniform(0, 1, &eng);
            if(r<p){param.alphavector.at(unit)=newalpha.at(unit);logl=newlogl;}
            /** End update a **/
            //cout <<"unitb="<<unit<<endl;
            /**update param.b **/
            newbeta.at(unit)=abs(param.betavector.at(unit)+samplefromnormal(0,sigmab,&eng));
            
            for(ward=0;ward<numberofwards;ward++)
            {
              paramnew.alphavector.at(ward)=param.alphavector.at(ward);
              paramnew.betavector.at(ward)=param.betavector.at(ward);
              paramnew.gammavector.at(ward)=param.gammavector.at(ward);
              paramnew.nuvector.at(ward)=param.nuvector.at(ward);
              paramnew.muvector.at(ward)=param.muvector.at(ward);
            }
            paramnew.betavector.at(unit)=newbeta.at(unit);
            paramnew.frequencydependence= param.frequencydependence;
            // ADDED Environmental contamination E
            newlogl=loglikelihood(acquisitionday,acquisition,culturearray,N,I,numberofacquisitions,admissionstate,paramnew,startperiod,endperiod,admissionday,multipleparameters,E,addingenvironment,bedwise,previousCol,present);
            
            p=min(1.,exp(newlogl-logl)*pdf(priorbeta,newb)/pdf(priorbeta,param.b));
            r=samplefromrealuniform(0, 1, &eng);
            if(r<p){param.betavector.at(unit)=newbeta.at(unit);logl=newlogl;}
            /** End update b **/
            
            // ADDED
            /**update param.c **/
            newgamma.at(unit)=abs(param.gammavector.at(unit)+samplefromnormal(0,sigmac,&eng));
            
            for(ward=0;ward<numberofwards;ward++)
            {
              paramnew.alphavector.at(ward)=param.alphavector.at(ward);
              paramnew.betavector.at(ward)=param.betavector.at(ward);
              paramnew.gammavector.at(ward)=param.gammavector.at(ward);
              paramnew.nuvector.at(ward)=param.nuvector.at(ward);
              paramnew.muvector.at(ward)=param.muvector.at(ward);
            }
            paramnew.gammavector.at(unit)=newgamma.at(unit);
            paramnew.frequencydependence= param.frequencydependence;
            // ADDED Environmental contamination E
            newlogl=loglikelihood(acquisitionday,acquisition,culturearray,N,I,numberofacquisitions,admissionstate,paramnew,startperiod,endperiod,admissionday,multipleparameters,E,addingenvironment,bedwise,previousCol,present);
            
            p=min(1.,exp(newlogl-logl)*pdf(priorgamma,newc)/pdf(priorgamma,param.c));
            r=samplefromrealuniform(0, 1, &eng);
            if(r<p){param.gammavector.at(unit)=newgamma.at(unit);logl=newlogl;}
            /** End update c **/
            
            // No update for nu, since param.nu = 1
            
            /**update param.mu **/
            newmuvector.at(unit)=abs(param.muvector.at(unit)+samplefromnormal(0,sigmamu,&eng));
            
            for(ward=0;ward<numberofwards;ward++)
            {
              paramnew.alphavector.at(ward)=param.alphavector.at(ward);
              paramnew.betavector.at(ward)=param.betavector.at(ward);
              paramnew.gammavector.at(ward)=param.gammavector.at(ward);
              paramnew.nuvector.at(ward)=param.nuvector.at(ward);
              paramnew.muvector.at(ward)=param.muvector.at(ward);
            }
            paramnew.muvector.at(unit)=newmuvector.at(unit);
            paramnew.frequencydependence= param.frequencydependence;
            // ADDED Environmental contamination E
            newlogl=loglikelihood(acquisitionday,acquisition,culturearray,N,I,numberofacquisitions,admissionstate,paramnew,startperiod,endperiod,admissionday,multipleparameters,E,addingenvironment,bedwise,previousCol,present);
            
            p=min(1.,exp(newlogl-logl)*pdf(priormu,newmu)/pdf(priormu,param.mu));
            r=samplefromrealuniform(0, 1, &eng);
            if(r<p){param.muvector.at(unit)=newmuvector.at(unit);logl=newlogl;}
            /** End update mu **/
            // END ADDED
            
          }
        }// End update multiple parameters
        
        /**calculate number of false negative tests,
        could be changed in changeloglandI as well, to avoid running over all culures, but gain limited
        **/
        //truepositiveandfalsenegative(&truepositive,&falsenegative,culturearray,admissionday,colstatus,startperiod,endperiod);
        falsenegativefunction(&falsenegative,culturearray,admissionday,colstatus,startperiod,endperiod);
        //cout << "TP="<<truepositive<<" FN="<<falsenegative<<endl;
        
        /********************************************************************************
        Update param.phi
        Sample from beta-distribution, see web appendix AJE Worby
        Suppose X~Gamma(1,α) and Y~Gamma(1,β). Then Z=X/(X+Y) has distribution Beta(α,β)
        *********************************************************************************/
        x=samplefromgamma(aa+truepositive,1.,&eng);
        y=samplefromgamma(bb+falsenegative,1.,&eng);
        logl-=truepositive*log(param.phi)+falsenegative*log(1.-param.phi);
        param.phi=x/(x+y);
        logl+=truepositive*log(param.phi)+falsenegative*log(1.-param.phi);
        //cout << "Update param.phi" << endl;
        //cout << "logl = "<< logl << endl;
        /*End update phi*/
        /* NOTE: Code can be made slightly faster by updating falsenegative, posatadmission and negatadmission in file changeloglandI.cpp */
        //cout <<"unitd="<<unit<<endl;
        /*********************************************************************************
        Update param.f
        Sample from beta-distribution, see web appendix AJE Worby
        Suppose X~Gamma(1,α) and Y~Gamma(1,β). Then Z=X/(X+Y) has distribution Beta(α,β)
        *********************************************************************************/
        numberposatadmisison(&posatadmission,&negatadmission, admissionstate,admissionday,startperiod, endperiod);
        x=samplefromgamma(aa+posatadmission,1.,&eng);
        y=samplefromgamma(bb+negatadmission,1.,&eng);
        logl-=posatadmission*log(param.f)+negatadmission*log(1.-param.f);
        if(importation==1) param.f=x/(x+y);
        logl+=posatadmission*log(param.f)+negatadmission*log(1.-param.f);
        //cout << "Update param.f" << endl;
        //cout << "logl = "<< logl << endl;
        /*End update f*/
      }
    }
    
    /**********************************************************************************
    Update environmental contamination after parameter update
    **********************************************************************************/
    //if(addingenvironment ==1) updateenvironment(param,I,N,E,numberofwards,maxdate);
    if(addingenvironment ==1) updateenvironment2(param,I,N,admissionstate,acquisition,acquisitionday,admissionday,dischargeday,E,E_stay,E_dis,numberofwards,maxdate);
    
    // ADDED Environmental contamination E
    if(onlyaugmented==1) newlogl=loglikelihood(acquisitionday,acquisition,culturearray,N,I,numberofacquisitions,admissionstate,param,startperiod,endperiod,admissionday,multipleparameters,E,addingenvironment,bedwise,previousCol,present);
    
    /**********************************************************************************
    Update augmented data
    **********************************************************************************/
    for(dummy=0;dummy<factoraugmented;dummy++){
      donothing=0;
      event=samplefromintuniform(0,2,&eng);
      //cout <<"event="<<event<<endl;
      switch(event)
      {
      case 0:   /**move a colonization time*/
    thesize=(int)(listdefpos.size()+listtemppos.size());
        if (thesize>0)
        {
          index=samplefromintuniform(0,thesize-1,&eng);
          if(index<(int)listdefpos.size())
          {
            patient=listdefpos.at(index);
            //cout << "Case 0:Listedefinitepositive patient=" << patient << endl;
          }
          else
          {
            patient=listtemppos.at(index-listdefpos.size());
            //cout << "Case 0:Listtemppositive patient=" << patient << endl;
          }
          r=samplefromrealuniform(0, 1, &eng);
          if(r<w)
          {
            newadmissionstate=1;
            
            // ADDED
            newacquisition=0;
            newacquisitionday=0;
            //
          }
          else
          {/* negative at admission*/
    if(possibleacquistiondays.at(patient)<=0){donothing=1;}
    else{
      newadmissionstate=0;
      newacquisition=1;
      indexday=samplefromintuniform(0,possibleacquistiondays.at(patient)-1,&eng);
      newacquisitionday=indexday;
    }
          }
          q=1;
          //if(admissionstate.at(patient)==1&newadmissionstate==1)q=1.;
          if(admissionstate.at(patient)==0&newadmissionstate==1)q=(1.-w)/w/(double)(possibleacquistiondays.at(patient));
          //if(admissionstate.at(patient)==0&newadmissionstate==0)q=1.;
          if(admissionstate.at(patient)==1&newadmissionstate==0)q=w*(double)(possibleacquistiondays.at(patient))/(1-w);
        }
        else{donothing=1;}
        break;
      case 1:  /** add a colonization time */
    if(listneverpos.size()>0)
    {
      index=samplefromintuniform(0,(int)listneverpos.size()-1,&eng);
      patient=listneverpos.at(index);
      //cout << "Case 1: Listneverpositive patient=" << patient << endl;
      r=samplefromrealuniform(0, 1, &eng);
      if(r<w)
      {
        /*patient admitted positive */
        newadmissionstate=1;
        newacquisition=0;
        newacquisitionday=0;
        q=(double)listneverpos.size()/w/(double)(listtemppos.size()+1);
      }
      else
      {
        /* patient uncolonized upon admission */
        //cout <<"patient uncolonized upon admission"<<endl;
        newadmissionstate=0;
        newacquisition=1;
        //cout <<"choices="<<possibleacquistiondays.at(patient)-1;
        indexday=samplefromintuniform(0,possibleacquistiondays.at(patient)-1,&eng);
        //cout<<" indexday="<<indexday<<endl;
        newacquisitionday=indexday;
        q=(double)listneverpos.size()*possibleacquistiondays.at(patient)/(1.-w)/(double)(listtemppos.size()+1);
      } 
    }
    else{donothing=1;}
    //cout<<"end case 1"<<endl;
    break;
      case 2:   /**remove a colonzation time */
        if(listtemppos.size()>0)
        {
          index=samplefromintuniform(0,(int)listtemppos.size()-1,&eng);
          patient=listtemppos.at(index);
          //cout << "Case 2: Listtemppositive patient=" << patient << endl;
          newadmissionstate=0;
          newacquisition=0;
          newacquisitionday=0;
          if(admissionstate.at(patient)==1)
          {
            q=listtemppos.size()*w/(double)(listneverpos.size()+1);
          }
          else
          {
            q=listtemppos.size()*(1-w)/(double)(listneverpos.size()+1)/(possibleacquistiondays.at(patient));
          }
        }
        else{donothing=1;}
        //cout<<"end case 2"<<endl;
        break;
      default:donothing=1;break;
      
      }
      if(donothing==0)//something is changing!
      {
        oldacquisition= acquisition.at(patient);
        oldacquisitionday=acquisitionday.at(patient);
        oldadmissionstate=admissionstate.at(patient);
        oldlogl=logl;
        //cout << "oldlogl = " << oldlogl << endl;
        
        // ADDED
        // Change vectors acquisition to be used in function loglikelihood
        if(addingenvironment == 1){
          acquisition.at(patient)=newacquisition;
          acquisitionday.at(patient)=newacquisitionday;
          admissionstate.at(patient)=newadmissionstate;
          // Change I according to changes in acquisition, acquisitionday, admissionstate etc.
          numberofinfectiouspatients2(I,whereabouts,colstatus,admissionday,acquisitionday,admissionstate,acquisition,numberofacquisitions,maxdate,numberofwards);
          
          //updateenvironment(param,I,N,E,numberofwards,maxdate);
          updateenvironment2(param,I,N,admissionstate,acquisition,acquisitionday,admissionday,dischargeday,E,E_stay,E_dis,numberofwards,maxdate);
          // Include Environmenmtal Contamination E, use "normal" loglikelihood function for added environmental contamination
          logl = loglikelihood(acquisitionday,acquisition,culturearray,N,I,numberofacquisitions,admissionstate,param,startperiod,endperiod,admissionday,multipleparameters,E,addingenvironment,bedwise,previousCol,present);
        }else{
          if(bedwise==1){
            acquisition.at(patient)=newacquisition;
            acquisitionday.at(patient)=newacquisitionday;
            admissionstate.at(patient)=newadmissionstate;
            // Change I according to changes in acquisition, acquisitionday, admissionstate etc.
            numberofinfectiouspatients2(I,whereabouts,colstatus,admissionday,acquisitionday,admissionstate,acquisition,numberofacquisitions,maxdate,numberofwards);
            logl = loglikelihood(acquisitionday,acquisition,culturearray,N,I,numberofacquisitions,admissionstate,param,startperiod,endperiod,admissionday,multipleparameters,E,addingenvironment,bedwise,previousCol,present);
          }else{
            changeloglandI(patient,oldacquisitionday,oldadmissionstate,oldacquisition,newacquisitionday,newadmissionstate,newacquisition,
                           startperiod,endperiod,&logl,param,present,I,N,colstatus,culturearray,cultureperpatient,whereabouts,numberofacquisitions,admissionday,multipleparameters);
          }
        }
        // END ADDED
        
        
        //logp=min(0.,logl-oldlogl+log(q));
        //p=exp(logp);
        p=min(1.,exp(logl-oldlogl)*q);
        r=samplefromrealuniform(0, 1, &eng);
        
        if(r>p)/* no change, recover originalvalues */
        {
          //cout <<"no change, recover originalvalues"<<endl;
          // ADDED
          // Change back (recover original values)
          if(addingenvironment == 1){
            acquisition.at(patient)=oldacquisition;
            acquisitionday.at(patient)=oldacquisitionday;
            admissionstate.at(patient)=oldadmissionstate;
            numberofinfectiouspatients2(I,whereabouts,colstatus,admissionday,acquisitionday,admissionstate,acquisition,numberofacquisitions,maxdate,numberofwards);
            //updateenvironment(param,I,N,E,numberofwards,maxdate);
            updateenvironment2(param,I,N,admissionstate,acquisition,acquisitionday,admissionday,dischargeday,E,E_stay,E_dis,numberofwards,maxdate);
            logl = loglikelihood(acquisitionday,acquisition,culturearray,N,I,numberofacquisitions,admissionstate,param,startperiod,endperiod,admissionday,multipleparameters,E,addingenvironment,bedwise,previousCol,present);
          } else{
            if(bedwise==1){
              acquisition.at(patient)=oldacquisition;
              acquisitionday.at(patient)=oldacquisitionday;
              admissionstate.at(patient)=oldadmissionstate;
              numberofinfectiouspatients2(I,whereabouts,colstatus,admissionday,acquisitionday,admissionstate,acquisition,numberofacquisitions,maxdate,numberofwards);
              logl = loglikelihood(acquisitionday,acquisition,culturearray,N,I,numberofacquisitions,admissionstate,param,startperiod,endperiod,admissionday,multipleparameters,E,addingenvironment,bedwise,previousCol,present);
            }else{
              changeloglandI(patient,oldacquisitionday,oldadmissionstate,oldacquisition,newacquisitionday,newadmissionstate,newacquisition,
                             startperiod,endperiod,&logl,param,present,I,N,colstatus,culturearray,cultureperpatient,whereabouts,numberofacquisitions,admissionday,multipleparameters);
            }
          }
          // END ADDED
          
        }
        else
        {
          if(bedwise==1){
            updatepreviouscol(patient,previouspatarray,previousCol,acquisition,admissionstate);
          }
          // if(addingenvironment == 0){
          //   acquisition.at(patient)=newacquisition;
          //   acquisitionday.at(patient)=newacquisitionday;
          //   admissionstate.at(patient)=newadmissionstate;
          //   numberofinfectiouspatients2(I,whereabouts,colstatus,admissionday,acquisitionday,admissionstate,acquisition,numberofacquisitions,maxdate,numberofwards);
          //   
          //   // Update previousCol array with patient where acquisition or admission might have changed 
          // 
          // }
          
          if(event==1&&listneverpos.size()>0)
          {
            addelementtovector(listtemppos,patient);
            deleteindexelementvector(listneverpos,index);
          }
          if(event==2&&listtemppos.size()>0)
          {
            addelementtovector(listneverpos,patient);
            deleteindexelementvector(listtemppos,index);
          }
        }
      }
      // cout<< "logl in augmented data = "<<logl<<"\n";
      

    }// End of for loop for data augmentation
    
    // Compute mean prevalence 
    int itotal=0, ntotal=0;
    for(unit=0;unit<numberofwards;unit++){
      for(day=0;day<maxdate;day++){
        itotal += I.at(unit).at(day); 
        ntotal += N.at(unit).at(day);
      }
    }
    meanPrev = (double)itotal/(double)ntotal;
  
    
    
    if(i/thinning==(i+thinning-1)/thinning)
    {
      /** Computing contribution of each transmission route *************/
      /** ONLY FOR FREQUENCY DEPENDENCE *********************************/
      numberofinfectiouspatients2(I,whereabouts,colstatus,admissionday,acquisitionday,admissionstate,acquisition,numberofacquisitions,maxdate,numberofwards);
      contralpha=0.0; contrbeta=0.0; contrE=0.0; contrCrossT=0.0; contrEnv=0.0, contrprev=0.0;
      double betaPart=0.0, envPart=0.0, prevPart=0.0;
      numberimportations = 0;
      int previous = 0;
      double delta_t=0;
      if(i==200 || i==100) mypatients << "Force of infection, i=200:\n";
      for(unit=0;unit<numberofwards;unit++){
        for(patient=0; patient<numberofpatients; patient++){
          if(admissionstate.at(patient)==1) numberimportations += 1;
          if(acquisition.at(patient)==1){
            acqday = admissionday.at(patient)+acquisitionday.at(patient);
            delta_t <- acquisitionday.at(patient);
            if(bedwise==1) previous = previousCol.at(patient);
            try{
              foi1(param,I.at(unit).at(acqday),N.at(unit).at(acqday),E.at(unit).at(acqday),addingenvironment,forceofinfection,bedwise,previousCol.at(patient),delta_t);
            } catch(const std::out_of_range& e1) {
              cout <<"Out of range in forceofinfection: "<<"Patient "<<patient<<", admissionday="<< admissionday.at(patient)<<", dischargeday="<<dischargeday.at(patient)<<", acquisitionday="<<acqday<< endl;
            }
            contralpha+=param.a/forceofinfection;
            try{
              betaPart=param.b*((double)I.at(unit).at(acqday)/(double)N.at(unit).at(acqday))/forceofinfection;
              contrbeta+=betaPart;
            } catch(const std::out_of_range& e2) {
              cout<<"Out of range in contrbeta: "<<"Patient "<<patient<<", admissionday="<<admissionday.at(patient)<<", dischargeday="<<dischargeday.at(patient)<<", acquisitionday="<<acqday<<endl;
            }
            prevPart=(double)bedwise*(double)previous*param.p*exp(-delta_t);
            contrprev+=prevPart/forceofinfection;
            try{
              contrE+=(double)addingenvironment*param.c*E.at(unit).at(acqday)/forceofinfection;
            } catch(const std::out_of_range& e3) {
              cout<<"Out of range in contrE: "<<"Patient "<<patient<<", admissionday="<<admissionday.at(patient)<<", dischargeday="<<dischargeday.at(patient)<<", acquisitionday="<<acqday<<endl;
            }
            try{
              contrCrossT+=betaPart+(double)addingenvironment*param.c*(E_stay.at(unit).at(acqday))/forceofinfection;
            } catch(const std::out_of_range& e4) {
              cout<<"Out of range in contrCrossT."<<"Patient "<<patient<<", admissionday="<<admissionday.at(patient)<<", dischargeday="<<dischargeday.at(patient)<<", acquisitionday="<<acqday<<endl;
            }
            try{
              envPart = (double)addingenvironment*param.c*(E.at(unit).at(acqday)-E_stay.at(unit).at(acqday))/forceofinfection;
              contrEnv+= envPart;
            } catch(const std::out_of_range& e5) {
              cout<<"Out of range in contrEnv."<<"Patient "<<patient<<", admissionday="<<admissionday.at(patient)<<", dischargeday="<<dischargeday.at(patient)<<", acquisitionday="<<acqday<<endl;
            }
            //if(envPart<0) cout<<"EnvPart negative. EnvPart="<<envPart<<", Patient "<<patient<<", admissionday="<<admissionday.at(patient)<<", dischargeday="<<dischargeday.at(patient)<<", acquisitionday="<<acqday<<endl;
            if(i==200 || i==100){
              mypatients<<forceofinfection <<",";
            }
          }
        }
      }
      
      totaldaysatrisk(daysatrisk,totalacquisitions,numberofacquisitions,N,I,startperiod,endperiod);
      myfile <<i<<",";
      myfile <<param.f<<",";
      myfile <<param.phi<<",";
      myfile <<logl<<",";
      myfile <<numberimportations<<",";
      myfile <<meanPrev<<",";
      myfile <<param.a<<",";
      myfile <<param.b<<",";
      myfile <<param.p<<",";
      myfile <<param.c<<",";
      myfile <<param.nu<<",";
      myfile <<param.mu<<",";
      myfile <<param.E0<<",";
      
      myfile <<truepositive<<","<<falsenegative<<",";
      myfile <<(double)aaccepted/(double)inew<<","<<(double)baccepted/(double)inew<<",";
      myfile <<(double)paccepted/(double)inew<<",";
      myfile <<(double)caccepted/(double)inew<<","<<(double)muaccepted/(double)inew<<",";
      myfile <<(double)sigmaa<<","<<(double)sigmab<<","<<(double)sigmap<<","<<(double)sigmac<<","<<(double)sigmamu<<",";
      
      for(unit=0;unit<numberofwards;unit++)
      {
        myfile <<"acq["<<unit<<"]="<<totalacquisitions.at(unit);
        myfile <<", daysatr["<<unit<<"]="<<daysatrisk.at(unit)<<",";
      }
      
      myfile <<contralpha<<",";
      myfile <<contrbeta<<",";
      myfile <<contrprev<<",";
      myfile <<contrE<<",";
      myfile <<contrCrossT<<",";
      myfile <<contrEnv<<"\n";
      
      if(i==100 || i==200){
        mypatients << endl;
        mypatients << "i=" << i << endl;
        mypatients << "logl="<<loglikelihood(acquisitionday,acquisition,culturearray,N,I,numberofacquisitions,admissionstate,param,startperiod,endperiod,admissionday,multipleparameters,E,addingenvironment,bedwise,previousCol,present)<<endl;
        truepositiveandfalsenegative(&truepositive,&falsenegative,culturearray,admissionday,colstatus,startperiod,endperiod);
        mypatients << "TP="<< truepositive << "; FN=" << falsenegative << endl;
        mypatients << param.a << ", " << param.b << ", " << param.p << ", " << param.f << ", " << param.phi << endl;
        mypatients << "I per day:" << endl;
        for(int d=startperiod;d<endperiod-1;d++){
          mypatients << I.at(0).at(d)<< ",";
        }
        mypatients << I.at(0).at(endperiod-1)<< "\n";
        mypatients << "N per day:" << endl;
        for(int d=startperiod;d<endperiod;d++){
          mypatients << N.at(0).at(d)<< ",";
        }
        mypatients<<N.at(0).at(endperiod-1)<< "\n";

        mypatients << "Previous colonized bed opccupant?" << endl;
        for(pat=0;pat<numberofpatients-1;pat++){
          mypatients << previousCol.at(pat) << ",";
        }
        mypatients << previousCol.at(numberofpatients-1) << "\n";

        mypatients << "Importation?" << endl;
        for(pat=0;pat<numberofpatients-1;pat++){
          mypatients << admissionstate.at(pat) << ",";
        }
        mypatients << admissionstate.at(numberofpatients-1) << "\n";

        mypatients << "Acquisition time:" << endl;
        for(pat=0;pat<numberofpatients-1;pat++){
          if(acquisition.at(pat)==1) mypatients << admissionday.at(pat) + acquisitionday.at(pat) <<",";
          else mypatients << "-1,";
        }
        if(acquisition.at(numberofpatients-1)==1) mypatients<<admissionday.at(numberofpatients-1)+acquisitionday.at(numberofpatients-1)<<"\n";
        else mypatients << "-1\n";

        mypatients << "Numberofacquisitions per day:" << endl;
        for(int d=startperiod;d<endperiod-1;d++){
          mypatients << numberofacquisitions.at(0).at(d)<< ",";
        }
        mypatients << numberofacquisitions.at(0).at(endperiod-1)<< "\n";
        mypatients << endl;
        
        long double logl1=0;
        double loglacq=0,loglnoacq=0,logladm=0,loglsens=0;
        double foipr=0;
        int patid,sumPat=0,numnoacq=0;
        pat=0;
        unsigned int culture;
        ntotal=0;
        itotal=0;
        if(multipleparameters==0){
          for(day=startperiod;day<endperiod;day++){
            for(unit=0;unit<N.size();unit++){
              itotal+=I.at(unit).at(day);
              ntotal+=N.at(unit).at(day);
              if(N.at(unit).at(day)-I.at(unit).at(day)>0){
                sumPat = present.at(unit).at(day).size();
                for(int j=0;j<sumPat;j++){
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
              }
            }
          }
          logl1 += loglnoacq + loglacq;
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
                logl1+= logpnoacquisitionward(param,I.at(unit).at(day),N.at(unit).at(day),unit,E.at(unit).at(day),addingenvironment)*(N.at(unit).at(day)-I.at(unit).at(day)-numberofacquisitions.at(unit).at(day));
                //cout<<"b"<<endl;
                /** day added */  
                logl1+= logpacquisitionward(param,I.at(unit).at(day),N.at(unit).at(day),unit,E.at(unit).at(day),addingenvironment)*numberofacquisitions.at(unit).at(day);
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
            logl1+=log(param.f);logladm+=log(param.f);
          }
          else
          {
            logl1+=log(1.-param.f);logladm+=log(1.-param.f);
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
                logl1+=log(1.0-param.phi);//false negative culture
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
                logl1+=log(param.phi);//positive culture for colonized patient
                loglsens+=log(param.phi);
                //cout << "True positive loglsens=" << loglsens << endl;
              }
              else
              {
                logl1-=1000000.0; //culture positive and patient uncolonized
                cout << "case logl-=1000000.0: logl=" << logl << endl;
                //loglsens-=1000000.0;
              }
            }
          }
          
        }
        mypatients<<"logl1="<<logl1<<", loglacq="<<loglacq<<", loglnoacq="<<loglnoacq<<"logladm="<<logladm<<", loglsens="<<loglsens<<endl;
      }
      
      /*
      for(unit=0;unit<numberofwards;unit++){
      for(day=0;day<maxdate-1;day++){
      myprevfile <<I.at(unit).at(day)<<",";
      myNfile <<N.at(unit).at(day)<<",";
      myacqfile << numberofacquisitions.at(unit).at(day) << ",";
      }
      myprevfile <<I.at(unit).at(maxdate-1)<<endl;
      myNfile <<N.at(unit).at(maxdate-1)<<endl;
      myacqfile << numberofacquisitions.at(unit).at(maxdate-1) << endl;
      }
      
      for(patient=0; patient<numberofpatients-1; patient++){
      myimpfile << admissionstate.at(patient) << ",";
      colstat << previousCol.at(patient) << ",";
      if(acquisition.at(patient)==1) mypatients << admissionday.at(patient) + acquisitionday.at(patient) <<",";
      else mypatients << "-1,";
      }
      if(acquisition.at(numberofpatients-1)==1) mypatients << admissionday.at(numberofpatients-1) + acquisitionday.at(numberofpatients-1)<<endl;
      else mypatients << "-1"<<endl;
      myimpfile << admissionstate.at(numberofpatients-1) << endl;
      colstat << previousCol.at(numberofpatients-1) << endl;
      */
      
      /*
      if(i==90000){
      for(unit=0;unit<numberofwards;unit++){
      for(day=0;day<maxdate-1;day++){
      myfoifile <<(double)endogenous*param.a+(double)exogenous*param.b*((double)I.at(unit).at(day)/(double)N.at(unit).at(day))+(double)addingenvironment*param.c*E.at(unit).at(day)<<",";
      myprevfile <<I.at(unit).at(day)<<",";
      myNfile <<N.at(unit).at(day)<<",";
      myacqfile << numberofacquisitions.at(unit).at(day) << ",";
      
      myenvfile <<E.at(unit).at(day)<<",";
      myenvstayfile <<E_stay.at(unit).at(day)<<",";
      myenvdisfile <<E_dis.at(unit).at(day)<<",";
      }
      
      myfoifile <<param.a+param.b*((double)I.at(unit).at(maxdate-1)/(double)N.at(unit).at(maxdate-1))+(double)addingenvironment*param.c*E.at(unit).at(maxdate-1)<<endl;
      myprevfile <<I.at(unit).at(maxdate-1)<<endl;
      myNfile <<N.at(unit).at(maxdate-1)<<endl;
      myacqfile << numberofacquisitions.at(unit).at(maxdate-1) << endl;
      
      
      myenvfile <<E.at(unit).at(maxdate-1)<<endl;
      myenvstayfile <<E_stay.at(unit).at(maxdate-1)<<endl;
      myenvdisfile <<E_dis.at(unit).at(maxdate-1)<<endl;
      }
      }
      */
      
      /*
      for(unit=0;unit<numberofwards;unit++){
      for(patient=0; patient<numberofpatients-1; patient++){
      if(acquisition.at(patient)==1) myacqfile << admissionday.at(patient) + acquisitionday.at(patient) <<", ";
      else myacqfile << "-1, ";
      }
      if(acquisition.at(numberofpatients-1)==1) myacqfile << admissionday.at(numberofpatients-1) + acquisitionday.at(numberofpatients-1);
      else myacqfile << "-1, ";
      myacqfile << endl;
      }
      * **/
      
      /**
      if(i==numberofupdates-thinning){
      // Save data set into file: Patient nr, admission state,
      // acquisition day, discharge day
      mypatients << "Patient, " << "admstate, " <<"acqstate, "<< "acqday, " << "admday, " << "disday" << endl;
      for(int p=0;p<dischargeday.size();p++){
      mypatients << p << ", "<<admissionstate.at(p)<<", " <<acquisition.at(p)<<", "<< admissionday.at(p)+acquisitionday.at(p) << ", "<< admissionday.at(p) << ", " << dischargeday.at(p) << endl;
      }
      mypatients<<"Number of acquisitions per day" << endl;
      for(unit=0;unit<numberofwards;unit++){
      for(day=0;day<maxdate-1;day++){
      mypatients<<numberofacquisitions.at(unit).at(day)<<", ";
      }
      mypatients<<numberofacquisitions.at(unit).at(maxdate-1)<<endl;
      }
      }
      * **/
    }
    
    
    /*
    if(i/thinningfoi==(i+thinningfoi-1)/thinningfoi){
    myfoifile << i << ", ";
    for(unit=0;unit<numberofwards;unit++){
    for(day=0;day<maxdate-1;day++){
    myfoifile <<(double)endogenous*param.a+(double)exogenous*param.b*((double)I.at(unit).at(day)/(double)N.at(unit).at(day))+(double)addingenvironment*param.c*E.at(unit).at(day)<<",";
    myprevfile <<I.at(unit).at(day)<<",";
    myNfile <<N.at(unit).at(day)<<",";
    myacqfile << numberofacquisitions.at(unit).at(day) << ",";
    
    myenvfile <<E.at(unit).at(day)<<",";
    myenvstayfile <<E_stay.at(unit).at(day)<<",";
    myenvdisfile <<E_dis.at(unit).at(day)<<",";
    }
    
    myfoifile <<param.a+param.b*((double)I.at(unit).at(maxdate-1)/(double)N.at(unit).at(maxdate-1))+(double)addingenvironment*param.c*E.at(unit).at(maxdate-1)<<endl;
    myprevfile <<I.at(unit).at(maxdate-1)<<endl;
    myNfile <<N.at(unit).at(maxdate-1)<<endl;
    myacqfile << numberofacquisitions.at(unit).at(maxdate-1) << endl;
    
    
    myenvfile <<E.at(unit).at(maxdate-1)<<endl;
    myenvstayfile <<E_stay.at(unit).at(maxdate-1)<<endl;
    myenvdisfile <<E_dis.at(unit).at(maxdate-1)<<endl;
    }
    }
    */
    
    
    if(i/thinningscreen==(i+thinningscreen-1)/thinningscreen)
    {
      cout<<"i="<<i<<", ";
      cout <<"a="<<param.a<<", ";
      cout <<"b="<<param.b<<", ";
      cout <<"f="<<param.f<<", ";
      // ADDED
      cout <<"p="<<param.p<<", ";
      cout <<"c="<<param.c<<", ";
      cout <<"nu="<<param.nu<<", ";
      cout <<"mu="<<param.mu<<", ";
      cout <<"E0="<<param.E0<<", ";
      // END ADDED
      cout <<"phi="<<param.phi<<", ";
      cout <<"TP="<<truepositive<<", FN="<<falsenegative<<", ";
      cout <<"aa="<<(double)aaccepted/(double)inew<<", ba="<<(double)baccepted/(double)inew<<", ";
      cout <<"ca="<<(double)caccepted/(double)inew<<", nua="<<(double)nuaccepted/(double)inew<<", mua="<<(double)muaccepted/(double)inew<<", ";
      // ADDED sigmac, sigmanu, sigmamu
      cout <<"sa="<<(double)sigmaa<<", sb="<<(double)sigmab<<", sc="<<(double)sigmac<<",snu="<<(double)sigmanu<<", smu="<<(double)sigmamu<<", ";
      
      totaldaysatrisk(daysatrisk,totalacquisitions,numberofacquisitions,N,I,startperiod,endperiod);
      for(unit=0;unit<numberofwards;unit++)
      {
        cout <<"acq["<<unit<<"]="<<totalacquisitions.at(unit)<<", ";
        cout <<"daysatr["<<unit<<"]="<<daysatrisk.at(unit)<<", ";
      }
      cout <<"numImp="<<numberimportations<<", ";
      cout <<"meanPrev="<<meanPrev<<", ";
      cout <<"logl="<<logl<<"\n ";
    }
    i+=1;
    }
  
  cout<<"Number of changes (sigma) = "<<nchange<<endl;
  mystartfile<<nchange*(double)updatestochangesigma/10+1<<endl;
  cout<<"Line = "<< nchange*(double)updatestochangesigma/10 + 1<<endl;
  /**********************************************************************
  * Compute running time
  * *******************************************************************/
  time_t timer2 = time(NULL);
  double seconds = difftime(timer2, timer);
  cout<< "Time to compute = "<<seconds<<"s\n";
  mystartfile<<seconds<<endl;
  /*********************************************************************/
  
  /// CLOSE STREAMS
  mystartfile.close();
  myfile.close();
  
  myfoifile.close();
  myacqfile.close();
  myimpfile.close();
  myprevfile.close();
  myNfile.close();
  
  myenvfile.close();
  myenvstayfile.close();
  myenvdisfile.close();
  
  mypatients.close();
  colstat.close();
  
  
  return 0;
  }