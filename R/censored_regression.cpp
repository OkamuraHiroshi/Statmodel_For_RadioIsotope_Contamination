// Censored Regression

#include <TMB.hpp>
#include <iostream>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // DATA //
  DATA_VECTOR(Cs_134);
  DATA_VECTOR(Cs_137);
  DATA_IVECTOR(Censor_134);
  DATA_IVECTOR(Censor_137);
  DATA_VECTOR(Time);
  DATA_INTEGER(mis);
  DATA_VECTOR(half_life);
  DATA_INTEGER(one_component);

  // PARAMETER //
  PARAMETER_VECTOR(aa);
  PARAMETER_VECTOR(bb);
  PARAMETER_VECTOR(cc);
  PARAMETER_VECTOR(dd);
  PARAMETER_VECTOR(log_sigma);
  PARAMETER(log_sigma_RE);
  PARAMETER_VECTOR(z);
  
  // PARAMETER TRANSFORMATION
  vector<Type> Q = exp(aa);
  vector<Type> P = 1/(1+exp(-bb));  
  vector<Type> K1 = exp(cc);
  vector<Type> K2 = exp(dd);
  vector<Type> HL = log(0.5)/(half_life*365);
  vector<Type> sigma = exp(log_sigma);
  Type sigma_RE = exp(log_sigma_RE);
  
  int N = Cs_134.size();
    
  vector<Type> log_Mean_134(N);
  vector<Type> log_Mean_137(N);
    
  Type nll=0.0;
  
  for (int i=0;i<N;i++){
    nll += -dnorm(z(i),Type(0.0),sigma_RE,true);
  }

  for (int i=0;i<N;i++){
    if (one_component==1){
      log_Mean_134(i) = z(i)+aa(0)+HL(0)*Time(i)-K1(0)*Time(i);
      log_Mean_137(i) = z(i)+aa(1)+HL(1)*Time(i)-K1(1)*Time(i);  
    } else{
      log_Mean_134(i) = z(i)+aa(0)+HL(0)*Time(i)+log(P(0)*exp(-K1(0)*Time(i))+(1-P(0))*exp(-K2(0)*Time(i)));
      log_Mean_137(i) = z(i)+aa(1)+HL(1)*Time(i)+log(P(1)*exp(-K1(1)*Time(i))+(1-P(1))*exp(-K2(1)*Time(i)));
    }
    
    if (Censor_134(i) == 0) nll += -dnorm(log(Cs_134(i)),log_Mean_134(i),sigma(0),true);
    if (Censor_134(i) > 0 && mis==0) nll += -log(pnorm(log(Cs_134(i)),log_Mean_134(i),sigma(0)));
    if (Censor_137(i) == 0) nll += -dnorm(log(Cs_137(i)),log_Mean_137(i),sigma(1),true);
    if (Censor_137(i) > 0 && mis==0) nll += -log(pnorm(log(Cs_137(i)),log_Mean_137(i),sigma(1)));
  }
  
  
  
  // REPORT
  
  REPORT(HL);
  REPORT(aa);
  REPORT(bb);
  REPORT(cc);
  REPORT(dd);
  REPORT(sigma);
  REPORT(sigma_RE);
  REPORT(Q);
  REPORT(P);
  REPORT(K1);
  REPORT(K2);
  ADREPORT(Q);
  ADREPORT(P);
  ADREPORT(K1);
  ADREPORT(K2);
  ADREPORT(sigma);
  ADREPORT(sigma_RE);
          
  return nll;
}
