// Risk Evaluation

#include <TMB.hpp>
#include <iostream>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // DATA //
  DATA_VECTOR(node);
  DATA_VECTOR(wt);
  DATA_SCALAR(Time);
  DATA_SCALAR(D);
  DATA_INTEGER(func_type);
  DATA_VECTOR(half_life);
  DATA_INTEGER(one_component);
  
  // PARAMETER //
  PARAMETER_VECTOR(aa);
  PARAMETER_VECTOR(bb);
  PARAMETER_VECTOR(cc);
  PARAMETER_VECTOR(dd);
  PARAMETER_VECTOR(log_sigma);
  PARAMETER(log_sigma_RE);
  PARAMETER(z);
  
  // PARAMETER TRANSFORMATION
  vector<Type> Q = exp(aa);
  vector<Type> P = 1/(1+exp(-bb));  
  vector<Type> K1 = exp(cc);
  vector<Type> K2 = exp(dd);
  vector<Type> HL = log(0.5)/(half_life*365);
  vector<Type> sigma = exp(log_sigma);
  Type sigma_RE = exp(log_sigma_RE);
  
  int N = node.size();
  
  Type log_Mean_134=0.0;
  Type log_Mean_137=0.0;
  
  Type nll=0.0;
  if (one_component==1){
    log_Mean_134 += z+aa(0)+HL(0)*Time-K1(0)*Time;
    log_Mean_137 += z+aa(1)+HL(1)*Time-K1(1)*Time;
  } else{
    log_Mean_134 += z+aa(0)+HL(0)*Time+log(P(0)*exp(-K1(0)*Time)+(1-P(0))*exp(-K2(0)*Time));
    log_Mean_137 += z+aa(1)+HL(1)*Time+log(P(1)*exp(-K1(1)*Time)+(1-P(1))*exp(-K2(1)*Time));
  }
  
  for (int i=0;i<N;i++){
    if (func_type==0 || func_type==1) nll += wt(i)*exp(dnorm(log(D)-node(i),log_Mean_134,sigma(0),true)+log(pnorm(log(D)+log(1-exp(-node(i))),log_Mean_137,sigma(1)))+dnorm(z,Type(0.0),sigma_RE,true));
  }
  
  if (func_type==0) nll = -log(nll);
  
  return nll;
}
