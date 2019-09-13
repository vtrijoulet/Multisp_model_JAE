/////////// Functions needed in the cpp code MS_SM.cpp ///////////////////

// Function that does the inverse logit transformation
template <class Type> 
Type inverse_logit(Type logit_x, Type lower, Type upper){
  Type x = upper-lower/(1+exp(-logit_x));
  return x;
}



// Function to force an argument to be positive
template<class Type>
Type posfun(Type x, Type eps, Type &pen){
  pen += CppAD::CondExpLt(x, eps, Type(0.01) * pow(x-eps,2), Type(0));
  return CppAD::CondExpGe(x, eps, x, eps/(Type(2)-x/eps));
}

// Function that square x
template <class Type> 
Type square(Type x){return x*x;}


// Function with options for age composition data
template<class Type>
Type get_acomp_ll(int year, int n_ages, Type Neff, int age_comp_model, vector<Type> paa_obs, vector<Type> paa_pred, vector<Type> age_comp_pars, int aref)
{
  Type zero = Type(0);
  Type one = Type(1);
  vector<Type> temp_n(n_ages);
  Type temp_Neff = zero, ll = zero;
  if(age_comp_model == 1) //multinomial
  {
    temp_Neff = Neff * exp(age_comp_pars(0));
    temp_n = temp_Neff * paa_obs;
    ll = lgamma(temp_Neff + one);
    for(int a = 0; a < n_ages; a++) ll += -lgamma(temp_n(a) + one) + temp_n(a) * log(paa_pred(a) + Type(1e-15));
  }
  return ll;
}


// Function to simulate age composition data
template<class Type>
vector<Type> sim_acomp(int year, int n_ages, Type Neff, int age_comp_model, vector<Type> paa_obs, vector<Type> paa_pred, vector<Type> age_comp_pars, int aref)
{
  Type one = Type(1);
  vector<Type> obs(n_ages);
  obs.setZero();
  if(age_comp_model == 1) //multinomial generated from condition binomials
  { 
    Type N = Neff * exp(age_comp_pars(0));
    obs(0) = rbinom(N, paa_pred(0));
    for(int a = 1; a < n_ages-1; a++) 
    {
      Type denom = one - paa_pred.head(a).sum();
      Type cond_N = N-obs.head(a).sum();
      if(denom > Type(1.0e-15)) if(cond_N > Type(1.0e-15))
      {
        Type cond_p = paa_pred(a)/denom; //.head first a components
        if(one - cond_p > Type(1.0e-15)) obs(a) = rbinom(cond_N,cond_p);
        else obs(a) = cond_N; //p pretty much 1
      }
    }
    obs(n_ages-1) = N - obs.sum();
    obs = obs/obs.sum();// proportions
  }
  return obs;
}



template <class Type>
Type ddirichlet(vector<Type> obs, vector<Type>p, Type phi, int do_log) 
{
  int n = obs.size();
  Type ll = lgamma(phi);
  vector<Type> alphas = (p)*phi;
  for(int i = 0; i < n; i++) ll +=  -lgamma(alphas(i)) + (alphas(i) - Type(1.0)) * log(obs(i));
  if(do_log == 1) return(ll);
  else return(exp(ll));
}



template<class Type>
vector<Type> rdirichlet(vector<Type> p, Type phi)
{
  int n = p.size();
  vector<Type> alpha = p * phi;
  vector<Type> obs(n);
  for(int i = 0; i < n; i++) obs(i) = rgamma(alpha(i),Type(1.0));
  if (obs.sum()!=0) obs = obs/obs.sum();
  return (obs);
}
