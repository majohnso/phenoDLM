# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

//////////////////////////////////////////////////////////////////////
// Kalman Filter/FFBS Component functions
//////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
arma::mat at(arma::mat m, arma::mat G) {
  return((G*m.t()).t());
}

// [[Rcpp::export]]
arma::mat Rt(arma::mat G, arma::mat C, arma::mat W) {
  return(G*C*G.t() + W);
}

// [[Rcpp::export]]
arma::mat ft(arma::mat F, arma::mat a) {
  return((F*a.t()).t());
}

// [[Rcpp::export]]
arma::mat Qt(arma::mat F, arma::mat R, arma::mat V) {
  return(F*R*F.t() + V);
}

// [[Rcpp::export]]
arma::mat mt(double Y, arma::mat F, arma::mat a, arma::mat R, arma::mat f, arma::mat Q) {
  return((a.t() + R*F.t()*inv(Q)*(Y - f)).t());
}

// [[Rcpp::export]]
arma::mat Ct(arma::mat F, arma::mat R, arma::mat Q) {
  return(R - R*F.t()*inv(Q)*F*R);
}

// [[Rcpp::export]]
arma::mat ht(arma::mat G, arma::mat C, arma::mat R, arma::mat theta, arma::mat a, arma::mat m) {
  return((m.t() + C*G.t()*inv(R)*(theta - a).t()).t());
}

// [[Rcpp::export]]
arma::mat Ht(arma::mat C, arma::mat G, arma::mat R) {
  return(C - C*G.t()*inv(R)*G*C);
}

// [[Rcpp::export]]
arma::mat ht2(arma::mat G, arma::mat C, arma::mat W, arma::mat theta, arma::mat H, arma::mat m) {
  return((H*(inv(C)*m.t() + G.t()*inv(W)*theta.t())).t());
}

// [[Rcpp::export]]
arma::mat Ht2(arma::mat C, arma::mat G, arma::mat W) {
  return(inv(inv(C) + G.t()*inv(W)*G));
}

//////////////////////////////////////////////////////////////////////
// Drawing the Multivariate Normal
//////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
arma::mat mvrnormChol(int n, arma::mat mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = randn(n, ncols);
  return repmat(mu, 1, n) + Y * chol(sigma);
}

// [[Rcpp::export]]
arma::mat mvrnormSVD(int n, arma::mat mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = randn(n, ncols);
  arma::mat U;
  arma::vec d;
  arma::mat V;
  svd(U, d, V, sigma);
  // return repmat(mu, 1, n) + Y * U*diagmat(sqrt(d))*V;
  return (repmat(mu.t(), 1, n) + U*diagmat(sqrt(d))*Y.t()).t();
}

// [[Rcpp::export]]
double normden(double x, double mu, double sd) {
  return (1/(pow(2*M_PI,0.5)*sd))*exp(-0.5*pow((x-mu)/sd, 2));
}

//////////////////////////////////////////////////////////////////////
// Kalman Filter
//////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
List KalmanFilter(arma::vec Y, arma::mat V, arma::mat W, List model){
  // Set initial values
  arma::mat F = as<arma::mat>(model["F"]);
  arma::mat G = as<arma::mat>(model["G"]);
  arma::vec m0 = as<arma::vec>(model["m0"]);
  arma::mat C0 = as<arma::mat>(model["C0"]);
  
  int T = Y.n_rows;
  int p = G.n_rows;
  
  // storing output
  arma::mat mt_keep(T+1,p);
  arma::mat at_keep(T,p);
  arma::mat ft_keep(T,1);
  List Rt_keep(T);
  List Qt_keep(T);
  List Ct_keep(T+1);
  
  // reused during kalman filter
  arma::mat R;
  arma::mat C;
  arma::mat Q;
  
  mt_keep.row(0) = m0.t();
  Ct_keep[0] = C0;
  C = C0;
  
  for (int i=0; i<T; i++){
    at_keep.row(i) = at(mt_keep.row(i), G);
    R = Rt(G, C, W);
    
    ft_keep.row(i) = ft(F, at_keep.row(i));
    Q = Qt(F, R, V);
    
    if (R_IsNA(Y[i])) {
      mt_keep.row(i+1) = at_keep.row(i);
      C = R;
    } else {
      mt_keep.row(i+1) = mt(Y[i], F, at_keep.row(i), R, ft_keep.row(i), Q);
      C = Ct(F, R, Q);
    }
    Rt_keep[i] = R;
    Qt_keep[i] = Q;
    Ct_keep[i+1] = C;
  }
  return List::create(
    Named("at") = at_keep,
    Named("ft") = ft_keep,
    Named("mt") = mt_keep,
    Named("Rt") = Rt_keep,
    Named("Qt") = Qt_keep,
    Named("Ct") = Ct_keep);
}

//////////////////////////////////////////////////////////////////////
// Fast Forward Backward Sampling (FFBS)
//////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
List FFBS(arma::mat G, arma::mat W, List Kalman, bool SVD = false){
  arma::mat at = as<arma::mat>(Kalman["at"]);
  arma::mat mt = as<arma::mat>(Kalman["mt"]);
  List Rt = as<List>(Kalman["Rt"]);
  List Ct = as<List>(Kalman["Ct"]);
  
  int nrow = at.n_rows;
  int T = nrow - 1;
  int p = G.n_rows;
  
  // storing output
  arma::mat ht_keep(nrow,p);
  arma::mat theta_keep(nrow+1,p);
  List Ht_keep(nrow);
  
  // Draw theta_T from filtering distribution
  if(SVD) {
    theta_keep.row(T+1) = mvrnormSVD(1, mt.row(T+1), Ct[T+1]);
  } else {
    theta_keep.row(T+1) = mvrnormChol(1, mt.row(T+1), Ct[T+1]);
  }
  
  for (int i=1; i<=nrow; i++){
    Ht_keep[nrow-i] = Ht2(Ct[nrow-i], G, W);
    ht_keep.row(nrow-i) = ht2(G, Ct[nrow-i], W, theta_keep.row(nrow-i+1), Ht_keep[nrow-i], mt.row(nrow-i));
    
    // Draw theta from posterior
    if (SVD) {
      theta_keep.row(nrow-i) = mvrnormSVD(1, ht_keep.row(nrow-i), Ht_keep[nrow-i]);
    } else {
      theta_keep.row(nrow-i) = mvrnormChol(1, ht_keep.row(nrow-i), Ht_keep[nrow-i]);
    }
  }
  return List::create(
    Named("theta") = theta_keep,
    Named("ht") = ht_keep,
    Named("Ht") = Ht_keep);
}

//////////////////////////////////////////////////////////////////////
// Calculating Sums of Squares
//////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
arma::mat calc_SSE_error(arma::mat y, arma::mat theta, arma::mat F) {
  int T = y.n_rows;
  arma::mat SSE(1,1);
  SSE.zeros();
  for (int i=0; i<T; i++) {
    if (R_IsNA(y[i])) {
      SSE += 0.0;
    } else {
      SSE += (y[i]-F*theta.row(i+1).t())*(y[i]-F*theta.row(i+1).t());
    }
  }
  return SSE;
}

// [[Rcpp::export]]
arma::mat calc_SSE_theta(arma::mat theta, arma::mat G) {
  int T = theta.n_rows;
  arma::mat SSE(1,1);
  SSE.zeros();
  for (int i=1; i<T; i++) 
    SSE += ((theta.row(i).t()-G*theta.row(i-1).t()).t()*(theta.row(i).t()-G*theta.row(i-1).t()));
  return SSE;
}

//////////////////////////////////////////////////////////////////////
// Sampling Functions
//////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
double cauchy_MH(double par_current, double shape, double rate, double c) {
  double par_proposed = 1/R::rgamma(shape, 1/rate);
  double log_rho = log(1+par_current/(c*c)) - log(1+par_proposed/(c*c));
  return ifelse(log(runif(1)) < log_rho, par_proposed, par_current)[0];
}

//////////////////////////////////////////////////////////////////////
// MCMC
//////////////////////////////////////////////////////////////////////

// Cauchy priors, MH - IG proposal
// [[Rcpp::export]]
List mcmc_seas_IG_NA(
    int n_reps, 
    arma::mat dat, 
    List initial_values,
    NumericVector c,
    List model,
    double T_star,
    bool svd = false,
    bool save_theta = true){
  
  // Set initial values
  double sig2_e = as<double>(initial_values["sig2_e"]);
  double sig2_w = as<double>(initial_values["sig2_w"]);
  
  // Model things
  arma::mat F = as<arma::mat>(model["F"]);
  arma::mat G = as<arma::mat>(model["G"]);
  arma::mat W = as<arma::mat>(model["W"]);
  arma::mat V = as<arma::mat>(model["V"]);
  

  List keep_theta(n_reps);
  NumericVector keep_sigma_e(n_reps);
  NumericVector keep_sigma_w(n_reps);
  
  int T = dat.n_rows;
  int p = G.n_cols;
  
  //reused during sampling
  List kalman;
  arma::mat rate1, rate2;
  double sh1 = (T_star-1.0)/2.0, sh2 = (T*p-1.0)/2.0;
  
  for (int i=0; i<n_reps; i++) {
    // draw theta_0:T
    kalman = KalmanFilter(dat, V, W, model); // Kalman Filter
    arma::mat theta = FFBS(G, W, kalman, svd)[0]; 
    
    // draw sigma_e^2
    rate1 = calc_SSE_error(dat, theta, F)/2.0;
    double rate1b = rate1[0]; 
    sig2_e = cauchy_MH(sig2_e, sh1, rate1b, c[0]); 
    
    // draw sigma_w^2
    rate2 = calc_SSE_theta(theta, G)/2.0;
    double rate2b = rate2[0]; 
    sig2_w = cauchy_MH(sig2_w, sh2, rate2b, c[1]); 
    
    // update V    
    V = sig2_e;
    
    // update W
    for (int j=0; j < p; j++){
      W(j,j) = sig2_w;
    }
    
    // Update storage
    if(save_theta){
      keep_theta[i] = theta;
    }
    keep_sigma_e[i] = sig2_e;
    keep_sigma_w[i] = sig2_w;  
  }
  if(save_theta){
    return List::create(
      Named("sigma2_e") = keep_sigma_e,
      Named("sigma2_w") = keep_sigma_w,
      Named("theta") = keep_theta);
  }else{
    return List::create(
      Named("sigma2_e") = keep_sigma_e,
      Named("sigma2_w") = keep_sigma_w);
  }
}

// Conjugate IG priors
// [[Rcpp::export]]
List mcmc_local_seas_NA(
    int n_reps, 
    arma::mat dat, 
    List initial_values,
    List prior,
    List model,
    double T_star,
    bool svd = false){
  
  // Set initial values
  double psi1 = as<double>(initial_values["psi1"]);
  double psi2 = as<double>(initial_values["psi2"]);

  // Set prior values
  double a1 = as<double>(prior["a1"]); 
  double a2 = as<double>(prior["a2"]);
  double b1 = as<double>(prior["b1"]);
  double b2 = as<double>(prior["b2"]); 

  // Model things
  arma::mat F = as<arma::mat>(model["F"]);
  arma::mat G = as<arma::mat>(model["G"]);
  arma::mat W = as<arma::mat>(model["W"]);
  arma::mat V = as<arma::mat>(model["V"]);
  
  List keep_theta(n_reps);
  NumericVector keep_sigma_e(n_reps);
  NumericVector keep_sigma_w(n_reps);

  int T = dat.n_rows;
  int p = G.n_cols;
  
  //reused during sampling
  List kalman;
  arma::mat rate1, rate2;
  
  double sh1 = a1 + T_star/2.0, sh2 = a2 + T*p/2.0;
  
  for (int i=0; i<n_reps; i++) {
    // draw theta_0:T
    kalman = KalmanFilter(dat, V, W, model); // Kalman Filter
    arma::mat theta = FFBS(G, W, kalman, svd)[0]; 
    
    // draw 1/sigma_e^2
    rate1 = b1 + calc_SSE_error(dat, theta, F)/2.0;
    double rate1b = rate1[0]; 
    psi1 = rgamma(1, sh1, 1/rate1b)[0];
    
    // draw 1/sigma_w^2
    rate2 = b2 + calc_SSE_theta(theta, G)/2.0;
    double rate2b = rate2[0]; 
    psi2 = rgamma(1, sh2, 1/rate2b)[0];
    
    // update V    
    V = 1/psi1;
    
    // update W
    for (int j=0; j < p; j++){
      W(j,j) = 1/psi2;
    }
    
    // Update storage
    keep_theta[i] = theta;
    keep_sigma_e[i] = 1/psi1;
    keep_sigma_w[i] = 1/psi2;  
  }
  
  return List::create(
    Named("theta") = keep_theta,
    Named("sigma2_e") = keep_sigma_e,
    Named("sigma2_w") = keep_sigma_w);
}


