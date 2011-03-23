
// includes from the plugin
#include <RcppArmadillo.h>
#include <Rcpp.h>


#ifndef BEGIN_RCPP
#define BEGIN_RCPP
#endif

#ifndef END_RCPP
#define END_RCPP
#endif

using namespace Rcpp;


// user includes

typedef Rcpp::NumericMatrix::iterator mat_iterator;
typedef Rcpp::NumericVector::iterator vec_iterator;

// declarations
extern "C" {
SEXP updateScoresC( SEXP Z_, SEXP A_, SEXP F_) ;
SEXP updateLoadingsVarC( SEXP tauinv_, SEXP A_, SEXP taua_, SEXP taub_) ;
SEXP updateLoadingsVarJC( SEXP tauinv_, SEXP A_, SEXP taua_, SEXP taub_) ;
SEXP updateRho( SEXP rho_, SEXP A_, SEXP rhoa_, SEXP rhob_) ;
SEXP updateSparseLoadings( SEXP Z_, SEXP A_, SEXP F_, SEXP tauinv_, SEXP rho_) ;
SEXP updateSparseLoadingsJ( SEXP Z_, SEXP A_, SEXP F_, SEXP tauinv_, SEXP rho_) ;
//SEXP MCMCstep( SEXP Z_, SEXP A_, SEXP F_, SEXP tauinv_, SEXP rho_ );
SEXP MCMCstep( SEXP Z_, SEXP A_, SEXP F_, SEXP tauinv_, SEXP rho_ , SEXP Ra_, SEXP maxes_, SEXP argsorts_, SEXP priors_, SEXP nsim_, SEXP nburn_, SEXP thin_, SEXP printstatus_, SEXP keepscores_, SEXP keeploadings_);
}

// definition

void sampleScores(arma::mat& Z, arma::mat& A, arma::mat& F, int n, int k) {
  arma::mat U = arma::randn(k, n);
  arma::mat fcov = arma::inv(arma::eye(k,k) + arma::trans(A)*A);
  arma::mat R = arma::chol(fcov);
  F = R*U + fcov*arma::trans(A)*Z;
}

SEXP updateScoresC( SEXP Z_, SEXP A_, SEXP F_ ){
BEGIN_RCPP

Rcpp::NumericMatrix Fr(F_);			// creates Rcpp matrix from SEXP
Rcpp::NumericMatrix Ar(A_);			// creates Rcpp matrix from SEXP
Rcpp::NumericMatrix Zr(Z_);			// creates Rcpp matrix from SEXP

int n = Zr.ncol(); 
int k = Ar.ncol(); 
int p = Ar.nrow();

arma::mat Z(Zr.begin(), p, n, false);   	// reuses memory and avoids extra copy
arma::mat A(Ar.begin(), p, k, false);   	// reuses memory and avoids extra copy
arma::mat F(Fr.begin(), k, n, false);   	// reuses memory and avoids extra copy

//arma::mat U = arma::randn(k, n);

//arma::mat fcov = arma::inv(arma::eye(k,k) + arma::trans(A)*A);
//arma::mat R = arma::chol(fcov);
//F = R*U + fcov*arma::trans(A)*Z;

sampleScores(Z, A, F, n, k);

END_RCPP
}

//todo: sampling fcn for below

SEXP updateLoadingsVarC( SEXP tauinv_, SEXP A_, SEXP taua_, SEXP taub_ ){
BEGIN_RCPP

Rcpp::NumericMatrix A(A_);	
Rcpp::NumericVector tauinv(tauinv_);	
double taua = as<double>(taua_);
double taub = as<double>(taub_);

int k = A.ncol(); 
int p = A.nrow();

for (int i=0; i<k; i++) {
	NumericMatrix::Column col = A.column(i) ;
	double ssA = std::inner_product( col.begin(), col.end(), col.begin(), 0.0 );
	double tauafc = taua + 0.5*ssA;
	double taubfc = taub + 0.5*(p);
	tauinv(i) = Rf_rgamma(tauafc, 1.0)/taubfc;
}


END_RCPP
}

void sampleLoadingsVarJ(Rcpp::NumericVector& tauinv, Rcpp::NumericMatrix& A, double& taua, double& taub) {
	int k = A.ncol(); 
	int p = A.nrow();
	for (int i=0; i<p; i++) {
		NumericMatrix::Row row = A.row(i) ;
		double ssA = std::inner_product( row.begin(), row.end(), row.begin(), 0.0 );
		double tauafc = taua + 0.5*ssA;
		double taubfc = taub + 0.5*(k);
		tauinv(i) = Rf_rgamma(tauafc, 1.0)/taubfc;
	}
}

SEXP updateLoadingsVarJC( SEXP tauinv_, SEXP A_, SEXP taua_, SEXP taub_ ){
BEGIN_RCPP

Rcpp::NumericMatrix A(A_);  
Rcpp::NumericVector tauinv(tauinv_);	
double taua = as<double>(taua_);
double taub = as<double>(taub_);

//int k = A.ncol(); 
//int p = A.nrow();

//for (int i=0; i<p; i++) {
//	NumericMatrix::Row row = A.row(i) ;
//	double ssA = std::inner_product( row.begin(), row.end(), row.begin(), 0.0 );
//	double tauafc = taua + 0.5*ssA;
//	double taubfc = taub + 0.5*(k);
//	tauinv(i) = Rf_rgamma(tauafc, 1.0)/taubfc;
//}
sampleLoadingsVarJ(tauinv, A, taua, taub);

END_RCPP
}

void sampleRho(Rcpp::NumericVector& rho, Rcpp::NumericMatrix& A, double rhoa, double rhob){
	int k = A.ncol(); 
	int p = A.nrow();
	for (int i=0; i<k; i++) {
	  NumericMatrix::Column col = A.column(i) ;
	  NumericVector As = ifelse(abs(col)<1e-10, 0.0, 1.0);
	  double nnz = std::accumulate(As.begin(), As.end(), 0.0);
	  double maxnnz = p;
	  rho(i) = Rf_rbeta(rhoa + nnz, rhob + maxnnz-nnz);
	}
}

SEXP updateRho( SEXP rho_, SEXP A_, SEXP rhoa_, SEXP rhob_ ){
BEGIN_RCPP

Rcpp::NumericVector rho(rho_);  
Rcpp::NumericMatrix A(A_);	
double rhoa = as<double>(rhoa_);
double rhob = as<double>(rhob_);

//int k = A.ncol(); 
//int p = A.nrow();

//for (int i=0; i<k; i++) {
//  NumericMatrix::Column col = A.column(i) ;
//  NumericVector As = ifelse(abs(col)<1e-10, 0.0, 1.0);
//  double nnz = std::accumulate(As.begin(), As.end(), 0.0);
//  double maxnnz = p;
//  rho(i) = Rf_rbeta(rhoa + nnz, rhob + maxnnz-nnz);
//}

sampleRho(rho, A, rhoa, rhob);

END_RCPP
}

//todo: sampling fcn for below

SEXP updateSparseLoadings( SEXP Z_, SEXP A_, SEXP F_, SEXP tauinv_, SEXP rho_ ){
BEGIN_RCPP

Rcpp::NumericMatrix Fr(F_);  		// creates Rcpp matrix from SEXP
Rcpp::NumericMatrix Ar(A_);			// creates Rcpp matrix from SEXP
Rcpp::NumericMatrix Zr(Z_);			// creates Rcpp matrix from SEXP
Rcpp::NumericVector tauinv(tauinv_);			// creates Rcpp matrix from SEXP
Rcpp::NumericVector rho(rho_);

GetRNGstate();

int n = Zr.ncol(); 
int k = Ar.ncol(); 
int p = Ar.nrow();

double out = 0.0;

arma::mat Z(Zr.begin(), p, n, false);   	// reuses memory and avoids extra copy
arma::mat A(Ar.begin(), p, k, false);   	// reuses memory and avoids extra copy
arma::mat F(Fr.begin(), k, n, false);   	// reuses memory and avoids extra copy

arma::mat t;
arma::mat uvector;
arma::vec sumf2 = arma::sum(arma::pow(F, 2), 1);

double u=0.0, uu=0.0, v=0.0, loglog=0.0, logc=0.0;
for (int i=0; i<p; i++) {
	for (int j=0; j<k; j++){
		t = Z.row(i) - A.row(i)*F + A(i,j)*F.row(j);
		u = arma::accu(F.row(j)%t);
		v = sumf2(j) + tauinv[j];
		
		loglog = log(u*u) - log(2.0*v);
        
        logc = log(rho[j]) - log(1-rho[j]) - 0.5*log(v) + 0.5*log(tauinv[j]) + exp(loglog);
        
        uu = Rf_runif(0.0, 1.0);

        log(1/uu-1) > logc ? A(i,j) = 0.0 : A(i,j) = Rf_rnorm(u/v, 1/sqrt(v));
		}
	}

PutRNGstate();


END_RCPP
}

void sampleSparseLoadingsJ(arma::mat& Z, arma::mat& A, arma::mat& F, Rcpp::NumericVector& tauinv, Rcpp::NumericVector& rho, int n, int p, int k ){

arma::mat t;
arma::mat uvector;
arma::vec sumf2 = arma::sum(arma::pow(F, 2), 1);

double u=0.0, uu=0.0, v=0.0, loglog=0.0, logc=0.0;
for (int i=0; i<p; i++) {
	for (int j=0; j<k; j++){
		t = Z.row(i) - A.row(i)*F + A(i,j)*F.row(j);
		u = arma::accu(F.row(j)%t);
		v = sumf2(j) + tauinv[i];
		
		loglog = log(u*u) - log(2.0*v);
        
        logc = log(rho[j]) - log(1-rho[j]) - 0.5*log(v) + 0.5*log(tauinv[i]) + exp(loglog);
        
        uu = Rf_runif(0.0, 1.0);

        log(1/uu-1) > logc ? A(i,j) = 0.0 : A(i,j) = Rf_rnorm(u/v, 1/sqrt(v));
		}
	}
}

SEXP updateSparseLoadingsJ( SEXP Z_, SEXP A_, SEXP F_, SEXP tauinv_, SEXP rho_ ){
BEGIN_RCPP

Rcpp::NumericMatrix Fr(F_);  		// creates Rcpp matrix from SEXP
Rcpp::NumericMatrix Ar(A_);			// creates Rcpp matrix from SEXP
Rcpp::NumericMatrix Zr(Z_);			// creates Rcpp matrix from SEXP
Rcpp::NumericVector tauinv(tauinv_);			// creates Rcpp matrix from SEXP
Rcpp::NumericVector rho(rho_);

GetRNGstate();

int n = Zr.ncol(); 
int k = Ar.ncol(); 
int p = Ar.nrow();

double out = 0.0;

arma::mat Z(Zr.begin(), p, n, false);   	// reuses memory and avoids extra copy
arma::mat A(Ar.begin(), p, k, false);   	// reuses memory and avoids extra copy
arma::mat F(Fr.begin(), k, n, false);   	// reuses memory and avoids extra copy

sampleSparseLoadingsJ(Z, A, F, tauinv, rho, n, p, k );

//arma::mat t;
//arma::mat uvector;
//arma::vec sumf2 = arma::sum(arma::pow(F, 2), 1);

//double u=0.0, uu=0.0, v=0.0, loglog=0.0, logc=0.0;
//for (int i=0; i<p; i++) {
//	for (int j=0; j<k; j++){
//		t = Z.row(i) - A.row(i)*F + A(i,j)*F.row(j);
//		u = arma::accu(F.row(j)%t);
//		v = sumf2(j) + tauinv[i];
//		
//		loglog = log(u*u) - log(2.0*v);
//        
//        logc = log(rho[j]) - log(1-rho[j]) - 0.5*log(v) + 0.5*log(tauinv[i]) + exp(loglog);
        
//        uu = Rf_runif(0.0, 1.0);

//        log(1/uu-1) > logc ? A(i,j) = 0.0 : A(i,j) = Rf_rnorm(u/v, 1/sqrt(v));
//		}
//	}

PutRNGstate();


END_RCPP
}

double rtnorm_1(double mu, double sigma, double a, double b){
	//draw from mean mu, unit variance truncated normal 
	//sigma argument is a farce
	double lo = Rcpp::stats::pnorm_1(a, mu, true, false);
	double hi = Rcpp::stats::pnorm_1(b, mu, true, false);
    double u = Rf_runif(lo, hi);
	return Rcpp::stats::qnorm_1(u, mu, true, false);
}


void sampleZ( Rcpp::NumericMatrix& Z , Rcpp::IntegerMatrix& Ra, Rcpp::NumericMatrix& maxes, Rcpp::IntegerMatrix& argsorts, arma::mat& A, arma::mat& F){

int cn;
arma::mat mu = A*F;
//Rcpp::NumericMatrix mu(mu0);

int n = Ra.ncol(), p = Ra.nrow();

double etol = 1e-13;
double cur_max, cur_min, next_min;
int k, first_rank, cur_rank;
double cm = -7.5;

//impute z
for (int i = 0; i < p; i++) {
	first_rank = 1;
	cur_min = -7.5;
	next_min = -7.5;
	cur_max = maxes(i, 0);
	cur_rank = 0;
	for (int j = 0; j < n; j++) {
		k = argsorts(i,j);
		if (Ra(i,k) < 0) {
			Z(i, k) = Rf_rnorm(mu(i,k), 1.0);
			cn = 1;
		} else if (Ra(i,k)==0 && first_rank==1) {
			Z(i,k) = rtnorm_1(mu(i,k), 1.0, cur_min, cur_max);
			first_rank = 0;
			next_min = Z(i,k);
			cn = 2;
		} else if (cur_rank==Ra(i,k) and cur_rank==0){
			Z(i,k) = rtnorm_1(mu(i,k), 1.0, cur_min, cur_max);
			next_min = std::max(Z(i,k), next_min);
			cn = 3;
		} else if (cur_rank==Ra(i,k)) {
			Z(i,k) = rtnorm_1(mu(i,k), 1.0, cur_min, cur_max);
			next_min = std::max(Z(i,k), next_min);
			maxes(i, (cur_rank-1) ) = std::min(Z(i,k), maxes(i, (cur_rank-1) ));
			cn = 4;
		} else {
			cur_rank += 1;
			cur_min = next_min;
			cur_max = maxes(i, cur_rank);
			Z(i,k) = rtnorm_1(mu(i,k), 1.0, cur_min, cur_max);
			next_min = std::max(Z(i,k), next_min);
			maxes(i, (cur_rank-1) ) = Z(i,k);
			cn = 5;
		}
	}
}

}


SEXP MCMCstep( SEXP Z_, SEXP A_, SEXP F_, SEXP tauinv_, SEXP rho_ , SEXP Ra_, SEXP maxes_, 
			   SEXP argsorts_, SEXP priors_, SEXP nsim_, SEXP nburn_, SEXP thin_, SEXP printstatus_, 
			   SEXP keepscores_, SEXP keeploadings_){
BEGIN_RCPP
	
//Rcpp::List model(model_) ;

GetRNGstate();

bool keepscores   = as<bool>(keepscores_);
bool keeploadings = as<bool>(keeploadings_);

int printstatus = as<int>(printstatus_);

Rcpp::IntegerMatrix Ra(Ra_);
Rcpp::IntegerMatrix argsorts(argsorts_);
Rcpp::NumericMatrix maxes(maxes_);
Rcpp::NumericMatrix Fr(F_);			
Rcpp::NumericMatrix Ar(A_);		
Rcpp::NumericMatrix Zr(Z_);	
Rcpp::NumericVector tauinv(tauinv_);
Rcpp::NumericVector rho(rho_);

int nsim = as<int>(nsim_);
int nburn = as<int>(nburn_);
int thin = as<int>(thin_);

Rcpp::List priors(priors_);
double taua = as<double>(priors["taua"]);
double taub = as<double>(priors["taub"]);
double rhoa = as<double>(priors["rhoa"]);
double rhob = as<double>(priors["rhob"]);
//fix for debugging
//priors = list(taua=1, taub=1, rhoa=100000, rhob=0.1)
//double taua=1.0;
//double taub=1.0;
//double rhoa=100000.0;
//double rhob=0.1;

int n = Zr.ncol(); 
int k = Ar.ncol(); 
int p = Ar.nrow();

arma::mat Z(Zr.begin(), p, n, false);   	// reuses memory and avoids extra copy
arma::mat A(Ar.begin(), p, k, false);   	// reuses memory and avoids extra copy
arma::mat F(Fr.begin(), k, n, false);   	// reuses memory and avoids extra copy

arma::colvec Aravel(Ar.begin(), p*k, false);
arma::running_stat_vec<double> Apstat;
arma::colvec Fravel(Fr.begin(), k*n, false);
arma::running_stat_vec<double> Fpstat;

int samples = nsim/thin;
int isamp = 0;

mat_iterator iF = F.begin();
mat_iterator iA = A.begin();

int score_size = keepscores ? samples*k*n : 1;
Rcpp::NumericVector Fp(score_size);
vec_iterator iFp = Fp.begin();

int load_size = keeploadings ? samples*k*p : 1;
Rcpp::NumericVector Ap(load_size);
vec_iterator iAp = Ap.begin();

Rcpp::NumericVector scAp(k*p);
vec_iterator iscAp = scAp.begin();
arma::colvec armascAp(scAp.begin(), p*k, false);

Rprintf("Beginning MCMC...");

for (int i=0; i<nburn+nsim; i++) {
	sampleSparseLoadingsJ(Z, A, F, tauinv, rho, n, p, k );
	sampleScores(Z, A, F, n, k) ;
	sampleLoadingsVarJ(tauinv, Ar, taua, taub) ;
	sampleRho(rho, Ar, rhoa, rhob);
	sampleZ(Zr, Ra, maxes, argsorts, A, F);
	
	if (i%printstatus==0) {
		Rprintf("iteration %d\n",i);
	}
	
	if (i>=nburn && i%thin==0){
		
		if (keepscores) {
			for (int j = 0; j<k*n; j++) {
				iFp[isamp*k*n + j] = iF[j];
			}
		}
		Fpstat(Fravel);
		
		//cacluate scaling matrix
		arma::colvec b = 1 + arma::sum(arma::pow(A,2), 1);

		//save scaled loadings
		if (keeploadings) {
			for (int j = 0; j<k*p; j++) {
				iAp[isamp*k*p + j] = iA[j]/sqrt(b[j%p]);
				iscAp[j] = iA[j]/sqrt(b[j%p]);
			}
		} else {
			for (int j = 0; j<k*p; j++) {
				iscAp[j] = iA[j]/sqrt(b[j%p]);
			}
		}
		Apstat(armascAp);
		isamp++;
	}
}

if (!keepscores)   { Fp[0] = NA_REAL; }
if (!keeploadings) { Ap[0] = NA_REAL; }

return List::create(
	_["Ap"] = Ap,
	_["Ap.mean"] = Apstat.mean(),
	_["Ap.var"] = Apstat.var(),
	_["Fp"] = Fp,
	_["Fp.mean"] = Fpstat.mean(),
	_["Fp.var"] = Fpstat.var());

PutRNGstate();

END_RCPP
}


