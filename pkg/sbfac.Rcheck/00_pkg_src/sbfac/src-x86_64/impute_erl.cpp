
// includes from the plugin

#include <Rcpp.h>


#ifndef BEGIN_RCPP
#define BEGIN_RCPP
#endif

#ifndef END_RCPP
#define END_RCPP
#endif

using namespace Rcpp;


// user includes

double rtnorm(double mu, double sigma, double ain, double bin) {
	Rcpp::NumericVector a(1);
	a[0] = ain;
	Rcpp::NumericVector b(1);
	b[0] = bin;
	Rcpp::NumericVector lo = pnorm(a, mu, sigma);
	Rcpp::NumericVector hi = pnorm(b, mu, sigma);
	Rcpp::NumericVector u = runif(1, lo[0], hi[0]);
	Rcpp::NumericVector z = qnorm(u, mu, sigma);
	return Rcpp::as<double>(z);
}

double rtnorm1(double mu, double sigma, double a, double b){
	//draw from mean mu, unit variance truncated normal 
	//sigma argument is a farce
	//cdef double lo = F(a - mu)
	//cdef double hi = F(b - mu) 
	//return mu + Finv(self.runif(lo, hi))
	double lo = Rcpp::stats::pnorm_1(a, mu, true, false);
	double hi = Rcpp::stats::pnorm_1(b, mu, true, false);
    double u = Rf_runif(lo, hi);
	return Rcpp::stats::qnorm_1(u, mu, true, false);
}


// declarations
extern "C" {
SEXP file1c06dac8( SEXP Z0, SEXP Ra0, SEXP maxes0, SEXP argsorts0, SEXP mu0) ;
}

// definition

SEXP file1c06dac8( SEXP Z0, SEXP Ra0, SEXP maxes0, SEXP argsorts0, SEXP mu0 ){
BEGIN_RCPP
int cn;

Rcpp::IntegerMatrix Ra(Ra0);
Rcpp::IntegerMatrix argsorts(argsorts0);
Rcpp::NumericMatrix maxes(maxes0);
Rcpp::NumericMatrix Z(Z0);
Rcpp::NumericMatrix mu(mu0);

int n = Ra.ncol(), p = Ra.nrow();

double etol = 1e-13;
double cur_max, cur_min, next_min;
int k, first_rank, cur_rank;
double cm = -7.5;

GetRNGstate();

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
			Z(i,k) = rtnorm1(mu(i,k), 1.0, cur_min, cur_max);
			first_rank = 0;
			next_min = Z(i,k);
			cn = 2;
		} else if (cur_rank==Ra(i,k) and cur_rank==0){
			Z(i,k) = rtnorm1(mu(i,k), 1.0, cur_min, cur_max);
			next_min = std::max(Z(i,k), next_min);
			cn = 3;
		} else if (cur_rank==Ra(i,k)) {
			Z(i,k) = rtnorm1(mu(i,k), 1.0, cur_min, cur_max);
			next_min = std::max(Z(i,k), next_min);
			maxes(i, (cur_rank-1) ) = std::min(Z(i,k), maxes(i, (cur_rank-1) ));
			cn = 4;
		} else {
			cur_rank += 1;
			cur_min = next_min;
			cur_max = maxes(i, cur_rank);
			Z(i,k) = rtnorm1(mu(i,k), 1.0, cur_min, cur_max);
			next_min = std::max(Z(i,k), next_min);
			maxes(i, (cur_rank-1) ) = Z(i,k);
			cn = 5;
		}
		//if (Z(i,k)>10) {
		//	std::cout << "troubles at " << i << " , " << k <<" ";
		//}
	}
}

PutRNGstate();
END_RCPP
}



