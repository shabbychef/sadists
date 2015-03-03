#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector ddnf(NumericVector x, double df1, double df2, double ncp1, double ncp2) {
	// based on algorithm 10.8 of Paolella
	int n = x.size();
	NumericVector pdf(n);



	return pdf;
}


//for vim modeline: (do not edit)
// vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=cpp:ft=cpp
