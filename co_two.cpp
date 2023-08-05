#include <Rcpp.h>
#include<algorithm>
using namespace Rcpp;
using namespace std;

inline Rcpp::IntegerVector from_std_vec(std::vector<int> v) {
  return Rcpp::IntegerVector(v.begin(), v.end());
}


// [[Rcpp::export]]
NumericVector two_co(IntegerVector x,IntegerVector y) {
  int xm = Rcpp::max(x);
  int ym = Rcpp::max(y);
  
  NumericVector value((ym+1));
  
  NumericMatrix re_cod((xm+1),(ym+1));
  
  for (int j = 0; j < (xm+1); j++){
    IntegerVector temp;
    for (size_t i = 0; i < x.size(); ++i ){
      if (x[i]==(j)){
        temp.push_back(y[i]);
      }
    }
    
    int n = temp.size();
    NumericVector counts(ym);
    for (int i = 0; i < n; i++) {
      counts[(temp[i])]++;
    }
    
    for (int i = 0; i < value.size(); i++ ){
      value[i] = counts[i]/n;
    }
    
    re_cod.row(j) = value;
    
  }
  
  NumericVector max_ra((ym+1));
  for (int i=0 ; i < (ym+1); i++){
    max_ra[i] = Rcpp::max(re_cod.column(i));
  }
  
  return max_ra;
}



/*** R
two_co(c(0,0,1,1,0),c(0,0,1,1,1))
*/
