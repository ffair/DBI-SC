#include <Rcpp.h>
#include<algorithm>
using namespace Rcpp;
using namespace std;

inline NumericVector two_co(IntegerVector x,IntegerVector y) {
  int xm = 2;
  int ym = 2;
  
  NumericVector value(ym);
  
  NumericMatrix re_cod(xm,ym);
  
  for (int j = 0; j < xm; j++){
    IntegerVector temp;
    for (size_t i = 0; i < x.size(); ++i ){
      if (x[i]==j){
        temp.push_back(y[i]);
      }
    }
    
    int n = temp.size();
    NumericVector counts(ym);
    for (int i = 0; i < n; i++) {
      counts[temp[i]]++;
    }
    
    for (int i = 0; i < value.size(); i++ ){
      value[i] = counts[i]/n;
    }
    
    re_cod.row(j) = value;
    
  }
  
  NumericVector max_ra(ym);
  for (int i=0 ; i < ym; i++){
    max_ra[i] = Rcpp::max(re_cod.column(i));
  }
  
  return max_ra;
}

// [[Rcpp::export]]
List dist_co(IntegerMatrix m) {
  int num = m.ncol();
  NumericMatrix result(num,num);
  for (int i = 0; i < (num-1); i++){
    for (int j = (i+1); j < num; j++){
      NumericVector temp_vec;
      temp_vec = two_co(m.column(i),m.column(j));
      result(i,j) = sum(temp_vec);
      result(j,i) = result(i,j);
    }
  }
  
  NumericVector dpd(num);
  for (int i = 0; i < num; i++){
    dpd[i] = sum(result.row(i)) - num + 1;
  }
  
  int n = m.nrow();
  NumericMatrix dist(n,n);
  for (int i = 0; i < (n-1); i++) {
    for (int j = (i+1); j < n; j++) {
      NumericVector aka = abs(m.row(i)-m.row(j));
      dist(i,j) = sum(aka * dpd);
      dist(j,i) = dist(i,j);
    }
  }
  
  List ret;
  ret["d"] = dist;
  ret["v"] = dpd;
  ret["vm"] = result;

  
  return ret;
}





