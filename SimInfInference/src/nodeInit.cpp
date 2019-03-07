#include <Rcpp.h>
#include <math.h>

using namespace Rcpp;

//' Init the nodes in u0.
//' p proportion should be infected.
//' @export
// [[Rcpp::export]]
IntegerVector nodeIniter(IntegerVector s0, double p) {
    double carryflag = 0.0;
    double ipreFlag = 0.0;

    int len = s0.size();
    IntegerVector i = IntegerVector(len);

    for(int j = 0; j < len; j++) {
        ipreFlag = p * s0[j] + carryflag;
        i[j] = round(ipreFlag);
        carryflag = ipreFlag - i[j];
    }

    return(i);
}
