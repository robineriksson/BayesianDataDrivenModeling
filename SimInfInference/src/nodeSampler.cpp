//#include <Rcpp.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <cmath>
#include <algorithm>

using namespace Rcpp;

//' set the seed
//' @export
// [[Rcpp::export]]
void set_seed(unsigned int seed) {
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(seed);
}



//' draw random sample
//' @export
// [[Rcpp::export]]
NumericVector csample_num( NumericVector x,
                           int size,
                           bool replace,
                           NumericVector prob = NumericVector::create()
                           ) {
  NumericVector ret = RcppArmadillo::sample(x, size, replace, prob);
  return ret;
}

//' a vector shuffle
void myShuffle(int *array, int stop) {

    NumericVector x(stop);
    for(int i = 0; i < stop; i++)
        x(i) = i;

    NumericVector order = csample_num(x, stop, false);

    int shuffledArray[stop];

    int j;
    for(int i = 0; i < stop; i++) {
        j = (int)order(i);
        // Rprintf("j: %d\n",j);
        // Rprintf("*(array + j) %d\n",*(array + j));
        shuffledArray[i] = *(array + j);
    }

    for(int i = 0; i < stop; i++) {
        *(array + i) = *(shuffledArray + i);
        // Rprintf("%d ", array[i]);
    }
}

//' a tester for the shuffle
//' @export
// [[Rcpp::export]]
IntegerVector testShuffle(IntegerVector array) {
    int size = array.size();
    int inArray[size];
    IntegerVector outArray(size);

    // from IntegerVector to *int
    for(int i = 0; i < size; i++)
        inArray[i] = array(i);


    myShuffle(inArray, size);

    // from *int to IntegerVector
    for(int i = 0; i < size; i++)
        outArray(i) = inArray[i];

    return(outArray);
}

//' Sample the pools
//' @export
// [[Rcpp::export]]
IntegerVector sample_pools_Rcpp(int S, int I){
    if((S + I) < 1)
        return 0;
    int cols = 3;

    while( (S + I) % cols != 0)
        S += 1;

    int rows = (S+I)/cols;

    int siVec[S+I];
    for (int i = 0; i < S; i++)
        siVec[i] = 0;
    for (int i = S; i < S+I; i++)
        siVec[i] = 1;


    // randomly shuffly the vector
    // must be 1 out of the array.
    // std::random_shuffle(&siVec[0], &siVec[S+I]);
    int stop = S+I;
    myShuffle(siVec, stop);

    // Create pools 3 and 3 from the vector of individuals. Create
    // a matrix nrow x 3 where each row is one pool. Ignore that
    // the number of individuals are not always divisble by three.
    IntegerMatrix mMat(rows,cols);

    IntegerVector rowSum(rows);

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            mMat[i,j] = siVec[cols*i+j];
            rowSum[i] += mMat[i,j];
        }
    }

    // Randomly determine the status of each pool based on the number
    // of positive samples in each pool and the sensitivity.
    NumericVector rands = runif(rows);
    float sensitivity [4] = {0, 0.2775461, 0.8552848, 0.9891212};
    IntegerVector status(S+I);
    int pos = 0;
    for (int i = 0; i < rows; i++) {
        //Rprintf("rowSum:%d,\t", rowSum[i]);
        //Rprintf("sensitivity:%f,\t", sensitivity[rowSum[i]]);
        //Rprintf("rand:%f,\t", rands[i]);
        status[i] = rands[i] < sensitivity[rowSum[i]];
        //Rprintf("status:%d\n", status[i]);
        pos += status[i];
    }

    int n = rows;

    IntegerVector out = IntegerVector::create(pos,
                                              n);
    return out;
}

//' Predict the outcome of the env. sample
//' @export
// [[Rcpp::export]]
bool predict_env_sample_SISe_Rcpp(float whpp) {
    float x0 = -382.32431;
    float x1 = 91.49417;
    float a = exp(x0 + x1*whpp);
    float prob = a / (1+a);

    if(std::isnan(prob))
        prob = 1.0;

    float u = R::runif(0,1);

    return u < prob;
}

//' Sample the herd
//' @export
// [[Rcpp::export]]
bool sample_herd_Rcpp(int s, int i) {
    // Sample pools
    IntegerVector pools = sample_pools_Rcpp(s,i);

    // Determine within-group pool
    float denom = 0.0;
    if(pools[1] > 0)
        denom = pools[1];
    else
        denom = 1;

    float whpp = 100 * pools[0] / denom;


    // predict the outcome of the sample
    bool p = predict_env_sample_SISe_Rcpp(whpp);

    return p;
}

//' Sample a node
//' @export
// [[Rcpp::export]]
LogicalVector sample_single_node_Rcpp(IntegerVector svec, IntegerVector ivec) {
    int s,i;
    int size = svec.size();
    LogicalVector out(size);
    for(int j = 0; j < svec.size(); j++) {
        s = svec[j];
        i = ivec[j];
        out[j] = sample_herd_Rcpp(s,i);
    }
    return out;
}

//' Sample multiple nodes
//' @export
// [[Rcpp::export]]
LogicalVector sample_nodes_Rcpp(IntegerMatrix resMat, int seed) { //, int seed){
    //std::default_random_engine generator(seed);
    set_seed(seed);
    //RNGScope scope;
    int numComp = 2;
    int ncol = resMat.ncol();
    int nrow = resMat.nrow();
    //Rprintf("col:%d\n", ncol);
    //Rprintf("row:%d\n", nrow);

    //IntegerVector svec(ncol);
    //IntegerVector ivec(ncol);
    int ss, ii;
    //LogicalVector sampleVecNode(ncol);
    LogicalVector sampleOut(nrow/numComp*ncol);

    for(int i = 0; i < nrow; i+=numComp) {
        for(int j = 0; j < ncol; j++) {
            //svec(j) = resMat(i,j);
            //ivec(j) = resMat(i+1,j);
            ss = resMat(i,j);
            //Rprintf("ss: %d\t", ss);
            ii = resMat((i+1),j);
            //Rprintf("ii: %d\t", ii);
            bool foo = sample_herd_Rcpp(ss,ii);
            //Rprintf("%d\n",foo);
            sampleOut(i/numComp*ncol + j) = foo;
        }
        //sampleVecNode = sample_single_node_Rcpp(svec, ivec);
        //for(int j = 0; j < ncol; j++) {
        //    sampleOut(i/numComp*ncol + j) = sampleVecNode(j);
        //}
    }

    return sampleOut;

}


//' Helper function for timeMatch.
//' Find the last entry in the vector of the same entry as
//' in the "starting" position.
//' @export
// [[Rcpp::export]]
int findLast(IntegerMatrix obs, int start, int node){
    bool foundLast = false;
    int len = obs.nrow();
    int last = len;

    int i = start;
// obs is: [node, time, sample]
    while(!foundLast){
        if(obs(i,0) != node){
            foundLast = true;
            last = i;
        }

        i += 1;
        if(i == len)
            foundLast = true;
    }
    return last;
}


//' match element that we should extract.
//' @export
// [[Rcpp::export]]
IntegerVector timeMatch(IntegerMatrix sim, IntegerMatrix obs, IntegerVector nodes){
    int numNodes = nodes.size();
    int lenTime = sim.nrow()/numNodes; // equal length of data for all nodes.
    int lenObsTime = obs.ncol();

    int start = 0;
    int stop = lenTime;

    IntegerVector elemMatch(lenObsTime);
    IntegerVector simNode(lenTime);
    IntegerVector simTime(lenTime);
    // ["time", "node", ...]
    for(int i = 0; i < numNodes; i++) {
        int node = nodes[i];

        // Get the elements for which the node is nObs[i] <-- IN ORDER!
        // sim <- traj.df[traj.df$node == nObs[i],]
        for(int j = 0; j < lenTime; j++) {
            simNode[j] = sim(j, 0);
        }

        // Get the elements for which the node is nObs[i] <-- un-regular length
        // obs <- obsDates[obsDates$node == nObs[i],]
        stop = findLast(obs, start, node);
        IntegerVector obsTime(stop-start);
        for(int j = 0; j < obsTime.size(); j++)
            obsTime[j] = obs[start + j];
        start = stop + 1;


    }


    // jOld <- 1
    // jNew <- 0
    // for(i in 1:length(nObs)){
    //     sim <- traj.df[traj.df$node == nObs[i],]
    //     obs <- obsDates[obsDates$node == nObs[i],]
    //     found <- which(sim$time %in% obs$numTime)
    //     jNew <- jOld + length(found)
    //     ##if(i == 1) jOld <- 1;
    //     obsTimeElem[jOld:(jNew-1)] <- length(tObs)*(i-1) + found
    //     jOld <- jNew
    // }
    return 0;
}
