#include<Rcpp.h>
#include<RcppEigen.h>
#include<Eigen/Core>
#ifdef RASTERFLOW_USECL
#define __CL_ENABLE_EXCEPTIONS
#include<CL/cl.hpp>
#endif 

#ifndef RASTERFLOW_DENSERASTERFLOWSIMULATOR
#define RASTERFLOW_DENSERASTERFLOWSIMULATOR
class DenseRasterFlowSimulator
{
    public:
        DenseRasterFlowSimulator(SEXP heightMatrix,
                                 SEXP membershipMatrix,
                                 SEXP nLevels,
                                 SEXP seed);
        void simulate(int n);
        void simulateOCL(int n);
        Rcpp::IntegerMatrix getContactMatrix(); // square matrix
        Rcpp::IntegerMatrix getDensityMatrix(); // dimension of original raster
        Rcpp::IntegerMatrix getDensityMatrix2(SEXP inIndices);
        ~DenseRasterFlowSimulator();
    private:
        int getDescentDirection(int Xidx, int Yidx);
        int membershipLevels;
        std::mt19937* generator;
        std::vector<std::uniform_int_distribution<int> > dirDists;
        Rcpp::IntegerMatrix directionMatrix;
        Rcpp::IntegerMatrix currentMembership; 
        Rcpp::IntegerMatrix heightMatrix; // from raster object
        Rcpp::IntegerMatrix originalHeightMatrix;
        Rcpp::IntegerMatrix membershipMatrix; // ID matrix

};
#endif
