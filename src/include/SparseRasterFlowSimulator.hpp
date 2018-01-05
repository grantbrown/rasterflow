#include<Rcpp.h>
#include<RcppEigen.h>
#include<Eigen/Core>
#ifdef RASTERFLOW_USECL
#define __CL_ENABLE_EXCEPTIONS
#include<CL/cl.hpp>
#endif

#ifndef RASTERFLOW_SPARSERASTERFLOWSIMULATOR
#define RASTERFLOW_SPARSERASTERFLOWSIMULATOR
class SparseRasterFlowSimulator
{
    public:
        SparseRasterFlowSimulator(SEXP heightMatrix,
                                  SEXP membershipMatrix,
                                  SEXP startParticles,
                                  SEXP nLevels,
                                  SEXP seed);
        void simulate(int n);
        void simulateOCL(int n);
        void addParticles(SEXP particles);
        Rcpp::IntegerMatrix getContactMatrix(); // square matrix
        Rcpp::IntegerMatrix getDensityMatrix(); // dimension of original raster
        Rcpp::IntegerMatrix getDensityMatrix2(SEXP inIndices);
        ~SparseRasterFlowSimulator();
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
