#include <SparseRasterFlowSimulator.hpp>
#include <DenseRasterFlowSimulator.hpp>
#include <random>


static int dirRowIdx[9]= {0,-1,-1,-1,0,1,1,1,0};
static int dirColIdx[9] = {0,-1,0,1,1,1,0,-1,-1};




SparseRasterFlowSimulator::SparseRasterFlowSimulator(SEXP inputHeightMatrix,
                                                     SEXP inputMembershipMatrix,
                                                     SEXP startParticles,
                                                     SEXP inLocations,
                                                     SEXP seed)
{
    // Store image data    
    int idx, i, j;
    int random_seed = Rcpp::IntegerVector(seed)(0);
    heightMatrix = Rcpp::IntegerMatrix(inputHeightMatrix); 
    originalHeightMatrix = Rcpp::clone(heightMatrix);
    membershipMatrix = Rcpp::IntegerMatrix(inputMembershipMatrix); 
    membershipLevels = (Rcpp::IntegerVector(inLocations))(0);

    Rcpp::IntegerMatrix particleMatrix = Rcpp::IntegerMatrix(startParticles);

    // Cache random number generators which we'll use later
    dirDists = std::vector<std::uniform_int_distribution<int> >();
    for (i = 0; i < 9; i++)
    {
        dirDists.push_back(std::uniform_int_distribution<int>(0,i));
    }



    if (heightMatrix.ncol() != membershipMatrix.ncol() || 
        heightMatrix.nrow() != membershipMatrix.nrow())
    {
        Rcpp::stop("Height matrix and membership matrix must be the same dimension.");
    }

    directionMatrix = Rcpp::IntegerMatrix(heightMatrix.nrow(), heightMatrix.ncol());
    // Columns: 
    // 1. starting row
    // 2. starting col
    // 3. current row
    // 4. current col
    // 5. group number
    currentMembership = Rcpp::IntegerMatrix(particleMatrix.nrow(), 
                                            5);    
    for (i = 0; i < particleMatrix.nrow(); i++)
    {
        currentMembership(i, 0) = particleMatrix(i,0);
        currentMembership(i, 1) = particleMatrix(i,1);
        currentMembership(i, 2) = particleMatrix(i,0);
        currentMembership(i, 3) = particleMatrix(i,1);
        currentMembership(i, 4) = membershipMatrix(particleMatrix(i,0),
                                                   particleMatrix(i,1));
    }

    // Set up random number generator
    std::minstd_rand0 lc_generator(random_seed);
    std::uint_least32_t seed_data[std::mt19937::state_size];
    std::generate_n(seed_data, std::mt19937::state_size, std::ref(lc_generator));
    std::seed_seq q(std::begin(seed_data), std::end(seed_data));
    generator = new std::mt19937{q};   
}

SparseRasterFlowSimulator::~SparseRasterFlowSimulator()
{
    delete generator;
}

int SparseRasterFlowSimulator::getDescentDirection(int Xidx, int Yidx)
{
    int dir,l;
    // Make sure we're on the interior

    if (Xidx == 0 
            || Yidx == 0 
            || Xidx == heightMatrix.nrow() - 1
            || Yidx == heightMatrix.ncol() - 1)
    {
        return(0);
    }
    
    // Water takes up space.
    std::vector<int> descentDirections;

    // Get direction and set new vals
    descentDirections.clear();
    dir = 0;
    int maxDrop = -1;
    int currentHt = heightMatrix(Xidx,Yidx) + originalHeightMatrix(Xidx,Yidx);
    int htDrops[9] = {
        0,
        currentHt - (heightMatrix(Xidx-1, Yidx-1) + originalHeightMatrix(Xidx - 1, Yidx - 1)),    // up, left
        currentHt - (heightMatrix(Xidx-1, Yidx) + originalHeightMatrix(Xidx - 1, Yidx)),      // up
        currentHt - (heightMatrix(Xidx-1, Yidx+1) + originalHeightMatrix(Xidx - 1, Yidx + 1)),    // up, right
        currentHt - (heightMatrix(Xidx  , Yidx+1) + originalHeightMatrix(Xidx, Yidx + 1)),    // right
        currentHt - (heightMatrix(Xidx+1, Yidx+1) + originalHeightMatrix(Xidx + 1, Yidx + 1)),    // down, right
        currentHt - (heightMatrix(Xidx+1, Yidx) + originalHeightMatrix(Xidx + 1, Yidx )),      // down
        currentHt - (heightMatrix(Xidx+1, Yidx-1) + originalHeightMatrix(Xidx + 1, Yidx - 1)),    // down, left
        currentHt - (heightMatrix(Xidx  , Yidx-1) + originalHeightMatrix(Xidx, Yidx - 1))};   // left
    for (l = 0; l < 9; l++)
    {
        maxDrop = std::max(maxDrop, htDrops[l]);
    }
    for (l = 0; l < 9; l++)
    {
        if (htDrops[l] == maxDrop)
        {
            descentDirections.push_back(l);
        }
    }
    if (descentDirections.size() == 0)
    {
        Rcpp::stop("No valid descent directons found. This shouldn't happen!\n");
    }

    dir = descentDirections[dirDists[descentDirections.size() - 1](*generator)];
    return(dir);
}

void SparseRasterFlowSimulator::addParticles(SEXP particles)
{
    Rcpp::IntegerMatrix particleMatrix = Rcpp::IntegerMatrix(particles);
    Rcpp::IntegerMatrix newCurrentMembership = Rcpp::IntegerMatrix(
                                                    currentMembership.nrow() 
                                                  + particleMatrix.nrow(), 
                                                  5);    
    int idx = 0;
    int i,j;
    for (i = 0; i < currentMembership.nrow(); i++)
    {
        for (j = 0; j < currentMembership.ncol(); j++)
        {
            newCurrentMembership(i,j) = currentMembership(i,j);
        }
        idx++;
    }
    for (i = 0; i < particleMatrix.nrow(); i++)
    {
        newCurrentMembership(idx, 0) = particleMatrix(i,0);
        newCurrentMembership(idx, 1) = particleMatrix(i,1);
        newCurrentMembership(idx, 2) = particleMatrix(i,0);
        newCurrentMembership(idx, 3) = particleMatrix(i,1);
        newCurrentMembership(idx, 4) = membershipMatrix(particleMatrix(i,0),
                                                     particleMatrix(i,1));
        idx++;
    }
    currentMembership = newCurrentMembership;
}

void SparseRasterFlowSimulator::simulate(int n)
{
    int idx, itr, i, j, dir, Xidx, Yidx;
    // Get starting height matrix
    heightMatrix = getDensityMatrix();
    for (itr = 0; itr < n; itr++)
    {
        for (idx = 0; idx < currentMembership.nrow(); idx++)
        {
            if (currentMembership(idx,4) != 1){
                Xidx = currentMembership(idx,2);
                Yidx = currentMembership(idx,3);
                dir = getDescentDirection(Xidx,Yidx);

                currentMembership(idx,2) += dirRowIdx[dir];
                currentMembership(idx,3) += dirColIdx[dir];
                heightMatrix(Xidx, Yidx) -= 1;
                // Might cause water to pile up at shoreline; shouldn't be an 
                // issue for Cholera models. 
                // Better would be to check if next location is water, like
                // getDensityMatrix does, but this should be good enough. 
                heightMatrix(currentMembership(idx, 2), 
                        currentMembership(idx, 3)) += 1; 
            }
        }
    }
}

void SparseRasterFlowSimulator::simulateOCL(int n)
{
    // Prepare OpenCL Environment
    std::vector<cl::Platform> platformVector;
    cl::Platform::get(&platformVector);

    try
    {
        std::vector<std::vector<cl::Device> > platformCPUDevices;
        std::vector<std::vector<cl::Device> > platformGPUDevices;
        std::vector<cl::Device> tmpCPUDeviceVector;
        std::vector<cl::Device> tmpGPUDeviceVector;

        for (int i = 0; i < platformVector.size(); i++)
        {
            tmpGPUDeviceVector = std::vector<cl::Device>();
            tmpCPUDeviceVector = std::vector<cl::Device>();
            try
            {
                platformVector[i].getDevices(CL_DEVICE_TYPE_CPU, 
                        &tmpCPUDeviceVector);
                platformCPUDevices.push_back(tmpCPUDeviceVector);
            }
            catch(cl::Error e)
            {
                if (e.err() == -1)
                {
                    Rcpp::Rcout << "Note: no supported OpenCL CPU devices found on system.\n";
                }
                else
                {
                    throw(e);
                }
            }
            try
            {  
                platformVector[i].getDevices(CL_DEVICE_TYPE_GPU, 
                        &tmpGPUDeviceVector);
                platformGPUDevices.push_back(tmpGPUDeviceVector);

            }
            catch(cl::Error e)
            {
                if (e.err() == -1)
                {
                    Rcpp::Rcout << "Note: no supported OpenCL GPU devices found on system.\n";
                }
                else
                {
                    throw(e);
                }
            }

        }
        if (platformVector.size() == 0)
        {
            Rcpp::Rcout << "Note: no supported OpenCL platforms detected.\n";
        }

        std::string platformName;
        for (int i = 0; i < platformVector.size(); i++)
        {
            Rcpp::Rcout << "Platform " << i << ": ";
            platformVector[i].getInfo(CL_PLATFORM_NAME, &platformName);
            Rcpp::Rcout << platformName << "\n";

        }

    }
    catch(cl::Error e)
    {
        Rcpp::Rcout << "CL Error Encountered\n";
        Rcpp::Rcout << e.what() << "\n";
        Rcpp::Rcout << e.err() << "\n";
        Rcpp::stop("exiting...");
    }
    Rcpp::stop("Hey, I haven't finished the OpenCL version yet. Hold your horses.");
}


Rcpp::IntegerMatrix SparseRasterFlowSimulator::getContactMatrix()
{
    Rcpp::IntegerMatrix output(membershipLevels, membershipLevels);
    int i;
    // Columns: 
    // 1. starting row
    // 2. starting col
    // 3. current row
    // 4. current col
    // 5. group number
    int sourceMemberNum;
    int targetMemberNum;
    int Xidx, Yidx, cidx;
    for (i = 0; i < currentMembership.nrow(); i++)
    {
        Xidx = currentMembership(i,2);
        Yidx = currentMembership(i,3);
       
        sourceMemberNum = currentMembership(i,4);    

        targetMemberNum = membershipMatrix(Xidx, Yidx);
        output(sourceMemberNum - 1, targetMemberNum -1) += 1;
    }
    return(output);
}

Rcpp::IntegerMatrix SparseRasterFlowSimulator::getDensityMatrix()
{
    Rcpp::Rcout << "Calling gdm1\n";
    Rcpp::IntegerMatrix outputMatrix(heightMatrix.nrow(), heightMatrix.ncol());
    int i, Xidx, Yidx;
    for (i = 0; i < currentMembership.nrow(); i++)
    {
        Xidx = currentMembership(i,2);
        Yidx = currentMembership(i,3);
        // Ocean is always zero density
        if (currentMembership(i,4) != 1)
        {
            outputMatrix(Xidx,Yidx) += 1;
        }
    }
    return(outputMatrix);
}

Rcpp::IntegerMatrix SparseRasterFlowSimulator::getDensityMatrix2(SEXP inIndices)
{
    Rcpp::IntegerVector indices(inIndices);
    Rcpp::Rcout << "Calling gdm2\n";
    Rcpp::Rcout << indices.size() << "\n";
    Rcpp::Rcout << indices.length() << "\n";
    Rcpp::IntegerMatrix outputMatrix(heightMatrix.nrow(), heightMatrix.ncol());
    int i, Xidx, Yidx;
    for (i = 0; i < indices.length(); i++)
    {
        Xidx = currentMembership(indices(i),2);
        Yidx = currentMembership(indices(i),3);
        outputMatrix(Xidx,Yidx) += 1;
    }
    return(outputMatrix);

}



RCPP_MODULE(SparseRasterFlowSimulator){
    using namespace Rcpp ;
    class_<SparseRasterFlowSimulator>( "SparseRasterFlowSimulator" )
    .constructor<SEXP, SEXP, SEXP, SEXP, SEXP>()
    .method("simulate", &SparseRasterFlowSimulator::simulate)
    .method("simulateOCL", &SparseRasterFlowSimulator::simulate)
    .method( "getContactMatrix", &SparseRasterFlowSimulator::getContactMatrix )
    .method( "getDensityMatrix", &SparseRasterFlowSimulator::getDensityMatrix2 )
    .method( "getDensityMatrix", &SparseRasterFlowSimulator::getDensityMatrix )
    .method( "addParticles", &SparseRasterFlowSimulator::addParticles );


}
