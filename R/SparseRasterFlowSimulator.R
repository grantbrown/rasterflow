sparseRasterFlowSimulator = function(heightMatrix, membershipMatrix, n,
                                     particleAdditions,
                                     keepIndices = NA)
{

    if (length(keepIndices) == 1 && all(is.na(keepIndices)))
    {
        print("No keepIndices included")
    }

    uq = unique(as.integer(membershipMatrix))
    if (!(all(uq[order(uq)] == 1:length(uq))))
    {
        stop(paste("The membership matrix should be composed of",
                   "integers numbered from 1 to the number of groups.\n", 
                   sep = " "))
    }
    if (length(particleAdditions) != length(n))
    {
        stop("n must be the same length as the list of particle addition matrices.")
    }
    if (class(particleAdditions) != "list")
    {
        stop("particleAdditions must be a list")
    }
    if (!all(sapply(particleAdditions, class) == "matrix"))
    {
        stop("Each member of particleAdditions shoudl be a matrix")
    }
    len = length(unique(as.integer(membershipMatrix)))
    # Draw random seed using current RNG state
    # Our C++ code is robust to bad seeds. 
    seed = abs(floor(rcauchy(1,location = 1000000, scale = 100000)))
    print(paste("Number of unique values: ", len, sep = ""))
    SRFS = new(SparseRasterFlowSimulator, heightMatrix, membershipMatrix, 
               particleAdditions[[1]],
               len, seed)

    if (length(n) == 1)
    {
        SRFS$simulate(n);
        if (length(keepIndices) == 1 && all(is.na(keepIndices)))
        {
             output1 = SRFS$getDensityMatrix();       
        }
        else
        {
             output1 = SRFS$getDensityMatrix(keepIndices);       
        }

        output2 = SRFS$getContactMatrix();
        output = list(density=output1, 
                      contact = output2)
    }
    else
    {
        itr = 0
        output = lapply(1:length(n), function(x){x})
        for (i in 1:length(n))
        {
            cat(paste("Iteration ", itr, "\n", sep = ""))
            SRFS$simulate(n[i]);
            itr = itr + n[i]
            if (length(keepIndices) == 1 && all(is.na(keepIndices)))
            {
                output[[i]] = list(density=SRFS$getDensityMatrix(),
                                   contact=SRFS$getContactMatrix()); 
            }
            else
            {
                output[[i]] = list(density=SRFS$getDensityMatrix(keepIndices), 
                                   contact=SRFS$getContactMatrix()); 

            }
            if (i != length(n))
            {
                SRFS$addParticles(particleAdditions[[i+1]])
            }
        }
    }
    rm(SRFS)
    return(output)
}
