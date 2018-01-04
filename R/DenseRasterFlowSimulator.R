denseRasterFlowSimulator = function(heightMatrix, membershipMatrix, n, 
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
    len = length(unique(as.integer(membershipMatrix)))
    # Draw random seed using current RNG state
    # Our C++ code is robust to bad seeds. 
    seed = abs(floor(rcauchy(1,location = 1000000, scale = 100000)))
    print(paste("Number of unique values: ", len, sep = ""))
    DRFS = new(DenseRasterFlowSimulator, heightMatrix, membershipMatrix, 
               len, seed)

    if (length(n) == 1)
    {
        DRFS$simulate(n);
        if (length(keepIndices) == 1 && all(is.na(keepIndices)))
        {
             output1 = DRFS$getDensityMatrix();       
        }
        else
        {
             output1 = DRFS$getDensityMatrix(keepIndices);       
        }

        output2 = DRFS$getContactMatrix();
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
            DRFS$simulate(n[i]);
            itr = itr + n[i]
            if (length(keepIndices) == 1 && all(is.na(keepIndices)))
            {
                output[[i]] = list(density=DRFS$getDensityMatrix(),
                                   contact=DRFS$getContactMatrix()); 
            }
            else
            {
                output[[i]] = list(density=DRFS$getDensityMatrix(keepIndices), 
                                   contact=DRFS$getContactMatrix()); 

            }


        }
    }
    rm(DRFS)
    return(output)
}
