using Plots
using Distributions


## set up

    ## individuals 
    n = 2
    ## traits
    traits = [1, 2, 3, 4, 5]
    ## variants
    variants = [1, 2, 3, 4]

    ## store information in matrix
    pop = zeros(Int64, n, length(traits))
    ## randomly asign variants to traits for n individuals
    for i in 1:n 

        pop[i,:] = rand(variants, 5)
    
    end

    


