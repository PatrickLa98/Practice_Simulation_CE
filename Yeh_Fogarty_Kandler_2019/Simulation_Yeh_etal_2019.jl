#=
1. Initialize population with random variants at all traits and no links
2. Iterate through timesteps
    * individuals simultanously chose interaction partner 
    * randomly pick a trait to copy from this partner 
    * if trait is linked the focal individuals copies links with propability 1 - b
    * copying successful with propability c (unsuccessful = id keeps own variant)
    * after transmission new links form at rate of association a (between any unlinked traits)
    * each trait can be innovated with probability μ (switch to another variant)
3. Burn in period of 5000 timesteps with unbiased transmission
4. additional 2000 timesteps with:
    * unbiased transmission
    * pay- off biased transmission
    * conformity-biased (only one trait modeled no burn in needed)
        
record variant frequencies at each timestep

=#



using Plots
using Distributions


## 1. Initialize population

## individuals 
n = 10
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

## 2. Iterate through timesteps


        








