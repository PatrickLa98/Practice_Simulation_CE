#=
1. Initiate population with random variants at all traits and no links
2. Initiate array, or list of matrices for links (0 1 storage)
3. Iterate through timesteps
    3.1 individuals simultanously chose interaction partner 
    3.2 randomly pick a trait to copy from this partner 
        3.2.1 copying successful with propability c (unsuccessful = id keeps own variant)
        3.2.2 if trait is linked the focal individuals copies links with propability 1 - b
    3.3 after transmission new links form at rate of association a (between any unlinked traits)
    3.4 each trait can be innovated with probability μ (switch to another variant)
4. Burn in period of 5000 timesteps with unbiased transmission
5. additional 2000 timesteps with:
    * unbiased transmission
    * pay-off biased transmission
    * conformity-biased (only one trait modeled no burn in needed)
        
record variant frequencies at each timestep

Questions: can links only break during the transmission process? (3.2.2)
=#

using Plots
using Distributions

## parameter settings 

## individuals 
n = 10
## number of traits
traits = 5
## variants
variants = [1, 2, 3, 4]
## timesteps
timesteps = 100
## probability of successful copying of variant
c = 0.99



## 1. Initiate population

## store information in matrix
pop = zeros(Int64, n, traits)
## randomly asign variants to traits for n individuals
for i in 1:n 

    pop[i,:] = rand(variants, 5)
        
end

## 2. Initiate array for links

links_array = zeros(Int64, traits, traits, n) ## rows and columns are traits, 3rd dimension are individuals

## 3. Iterate through timesteps

#for i in 1:timesteps

    ids_seq = collect(1:n)
    interacts_with = Int[]

    new_pop = deepcopy(pop)

    #for id in 1:n

        ## 3.1 choose interaction partner
        potential_interaction_partners = setdiff(ids_seq, id) # exclude possibility to interact with oneselve
        push!(interacts_with, rand(potential_interaction_partners, 1)[1])
        
        ## 3.2 choose trait and copy respective variant 
        traits_of_interacting_individual = pop[interacts_with[id],:]
        chosen_trait = rand(1:traits, 1)[1]    
        variant_of_chosen_trait_of_interactant = traits_of_interacting_individual[chosen_trait]
        variant_of_chosen_trait_of_focalid = pop[id, chosen_trait]

            ## 3.2.1 copying either successful (adapt interactants variant) or unsuccesful (keep own variant)
            new_variant = [variant_of_chosen_trait_of_interactant, variant_of_chosen_trait_of_focalid][rand(Categorical([c, (1 - c)]), 1)][1]
            
            ## store new variant
            new_pop[id, chosen_trait] = new_variant
            
            ## 3.2.2 check if variant is linked and copy link with prob 1 - b ( + copy variant of the linked traits)

            links = links_array[chosen_trait,:,interacts_with[id]]

            for l in 1:traits
                
                if links[l] == 1
                    



            end

    #end

    pop = new_pop






#end



        








