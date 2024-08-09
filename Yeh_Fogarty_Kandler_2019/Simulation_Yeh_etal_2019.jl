#=
1. Initiate population with random variants at all traits and no links
2. Initiate array, or list of matrices for links (0 1 storage)
3. Iterate through timesteps
    3.1 individuals simultanously chose interaction partner 
    3.2 randomly pick a trait to copy from this partner 
        3.2.1 copying successful with propability c (unsuccessful = id keeps own variant)
        3.2.2 breaks all links associated with trait other than links that are also present in interaction partner
        3.2.3 if trait is linked the focal individuals copies links with propability 1 - b
    3.3 after transmission new links form at rate of association a (between any unlinked traits)
    3.4 each trait can be innovated with probability μ (switch to another variant)
4. Burn in period of 5000 timesteps with unbiased transmission
5. additional 2000 timesteps with:
    * unbiased transmission
    * pay-off biased transmission
    * conformity-biased (only one trait modeled no burn in needed)
        
record variant frequencies at each timestep

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
## probability that link breaks during transmission
b = 0.2
## rate of association (form a new link)
a = 0.1



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
    new_links_array = deepcopy(links_array)
   # for id in 1:n

        ## 3.1 choose interaction partner
        potential_interaction_partners = setdiff(ids_seq, id) # exclude possibility to interact with oneselve
        push!(interacts_with, rand(potential_interaction_partners, 1)[1])
        
        ## 3.2 choose focal trait and copy package
        traits_of_interacting_individual = pop[interacts_with[id],:]
        chosen_trait = rand(1:traits, 1)[1]    
        variant_of_chosen_trait_of_interactant = traits_of_interacting_individual[chosen_trait]
        variant_of_chosen_trait_of_focalid = pop[id, chosen_trait]

            links_interactant = links_array[chosen_trait,:,interacts_with[id]]
            links_focal = links_array[chosen_trait,:, id]


            ## 3.2.1 copying either successful (adapt interactants variant) or unsuccesful (keep own variant)
            new_variant = [variant_of_chosen_trait_of_interactant, variant_of_chosen_trait_of_focalid][rand(Categorical([c, (1 - c)]), 1)][1]
            
            ## store new variant
            new_pop[id, chosen_trait] = new_variant
            
            #for t in 1:traits
              
                ## 3.2.2 if the interactant has a link on the trait to be copied AND the focal id has a link to the trait t AND the interactant   
                ## doesnt have the link to the trait t (but to another trait) --> break the link of the focal individual
                if sum(links_interactant) != 0 && links_focal[t] == 1 && links_interactant[t] == 0

                    ## store information about new link structure
                    new_links_array[chosen_trait, t, id] = 0
                    new_links_array[t, chosen_trait, id] = 0
                end

                ## 3.2.3 check if variant is linked and copy link with prob 1 - b (+ copy variant of the linked traits)
                if links_interactant[t] == 1
                 
                    variant_of_linked_trait_of_interactant = traits_of_interacting_individual[t]
                    variant_of_linked_trait_of_focalid = pop[id, t]

                    ## decide if linked is copied
                    link_copied =  [1, 0][rand(Categorical([(1 - b), b]), 1)][1]
                    ## store information about new link structure
                    new_links_array[chosen_trait, t, id] = link_copied
                    new_links_array[t, chosen_trait, id] = link_copied

                    ## decide if copying is successful (repeat 3.2.1; provided that link doesnt break)
                    if link_copied == 1
                          
                           new_variant = [variant_of_linked_trait_of_interactant, variant_of_linked_trait_of_focalid][rand(Categorical([c, (1 - c)]), 1)][1]
                           new_pop[id, t] = new_variant
                    end

                end

                ## 3.3 new links form with rate of association a (similar to probability)

                links_focal_after_transmission = new_links_array[chosen_trait,:, id]

                for t in 1:traits

                    new_link = 0
                    ## if there exists no prior link, decide if new link is formed
                    if links_focal_after_transmission[t] == 0

                        new_link = [0, 1][rand(Categorical([(1 - a), a]), 1)][1]
                    end

                    ## if link was formed, decide to which trait
                    if new_link == 1

                        newly_linked_trait = rand(1:5, 1)[1]
                        new_links_array[newly_linked_trait, t, id]

                    end

                end

## NOTE :New link formation doesnt work yet I think, potentially issue with copy vs deepcopy and setting vectors locally (in loop) vs globally (outside loop) 

                ## 3.4 each trait can be innovated (switch to random other variant) with probability μ


            #end

   # end

    pop = new_pop
    links_array = new_links_array






#end



        








