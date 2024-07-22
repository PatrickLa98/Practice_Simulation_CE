## Goal: Model the haploid version of the hitch-hiking effect of a favourable gene
## 1. numerical deterministic simulation and 2. stochastic abm 
#= 
allele b subsitituted by favourable B
neighbouring locus neutral alleles a and A 
closer neighbouring loci are more likely to be "hitchhiked" denoted as distance/recombination fraction c (smaller values mean stronger hitchhiking effects)
B individuals descended from single mutant aB in initial generation
    -> thus the fraction of individuals with AB in the initial generation is Q0 = 0

 The hitchhiking effect acts on allele a since the initial favourable mutation B occurs in a chromosome with alleles aB at neighbouring loci 
 By tracking A we record the decrease in frequency due to the hitchhiking effect
=#

using Plots
using Distributions

cd("/Users/patricklauer/Documents/GitHub/Practice_Simulation_CE/Maynard_Smith&Haigh_1974")
include("Equations_Maynard_Smith&Haigh_1974.jl")

## parameter settings
    ## Population size = n
    N = 10^2
    ## initial frequency of favourable allele B in population of N individuals = p0
    p0 = 1/N
    ## initial proportion of A in chromosomes containing B = Q
    Q0 = 0
    ## initial proportions of A in chromosomes containing b = R
    R0 = 0.5
    ## selection coefficient = s
    s = 0.02
    ## recombination fraction between the two loci = c
    #c = 0.001 


## 1. Numerical simulation to replicate Figure 1
## iterate through different recombination fraction (c) values
recombination_fractions = [c for c in 0:0.0001:0.01]
hitchhiking_effect = []
approx_hitchhiking_effect = []

for c in recombination_fractions

    push!(hitchhiking_effect, Equation8(c, p0, s, R0))
    
    push!(approx_hitchhiking_effect, Equation14(c, s, p0))
    
end

plot(recombination_fractions, hitchhiking_effect,
    ylims = (0,1),
    label = "'exact' approach",
    xlabel = "recombination fraction",
    ylabel = "hitch-hiking effect")
plot!(recombination_fractions, approx_hitchhiking_effect,
    label = "approximation")

## 2. ABM to replicate Figure 1
## transform deterministic system to stochastic system 

#=
* notation: AB = 1, Ab = 2, aB = 3, ab = 4
* starting population: N - 1 individuls sampled as 2 (Ab) or 4 (ab) (proportion dependent on R0 value)
                     0 individuals with genotype 1 (AB)
                     1 indiviudal with genotype 3 (aB)
* individual with genotype 3 will spread and depending on strength selection and of recombination fraction c 1 will also increase
* run till all individuals are either 1 or 3. 
* proportion of predefined R0 to sum 3 individuals at the end is Qinf/R0
=#

hitchhicking_ABM = function (N, R0, c, s)    

    ## number of generations 
    generations = 100
    ## measure if favourable allele B becomes fixed
    fixed = 0
    ## measure the number of replicates need for B to become fixed
    replicate = 0
    ## define result fixed pop outside of while loop (because julia stores it as local object otherwise)
    fixed_pop = []

    ## parent population reproduces by randomly sampling through parent generation and copying genotype (no further mutation)
    ## the neighbouring loci may recombine with recombination fraction c
    ## offspring survives 100 % if B and survives 1 - s if b
    ## if offspring dies a new random parent is sampled


    ## REPLICATE SIMULATION TILL B FIXED
    while fixed == 0

        ## create initial population
        ## asign types 2 (Ab) and 4 (ab) to N-1 individuals proportionally to R0 and 1 - R0
        pop = [2, 4][rand(Categorical([R0, (1 - R0)]), N -1)]
        ## add mutated individual with type 3 (aB)
        push!(pop, 3)

        ## LOOP THROUGH PREDEFINED GENERATIONS (TO ENSURE ABSORPTION STATE IS REACHED) 
        for n in 1:generations

            i = 1
            new_pop = []
        
            ## LOOP THROUGH PARENT-POPULATION TILL THE SAME NUMBER OF OFFSPRING SURVIVES AS PRESENT IN PARENT-POPULATION
            while i <= length(pop)

                offspring = rand(pop) # randomly sample parent from previous pop to reproduce

                ## RECOMBINATION
                ## offsprings genotype can recombine depending on recombination fraction c
                ## for simplicity and since there is selection for B, only allow recombination of A and a (hence + or - 2 jumps)
                if offspring in [1, 2]

                    offspring = [offspring, offspring + 2][rand(Categorical([(1 - c), c]), 1)][1]

                elseif offspring in [3, 4]
                    
                    offspring = [offspring, offspring - 2][rand(Categorical([(1 - c), c]), 1)][1]

                end

                ## SELECTION
                if offspring in [2, 4]  # offspring has change of dying if unfavourable allele b

                    survives = [0, 1][rand(Categorical([s, (1 - s)]), 1)]

                    if survives == [1]
                        
                        push!(new_pop, offspring)
                        i += 1
                    elseif survives == [0]
                        continue
                    end
                elseif offspring in [1, 3]

                    push!(new_pop, offspring)
                    i += 1
                end
            end

            ## replace parent gen
            pop = new_pop

        end

        replicate += 1
        ## check if favourable allele B is fixed in population
        if sum(pop .== 1) + sum(pop .== 3) == N

            fixed = 1
            fixed_pop = push!(fixed_pop, pop)
        end
    end
    return fixed_pop
end


test = hitchhicking_ABM(N, R0, c , s)

sum(test .== 1)