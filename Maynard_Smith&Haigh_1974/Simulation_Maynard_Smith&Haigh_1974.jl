## Goal: Model the haploid version of the hitch-hiking effect of a favourable gene
## numerical simulation and abm 
#= 
allele b subsitituted by favourable B
neighbouring locus neutral alleles a and A 
closer neighbouring loci are more likely to be "hitchhiked" denoted as distance/recombination fraction c (smaller values mean stronger hitchhiking effects)
B individuals descended from single mutant aB in initial generation
    -> thus the fraction of individuals with AB in the initial generation is Q0 = 0
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


## Numerical simulation to replicate Figure 1
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

## ABM to replicate Figure 1

#=
* notation: AB = 1, Ab = 2, aB = 3, ab = 4
* starting population: N - 1 individuls sampled as 2 or 4 (proportion dependent on R0 value)
                     0 individuals 1
                     1 indiviudal 3
* individual with genotype 3 will spread and depending on strength of recombination fraction c 1 will also increase
* run till all individuals are either 1 or 3. 
* proportion of predefined R0 to sum 3 individuals at the end is Qinf/R0
=#

## population of N - 1 individuals with b allele(denotetd as 0)        
pop = repeat([0], N - 1)
## 1 individual mutates and recieves B allele
push!(pop, 1)


