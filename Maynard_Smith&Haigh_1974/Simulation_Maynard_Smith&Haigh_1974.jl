## Goal: Model the hitch-hiking effect of a favourable gene
## numerical simulation and abm 
#= 
for now only do the haploid Model
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
    N = 10^4
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
recombination_fractions = [c for c in 0:0.001:0.01]

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


        



