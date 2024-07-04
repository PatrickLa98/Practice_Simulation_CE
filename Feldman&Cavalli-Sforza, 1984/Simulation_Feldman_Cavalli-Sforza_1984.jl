using Plots
using Distributions

cd("/Users/patricklauer/Documents/GitHub/Practice_Simulation_CE/Feldman&Cavalli-Sforza, 1984/")

include("Equations_Feldman_Cavalli-Sforza_1984.jl")


    ## Iterating through generations and recording cultural variant frequencies and gene-culture disequilibrium

    μ = 0.01
    v = 0.01
    
    β1 = 0.5    ## A more likely to transmit 2
    β2 = 0.9        
    β3 = 0.9
    β4 = 0.9
    ## gamma values irrelevant for vertical transmission in haploids with very low mutation rate
    γ1 = 0.9 
    γ2 = 0.9
    γ3 = 0.9
    γ4 = 0.9
    x1 = 0.25
    x2 = 0.25
    x3 = 0.25
    x4 = 0.25

    n_gen = 100
    ## initial frequency of cultural variants
    y_num = [x1 + x3]
    ## initial gene-culture disequilibrium
    D_num = [x1*x4 - x2*x3]

    for i in 1:n_gen

        ## set geno-pheontype frequencies to next generation
        x1, x2, x3, x4 = Equation1(μ, v, β1, β2, β3, β4, γ1, γ2, γ3, γ4, x1, x2, x3, x4) 
        
        ## record cultural variant frequencies 
        push!(y_num, x1 + x3)
        ## record gene-culture disequilibrium
        push!(D_num, x1*x4 - x2*x3)
    end


  ## compare the equilibrium frequencies/disequilibrium found by iteration 
  ## with the equilibrium frequencies in postulated mathematical models

  ## set geno-pheotype frequencies back to initial state as input for mathematical models
  x1 = 0.25
  x2 = 0.25
  x3 = 0.25
  x4 = 0.25
    ## frequecy of cultural variant 1
    Equation6a(μ, v, β1, β2, β3, β4, γ1, γ2, γ3, γ4, x1, x2, x3, x4)
    ## only works with neglectable mutation
    Equation7(μ, v, β1, β2, β3, β4, γ1, γ2, γ3, γ4, x1, x2, x3, x4) 
    ## doesnt work
    EquationA7(μ, v, β1, β2, β3, β4, γ1, γ2, γ3, γ4, x1, x2, x3, x4)
    ## works best 
    
    ## gene-culture disequilibrium
    Equation6c(μ, v, β1, β2, β3, β4, γ1, γ2, γ3, γ4, x1, x2, x3, x4)
    ## doesn't match
    EquationA8(μ, v, β1, β2, β3, β4, γ1, γ2, γ3, γ4, x1, x2, x3, x4)
    ## matches

## ABM of model 1

# number of individuals
n_pop = 10000
## Genotype-Culture Combinations
genotype_culture = [1, 2, 3, 4]
## relative frequencies of genotype culture combinations x1, x2, x3, x4
freq = [0.25, 0.25, 0.25, 0.25]
## population according to frequency of genotype culture combinations
pop = genotype_culture[rand(Categorical(freq), n_pop)]
## transmission coefficients
## A1 to A1
β1 = 0.5
## A2 to A2
β2 = 0.9 
## a1 to a1
β3 = 0.9
## a2 to a2
β4 = 0.9

## mutation rates

μ = 0.001
v = 0.001

n_generations = 100
## initial frequency of cultural variant 1
x1 = sum(pop .== 1)/n_pop
x3 = sum(pop .== 3)/n_pop
y = [x1 + x3]

for g in 1:n_generations

    for i in 1:length(pop)

     ## replace parent genotypes according to mutation rates μ or values
        ## genotype A to a
        if pop[i] in [1, 2]

            choice = [pop[i], pop[i] + 2] # either remain same genotype or convert genotype (culture remains the same)
            prob = [1 - μ, μ]
            pop[i] = choice[rand(Categorical(prob))]    #convert offspring to genotype a with probability μ
        ## genotype a to A
        elseif pop[i] in [3, 4]

            choice = [pop[i], pop[i] + 2] # either remain same genotype or convert genotype (culture remains the same)
            prob = [1 - v, v]
            pop[i] = choice[rand(Categorical(prob))]    # convert offspring to genotype a with probability μ
        end



     ## replace or keep parents culture according to transmission coefficients
        if pop[i] == 1
            ## replace parent individual with offspring individual that has same culture with prob β1 
            choice = [1, 2]
            prob = [β1, (1 - β1)]

            pop[i] = choice[rand(Categorical(prob))]
        elseif pop[i] == 2
            ## replace parent individual with offspring individual that has same culture with prob β2 
            choice = [2, 1]
            prob = [β2, (1 - β2)]

            pop[i] = choice[rand(Categorical(prob))]
        elseif pop[i] == 3
            ## replace parent individual with offspring individual that has same culture with prob β3 
            choice = [3, 4]
            prob = [β3, (1 - β3)]

            pop[i] = choice[rand(Categorical(prob))]
        elseif pop[i] == 4
            ## replace parent individual with offspring individual that has same culture with prob β3 
            choice = [4, 3]
            prob = [β4, (1 - β4)]

            pop[i] = choice[rand(Categorical(prob))]
        end
    end
    ## keep track of frequency of cultural variants over time
    x1 = sum(pop .== 1)/n_pop
    x3 = sum(pop .== 3)/n_pop
    push!(y, x1 + x3)

end


 plot(
    collect(0:n_gen), y,
    line = (1, :red),
    marker = (:circle, 4, :red),
    xlabel = "time in generations",
    ylabel = "frequency of cultural variant1",
    label = "ABM"  
)
plot!(
    collect(0:n_gen), y_num,
    line = (1, :black),
    marker = (:circle, 4, :black),
    label = "Numeric"  
)
