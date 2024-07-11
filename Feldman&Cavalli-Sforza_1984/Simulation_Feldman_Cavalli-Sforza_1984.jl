## Goal: Model a dichotomous (2 variants)trait, genetic transmission diallellc haploid gene (A and a)
## Cultural transmission with 2 cultural states (indicated by subscripts 1 and 2):
    ## 1. vertical
    ## 2. oblique
    ## 3. vertical with selection


    ## pheno-genotype combinations A1, A2, a1, a2 with denoted frequencies x1, x2, x3, x4
    ## combinations lead to 8 transmission coefficients from parent to child (see Table 1 in paper) denoted as:
        ## β (with respective subscript) for the transmission probabilities of same genotype parent and child,
        ## γ (with respective subscript) for the transmission probabilities of different genotype parent and child
    ## mutations are denoted as μ for A to a and v for a to A

using Plots
using Distributions

## source the equations from paper to this script
cd("/Users/patricklauer/Documents/GitHub/Practice_Simulation_CE/Feldman&Cavalli-Sforza_1984/")
include("Equations_Feldman_Cavalli-Sforza_1984.jl")

## parameter settings

## mutation rates
    ## A to a
    μ = 0.01
    ## a to A
    v = 0.001
## transmission coefficients
    ## A1 to A1
    β1 = 0.95
    ## A2 to A2
    β2 = 0.9 
    ## a1 to a1
    β3 = 0.6
    ## a2 to a2
    β4 = 0.95
    ## A1 to a1
    γ1 = 0.7 
    ## A2 to a2
    γ2 = 0.9
    ## a1 to A1
    γ3 = 0.9
    ## a2 to A2
    γ4 = 0.9
## starting frequency of geno-culturetype
    ## A1
    x1 = 0.25
    ## A2
    x2 = 0.25
    ## a1
    x3 = 0.25
    ## a2
    x4 = 0.25
## number of generations   
    n_gen = 500

##################################################################
###### MODEL 1: VERTICAL TRANSMISSION WITH NO SELECTION ##########
##################################################################

## Numerical Simulation
## Iterating through generations and recording cultural variant frequencies and gene-culture disequilibrium
    
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


  #########################################################################
  ## NOT PART OF NUMERICAL SIMULATION, JUST COMPARE GIVEN EQUATIONS WITH RESULTS OF ITERATING THROUGH GENERATIONS NUMERICALLY
    ## compare the equilibrium frequencies/disequilibrium found by iteration 
    ## with the equilibrium frequencies in postulated mathematical models
    ## set geno-pheotype frequencies back to initial state as input for mathematical models (need to be reset as frequencies update through iterations)
     x1 = 0.25
     x2 = 0.25
     x3 = 0.25
     x4 = 0.25
        ## frequecy of cultural variant 1,
            ## if no mutation is included
            Equation6a(μ, v, β1, β2, β3, β4, γ1, γ2, γ3, γ4, x1, x2, x3, x4)
            ## dont fully understand under which condition this works
            Equation7(μ, v, β1, β2, β3, β4, γ1, γ2, γ3, γ4, x1, x2, x3, x4) 
            ## works best! if neglecting quadratic terms in μ and v (also not sure what that means)
            EquationA7(μ, v, β1, β2, β3, β4, γ1, γ2, γ3, γ4, x1, x2, x3, x4)
        
        ## gene-culture disequilibrium
            ## if no mutation is included
            Equation6c(μ, v, β1, β2, β3, β4, γ1, γ2, γ3, γ4, x1, x2, x3, x4)
            ## works best! if neglecting quadratic terms in μ and v (also not sure what that means)
            EquationA8(μ, v, β1, β2, β3, β4, γ1, γ2, γ3, γ4, x1, x2, x3, x4)
 #############################################################

## ABM of model 1

## number of individuals
n_pop = 1000
## Genotype-Culture Combinations
genotype_culture = [1, 2, 3, 4]
## relative frequencies of genotype culture combinations x1, x2, x3, x4
freq = [0.25, 0.25, 0.25, 0.25]
## population according to frequency of genotype culture combinations
pop = genotype_culture[rand(Categorical(freq), n_pop)]
## initial frequency of cultural variant 1
x1 = sum(pop .== 1)/n_pop
x3 = sum(pop .== 3)/n_pop
y_abm = [x1 + x3]
## intital gene-culture disequilibrium
x2 = sum(pop .== 2)/n_pop
x4 = sum(pop .== 4)/n_pop
D_abm = [x1*x4 - x2*x3]

for g in 1:n_gen

    for i in 1:length(pop)

     ## replace parent genotypes according to mutation rates μ or v
        ## genotype A to a
        if pop[i] in [1, 2]

            choice = [pop[i], pop[i] + 2] # either remain same genotype or convert genotype (culture remains the same)
            prob = [1 - μ, μ]
            pop[i] = choice[rand(Categorical(prob))]    # convert offspring to genotype a with probability μ
        ## genotype a to A
        elseif pop[i] in [3, 4]

            choice = [pop[i], pop[i] - 2] # either remain same genotype or convert genotype (culture remains the same)
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
    ## keep track of frequency of cultural variants and disequilibrium  over time
    x1 = sum(pop .== 1)/n_pop
    x2 = sum(pop .== 2)/n_pop
    x3 = sum(pop .== 3)/n_pop
    x4 = sum(pop .== 4)/n_pop

    push!(y_abm, x1 + x3)

    push!(D_abm, x1*x4 - x2*x3)
end

## PLOT

## Frequency of cultural variant 1 over time

plot(
    collect(0:n_gen), y_abm,
    line = (1, :red),
    marker = (:circle, 4, :red),
    xlabel = "time in generations",
    ylabel = "frequency of cultural variant1",
    label = "ABM",
    title = "Model 1 - Cultural variant frequencies"   
)
plot!(
    collect(0:n_gen), y_num,
    line = (1, :blue),
    marker = (:circle, 4, :blue),
    label = "Numeric")

## Gene-culture disequilibrium over time
plot(
    collect(0:n_gen), D_abm,
    line = (1, :red),
    marker = (:circle, 4, :red),
    xlabel = "time in generations",
    ylabel = "gene-culture disequilibrium",
    label = "ABM",
    title = "Model 1 - Gene-culture disequilibrium"  
)
plot!(
    collect(0:n_gen), D_num,
    line = (1, :blue),
    marker = (:circle, 4, :blue),
    label = "Numeric")

##################################################################
###### MODEL 2: OBLIQUE TRANSMISSION WITH NO SELECTION ###########
##################################################################

## Numerical Simulation
## Iterating through generations and recording cultural variant frequencies and gene-culture disequilibrium
    ## set geno-pheotype frequencies back to initial state  (need to be reset as frequencies updated by prior iterations)
    x1 = 0.25
    x2 = 0.25
    x3 = 0.25
    x4 = 0.25
    ## initial frequency of cultural variants
    y_num2 = [x1 + x3]
    ## initial gene-culture disequilibrium
    D_num2 = [x1*x4 - x2*x3]

    for i in 1:n_gen

        ## set geno-pheontype frequencies to next generation
        x1, x2, x3, x4 = Equation8(μ, v, β1, β2, β3, β4, γ1, γ2, γ3, γ4, x1, x2, x3, x4) 
        
        ## record cultural variant frequencies 
        push!(y_num2, x1 + x3)
        ## record gene-culture disequilibrium
        push!(D_num2, x1*x4 - x2*x3)
    end

   ###################################################################################
   ## NOT PART OF NUMERICAL SIMULATION, JUST COMPARE GIVEN EQUATIONS WITH RESULTS OF ITERATING THROUGH GENERATIONS NUMERICALLY
    ## compare the equilibrium frequencies/disequilibrium found by iteration 
    ## with the equilibrium frequencies in postulated mathematical models
    ## set geno-pheotype frequencies back to initial state as input for mathematical models (need to be reset as frequencies update through iterations)
    x1 = 0.25
    x2 = 0.25
    x3 = 0.25
    x4 = 0.25
    ## frequecy of cultural variant 1 if βi = γi, equil_D = 0
    Equation11(μ, v, β1, β2, β3, β4, γ1, γ2, γ3, γ4, x1, x2, x3, x4)

    ## gene-culture disequilibrium if no genetic differences in teaching abilities within each cultural state
    Equation12(μ, v, β1, β2, β3, β4, γ1, γ2, γ3, γ4, x1, x2, x3, x4)
   #########################################################################

## ABM of model 2

## Create a Dictionary for transmission coefficients ("parent-genophenotype offspring-genophenotype")
Transmission = Dict{String, Float64}("1 1" => β1, 
                      "2 2" => β2,
                      "3 3" => β3, 
                      "4 4" => β4, 
                      "1 3" => γ1, 
                      "2 4" => γ2, 
                      "3 1" => γ3, 
                      "4 2" => γ4,
                      "1 2" => 1 - β1, 
                      "2 1" => 1 - β2,
                      "3 4" => 1 - β3, 
                      "4 3" => 1 - β4, 
                      "1 4" => 1 - γ1, 
                      "2 3" => 1 - γ2, 
                      "3 2" => 1 - γ3, 
                      "4 1" => 1 - γ4)
## number of individuals
n_pop = 1000
## Genotype-Culture Combinations
genotype_culture = [1, 2, 3, 4]
## relative frequencies of genotype culture combinations x1, x2, x3, x4
freq = [0.25, 0.25, 0.25, 0.25]
## population according to frequency of genotype culture combinations
pop = genotype_culture[rand(Categorical(freq), n_pop)]
## initial frequency of cultural variant 1
x1 = sum(pop .== 1)/n_pop
x3 = sum(pop .== 3)/n_pop
y_abm2 = [x1 + x3]
## intital gene-culture disequilibrium
x2 = sum(pop .== 2)/n_pop
x4 = sum(pop .== 4)/n_pop
D_abm2 = [x1*x4 - x2*x3]

## separate genotype and culture transmission, 2 step process:
## 1. genotype frequency of gen i similar to gen i - 1, differences due to mutation rates

for g in 1:n_gen

    pop_new_gen = zeros(Int64, length(pop))

    for i in 1:length(pop)
        ## replace "parent" generation genotypes according to mutation rates μ or v
        ## genotype A to a
        if pop[i] in [1, 2]

            choice = ["A","a"] # either remain same genotype or convert genotype (culture remains the same)
            prob = [1 - μ, μ]
            genotype = choice[rand(Categorical(prob))]    # convert offspring to genotype a with probability μ
        ## genotype a to A
        elseif pop[i] in [3, 4]

            choice = ["a", "A"] # either remain same genotype or convert genotype (culture remains the same)
            prob = [1 - v, v]
            genotype = choice[rand(Categorical(prob))]    # convert offspring to genotype a with probability μ
        end

     ## add culture to the produced genotype, by randomly selecting an individual from "parent" generation
        ## sample from parent gen
        model = rand(pop)
            
        if genotype == "A"
          ## if id has genotype A it can adapt genophenotype 1 or 2 (copies from model according to transmission coefficents)
          choices =  [string(model, " ", 1), string(model, " ", 2)]
          prob = [Transmission[choices[1]], Transmission[choices[2]]]
          
          pop_new_gen[i] = [1, 2][rand(Categorical(prob))]   

        elseif genotype == "a"
            ## if id has genotype a it can adapt genophenotype 3 or 4 (copies from model according to transmission coefficents)
            choices =  [string(model, " ", 3), string(model, " ", 4)]
            prob = [Transmission[choices[1]], Transmission[choices[2]]]
            
            pop_new_gen[i] = [3, 4][rand(Categorical(prob))]
        end    

    end

    ## update generation after all ids are done learning from parent generation
    pop = pop_new_gen
    ## keep track of frequency of cultural variants and disequilibrium  over time
    x1 = sum(pop .== 1)/n_pop
    x2 = sum(pop .== 2)/n_pop
    x3 = sum(pop .== 3)/n_pop
    x4 = sum(pop .== 4)/n_pop

    push!(y_abm2, x1 + x3)

    push!(D_abm2, x1*x4 - x2*x3)
end

## PLOT

## Frequency of cultural variant 1 over time
plot(
    collect(0:n_gen), y_abm2,
    line = (1, :red),
    marker = (:circle, 4, :red),
    xlabel = "time in generations",
    ylabel = "frequency of cultural variant1",
    label = "ABM",
    title = "Model 2 - Cultural variant frequencies" 
)
plot!(
    collect(0:n_gen), y_num2,
    line = (1, :blue),
    marker = (:circle, 4, :blue),
    label = "Numeric")

## Gene-culture disequilibrium over time
plot(
    collect(0:n_gen), D_abm2,
    line = (1, :red),
    marker = (:circle, 4, :red),
    xlabel = "time in generations",
    ylabel = "gene-culture disequilibrium",
    label = "ABM",  
    title = "Model 2 - Gene-culture disequilibrium"
)
plot!(
    collect(0:n_gen), D_num2,
    line = (1, :blue),
    marker = (:circle, 4, :blue),
    label = "Numeric")


##################################################################
###### MODEL 3: VERTICAL TRANSMISSION WITH SELECTION #############
##################################################################

## no mutation, γ values irrelevant

## selection coefficient
s = 0.1

## reset gene-culturetypes
x1 = 0.25
x2 = 0.25
x3 = 0.25
x4 = 0.25
## Numerical Simulation
## Iterating through generations and recording cultural variant frequencies and gene-culture disequilibrium
    
    ## initial frequency of cultural variants
    y_num3 = [x1 + x3]
    ## initial gene-culture disequilibrium
    D_num3 = [x1*x4 - x2*x3]

    for i in 1:n_gen

        ## set geno-pheontype frequencies to next generation
        x1, x2, x3, x4 = Equation14(β1, β2, β3, β4, x1, x2, x3, x4, s) 
        
        ## record cultural variant frequencies 
        push!(y_num3, x1 + x3)
        ## record gene-culture disequilibrium
        push!(D_num3, x1*x4 - x2*x3)
    end

## ABM of model 3
#=
1. randomly sampled parent repoduces asexually, offspring always adopts its genotype (no mutation). 
2. The cultural state of the parent is transmitted according to transmition coefficents.
3. Offspring has chance to die if there is selection against the adopted cultural variant (survival probability 1 - s)
If offspring dies step 1 - 3 repeat
4. if offspring survives it enters the next_population
5. iteration i moves to i + 1
6. as soon as the offspring population (next_population) reaches the number of individuals of parent population 
   the offspring population replaces the parent population
=#


## number of individuals
n_pop = 1000
## Genotype-Culture Combinations
genotype_culture = [1, 2, 3, 4]
## relative frequencies of genotype culture combinations x1, x2, x3, x4
freq = [0.25, 0.25, 0.25, 0.25]
## population according to frequency of genotype culture combinations
pop = genotype_culture[rand(Categorical(freq), n_pop)]
## initial frequency of cultural variant 1
x1 = sum(pop .== 1)/n_pop
x3 = sum(pop .== 3)/n_pop
y_abm3 = [x1 + x3]
## intital gene-culture disequilibrium
x2 = sum(pop .== 2)/n_pop
x4 = sum(pop .== 4)/n_pop
D_abm3 = [x1*x4 - x2*x3]
# vector containing offspring individuals
next_pop = zeros(Int64, length(pop))

for g in 1:n_gen

    i = 1

    ## while loop because parent id can be "reproduce" multiply times if offspring dies 
    while i <= length(pop)

        parent = rand(pop)
        
        
        ## replace or keep parents culture according to transmission coefficients
        if parent == 1
            ## replace parent individual with offspring individual that has same culture with prob β1 
            choice = [1, 2]
            prob = [β1, (1 - β1)]

            offspring = choice[rand(Categorical(prob))]
        elseif parent == 2
            ## replace parent individual with offspring individual that has same culture with prob β2 
            choice = [2, 1]
            prob = [β2, (1 - β2)]

            offspring = choice[rand(Categorical(prob))]
        elseif parent == 3
            ## replace parent individual with offspring individual that has same culture with prob β3 
            choice = [3, 4]
            prob = [β3, (1 - β3)]

            offspring = choice[rand(Categorical(prob))]
        elseif parent == 4
            ## replace parent individual with offspring individual that has same culture with prob β3 
            choice = [4, 3]
            prob = [β4, (1 - β4)]

            offspring = choice[rand(Categorical(prob))]
        end

        ## assume that survival probability of ids with cultural variant1 is 1 and of ids with cultural variant2 is 1 - s
        ## if an individual dies replace it with a new offspring (from randomly sampled parent)

        if offspring in [2, 4]
            
          death = rand(Categorical(s, 1 - s)) # 1 = dies 2 = survives
          
          ## if offspring doesnt die replace parent with it and move to next parent individual
            if death == 2  
            
                next_pop[i] = offspring
                i += 1
            else
                continue
        
            end
          ## if offspring has no maladaptive cultural variant it survives, replaces parent in population and loop goes to next parent individual  
        elseif offspring in [1, 3]

            next_pop[i] = offspring
            i += 1
        end

    end
    ## replace parent with offspring generation
    pop = next_pop

    ## keep track of frequency of cultural variants and disequilibrium  over time
    x1 = sum(pop .== 1)/n_pop
    x2 = sum(pop .== 2)/n_pop
    x3 = sum(pop .== 3)/n_pop
    x4 = sum(pop .== 4)/n_pop

    push!(y_abm3, x1 + x3)

    push!(D_abm3, x1*x4 - x2*x3)
end


## PLOT

## Frequency of cultural variant 1 over time
plot(
    collect(0:n_gen), y_abm3,
    line = (1, :red),
    marker = (:circle, 4, :red),
    xlabel = "time in generations",
    ylabel = "frequency of cultural variant1",
    label = "ABM",
    title = "Model 3 - Cultural variant frequencies" 
)
plot!(
    collect(0:n_gen), y_num3,
    line = (1, :blue),
    marker = (:circle, 4, :blue),
    label = "Numeric")

## Gene-culture disequilibrium over time
plot(
    collect(0:n_gen), D_abm3,
    line = (1, :red),
    marker = (:circle, 4, :red),
    xlabel = "time in generations",
    ylabel = "gene-culture disequilibrium",
    label = "ABM",  
    title = "Model 3 - Gene-culture disequilibrium"
)
plot!(
    collect(0:n_gen), D_num3,
    line = (1, :blue),
    marker = (:circle, 4, :blue),
    label = "Numeric")


    ## to avoid exporting plots in weird .svg format use the following command (store plot as p first)
    # savefig(p, "plot.png")




    










