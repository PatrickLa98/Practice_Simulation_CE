using Plots

#=
Goal: Create an agent-based version of the Wright-Fisher model
       Many alleles version , with mutation

    Track the changes of Allele 1 over
    generations in a population with k 
    different Alleles and a constant mutation rate μ 
=#

## Function that creates Matrix containing individual alleles over generations (Many Alleles version)
function wf_sampling_ma(N, A0, k, μ, Generations) 

 ## Individuals with allele 1 
 n_1 = floor(Int, N * A0) 
 Ind_1 = fill(1,n_1)

 ## Fill up the remaining population randomly with k variants
 n_!1 = N - n_1
 variants = collect(2:k)
 Ind_!1 = rand(variants, n_!1)

 ## Combine population of individuals with variant 1 and individuals without variant 1
 population = [Ind_1; Ind_!1]

    results = zeros(Int64, Generations + 1, N)
    results[1, :] = population # first row in matrix is the starting composition of population
  
    for i in 1:Generations

        population = rand(population, N) # randomly draw alleles from previous population with replacement
        
        ### the newly assigned variant can mutate with mutationrate \mu and become the new variant k_max + 1
        for j in 1:length(population)
                     
            if rand() < μ
            population[j] = k + 1
            end

            k = maximum(population) # need to update variant number after potential mutation event

        end

        results[i+1,:] = population

    end
    return results
end 

## Function that transforms individual alleles to allele proportions (Many Alleles Version)
function allele_proportions_ma(N, A0, k, μ, Generations)

    wf_matrix = wf_sampling_ma(N, A0, k, μ, Generations) 
    Allele_Proportion = zeros(Generations +1)
    for i in 1:Generations+1

        Allele_Proportion[i] = count(wf_matrix[i,:] .== 1) / N

    end

    return Allele_Proportion
end

## Starting Conditions   
N = 200 # population size
A0 = 0.2 # starting frequency of Allele 1
k = 10 # number of variants
μ = 0 # mutation rate
Generations = 2000 # number of generations

wf_sampling_ma(N, A0, k, μ, Generations)
allele_proportions_ma(N, A0, k, μ, Generations)

## Plotting allele frequencies over time with r replicates

r = 30

wf_plot = plot(allele_proportions_ma(N, A0, k, μ, Generations), 
legend = false,
xlabel = "# of Generations",
ylabel = "Frequency Allele A",
color = "black"
)

for r in 1:r 
wf_plot = plot!(wf_plot,
 allele_proportions_ma(N, A0, k, μ, Generations),
 color = "black")
end

wf_plot