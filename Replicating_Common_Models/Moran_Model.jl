using Plots

#=
Goal: Create a simulation of the Moran model
Constant finite population with size N
allele with k variants
1 birth and 1 death per timestep 
tracking the frequency of allele 1 that has the starting frequency A0
=#



function moran_model(N, k, A0, μ, generations)

## setting starting population according to input allele frequency
allele_1 = fill(1, floor(Int64, N*A0))
allele_!1= rand(collect(2:k), N - length(allele_1))
population = [allele_1; allele_!1]


### iterating through generations

mmm = zeros(generations + 1, N) # mmm = moran model matrix
mmm[1,:] = population # first row is the starting allele frequency

for i in 1:generations

birth = rand(population, 1) # sampling birthing individual
## individual with new variant is instead born proportional to mutation rate μ
if μ > rand()
    birth = maximum(population) + 1
end
death = rand(population, 1) # sampling dieing individual

population = [population; birth] # update population with 1 birth event
deleteat!(population, findfirst(population .== death)) # update population with 1 death event 

mmm[i+1,:] = population

end

return mmm

end

function moran_model_proportions(N, k, A0, μ, generations)
    
mmm = moran_model(N, k, A0, generations)  #accessing previous defined function to get individual based matrix of moran model

prop_mm = zeros(generations + 1)  #set up vector to store proportions

## iterate through generations to compress individual based alleles to population level proportions
for i in 1: generations + 1
prop_mm[i] = count(mmm[i,:] .== 1) / N
end

return prop_mm

end

## inputs
N = 100
k = 5
A0 = 0.6
μ = 0.01
generations = 10000

 #mmm = moran_model(N,k,A0, μ, generations)

## Plotting allele frequencies over time with r replicates
r = 5

mm_plot = plot(moran_model_proportions(N, k, A0, generations), 
legend = false,
xlabel = "# of Generations",
ylabel = "Frequency Allele A",
color = "black"
)

for r in 1:r 
mm_plot = plot!(mm_plot,
moran_model_proportions(N, k, A0, generations),
 color = "black")
end

mm_plot



