## Numerical Simulation based on mathematical models in Cultural and biological evolutionary processes: Gene-culture disquilibrium, Feldman and Cavalli-Sforza, 1984

## Goal model a dichotomous (2 variants)trait, genetic transmission diallellc haploid gene (A and a)
## Cultural transmission with 2 cultural states (indicated by subscripts 1 and 2):
    ## 1. vertical
    ## 2. oblique
    ## 3. vertical with selection


## 1. Vertical transmission with no selection

    ## pheno-genotype combinations A1, A2, a1, a2 with denoted frequencies x1, x2, x3, x4
    ## combinations lead to 8 transmission coefficients from parent to child (see Table 1 in paper) denoted as:
    ## β (with respective subscript) for the transmission probabilities of same genotype parent and child,
    ## γ (with respective subscript) for the transmission probabilities of different genotype parent and child
    ## mutations are denoted as μ for A to a and v for a to A

    ## parameter settings
    μ = 0.001
    v = 0.001
    β1 = 0.9
    β2 = 0.9
    β3 = 0.9
    β4 = 0.9
    γ1 = 0.9
    γ2 = 0.9
    γ3 = 0.9
    γ4 = 0.9
    x1 = 0.25
    x2 = 0.25
    x3 = 0.25
    x4 = 0.25

    ## next generation frequencies
    Equation1 = function (μ, v, β1, β2, β3, β4, γ1, γ2, γ3, γ4, x1, x2, x3, x4)

        next_x1 = (1 - μ)*(β1*x1 + (1 - β2)*x2) + v*(γ3*x3 + (1 - γ4)*x4)
        next_x2 = (1 - μ)*((1 - β1)*x1 + β2*x2) + v*((1-γ3)*x3 + γ4*x4)
        next_x3 = μ*(γ1*x1 + (1 - γ2)*x2) + (1 - v)*(β3*x3 + (1 - β4)*x4)
        next_x4 = μ*((1 - γ1)*x1 + γ2*x2) + (1 - v)*((1 - β3)*x3 + β4*x4)
        
        return next_x1, next_x2, next_x3, next_x4
    end
    
    Equation1(
        0.001, 0.001,
        0.9, 0.9, 0.9, 0.9,
        0.9, 0.9, 0.9, 0.9,
        0.25, 0.25, 0.25, 0.25
    )

    ## next generation Allele A frequency
    Equation2a = function(x1, x2, μ, v)
        ## frequency of allele A is 
        p = x1 + x2
        ## frequency of next generation allele A is 
        next_p = (1 - μ)*p + v*(1-p)
    
        return next_p
    end

    Equation2a(0.25, 0.25, 0.001, 0.001)

    ## next generation culutural state 1 frequency
    Equation2b = function(μ, v, β1, β2, β3, β4, γ1, γ2, γ3, γ4, x1, x2, x3, x4)
        ## frequency of cultural state 1
        y = x1 + x3
        ## frequency of allele A is (used to calculte Coefficients)
        p = x1 + x2
        ## Coefficient (described in Appendix 1)
        L = p*((1 - μ)*(β1 + β2 - 1) + μ*(γ1 + γ2 - 1)) +
            (1 - p)*(v*(γ3 + γ4 - 1) + (1 - v)*(β3 + β4 -1))
        ## Coefficient (described in Appendix 2)
        M = (1 - μ) * (β1 + β2 - 1) + μ*(γ1 + γ2 - 1) +
            (1 - v) * (1 - β3 - β4) + v*(1 - γ3 - γ4)
        ## Gene-Culture Disequilibrium
        D = x1*x4 - x2*x3
        ## Coefficient (described in Appendix 3)
        N = p*((1 - μ)*(1 - β2) + μ*(1 - γ2)) + (1 - p)*(v*(1 - γ4) +
            (1 - v)*(1 - β4))
        ## next frequency of cultural state 1
        next_y = L*y + M*D + N

        return next_y
    end

    Equation2b(
        0.001, 0.001,
        0.9, 0.9, 0.9, 0.9,
        0.9, 0.9, 0.9, 0.9,
        0.25, 0.25, 0.25, 0.25
    )

    ## next generation gene culture Disequilibrium
    Equation2c = function(μ, v, β1, β2, β3, β4, γ1, γ2, γ3, γ4, x1, x2, x3, x4)
      ## frequency of allele A is (used to calculte Coefficients)
      p = x1 + x2
      ## Coefficient (described in Appendix 4)
      R = μ*(1 -μ)*(p^2)*(β1 + β2 - γ1 - γ2) +
          v*(1 - v)*((1 - p)^2)*(γ3 + γ4 - β3 - β4) +
          p*(1 - p)*((1 - μ)*(1 - v)*(β1 + β2 - β3 - β4) +
          μ*v*(γ3 + γ4 - γ1 - γ2))
      ## frequency of cultural state 1
      y = x1 + x3
      ## Coefficient (described in Appendix 5)
      S = μ*(1 -μ)*p*(β1 + β2 - γ1 - γ2) + v*(1 - v)*(1 - p)*(β3 + β4 - γ3 - γ4) +
          (1 - μ)*(1 - v)*(p*(β3 + β4 - 1) +
          (1 - p)*(β1 + β2 - 1)) + μ*v*(p*(1 - γ3 - γ4) +
          (1 - p)*(1 - γ1 - γ2))
      ## Gene-Culture Disequilibrium
      D = x1*x4 - x2*x3
      ## Coefficient (described in Appendix 6)
      T =  μ*(1 - μ)*(p^2)*(γ2 - β2) + v*(1 - v)*((1 - p)^2)*(β4 - γ4) +
      p*(1 - p) * ((1 - μ)*(1 - v)*(β4 - β2) + μ*v*(γ2 - γ4))
      ## Next Generation Gene-Culture Disequilibrium
      next_D = R*y + S*D + T

      return next_D
    end

    Equation2c(
        0.001, 0.001,
        0.9, 0.9, 0.9, 0.9,
        0.9, 0.9, 0.9, 0.9,
        0.25, 0.25, 0.25, 0.25
    )

    ## Convergence properties of y (frequency of cultural variant 1) and D (gene cultural Disequilibrium) over time 
    ## can be represented by matrix U, given p (frequency of Allele A) is in equilibrium state
    Equation3 = function(μ, v, β1, β2, β3, β4, γ1, γ2, γ3, γ4, x1, x2, x3, x4)
        ## Coefficient (described in Appendix 1)
        L = p*((1 - μ)*(β1 + β2 - 1) + μ*(γ1 + γ2 - 1)) +
        (1 - p)*(v*(γ3 + γ4 - 1) + (1 - v)*(β3 + β4 - 1))
        ## Coefficient (described in Appendix 2)
        M = (1 - μ) * (β1 + β2 - 1) + μ*(γ1 + γ2 - 1) +
        (1 - v) * (1 - β3 - β4) + v*(1 - γ3 - γ4)
        ## Coefficient (described in Appendix 4)
        R = μ*(1 - μ)*(p^2)*(β1 + β2 - γ1 - γ2) +
        v*(1 - v)*((1 - p)^2)*(γ3 + γ4 - β3 - β4) +
        p*(1 - p)*((1 - μ)*(1 - v)*(β1 + β2 - β3 - β4) +
        μ*v*(γ3 + γ4 - γ1 - γ2))   
        ## Coefficient (described in Appendix 5)
        S = μ*(1 - μ)*p*(β1 + β2 - γ1 - γ2) + v*(1 - v)*(1 - p)*(β3 + β4 - γ3 - γ4) +
        (1 - μ)*(1 - v)*(p*(β3 + β4 - 1) +
        (1 - p)*(β1 + β2 - 1)) + μ*v*(p*(1 - γ3 - γ4) +
        (1 - p)*(1 - γ1 - γ2))
        ## Convergence Matrix U
        U = [L M
            R S]

        return U
    end


   U =  Equation3(
        0.001, 0.001,
        0.9, 0.9, 0.9, 0.9,
        0.9, 0.9, 0.9, 0.9,
        0.25, 0.25, 0.25, 0.25
    )

    using LinearAlgebra
    eigen(U)

    ## The eigenvalues of U can also be represented by the roots of the quadratic Q(λ) = 0

    ## Equation 5
    λ12 = β1 + β2 - 1
    λ34 = β3 + β4 - 1
    γ12 = γ1 + γ2 - 1
    γ34 = γ3 + γ4 - 1

    ## Equation 4
    Q(λ) = (λ^2) - λ*(λ12*(1 - μ) + λ34*(1 - v)) +
            λ12*λ34*(1 - v)*(1 - μ) - γ12*γ34*μ*v

    λ_seq = collect(0:0.1:1)        
    Qs = [Q(λ) for λ in λ_seq]
    using Plots
    plot(λ_seq, Qs, legend = false)


    ## Equilibrium y if no mutation occurs
    Equation6a = function(μ, v, β1, β2, β3, β4, γ1, γ2, γ3, γ4, x1, x2, x3, x4)
        ## if no mutation and no selection we assume p0 is p_equilibrium
        p0 = x1 + x2
        ## Coefficient H (Equation 6b)
        H = (1 - β2)*(1 - β3) - (1 - β1)*(1 - β4)

        equil_y = ((1 - λ12)^-1)*((1 - λ34)^-1)*((1 - β4)*(1 - λ12) + p0*H)
        
        return equil_y
    end
    
    Equation6a(
        0.001, 0.001,
        0.9, 0.9, 0.9, 0.9,
        0.9, 0.9, 0.9, 0.9,
        0.25, 0.25, 0.25, 0.25
    )

    ## Equilibrium D if no mutation occurs
    Equation6c = function(μ, v, β1, β2, β3, β4, γ1, γ2, γ3, γ4, x1, x2, x3, x4)
        ## if no mutation and no selection we assume p0 is p_equilibrium
        p0 = x1 + x2
        ## Coefficient H (Equation 6b)
        H = (1 - β2)*(1 - β3) - (1 - β1)*(1 - β4)

        equil_D = p0*(1 - p0)*H/(1 - λ12)*(1 - λ34)

        return equil_D
    end

    Equation6c(
        0.001, 0.001,
        0.9, 0.9, 0.9, 0.9,
        0.9, 0.9, 0.9, 0.9,
        0.25, 0.25, 0.25, 0.25
    )

    



