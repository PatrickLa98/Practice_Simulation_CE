## Goal: Model the hitch-hiking effect of a favourable gene
#= 
for now only do the haploid Model
allele b subsitituted by favourable B
neighbouring locus neutral alleles a and A 
closer neighbouring loci are more likely to be "hitchhiked" denoted as distance/recombination fraction c (smaller values mean stronger hitchhiking effects)
B individuals descended from single mutant aB in initial generation
    -> thus the fraction of individuals with AB in the initial generation is Q0 = 0
=#

## parameter settings

    ## initial frequency of favourable allele B = p
        p0 = 0.01
    ## initial proportion of A in chromosomes containing B = Q
        Q0 = 0
    ## initial proportions of A in chromosomes containing b = R
        R0 = 0.4
    ## selection coefficient = s
        s = 0.2
    ## recombination fraction between the two loci = c
        c = 0.001
    ## number of generations = n
        n = 10
    
## --> frequency of genotype = freqAB 
        freqAB= p0*Q0

## Equations

## Equation 0: frequency of AB in generation n + 1
Equation0 = function(p, Q, R, s, c)
    
    next_freqAB = (((1 + s)^2)*(p^2)*(Q^2) + 
                  ((1 + s)^2)*(p^2)*Q*(1 - Q) +
                  (1 + s)*p*(1 - p)*Q*R + 
                  (1 + s)*p*(1 - p)*Q*(1 - R)*(1 - c) +
                  (1 + s)*p*(1 - p)*(1 - Q)*R*c)/
                  (1 + p*s)^2

                  return next_freqAB
end

Equation0(p0, Q0, R0, s, c)

## Equation 1: frequency of AB in generation n + 1 (simplified form of Equation 0)
Equation1 = function(p, Q, R, s, c)

    next_freqAB = (p*(1 + s)*(Q*(1 + p*s) + c*(1 - p)*(R - Q)))/(1 + p*s)^2

    return next_freqAB
end

Equation1(p0, Q0, R0, s, c)

## Equation 2: frequency of B in generation n + 1
Equation2 = function(p, Q, R, s, c)

    next_p = (p*(1 + s)*(1 + p*s))/(1 + p*s)^2 

    return next_p
end

Equation2(p0, Q0, R0, s, c)

## Equation 3: proportion of A in chromosomes containing B in generation n + 1
Equation3 = function(p, Q, R, s, c)

    next_Q = ((1 + p*s)*Q + c*(1 - p)*(R - Q))/(1 + p*s)

    return next_Q
end

Equation3(p0, Q0, R0, s, c)


## Sidenote: Equation 2 times Equation 3 should be equal to output in Equation 1 since  it represents p*Q 
#Equation2(p0, Q0, R0, s, c) * Equation3(p0, Q0, R0, s, c)

## Equation 4: proportion of A in chromosomes containing b in generation n + 1
Equation4 = function(p, Q, R, s, c)

    next_R = ((1 + p*s)*R + c*(1 + s)*p*(Q - R))/(1 + p*s)

    return next_R
end

Equation4(p0, Q0, R0, s, c)

## Equation 6: frequency of B in generation n (Equation 2 solved explicitly)
Equation6 = function(p0, s, c, n)

    p_n = p0*((1 + s)^n)/(1 - p0 + p0*((1 + s)^n))

    return p_n
end

Equation6(p0, s, c, n)

## Equation 7: proportion of A in chromosomes containing B in generation n + 1 (recurrence relation)
Equation7 = function(p0, Q, R0, s, c, n)

    next_Q = Q + c*R0*(1 - p0)*((1 - c)^n)/(1 - p0 + p0*(1 + s)^n + 1)

    return next_Q
end

Equation7(p0, 0.5, R0, s, c, n)


## Equation 8: Exact calculation of proporional frequency reduction of A after fixation of a linked favourable allele B
Equation8 = function(c, p0, s, R0) 
    Eq8_sum = 0
    for n in 0:10^5  ## 10^5 is arbitrary, used instead of infinity
        Eq8_sum += ((1 - c)^n)/(1 - p0 + p0*(1 + s)^(n+1))
    end
    Q_inf = c*R0*(1 - p0)*Eq8_sum
    ## divide Qinf by R0 to be consistent with Equation 14
    hitchhiking_effect = Q_inf / R0  # not sure but I think this is the wanted measure of the hitchhiking effect ??

    return hitchhiking_effect
end

Equation8(c, p0, s, R0)

### Equation 10: differentiation of allele B frequency with respect to time (in generations n)
Equation10 = function(s, p)

    p_dot = s*p*(1 - p)/(1 + s*p)

    return p_dot
end

Equation10(s, p0)

### Equation 11: differentiation of proportion of A in chromosomes containing B with respect to time (in generations n)
Equation11 = function(c, p, R, Q, s)

    Q_dot = c*(1 - p)*(R - Q)/(1 + s*p)

    return Q_dot
end

Equation11(c, p0, R0, Q0, s)

### Equation 12: differentiation of proportion of A in chromosomes containing b with respect to time (in generations n)
Equation12 = function(c, p, R, Q, s)

    R_dot = -c*p*(1 + s)*(R - Q)/(1 + s*p)

    return R_dot
end

Equation12(c, p0, R0, Q0, s)

## Equation 14: Approximation for the proporional frequency reduction of A after fixation of a linked favourable allele B - Note "exact" calculation in Equation 8 
Equation14 = function(c, s , p0)
    
    hitchhiking_effect = (c/s)*log(1/p0)
    
    return hitchhiking_effect
end

Equation14(c, s, p0)

