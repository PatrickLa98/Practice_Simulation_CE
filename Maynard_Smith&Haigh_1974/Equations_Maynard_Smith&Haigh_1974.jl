## Goal: Model the hitch-hiking effect of a favourable gene

## parameter settings

    ## frequency of favourable allele B = p
        p = 0.5
    ## proportion of A in chromosomes containing B = Q
        Q = 0.8
    ## proportions of A in chromosomes containing b = R
        R = 0.4
    ## selection coefficient = s
        s = 0.2
    ## recombination fraction between the two loci = c
        c = 0.1
    
## --> frequency of genotype = freqAB 
        freqAB= p*Q

## Equations

## frequency of AB in generation n + 1
Equation0 = function(p, Q, R, s, c)
    
    next_freqAB = (((1 + s)^2)*(p^2)*(Q^2) + 
                  ((1 + s)^2)*(p^2)*Q*(1 - Q) +
                  (1 + s)*p*(1 - p)*Q*R + 
                  (1 + s)*p*(1 - p)*Q*(1 - R)*(1 - c) +
                  (1 + s)*p*(1 - p)*(1 - Q)*R*c)/
                  (1 + p*s)^2

                  return next_freqAB
end

Equation0(p, Q, R, s, c)

Equation1 = function(p, Q, R, s, c)

    next_freqAB = (p*(1 + s)*(Q*(1 + p*s) + c*(1 - p)*(R - Q)))/(1 + p*s)^2

    return next_freqAB
end

Equation1(p, Q, R, s, c)

## frequency of B in generation n + 1
Equation2 = function(p, Q, R, s, c)

    next_p = (p*(1 + s)*(1 + p*s))/(1 + p*s)^2 

    return next_p
end

Equation2(p, Q, R, s, c)

##  proportion of A in chromosomes containing B in generation n + 1
Equation3 = function(p, Q, R, s, c)

    next_Q = ((1 + p*s)*Q + c*(1 - p)*(R - Q))/(1 + p*s)

    return next_Q
end

Equation3(p, Q, R, s, c)


## Sidenote: Equation 2 times Equation 3 should be equal to output in Equation 1 since  it represents p*Q 
Equation2(p, Q, R, s, c) * Equation3(p, Q, R, s, c)

## proportion of A in chromosomes containing b in generation n + 1
Equation4 = function(p, Q, R, s, c)

    next_R = ((1 + p*s)*R + c*(1 + s)*p*(Q - R))/(1 + p*s)

    return next_R
end

Equation4(p, Q, R, s, c)


