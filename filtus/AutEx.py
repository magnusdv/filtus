from math import exp, log, expm1, log1p
 
#############################################
### Fra Kristinas hmm.py ###################
#############################################

### Endringer av MDV:
# 1) Obsvec, posvec og Bfreqvec er naa argumenter i funksjonen (for aa unngaa unodv lesing av fil)
# 2) Fjernet argument 'logs' (siden dette ikke er optional i resten av koden)
# 3) Forenklet output til en enkel vektor: posteriori-sannsynlighetene for IBD
# 4) Mer noyaktig og raskere logaritme-beregninger i transition() (se eget notat)
# 5) Endret koding av obsvektoren, slik at 0=homREF, 1=heteroz, 2=homALT

def logsum(v):
    """ Calculating the sum of a vector of logarithmic numbers."""
    a = max(v)
    return a + log(sum(exp(x-a) for x in v))
    
def transition(dist, a, f, logspace=0):
    """
    Compute transition probabilities for a HMM. 
    
    to compute transition probabilities between hidden-states,
    when moving from time t to t+1,
    the genetic distance (cM) between the two markers are required.
    
    Assuming known parameters a and f.
    lf logspace = 1, calculations are log-transformed.
    
    Key in dictionary: 0 = not-IBD, 1 = IBD.
    """
    if logspace == 0:
        qk = exp(-a*dist)

        T = { # 0 = not-IBD, 1 = IBD
            1: {1: (1-qk)*f + qk, 0: (1-qk)*(1-f)},
            0: {1: (1-qk)*f, 0: (1-qk)*(1-f) + qk}
            }
        
    else:
        if dist == 0:
            dist = 1e-06

        ff = 1-f
        ad = a*dist
        A = expm1(ad)
        AA = -expm1(-ad)
        T = { # 0 = not-IBD, 1 = IBD
            1: {1: log1p(A*f)-ad, 0: log(AA*ff)},
            0: {1: log(AA*f), 0: log1p(A*ff)-ad}}
    return T

def initiation(f, logspace=0):
    """
    Compute initiation probabilities for a HMM.

    Assuming that the probability for the first state being IBD or not-IBD
    is equivalent with the inbreeding-coefficient (f).

    lf logspace = 1, calculations are log-transfermed.
    
    Key in dictionary: 0 = not-IBD, 1 = IBD.
    """

    if logspace==0:
        I = {1: f, 0: 1-f} 
    else:
        I = {1: log(f), 0: log(1-f)}
    
    return I

def emission(q, error, logspace=0):
    """
    Compute emission probabilities for a HMM.
    
    This function takes as input three parameters:
    q = the minor allele frequency
    error = error-rate
    logspace = if 1, calculations are log-transfermed

    The output is an emission-matrix of probabilities between a given hidden-state
    and the observed event. 
    
    Key in dictionary: 0 = not-IBD, 1 = IBD.
    #Key in inner dictionary: 0 = heterozygous, 1 = hom ref allele, 2=hom alt allele.
    Key in inner dictionary: 1 = heterozygous, 0 = hom ref allele, 2=hom alt allele.
    """
    
    if logspace==0:
        p = 1-q
    
        E = {
            1: {0: (1-error)*p + error*p**2, 1: error*2*p*q, 
                    2: (1-error)*q + error*q**2 },
            0: {0: p**2, 1: 2*p*q, 2: q**2}
            }
          
    else:
        #q = min(0.999999, max(q, 1e-6))
        #error = min(0.999999, max(error, 1e-6))
        
        lp = log(1-q)
        lq = log(q)
        le = log(error)
        lne = log(1-error)

        E = {
            1: {0: logsum([lne + lp , le + 2*lp]), 1: le + log(2) + lp + lq, 
                    2: logsum([lne + lq , le + 2*lq])},
            0: {0: 2*lp, 1: log(2) + lp + lq, 2: 2*lq}
            }
       
    return E


#########################################################
### From Kristinas fwdbwd.py ###################
#########################################################

    
def FwdBwd(obsvec, posvec, Bfreqvec, error, a, f):
    """ Forward-Backward algorithm to compute posteriori probabilities. 
    
        A local optimization to find the most probable hidden state at a 
        given position t based on a observation sequence. 
        Assuming known initiation- , transition- and emission probabilities.
        
        Input to the function:
        obsvec = list of 0 = heterozygous, 1 = hom ref allele, 2 = hom alt allele
        posvec = CM marker positions
        Bfreqvec = population frequencies of second alleles

        error = the genotyping error-rate
        a = the distance between two states, exponential distributed
        f = inbreeding coefficient
        
        Returning the posteriori-probabilites for IBD at each position.        
    """

    states = (1, 0) # 1 = IBD, 0 = not-IBD
    n = len(posvec)
    distvec = [y-x for x,y in zip(posvec[:-1], posvec[1:])]
    
    # Probability matrices
    I = initiation(f, logspace=1)
    Evec = [emission(q, error, logspace=1) for q in Bfreqvec]
    Tvec = [transition(d, a, f, logspace=1) for  d in distvec]
    
    # FORWARD PART #
    lfwd = [{}]

    # Initiation
    for j in states:
        lfwd[0][j] = I[j] + Evec[0][j][obsvec[0]]

    # Iteration
    for t in range(1,n):
        lfwd.append({})
        for j in states:
           lfwd[t][j] = logsum([lfwd[t-1][i] + Tvec[t-1][i][j] + Evec[t][j][obsvec[t]] for i in states])


    # Termination
    lsum_fwd = logsum([lfwd[-1][j] for j in states])
    
    # BACKWARD PART #
    lbwd = [{} for x in range(n)]

    # Initialisation
    for j in states:
        lbwd[n-1][j] = log(1)

    # Iteration
    for t in reversed(range(n-1)): 
        for j in states:                                    
            lbwd[t][j] = logsum([lbwd[t+1][i] + Tvec[t][j][i] + Evec[t+1][i][obsvec[t+1]] for i in states])
                                     
    # POSTERIORI PART #
    # Posterior probabilities in log-space
    lprob = [{j : lfwd[t][j] + lbwd[t][j] - lsum_fwd for j in states} for t in range(n)]

    # Posterior probabilities exp() transformed   
    return [exp(lp[1]) for lp in lprob]
 

def read_HMMdata(infile):
    """ Reading a csv.file."""
    
    def gt(a,b):
        if a==b:  
            return 1 if a == '0' else 2
        return 0
        
    with open(infile, 'r') as fil: 
        next(fil) # skip header
        data = [ln.split() for ln in fil]	
    
    return zip(*[(float(pos), gt(a,b), float(bfr)) for pos,a,b,bfr in data])


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 6:
        raise Error("Syntax error. Too few arguments")
    infile = sys.argv[1]
    error = float(sys.argv[2])
    a = float(sys.argv[3])
    f = float(sys.argv[4])
    outfile = sys.argv[5]
    posvec, obsvec, Bfreqs = read_HMMdata(infile)
    postprobs = FwdBwd(obsvec, posvec, Bfreqs, error, a, f)  
    with open(outfile, 'w') as out:
        out.write('\n'.join(map(str, postprobs)))
    