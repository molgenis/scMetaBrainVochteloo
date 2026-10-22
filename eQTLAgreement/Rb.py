import math as m
import statistics as s
from scipy.stats.distributions import chi2
# import concurrent.futures

def squareVal(val):
    return val*val

def squareList(list):
    return [i ** 2 for i in list]

def cov(a,b):
    meana = s.mean(a)
    meanb = s.mean(b)
    return covWMean(a,b,meana,meanb)

def covWMean(a, b, meana, meanb):
    sigma = 0
    for i in range(0,len(a)):
        sigma += (a[i]-meana)*(b[i]-meanb)
    return sigma / (len(a)-1)

def jacknife(b1,se1,b2,se2,theta,k):
    b1_jack = b1.copy()
    b2_jack = b2.copy()
    se1_jack = se1.copy()
    se2_jack = se2.copy()
    theta_jack = None
    if theta is not None:
        theta_jack = theta.copy()
        theta_jack.pop(k)

    # print(f"{k} - {n}")
    b1_jack.pop(k)
    b2_jack.pop(k)
    se1_jack.pop(k)
    se2_jack.pop(k)

    se1_jack_sq = squareList(se1_jack)
    se1_jack_sq_mean = s.mean(se1_jack_sq)

    se2_jack_sq = squareList(se2_jack)
    se2_jack_sq_mean = s.mean(se2_jack_sq)
    
    var_b1_jack = s.variance(b1_jack) - se1_jack_sq_mean
    var_b2_jack = s.variance(b2_jack) - se2_jack_sq_mean
    
    if var_b1_jack < 0:
        var_b1_jack = s.variance(b1_jack)
    
    if var_b2_jack < 0:
        var_b2_jack = s.variance(b2_jack)
    
    cov_e1_jack_e2_jack = 0
    if theta is not None:
        cov_e1_jack_e2_jack = s.mean(theta_jack) * m.sqrt( se1_jack_sq_mean * se2_jack_sq_mean )

    cov_b1_b2_jack = cov(b1_jack, b2_jack) - cov_e1_jack_e2_jack       
    r_tmp = cov_b1_b2_jack / m.sqrt(var_b1_jack * var_b2_jack)

    return r_tmp

# This function is slightly adapted from the script shared by Ting Qi.
def calcu_cor_true(b1, se1, b2, se2, theta, threads = 1, skip_jacknife=False):
    if len(b1) != len(b2) or len(b1) != len(se1) or len(b2) != len(se2):
        print("Error Rb in calculation: vectors not the same length")
        exit()
    if len(b1) == 0:
        return [None, None, None]
    # remove infinite, none and nan
    idx = []
    for i in range(0,len(b1)):
        if b1[i] is None or b2[i] is None or se1[i] is None or se2[i] is None:
            idx.append(i)
        if m.isnan(b1[i]) or m.isnan(b2[i]) or m.isnan(se1[i]) or m.isnan(se2[i]):
            idx.append(i)
        elif m.isinf(b1[i]) or m.isinf(b2[i]) or m.isinf(se1[i]) or m.isinf(se2[i]):
            idx.append(i)

    if len(idx) > 0:
        idx.sort(reverse=True)
        for i in idx:
            b1.pop(i)
            b2.pop(i)
            se1.pop(i)
            se2.pop(i)
            if theta is not None:
                theta.pop(i)
    if len(b1) < 1:
        return [None, None, None]

    se1sq = squareList(se1)
    se2sq = squareList(se2)
    se1sqmean = s.mean(se1sq)
    se2sqmean = s.mean(se2sq)
    var_b1 = s.variance(b1) - se1sqmean
    var_b2 = s.variance(b2) - se2sqmean

    print(var_b1)
    print(var_b2)
    if var_b1<0:
        var_b1 = s.variance(b1)
    if var_b2<0:
        var_b2 = s.variance(b2)
    print(var_b1)
    print(var_b2)
    
    cov_b1_b2 = 0
    if theta is not None:

        cov_b1_b2 = cov(b1,b2) - s.mean(theta) *  m.sqrt(  se1sqmean * se2sqmean )
    else:
        cov_b1_b2 = cov(b1,b2)
    r = cov_b1_b2 / m.sqrt(var_b1 * var_b2)

    # jacknife
    if skip_jacknife:
        return [r,-1,-1]
    else:
        r_jack = []
        n = len(b1)
        # futures = []
        # print(f"nrthreads: {threads}")
        # the below can be multithreaded, but is still quite slow in python for some reason
        for k in range(0,n):
            r_tmp = jacknife(b1,se1,b2,se2,theta,k)
            if not m.isnan(r_tmp):
                r_jack.append(r_tmp)
        
        r_jack_mean = s.mean(r_jack)
        # print(f"rmeanjack: {r_jack_mean}")
        r_jack_sumsq = 0
        for r_jack_v in r_jack:
            delta = (r_jack_v - r_jack_mean)
            r_jack_sumsq += delta*delta
        # print(f"rjacksumsq: {r_jack_sumsq}")
        se_r = m.sqrt( (n - 1) / n * r_jack_sumsq )
        
        chi = r / se_r
        chisq = chi * chi
        p = chi2.sf(chisq, 1)
        if r > 1:
            r=1
        if r < -1:
            r=-1
        return [r,se_r,p]

    



