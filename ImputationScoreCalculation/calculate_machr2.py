# Adapted from: SNPSummaryComputation.cpp
# part of the Impute software?

import numpy as np
import gzip
import sys, os, math

input_file = sys.argv[1] # path to input VCF
output_file = sys.argv[2] # out file name
snp_file_path = sys.argv[3] # out file name for a SNP-R2 score overview

info = 0.0
impinfo = 0.0
levels = [0.0, 1.0, 2.0]
levelsSquared = [0.0, 1.0, 4.0]


def colsum(matrix, colid):
    sum = 0
    for i in range(0,len(matrix)):
        sum += float(matrix[i][colid])
    return sum

def calcsum(vector):
    sum = 0
    for val in vector:
        sum += float(val)
    return sum

def compVar(vector):
    mean = 0
    for i in range(0,len(levels)):
        mean += vector[i] * levels[i]
    variance = 0
    for i in range(0,len(levels)):
        variance += vector[i] * levelsSquared[i]
    variance -= (mean*mean)
    return variance

def compExpVar(matrix, dist):
    result = 0
    for i in range(0,len(matrix)):
        vect = matrix[i]
        colSum = calcsum(vect)
        indProbs = [
            float(vect[0]) + (1-colSum) * dist[0],
            float(vect[1]) + (1-colSum) * dist[1],
            float(vect[2]) + (1-colSum) * dist[2]
        ]
        v = compVar(indProbs)
        result += v
    
    output = result / len(matrix)
    if np.isnan(output):
        return 1
    return output

def computeAutosomal(variant): 
    probs = variant
    # variant has structure: [ind][allele], with 3 columns with probabilities (homAA, het, homBB) for each ind
    sumAA = colsum(variant,0)
    sumAB = colsum(variant,1)
    sumBB = colsum(variant,2)
    totalSum = sumAA + sumAB + sumBB
    theta_mle = (sumAB + 2.0 * sumBB) / (2*totalSum)

    if theta_mle == 0 or theta_mle == 1:
        return [None, None]

    imputeFallBackDist = [0.0,0.0,0.0]
    fallBackDist = [   ((1-theta_mle) * (1-theta_mle)),
                        (2 * theta_mle * (1 - theta_mle)),
                        (theta_mle * theta_mle)]
    infoscore = 1 - compExpVar(probs, fallBackDist) / (2 * theta_mle * (1 - theta_mle))
    infoscoreimputation = 1 - compExpVar(probs, imputeFallBackDist) / (2 * theta_mle * (1 - theta_mle))

    return [infoscore, infoscoreimputation]


line_counter = 0

with gzip.open(input_file, 'rt') as vcf_file, \
    gzip.open(output_file, 'wt') as out_file, \
    gzip.open(snp_file_path, 'wt') as snp_file:

    snp_file.write('snp\tinfo_score\tinfo_imputation_score\n')

    for line in vcf_file:
        line = line.strip()

        if line.startswith('#'):
            out_file.write(line + '\n')
            continue

        elems = line.split('\t')
        snp = elems[2]

        format_field = elems[8]

        if not 'GP' in format_field.split(':'):
            continue 

        gp_col_index = format_field.split(':').index('GP')
        
        gt_values = elems[9:]

        gt_probs = [gt.split(':')[gp_col_index].split(',') for gt in gt_values]
        gt_probs = [g for g in gt_probs if g != ['.']]

        matrix = np.array(np.array(gt_probs))

        info_score, info_score_imputation = computeAutosomal(matrix)
        
        # Skip variant if info score is None
        # This could happen when it is monomorphic
        if info_score is None:
            continue

        snp_file.write(f'{snp}\t{info_score}\t{info_score_imputation}\n')

        if info_score > 0.8 and info_score_imputation > 0.8:
            out_file.write(line + '\n')

        line_counter += 1

        if line_counter % 100 == 0:
            print(f'{line_counter} lines processed')

