import gzip
import sys
import statsmodels.api as sm 
import pandas as pd    
from scipy import stats    
import numpy as np

if len(sys.argv) < 5:
    print("Usage: summaryfile mashrfile output suffix")
    exit()

print("merging in mashr info")

f1 = sys.argv[2]
f2 = sys.argv[1]
f3 = sys.argv[3]
suffix = sys.argv[4]


mashrdata = {}
fh = open(f1,'r')
header = fh.readline().rstrip().split("\t")
for line in fh:
    elems = line.rstrip().split("\t")
    ct2 = elems[0]
    for i in range(1,len(elems)):
        ct1 = header[i]
        pair = ct1+"_"+ct2
        mashrdata[pair] = elems[i]

fh.close()

fh = open(f2,'r')
fho = open(f3,'w')
origheader = fh.readline().rstrip()
header = origheader+f"\tMASHR_{suffix}\n"
fho.write(header)

print(f2)

colnames = origheader.split("\t")
print(colnames)


alldata = {}
numeric = set()
for line in fh:
    elems = line.strip().split("\t")
    
    ct1 = elems[0].strip()
    ct2 = elems[1].strip()

    lnout = line.rstrip()
    pair = ct1+"_"+ct2
    other = mashrdata.get(pair)
    lnout += f"\t{other}\n"
    fho.write(lnout)

    
    pairdata = {}
    pairdata["mashr"] = float(other)
    numeric.add("mashr")

    
    for i in range(2,len(elems)):
        try:
            colname = colnames[i]
            
            pairdata[colname] = float(elems[i])
            numeric.add(colname)
        except Exception as e:
            print(colname)
            pass

    alldata[pair] = pairdata
fho.close()
fh.close()

fho = open(f3+"_correlated.txt",'w')
numeric = list(numeric)
numeric.sort()
header = "-"
for statA in numeric:
    header += "\t"+statA

fho.write(header+"\n")

# correlate the stats
for statA in numeric:
    correls = []
    outln = statA
    for statB in numeric:
        vals1 = []
        vals2 = []
        for pair in alldata.keys():
            cts = pair.split("_")
            data = alldata.get(pair)
            if cts[0] != cts[1]:
                v1 = data[statA]
                v2 = data[statB]
                vals1.append(v1)
                vals2.append(v2)
        # correlate

        # print(f"{len(vals1)}  - {len(vals2)}")

        correl = stats.pearsonr(vals1,vals2)

        # OLS
        x = vals2
        x = sm.add_constant(x)
        result = sm.OLS(vals1, x).fit()
        # result = "OK"
        print(f"{statA} - {statB} - {correl.statistic} - {result.params}")
        outln += f"\t{correl.statistic}"
    fho.write(outln+"\n")
fho.close()
    
