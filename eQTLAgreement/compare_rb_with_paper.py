
import glob
files = glob.glob("./output/*-summary.txt")

data={}

for file in files:
    fh = open(file,'r')
    fh.readline()
    elems = fh.readline().strip().split("\t")
    cta = elems[0]
    ctb = elems[1]
    rb = elems[5]
    n = elems[3]

    data[cta+"-"+ctb] = [rb,n]
    # CT1	CT2	NrSignCT1	NrOverlapCT2	AC	Rb	Rb_se
    fh.close()

fho = open("ComparisonWithPaper.txt",'w')
fh = open("/path/to/scMetaBrain-Vochteloo-Freeze2-Rb.txt",'r')
# Discovery cell type	Replication cell type	Number of overlapping eQTL	AC	Rb	Rb standard error
fh.readline()
fho.write("ctA\tctB\tRbPaper\tRbNew\tnPaper\tnNew\n")
for line in fh:
    elems = line.rstrip().split("\t")
    cta = elems[0]
    ctb = elems[1]
    rb = elems[-2]
    n = elems[2]
    other = data.get(cta+"-"+ctb)
    if other is not None:
        outln = f"{cta}\t{ctb}\t{rb}\t{other[0]}\t{n}\t{other[1]}\n"
        fho.write(outln)
fho.close()