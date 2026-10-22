import gzip
import sys

if len(sys.argv) < 6:
    print("Usage: f1 label1 f2 label2 out")
    exit()

f1 = sys.argv[1]
label1 = sys.argv[2]
if label1 == "None":
    label1 = None
f2 = sys.argv[3]
label2 = sys.argv[4]
if label2 == "None":
    label2 = None
f3 = sys.argv[5]

print()
print("Concatenating summaries")
print(f1)
print(f2)
print(f3)

print(f1)
fh = open(f1,'r')
header1 = fh.readline().rstrip()
data = {}
for line in fh:
    elems = line.strip().split("\t")
    
    ct1 = elems[0].strip()
    ct2 = elems[1].strip()
    pair = ct1+"_"+ct2
    data[pair] = line.strip()

fh.close()


fho = open(f3,'w')



headerout = []
helems1 = header1.strip().split("\t")
for e in helems1:
    if label1 is None:
        headerout.append(e)
    else:
        headerout.append(label1+"_"+e)

fh2 = open(f2,'r')
header2 = fh2.readline().rstrip()
helems2 = header2.strip().split("\t")
for e in helems2:
    if label2 is None:
        headerout.append(e)
    else:
        headerout.append(label2+"_"+e)




headerout = "\t".join(headerout)+"\n"
print(headerout)
fho.write(headerout)


for line in fh2:
    
    elems = line.strip().split("\t")
    
    ct1 = elems[0].strip()
    ct2 = elems[1].strip()
    
    pair = ct1+"_"+ct2
    # print("f2:\t"+pair)
    otherline = data.get(pair)
   
    if otherline is not None:
        # print(otherline)
        outln = otherline+"\t"+line
        fho.write(outln)
fh.close()
fho.close()