import glob
import sys

indir=sys.argv[1]
outfile=sys.argv[2]


files = glob.glob(indir+"/*-summary.txt")

print(f"{indir} - {len(files)} files")
fctr = 0
fo = open(outfile,'w')


disc = ["EX", "IN", "OLI", "OPC", "MIC", "AST", "END", "PER"]
discset = set()
for ct in disc:
	discset.add(ct)


for file in files:
	fname = file.split("/")[-1]
	fname = fname.split("-vs-")[0]
	if fname in  discset:
		try:
			fh = open(file,'r')
			if fctr == 0:
				fo.write(fh.readline())
			else:
				fh.readline()
			for line in fh:
				if len(line.strip()) > 0:
					fo.write(line)
			fctr += 1
		except Exception as e:
			print("Issue with file: "+fname)
			print(repr(e))
fo.close()

