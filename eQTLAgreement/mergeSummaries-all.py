import glob
import sys

indir=sys.argv[1]
outfile=sys.argv[2]


files = glob.glob(indir+"/*-summary.txt")

print(f"{indir} - {len(files)} files")
fctr = 0
fo = open(outfile,'w')



for file in files:
	fname = file.split("/")[-1]
	fname = fname.split("-vs-")[0]
	fh = open(file,'r')
	if fctr == 0:
		fo.write(fh.readline())
	else:
		fh.readline()
	for line in fh:
		if len(line.strip()) > 0:
			fo.write(line)
	fctr += 1
fo.close()

