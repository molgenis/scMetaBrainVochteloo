import sys
import pandas as pd


if len(sys.argv) < 2:
    print("Usage: input.txt output.txt")
    exit()

inf = sys.argv[1]
ouf = sys.argv[2]

df = pd.read_csv(inf,sep='\t',header=0)
df = df.sort_values('PPH4',ascending=False)
df["deltaPPH4"] = 1.0-df["PPH4"]
rctr = 1
sum = 0
lowestpph4belowcutoff = None
nrsignificant = 0
fdrs = []
signf = []
for index, row in df.iterrows():
    sum += row["deltaPPH4"]
    mean = sum / rctr
    sign = False
    if mean < 0.05:
        lowestpph4belowcutoff = row["PPH4"]
        nrsignificant += 1
        sign = True
    fdrs.append(mean)
    signf.append(sign)
    # print(f"{row["PPH4"]}\t{row["deltaPPH4"]}\t{rctr}\t{mean}\t{sign}")

    rctr += 1

df["FDR"] = fdrs
df["Significant"] = signf

print(f"Lowest PPH4 with FDR<0.05: {lowestpph4belowcutoff}; nr significant: {nrsignificant}")

df.drop(columns=["deltaPPH4"], inplace=True)
df.to_csv(ouf,index=False,sep='\t')