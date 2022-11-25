import pandas as pd
import numpy as np
import re
df = pd.read_csv('Var.csv')
name = pd.read_csv('names.csv')
name2=name["Name"].to_list()
n_names=df.index[df["AbsoluteFrequency"]=="AbsoluteFrequency"]
NAME=[]

for x in range(n_names[0]):
    NAME.append(name2[0])

for i in range(len(n_names)-1):
    x1=n_names[i]
    x2=n_names[i+1]
    for ii in range(x1,x2):
        NAME.append(name2[i+1])

for i in range(len(df)-x2):
    NAME.append(name2[len(name2)-1])

df["Name"]=NAME
df.dropna(subset = ["AbsoluteFrequency"], inplace=True)
df=df[df.AbsoluteFrequency!="AbsoluteFrequency"]

Fish=[]
for i in range(0,len(df)):
    if bool(re.match(r'KI[0-9]+[ZW][E][G]*-',df.Name.iloc[i])):
        a=re.split(r'KI[0-9]+[ZW][E][G]*-',df.Name.iloc[i])
        Fish.append(a[1])
    else:
        Fish.append(float("NAN"))
df["Fish"] = Fish
df=df.dropna()

df["Group"]=df['Name'].str.extract(r'KI[0-9]+([ZW][E][G]*)*')
df.NrOfReads=df['NrOfReads'].astype('int')
df = df.drop(df[df.NrOfReads < 100].index)
df_WE=df[df.Group == "WE"]
df_ZEG=df[df.Group == "ZEG"]
ZEG_U=df_ZEG.Fish.unique()
WE_U=df_WE.Fish.unique()
AF=np.hstack((ZEG_U,WE_U))
BF=[x for x in AF if list(AF.flatten()).count(x)==1]
df=df[~df['Fish'].isin(BF)]
df["ID"] = df["Fish"] + df["TypeLengthFlankingSequence (+ strand)"]


df=df.rename(columns={"TypeLengthFlankingSequence (+ strand)": "Seq"})
df=df.rename(columns={"NrOfReads": "NofReads"})
df=df.rename(columns={"AbsoluteFrequency": "Abs"})
df.NofReads=df['NofReads'].astype('int')
df.Abs=df['Abs'].astype('int')
df["Rel"] =df.Abs/df.NofReads

df.to_csv('df.csv')

df_WE=df[df["Group"]=="WE"]
US=pd.unique(df_WE.Seq)
Rus=[]
for i in range(0,len(US)):
    Rus.append(sum(df_WE.Rel[df_WE.Seq == US[i]])/len(df.Fish.unique()))

FUS=pd.DataFrame({"Seq":US,"Rel":Rus})
FUS = FUS.sort_values(by = "Rel",ascending=False)
FUS =FUS.reindex()
FUS =FUS.iloc[0:10]
df2=df[df['Seq'].isin(FUS.Seq)]

U=pd.unique(df2.ID)
fdf = pd.DataFrame({'ID' : [], 'Fish' : [], "WE_Abs" : [], "ZEG_Abs" : [], "WE_Rel" : [], "ZEG_Rel" : [],"WE_NR" : [], "ZEG_NR" : []})
for i in range(0,len(U)):
    T=df2[df2.ID == U[i]]
    if len(T) == 2:
        d = {'ID': T.ID.iloc[0], 'Fish': T.Fish.iloc[0], "WE_Abs": int(sum(T.Abs[T.Group == "WE"])), "ZEG_Abs": int(sum(T.Abs[T.Group == "ZEG"])), "WE_Rel": float(sum(T.Rel[T.Group == "WE"])), "ZEG_Rel": float(sum(T.Rel[T.Group == "ZEG"])), "WE_NR": int(sum(T.NofReads[T.Group == "WE"])), "ZEG_NR": int(sum(T.NofReads[T.Group == "ZEG"]))}
    elif T.Group.iloc[0] == "WE":
        V = df[df.Fish == T.Fish.iloc[0]]
        V=V[V.Group == "ZEG"]
        d = {'ID': T.ID.iloc[0], 'Fish': T.Fish.iloc[0], "WE_Abs": int(sum(T.Abs[T.Group == "WE"])), "ZEG_Abs": int(sum(T.Abs[T.Group == "ZEG"])), "WE_Rel": float(sum(T.Rel[T.Group == "WE"])),"ZEG_Rel": float(sum(T.Rel[T.Group == "ZEG"])), "WE_NR": int(sum(T.NofReads[T.Group == "WE"])),"ZEG_NR": int(V.NofReads.iloc[0])}
    else:
        V = df[df.Fish == T.Fish.iloc[0]]
        V = V[V.Group == "WE"]
        d = {'ID': T.ID.iloc[0], 'Fish': T.Fish.iloc[0], "WE_Abs": int(sum(T.Abs[T.Group == "WE"])), "ZEG_Abs": int(sum(T.Abs[T.Group == "ZEG"])), "WE_Rel": float(sum(T.Rel[T.Group == "WE"])),"ZEG_Rel": float(sum(T.Rel[T.Group == "ZEG"])), "WE_NR": int(V.NofReads.iloc[0]),"ZEG_NR": int(sum(T.NofReads[T.Group == "ZEG"]))}
    fdf = fdf.append(d, ignore_index=True)
fdf["Seq"]=fdf['ID'].str.extract(r'[ssODN]*[1-9]*[-]*[0-9]+([DI][EN][LS][0-9]+[A\]\[CTG]+)')



fdf.to_csv('fdf.csv')

