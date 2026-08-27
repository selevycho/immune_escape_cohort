#!/usr/bin/env python3
"""
Is expression LOWER in mutated samples than in non-mutated samples of the SAME gene?
This separates real silencing from 'mutations land in low-expressed genes'.
"""
import sys, pandas as pd, numpy as np
from scipy import stats
TCGA_DIR, OUT_DIR, COHORT = sys.argv[1], sys.argv[2], sys.argv[3]

def norm_id(s):
    p=str(s).split("-"); return "-".join(p[:4])[:15] if len(p)>=4 else str(s)

expr = pd.read_csv(f"{TCGA_DIR}/data_mrna_seq_v2_rsem.txt", sep="\t", low_memory=False)
expr = expr[expr.Hugo_Symbol.notna()].drop_duplicates("Hugo_Symbol").set_index("Hugo_Symbol")
expr = expr.drop(columns=["Entrez_Gene_Id"], errors="ignore")
expr.columns = [norm_id(c) for c in expr.columns]
expr = expr.loc[:, ~expr.columns.duplicated()]

maf = pd.read_csv(f"{TCGA_DIR}/data_mutations.txt", sep="\t", low_memory=False,
                  usecols=["Hugo_Symbol","Tumor_Sample_Barcode","Variant_Classification"])
NONSYN={"Missense_Mutation","Nonsense_Mutation","Frame_Shift_Del","Frame_Shift_Ins",
        "In_Frame_Del","In_Frame_Ins","Splice_Site","Nonstop_Mutation","Translation_Start_Site"}
maf = maf[maf.Variant_Classification.isin(NONSYN)]
maf["sample"] = maf.Tumor_Sample_Barcode.map(norm_id)

rows=[]
counts = maf.Hugo_Symbol.value_counts()
for g in counts[counts>=20].index:
    if g not in expr.index: continue
    mut = set(maf.loc[maf.Hugo_Symbol==g,"sample"]) & set(expr.columns)
    non = set(expr.columns) - mut
    if len(mut)<10 or len(non)<10: continue
    a = np.log2(expr.loc[g, list(mut)].clip(lower=0)+1)
    b = np.log2(expr.loc[g, list(non)].clip(lower=0)+1)
    u,p = stats.mannwhitneyu(a,b,alternative="two-sided")
    pool = np.sqrt((a.var()+b.var())/2)
    rows.append({"gene":g,"n_mut":len(mut),"n_non":len(non),
                 "mean_mut":a.mean(),"mean_non":b.mean(),
                 "diff":a.mean()-b.mean(),
                 "cohens_d":(a.mean()-b.mean())/pool if pool>0 else np.nan,"p":p})

res=pd.DataFrame(rows).sort_values("cohens_d")
res["p_bonf"]=(res.p*len(res)).clip(upper=1.0)
res["cohort"]=COHORT
res.to_csv(f"{OUT_DIR}/{COHORT}_line1_control.tsv",sep="\t",index=False)

print(f"=== {COHORT.upper()} — mutated vs non-mutated, same gene ===")
print(f"genes tested: {len(res)}")
print(f"median cohens_d: {res.cohens_d.median():.3f}")
n_low = int(((res.p_bonf<0.05) & (res.loc[:,"diff"]<0)).sum())
print("genes with significantly LOWER expression when mutated:", n_low)
n_high = int(((res.p_bonf<0.05) & (res.loc[:,"diff"]>0)).sum())
print("genes with significantly HIGHER expression when mutated:", n_high)
print()
print("most down-regulated when mutated:")
print(res.head(10).to_string(index=False,float_format=lambda x:f"{x:.3f}"))
