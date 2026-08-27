import sys, pandas as pd, numpy as np
TCGA_DIR, COHORT = sys.argv[1], sys.argv[2]
def norm_id(s):
    p=str(s).split("-"); return "-".join(p[:4])[:15] if len(p)>=4 else str(s)
cna = pd.read_csv(f"{TCGA_DIR}/data_cna.txt", sep="\t", low_memory=False)
cna = cna[cna.Hugo_Symbol.notna()].drop_duplicates("Hugo_Symbol").set_index("Hugo_Symbol")
cna = cna.drop(columns=["Entrez_Gene_Id"], errors="ignore")
# neighbours on 15q (B2M) and 16q (NLRC5)
sets = {"B2M_15q": ["B2M","TRIM69","RAB27A","SPG11","PATL2","TGM5"],
        "NLRC5_16q":["NLRC5","CYLD","CBLN1","NKD1","SNX20","CIITA"]}
for name, genes in sets.items():
    g = [x for x in genes if x in cna.index]
    if len(g) < 2: continue
    loss = (cna.loc[g] < 0).mean(axis=1)
    print(f"[{COHORT}] {name}")
    print(loss.to_string(float_format=lambda x: f"{x:.3f}"))
    print()
