import pandas as pd

gistic = (
    pd.read_table("gistic_gene_matrix.tsv")      # 1️⃣ your file
      .set_index("Hugo_Symbol")                  # rows = genes
      .sort_index()
)

# 2️⃣ Which samples?  (values EXACTLY == -2)
brca2_row = gistic.loc["BRCA2"]
homdel_samples = brca2_row[brca2_row == -2].index.tolist()

print(f"{len(homdel_samples)} samples with BRCA2 -2")

from mygene import MyGeneInfo
mg = MyGeneInfo()

# Cache to avoid repeated API calls
coord_cache = "chr13_coords.parquet"

try:
    coords = pd.read_parquet(coord_cache)

except FileNotFoundError:
    anno = mg.querymany(
        gistic.index.tolist(), scopes="symbol", fields="chrom,start,end", species="human", as_dataframe=True
    )
    # Keep only GRCh38 chr13 hits; drop duplicates, sort by position
    coords = (
        anno.query('chrom == "13" & not start.isna()')
             .reset_index()
             .groupby("symbol", as_index=False)
             .agg(start=("start", "min"), end=("end", "max"))
             .sort_values("start")
             .rename(columns={"symbol": "Hugo_Symbol"})
             .set_index("Hugo_Symbol")
    )
    coords.to_parquet(coord_cache)

# Merge coordinates with the GISTIC values for chr13 genes only
chr13 = coords.join(gistic, how="left")


results = []

# Absolute BRCA2 position for quick lookup
brca2_pos = chr13.loc["BRCA2", "start"]

for sample in homdel_samples:
    calls = chr13[sample]          # pandas Series (index = Hugo, values = CNV)
    is_del = calls == -2

    # Index of BRCA2 row
    brca2_idx = chr13.index.get_loc("BRCA2")

    # Walk left
    left_idx = brca2_idx
    while left_idx > 0 and is_del.iat[left_idx - 1]:
        left_idx -= 1

    # Walk right
    right_idx = brca2_idx
    while right_idx < len(chr13) - 1 and is_del.iat[right_idx + 1]:
        right_idx += 1

    left_gene   = chr13.index[left_idx]
    right_gene  = chr13.index[right_idx]
    left_pos    = chr13.iloc[left_idx]["start"]
    right_pos   = chr13.iloc[right_idx]["end"]

    span_bp     = right_pos - left_pos

    results.append(
        {
            "sample": sample,
            "left_gene": left_gene,
            "right_gene": right_gene,
            "gene_count": right_idx - left_idx + 1,
            "span_bp": span_bp,
        }
    )

summary = pd.DataFrame(results).sort_values("span_bp", ascending=False)


### Plotting

import matplotlib.pyplot as plt

plt.figure(figsize=(6,4))
plt.hist(summary["span_bp"]/1e6, bins=20)
plt.xlabel("BRCA2 homozygous-deletion span (Mb)")
plt.ylabel("# samples")
plt.title("Extent of chr13 -2 blocks containing BRCA2")
plt.tight_layout()
plt.show()

