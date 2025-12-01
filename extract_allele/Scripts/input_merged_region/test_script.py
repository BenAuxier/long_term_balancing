import pandas as pd

# 读入数据
df = pd.read_excel("yourfile.xlsx")  # 或 read_csv

# 除 assembly_id 外作为分组键的列
keys = ["genomic_region", "assembly_label", "gene",
        "domain_name", "ipr_name", "go_terms"]

# 分组合并并统计数量
df2 = (
    df.groupby(keys, dropna=False)
      .agg({
          "assembly_id": list,
      })
      .reset_index()
)

# 新增统计列：合并后 assembly 数量
df2["assembly_number"] = df2["assembly_id"].apply(len)

# 导出
df2.to_excel("merged_with_count.xlsx", index=False)
