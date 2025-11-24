import pandas as pd

data = []
current_title = None

with open("your_file.txt", "r") as f:
    for line in f:
        line = line.strip()

        # 标题行（以 "vs" 识别）
        if " vs " in line:
            current_title = line
            continue

        # 跳过分隔线和空行
        if not line or line.startswith("-") or line.startswith("Query"):
            continue

        # 解析表格行
        parts = line.split()
        if len(parts) == 4:
            query, target, identity, similarity = parts
            data.append({
                "Comparison": current_title,
                "Query": query,
                "Target": target,
                "Identity": float(identity),
                "Similarity": float(similarity)
            })

# 转为 DataFrame
df = pd.DataFrame(data)
print(df)
