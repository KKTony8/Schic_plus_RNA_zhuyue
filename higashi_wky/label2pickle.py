import pandas as pd
import pickle

# labels.txt 文件路径
txt_file = "/data5/GPT/Wuky/Higashi/Higashi/WKY_data/labels.txt"

# 输出 pickle 文件
pickle_file = "/data5/GPT/Wuky/Higashi/Higashi/WKY_data/label_info.pickle"

# 读取 txt
df = pd.read_csv(txt_file, sep="\t")

# 修改 key 名，保持和 demo 一致
label_dict = {
    "cell type": df["group"].tolist(),  # 把 group 改成 cell type
    "batch": ["batch1"] * len(df)       # 如果没有 batch 信息，全部设为 batch1
}

# 保存为 pickle
with open(pickle_file, "wb") as f:
    pickle.dump(label_dict, f)

print(f"✅ 已生成 pickle 文件: {pickle_file}")

# 测试读取
with open(pickle_file, "rb") as f:
    data = pickle.load(f)

print("\n字典 keys:", data.keys())
print("前 5 个 cell type:", data['cell type'][:5])
print("前 5 个 batch:", data['batch'][:5])

