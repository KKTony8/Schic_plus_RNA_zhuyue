import pandas as pd
import pickle

# labels.txt 文件路径
txt_file = "/data5/GPT/Wuky/Higashi/Higashi/Monocyte_data/labels.txt"

# 输出 pickle 文件
pickle_file = "/data5/GPT/Wuky/Higashi/Higashi/Monocyte_data/label_info.pickle"

# 读取 txt
# 假设 labels.txt 是两列 tab 分隔: cell_id group
df = pd.read_csv(txt_file, sep="\t")

# 将 DataFrame 转换为 dict，每列对应一个 key，值是 list
label_dict = {col: df[col].tolist() for col in df.columns}

# 保存为 pickle
with open(pickle_file, "wb") as f:
    pickle.dump(label_dict, f)

print(f"✅ 已生成 pickle 文件: {pickle_file}")

# 测试读取
with open(pickle_file, "rb") as f:
    data = pickle.load(f)
    
print("\n字典 keys:", data.keys())
print("前 5 个 cell_id:", data['cell_id'][:5])
print("前 5 个 group:", data['group'][:5])

