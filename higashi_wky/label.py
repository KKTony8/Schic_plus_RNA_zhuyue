import os
import glob
import pandas as pd

# allValidPairs 文件所在目录
input_dir = "/data5/GPT/Wuky/Higashi/Higashi/demo_data/WKYmonocyteQtxt/"
# 输出 label 文件
output_file = "/data5/GPT/Wuky/Higashi/Higashi/Monocyte_data/labels.txt"

# 定义 group 对应关系
group_map = {
    "M": "middle",  # 中年
    "O": "old",     # 老年
    "Y": "young"    # 年轻
}

# 获取所有 allValidPairs 文件
all_files = glob.glob(os.path.join(input_dir, "*.allValidPairs"))

# 按固定顺序排序：Y -> M -> O
def sort_key(f):
    name = os.path.basename(f)
    prefix = name[0]
    order = {'Y': 0, 'M': 1, 'O': 2}
    return order.get(prefix, 99), name

all_files_sorted = sorted(all_files, key=sort_key)

groups = []

for f in all_files_sorted:
    file_name = os.path.basename(f)
    cell_id = file_name.replace(".allValidPairs", "")
    prefix = cell_id[0]
    group = group_map.get(prefix, "unknown")
    groups.append(group)

# 构建 DataFrame：cell_id 为连续数字索引
df_label = pd.DataFrame({
    "cell_id": range(len(groups)),
    "group": groups
})

# 保存为 tab 分隔 txt
df_label.to_csv(output_file, sep="\t", index=False)
print(f"✅ labels.txt 已生成: {output_file}")

