cd /data5/GPT/Wuky/
mkdir Higashi
git clone https://github.com/ma-compbio/Higashi.git
cd Higashi
conda create -n higashi python=3.8 -y
conda activate higashi
# 安装必要的库
pip install numpy==1.19.2 pandas==1.1.3 h5py==2.10.0 scikit-learn==0.23.2 fbpca==1.0 tqdm==4.50.2 matplotlib seaborn umap-learn cooler
pip install bokeh==2.1.1 pillow==7.2.0 cachetools cmocean
# 安装pytorch
nvidia smi
conda install pytorch torchvision torchaudio pytorch-cuda=11.7 -c pytorch -c nvidia

# jupyter的用法
conda activate higashi
tmux new -s jupyter
jupyter notebook --ip=服务器ip --port=1888 --no-browser --NotebookApp.token='' --NotebookApp.password=''

tmux a -t jupiter

# 在自己浏览器输入url
# 在浏览器便可以打开.ipy文件进行操作
# 下载好示例data,完成好config即可跑示例脚本
"/data5/GPT/Wuky/Higashi/Higashi/tutorials/4DN_sci-Hi-C_Kim et al.ipynb"

