conda create -n alphaquant python=3.8 -y
conda activate alphaquant
pip install -e '../.[development]'
alphaquant
conda deactivate
