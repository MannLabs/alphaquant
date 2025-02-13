conda create -n alphaquant python=3.11 -y
conda activate alphaquant
pip install -e '../.[development, gui]'
alphaquant
conda deactivate
