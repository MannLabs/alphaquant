conda create -n alphaquant python=3.9 -y
conda activate alphaquant
pip install -e '../.[development, gui]'
alphaquant
conda deactivate
