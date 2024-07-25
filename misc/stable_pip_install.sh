conda create -n alphaquant python=3.9 -y
conda activate alphaquant
pip install -e '../.[stable,development-stable]'
alphaquant
conda deactivate
