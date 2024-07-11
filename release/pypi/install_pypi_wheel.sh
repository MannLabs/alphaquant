conda create -n alphaquant_pip_test python=3.10 -y
conda activate alphaquant_pip_test
pip install "alphaquant[stable]"
alphaquant
conda deactivate
