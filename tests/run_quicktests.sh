conda activate alphaquant
pip install pytest
pip install nbmake==1.5.3
python download_testfiles.py quicktests
echo "Running quicktests"
pytest --nbmake quicktests/mixed_species.ipynb
pytest --nbmake quicktests/phospho.ipynb
conda deactivate
