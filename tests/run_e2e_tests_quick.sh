#!/bin/bash
conda activate alphaquant
pip install pytest
pip install nbmake==1.5.3
echo "Running e2e tests quick"
pytest --nbmake quicktests/mixed_species.ipynb
pytest --nbmake quicktests/phospho.ipynb
conda deactivate
