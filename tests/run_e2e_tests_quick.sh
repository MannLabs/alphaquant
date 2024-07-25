#!/bin/bash
conda activate alphaquant
pip install pytest
pip install nbmake==1.5.3
echo "Running e2e tests quick"
pytest --nbmake e2e_tests_small/mixed_species.ipynb
pytest --nbmake e2e_tests_small/phospho.ipynb
conda deactivate
