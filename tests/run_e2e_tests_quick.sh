#!/bin/bash
conda activate alphaquant
echo "Running e2e tests quick"
pytest --nbmake e2e_tests_small/mixed_species.ipynb
pytest --nbmake e2e_tests_small/phospho.ipynb
conda deactivate
