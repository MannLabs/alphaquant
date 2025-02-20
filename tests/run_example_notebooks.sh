#!/bin/bash
conda activate alphaquant
echo "Running example notebooks"
pytest --nbmake example_nbs/*ipynb
conda deactivate
