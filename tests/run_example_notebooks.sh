#!/bin/bash
conda activate alphaquant
echo "Running example notebooks"
pytest --nbmake example_nbs/
conda deactivate
