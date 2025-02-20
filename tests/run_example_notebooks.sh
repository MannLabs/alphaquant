#!/bin/bash
conda activate alphaquant
cd ../example_nbs
echo "Running example notebooks"
pytest --nbmake .
conda deactivate
