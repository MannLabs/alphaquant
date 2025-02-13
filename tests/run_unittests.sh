#!/bin/bash
conda activate alphaquant
cd unit_tests
python -m pytest -v
cd ..
conda deactivate
