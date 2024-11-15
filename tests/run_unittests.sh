#!/bin/bash
conda activate alphaquant
pip install pytest
cd unit_tests
python -m pytest -v
cd ..
conda deactivate