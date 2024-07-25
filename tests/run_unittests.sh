#!/bin/bash
conda activate alphaquant
pip install pytest
cd unit_tests
pytest test_*.py -v
cd ..
conda deactivate