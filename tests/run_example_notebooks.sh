#!/bin/bash
conda activate alphaquant
cd ../example_nbs
echo "Running example notebooks"
pytest --nbmake . --ignore=tree_example.ipynb
cd -

TEST_NBS=$(find ./nbs -name "*.ipynb")

pytest --nbmake $(echo $TEST_NBS)

conda deactivate
