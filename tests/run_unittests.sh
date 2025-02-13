conda activate alphaquant

coverage run --source=../alphaquant -m pytest -k 'not slow'
conda deactivate
