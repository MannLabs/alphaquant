set -e -u

INSTALL_TYPE=$1 # stable, loose, etc..
ENV_NAME=${2:-alphaquant}
PYTHON_VERSION=${3:-3.11}
INSTALL_MONO=${4:-false}


if [ "$INSTALL_MONO" = "true" ]; then
  conda create -n $ENV_NAME python=$PYTHON_VERSION mono -y
else
  conda create -n $ENV_NAME python=$PYTHON_VERSION -y
fi

if [ "$INSTALL_TYPE" = "loose" ]; then
  INSTALL_STRING=""
else
  INSTALL_STRING="[${INSTALL_TYPE}]"
fi

# print pip environment for reproducibility
conda run -n $ENV_NAME --no-capture-output pip freeze

# conda 'run' vs. 'activate', cf. https://stackoverflow.com/a/72395091
conda run -n $ENV_NAME --no-capture-output pip install -e "../.$INSTALL_STRING"
conda run -n $ENV_NAME --no-capture-output alphaquant -v
