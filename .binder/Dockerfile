# adapted from https://github.com/binder-examples/minimal-dockerfile
FROM python:3.11-slim
# install the notebook package
RUN pip install --no-cache --upgrade pip && \
    pip install --no-cache notebook jupyterlab

# create user with a home directory
ARG NB_USER=x
ARG NB_UID=1000
ENV USER ${NB_USER}
ENV HOME /home/${NB_USER}

RUN adduser --disabled-password \
    --gecos "Default user" \
    --uid ${NB_UID} \
    ${NB_USER}
WORKDIR ${HOME}
USER ${USER}


COPY requirements requirements
RUN pip install -r requirements/requirements.txt
RUN pip install -r requirements/requirements_gui.txt

COPY alphaquant alphaquant
COPY MANIFEST.in MANIFEST.in
COPY LICENSE LICENSE
COPY README.md README.md
COPY pyproject.toml pyproject.toml

RUN pip install ".[stable,gui-stable]"
