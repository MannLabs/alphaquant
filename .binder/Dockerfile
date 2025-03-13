# adapted from https://github.com/binder-examples/minimal-dockerfile
FROM python:3.11-slim
# install the notebook package
RUN pip install --no-cache --upgrade pip && \
    pip install --no-cache notebook jupyterlab

# create user with a home directory
ARG NB_USER
ARG NB_UID
ENV USER ${NB_USER}
ENV HOME /home/${NB_USER}

RUN adduser --disabled-password \
    --gecos "Default user" \
    --uid ${NB_UID} \
    ${NB_USER}
WORKDIR ${HOME}

# copying requirements first to leverage caching
COPY requirements requirements
RUN pip install -r requirements/requirements.txt

COPY alphaquant alphaquant
COPY MANIFEST.in MANIFEST.in
COPY LICENSE LICENSE
COPY README.md README.md
COPY pyproject.toml pyproject.toml

RUN pip install ".[stable]"

# change owner to be able to write to that folder
COPY --chown=${USER}:${USER} example_nbs example_nbs

# build & debug:
# docker build --progress=plain --build-arg="NB_USER=user" --build-arg="NB_UID=1000" -f .binder/Dockerfile -t alphaquant-binder .
# docker run -it alphaquant-binder bash
