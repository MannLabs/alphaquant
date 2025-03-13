# syntax=docker/dockerfile:1

FROM python:3.11

# Prevents Python from writing pyc files.
ENV PYTHONDONTWRITEBYTECODE=1
# Keeps Python from buffering stdout and stderr to avoid situations where
# the application crashes without emitting any logs due to buffering.
ENV PYTHONUNBUFFERED=1

WORKDIR /app

COPY requirements requirements
RUN pip install -r requirements/requirements.txt
RUN pip install -r requirements/requirements_gui.txt

COPY alphaquant alphaquant
COPY MANIFEST.in MANIFEST.in
COPY LICENSE LICENSE
COPY README.md README.md
COPY pyproject.toml pyproject.toml

RUN pip install ".[stable,gui-stable]"

ENV PORT=41215
EXPOSE 41215

CMD ["/usr/local/bin/alphaquant", "gui"]

# build & run:
# docker build --progress=plain -t alphaquant .
# DATA_FOLDER=/path/to/local/data
# docker run -p 41215:41215 -v $DATA_FOLDER:/app/data/ -t alphaquant