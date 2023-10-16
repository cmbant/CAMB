
#Dockerfile for running pycamb notebooks with binder
#https://mybinder.org/v2/gh/cmbant/camb/master?filepath=docs%2FCAMBdemo.ipynb

FROM cmbant/cosmobox:latest


RUN pip install --no-cache --upgrade pip && \
    pip install --no-cache notebook

ARG NB_USER
ARG NB_UID
ENV NB_USER jovyan
ENV NB_UID 1000
ENV USER ${NB_USER}
ENV HOME /home/${NB_USER}

RUN adduser --disabled-password \
    --gecos "Default user" \
    --uid ${NB_UID} \
    ${NB_USER}

COPY . ${HOME}
USER root
RUN chown -R ${NB_UID} ${HOME}
USER ${NB_USER}

WORKDIR ${HOME}

RUN pip install --no-cache-dir -e .

WORKDIR ${HOME}
USER ${USER}
