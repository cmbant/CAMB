#Dockerfile for running pycamb notebooks with binder
#https://mybinder.org/v2/gh/cmbant/camb/master?filepath=pycamb%2Fdocs%2FCAMBdemo.ipynb

FROM cmbant/cosmobox:python3

ENV NB_USER jovyan
ENV NB_UID 1000
ENV HOME /home/${NB_USER}

RUN adduser --disabled-password \
    --gecos "Default user" \
    --uid ${NB_UID} \
    ${NB_USER}
    
COPY . ${HOME}
USER root
RUN chown -R ${NB_UID} ${HOME}
USER ${NB_USER}

WORKDIR ${HOME}/pycamb
RUN python setup.py build

WORKDIR ${HOME}

