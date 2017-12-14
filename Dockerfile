#Dockerfile for running pycamb notebooks with binder

FROM cmbant/cosmobox:gcc6

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

RUN cd pycamb; python setup.py install; cd ..
CMD jupyter-notebook --ip=* --no-browser
