
FROM cmbant/docker-gcc-build

MAINTAINER "Antony Lewis" <antony@cosmologist.info>, "Andre Vehreschild" <vehre@gmx.de>

## Create an unprivileged user.
RUN useradd camb \
  && mkdir /home/camb \
  && chown camb:camb /home/camb

RUN apt-get update \ 
  && apt-get install -y --no-install-recommends \
    python \
  && rm -rf /var/lib/apt/lists/*

## Configure default locale
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen \
  && locale-gen en_US.utf8 \
  && /usr/sbin/update-locale LANG=en_US.UTF-8

ENV LC_ALL en_US.UTF-8

### Run as user camb
USER camb

WORKDIR /home/camb

### Get forutils and CAMB from github
CMD git clone https://github.com/cmbant/forutils.git forutils \
  && git clone -b vehre_formatted https://github.com/cmbant/CAMB.git CAMB \
  && git clone https://github.com/vehre/CAMB_expected.git CAMB/Expected_Results \
  && cd forutils && make Release \
  && cd /home/camb/CAMB && make Release \
  && python python/CAMB_test_files.py testfiles --make_ini \
  && cd testfiles \
  && ln -s ../Expected_Results/test_outputs expected_outputs \
  && cd .. \
  && python python/CAMB_test_files.py testfiles --diff_to expected_outputs --verbose

   

