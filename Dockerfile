FROM mambaorg/micromamba:1.4.2
USER root

RUN apt-get update && apt-get install gcc g++ bowtie2 samtools \
  -y --no-install-recommends \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/* \
  && rm -rf /usr/share/man/* \
  && rm -rf /usr/share/doc/* \
  && micromamba install -c conda-forge -n base -y ncurses \
  && micromamba install -c defaults -c conda-forge -c bioconda -y -n base --debug -c bioconda \
      trimmomatic flash numpy cython jinja2 scipy matplotlib pandas plotly \
      crispresso2 biopython cutadapt cas-offinder samtools=1.17 \
  && micromamba clean --all --yes

#install ms fonts
RUN echo "deb http://httpredir.debian.org/debian buster main contrib" > /etc/apt/sources.list \
  && echo "deb http://security.debian.org/ buster/updates main contrib" >> /etc/apt/sources.list \
  && echo "ttf-mscorefonts-installer msttcorefonts/accepted-mscorefonts-eula select true" | debconf-set-selections \
  && apt-get update \
  && apt-get install -y ttf-mscorefonts-installer \
  && apt-get clean \
  && apt-get autoremove -y \
  && rm -rf /var/lib/apt/lists/* \
  && rm -rf /usr/share/man/* \
  && rm -rf /usr/share/doc/* \
  && rm -rf /usr/share/zoneinfo

# install crispresso
COPY CRISPRlungo.py /CRISPRlungo/CRISPRlungo.py
COPY lib /CRISPRlungo/lib
WORKDIR /CRISPRlungo

ENTRYPOINT ["python","/CRISPRlungo/CRISPRlungo.py"]
