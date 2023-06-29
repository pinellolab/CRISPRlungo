FROM mambaorg/micromamba:1.4.2
USER root

ARG MAMBA_DOCKERFILE_ACTIVATE=1

RUN apt-get update && apt-get install gcc g++ bowtie2 samtools \
  -y --no-install-recommends \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/* \
  && rm -rf /usr/share/man/* \
  && rm -rf /usr/share/doc/* 

RUN micromamba install -c conda-forge -c bioconda -n base -y ncurses python=3.10.11 trimmomatic flash numpy scipy cython\
  && micromamba clean --all --yes

RUN micromamba install -c defaults -c conda-forge -c bioconda -y -n base \
      jinja2 matplotlib pandas crispresso2=2.2.12 \
  && micromamba clean --all --yes


RUN micromamba install -c conda-forge -c bioconda -y -n base --debug -c bioconda \
     biopython cas-offinder plotly samtools=1.17 \
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

# install crisplungo
COPY . /CRISPRlungo/
WORKDIR /CRISPRlungo
RUN pip install .

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "python","/CRISPRlungo/src/CRISPRlungo/cli.py"]
