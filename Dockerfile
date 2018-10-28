FROM perl:latest

# install dependencies
RUN apt-get update && apt-get update && apt-get install -y \
  rsync        \
  zlib1g-dev   \
  libbz2-dev   \
  liblzma-dev

# install htslib, which is required for bcftools
WORKDIR /app
RUN git clone https://github.com/samtools/htslib.git
WORKDIR /app/htslib
RUN autoheader  && \
  autoconf      && \
  ./configure   && \
  make          && \
  make install

# install bcftools
WORKDIR /app
RUN git clone https://github.com/samtools/bcftools.git
WORKDIR /app/bcftools
RUN make install

# install GenPro and dependencies
WORKDIR /app
RUN curl -L http://cpanmin.us | perl - App::cpanminus
RUN cachebuster=1 git clone https://github.com/wingolab-org/GenPro.git
WORKDIR /app/GenPro
RUN cpanm GenPro.tar.gz
