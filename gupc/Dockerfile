
FROM ubuntu:18.04

MAINTAINER Jason DeVinney <jgdevin@super.org>

# Build-arg: NPROCS for parallel make
ARG NPROCS=8

ENV PREFIX=/usr/local
RUN mkdir -p "${PREFIX}"
RUN mkdir -p "${PREFIX}/src"

# This will prevent questions from being asked during the install
ENV DEBIAN_FRONTEND noninteractive

# Get certs so HTTPS requests (e.g., git clone) works
# do we need to do this?
RUN apt-get -y update \
    && apt-get install -y --no-install-recommends \
        ca-certificates \
    && apt-get clean

# Install various useful/prereq packages 
RUN apt-get -y update \
    && apt-get install -y --no-install-recommends \
        autoconf \
        automake \
        g++ \
        gcc \
        gfortran \
        git \
	libgmp10 \
	libmpc3 \
	libmpfr6 \
	libnuma-dev \
        libtool \
        m4 \
        make \
	openssh-client \
        patch \
	pkg-config \
        wget \
        zlib1g-dev \
    && apt-get clean

RUN export GUPC_URL=https://www.gccupc.org/gupc-5201-1/30-gupc-5201-x8664-ubuntu-1204/file
RUN export GUPC_ARCHIVE=${PREFIX}/src/gupc.tar.gz
RUN export GUPC_SRC_DIR=${PREFIX}/src/gupc/
RUN export GUPC_INSTALL_DIR=${PREFIX}
RUN wget --quiet ${GUPC_URL} --output-document=${GUPC_ARCHIVE}
RUN tar -xzf ${GUPC_ARCHIVE} -C ${GUPC_SRC_DIR} --strip-components=1
WORKDIR ${GUPC_SRC_DIR}
RUN ./configure --prefix=${GUPC_INSTALL_DIR} --enable-languages=c,c++
RUN make -j ${NPROCS}
RUN make install
RUN rm -rf ${GUPC_ARCHIVE}

ENV PATH=${PREFIX}/bin:$PATH
RUN adduser --disable-password --gecos "" bale_user
CMD ["/bin/bash"]