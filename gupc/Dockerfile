
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
	flex \
        g++ \
        gcc \
	gcc-multilib \
        gfortran \
        git \
	libgmp3-dev \
	libmpc-dev \
	libmpfr-dev \
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

ENV GUPC_URL=https://www.gccupc.org/gupc-5201-1/32-gupc-5-2-0-1-source-release/file
#ENV GUPC_URL=https://www.gccupc.org/gupc-5201-1/30-gupc-5201-x8664-ubuntu-1204/file
ENV GUPC_ARCHIVE=${PREFIX}/src/gupc.tar.bz2
ENV GUPC_SRC_DIR=${PREFIX}/src/gupc/
ENV GUPC_BUILD_DIR=${PREFIX}/build/gupc/
ENV GUPC_INSTALL_DIR=${PREFIX}
RUN mkdir -p ${GUPC_SRC_DIR}
RUN mkdir -p ${GUPC_BUILD_DIR}
RUN mkdir -p ${GUPC_INSTALL_DIR}
RUN wget --quiet ${GUPC_URL} --output-document=${GUPC_ARCHIVE}
RUN tar -jxf ${GUPC_ARCHIVE} -C ${GUPC_SRC_DIR} --strip-components=1
WORKDIR ${GUPC_BUILD_DIR}
RUN ${GUPC_SRC_DIR}/configure --prefix=${GUPC_INSTALL_DIR} --enable-languages=c \
    --with-mpfr-include=/usr/include --with-mpfr-lib=/usr/lib/x86_64-linux-gnu \
    --with-gmp-include=/usr/include --with-gmp-lib=/usr/lib/x86_64-linux-gnu \
    --with-mpc-include=/usr/include --with-mpc-lib=/usr/lib/x86_64-linux-gnu\
    --disable-multilib --disable-bootstrap
RUN make -j ${NPROCS}
RUN make install
RUN rm -rf ${GUPC_ARCHIVE}

ENV PATH=${PREFIX}/bin:$PATH
RUN adduser --disabled-password --gecos "" bale_user
CMD ["/bin/bash"]
