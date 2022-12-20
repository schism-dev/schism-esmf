# SPDX-FileCopyrightText: 2022 Helmholtz-Zentrum Hereon
# SPDX-License-Identifier: CC0-1.0
# SPDX-FileContributor: Carsten Lemmen <carsten.lemmen@hereon.de

#FROM platipodium/esmf:openmpi-v8.4.0
FROM registry.hzdr.de/schism/esmf-docker/esmf:v8.4.0-openmpi

LABEL description="SCHISM-ESMF Docker environment based on Ubuntu"
LABEL author="Carsten Lemmen <carsten.lemmen@hereon.de>"
LABEL license="CC0-1.0"
LABEL copyright="2022 Helmholtz-Zentrum Hereon"

# Arguments can be passed via the --build-arg key=value command to the 
# docker build command.  The default values are set below to superbee 
ARG TVD_LIM="SB"
ARG COMMUNICATOR="openmpi"

# Make ARG variables available within the container
ENV TVD_LIM ${TVD_LIM}
ENV COMMUNICATOR ${COMMUNICATOR}

RUN apt-get update && apt-get -qy upgrade
RUN apt-get update && apt-get -qy \
    install make lib${COMMUNICATOR}-dev \
    cmake wget python3 python3-pip \
    python-is-python3 libnetcdf-dev \
    libnetcdff-dev libxerces-c-dev liblapack-dev libyaml-cpp-dev \
    libparmetis-dev libmetis-dev subversion cvs git \
    gcc-11 gfortran-11 g++-11

# Remove all mpich related packages if communicator is not mpich, or
# do the same for openmpi
RUN if [ "x${COMMUNICATOR}" != "xmpich" ] ; then apt-get remove -qy *mpich* ; fi
RUN if [ "x${COMMUNICATOR}" != "xopenmpi" ] ; then apt-get remove -qy *openmpi* ; fi

ENV PATH="/usr/lib64/${COMMUNICATOR}/bin:${PATH}"

WORKDIR /usr/src

ENV SCHISM_DIR="/usr/src/schism"
ENV SCHISM_BUILD_DIR="/usr/src/schism/build"
ENV SCHISM_ESMF_DIR="/usr/src/schism-esmf"

ENV OLDIO="ON"

RUN git clone  --branch master --depth 1 https://github.com/schism-dev/schism ${SCHISM_DIR}
RUN mkdir -p ${SCHISM_BUILD_DIR}
RUN cmake -S ${SCHISM_DIR}/src -B ${SCHISM_BUILD_DIR} -DOLDIO=${OLDIO} -DTVD_LIM=${TVD_LIM}
RUN make -C /usr/src/schism/build -j8 pschism

RUN git clone  --branch master --depth 1 https://github.com/schism-dev/schism-esmf ${SCHISM_ESMF_DIR}
RUN make -C ${SCHISM_ESMF_DIR} -j8 all

