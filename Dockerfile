# SPDX-FileCopyrightText: 2022 Helmholtz-Zentrum Hereon
# SPDX-License-Identifier: CC0-1.0
# SPDX-FileContributor: Carsten Lemmen <carsten.lemmen@hereon.de

#FROM phusion/baseimage:jammy-1.0.1
FROM platipodium/esmf:mpich-v8.4.0

LABEL description="SCHISM-ESMF Docker environment based on Ubuntu"
LABEL author="Carsten Lemmen <carsten.lemmen@hereon.de>"
LABEL license="CC0-1.0"
LABEL copyright="2022 Helmholtz-Zentrum Hereon"

# Arguments can be passed via the --build-arg key=value command to the 
# docker build command.  The default values are set below to superbee 
ARG TVD_LIM="SB"
ARG COMMUNICATOR="mpich"

#RUN apt update && apt -qy install cmake gcc-11 python3 \
#    python-is-python3 lib${COMMUNICATOR}-dev libmetis-dev libnetcdf-dev \
#    libnetcdff-dev libparmetis-dev git
RUN apt update

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

