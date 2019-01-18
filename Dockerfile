# ================================================================
#
#   Dockerfile for RELI (public release)
#   Version:  0.0.2
#   Author:   Kevin Ernst <kevin.ernst -at- cchmc.org>
#   Date:     25 July 2018
#
#   Build:    docker build --rm -t weirauchlab/reli:v0.0.2 -f Dockerfile .
#   Pull:     docker pull weirauchlab/reli:v0.0.2
#   Run:      docker run --rm -ti weirauchlab/reli:v0.0.2
#
# ================================================================

# start with https://store.docker.com/images/centos
FROM alpine

# install prerequisite development packages
RUN apk update
RUN apk add g++ make gsl gsl-dev bzip2 bash which

RUN mkdir -p /reli/src
WORKDIR /reli
COPY src src/
COPY Makefile .

# build RELI from source
RUN make
RUN ln -s RELI reli
ENV PATH $PATH:/reli

# Remove unnecessary packages pulled in as dependencies of g++
RUN apk del g++
RUN apk add libstdc++ libgcc
