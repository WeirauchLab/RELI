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

WORKDIR /reli
ENV PATH $PATH:/reli

COPY src src/
COPY Makefile .

# install prerequisite development packages
# build RELI from source
# Remove unnecessary packages pulled in as dependencies of g++
RUN apk update &&\
    apk add g++ make gsl gsl-dev bzip2 bash which libstdc++ libgcc &&\
    make &&\
    ln -s RELI reli &&\
    apk del g++