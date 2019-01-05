# ================================================================
#
#   Dockerfile for RELI (public release)
#
#   Author:   Kevin Ernst <kevin.ernst -at- cchmc.org>
#   Date:     25 July 2018
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
