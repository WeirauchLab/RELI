# ================================================================
#
#   Dockerfile for RELI (public release)
#
#   Author:   Kevin Ernst <kevin.ernst -at- cchmc.org>
#   Date:     25 July 2018
#
# ================================================================

# start with https://store.docker.com/images/centos
FROM centos:7

# install prerequisite development packages
RUN yum update -y
RUN yum install -y gcc gcc-c++ make gsl gsl-devel bzip2 which

RUN mkdir -p /reli/src /reli/example
WORKDIR /reli
COPY src src/
COPY example example/
COPY Makefile .

# build RELI from source
RUN make
RUN ln -s RELI reli
ENV PATH $PATH:/reli
