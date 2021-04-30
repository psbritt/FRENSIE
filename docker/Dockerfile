FROM ubuntu:18.04

RUN apt-get update; \
    apt-get install -y \
        git \
        wget \
        doxygen \
        libpcre3 libpcre3-dev \
        gfortran \
        gcc \
        libssl-dev \
        python-pip \
        python-dev \
        libeigen3-dev \
        autoconf; \
    pip install numpy

RUN mkdir ~/temp; \
    cd ~/temp;\
    wget https://cmake.org/files/v3.17/cmake-3.17.1.tar.gz;\
    tar -xvf cmake-3.17.1.tar.gz;\
    cd cmake-3.17.1;\
    mkdir build ;\
    cd build ;\
    ../configure --prefix="/usr/local";\
    make -j16;\
    make install
