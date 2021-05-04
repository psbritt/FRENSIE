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
        libblas-dev \
        liblapack-dev \
        python-pip \
        python-dev \
        libeigen3-dev \
        autoconf \
        openmpi-bin \
        libopenmpi-dev; \
    pip install numpy

RUN mkdir /root/dependencies; \
    cd /root/dependencies; \
    mkdir /root/dependencies/cmake; \
    mkdir /root/temp; \
    cd /root/temp; \
    wget https://cmake.org/files/v3.17/cmake-3.17.1.tar.gz; \
    tar -xvf cmake-3.17.1.tar.gz; \
    cd cmake-3.17.1; \
    mkdir build; \
    cd build; \
    ../configure --prefix="/root/dependencies/cmake"; \
    make -j16; \
    make install; \
    cd /root/; \
    echo 'export PATH=/root/dependencies/cmake/bin:$PATH' >> /root/.bashrc 

ENV PATH /root/dependencies/cmake/bin:$PATH 

RUN mkdir /root/dependencies/hdf5; \
    cd /root/temp; \
    wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.13/src/hdf5-1.8.13.tar.gz; \
    tar -xvf hdf5-1.8.13.tar.gz; \
    cd hdf5-1.8.13; \
    mkdir build; \
    cd build; \
    ../configure --prefix="/root/dependencies/hdf5" --enable-optimized --enable-shared --enable-cxx --enable-hl --disable-debug; \
    make -j 16; \
    make install; \ 
    echo 'export PATH=/root/dependencies/hdf5/bin:$PATH' >> /root/.bashrc; \
    echo 'export LD_LIBRARY_PATH=/root/dependencies/hdf5/lib:$LD_LIBRARY_PATH' >> /root/.bashrc; \
    exec bash 

ENV PATH /root/dependencies/hdf5/bin:$PATH 
ENV LD_LIBRARY_PATH /root/dependencies/hdf5/lib:$LD_LIBRARY_PATH

RUN mkdir /root/dependencies/swig; \
    cd /root/temp; \
    wget -O swig-4.0.0.tar.gz https://sourceforge.net/projects/swig/files/swig/swig-4.0.0/swig-4.0.0.tar.gz/download; \
    tar -xvf swig-4.0.0.tar.gz; \
    cd swig-4.0.0; \
    mkdir build; \
    cd build; \
    ../configure --prefix="/root/dependencies/swig"; \
    make -j 16; \
    make install; \
    echo 'export PATH=/root/dependencies/swig/bin:$PATH' >> /root/.bashrc; \
    exec bash

ENV PATH /root/dependencies/swig/bin:$PATH

RUN mkdir /root/dependencies/boost; \
    cd /root/temp; \
    wget -O boost_1_72_0.tar.gz https://sourceforge.net/projects/boost/files/boost/1.72.0/boost_1_72_0.tar.gz/download; \
    tar -xvf boost_1_72_0.tar.gz; \
    cd boost_1_72_0; \
    ./bootstrap.sh --prefix=/root/dependencies/boost; \
    sed -i "$ a using mpi ;" project-config.jam; \
    ./b2 -j16 --prefix=/root/dependencies/boost -s NO_BZIP2=1 link=shared runtime-link=shared install; \
    echo 'export LD_LIBRARY_PATH=/root/dependencies/boost/lib:$LD_LIBRARY_PATH' >> /root/.bashrc; \
    exec bash

ENV LD_LIBRARY_PATH /root/dependencies/boost/lib:$LD_LIBRARY_PATH

RUN mkdir /root/dependencies/moab; \
    cd /root/temp; \
    git clone --single-branch -b Version5.1.0 https://bitbucket.org/fathomteam/moab; \
    cd moab; \
    autoreconf -fi; \
    mkdir build; \
    cd build; \
    ../configure --enable-optimize --enable-shared --disable-debug --with-hdf5=/root/dependencies/hdf5 --prefix=/root/dependencies/moab/; \
    make -j 16; \
    make install; \
    echo 'export PATH=/root/dependencies/moab/bin:$PATH' >> /root/.bashrc; \
    echo 'export LD_LIBRARY_PATH=/root/dependencies/moab/lib:$LD_LIBRARY_PATH' >> /root/.bashrc; \
    exec bash

ENV PATH root/dependencies/moab/bin:$PATH
ENV LD_LIBRARY_PATH /root/dependencies/moab/lib:$LD_LIBRARY_PATH

RUN mkdir /root/dependencies/dagmc; \   
    cd /root/temp; \
    git clone -b develop --single-branch https://github.com/svalinn/DAGMC.git; \
    cd DAGMC; \
    mkdir build; \
    cd build; \
    env HDF5_ROOT=/root/dependencies/hdf5; \ 
    cmake .. -DCMAKE_INSTALL_PREFIX=/root/dependencies/dagmc -DCMAKE_BUILD_TYPE:STRING=Release -DMOAB_DIR=/root/dependencies/moab; \
    make -j 16; \
    make install; \
    echo 'export PATH=/root/dependencies/dagmc/bin:$PATH' >> /root/.bashrc; \
    echo 'export LD_LIBRARY_PATH=/root/dependencies/dagmc/lib:$LD_LIBRARY_PATH' >> /root/.bashrc; \
    exec bash

ENV PATH /root/dependencies/dagmc/bin:$PATH
ENV LD_LIBRARY_PATH /root/dependencies/dagmc/lib:$LD_LIBRARY_PATH

# Remove all the extra dependency source code and build files
#RUN cd ~; \
#    rm -rf temp

# Create directory that will hold cross section data
RUN cd ~; \
    mkdir data; \
    cd data; \
    mkdir native

COPY database.xml /root/data/database.xml

COPY native /root/data/native

RUN mkdir /root/FRENSIE; \
    cd /root/FRENSIE; \
    git clone -b HPC_docker --single-branch https://github.com/psbritt/FRENSIE.git; \
    ln -s FRENSIE src; \
    mkdir build; \
    cd build; \
    source ~/.bashrc; \
    cmake -D CMAKE_INSTALL_PREFIX:PATH=/root/FRENSIE \
          -D CMAKE_BUILD_TYPE:STRING=DEBUG \
          -D CMAKE_VERBOSE_CONFIGURE:BOOL=OFF \
          -D CMAKE_VERBOSE_MAKEFILE:BOOL=ON \
          -D FRENSIE_ENABLE_DBC:BOOL=ON \
          -D FRENSIE_ENABLE_COLOR_OUTPUT:BOOL=OFF \
          -D FRENSIE_ENABLE_OPENMP:BOOL=ON \
          -D FRENSIE_ENABLE_MPI:BOOL=ON \
          -D FRENSIE_ENABLE_MOAB:BOOL=ON \
          -D FRENSIE_ENABLE_DAGMC:BOOL=ON \
          -D FRENSIE_ENABLE_ROOT:BOOL=OFF \
          -D FRENSIE_ENABLE_PROFILING:BOOL=OFF \
          -D FRENSIE_ENABLE_COVERAGE:BOOL=OFF \
          -D FRENSIE_ENABLE_EXPLICIT_TEMPLATE_INST:BOOL=ON \
          -D MOAB_PREFIX:PATH=/root/dependencies/moab \
          -D DAGMC_PREFIX:PATH=/root/dependencies/dagmc \
          -D HDF5_PREFIX:PATH=/root/dependencies/hdf5 \
          -D BOOST_PREFIX:PATH=/root/dependencies/boost \
          -D SWIG_PREFIX:PATH=/root/dependencies/swig \
          -D FRENSIE_ENABLE_DASHBOARD_CLIENT:BOOL=ON \
          -D XSDIR:PATH=/root/data \
          $@ \
          /root/FRENSIE/src; \
    make -j14; \
    make install