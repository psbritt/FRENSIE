FROM ubuntu:18.04

RUN useradd -u 8877 philip

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

RUN mkdir /home/philip/dependencies; \
    cd /home/philip/dependencies; \
    mkdir /home/philip/dependencies/cmake; \
    mkdir /home/philip/temp; \
    cd /home/philip/temp; \
    wget https://cmake.org/files/v3.17/cmake-3.17.1.tar.gz; \
    tar -xvf cmake-3.17.1.tar.gz; \
    cd cmake-3.17.1; \
    mkdir build; \
    cd build; \
    ../configure --prefix="/home/philip/dependencies/cmake"; \
    make -j16; \
    make install; \
    cd /home/philip/; \
    echo 'export PATH=/home/philip/dependencies/cmake/bin:$PATH' >> /home/philip/.bashrc 

ENV PATH /home/philip/dependencies/cmake/bin:$PATH 

RUN mkdir /home/philip/dependencies/hdf5; \
    cd /home/philip/temp; \
    wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.13/src/hdf5-1.8.13.tar.gz; \
    tar -xvf hdf5-1.8.13.tar.gz; \
    cd hdf5-1.8.13; \
    mkdir build; \
    cd build; \
    ../configure --prefix="/home/philip/dependencies/hdf5" --enable-optimized --enable-shared --enable-cxx --enable-hl --disable-debug; \
    make -j 16; \
    make install; \ 
    echo 'export PATH=/home/philip/dependencies/hdf5/bin:$PATH' >> /home/philip/.bashrc; \
    echo 'export LD_LIBRARY_PATH=/home/philip/dependencies/hdf5/lib:$LD_LIBRARY_PATH' >> /home/philip/.bashrc; \
    exec bash 

ENV PATH /home/philip/dependencies/hdf5/bin:$PATH 
ENV LD_LIBRARY_PATH /home/philip/dependencies/hdf5/lib:$LD_LIBRARY_PATH

RUN mkdir /home/philip/dependencies/swig; \
    cd /home/philip/temp; \
    wget -O swig-4.0.0.tar.gz https://sourceforge.net/projects/swig/files/swig/swig-4.0.0/swig-4.0.0.tar.gz/download; \
    tar -xvf swig-4.0.0.tar.gz; \
    cd swig-4.0.0; \
    mkdir build; \
    cd build; \
    ../configure --prefix="/home/philip/dependencies/swig"; \
    make -j 16; \
    make install; \
    echo 'export PATH=/home/philip/dependencies/swig/bin:$PATH' >> /home/philip/.bashrc; \
    exec bash

ENV PATH /home/philip/dependencies/swig/bin:$PATH

RUN mkdir /home/philip/dependencies/boost; \
    cd /home/philip/temp; \
    wget -O boost_1_72_0.tar.gz https://sourceforge.net/projects/boost/files/boost/1.72.0/boost_1_72_0.tar.gz/download; \
    tar -xvf boost_1_72_0.tar.gz; \
    cd boost_1_72_0; \
    ./bootstrap.sh --prefix=/home/philip/dependencies/boost; \
    sed -i "$ a using mpi ;" project-config.jam; \
    ./b2 -j16 --prefix=/home/philip/dependencies/boost -s NO_BZIP2=1 link=shared runtime-link=shared install; \
    echo 'export LD_LIBRARY_PATH=/home/philip/dependencies/boost/lib:$LD_LIBRARY_PATH' >> /home/philip/.bashrc; \
    exec bash

ENV LD_LIBRARY_PATH /home/philip/dependencies/boost/lib:$LD_LIBRARY_PATH

RUN mkdir /home/philip/dependencies/moab; \
    cd /home/philip/temp; \
    git clone --single-branch -b Version5.1.0 https://bitbucket.org/fathomteam/moab; \
    cd moab; \
    autoreconf -fi; \
    mkdir build; \
    cd build; \
    ../configure --enable-optimize --enable-shared --disable-debug --with-hdf5=/home/philip/dependencies/hdf5 --prefix=/home/philip/dependencies/moab/; \
    make -j 16; \
    make install; \
    echo 'export PATH=/home/philip/dependencies/moab/bin:$PATH' >> /home/philip/.bashrc; \
    echo 'export LD_LIBRARY_PATH=/home/philip/dependencies/moab/lib:$LD_LIBRARY_PATH' >> /home/philip/.bashrc; \
    exec bash

ENV PATH /home/philip/dependencies/moab/bin:$PATH
ENV LD_LIBRARY_PATH /home/philip/dependencies/moab/lib:$LD_LIBRARY_PATH

RUN mkdir /home/philip/dependencies/dagmc; \   
    cd /home/philip/temp; \
    git clone -b develop --single-branch https://github.com/svalinn/DAGMC.git; \
    cd DAGMC; \
    mkdir build; \
    cd build; \
    env HDF5_/home/philip=/home/philip/dependencies/hdf5; \ 
    cmake .. -DCMAKE_INSTALL_PREFIX=/home/philip/dependencies/dagmc -DCMAKE_BUILD_TYPE:STRING=Release -DMOAB_DIR=/home/philip/dependencies/moab; \
    make -j 16; \
    make install; \
    echo 'export PATH=/home/philip/dependencies/dagmc/bin:$PATH' >> /home/philip/.bashrc; \
    echo 'export LD_LIBRARY_PATH=/home/philip/dependencies/dagmc/lib:$LD_LIBRARY_PATH' >> /home/philip/.bashrc; \
    exec bash

ENV PATH /home/philip/dependencies/dagmc/bin:$PATH
ENV LD_LIBRARY_PATH /home/philip/dependencies/dagmc/lib:$LD_LIBRARY_PATH

COPY database.xml /home/philip/data/database.xml

COPY native /home/philip/data/native

RUN mkdir /home/philip/FRENSIE; \
    cd /home/philip/FRENSIE; \
    git clone -b HPC_docker --single-branch https://github.com/psbritt/FRENSIE.git; \
    ln -s FRENSIE src; \
    mkdir build; \
    cd build; \
    source ~/.bashrc; \
    cmake -D CMAKE_INSTALL_PREFIX:PATH=/home/philip/FRENSIE \
          -D CMAKE_BUILD_TYPE:STRING=DEBUG \
          -D CMAKE_VERBOSE_CONFIGURE:BOOL=OFF \
          -D CMAKE_VERBOSE_MAKEFILE:BOOL=ON \
          -D FRENSIE_ENABLE_DBC:BOOL=ON \
          -D FRENSIE_ENABLE_COLOR_OUTPUT:BOOL=OFF \
          -D FRENSIE_ENABLE_OPENMP:BOOL=ON \
          -D FRENSIE_ENABLE_MPI:BOOL=ON \
          -D FRENSIE_ENABLE_MOAB:BOOL=ON \
          -D FRENSIE_ENABLE_DAGMC:BOOL=ON \
          -D FRENSIE_ENABLE_/home/philip:BOOL=OFF \
          -D FRENSIE_ENABLE_PROFILING:BOOL=OFF \
          -D FRENSIE_ENABLE_COVERAGE:BOOL=OFF \
          -D FRENSIE_ENABLE_EXPLICIT_TEMPLATE_INST:BOOL=ON \
          -D MOAB_PREFIX:PATH=/home/philip/dependencies/moab \
          -D DAGMC_PREFIX:PATH=/home/philip/dependencies/dagmc \
          -D HDF5_PREFIX:PATH=/home/philip/dependencies/hdf5 \
          -D BOOST_PREFIX:PATH=/home/philip/dependencies/boost \
          -D SWIG_PREFIX:PATH=/home/philip/dependencies/swig \
          -D FRENSIE_ENABLE_DASHBOARD_CLIENT:BOOL=ON \
          -D XSDIR:PATH=/home/philip/data \
          $@ \
          /home/philip/FRENSIE/src; \
    make -j14; \
    make install

RUN rm -rf /home/philip/temp

RUN chown -R philip /home/philip