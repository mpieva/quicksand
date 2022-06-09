FROM ubuntu:20.04

ENV PATH=$PATH:/opt/ghc/bin:/root/.cabal/bin/:/root/.cargo/bin:/usr/local/krakenuniq/:/usr/local/kraken/

COPY pkgs/ /home/root/

WORKDIR /home/root/

#get the correct c-compile-library to compile bwa
RUN echo "deb http://dk.archive.ubuntu.com/ubuntu xenial main" >> /etc/apt/sources.list && \
    echo "deb http://dk.archive.ubuntu.com/ubuntu xenial universe" >> /etc/apt/sources.list && \
    apt update && apt install -y g++-4.8 gcc-4.8 && \
    DEBIAN_FRONTEND=noninteractive apt install -y tzdata build-essential && \
    update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-4.8 1 && \
    update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-4.8 1

#install bwa
RUN apt-get install -y libzmq3-dev pkg-config zlib1g-dev git && \
    git clone -b 0.5.10-evan.10 --depth 1 https://github.com/mpieva/network-aware-bwa && \
    cd network-aware-bwa && \
    make && \
    mv bwa /bin/ && \
    cd .. && \
    rm -r network-aware-bwa

#install ghc
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y software-properties-common && \
    add-apt-repository -y ppa:hvr/ghc && \
    apt-get update && \
    apt-get install -y ghc-8.0.2 libjudy-dev wget

#install cabal and cabal-install (haskell library-installer)
RUN update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-9 2 && \
    update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-9 2 && \
    update-alternatives --set g++ /usr/bin/g++-9 && \
    update-alternatives --set gcc /usr/bin/gcc-9 && \
    wget https://downloads.haskell.org/~cabal/Cabal-1.24.2.0/Cabal-1.24.2.0.tar.gz && \
    tar xzf Cabal-1.24.2.0.tar.gz && \
    rm Cabal-1.24.2.0.tar.gz && \
    cd Cabal-1.24.2.0 && \
    ghc -threaded --make Setup && \
    ./Setup configure && \
    ./Setup build && \
    ./Setup install && \
    cd .. && \
    rm -r Cabal-1.24.2.0

RUN wget https://hackage.haskell.org/package/cabal-install-1.24.0.2/cabal-install-1.24.0.2.tar.gz && \
    tar xzf cabal-install-1.24.0.2.tar.gz && \
    rm cabal-install-1.24.0.2.tar.gz && \
    cd cabal-install-1.24.0.2 && \
    EXTRA_CONFIGURE_OPTS="" ./bootstrap.sh && \
    cabal update && \
    cd .. && \
    rm -r cabal-install-1.24.0.2

# install Biohazard
RUN cabal install --global base-prelude-1.2.0.1 && \
    git clone https://github.com/mpieva/biohazard && \
    cd biohazard && \
    git checkout 0.6.15 && \
    cabal install --global . && \
    cd .. && \
    rm -r biohazard

# and Biohazard-tools
RUN apt-get install -y libsnappy-dev && \
    git clone https://github.com/mpieva/biohazard-tools && \
    cd biohazard-tools && \
    git checkout 2231874 && \
    cabal install --global . && \
    cd .. && \
    rm -r biohazard-tools

# get the additional software
RUN apt-get install -y samtools bedtools python3-pip && \
    python3 -m pip install pysam && \
    python3 -m pip install pandas &&\
    python3 -m pip install scipy

# install Kraken
RUN git clone https://github.com/DerrickWood/kraken && \
    cd kraken && \
    ./install_kraken.sh /usr/local/kraken 

# install KrakenUniq
RUN git clone https://github.com/fbreitwieser/krakenuniq && \
    cd krakenuniq && \
    ./install_krakenuniq.sh /usr/local/krakenuniq 

# install jellyfish 1.1.12
RUN wget https://github.com/gmarcais/Jellyfish/releases/download/v1.1.12/jellyfish-1.1.12.tar.gz && \
    tar -zxf jellyfish-1.1.12.tar.gz && \
    rm jellyfish-1.1.12.tar.gz && \
    cd jellyfish-1.1.12 && \
    ./configure --prefix=/usr/local/ && \
    make -j 4 && \
    make install


# install splitbam
RUN tar xzf splitbam-0.1.6a4.tar.gz && \
    python3 -m pip install splitbam-0.1.6a4/ && \
    rm splitbam-0.1.6a4.tar.gz && \
    rm -r splitbam-0.1.6a4

#install the rust-packages
RUN apt-get install -y cargo && \
    apt-get install -y clang && \
    tar xzf bamfilter-0.2.9.crate && \
    cargo install --path bamfilter-0.2.9/ --root /usr/local/ && \
    rm bamfilter-0.2.9.crate && \
    rm -r bamfilter-0.2.9 && \
    tar xzf bam-lengthfilter-0.1.1.crate && \
    cargo install --path bam-lengthfilter-0.1.1 --root /usr/local/ && \
    rm bam-lengthfilter-0.1.1.crate && \
    rm -r bam-lengthfilter-0.1.1
    
#install packages required for the datastructure
RUN apt-get install -y rsync ncbi-blast+ 
RUN python3 -m pip install biopython

