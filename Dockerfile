FROM ubuntu:xenial
ENV DEBIAN_FRONTEND noninteractive
ENV BOOST_VERSION 1.63.0
ENV BOOST_BUILD_SRC="/opt/boost-src"
ENV BOOST_BUILD_ROOT="/usr/local/share/boost-build/"

# Boost build section is based on:
# https://hub.docker.com/r/teego/mingw-boost/

RUN apt-get update \
    && apt-get install -y software-properties-common \
    && apt-get update \
    && apt-get install -y \
       build-essential \
       openssl \
       autotools-dev \
       libicu-dev \
       libbz2-dev \
       libssl-dev \
       libreadline6-dev \
       git-core \
       python3.5 \
       python3.5-dev \
       libpython3.5-dev \ 
#       libboost1.58-all-dev \
       clang-3.6 \
       emacs-nox \
       wget \
    && rm -rf /var/lib/apt/lists/* \
    && ln -s /usr/bin/clang++-3.6 /usr/bin/clang++ \
    && mkdir -p $BOOST_BUILD_SRC \
    && cd $BOOST_BUILD_SRC \
    && wget --content-disposition "https://downloads.sourceforge.net/project/boost/boost/1.63.0/boost_1_63_0.tar.bz2?use_mirror=netix" \
    && tar -xjf boost_1_63_0.tar.bz2 \
    && rm boost_1_63_0.tar.bz2 \
    && cd boost_1_63_0 \
    && ./bootstrap.sh \
            --with-toolset=gcc \
            --without-icu \
            --with-python=/usr/bin/python3.5 \
    && ( \
       ./b2  \
            --with-python \
              toolset=clang \
              variant=debug \
              install \
              || /bin/true \
              ) \
    && cd tools/build/ \
    && ./bootstrap.sh \
    && ./b2 install \
    && cd ../.. \
    && git clone https://github.com/michaelpeeri/hmmdsl.git --branch v1 /opt/hmmdsl \
    && cd /opt/hmmdsl/ \
    && bjam hmmdsl_py toolset=clang variant=debug \
    && bjam hmmdsl_py toolset=clang variant=release
CMD ["bash"]
