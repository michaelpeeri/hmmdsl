#!/bin/bash

export BOOST_ROOT=/home/michael/src/boost_1_60_0
export BOOST_BUILD_PATH=/home/michael/src/boost_1_60_0

export PYTHONPATH=/usr/local/lib/python3.4/dist-packages:/usr/lib/python3/dist-packages:/usr/lib/python3.4/dist-packages:/home/michael/signal/hmmdsl11/bin/clang-linux-3.6/debug/

export LD_LIBRARY_PATH=/home/michael/src/boost_1_60_0/bin.v2/libs/python/build/clang-linux-3.6/debug/:/home/michael/src/boost_1_60_0/bin.v2/libs/random/build/clang-linux-3.6/debug/threading-multi/:/home/michael/src/boost_1_60_0/bin.v2/libs/system/build/clang-linux-3.6/debug/threading-multi/:/home/michael/signal/hmmdsl11/bin/clang-linux-3.6/debug/


export PATH=$PATH:$BOOST_ROOT

alias python3=/home/michael/src/Python-3.4.1/python

alias clang++=clang++-3.6

echo Now run ldconfig!
