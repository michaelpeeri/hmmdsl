#!/bin/bash

export BOOST_ROOT=/home/michael/src/boost_1_57_0
export BOOST_BUILD_PATH=/home/michael/src/boost_1_57_0

export PYTHONPATH=/usr/local/lib/python3.4/dist-packages:/usr/lib/python3/dist-packages:/usr/lib/python3.4/dist-packages:/home/michael/signal/hmmdsl11/bin/clang-linux-3.5/release/

export LD_LIBRARY_PATH=/home/michael/src/boost_1_57_0/bin.v2/libs/python/build/clang-linux-3.5/release/:/home/michael/src/boost_1_57_0/bin.v2/libs/random/build/clang-linux-3.5/release/threading-multi/:/home/michael/src/boost_1_57_0/bin.v2/libs/system/build/clang-linux-3.5/release/threading-multi/:/home/michael/signal/hmmdsl11/bin/clang-linux-3.5/release/

export PATH=$PATH:$BOOST_ROOT

echo Now run ldconfig!
