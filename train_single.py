#-------------------------------------------------------------------------------------
#  Copyright 2014 Michael Peeri
#
#  This file is part of hmmdsl.
#  hmmdsl is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  hmmdsl is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with hmmdsl.  If not, see <http://www.gnu.org/licenses/>.
#-------------------------------------------------------------------------------------
import sys
import getopt
import collections # for deque
from math import exp
import hmmdsl_py
import dotio
import training


# Simple nullary function returning RNG-generated numbers, seeded using /dev/random
class _rng:
    def __init__(self):
        from random import random, Random, SystemRandom
        init = SystemRandom()
        self.rng = Random(init.random())
    def __call__(self):
        return self.rng.random()

#a = UtilityMaximiser(3)
#a.SetUtility(0, 5)
#a.SetUtility(1, 2)
#a.SetUtility(2, 3)

def main_test(model_file, fasta_file, args = None):

    print( " *** Training model ***" )
    model = dotio.load( model_file )
    print(type(model))

    training.train(fasta_file, model, [], 0.05)
        
    #model.debug_print()
    dotio.WriteDOT(model, "{}_train_vs_effectors.faa.dot".format(model_file) )
    dotio.store(model, "{}_train_vs_effectors.faa.pickle".format(model_file) )


class UsageError(Exception):
    def __init__(self, msg):
        self.msg = msg


def main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "h", ["help"])
        except getopt.GetoptError as err:
            raise UsageError(err)

        if( len(argv) < 2 ):
            raise UsageError("No model specified");

        if( len(argv) < 3 ):
            raise UsageError("No sequences specified");

        main_test(argv[1], argv[2], argv)

    except UsageError as err:
        print( err.msg )
        return 2;



if __name__ == "__main__":
    sys.exit(main())

