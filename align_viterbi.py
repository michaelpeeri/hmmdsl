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
from math import exp
import hmmdsl_py
import dotio

def main_test(model_file, fasta_file, args = None):
    model = dotio.load(model_file)

    #model.debug_print()
    if isinstance(model, hmmdsl_py.Model):
        hmmdsl_py.align_viterbi(fasta_file, model, [], 0 )
    else:
        hmmdsl_py.align_viterbi_hmm(fasta_file, model, [], 0 )


class UsageError(Exception):
    def __init__(self, msg):
        self.msg = msg

    def __init__(self):
        self.msg = "Usage: score_sequences <model.pickle> <sequences.fa>";


def main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "h", ["help"])
        except getopt.GetoptError as msg:
            raise UsageError(msg)

        if( len(argv) < 3 ):
            raise UsageError();

        main_test(argv[1], argv[2], argv)

    except UsageError as err:
        print( err.msg )
        return 2;



if __name__ == "__main__":
    sys.exit(main())

