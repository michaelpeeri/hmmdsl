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
import hmmdsl_py
import dotio


#model_file = "/groups/pupko/mich1/signal/hmmdsl/effectors.283.clean/perturb1001.283.clean.train_vs_all.pickle"
#model_file = "/groups/pupko/mich1/signal/hmmdsl/effectors.283.clean/perturb1001.283.clean.train_vs_all.pickle"
#model_file = "/groups/pupko/mich1/signal/hmmdsl/params_hsmm22_35c/perturb356.pickle_train_vs_effectors.faa.pickle"
#model_file = "/groups/pupko/mich1/signal/hmmdsl/params_hsmm14_35c/perturb287.pickle_train_vs_effectors.faa.pickle"
#model_file = "/groups/pupko/mich1/signal/hmmdsl/params_hsmm18_35c/perturb53.pickle_train_vs_effectors.faa.pickle"

def emit_sequences_impl(model_file, outfile, print_progress):
    model = dotio.load(model_file)
    e = hmmdsl_py.Emitter(model)

    num_seqs = 1000
    for i in range(num_seqs):
        while(True):
            seq = e.emit()
            if( True or len(seq) == 35  ):
                print( ">test_{}".format(i), file=outfile)
                print( seq, file=outfile )
                break

        if( print_progress ):
            if( (i%10000==0) and (i > 0) and (i < num_seqs) ):
                print("{:2.0f}% done".format(float(i)/num_seqs*100))

def emit_sequences(model_file, outfile=None):
    if( outfile == None ):
        emit_sequences_impl(model_file, sys.stdout, False )
    else:
        emit_sequences_impl(model_file, open(outfile, "w"), True )



class UsageError(Exception):
    def __init__(self, msg):
        self.msg = msg

def main(argv=None):
    if argv is None:
        argv = sys.argv
        try:
            try:
                opts, args = getopt.getopt(argv[1:], "h", ["help"])
            except getopt.GetoptError as msg:
                raise UsageError(msg)
            
            if( len(argv) < 2 ):
                raise UsageError("Model file not specified")

        except UsageError as err:
            print( err.msg )
            return 2;

        if( len(argv) >= 3 ):
            emit_sequences(argv[1], argv[2])
        else:
            emit_sequences(argv[1])


if __name__ == "__main__":
    sys.exit(main())
