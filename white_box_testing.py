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
from math import exp, log, isnan
from tempfile import NamedTemporaryFile
from textwrap import wrap
from functools import reduce
import subprocess 
import hmmdsl_py
import dotio
import training
import make_hsmm_validation_prototype_basic01


def main(argv=None):

    # Use temporary fasta file to transmit test sequence to impl.
    with NamedTemporaryFile() as seqsFile:
        # Obtain testing model
        #target_model = make_hsmm_validation_prototype_basic01.create_target_model()
        test_model = make_hsmm_validation_prototype_basic01.create_test_model()

        hmmdsl_py.print_model_statistics(test_model)

        # Use single test-sequence (for now) 
        sequence = "AAAAACCC"

        # Write problem to text file, and invoke reference (Haskell) impl.
        dotio.WriteHaskellTextFormat(test_model, sequence, "/tmp/text.dat")

        # Write training sequence(s)
        seqsFile.write( bytes(">test_{}\n".format(1), "ASCII") )
        seqsFile.write( bytes("\n".join(wrap(sequence)), "ASCII") )
        seqsFile.write( bytes("\n", "ASCII") )
        seqsFile.flush()

        # Run model on sequence(s), to calculate matrices from core algos
        elems = []
        em = hmmdsl_py.EM(test_model, seqsFile.name, elems)
        #f1 = lambda em: em.reestimate_b()
        #f2 = lambda em: em.reestimate_a()
        #f3 = lambda em: em.reestimate_scale()
        #f4 = lambda em: em.reestimate_shape()

        test_model.debug_print()

        # TODO - get matrices from impl.
        def readMatrix( seqid, rows, cols, getval ):
            mat = []
            for i in range(rows):
                row = list(map( lambda t: getval(seqid, i, t), range(cols) ))
                mat.append(row)
            return mat


        safelog = lambda x: log(x) if x>1e-50 else float('-inf')
        
        em.reestimate_b()
        coreMats = {}
        coreMats["forward"] = readMatrix( 0, 10, 9, lambda n,i,t: em.forward(n, t, i) )
        coreMats["forward_begin"] = readMatrix( 0, 10, 9,  lambda n,i,t: em.forward_begin(n, t, i) )
        coreMats["backward"] = readMatrix( 0, 10, 9, lambda n,i,t: em.backward(n, t, i) )
        coreMats["backward_begin"] = readMatrix( 0, 10, 9, lambda n,i,t: em.backward_begin(n, t, i) )
        coreMats["gamma"] = readMatrix( 0, 10, 9, lambda n,i,t: safelog(em.gamma(n, t, i)) )
        coreMats["xi"] = readMatrix( 0, 9, 9,  lambda n,i,j: safelog(em.sigma_t_xi(n, j, i)) )

        print(coreMats["forward"])

        # Write problem to text file, and invoke reference (Haskell) impl.

        # Invoke the reference impl.
        ret = subprocess.call(('./test_hsmm', '/tmp/text.dat'))
        if( ret != 0 ): raise Exception("Got return value {0}!".format(ret))
        
        # Read ref impl. output values
        #def convertInf(x):
        #    if( x == float('inf') ):
        #        return 
        hsMats = dotio.ReadHaskellTextFormat("results.dat" )

        assert("forward" in hsMats)
        assert("forward_begin" in hsMats)
        assert("backward" in hsMats)
        assert("backward_begin" in hsMats)
        assert("gamma" in hsMats)
        assert("xi" in hsMats)

        # Write the comparison report
        report = open("report.html", "w")
        report.write("<html>")
        def diff(mat1, mat2, combinator_func, title="matrix", colorfunc=lambda x: "rgb(0,0,0)", backcolorfunc=lambda x: "white" ):
            out = []
            out.append( "<div>{0}</div>\n".format(title) )
            out.append("<table>")
            for row1, row2 in zip(mat1, mat2):
                out.append("<tr>")
                for v1, v2 in zip(row1, row2):
                    val = combinator_func(v1,v2)
                    color = colorfunc(val)
                    back = backcolorfunc(val)
                    out.append("""<td><div style="color: {1}; background-color: {4}" title="{2:.3} // {3:.3}">{0:.5}</div></td>""".format(val[0], color, val[1], val[2], back))
                out.append("</tr>\n")
            out.append("</table>\n")
            return "".join(out)
            #return "<div>{0}</div>".format(title) + "<table>" + "".join(reduce( lambda x,y:x+y, map(
            #        lambda row1, row2: ["<tr>"] +
            #        list(map(lambda x1, x2: "<td>{0:.4}</td>".format(combinator_func(x1, x2)),
            #                                        row1, row2))
            #        + ["</tr>"]
            #        , mat1, mat2 ))) + "</table>\n"
        
        
        def redorgreen(x, red="rgb(225, 35, 42)", green="rgb(100, 230, 120)"):
            return red if abs(x[0])>1e-4 else green

        def redtowhite(x):
            err = max(0, min(100, abs(x[0])))
            #print("{0} -> {1} ->".format(x,err))
            if( err < 1e-4 ):
                return "white"

            mag = int(max(0, min(200, 20*log(err+1))))
            
            return "rgb({0},{1},{2})".format(255,255,255-mag)
        


        def difference(x,y):
            # Handle +numeric_limits<double>::max, +inf
            x = min( x, 1e+300 )
            y = min( y, 1e+300 )

            # Handle -numeric_limits<double>::max, -inf
            x = max( x, -1e+300 )
            y = max( y, -1e+300 )
            
            # Handle nan
            if( isnan(x) ):
                if( isnan(y) ):
                    return nan
                else:
                    return 1e+300
            
            return y-x
            

        # TODO - compare matrices
        for s in ["forward", "forward_begin", "backward", "backward_begin", "gamma", "xi"]:
            report.write(diff(coreMats[s], hsMats[s], lambda x,y: (difference(x,y),x,y), title="{0} [delta]".format(s), colorfunc=redorgreen, backcolorfunc=redtowhite))

        report.write("""<div style="position: relative;">""")

        report.write("""<div style="width: 50%; position: relative; top: 0px; left: 0%; ">""")
        for s in ["forward", "forward_begin", "backward", "backward_begin", "gamma", "xi"]:
            report.write(diff(coreMats[s], hsMats[s], lambda x,y: (y,x,y), title="[ref] {0}".format(s)))
        report.write("""</div>""")

        report.write("""<div style="width: 50%; position: absolute; top: 0px; left: 50%; ">""")
        for s in ["forward", "forward_begin", "backward", "backward_begin", "gamma", "xi"]:
            report.write(diff(coreMats[s], hsMats[s], lambda x,y: (x,x,y), title="[core] {0}".format(s)))
        report.write("""</div>""")

        report.write("""</div>""")

        #report.write(diff(coreMats["forward"], hsMats["forward"], title="forward"))
        #report.write(diff(coreMats["forward_begin"], hsMats["forward_begin"], title="forward-begin"))
        #report.write(diff(coreMats["backward"], hsMats["backward"], title="backward"))
        #report.write(diff(coreMats["backward_begin"], hsMats["backward_begin"], title="backward-begin"))
        #report.write(diff(coreMats["gamma"], hsMats["gamma"], title="gamma"))
        #report.write(diff(coreMats["xi"], hsMats["xi"], title="xi"))
        report.write("</html>")
        report.close()
        
        
        


if __name__ == "__main__":
    sys.exit(main())
