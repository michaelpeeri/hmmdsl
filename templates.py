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
import re


_symbol = re.compile("<<(\S+)>>")

class Evaluator:
    def __init__(self, globals, locals):
        self._globals = globals
        self._locals = locals
    
    def __call__(self, match):
        expression = match.group(1)

        try:
            return str(eval(expression, self._globals, self._locals))
        except Exception as err:
            print("Error in expression <<{}>>".format(expression))
            raise err


def generate(filein, fileout, globals, locals):
    fin = open(filein, "r")
    fout = open(fileout, "w")

    repl = Evaluator(globals, locals)
        
    for line in fin.readlines():
        line = _symbol.sub(repl, line)
        fout.write(line)



