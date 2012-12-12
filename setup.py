# pmx  Copyright Notice
# ============================
#
# The pmx source code is copyrighted, but you can freely use and
# copy it as long as you don't change or remove any of the copyright
# notices.
#
# ----------------------------------------------------------------------
# pmx is Copyright (C) 2006-2013 by Daniel Seeliger
#
#                        All Rights Reserved
#
# Permission to use, copy, modify, distribute, and distribute modified
# versions of this software and its documentation for any purpose and
# without fee is hereby granted, provided that the above copyright
# notice appear in all copies and that both the copyright notice and
# this permission notice appear in supporting documentation, and that
# the name of Daniel Seeliger not be used in advertising or publicity
# pertaining to distribution of the software without specific, written
# prior permission.
#
# DANIEL SEELIGER DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS
# SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND
# FITNESS.  IN NO EVENT SHALL DANIEL SEELIGER BE LIABLE FOR ANY
# SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER
# RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF
# CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
# CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
# ----------------------------------------------------------------------

import os,sys
from distutils.core import setup, Extension
from distutils.command.install_data import install_data

class install_data_files(install_data):

    def run(self):
        install_cmd = self.get_finalized_command('install')
        self.install_dir = getattr(install_cmd, 'install_lib')
        return install_data.run(self)


pmx = Extension('pmx/_pmx',
                   libraries = ['m'],
                   include_dirs = ['src'],
                   sources = ['src/Geometry.c','src/wrap_Geometry.c',
                              'src/init.c', 'src/Energy.c']
                   
                 )

## cpp_test = Extension('pmx/_cpp_test',
##                    libraries = ['m'],
##                    include_dirs = ['src'],
##                    sources = ['src/cpp_test.cpp']
                   
##                  )



setup (name = 'pmx',
       version = '1.0.0',
       description = 'Python Toolbox structure file editing and writing simulation setup/analysis tools',
       author = 'Daniel Seeliger',
       author_email = 'seeliger.biosoft@gmail.de',
       url = 'http://code.google.com/p/pmx/',
       long_description = '''Tools to play with structure files, topologies, index files, xtc files, etc....''',
       packages = ['pmx'],
       data_files = [('pmx/data',['data/bbdep.pkl']),
                     ('pmx/data',['data/bp.pkl']),
                     #('pmx/data',['data/fragments.pkl']),
                     ('pmx/data',['data/ffamber99sb.rtp']),
                     ('pmx/data',['data/ffamber99sbbon.itp']),
                     ('pmx/data',['data/ffamber99sbnb.itp']),
                     ('pmx/data',['data/blosum62_new.mat'])
                     ],
       ext_modules = [pmx],
       cmdclass = {'install_data': install_data_files
                   },

       )

