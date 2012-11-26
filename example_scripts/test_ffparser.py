import sys, os

from pmx.forcefield2 import *

top = ITPFile( sys.argv[1])
## for atom in top.atoms:
##     atom.atomtypeB = 'DUM'
##     atom.qB = 0
##     atom.mB = atom.m
top.write('new.itp')


## from pmx.ffparser import *
## ## from pmx.library import _aacids_dic
## rtp = RTPParser( 'ffamber99sb.rtp')
## print rtp['ALA']

## for name, entry in rtp:
##     print name, entry['atoms'][0]
## print 'XXX'
## for name, entry in rtp:
##     print name, entry['atoms'][0]
## ala = rtp['ALA']
## print _aacids_dic
## for aa in _aacids_dic.values():
##     try:
##         del rtp[aa]
##     except:
##         pass
## rtp.add_entry( "ala", ala )
## print "ala" in rtp
#rtp.write(sys.stdout)
#print rtp["ALA"]['bonds']

#nb = NBParser("ffamber99sbnb.itp")
#print nb

