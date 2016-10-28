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
__doc__="""
Some functions to write parsers for different data formats.

Usage:

>>> lst = open(infile,'r').readlines()     # read file
>>> lst = kickOutComments(lst,'#')       # remove comments with # or ;
>>> lst = parseList('ifs',lst)            # parse each line and return
[integer, float, string] triples


>>> subl = readSection(lst,'[ begin ]','[ end ]') # get all lines between
[ begin ] and [ end ]


"""

import sys
from collections import OrderedDict
import logging

logger = logging.getLogger()




def kickOutComments( lines, comment = '#'):
    ret = []
    for line in lines:
        if comment in line:
            idx = line.index(comment)
            new_line = line[:idx].strip()
            if new_line:
                ret.append( new_line )
        else:
            new_line = line.strip()
            if new_line:
                ret.append( new_line )
    return ret

def readSection(lines, begin, end):
    ret = []
    for i, line in enumerate(lines):
        if line.strip() == begin:
            for line2 in lines[i+1:]:
                if not end in line2:
                    ret.append( line2 )
                else:
                    return ret
    return ret

def __parse_error(msg, line):
	logger.error('Error while parsing line: %s' % line)
	logger.error(msg)
	sys.exit(1)

def __parse_entry(entr, tp):
    new = None
    if tp == 's':
        new = entr
    elif tp == 'i':
        try:
            new = int(entr) 
        except:
            __parse_error("Integer conversion failed", entr)
    elif tp == 'f':
        try:
            new = float(entr) 
        except:
            __parse_error("Float conversion failed", entr)
    return new
    

def parseList(format_string, lst, ignore_missing = False):

	""" Read a list of lines into a defined format.
	The format_string can contain i, f and s (for integer, float and string)
	e.g. 'ssifff' would expect 6 entries and these would be imported as string, string,
	integer, float, float. By default the function throws an error if data is missing. This can be
	switched off by setting ignore_missing to True.

	:param format_string: e.g. 'sffiiff'
	:param lst: list of lines
	:param ignore_missing: Ignore if a line does not match the format string
	:returns: formatted list of lists
	"""

    ret = []
    format_list = [ a for a in format_string ]
    for line in lst:
        entr = line.split()
        if len(entr) != len(format_list ):
            if not ignore_missing or len(entr) < len(format_list ):
                __parse_error("Cannot convert line into format: %s" % format_string, line)
        new_list = []
        for i, tp in enumerate( format_list ):
            new_list.append( __parse_entry( entr[i], tp ) )
        ret.append( new_list )
    return ret
            

def read_and_format(filename, format_string, comment = '#', ignore_missing = False):
	""" A shortcut function to read a file into a defined format.

	:param filename: filename
	:param format_string: same is in parseList, e.g. 'siffiff'
	:param comment: content after this character is ignored, e.g. the line '1 2 3 # XX' will be parsed as [1,2,3]
	:param ignore_missing:  Ignore if a line does not match format string 
	:returns: formatted list of lists

	"""

    l = open(filename).readlines()
    if comment is not None:
        l = kickOutComments(l, comment)
    n = parseList( format_string, l, ignore_missing)
    return n



#=================================================
# some file format parsers frequently needed

def read_fasta( fn ):
    l = open(fn).readlines()
    l = kickOutComments(l)
    dic = OrderedDict()
    name = None
    for line in l:
        if line.startswith('>'):
            if name is not None:
                dic[name] = new_s
            new_s = ''
            name = line[1:].strip()
        else:
            seq = line.strip().upper()
            new_s+=seq
    dic[name] = new_s
    return dic


def read_xvg( fn,  style='xy'):
    l = open(fn).readlines()
    l = kickOutComments(l,'@')
    l = kickOutComments(l,'#')
    l = kickOutComments(l,'&')
    res = parseList('ff', l)
    if style == 'list':
        return res
    else:
        x = map(lambda a: a[0],  res)        
        y = map(lambda a: a[1],  res)
        return x,  y



