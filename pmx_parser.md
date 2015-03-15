This module contains some functions that speed up file parsing

## Functions: ##

```
def kickOutComments( lines, comment='#'): # remove all comments and empty lines

def readSection( lines, beg, end ):  # read entries between two keywords
>>> l = readSection( lines, "[ atomtypes ]", "[")  # starts reading at "atomtypes" and stops at the first occurence of "["

def parseList( how, lines ):     # converts a list of lines into a defined format
>>> l = parseList("ssiff", lines )     # returns a list of [string, string, int, float, float] lists

def read_and_format(filename, format_string, comment = "#")  # commonly use combination of the above functions
>>> l = read_and_format( "datafile.dat", "ssiff")
```