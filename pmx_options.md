This module contains classes that can be used to write programs with fancy command line options and program descriptions.
It is adapted from the Gromacs command line parsing functionality.

The commandline parsing system distinguishes two kinds of options, options that accept a value from the commandline and options that take one or more filenames as arguments.
The first kind of options is initialized by the "Option" class.
an Option instance is initialized by the commandline flag, the data type of the argument, the default value and a description of the option. Multiple options can be initialized and stored in a list of options as indicated below.

```
options = [
Option( "-s", "string", "default_string", "Description of what this string options is doing"),
Option("-f", "float", 3.1415,"Description of what this float option is doing"),
Option("-i","int",42,"Description of what this int option is doing"),
Option("-b","bool",True,"Description of what this bool option is doing"),
Option("-vec1","rvec", [2.5, 3.8, -1.1], "Description of what this real vector option is doing"),
Option("-vec2","ivec", [2, 3, 0], "Description of what this int vector option is doing"),
Option("-vec3","svec", ["a", "b", "c"], "Description of what this string vector option is doing")
]
```
The file options are similarly initialized
```
file_options = [
FileOption("-pdb", "r",["pdb"], "protein.pdb", input pdb file"),
FileOption("-opdb", "w",["pdb","gro"], "out.pdb", output pdb or gro file"),
FileOption("-mpdb", "r/m",["pdb","gro"], "one_of_many.pdb", several pdb files"),
  ]
```

the first argument is the commandline flag, the second indicates whether we would like to read or write this file and if we would like to read more than one file (r/m). The third argument determines the file type(s) and again the last argument is a short description.

In addition you can provide a more exhaustive description of your program, simply as a list of strings.
```
help_text = [ "This program does useful things",
                  "as long as you use option a,b and c",
                  "but not in combination with d and e"]
```

Having all options initialized you can call the commandline parser class

```
cmdl = Commandline( sys.argv, options = options, fileoptions = file_options, program_desc = help_text, version = "2.3")
```

The option "-h" is automatically added and triggers the display of the help message when executing the script.

```
$ python test_opt.py -h
HELP TEXT for "test_opt.py"
-----------------------------------------------------------------------------------------

This program does useful things
as long as you use option a,b and c
but not in combination with d and e
-----------------------------------------------------------------------------------------
 Program: test_opt.py ( v. 2.3)    | pymacs version 0.6.0
-----------------------------------------------------------------------------------------
     File Options          | Type(s)|Mode    | File(s)                | Description
-----------------------------------------------------------------------------------------
     -pdb                  | pdb|r           | protein.pdb            | input pdb file
     -opdb                 | pdb,gro|w       | out.pdb                | output pdb or gro file
     -mpdb                 | pdb,gro|r/m     | one_of_many.pdb        | several pdb files
-----------------------------------------------------------------------------------------
     Options               | Type            | Value                  | Description
-----------------------------------------------------------------------------------------
     -[no]h !              | bool            | True                   | Show help message and quit 
     -s                    | string          | default_string         | Description of what this 
                           |                 |                        | string options is doing
     -f                    | float           | 3.1415                 | Description of what this 
                           |                 |                        | float option is doing
     -i                    | int             | 42                     | Description of what this 
                           |                 |                        | int option is doing
     -[no]b                | bool            | True                   | Description of what this 
                           |                 |                        | bool option is doing
     -vec1                 | rvec            | [2.5, 3.0, -1.0]       | Description of what this 
                           |                 |                        | real vector option is doing
     -vec2                 | ivec            | [2, 3, 0]              | Description of what this 
                           |                 |                        | int vector option is doing
     -vec3                 | svec            | ['a', 'b', 'c']        | Description of what this 
                           |                 |                        | string vector option is doing
-----------------------------------------------------------------------------------------

```

Once this is executed you can access the parsed values by the commandline flag

```
string_opt = cmdl['-s']
input_pdb = cmdl['-pdb']
```