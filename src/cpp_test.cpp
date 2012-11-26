// pmx  Copyright Notice
// ============================
//
// The pmx source code is copyrighted, but you can freely use and
// copy it as long as you don't change or remove any of the copyright
// notices.
//
// ----------------------------------------------------------------------
// pmx is Copyright (C) 2006-2011 by Daniel Seeliger
//
//                        All Rights Reserved
//
// Permission to use, copy, modify, distribute, and distribute modified
// versions of this software and its documentation for any purpose and
// without fee is hereby granted, provided that the above copyright
// notice appear in all copies and that both the copyright notice and
// this permission notice appear in supporting documentation, and that
// the name of Daniel Seeliger not be used in advertising or publicity
// pertaining to distribution of the software without specific, written
// prior permission.
//
// DANIEL SEELIGER DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS
// SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND
// FITNESS.  IN NO EVENT SHALL DANIEL SEELIGER BE LIABLE FOR ANY
// SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER
// RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF
// CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
// CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
// ----------------------------------------------------------------------
#include <Python.h>
#include <string>

using namespace std;

string string_add(const string& s1, const string& s2)
{
  return s1+s2;
}

PyObject *wrap_string_add( PyObject *self, PyObject *args )
{
  char *s1;
  char *s2;
  
  if(!PyArg_ParseTuple(args,"ss",&s1,&s2))
    return NULL;
  string ss1(s1);
  string ss2(s2);
  string s = string_add(ss1,ss2);
  return Py_BuildValue("s", s.c_str());
}

extern "C" {
  
  static PyMethodDef testMethods[] = {
    {"string_add", wrap_string_add, METH_VARARGS },
    {NULL, NULL}
  };
  
  void init_cpp_test()
  {    
    (void) Py_InitModule3("_cpp_test", testMethods, NULL);
  }
};



