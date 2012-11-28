/* Header to test of C modules for arrays for Python: C_test.c */

typedef struct XDR_HEAD XDR_HEAD;

/* ==== Prototypes =================================== */
//init function for python
void initNGMX_XTC(void);

// .... Python callable Vector functions ..................
int skip_frame(XDRFILE* xfp, XDR_HEAD* header, long file_byte_len);
int read_headder(XDRFILE* xfp, XDR_HEAD* header, long file_byte_len);
static PyObject *xdrfile_NT(PyObject *self, PyObject *args);
static PyObject *xdrfile_read(PyObject *self, PyObject *args);
static PyObject *xdrfile_write(PyObject *self, PyObject *args);
//static PyObject *xdrfile_readI(PyObject *self, PyObject *args);

