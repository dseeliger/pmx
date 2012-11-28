#include "Python.h"
#include "arrayobject.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "NGMX_geometry.h"
#include "lap.h"			//for linear assignment problem


/* ==== Set up the methods table ====================== */
static PyMethodDef NGMX_geometryMethods[] = {
	{"test_fct", 		test_fct, 		METH_VARARGS},
	{"linAssignment", 	linAssignment, 	METH_VARARGS},
   // {"fit_rotmatrix", 	fit_rotmatrix, 	METH_VARARGS},
    {NULL, NULL}     /* Sentinel - marks the end of this structure */
};

/* ==== Initialize the C_test functions ====================== */
// Module name must be _C_arraytest in compile and linked
void initNGMX_geometry()  {
    (void) Py_InitModule("NGMX_geometry", NGMX_geometryMethods);
    import_array();  // Must be present for NumPy.  Called first after above line.
}

static PyObject *test_fct(PyObject *self, PyObject *args)
{
	return Py_BuildValue("i", 1);
}

static PyObject *linAssignment(PyObject *self, PyObject *args)
{
	PyArrayObject *pcost, *rowsol;
	double **ccost;
	int i, dim, *colsol;
	double *u, *v, lapcost;
	npy_intp pdim[0];

	if (!PyArg_ParseTuple(args, "O!:linAssignment",
		&PyArray_Type, &pcost
        ))  return Py_BuildValue("i", 0);


	dim     = PyArray_DIM(pcost, 0);
	pdim[0] = dim;
	rowsol  = PyArray_ZEROS(1, pdim, PyArray_INT, 0);  // vector that is to be given back to Python with the result
	ccost   = calloc(dim, sizeof *ccost);			   // converting python cost matrix --> **cost (pointer on pointer)
	colsol  = calloc(dim, sizeof(int));			       // initializing vectors needed for the routine
	u       = calloc(dim, sizeof(double));
	v       = calloc(dim, sizeof(double));

	for(i=0; i<dim; i++)
	{
		ccost[i] = PyArray_GETPTR2(pcost, 0, i );
	};

	// Solve the LAP //
    lapcost = lap(dim, ccost, PyArray_GETPTR1(rowsol,0), colsol, u, v);
    // check the computed assignment for consistency //
    //TODO make errormasseges for python instead of the c halt stuff in the next function!
    checklap(dim, ccost, PyArray_GETPTR1(rowsol,0), colsol, u, v);

	free(ccost);
	free(colsol);
	free(u);
	free(v);

    return rowsol;

};

/*static PyObject *fit_rotmatrix(PyObject *self, PyObject *args)
{
	return 9;
}//*/
