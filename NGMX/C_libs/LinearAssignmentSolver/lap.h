/************************************************************************
*
*  lap.h
   version 1.0 - 21 june 1996
   author  Roy Jonker, MagicLogic Optimization Inc.
   
   header file for LAP
*
**************************************************************************/

//#include <typedefs.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#ifndef __LAP_H
#define __LAP_H

#ifdef CPP
#define EXPORT extern "C"
#endif

#ifndef CPP
#define EXPORT extern
#endif

/*************** CONSTANTS  *******************/

  #define BIG 100000

/*************** TYPES      *******************/

  typedef int row;
  typedef int col;
  typedef double cost;

/*************** FUNCTIONS  *******************/

//EXPORT         cost lap(int dim, cost **assigncost, int *rowsol, int *colsol, cost *u, cost *v);
cost lap(int dim, cost **assigncost, int *rowsol, int *colsol, cost *u, cost *v);

//EXPORT         void checklap(int dim, cost **assigncost, int *rowsol, int *colsol, cost *u, cost *v);
void checklap(int dim, cost **assigncost, int *rowsol, int *colsol, cost *u, cost *v);

#endif /* __LAP_H */
