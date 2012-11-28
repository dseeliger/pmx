/*
 * This file contains some functions to simplify memory handling
 * of matrices and vectors. It was taken from   
 * the Jonker / Volgenant LAP code and adapted slightly
 * see http://www.magiclogic.com/assignment.html for details
 */

#include <memory.h>
#include "lap.h"

#include <stdio.h>
#include <stdlib.h>

void terminate(void);
int *alloc_ints(int len);
void free_ints(int *ints);
cost *alloc_costs(int len);
void free_costs(cost *costs);
cost   **alloc_costmatrix(int width, int height);
void free_costmatrix(cost **cmat, int width, int height);

