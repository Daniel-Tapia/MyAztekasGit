/**
 * @file boundaries.c
 *
 * @brief Standard boundary conditions. Outflow, Periodic and Reflection.
 *
 * @author Alejandro Aguayo-Ortiz
 *
 * @date 04-10-2020 - 13:09:44
 *
 * E-mail: aaguayoo92@ciencias.unam.mx
 *
 * Created on: 20-10-19
 */

#include"main.h"

/**
 * Call the default boundary conditions implemented in \a aztekas:
 * - Reflection.
 * - Periodic.
 * - Outflow.
 *
 * Input/Output: double *B
 */
void Boundaries(double *B)
{
#if reflective_x1max == TRUE \
 || reflective_x1min == TRUE \
 || reflective_x2max == TRUE \
 || reflective_x2min == TRUE \
 || reflective_x3max == TRUE \
 || reflective_x3min == TRUE 

   Reflection(B);

#endif

#if periodic_x1 == TRUE \
 || periodic_x2 == TRUE \
 || periodic_x3 == TRUE

   Periodic(B);

#endif

#if outflow_x1max == TRUE \
 || outflow_x1min == TRUE \
 || outflow_x2max == TRUE \
 || outflow_x2min == TRUE \
 || outflow_x3max == TRUE \
 || outflow_x3min == TRUE 

   Outflow(B);

#endif

   User_Boundaries(B);
}
