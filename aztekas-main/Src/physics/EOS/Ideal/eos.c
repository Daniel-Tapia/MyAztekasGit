/**
 * @file eos.c
 *
 * @author Alejandro Aguayo-Ortiz
 *
 * @brief Equation of state an Ideal Gas
 *
 */

#include"main.h"

void EoS(eos_ *eos, double *u, gauge_ *local_grid)
{
   double rho, p;
   rho = u[RHO];
   p   = u[PRE];

   eos->e = p / (rho * (K - 1.0));

#if PHYSICS == HD
   eos->cs = sqrt(K * p / rho);
#elif PHYSICS == RHD
   eos->h  = 1.0 + eos->e + p/rho;
   #if INTEGRATION == PVRS
   eos->dhdrho = - K*p/(rho*rho*(K - 1.0));
   eos->dhdp   = K/(rho*(K - 1.0));
   #endif
   eos->cs = sqrt(K * p / (rho * eos->h));
#endif
}
