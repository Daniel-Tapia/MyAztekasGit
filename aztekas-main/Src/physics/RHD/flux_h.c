/**
 * @file /RHD/flux_h.c
 *
 * @author Alejandro Aguayo-Ortiz
 *
 */

#include"main.h"
    
void Prim2FluxH(double *f, double *v, double *u, gauge_ *local_grid)
{
   int i;
   double rho, p, v_cov[3], v_con[3];
   double D, U, tau, Lorentz, h, cs, VV;
   double W[3][3], S_cov[3], S_con[3], V[3];
   double gamma, beta, lapse, vel;
   eos_ eos;

   gamma = local_grid->gamma_con[2][2];
   beta  = local_grid->beta_con[2];
   lapse = local_grid->lapse;

   // Density and Pressure
   rho = u[RHO];
   p   = u[PRE];

   // Covariant components of the 3-velocity
#if DIM == 1
   v_cov[0] = u[VX1];
   v_cov[1] = 0.0;
   v_cov[2] = 0.0;
#elif DIM == 2
   v_cov[0] = u[VX1];
   v_cov[1] = u[VX2];
   v_cov[2] = 0.0;
#elif DIM == 3 || DIM == 4
   v_cov[0] = u[VX1];
   v_cov[1] = u[VX2];
   v_cov[2] = u[VX3];
#endif

   // Contravariant components of the 3-velocity
   v_con[0] = local_grid->gamma_con[0][0]*v_cov[0] + \
              local_grid->gamma_con[0][1]*v_cov[1] + \
              local_grid->gamma_con[0][2]*v_cov[2];
   v_con[1] = local_grid->gamma_con[1][0]*v_cov[0] + \
              local_grid->gamma_con[1][1]*v_cov[1] + \
              local_grid->gamma_con[1][2]*v_cov[2];
   v_con[2] = local_grid->gamma_con[2][0]*v_cov[0] + \
              local_grid->gamma_con[2][1]*v_cov[1] + \
              local_grid->gamma_con[2][2]*v_cov[2];

   // Contraction v_i v^i
   VV = v_con[0]*v_cov[0] + v_con[1]*v_cov[1] + v_con[2]*v_cov[2];

   if(VV > 1.0)
   {
      VV = 0.99;
   }

   // Lorentz Factor
   Lorentz = 1/sqrt(1 - VV);

   // Equation of State
   EoS(&eos,u,local_grid);
   h  = eos.h;
   cs = eos.cs;

   // Define conservative density D, variable U and conservative energy \tau
   D   = rho*Lorentz;
   U   = rho*h*Lorentz*Lorentz - p;
   tau = U - D;

   // Velocity affected by the coordinates
   V[2] = lapse*v_con[2] - beta;

   // Compute the covariant and contravariant components of the 3-momentum
   S_con[2] = rho*eos.h*Lorentz*Lorentz*v_con[2];

   for(i = 0; i < 3; i++)
   {
      S_cov[i] = rho*eos.h*Lorentz*Lorentz*v_cov[i];
   }

   // Compute useful 2-tensor W^i_j (see BHAC article)
   W[2][0] = S_con[2]*v_cov[0];
   W[2][1] = S_con[2]*v_cov[1];
   W[2][2] = S_con[2]*v_cov[2] + p;

   // Compute fluxes
   f[DEN] = D*V[2];
   f[ENE] = lapse*(S_con[2] - v_con[2]*D) - beta*tau;
   f[MX1] = lapse*W[2][0] - beta*S_cov[0];
   f[MX2] = lapse*W[2][1] - beta*S_cov[1];
   f[MX3] = lapse*W[2][2] - beta*S_cov[2];

   // Computed characteristic velocities
   vel   = v_con[2];
   double cs2  = cs*cs;
   double vel2 = vel*vel;

   v[0] = (lapse/(1.0 - VV*cs2))*(vel*(1.0 - cs2) + sqrt(cs2*(1.0 - VV)*(gamma*(1.0 - VV*cs2) - vel2*(1.0 - cs2)))) - beta;
   v[1] = (lapse/(1.0 - VV*cs2))*(vel*(1.0 - cs2) - sqrt(cs2*(1.0 - VV)*(gamma*(1.0 - VV*cs2) - vel2*(1.0 - cs2)))) - beta;
   v[2] = lapse*vel - beta;
}
