/*
 * File Name : initial.c
 * Description : aztekas initial module for Choked Accretion
 * Creation Date : 27-09-2019
 * Last Modified : 27-09-2019 09:57:48
 * Created By : Alejandro Aguayo-Ortiz
 */

#include"main.h"

void Initial()
{
   int n, i, j, k, cell;
   double unif1, unif2, Noise, Rmax, Pi;

   Rmax = 2147483647;
   Pi = 3.14159265358979323846;
   //Initialize grid.time
   grid.time = 0.0;

   //Initialize dt
   dt = 0.0;

   //Ex1 = (double *)malloc((Nx2+1)*sizeof(double));
   //Ex2 = (double *)malloc((Nx2+1)*sizeof(double));
   bounddenst = (double *)malloc((Nx2+1)*sizeof(double));
   boundpress = (double *)malloc((Nx2+1)*sizeof(double));

   for( j = 0; j <= Nx2; j++ )
   {
      unif1 = (rand()/Rmax);//pow(-1,rand() % 2)*(rand()/Rmax);
//      printf("unif1 en %d es:	%.6f\n", j, unif1);
      unif2 = (rand()/Rmax);//pow(-1,rand() % 2)*(rand()/Rmax);
//      printf("unif2 en %d es:	%.6f\n", j, unif2);
      Noise = densamp*sqrt(-2.*log(unif1))*cos(2*Pi*unif2);
      //Noise = 0;
//      printf("Noise en %d es:	%.6f\n", j, Noise);
      bounddenst[j] = density_0*(1.0 - (deltarho)*pow(cos(grid.X2[j]),2.0))*(1+Noise);
//      printf("bounddenst en %d es:	%.6f\n", j, bounddenst[j]);
      boundpress[j] = pow(bounddenst[j],K)/K;
//      printf("boundpress en %d es:	%.6f\n", j, boundpress[j]);
   }

#if DIM == 1 

   ///////////////////////////
   //-------Bondi-1D--------//
   ///////////////////////////

   for(i = 0; i <= Nx1; i++)
   {
      U(0,i) = density_0;
      U(1,i) = pressure_0;
      U(2,i) = velocity_0;
   }

#elif DIM == 2

   ///////////////////////////
   //-------Bondi-2D--------//
   ///////////////////////////

   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         U(0,i,j) = density_0;
         U(1,i,j) = pressure_0;
         U(2,i,j) = velocity_0;
         U(3,i,j) = 0.0;
      }
   }

#elif DIM == 4

   ///////////////////////////
   //-------Bondi-2D--------//
   ///////////////////////////

   for(i = 0; i <= Nx1; i++)
   {
      for(j = 0; j <= Nx2; j++)
      {
         U(0,i,j) = density_0;
         U(1,i,j) = pressure_0;
         U(2,i,j) = velocity_0;
         U(3,i,j) = 0.0;
         U(4,i,j) = 0.0;
      }
   }

#endif
}
