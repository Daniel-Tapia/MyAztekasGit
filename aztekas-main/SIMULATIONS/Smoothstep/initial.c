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

   float Smoothstep(float t);

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
//	  IMPLEMENTACIOON DE RUIDO A COS^2(THETA)

//      unif1 = (rand()/Rmax);//pow(-1,rand() % 2)*(rand()/Rmax);
//      printf("unif1 en %d es:	%.6f\n", j, unif1);
//      unif2 = (rand()/Rmax);//pow(-1,rand() % 2)*(rand()/Rmax);
//      printf("unif2 en %d es:	%.6f\n", j, unif2);
//      Noise = 0; //										SIN RUIDO
//      Noise = densamp*((2*unif1)-1); //							RUIDO UNIFORME
//      Noise = densamp*sqrt(-2.*log(unif1))*cos(2*Pi*unif2); // 				RUIDO GAUSSIANO
//      printf("Noise en %d es:	%.6f\n", j, Noise);
//      bounddenst[j] = density_0*(1.0 - (deltarho)*pow(cos(grid.X2[j]),2.0))*(1+Noise); //	DENSIDAD CON RUIDO
//      printf("bounddenst en %d es:      %.6f\n", j, bounddenst[j]);
//      boundpress[j] = pow(bounddenst[j],K)/K; //						PRESIOON
//      printf("boundpress en %d es:      %.6f\n", j, boundpress[j]);
                                                 
//        IMPLEMENTACIOON DE TRANSICIOON SIGMOIDE
      bounddenst[j] = Smoothstep(grid.X2[j]); //						DENSIDAD
      boundpress[j] = pow(bounddenst[j],K)/K; //						PRESIOON
//      printf("%.6f\n", theta_p);
//      printf("%.6f	%.6f	%.6f\n", grid.X2[j], bounddenst[j], boundpress[j]);
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

float Smoothstep(float t)
{
   double Pi = 3.14159265358979323846;
   //theta_p = Pi/2.0 - Pi/8.0 -Pi/4.0;//(theta_e+tr_width);

   if (t<polarrho)
   {
      return density_0*(1.0-deltarho);
   }
   else if (polarrho<=t && t<=polarrho+tr_width)
   {
      return density_0*(1.0-deltarho) + density_0*deltarho*(3.0*pow((t-polarrho)/tr_width, 2.0) - 2.0*pow((t-polarrho)/tr_width, 3.0));
   }
   else if (polarrho+tr_width<t)// || (fabs(t-Pi/2.0) < 0.0001))
   {
      return density_0;
   }
}
