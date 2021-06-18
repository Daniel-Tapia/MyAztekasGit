/*
 * File Name : boundaries.c
 * Description : aztekas boundaries module for Choked Accretion
 * Creation Date : 27-09-2019
 * Last Modified : 28-10-2019 18:20:53
 * Created By : Alejandro Aguayo-Ortiz
 */

#include"main.h"

void User_Boundaries(double *B)
{
#if DIM == 1

   for(int i = 0; i <= Nx1; i++)
   {
      if(i >= Nx1-gc)
      {
         B(RHO,i) = density_0;
         B(PRE,i) = pressure_0;
      }
   }

#elif DIM == 2
//   printf("Inside User_Boundaries Dim=2\n");
   for(int j = 0; j <= Nx2; j++)
   {
      for(int i = 0; i <= Nx1; i++)
      {
         if(i >= Nx1-gc)
         {
            //B(RHO,i,j) = density_0*gtheta(grid.X2[j]);
            //B(RHO,i,j) = density_0*gtheta2(grid.X2[j], Nx2);
               //B(RHO,i,j) = density_0*gtheta3(grid.X2[j], j);
	    //B(RHO,i,j) = density_0;
            B(RHO,i,j) = bounddenst[j];
            
	    //B(PRE,i,j) = pow(B(RHO,i,j),K)/K;
            B(PRE,i,j) = boundpress[j];
            //B(VX1,i,j) = 0.0;
            //B(VX2,i,j) = 0.0;
         }
      }
   }

#elif DIM == 4

   for(int j = 0; j <= Nx2; j++)
   {
      for(int i = 0; i <= Nx1; i++)
      {
         if(i >= Nx1-gc)
         {
            B(RHO,i,j) = density_0*gtheta(grid.X2[j]);
            B(PRE,i,j) = pow(B(RHO,i,j),K)/K;
         }
      }
   }

#endif
}
