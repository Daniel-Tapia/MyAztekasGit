/**
 * @file /common/print_values.c
 *
 * @authors Alejandro Aguayo-Oritz and Emilio Tejeda
 *
 * @brief Print values.
 */

#include"main.h"

void Print_Values(double *tprint, double *dtprint, int *itprint)
{
   /**
    * Print the Mass Accretion Rate \f[ \dot{M} \f]
    */
   #if MDOT == TRUE
   Mass_Accretion_Rate(U);
   #endif

   if(grid.time >= *tprint || CHECK_NAN == TRUE)
   {
      printf("Time = %e, dt = %e\n",grid.time,dt);
      if(graf == 1)
      {
         if(binary == 1)
         {
            Output1_bin(itprint);
         }
         else
         {
            Output1(itprint);
         }
      }
      else if(graf == 2)
      {
         if(binary == 1)
         {
            Output2_bin(itprint);
         }
         else
         {
            Output2(itprint);
         }
      }
      else if(graf == 3)
      {
         if(binary == 1)
         {
            Output3_bin(itprint);
         }
         else
         {
            Output3(itprint);
         }
      }

      /**
       * Increase the time of printing by dtprint and increase
       * the humber of the output file
       */
      *tprint = *tprint + *dtprint;
      ++*itprint;

      /**
       * Forcing termination
       */
      Alternative_Termination();
   }
}
