/**
 * @file functions.c
 * 
 * @brief Useful functions
 *
 * @author Alejandro Aguayo-Ortiz
 *
 * @date 28-07-2020 - 20:43:20
 *
 * E-mail: aaguayoo92@ciencias.unam.mx
 *
 * Created on: 28-04-2020
 */

#include"main.h"

void Mass_Accretion_Rate(double *B)
{
   FILE *file_time, *file_domain;
   char name_time[100], name_domain[100];
   char ext[50], dat[50];
   int num;

   gauge_ local_grid;
   #if COORDINATES == SPHERICAL
   double time, r, th, thm, phi;
   double y_rr, y_thth, y_phiphi;
   double y_rth, y_rphi, y_thphi;
   double beta_r, beta_th, beta_phi;
   double lapse, dety;
   double rho, pre, vr, vth, vphi;
   double Vr, Vth, Vphi;
   double VV, W, Ur;
   double Mdot, Min, Mej, Vol;
   double Mean[Nx1+1];
   double Mdot_Acc = 0.0;
   double Mdot_Mean_Time = 0.0;
   double Delta;
   double Mb = 4.0*M_PI*0.25*pow(2.0/(5.0 - 3.0*K),(5.0 - 3.0*K)/(2.0*(K - 1.0)));
   double rs = (5.0 - 3.0*K)/4.0;
   #endif

   local_grid.x[0] = grid.time;
   time            = local_grid.x[0];
   Mdot = 0.0;

   num = itprint;

   strcpy(ext,".dat");
   snprintf(dat,8,"%d",num);

   strcpy(name_time,outputdirectory);
   strcat(name_time,"Mdot/");
   strcat(name_time,"Mdot_Mean_Time");
   strcat(name_time,ext);

   strcpy(name_domain,outputdirectory);
   strcat(name_domain,"Mdot/");
   strcat(name_domain,"Mdot_Mean_Domain_");
   strcat(name_domain,dat);
   strcat(name_domain,ext);

   Mdot = 0.0;
   Min  = 0.0;
   Mej  = 0.0;

   for(int i = gc; i <= Nx1-gc; i++)
   {
      for(int j = gc+1; j <= Nx2-gc-0; j++)
      {
      #if DIM == 1                                                                    
                                                                                      
         local_grid.x[1] = grid.X1[i];                                                
         local_grid.x[2] = 0.0;                                                       
         local_grid.x[3] = 0.0;                                                       
         #if COORDINATES == SPHERICAL                                                 
         local_grid.x[2] = M_PI_2;                                                    
         #endif                                                                       
                                                                                      
      #elif DIM == 2 || DIM == 4                                                      
                                                                                      
         local_grid.x[1] = grid.X1[i];                                                
         local_grid.x[2] = grid.X2[j];                                                
         local_grid.x[3] = 0.0;                                                       
         #if POLAR == TRUE                                                            
         local_grid.x[2] = M_PI_2;                                                    
         #endif                                                                       
                                                                                      
      #elif DIM == 3                                                                  
                                                                                      
         local_grid.x[1] = grid.X1[i];                                                
         local_grid.x[2] = grid.X2[j];                                                
         local_grid.x[3] = grid.X3[k];                                                
                                                                                      
      #endif 

      #if PHYSICS == RHD
         Get_Metric_Components(&local_grid);
      #endif 

         r   = local_grid.x[1];
         th  = local_grid.x[2];
         phi = local_grid.x[3];

         y_rr     = 1.0;
         y_thth   = 1.0/pow(r,2.0);
         y_phiphi = 1.0;

         beta_r   = local_grid.beta_con[0];
         beta_th  = local_grid.beta_con[1];
         beta_phi = local_grid.beta_con[2];

         lapse = local_grid.lapse;
         dety  = r*r*sin(th);

         rho  = B(RHO,i,j);
         pre  = B(PRE,i,j);
         vr   = B(VX1,i,j);
         vth  = B(VX2,i,j);
         vphi = B(VX3,i,j);

         Vr   = y_rr*vr;
         Vth  = y_thth*vth;
         Vphi = y_phiphi*vphi;
         
         VV = Vr*vr + Vth*vth + Vphi*vphi;
         W  = 1.0/sqrt(1.0 - VV);

         Ur = Vr;

         if(j == gc+1 || j == Nx2-gc-1)
         {
            Mdot = Mdot - 4.0*M_PI*rho*Ur*dety*(grid.X2p[j] - grid.X2m[j])/2.0;
            if(i == Nx1-gc-1 && Ur <= 0.0)
            {
               Min = Min - 4.0*M_PI*rho*Ur*dety*((grid.X2p[j] - grid.X2m[j])/2.0);
            }
            else if(i == Nx1-gc-1 && Ur > 0.0)
            {
               Mej = Mej + 4.0*M_PI*rho*Ur*dety*((grid.X2p[j] - grid.X2m[j])/2.0);
            }
         }
         else
         {
            Mdot = Mdot - 4.0*M_PI*rho*Ur*dety*(grid.X2p[j] - grid.X2m[j]);
            if(i == Nx1-gc-1 && Ur <= 0.0)
            {
               Min = Min - 4.0*M_PI*rho*Ur*dety*((grid.X2p[j] - grid.X2m[j]));
            }
            else if(i == Nx1-gc-1 && Ur > 0.0)
            {
               Mej = Mej + 4.0*M_PI*rho*Ur*dety*((grid.X2p[j] - grid.X2m[j]));
            }
         }
      }
      Mean[i] = Mdot;
      Mdot    = 0.0;
      Vol     = 0.0;
   }

   // Clean the Mdot/Mdot_Mean_Time.dat file if t = 0
   file_time   = fopen(name_time,"a");
   if(time == 0.0)
   {
      char dum1[50] = "cat /dev/null > ";
      strcat(dum1," ");
      strcat(dum1,name_time);
      system(dum1);
      Mdot_Mean = 0.0;
   }

   // Computes both Mdot_Mean all over the domain
   // and Mdot at a certain radius R(4)
   for(int i = gc + 1; i <= Nx1-gc-1; i++)
   {
      local_grid.x[1] = grid.X1[i];                                                
      local_grid.x[2] = 0.0;                                                       
      local_grid.x[3] = 0.0;                                                       
      #if COORDINATES == SPHERICAL                                                 
      local_grid.x[2] = M_PI_2;                                                    
      #endif                                                                       

      #if PHYSICS == RHD
      Get_Metric_Components(&local_grid);
      #endif 

      r   = local_grid.x[1];
      th  = local_grid.x[2];
      phi = local_grid.x[3];

      y_rr     = local_grid.gamma_con[0][0];
      y_thth   = local_grid.gamma_con[1][1];
      y_phiphi = local_grid.gamma_con[2][2];
      y_rth    = local_grid.gamma_con[0][1];
      y_rphi   = local_grid.gamma_con[0][2];
      y_thphi  = local_grid.gamma_con[1][2];

      beta_r   = local_grid.beta_con[0];
      beta_th  = local_grid.beta_con[1];
      beta_phi = local_grid.beta_con[2];

      lapse = local_grid.lapse;
      dety  = r*r*sin(th);
                                                                                      
      Mdot    += Mean[i];//*dety*(grid.X1p[i] - grid.X1m[i]);

      if(local_grid.x[1] > rs && local_grid.x[1] < rs + (grid.X1p[i] - grid.X1m[i]))
         Mdot_Acc = Mean[i]; 
   }

   Vol = (Nx1-2);
   // Computes the difference between Mdot_mean(n) and Mdot_mean(n-1)
   Delta = Mdot/Vol - Mdot_0;

   if(Delta > 0 && plus == FALSE)
   {
      Mdot_Max = Mdot/Vol;
      plus     = TRUE;
      minus    = FALSE;
      count++;
   }

   if(Delta < 0 && minus == FALSE)
   {
      Mdot_Min = Mdot/Vol;
      plus     = FALSE;
      minus    = TRUE;
      count++;
   }

   #if MDOT_END == TRUE
   if(time == 0.0)
   {
      restart_file = FALSE;
   }

   if(restart_simulation == FALSE && fabs(1.0 - Mdot_Min/Mdot_Max) < 1.0e-02 && restart_file == FALSE)
   {
      tprint       = time;
      restart_file = TRUE;
      printf("\n");
      printf("\n");
      printf("Mass accretion rate oscillation error (ΔM/M) drops below %e. Restart file created.\n",1.0e-02);
      printf("Simulation continues until ΔM/M < %e.\n",MDOT_ERR);
      printf("\n");
   }
   else if(fabs(1.0 - Mdot_Min/Mdot_Max) < MDOT_ERR && count > 1000)
   {
      tmax   = time;
      tprint = tmax;
      Mdot_end  = TRUE;
      printf("\n");
      printf("\n");
      printf("Error below %e. Termination. Count = %d\n",MDOT_ERR,count);
      printf("\n");
   }
//   //else if(Black_Hole_Spin >= 0.0 && count > 1000)
//   else if(Black_Hole_Spin >= 0.0)
//   {
//      tmax   = time;
//      tprint = tmax;
//      Mdot_end  = TRUE;
//      printf("\n");
//      printf("\n");
//      printf("Error below %e. Termination. Count = %d\n",MDOT_ERR,count);
//      printf("\n");
//   }
   //if(fabs(1.0 - Mdot_Min/Mdot_Max) < MDOT_ERR && restart_file == FALSE)
//   if(Black_Hole_Spin >= 0.0 && fabs(1.0 - Mdot_Min/Mdot_Max) < MDOT_ERR && restart_file == FALSE)
//   {
//      tprint       = time;
//      //tprint       = time + tmax/5.0;
//      restart_file = TRUE;
//      tmax         = time + tmax/5.0;
//      printf("\n");
//      printf("\n");
//      printf("Error below %e. Restart file created.\n",MDOT_ERR);
//      printf("Simulation continues up to t = %e.\n",tmax);
//      printf("\n");
//   }
//   else if(Black_Hole_Spin >= 0.0 && restart_file == TRUE)
//   {
//      tprint = tmax;
//   }
//   else if(Black_Hole_Spin > 0.0 && restart_file == FALSE && time > 0.0)
//   {
//      tmax   = tmax/5.0;
//      tprint = tmax;
//      restart_file = TRUE;
//   }
//   else if(Mdot == Mdot_0 || Mdot_end == TRUE)
//   {
//      tmax   = time;
//      tprint = tmax;
//      Mdot_end  = TRUE;
//      printf("Mdot = Mdot_0\n");
//   }
//   else if(count > 20000 && fabs(1.0 - Mdot_Min/Mdot_Max) < 10*MDOT_ERR || Mdot_end == TRUE)
//   {
//      tmax   = time;
//      tprint = tmax;
//      Mdot_end  = TRUE;
//      printf("Number of points %d\n",count);
//   }
   #else
   Mdot_end = FALSE;
   #endif

   Mdot_Mean_Time = (Mdot_Max + Mdot_Min)/2.0;
   if(time == 0.0)
   {
      Mdot_Mean_Time = Mdot/Vol;
   }

   if(time >= Mdot_tprint)
   {
      fprintf(file_time,"%e %e %e %e %e %e %e\n",time,Mdot/Vol/Mb,Mdot_Acc/Mb,Mdot_Mean_Time/Mb,fabs(1.0 - Mdot_Min/Mdot_Max),Min/Mb,fabs(Mdot_Mean_Time - Min)/Mb);
      Mdot_tprint = Mdot_tprint + MDOT_TIME;
   }
   fclose(file_time);

   file_domain = fopen(name_domain,"w");
   if(time >= tprint)
   {
      fprintf(file_domain,"###############PARAM###############\n");
      fprintf(file_domain,"%e \n",time);
      fprintf(file_domain,"%d \n",Nx1-2*gc+1);
      fprintf(file_domain,"###################################\n");

      for(int i = gc; i <= Nx1-gc; i++)
      {
         fprintf(file_domain,"%e %e %e %e\n",grid.X1[i],Mb,Mean[i],(fabs(Mean[i]-Mb))/Mb);
      }

      printf("Computed the mass accretion rate: %s\n",name_domain); 
   }
   fclose(file_domain);

   Mdot_0 = Mdot/Vol;

//   if(time >= tmax)
//   {
//      FILE *results;
//      strcpy(name_time,outputdirectory);
//      strcat(name_time,"../");
//      strcat(name_time,"results_a");
//      strcat(name_time,ext);
//
//      results = fopen(name_time,"a");
//      fprintf(results,"%e %e %e %e %e\n",K,cs,x1max,Mdot_Mean_Time);
//      fclose(results);
//
//      strcpy(name_time,outputdirectory);
//      strcat(name_time,"../../");
//      strcat(name_time,"results");
//      strcat(name_time,ext);
//
//      results = fopen(name_time,"a");
//      fprintf(results,"%e %e %e %e %e\n",K,cs,x1max,Mdot_Mean_Time);
//      fclose(results);
//   }
}

double gtheta(double th)
{
   //float  rax, rnum;//sgn;
   
   //int    rnum;
   
   //double noise, epsilon;

   //double dum = (1.0 - rho_atm)*pow(sin(th),2.0) + rho_atm;
   //sgn	= pow(-1,rand() % 2);
   //rnum = rand() % 101; //11
   // noise= (double) rnum * 0.005; //0.01
   
   //rax = RAND_MAX;
   //epsilon = 0.01;
   //rnum = epsilon*(2*(rand()/rax)-1);

   //double dum = 1.0 - (0.01)*pow(cos(th),2.0);
   //double dum = (1.0 - (0.1)*pow(cos(th),2.0))+(sgn*noise*(1.0 - (0.1)*pow(cos(th),2.0)));

   double dum = 1.0 - (0.1)*pow(cos(th),2.0)+(RNUM*densamp*(1.0 - (0.1)*pow(cos(th),2.0)));
   return dum;
}

double gtheta2(double th, int n2)
{
   float randarr[] = {-0.736924, -0.0827, -0.562082, 0.357729, 0.869386, 0.038833, -0.930856, 0.0594, -0.984604, -0.866316, 0.373545, 0.860873, 0.053858, 0.307838, 0.402381, 0.524396, -0.905071, -0.343532, 0.512821, -0.269323, 0.965101, 0.506712, -0.854628, 0.769414, -0.127177, -0.044536, -0.450186, -0.666986, 0.795313, -0.878871, 0.009046, -0.361934, -0.012047, -0.818534, -0.852502, -0.231716, 0.827635, -0.071108, -0.899832, 0.540409, -0.749269, 0.376911, 0.259087, 0.450824, 0.777144, -0.387356, 0.026547, 0.691963, 0.683021, -0.169211, -0.064165, -0.643345, 0.14331, -0.933892, -0.00304, 0.496585, 0.781475, 0.684079, -0.574497, -0.739145, -0.450824, -0.171413, 0.419639, -0.520178, -0.364921, 0.304117, 0.362692, -0.224549, -0.704934, 0.691151, 0.910818, -0.703697, -0.182467, 0.129797, -0.022971, 0.92219, -0.600486, 0.258538, 0.302507, 0.606146, -0.047136, -0.593499, 0.803347, -0.715958, -0.179374, 0.771297, -0.675603, -0.269322, -0.729781, -0.089385, -0.0954, 0.863349, -0.569503, 0.817844, 0.72172, 0.011912, 0.635123, -0.07551, 0.265477, 0.649395, 0.404413, 0.90883, -0.421368, 0.028869, -0.171943, 0.753131, 0.459495, 0.431285, 0.413071, -0.961815, 0.049975, -0.869612, -0.022114, 0.364098, 0.833268, 0.780037, -0.72161, 0.978725, -0.107953, 0.029319, -0.120549, 0.6133, -0.576962, -0.692791, 0.232701, -0.998242, 0.45467, -0.164553, 0.361124};
   float densnoise[n2+1];
   int i;
   double dum;
   float theta;
   for( i = 0 ; i <= n2; i++ )
   {
      theta = (M_PI/(2*n2))*i;
      if (fabs(th-theta) < 0.00001)
      {
         densnoise[i] = (1.0 - (0.1)*pow(cos(th),2.0))+(randarr[i]*densamp*(1.0 - (0.1)*pow(cos(th),2.0)));
         dum = densnoise[i];
      }
   }
   
   return dum;
}

double gtheta3(double th, int j)
{
   double dum = 0.0; //1.0 - (0.1)*pow(cos(th),2.0)+(uninoise1[j]*densamp*(1.0 - (0.1)*pow(cos(th),2.0)));
   return dum;
}

