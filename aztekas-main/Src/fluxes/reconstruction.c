/**
 * @file /fluxes/reconstruction.c
 *
 * @author Alejandro Aguayo-Ortiz
 *
 * @brief Define the reconstruction vectors \f$ (Up, Um) \f$
 * at the interfaces \f$ i \pm 1/2 \f$, respectively.
 */

#include"main.h"

#if DIM == 1

void Primitive_Reconstruction()
{
   int reconst;
   double dup, dum, sx;

#ifdef _OPENMP
   #pragma omp parallel shared(U,U1p,U1m)
   #pragma omp for private(reconst,dup,dum,sx)
#endif
   for(int i = 1; i <= Nx1-1; i++)
   {
      reconst = RECONST;
#if x1min_exc == TRUE
      if (i <= 2*gc && RECONST == WENO5)
      {
         reconst = MC;
      }
      if (i <= 2*gc && RECONST != WENO5)
      {
         reconst = GODUNOV;
      }
#endif

      for(int n = 0; n < eq; n++)
      {
         if(reconst == WENO5)
         {
            U1m(n,i) = Weno5(U0(n,i+2),\
                             U0(n,i+1),\
                             U0(n, i ),\
                             U0(n,i-1),\
                             U0(n,i-2));
            U1p(n,i) = Weno5(U0(n,i-2),\
                             U0(n,i-1),\
                             U0(n, i ),\
                             U0(n,i+1),\
                             U0(n,i+2));
         }
         else
         {
            dum = U0(n, i ) - U0(n,i-1);
            dup = U0(n,i+1) - U0(n, i );

            sx  = Limiter(dup,dum,reconst);

            U1m(n,i) = U0(n,i) - 0.5*sx;
            U1p(n,i) = U0(n,i) + 0.5*sx;
         }
      }
   }
}

#elif DIM == 2 || DIM == 4

void Primitive_Reconstruction()
{
   int reconst;
   double dup, dum, sx;

#ifdef _OPENMP
   #pragma omp parallel shared(U,U1p,U1m,U2p,U2m)
   #pragma omp for private(reconst,dup,dum,sx) collapse(2)
#endif
   for(int j = 1; j <= Nx2-1; j++)
   {
      for(int i = 1; i <= Nx1-1; i++)
      {
         reconst = RECONST;
#if x1min_exc == TRUE
         if (i <= gc+1 && RECONST == WENO5)
         {
            reconst = MC;
         }
         if (i <= gc+1 && RECONST != WENO5)
         {
            reconst = GODUNOV;
         }
#endif
#if x1max_exc == TRUE
         if (i >= Nx1-gc-1 && RECONST == WENO5)
         {
            reconst = MC;
         }
         if (i >= Nx1-gc-1 && RECONST != WENO5)
         {
            reconst = GODUNOV;
         }
#endif
         for(int n = 0; n < eq; n++)
         {
            if(reconst == WENO5)
            {
               U1m(n,i,j) = Weno5(U0(n,i+2,j),\
                                  U0(n,i+1,j),\
                                  U0(n, i ,j),\
                                  U0(n,i-1,j),\
                                  U0(n,i-2,j));
               U1p(n,i,j) = Weno5(U0(n,i-2,j),\
                                  U0(n,i-1,j),\
                                  U0(n, i ,j),\
                                  U0(n,i+1,j),\
                                  U0(n,i+2,j));

               U2m(n,i,j) = Weno5(U0(n,i,j+2),\
                                  U0(n,i,j+1),\
                                  U0(n,i, j ),\
                                  U0(n,i,j-1),\
                                  U0(n,i,j-2));
               U2p(n,i,j) = Weno5(U0(n,i,j-2),\
                                  U0(n,i,j-1),\
                                  U0(n,i, j ),\
                                  U0(n,i,j+1),\
                                  U0(n,i,j+2));
            }
            else
            {
               dum = U0(n, i ,j) - U0(n,i-1,j);
               dup = U0(n,i+1,j) - U0(n, i ,j);

               sx  = Limiter(dup,dum,reconst);

               U1m(n,i,j) = U0(n,i,j) - 0.5*sx;
               U1p(n,i,j) = U0(n,i,j) + 0.5*sx;


               dum = U0(n,i, j ) - U0(n,i,j-1);
               dup = U0(n,i,j+1) - U0(n,i, j );

               sx  = Limiter(dup,dum,reconst);

               U2m(n,i,j) = U0(n,i,j) - 0.5*sx;
               U2p(n,i,j) = U0(n,i,j) + 0.5*sx;

            }
         }
      }
   }
}

#elif DIM == 3

int Reconst3D(double *u, lim_ *l, int *I)
{
   int n, i, j, k;
   int reconst;
   double dup2, dup1, dum1, dum2;

   i = I[0];
   j = I[1];
   k = I[2];
   reconst = RECONST;

   for(n = 0; n < eq; n++)
   {
      l->ux[n] = u(n,i,j,k);

      if(reconst == WENO5)
      {
         l->ux1m[0*eq + n] = Weno5(u(n,i-3,j,k),\
                                   u(n,i-2,j,k),\
                                   u(n,i-1,j,k),\
                                   u(n, i ,j,k),\
                                   u(n,i+1,j,k));
         l->ux1m[1*eq + n] = Weno5(u(n,i+2,j,k),\
                                   u(n,i+1,j,k),\
                                   u(n, i ,j,k),\
                                   u(n,i-1,j,k),\
                                   u(n,i-2,j,k));

         l->ux1p[0*eq + n] = Weno5(u(n,i-2,j,k),\
                                   u(n,i-1,j,k),\
                                   u(n, i ,j,k),\
                                   u(n,i+1,j,k),\
                                   u(n,i+2,j,k));
         l->ux1p[1*eq + n] = Weno5(u(n,i+3,j,k),\
                                   u(n,i+2,j,k),\
                                   u(n,i+1,j,k),\
                                   u(n, i ,j,k),\
                                   u(n,i-1,j,k));

         l->ux2m[0*eq + n] = Weno5(u(n,i,j-3,k),\
                                   u(n,i,j-2,k),\
                                   u(n,i,j-1,k),\
                                   u(n,i, j ,k),\
                                   u(n,i,j+1,k));
         l->ux2m[1*eq + n] = Weno5(u(n,i,j+2,k),\
                                   u(n,i,j+1,k),\
                                   u(n,i, j ,k),\
                                   u(n,i,j-1,k),\
                                   u(n,i,j-2,k));

         l->ux2p[0*eq + n] = Weno5(u(n,i,j-2,k),\
                                   u(n,i,j-1,k),\
                                   u(n,i, j ,k),\
                                   u(n,i,j+1,k),\
                                   u(n,i,j+2,k));
         l->ux2p[1*eq + n] = Weno5(u(n,i,j+3,k),\
                                   u(n,i,j+2,k),\
                                   u(n,i,j+1,k),\
                                   u(n,i, j ,k),\
                                   u(n,i,j-1,k));

         l->ux3m[0*eq + n] = Weno5(u(n,i,j,k-3),\
                                   u(n,i,j,k-2),\
                                   u(n,i,j,k-1),\
                                   u(n,i,j, k ),\
                                   u(n,i,j,k+1));
         l->ux3m[1*eq + n] = Weno5(u(n,i,j,k+2),\
                                   u(n,i,j,k+1),\
                                   u(n,i,j, k ),\
                                   u(n,i,j,k-1),\
                                   u(n,i,j,k-2));

         l->ux3p[0*eq + n] = Weno5(u(n,i,j,k-2),\
                                   u(n,i,j,k-1),\
                                   u(n,i,j, k ),\
                                   u(n,i,j,k+1),\
                                   u(n,i,j,k+2));
         l->ux3p[1*eq + n] = Weno5(u(n,i,j,k+3),\
                                   u(n,i,j,k+2),\
                                   u(n,i,j,k+1),\
                                   u(n,i,j, k ),\
                                   u(n,i,j,k-1));
      }
      else
      {
         dum2 = u(n,i-1,j,k) - u(n,i-2,j,k);
         dum1 = u(n,i  ,j,k) - u(n,i-1,j,k);
         dup1 = u(n,i+1,j,k) - u(n,i  ,j,k);
         dup2 = u(n,i+2,j,k) - u(n,i+1,j,k);

         l->sx1[0*eq + n] = Limiter(dum1,dum2,reconst);
         l->sx1[1*eq + n] = Limiter(dup1,dum1,reconst);
         l->sx1[2*eq + n] = Limiter(dup2,dup1,reconst);

         l->ux1m[0*eq + n] = u(n,i-1,j,k) + 0.5*l->sx1[0*eq + n];
         l->ux1m[1*eq + n] = u(n,i  ,j,k) - 0.5*l->sx1[1*eq + n];

         l->ux1p[0*eq + n] = u(n,i  ,j,k) + 0.5*l->sx1[1*eq + n];
         l->ux1p[1*eq + n] = u(n,i+1,j,k) - 0.5*l->sx1[2*eq + n];

         dum2 = u(n,i,j-1,k) - u(n,i,j-2,k);
         dum1 = u(n,i,j  ,k) - u(n,i,j-1,k);
         dup1 = u(n,i,j+1,k) - u(n,i,j  ,k);
         dup2 = u(n,i,j+2,k) - u(n,i,j+1,k);

         l->sx2[0*eq + n] = Limiter(dum1,dum2,reconst);
         l->sx2[1*eq + n] = Limiter(dup1,dum1,reconst);
         l->sx2[2*eq + n] = Limiter(dup2,dup1,reconst);

         l->ux2m[0*eq + n] = u(n,i,j-1,k) + 0.5*l->sx2[0*eq + n];
         l->ux2m[1*eq + n] = u(n,i,j  ,k) - 0.5*l->sx2[1*eq + n];

         l->ux2p[0*eq + n] = u(n,i,j  ,k) + 0.5*l->sx2[1*eq + n];
         l->ux2p[1*eq + n] = u(n,i,j+1,k) - 0.5*l->sx2[2*eq + n];

         dum2 = u(n,i,j,k-1) - u(n,i,j,k-2);
         dum1 = u(n,i,j,k  ) - u(n,i,j,k-1);
         dup1 = u(n,i,j,k+1) - u(n,i,j,k  );
         dup2 = u(n,i,j,k+2) - u(n,i,j,k+1);

         l->sx3[0*eq + n] = Limiter(dum1,dum2,reconst);
         l->sx3[1*eq + n] = Limiter(dup1,dum1,reconst);
         l->sx3[2*eq + n] = Limiter(dup2,dup1,reconst);

         l->ux3m[0*eq + n] = u(n,i,j,k-1) + 0.5*l->sx3[0*eq + n];
         l->ux3m[1*eq + n] = u(n,i,j,k  ) - 0.5*l->sx3[1*eq + n];

         l->ux3p[0*eq + n] = u(n,i,j,k  ) + 0.5*l->sx3[1*eq + n];
         l->ux3p[1*eq + n] = u(n,i,j,k+1) - 0.5*l->sx3[2*eq + n];
      }
   }

   return 0;
}

#endif
