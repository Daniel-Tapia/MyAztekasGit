/**
 * @file main.h
 *
 * @author Alejandro Aguayo-Ortiz
 *
 * @brief Main function, headers and variable declaration.
 */

#ifdef _OPENMP
   #include<omp.h>
   int MAX_NUM_THREADS;
   double start;
#else
   #include<time.h>
   clock_t start;
#endif

#include<unistd.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif
#ifndef M_PI_2
    #define M_PI_2 1.57079632679489661923
#endif
#ifndef RAND_MAX
    #define RAND_MAX 2147483647
#endif
#ifndef RNUM
    #define RNUM (2*(rand()/RAND_MAX)-1)
#endif


#include"mesh.h"
#include"physics.h"

#include"initial.h"
#include"integration.h"
#include"boundaries.h"
#include"limiters.h"
#include"flux.h"

#include"const.h"
#include"macros.h"
#include"io.h"
#include"user_param.h"
