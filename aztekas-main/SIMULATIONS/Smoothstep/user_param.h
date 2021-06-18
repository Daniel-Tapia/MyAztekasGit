/*
 * File Name : user_param.h
 * Description : aztekas user parameters header file for Choked Accretion
 * Creation Date : 27-09-2019
 * Last Modified : 28-10-2019 18:20:57
 * Created By : Alejandro Aguayo-Ortiz
 */

#include"macros.h"

#define outflow_x1max      TRUE
#define outflow_x1min      TRUE
#define reflective_x2max   TRUE
#define reflective_x2min   TRUE

#define RECONST            MC
#define FLUX               HLL
#define GRID               LOGMESH
#define lfac               1.0

#define MDOT               TRUE
#define MDOT_END           TRUE
#define MDOT_TIME          1.0
#define MDOT_ERR           1.0e-03

double density_0, pressure_0, velocity_0;
double deltarho, densamp;
double polarrho, tr_width;

char last[50];
int count, plus, minus, restart_file;
int Mdot_end,				count;
double Mdot_Max, Mdot_Min;
double Mdot_0, Mdot_tprint, Mdot_Mean;

void Mass_Accretion_Rate(double *B);
double gtheta(double th);
double gtheta2(double th, int n2);
double gtheta3(double th, int j);

//double *Ex1;
//double *Ex2;
double *bounddenst;
double *boundpress;
double prev;
//double prev;
