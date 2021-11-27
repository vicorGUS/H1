/*
 H1main.c
 
 Created by Anders Lindman on 2013-10-31.
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <float.h>
#include "H1lattice.h"
#include "H1potential.h"

void velocity_verlet(int n_timesteps, int nbr_atoms, double v[nbr_atoms][3], double pos[nbr_atoms][3], double a, double dt, double m, double T[n_timesteps+1], double P[n_timesteps+1], double *E_pot, double *E_kin, double tau_T, double tau_P, double T_eq, double P_eq, double q1[n_timesteps][2], double q2[n_timesteps][2], double q3[n_timesteps][2]);

void arange(double *array, double start, int len_t, double dt){
    for(int i = 0; i < len_t; i++){
    array[i] = start + i*dt;
    }
}

void write_to_file(char *fname, double *time_array,
           double *E, int n_points)
{
    FILE *fp = fopen(fname, "w");
    fprintf(fp, "Time, Energy\n");
    for(int i = 0; i < n_points; ++i){
        fprintf(fp, "%f,%f\n", time_array[i], E[i]);
    }
    fclose(fp);
}

void trajectiores(char *fname, int n_points, double q1[n_points][2], double q2[n_points][2], double q3[n_points][2])
{
    FILE *fp = fopen(fname, "w");
    for(int i = 0; i < n_points; ++i){
        fprintf(fp, "%f,%f,%f,%f,%f,%f\n", q1[i][0], q1[i][1], q2[i][0], q2[i][1], q3[i][0], q3[i][1]);
    }
    fclose(fp);
}

/* Main program */
int main()
{
    
    /*
     Code for generating a uniform random number between 0 and 1. srand should only
     be called once.
    */
    /*
     srand(time(NULL));
     double random_value;
     random_value = (double) rand() / (double) RAND_MAX;
    */
    
    /*
     Descriptions of the different functions in the files H1lattice.c and
     H1potential.c are listed below.
    */
    
    /*
     Function that generates a fcc lattice in units of [Å]. Nc is the number of
     primitive cells in each direction and a0 is the lattice parameter. The
     positions of all the atoms are stored in pos which should be a matrix of the
     size N x 3, where N is the number of atoms. The first, second and third column
     correspond to the x,y and z coordinate respectively.
    */
    /*
     init_fcc(pos, Nc, a0);
    */
    
    /*
     Function that calculates the potential energy in units of [eV]. pos should be
     a matrix containing the positions of all the atoms, L is the length of the
     supercell and N is the number of atoms.
    */
    /*
     double energy;
     energy = get_energy_AL(pos, L, N);
    */
    
    /*
     Function that calculates the virial in units of [eV]. pos should be a matrix
     containing the positions of all the atoms, L is the length of the supercell
     and N is the number of atoms.
    */
    /*
     double virial;
     virial = get_virial_AL(pos, L, N);
    */
    
    /*
     Function that calculates the forces on all atoms in units of [eV/Å]. the
     forces are stored in f which should be a matrix of size N x 3, where N is the
     number of atoms and column 1,2 and 3 correspond to the x,y and z component of
     the force resepctively . pos should be a matrix containing the positions of
     all the atoms, L is the length of the supercell and N is the number of atoms.
    */
    /*
     get_forces_AL(f,pos, L, N);
    */
    
    
    
    
    
    int Nc = 4, n_points=9;
    double pos[4*Nc*Nc*Nc][3];
    double a0[n_points], da0=(pow(68,1.0/3.0) - pow(64, 1.0/3.0)) / n_points;
    double E_pot[n_points];
    
    for (int i = 0; i < n_points; i++){
        a0[i] = 4 + da0 * i;
        init_fcc(pos, Nc, a0[i]);
        E_pot[i] = get_energy_AL(pos, Nc * a0[i], 256) / pow(Nc, 3);
    }

    FILE *fp = fopen("E_pot.csv", "w");
    fprintf(fp, "Volume(Å), Energy\n");
    for(int i = 0; i < n_points; ++i){
        fprintf(fp, "%f,%f\n", pow(a0[i],3), E_pot[i]);

    }
    fclose(fp);


    double a=4.0478;
    double m=0.002796;

    int timestep = 10000;
    double dt = 1e-2;

    double v[256][3];
    double pos_dev[256][3];
    double Ep[timestep+1];
    double Ek[timestep+1];
    
    double T[timestep+1];
    double P[timestep+1];
    
    double T_eq = 500.0 + 273.15;
    double P_eq = 1.0e-4;
    
    double tau_T = 2e2 * dt;
    double tau_P = 1e-1 * dt;
    
    double q1[timestep][2];
    double q2[timestep][2];
    double q3[timestep][2];

    init_fcc(pos_dev, Nc, a);

    srand((unsigned int)time(NULL));

    for (int j = 0; j < 256; j++){
        for(int d = 0; d < 3; d++){
            pos_dev[j][d] += -0.263107 + (double) rand() / ((double) RAND_MAX / (2 * 0.263107));
        }
    }
    velocity_verlet(timestep, 256, v, pos_dev, a, dt, m, T, P, Ep, Ek, tau_T, tau_P, T_eq, P_eq, q1, q2, q3);

    double T_avg = 0.0;
    double P_avg = 0.0;
    for (int i = 3000; i < timestep+1; i++){
        T_avg += T[i] / (timestep+1 - 1000.0);
        P_avg += P[i] / (timestep+1 - 1000.0);
    }

    double time_array[timestep+1];
    arange(time_array, 0, timestep, dt);

    write_to_file("Ep.csv", time_array, Ep, timestep);
    write_to_file("Ek.csv", time_array, Ek, timestep);
    
    write_to_file("T.csv", time_array, T, timestep);
    write_to_file("P.csv", time_array, P, timestep);
    
    trajectiores("q.csv", timestep, q1, q2, q3);

    printf("T_avg = %f T_eq = %f\n", T_avg, T_eq);
    printf("P_avg = %f P_eq = %f\n", P_avg, P_eq);

    return 0;
}

