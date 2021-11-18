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

void velocity_verlet(int n_timesteps, int nbr_atoms, double v[nbr_atoms][3], double pos[nbr_atoms][3], double L, double dt, double m, double *E_pot, double *E_kin);

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
    double dt = 0.0001;
    
    double v[256][3];
    double pos_dev[4*Nc*Nc*Nc][3];
    double Ep[timestep+1];
    double Ek[timestep+1];
    
    init_fcc(pos_dev, Nc, a);
    
    srand(time(NULL));
    
    for (int j = 0; j < 256; j++){
        for(int d=0; d<3; d++){
            pos_dev[j][d] += -0.263107 + (double) rand() / ((double) RAND_MAX / (2 * 0.263107));
        }
    }
    
    velocity_verlet(timestep, 256, v, pos_dev, a, dt, m, Ep, Ek);
    double Ek_avg=0;
    for (int i = 0; i < timestep+1; i++){
        if(fabs(Ek[i]) < 100.0){
            Ek_avg += Ek[i] / (timestep+1);
        }
        else{
            Ek[i]= 0.0/0.0;
        }
    }
    for (int i = 0; i < timestep+1; i++){
        if(Ep[i] > 0.0){
            Ep[i] = 0.0/0.0;
        }
    }

    double time_array[timestep+1];
    arange(time_array, 0, timestep, dt);

    //write_to_file("Ep.csv", time_array, Ep, timestep);
    //write_to_file("Ek.csv", time_array, Ek, timestep);
    
    double T;
    
    T = 2 / (3 * 256 * 8.6173 * pow(10, -5)) * (64 * Ek_avg);
    printf("%f\n", T);
    
    return 0;
}
