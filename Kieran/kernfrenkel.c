#include <stdio.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include "mt19937.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define NDIM 3
#define N 2000

/* Initialization variables */
const int    mc_steps      = 10000;
const int    output_steps  = 100;
const int    eq_steps      = 2000; // Equilibration steps
const double density       = 0.6;
const double delta         = 0.1; // Step size
const double beta          = 0.5;
const char*  init_filename = "../fcc/864.dat";

// Kern Frenkel Patchy particle model parameters
const double patchdistance = 1.5; // Lambda, multiple of the diameter to which the path extends
const double coshalfangle = 0.9; // ~cos(pi/8) as a first try, some nice value for now

/* Simulation variables */
int n_particles = 0;
double e_cut;
double r[N][NDIM];
double box[NDIM];

double energy = 0.0;


void run_simulation(void);

double particle_energy(int);
int    move_particle(void);

void   read_data(void);
void   write_data(int step);
void   set_density(void);

int main(int argc, char* argv[]){
    assert(delta > 0.0);
    
	size_t seed = time(NULL);
    dsfmt_seed(seed);
	printf("Seed: %lu\n", seed);

	run_simulation();

    return 0;
}

double particle_energy(int pid){
    double particle_energy = 0.0;
    int n, d;
	double dist2, min_d, temp;
    for(n = 0; n < n_particles; ++n){
        if(n == pid) continue;

//TODO Properly calculate the energy
//        dist2 = 0.0;
//        for(d = 0; d < NDIM; ++d){
//            min_d = r[pid][d] - r[n][d];
//            min_d -= (int)(2.0 * min_d / box[d]) * box[d];
//            dist2 += min_d * min_d;
//        }
//
//        if(dist2 <= r_cut * r_cut){
//            temp = 1.0 / (dist2 * dist2 * dist2);
//            particle_energy += 4.0 * temp * (temp - 1.0) - e_cut;
//        }
    }

    return particle_energy;
}

int move_particle(void){
    int rpid = n_particles * dsfmt_genrand();

    double dE = -particle_energy(rpid);

    double old_pos[NDIM];
    int d;
    for(d = 0; d < NDIM; ++d){
        old_pos[d] = r[rpid][d];
        r[rpid][d] += delta * (2.0 * dsfmt_genrand() - 1.0) + box[d];
        r[rpid][d] -= (int)(r[rpid][d] / box[d]) * box[d];
    }

    double dE += particle_energy(rpid);
    if(dE < 0.0 || dsfmt_genrand() < exp(-beta * dE)){
        energy += dE;
        return 1;
    }

    for(d = 0; d < NDIM; ++d) r[rpid][d] = old_pos[d];

    return 0;
}

void run_simulation(){
    int step, n, d;
    read_data();

	if(n_particles == 0){
        printf("Error: Number of particles, n_particles = 0.\n");
        return mean_values;
    }

    set_density();

    //for(d = 0; d < NDIM; ++d) assert(r_cut <= 0.5 * box[d]);

    for(n = 0; n < n_particles; ++n) energy += particle_energy(n);
    energy *= 0.5;

    int accepted = 0;
    for(step = -eq_steps; step < mc_steps; ++step){
        for(n = 0; n < n_particles; ++n){
            accepted += move_particle();
        }

        if(step % output_steps == 0){
//            printf("Step %d. Move acceptance: %f.\n",
//                step, (double)accepted / (n_particles * output_steps)
//            );
			
			if ((double)accepted / (n_particles * output_steps) < 0.4) delta *= 0.9;
			if ((double)accepted / (n_particles * output_steps) > 0.6) delta /= 0.9;
            
			accepted = 0;
            write_data(step);
        }
    }
}

void read_data(void){
    FILE* fp = fopen(init_filename, "r");
    int n, d;
    double dmin,dmax;
    fscanf(fp, "%d\n", &n_particles);
    for(d = 0; d < NDIM; ++d){
        fscanf(fp, "%lf %lf\n", &dmin, &dmax);
        box[d] = fabs(dmax-dmin);
    }
    for(n = 0; n < n_particles; ++n){
        for(d = 0; d < NDIM; ++d) fscanf(fp, "%lf\t", &r[n][d]);
        double diameter;
        fscanf(fp, "%lf\n", &diameter);
    }
    fclose(fp);
}

void write_data(int step){
    char buffer[128];
#if NDIM==2
	char extension[6] = ".patch";
#elif NDIM==3
	char extension[6] = "__.ptc"; // Not the nicest way, but prevents segfaults
#endif
    sprintf(buffer, "viscol/coords_step%07li%.6s", step, extension);
	FILE* fd = nice_fopen(buffer, "w");
	if (fd == NULL) return;
    int d, n;
    fprintf(fd, "%d\n", n_particles);
    for(d = 0; d < NDIM; d++) fprintf(fd, "%lf ", box[d]);
	fprintf(fd,"\n");
    for(n = 0; n < n_particles; n++){
		// <label> <x> <y> <z> <coreRadius> <cosHalfAngle> <capDiameter> <r00> ... <r22> [bondId ...]
		fprintf(fd, "a\t");
        for(d = 0; d < NDIM; d++) fprintf(fd, "%f\t", r[n][d][1]);
        fprintf(fd, "%lf\t%lf\t%lf\t", diameter/2.0, coshalfangle, diameter * patchdistance);
		
		// Rotation matrix, see https://en.wikipedia.org/wiki/Rotation_matrix?useskin=vector#In_three_dimensions
		// until we figure out a nice way to store rotation, use the identity matrix (== no rotation)
		fprintf(fd, "1.0\t0.0\t0.0\t0.0\t1.0\t0.0\t0.0\t0.0\t1.0\t");

		// Bonds
		// 3 patches -> 3 potential bonds
		// Eventually, this can become the ids of the particle it is bonded to I think
		// For now, -1 indicates no bond
		fprintf(fd, "-1\t-1\t-1\n");

    }
    fclose(fd);
	printf("got here\n");
}

void set_density(void){
    double volume = 1.0;
    int d, n;
    for(d = 0; d < NDIM; ++d) volume *= box[d];

    double target_volume = n_particles / density;
    double scale_factor = pow(target_volume / volume, 1.0 / NDIM);

    for(n = 0; n < n_particles; ++n){
        for(d = 0; d < NDIM; ++d) r[n][d] *= scale_factor;
    }
    for(d = 0; d < NDIM; ++d) box[d] *= scale_factor;
}

