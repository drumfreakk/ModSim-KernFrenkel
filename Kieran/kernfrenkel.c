#include <stdio.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include <pthread.h>
#include <signal.h>
#include "mt19937.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define NDIM 3
#define N 257
#define SAVEABLESTEPS 100000

#define NVE
//#define NVT

#define WRITE_DATA
#define MSD

/* Initialization variables */
const double max_time      = 10;
const int    output_steps  = 10;
const double delta_t       = 0.0001; // 0.001 should really be enough for MD
const double kBT           = 1.0;
const double density       = 0.6;
const char*  init_filename = "../fcc/256.dat";

// Andersen thermostat
const double heat_coupling = 0.01;
const int    temp_steps    = 5;

// Particle properties
const double mass     = 1.0;
const double e        = 1.0;
const double diameter = 1.0;
const double friction = 1.0;

const double r_cut    = 2.5;

// Multithreading 
#define NTHREADS 8
int          fraction_nos[N];
int          fractionsize;
pthread_t    threads[NTHREADS];


/* Simulation variables */
int          n_particles;
double       e_cut;
double       r[N][NDIM][3];
	// r[n][d][0] is the previous position
	// r[n][d][1] is the current position
	// r[n][d][2] is the new position
double       box[NDIM];
long         step;

FILE* output_fd = NULL;

#ifdef MSD
int equilibration_time = 1.0;
int lag_jump = 100;
double max_lag_time = 9.0;

// Store the positions as if there were no periodic boundary conditions, to avoid artificially capping the MSD
double       positions_unbounded[N][SAVEABLESTEPS][NDIM];

FILE* MSD_fd = NULL;
#endif

// Function declarations
void   set_density(void);
void   initialize_velocities(void);

void   run_simulation(void);
void   move_particles(void);
void*  move_particle_group(void* group);
void   andersen_thermostat(void);

void   get_force(int particle, double force[NDIM]);
double get_potential_energy(void);
double get_kinetic_energy(void);

FILE*  nice_fopen(const char* path, const char* mode);
FILE*  fopen_with_parameters(const char* prefix, const char* mode);
void   read_data(void);
void   write_data(void);
void   close_fds(int sig);

#ifdef MSD
void   get_MSD(void);
#endif


int main(int argc, char* argv[]){
	time_t seed = time(NULL);
    dsfmt_seed(seed);
	printf("Seed: %lu\n", seed);

    assert(delta_t > 0.0);

#ifdef MSD
	assert(SAVEABLESTEPS >= max_time/delta_t); 
#endif

	signal(SIGINT, *close_fds);

	printf("Running simulation...\n");
	run_simulation();
	
#ifdef MSD
	printf("Calculating mean squared displacement...\007\n");
	get_MSD();
#endif

	printf("Done!\007\n");

    return 0;
}

#ifdef MSD
void get_MSD(void){
	int n,d,t,lag_steps;
	double autocorrelation, displacement;
	
	int equilibration_steps = equilibration_time/delta_t;
	int nsteps = max_time/delta_t;

	MSD_fd = fopen_with_parameters("MSD/", "w");
	if (MSD_fd == NULL) return;

	for (lag_steps = lag_jump; lag_steps < (int)(max_lag_time/delta_t); lag_steps+=lag_jump){
		// For all lag times, calculate the MSD
		autocorrelation = 0.0;

		for (t = equilibration_steps; t < nsteps - lag_steps; t++){
			for (n = 0; n < n_particles; n++){
				for (d = 0; d < NDIM; d++){
					displacement = positions_unbounded[n][t+lag_steps][d] - positions_unbounded[n][t][d];
					autocorrelation += displacement * displacement;
				}
			}
		}

		autocorrelation /= n_particles * (nsteps - lag_steps - equilibration_steps) * diameter * diameter;
		fprintf(MSD_fd, "%lf\t%lf\n", lag_steps*delta_t, autocorrelation);
		
		if (lag_steps % (lag_jump * 100) == 0) printf("Lag time: %lf\n", lag_steps*delta_t);
	}

	fclose(MSD_fd);
}
#endif

double get_potential_energy(void){
	double energy = 0.0;
	int n,k,d;

	double min_d, distsq, temp;

	for (n = 0; n < n_particles; n++){
		for (k = 0; k < n; k++) {		
			distsq = 0.0;
    	    
			// Get their nearest neighbours
			for (d = 0; d < NDIM; d++){
    	        min_d = r[k][d][1] - r[n][d][1];
    	        min_d -= (int)(2.0 * min_d / box[d]) * box[d];
    	        distsq += min_d * min_d;
    	    }

			// Calculate the potential energy
    	    if (distsq <= r_cut * r_cut){
				temp = powf(diameter,6.0)/(distsq * distsq * distsq);
				energy += 4.0 * e * temp * (temp - 1.0) - e_cut;
			}
		}
	}
	
	return energy;
}

double get_kinetic_energy(){
	double ekin = 0.0;

	for (int n = 0; n < n_particles; n++){
		for (int d = 0; d < NDIM; d++) ekin += (r[n][d][0] - r[n][d][1]) * (r[n][d][0] - r[n][d][1]);
	}

	return 0.5 * mass * ekin / (delta_t * delta_t);
}

void get_force(int particle, double force[NDIM]){
	int n,d;

	for (d = 0; d < NDIM; d++) force[d] = 0.0;

	double distsq, temp, f_d;
	double director[NDIM];

	for (n = 0; n < n_particles; n++){
		if (n == particle) continue;
	
		distsq = 0.0;
        
		// Get their nearest neighbours
		for (d = 0; d < NDIM; d++){
            director[d] = r[particle][d][1] - r[n][d][1];
            director[d] -= (int)(2.0 * director[d] / box[d]) * box[d];
            distsq += director[d] * director[d];
        }

		// Calculate the force from the potential energy
        if (distsq <= r_cut * r_cut){
			temp = powf(diameter,6.0)/(distsq * distsq * distsq);
			f_d = 24.0 * e * temp * (2.0*temp - 1.0); // The force multiplied by the distance

        	for (d = 0; d < NDIM; d++){
				// Add the component of f in the right direction to the force
				// In one step, normalize the director and divide the force by the distance
				// Avoiding a relatively costly square root
				force[d] += f_d * director[d] / distsq;
			}
		}
	}
}

void* move_particle_group(void* group){
	double displacement, force[NDIM];
	int n,d;

	for (n = *(int*)group * fractionsize; n < (*(int*)group + 1)*fractionsize; n++){
		get_force(n, force);
		
		for (d = 0; d < NDIM; d++){
			// Calculate the new position and save it in [2]
			// Calculate the new old nearest image and save it in [0]
			// Keep [1] constant for now
			r[n][d][2] = -r[n][d][0] + 2*r[n][d][1] + force[d]/mass * delta_t*delta_t;

			// Avoid large velocities due to periodic boundary conditions
			displacement = r[n][d][1] - r[n][d][2];

#ifdef MSD
			positions_unbounded[n][step+1][d] = positions_unbounded[n][step][d] + displacement;
#endif

			if (r[n][d][2] < 0.0)    r[n][d][2] += box[d];
			if (r[n][d][2] > box[d]) r[n][d][2] -= box[d];
			
			// Set the old position, to avoid having to store it elsewhere
			r[n][d][0] = r[n][d][2] + displacement;
		}
	}
	return NULL;
}

void move_particles(void){
	int f,n,d;

	// Update all particle positions in NTHREADS threads, then wait for all threads to finish
	// There must be prettier ways to pass what particles to update to each thread, 
	// but I couldn't figure out a way to get the pointers working, so this will do.
	for (f = 0; f < NTHREADS; f++) pthread_create(&threads[f], NULL, move_particle_group, &fraction_nos[f]);
	for (f = 0; f < NTHREADS; f++) pthread_join(threads[f], NULL);
	
	// Move the new position to be the current position, now all forces have been calculated
	for (n = 0; n < n_particles; n++) for (d = 0; d < NDIM; d++) r[n][d][1] = r[n][d][2];
}

void andersen_thermostat(void){
#if NDIM != 3
	printf("Wrong NDIM\n");
	return;
#endif

	int n,d;

	double new_velocity[NDIM];
	double new_vel_sq;
	double U,V;

	double prefactor = pow(mass/(2.0*M_PI*kBT), 3.0/2.0) * 4.0 * M_PI;

	for (n = 0; n < n_particles; n++){
		if (dsfmt_genrand() >= temp_steps*delta_t*heat_coupling) continue;
		
		do {
			U = dsfmt_genrand();
			V = dsfmt_genrand();
			new_velocity[0] = sqrt(-2.0*log(U)) * cos(2.0*M_PI*V);
			new_velocity[1] = sqrt(-2.0*log(U)) * sin(2.0*M_PI*V);
			
			U = dsfmt_genrand();
			V = dsfmt_genrand();
			new_velocity[2] = sqrt(-2.0*log(U)) * cos(2.0*M_PI*V);
	
			new_vel_sq = 0.0;
			for (d = 0; d < NDIM; d++) new_vel_sq += new_velocity[d]*new_velocity[d];
		} while (dsfmt_genrand() >= prefactor * new_vel_sq * exp(-mass * new_vel_sq / (2.0 * kBT)));
		
		for (d = 0; d < NDIM; d++) r[n][d][0] = r[n][d][1] - new_velocity[d] * delta_t;
	}
}

void initialize_velocities(){
	int n,d;
	
	double v[N][NDIM];
	double vsq, scale;
	double avg_v   = 0.0;
	double avg_vsq = 0.0;

	for (n = 0; n < n_particles; n++){
		vsq = 0.0;
		for (d = 0; d < NDIM; d++){
			v[n][d] = 2.0 * dsfmt_genrand() - 1.0;
			vsq += v[n][d] * v[n][d];
		}
		avg_v   += sqrt(vsq);
		avg_vsq += vsq;
	}

	avg_v   /= n_particles;
	avg_vsq /= n_particles;
	scale    = sqrt(3*kBT/(mass*avg_vsq));

	// Set the previous position based on the scaled & shifted velocity
	for (n = 0; n < n_particles; n++){
		for (d = 0; d < NDIM; d++) r[n][d][0] = r[n][d][1] - delta_t * (v[n][d] - avg_v) * scale;
	}
}

void run_simulation(void){
    int n,d;
    double volume, ekin, epot;
	double avg_ekin = 0.0;
	double avg_epot = 0.0;
	

	read_data();

	if(n_particles == 0){
        printf("Error: Number of particles, n_particles = 0.\n");
    	return;
	}

    e_cut = 4.0 * e * (pow(diameter / r_cut, 12.0) - pow(diameter / r_cut, 6.0));

    set_density();

    for(d = 0; d < NDIM; d++) assert(r_cut <= 0.5 * box[d]);

	volume = 1.0;
    for(d = 0; d < NDIM; d++) volume *= box[d];
	
	initialize_velocities();

	fractionsize = n_particles/NTHREADS;
	for (n = 0; n < NTHREADS; n++) fraction_nos[n] = n;

	output_fd = fopen_with_parameters("output/", "w");
	if (output_fd == NULL) return;

	for (step = 0; step < (int)(max_time / delta_t); step++){
#ifdef NVT
		if (step % temp_steps == 0) andersen_thermostat();
#endif
		
		move_particles();

        if(step % output_steps == 0){
#ifdef WRITE_DATA
			write_data();
#endif

			ekin = get_kinetic_energy()/e;
			epot = get_potential_energy()/e;
			avg_ekin += ekin;
			avg_epot += epot;
			//fprintf(output_fd, "%lf\t%lf\t%lf\t%lf\n", step * delta_t, ekin+epot, ekin, epot);
			fprintf(output_fd, "%lf\t%lf\n", step * delta_t, ekin+epot);
    	}
	}
	
	fclose(output_fd);

	avg_ekin /= max_time / ((double)(output_steps) * delta_t);
	avg_epot /= max_time / ((double)(output_steps) * delta_t);

	printf("Average kinetic energy: %lf\t Average potential energy: %lf\n", avg_ekin, avg_epot);

	return;
}

FILE* fopen_with_parameters(const char* prefix, const char* mode){
	char buffer[128];

#ifdef NVE
	char type[3] = "NVE";
#endif
#ifdef NVT
	char type[3] = "NVT";
#endif

    sprintf(buffer, "%s%.3s_kBT_%.1lf_rho_%.2lf_dt_%lf.tsv",
	        prefix, type, kBT, density, delta_t);

	return nice_fopen(buffer, mode);
}

void read_data(void){
    FILE* fp = nice_fopen(init_filename, "r");
    if (fp == NULL) return;
	int n, d;
    double dmin,dmax;
	double temp;
    fscanf(fp, "%d\n", &n_particles);
    for(d = 0; d < NDIM; d++){
        fscanf(fp, "%lf %lf\n", &dmin, &dmax);
        box[d] = fabs(dmax-dmin);
    }
    for(n = 0; n < n_particles; n++){
        for(d = 0; d < NDIM; ++d) fscanf(fp, "%lf\t", &r[n][d][1]);
        fscanf(fp, "%lf\n",&temp);
    }
    fclose(fp);
}

void write_data(){
    char buffer[128];
    sprintf(buffer, "viscol/coords_step%07li.dat", step);
    FILE* fp = nice_fopen(buffer, "w");
	if (fp == NULL) return;
    int d, n;
    fprintf(fp, "%d\n", n_particles);
    for(d = 0; d < NDIM; d++){
        fprintf(fp, "%lf %lf\n",0.0,box[d]);
    }
    for(n = 0; n < n_particles; n++){
        for(d = 0; d < NDIM; ++d) fprintf(fp, "%f\t", r[n][d][1]);
        fprintf(fp, "%lf\n", 1.0);
    }
    fclose(fp);
}

void set_density(void){
    double volume = 1.0;
    int d, n;
    for(d = 0; d < NDIM; d++) volume *= box[d];

    double target_volume = n_particles / density;
    double scale_factor = pow(target_volume / volume, 1.0 / NDIM);

    for(n = 0; n < n_particles; n++){
        for(d = 0; d < NDIM; d++) r[n][d][1] *= scale_factor;
    }
    for(d = 0; d < NDIM; d++) box[d] *= scale_factor;
}

void close_fds(int sig){
	// Ensure the output files are still properly written to if the program is terminated
	if (output_fd != NULL) fclose(output_fd);

#ifdef MSD
	if (MSD_fd != NULL) fclose(MSD_fd);
#endif

	raise(SIGTERM);
}

FILE* nice_fopen(const char* path, const char* mode){
	FILE* fd = fopen(path, mode);
	if (fd == NULL)	printf("Unable to open file \"%s\"\n", path);
	return fd;
}
