#include <stdio.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include "mt19937.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define NDIM 3   // The actual number of dimensions used in the simulation
#define MAXDIM 3 // Maximum number of dimensions, to make the output files work nicely
#define N 2000

#define NPATCHES 3 

/* Initialization variables */
const int    mc_steps      = 10000;
const int    output_steps  = 100;
const double density       = 0.8;
double       delta_r       = 0.1; // Initial step size
double       delta_a       = 0.1; // Initial angle change size 
const double beta          = 50;
const char*  init_filename = "../fcc/32.dat";
// Fix the diameter at 1.0, add a variable to control that
const double sigma = 1.0;
const double epsilon = 1.0;

// Kern Frenkel Patchy particle model parameters
const double patchdistance = 1.5; // Lambda, multiple of the diameter to which the path extends
const double coshalfangle = 0.9; // ~cos(pi/8) as a first try, some nice value for now


/* Simulation variables */
int n_particles = 0;
double r[N][MAXDIM];           // Positions of the particles
double directors[N][4];        // Quaternion rotations of the particles. The source of truth for orientation
// Note that the quaternion is stored as:
// q[0] * i + q[1] * j + q[2] * k + q[3]

double rot[N][MAXDIM][MAXDIM]; // Rotation matrices of the particles. Stored since they are used a lot
double box[MAXDIM];

double energy = 0.0;

// bond book keeping, tracking the single bond per patch KF condition
// Example particle_bonds[0][1][2] = 60 means that the second patch of particle 0 is bonded to the third patch of particle 60 
int particle_bonds[N][NPATCHES][NPATCHES] = {{{-1}}}; // for each particle, stores -1 if patch unoccupied or the pid of the partner particle at the N_partner_patch location


void   run_simulation(void);

double particle_energy(int pid);
int    move_particle(void);
int    rotate_particle(void);

void   matrix_vector_product(double res[MAXDIM], const double matrix[MAXDIM][MAXDIM], const double vector[MAXDIM]);
double dot_product(const double a[MAXDIM], const double b[MAXDIM]);

void   generate_random_unit_quaternion(double rv[4]);
void   set_rotation_matrix_from_director(int pid, const double d[4]);

void   read_data(void);
void   write_data(long int step);
void   set_density(void);
void   init_rotations(void);

FILE*  nice_fopen(const char* path, const char* mode);

//TODO: generate bonds

int main(int argc, char* argv[]){
	size_t seed = time(NULL);
    dsfmt_seed(seed);
	printf("Seed: %lu\n", seed);

	assert(delta_r > 0.0);
    assert(delta_a > 0.0);

	run_simulation();

    return 0;
}

void init_rotations(void) {
	// Initialise all rotation matrices as identity matrices
	int n,i,j;
	for (n = 0; n < n_particles; n++){
		for (i = 0; i < MAXDIM; i++){
			for (j = 0; j < MAXDIM; j++) rot[n][i][j] = (i==j) ? 1.0 : 0.0;
		}
		for (i = 0; i < 3; i++) directors[n][i] = 0.0;
		directors[n][3] = 1.0;
	}

}

void matrix_vector_product(double res[MAXDIM], const double matrix[MAXDIM][MAXDIM], const double vector[MAXDIM]){
	int i,k;

	for (i = 0; i < MAXDIM; i++){
		res[i] = 0.0;
		for (k = 0; k < MAXDIM; k++) res[i] += matrix[i][k] * vector[k];
	}
}

double dot_product(const double a[MAXDIM], const double b[MAXDIM]){
	double product = 0.0;

	for (int d = 0; d < MAXDIM; d++) product += a[d] * b[d];

	return product;
}

double particle_energy(int pid){
    double particle_energy = 0.0;
    int n, d, i, j;
	double dist2;
	// if in lambda, determine orientation, otherwise infinity or 0

	const double original_patch_directions[NPATCHES][NDIM] = {{1.0,0.0,0.0},
												 {-0.5,0.0,sin(2.0*M_PI/3.0)},
											  	 {-0.5,0.0,-sin(2.0*M_PI/3.0)}};

    for(n = 0; n < n_particles; ++n){
        if(n == pid) continue;
		double r_ij[NDIM]; // replaced min_d with r_ij so we can store the distances for later
		dist2 = 0.0;
		for(d = 0; d < NDIM; ++d){
			r_ij[d] = r[pid][d] - r[n][d];
			r_ij[d] -= (int)(2.0 * r_ij[d] / box[d]) * box[d];
			dist2 += r_ij[d] * r_ij[d]; // dist2 is the r_ij vector dot r_ij, so square diagonal distance
		}
		double dist = sqrt(dist2);
		if (dist < sigma){
			// This should never occur with the new check in move_particle, but just to be safe.
			particle_energy = INFINITY; // set energy ot infinity if particles overlap
			// printf("ERROR: PARTICLE OVERLAP\n");
			return particle_energy; 
		}
		if (dist > sigma * patchdistance) continue; // 0 energy for non interacting particles

		// Get the patch directors properly rotated
		double patch_directions[2][NPATCHES][MAXDIM];
		for (i = 0; i < NPATCHES; i++) matrix_vector_product(patch_directions[0][i], rot[pid], original_patch_directions[i]);
		for (i = 0; i < NPATCHES; i++) matrix_vector_product(patch_directions[1][i], rot[n],   original_patch_directions[i]);
		
		double r_hat[NDIM]; // unit vector pointing from n to pid
		double dot_products[2] = {0.0,0.0}; // 0 for pid, 1 for nth particle
		for(d = 0; d < NDIM; ++d){
			r_hat[d] = r_ij[d]/dist;
	//		dot_product[0] += rot[pid][0][d] * r_hat[d];
	//		dot_product[1] += rot[n][0][d] * (-r_hat[d]); // r_hat points from n to pid, so -r_hat points from pid to n
		}
		// Iterate over all patch combinations
		for (i = 0; i < NPATCHES; i++){
			for (j = 0; j < NPATCHES; j++){

				// r_hat is pointing from n to pid, so opposite signs are needed, my mistake
				dot_products[0] = -dot_product(patch_directions[0][i], r_hat);
				dot_products[1] = dot_product(patch_directions[1][j], r_hat);
				if (dot_products[0] > coshalfangle && dot_products[1] > coshalfangle) {
					particle_energy -= epsilon; // attractive

					// I must obey detailed balance when asigning particle pairs, 
					// in the same way I must allow for the desctruction of the particle pairs.
					// The formation and destruction (simpler) must be deterministic 

					// One option is to: always assign the bond to the lowest numbered particle, 
					// and remove the bond if either particle rotates such that the patches are no longer aligned, 
					// or if the particles move too far apart.

					// possibly then, I need to check ALL the moved particle bonds each time and assign energy based on bond number
					// particle_energy = -epsilon * n_bonds
					// thus new bonds - old bonds is deciding for dE

					// TO DO
					// 
					// 1. Match bonds in a deterministic manner upon initialisation (maybe implement if init condition)
					// 2. Upon each particled moved, destroy pid bonds and reform them deterministcally
					// 3. Pass on the new energy based on new bonds
				}
			}
		}

		// I believe that this is correct, but below is my first intuition which I think is wrong but I'm leaving it here for now

		// if(dist2 <= patchdistance * patchdistance){
		// 	// check alignment of patches w dot product
		// 	double dot_product = 0.0;
		// 	for (int i = 0; i < MAXDIM; i++){
		// 		for (int j = 0; j < MAXDIM; j++) dot_product += rot[pid][i][j] * rot[n][i][j];
		// 	}
		// 	if (dot_product > coshalfangle) particle_energy -= 1.0; // attractive
		// }
// This is how far I've come
    }
    return particle_energy;
}

int move_particle(void){
    int rpid = n_particles * dsfmt_genrand();

    double dE = -particle_energy(rpid); // of course we also assume this is not overlapping

    double old_pos[NDIM];
    int n,d;
    for(d = 0; d < NDIM; d++){
        old_pos[d] = r[rpid][d];
        r[rpid][d] += delta_r * (2.0 * dsfmt_genrand() - 1.0) + box[d];
        r[rpid][d] -= (int)(r[rpid][d] / box[d]) * box[d];
    }

	// // Check for particle overlap, and reject the move if particles overlap.
	// double abs_distance_squared, axis_distance;
	// for (n = 0; n < n_particles; n++){
	// 	if (n == rpid) continue; // Don't compare to the original position of the particle

	// 	abs_distance_squared = 0.0;

	// 	for (d = 0; d < NDIM; d++){
	// 		axis_distance = r[rpid][d] - r[n][d];
	// 		axis_distance -= box[d] * round(axis_distance / box[d]);

	// 		abs_distance_squared += axis_distance * axis_distance;
	// 	}

	// 	if (abs_distance_squared < sigma * sigma){
    // 		for(d = 0; d < NDIM; d++) r[rpid][d] = old_pos[d];
	// 		return 0;
	// 	}
	// }

	if (particle_energy(rpid) == INFINITY) {return 0;} // autoreject overlapping moves

    dE += particle_energy(rpid);
    if(dE < 0.0 || dsfmt_genrand() < exp(-beta * dE)){
        energy += dE;
        return 1;
    }

    for(d = 0; d < NDIM; d++) r[rpid][d] = old_pos[d];

    return 0;
}

void generate_random_unit_quaternion(double rv[4]){
// Generate a random unit quaternion, according to the method from Marsaglia, 1972
	int n;
	double S[2];  // Lengths
	for (n = 0; n < 2; n++){
		do{
			rv[2*n] = 2.0*dsfmt_genrand() - 1.0;
			rv[2*n+1] = 2.0*dsfmt_genrand() - 1.0;
			S[n] = rv[2*n]*rv[2*n] + rv[2*n+1]*rv[2*n+1];
		} while (S[n] >= 1);
	}
//	assert(S[0] < 1.0);
//	assert(S[1] < 1.0);

	rv[2] *= sqrt((1-S[0])/S[1]);
	rv[3] *= sqrt((1-S[0])/S[1]);


/** TESTING CODE STARTS HERE **/
//	double length = 0.0;
//	
//	for (n = 0; n < 4; n++){
//		printf("%lf\t", rv[n]);
//		length += rv[n]*rv[n];
//	}
//	printf("\tL: %e\n", sqrt(length)-1.0);
//
//	// This assertion doesn't always pass,
//	// Sometimes (due to rounding errors), the length is around 1 +/- 1e-16
//	// Which I call close enough
//	assert(rv[0]*rv[0] + rv[1]*rv[1] + rv[2]*rv[2] + rv[3]*rv[3] == 1.0);
}

int rotate_particle(void){
    int rpid = n_particles * dsfmt_genrand();

    double dE = -particle_energy(rpid);

    double new_director[4];
	int n;
	
/* Based on Frenkel & Smit, 1996:
Specify the orientation of the particle with the unit quaternion orientation
Generate a random unit quaternion rv (Vesely, 1982)
Get a rotated orientation vector: rv = orientation + delta_a * rv
Normalise rv again
Get the rotation matrix associated with rv
*/

// Generate the unit quaternion
	double rv[4]; // Random unit 4-vector
	generate_random_unit_quaternion(rv);

	double length = 0.0;
	for (n = 0; n < 4; n++){
		new_director[n] = rv[n] * delta_a + directors[rpid][n];
		length += new_director[n] * new_director[n];
	}
	length = sqrt(length);
	for (n = 0; n < 4; n++) new_director[n] /= length;

	set_rotation_matrix_from_director(rpid, new_director);

    dE += particle_energy(rpid);
    if(dE < 0.0 || dsfmt_genrand() < exp(-beta * dE)){
        energy += dE;
		for (n = 0; n < 4; n++) directors[rpid][n] = new_director[n];
        return 1;
    }
	
	set_rotation_matrix_from_director(rpid, directors[rpid]);

	return 0;
}

// Set the (globally stored) rotation matrix of particle pid according to the rotation quaternion d
void set_rotation_matrix_from_director(int pid, const double d[4]){
	rot[pid][0][0] = d[0]*d[0] - d[1]*d[1] - d[2]*d[2] + d[3]*d[3];
	rot[pid][0][1] = 2.0 * (d[0] * d[1] - d[2] * d[3]);
	rot[pid][0][2] = 2.0 * (d[2] * d[0] + d[1] * d[3]);

	rot[pid][1][0] = 2.0 * (d[0] * d[1] + d[2] * d[3]);
	rot[pid][1][1] = d[1]*d[1] - d[2]*d[2] - d[0]*d[0] + d[3]*d[3];
	rot[pid][1][2] = 2.0 * (d[1] * d[2] - d[0] * d[3]);

	rot[pid][2][0] = 2.0 * (d[2] * d[0] - d[1] * d[3]);
	rot[pid][2][1] = 2.0 * (d[1] * d[2] + d[0] * d[3]);
	rot[pid][2][2] = d[2]*d[2] - d[0]*d[0] - d[1]*d[1] + d[3]*d[3];
}

void run_simulation(){
    long int step;
	int n;
    read_data();
	init_rotations();

	if(n_particles == 0){
        printf("Error: Number of particles, n_particles = 0.\n");
        return;
    }

    set_density();

    //for(d = 0; d < NDIM; ++d) assert(r_cut <= 0.5 * box[d]);

    for(n = 0; n < n_particles; ++n) energy += particle_energy(n);
    energy *= 0.5;

    int accepted_mov = 0;
	int accepted_rot = 0;
	int total_mov = 0;
	int total_rot = 0;
    for(step = 0; step < mc_steps; step++){
        for(n = 0; n < n_particles; n++){
			// Probabilistically choose whether to move or rotate a particle, to obey detailed balance
			if (dsfmt_genrand() < 0.5){
				accepted_mov += move_particle();
				total_mov++;
			} else {
				accepted_rot += rotate_particle();
				total_rot++;
			}
        }

        if(step % output_steps == 0){
//            printf("Step %ld. Move acceptance: %lf\tRotation acceptance: %lf\n",
//                step,
//				(double)accepted_mov / total_mov,
//				(double)accepted_rot / total_rot
//            );
//			printf("\tdelta_r: %lf\tdelta_a: %lf\n", delta_r, delta_a);
			
			if ((double)accepted_mov / total_mov < 0.45) delta_r *= 0.9;
			if ((double)accepted_mov / total_mov > 0.55) delta_r /= 0.9;
			if ((double)accepted_rot / total_rot < 0.45) delta_a *= 0.9;
			if ((double)accepted_rot / total_rot > 0.55) delta_a /= 0.9;
            
    		accepted_mov = 0;
			accepted_rot = 0;
			total_mov = 0;
			total_rot = 0;
            write_data(step);
        }
    }
    write_data(step);

	printf("Done simulating!\007\nFinal energy: %lf\n", energy);
}

void read_data(void){
    FILE* fp = nice_fopen(init_filename, "r");
    int n, d;
    double dmin,dmax, diameter;
    fscanf(fp, "%d\n", &n_particles);
    for(d = 0; d < NDIM; ++d){
        fscanf(fp, "%lf %lf\n", &dmin, &dmax);
        box[d] = fabs(dmax-dmin);
    }

	if (NDIM == 2) box[2] = 0.0; // Ensure the 2d case is happy

    for(n = 0; n < n_particles; ++n){
        for(d = 0; d < NDIM; ++d) fscanf(fp, "%lf\t", &r[n][d]);
		if (NDIM == 2) r[n][2] = 0.0; // Ensure the 2d case is happy
        fscanf(fp, "%lf\n", &diameter);
    }
    fclose(fp);
}

void write_data(long int step){
    char buffer[128];
#if NDIM==2
	char extension[6] = ".patch";
#elif NDIM==3
	char extension[6] = "__.ptc"; // Not the nicest way, but prevents segfaults
#endif
    sprintf(buffer, "viscol/coords_step%07li%.6s", step, extension);
	FILE* fd = nice_fopen(buffer, "w");
	if (fd == NULL) return;
    int d, n,i,j;
    fprintf(fd, "%d\n", n_particles);
    for(d = 0; d < MAXDIM; d++) fprintf(fd, "%lf ", box[d]);
	fprintf(fd,"\n");
    for(n = 0; n < n_particles; n++){
		// <label> <x> <y> <z> <coreRadius> <cosHalfAngle> <capDiameter> <r00> ... <r22> [bondId ...]
		fprintf(fd, "a\t");
        for(d = 0; d < MAXDIM; d++) fprintf(fd, "%f\t", r[n][d]);
        fprintf(fd, "%lf\t%lf\t%lf\t", 0.5, coshalfangle, patchdistance);
		
		// Rotation matrix, see https://en.wikipedia.org/wiki/Rotation_matrix?useskin=vector#In_three_dimensions
		for (i = 0; i < MAXDIM; i++){
			for (j = 0; j < MAXDIM; j++) fprintf(fd, "%lf\t", rot[n][i][j]);
		}

		// Bonds
		// 3 patches -> 3 potential bonds
		// Eventually, this can become the ids of the particle it is bonded to I think
		// For now, -1 indicates no bond
		for (i = 0; i < NPATCHES-1; i++) fprintf(fd, "-1\t");
		fprintf(fd,"-1\n");

    }
    fclose(fd);
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

FILE* nice_fopen(const char* path, const char* mode){
	FILE* fd = fopen(path, mode);
	if (fd == NULL)	printf("Unable to open file \"%s\"\n", path);
	return fd;
}
