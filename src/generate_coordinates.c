#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define D 1.0


// N is the amount of spheres per line on an axis
void cubic(double l, int N, FILE* output){

    fprintf(output, "%i\n", N*N*N);
    fprintf(output, "%lf\t%lf\n", 0.0, N * l);
    fprintf(output, "%lf\t%lf\n", 0.0, N * l);
    fprintf(output, "%lf\t%lf\n", 0.0, N * l);

    for (int x = 0; x < N; x++){
        for (int y = 0; y < N; y++){
            for (int z = 0; z < N; z++){
                fprintf(output, "%lf\t%lf\t%lf\t%lf\n", x*l, y*l, z*l, D);
            }
        }
    }
}

// N is the number of unit cells of 4 spheres per line on each axis
void fcc(double l, int N, FILE* output){

	double size = l * 2.0 * N;
   
	fprintf(output, "%i\n", 4*N*N*N);

	fprintf(output, "%lf\t%lf\n", 0.0, size);
    fprintf(output, "%lf\t%lf\n", 0.0, size);
    fprintf(output, "%lf\t%lf\n", 0.0, size);

	for (int x = 0; x < 2*N; x++){
		for (int y = 0; y < 2*N; y++){
			for (int z = 0; z < 2*N; z++){
				// Fancy math to automagically generate an FCC unit cell
				if ((x+y+z) % 2 == 0){
					fprintf(output, "%lf\t%lf\t%lf\t%lf\n", x * l, y * l, z * l, D);
				}
			}
		}
	}
}

void bcc(double l, int N, FILE* output){
    int count = 0;
	for (int x = 0; x < N; x++){
        for (int y = 0; y < N; y++){
            for (int z = 0; z < N; z++){
                count+=2;
			}
        }
    }

    fprintf(output, "%i\n", count);
    fprintf(output, "%lf\t%lf\n", 0.0, N * 2*l);
    fprintf(output, "%lf\t%lf\n", 0.0, N * 2*l);
    fprintf(output, "%lf\t%lf\n", 0.0, N * 2*l);

    for (int x = 0; x < N; x++){
        for (int y = 0; y < N; y++){
            for (int z = 0; z < N; z++){
                fprintf(output, "%lf\t%lf\t%lf\t%lf\n", x*2*l, y*2*l, z*2*l, D);
                //if (x != N-1 && y != N-1 && z != N-1){
					fprintf(output, "%lf\t%lf\t%lf\t%lf\n", x*2*l+l, y*2*l+l, z*2*l+l, D);
            	//}
			}
        }
    }
}

int main(){
    FILE* output;
	output = fopen("../bcc/bcc.dat", "w");
	bcc(1.0, 10, output);

//	output = fopen("cubic.xyz", "w");
//    cubic(1.0, 5, output);    
//    
//    output = fopen("fcc_6.xyz", "w");
//    fcc(1.0/sqrt(2.0), 4, output);
    
	fclose(output);

    return 0;
}
