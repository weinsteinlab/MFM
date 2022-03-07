#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define header_lines 11
#define title 2
#define trailor 5
#define PI 3.1415
#define eps_die 78.54
#define R 8.314
#define kbol 1.38
#define e_charge 1.6
#define eps_zero 8.85
#define T 300.0 
#define surf_z 40.0
#define a_lip 65.0
#define valency 4.0
#define f_memb 1648.457224946
#define f_syn 11398.40265387

main (int argc, char *argv[])

{
	int i, j, k, grid_points, nx, ny, nz, *nn_x_left, *nn_x_right, *nn_y_left, *nn_y_right, *nn_z_left, *nn_z_right;
	double *data, *diel, origin_x, origin_y, origin_z, counter, check_dist, gradient_square, radius, grad_x, grad_y, grad_z, hx, hy, hz, *mixing, *el_ene;
	double debye_inv_sq, bjerrum, integrand, const_term1, const_term2, pos_exp, neg_exp, grid[257][257], *fourth_term, fourth_term_memb;
	char headers[150];
	FILE *inptr, *paraptr, *outptr, *dieptr, *out2ptr, *in1ptr;

	
	grid_points = 0;

        if ((paraptr = fopen(argv[1], "r")) == NULL)
                printf ("\nERROR - Cannnot open the trajectory file\n");

        if ((inptr = fopen(argv[2], "r")) == NULL)
                printf ("\nERROR - Cannnot open the output file\n");

        if ((in1ptr = fopen(argv[3], "r")) == NULL)
                printf ("\nERROR - Cannnot open the output file\n");


	fscanf (paraptr, "%d\n", &nx);  		/*  grid points in x direction */		


	mixing = (double *) malloc(nx*(sizeof(double)));
	el_ene = (double *) malloc(nx*(sizeof(double)));
	fourth_term = (double *) malloc(nx*(sizeof(double)));

	for (i=0; i<nx; ++i) {
	        fscanf (inptr, "%lf %lf\n", &mixing[i], &fourth_term[i]);
	        fscanf (in1ptr, "%lf\n", &el_ene[i]);
		el_ene[i] = el_ene[i] / 2.494353;
		printf ("%d  %lf  %lf  %lf  %lf\n", i, mixing[i], el_ene[i], fourth_term[i], el_ene[i]+mixing[i]+fourth_term[i]-(f_memb/2.494353+f_syn/2.494353+fourth_term[0]));
	}
        


}	


