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

main (int argc, char *argv[])

{
	int i, j, k, grid_points, nx, ny, nz, *nn_x_left, *nn_x_right, *nn_y_left, *nn_y_right, *nn_z_left, *nn_z_right;
	double *data, *diel, origin_x, origin_y, origin_z, counter, check_dist, gradient_square, radius, grad_x, grad_y, grad_z, hx, hy, hz;
	double debye_inv_sq, bjerrum, integrand, const_term1, const_term2, pos_exp, neg_exp, grid[257][257];
	char headers[150];
	FILE *inptr, *paraptr, *outptr, *dieptr, *out2ptr, *out3ptr;

	
	grid_points = 0;

        if ((inptr = fopen(argv[1], "r")) == NULL)
                printf ("\nERROR - Cannnot open the trajectory file\n");

        if ((paraptr = fopen(argv[2], "r")) == NULL)
                printf ("\nERROR - Cannnot open the trajectory file\n");

        if ((outptr = fopen(argv[3], "w")) == NULL)
                printf ("\nERROR - Cannnot open the output file\n");



	fscanf (paraptr, "%d\n", &nx);  		/*  grid points in x direction */		
	fscanf (paraptr, "%d\n", &ny);			/*  grid points in y direction */
	fscanf (paraptr, "%d\n", &nz);  		/*  grid points in z direction */
	fscanf (paraptr, "%lf\n", &hx); 		/*  mesh size in x direction, in A */		
	fscanf (paraptr, "%lf\n", &hy);			/*  mesh size in y direction, in A */
	fscanf (paraptr, "%lf\n", &hz); 		/*  mesh size in z direction, in A */
	fscanf (paraptr, "%lf\n", &origin_x);  	/*  origin x coordinate */		
	fscanf (paraptr, "%lf\n", &origin_y);		/*  origin y coordinate */
	fscanf (paraptr, "%lf\n", &origin_z);  	/*  origin z coordinate */


	for (i=0; i<header_lines; ++i) {
	        fscanf (inptr, " %[^\n]", headers);
	}
        

	grid_points = nx*ny*nz;
	 
	data = (double *) malloc(grid_points*(sizeof(double)));

//	printf ("%d %d %d %d %lf %lf %lf\n", nx, ny, nz, grid_points, hx, hy, hz);	
//	printf ("%lf %lf %lf\n", origin_x, origin_y, origin_z);

	i = 0;
	
	do { 

              if (i==(grid_points-2))
//                   fscanf (inptr, "%lf\n", &data[i]);
                   fscanf (inptr, "%lf %lf\n", &data[i], &data[i+1]);
              else
                        fscanf (inptr, "%lf%lf%lf\n", &data[i], &data[i+1], &data[i+2]);
	

		i = i+3;
	
	} while (i < nx*ny*nz);


        for (i=0; i<nx; ++i) {
                for (j=0; j<ny; ++j) {
			grid[i][j] = 1.0;
		}
	}
	
	
			               fprintf (outptr, "%d\n", nx*ny);
/************************************** reading ends here *****************************************/	
		
        for (i=0; i<nx; ++i) {
                for (j=0; j<ny; ++j) {
                        for (k=0; k<nz; ++k) {
//                              if (data[i*ny*nz+j*nz+k] == 2.0) /* z direction profile */
                                if (origin_z+hz*k == surf_z) /* z direction profile */ {
			               fprintf (outptr, "%d  %d  %lf\n", i, j, -a_lip*data[i*ny*nz+j*nz+k]/valency/0.050375);
//			               fprintf (outptr, "%lf  %lf  %lf  %lf  %lf\n", origin_x+hx*i, origin_y+hy*j, origin_z+hz*k, data[i*ny*nz+j*nz+k], 0.0);
//				       if (-a_lip*data[i*ny*nz+j*nz+k]/valency/0.050375 > 1.08) printf ("%d  %d  %lf\n", i, j, -a_lip*data[i*ny*nz+j*nz+k]/valency/0.050375);
				       printf ("%d  %d  %lf\n", i, j, -a_lip*data[i*ny*nz+j*nz+k]/valency/0.050375);
				}

                        }
                }
        }

			               fprintf (outptr, "%d\n", 0);
				return 0;

					


                                       printf ("%d\n", nx*ny);
/************************************** reading ends here *****************************************/

        for (i=0; i<nx; ++i) {
                for (j=0; j<ny; ++j) {
                        for (k=0; k<nz; ++k) {
//                              if (data[i*ny*nz+j*nz+k] == 2.0) /* z direction profile */
                                if (origin_z+hz*k != surf_z && data[i*ny*nz+j*nz+k] != 0.0) /* z direction profile */ 
					grid[i][j] = 0.0;
//                                     fprintf (outptr, "%lf  %lf  %lf  %lf  %lf\n", origin_x+hx*i, origin_y+hy*j, origin_z+hz*k, data[i*ny*nz+j*nz+k], 0.0);

                        }
                }
        }


        for (i=0; i<nx; ++i) {
                for (j=0; j<ny; ++j) {
                        printf ("%d %d %lf\n", i, j, grid[i][j]);
                }
        }



                                       printf ("%d\n", 0);





	free(data);	
	fclose(inptr);
	fclose(paraptr);
	fclose(outptr);

}	


