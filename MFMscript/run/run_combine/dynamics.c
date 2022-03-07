#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#define header_lines 11
#define trailor 5
#define coor_number 4.0

void Randomize();
double Random(double);
double GRandom(double, double);
void ReadMap(FILE *inptr, double grid_points, double *data);
void WriteMap(FILE *outptr, double grid_points, double *data);

char headers[11][150], tails[5][150];


main (int argc, char *argv[])

{

	int i, j, k, nx, ny, nz, grid_points, **nn_x_left, **nn_x_right, **nn_y_up, **nn_y_down, tt;
	int **nnn_left_up_x, **nnn_left_down_x, **nnn_right_up_x, **nnn_right_down_x;
	int **nnn_left_up_y, **nnn_left_down_y, **nnn_right_up_y, **nnn_right_down_y;
	double noise, **niu1, **niu2, **order_parameter, *dataPotential, *dataCharge, *dataKappa, **SurfPotential, w, v, valency;
	double pos_exp, neg_exp, Iv_integrand, Ip_integrand, Is_integrand, const_rand, dt, laplacian_term, sigma_ion;
	double hx, hy, hz, origin_x, origin_y, origin_z, check_sum, average_op, average_sigma, surf_z, **dFdEta, up, down, sum_op, area_per_lipid, dFdEta_gradient_left, dFdEta_gradient_down, density_left, density_down, **flow_left, **flow_down;
	FILE *paraptr, *potptr, *chargeptr, *out3Dptr, *outptr, *kappaptr;

	Randomize();
	check_sum = 0;
	w = 1.0;
	v = 1.0-w;

    	if ((paraptr = fopen(argv[1], "r")) == NULL)
                printf ("\nERROR - Cannnot open parameter file\n");

        if ((potptr = fopen(argv[2], "r")) == NULL)
                printf ("\nERROR - Cannnot open potential map file\n");

        if ((chargeptr = fopen(argv[3], "r")) == NULL)
                printf ("\nERROR - Cannnot open charge map file\n");

        if ((kappaptr = fopen(argv[4], "r")) == NULL)
                printf ("\nERROR - Cannnot open kappa map file\n");

        if ((outptr = fopen(argv[5], "w")) == NULL)
                printf ("\nERROR - Cannnot open the output file\n");


/****** READ PARAMETER FILE ***/

        fscanf (paraptr, "%d\n", &nx);                  /*  grid points in x direction */
        fscanf (paraptr, "%d\n", &ny);                  /*  grid points in y direction */
        fscanf (paraptr, "%d\n", &nz);                  /*  grid points in z direction */
        fscanf (paraptr, "%lf\n", &hx);                 /*  mesh size in x direction, in A */
        fscanf (paraptr, "%lf\n", &hy);                 /*  mesh size in y direction, in A */
        fscanf (paraptr, "%lf\n", &hz);                 /*  mesh size in z direction, in A */
        fscanf (paraptr, "%lf\n", &origin_x); 		/*  origin x coordinate */
        fscanf (paraptr, "%lf\n", &origin_y);           /*  origin y coordinate */
        fscanf (paraptr, "%lf\n", &origin_z);   	/*  origin z coordinate */
        fscanf (paraptr, "%lf\n", &surf_z);   		/*  z coordinate of the surface */
        fscanf (paraptr, "%lf\n", &average_sigma);      /*  average charge density on the surface */
        fscanf (paraptr, "%lf\n", &dt);   		/*  dimensioneless timestep for the dynamics */
        fscanf (paraptr, "%lf\n", &valency);   		/*  valency of the lipid */
        fscanf (paraptr, "%lf\n", &area_per_lipid);     /*  area per lipid */

	const_rand = sqrt(2*dt);			 /*  constant term up front of the noise term in CH equation*/
	const_rand = 0.0;
	noise = 0.0;
	average_op = fabs(area_per_lipid*average_sigma/valency); /* this is average fraction of chargd lipids */
        grid_points = nx*ny*nz;

        dataPotential = (double *) malloc(grid_points*(sizeof(double)));
        dataCharge = (double *) malloc(grid_points*(sizeof(double)));
        dataKappa = (double *) malloc(grid_points*(sizeof(double)));

/*** READ APBS MAPS ****/

	ReadMap(potptr, grid_points, dataPotential);
	ReadMap(chargeptr, grid_points, dataCharge);


	order_parameter = (double **) malloc(nx * sizeof(double *));
        order_parameter[0] = malloc(nx * ny * sizeof(double));
        for(i = 1; i < nx; ++i)
                order_parameter[i] = order_parameter[0] + i * ny;

        SurfPotential = (double **) malloc(nx * sizeof(double *));
        SurfPotential[0] = malloc(nx * ny * sizeof(double));
        for(i = 1; i < nx; ++i)
                SurfPotential[i] = SurfPotential[0] + i * ny;

        niu1 = (double **) malloc(nx * sizeof(double *));
        niu1[0] = malloc(nx * ny * sizeof(double));
        for(i = 1; i < nx; ++i)
                niu1[i] = niu1[0] + i * ny;

        niu2 = (double **) malloc(nx * sizeof(double *));
        niu2[0] = malloc(nx * ny * sizeof(double));
        for(i = 1; i < nx; ++i)
                niu2[i] = niu2[0] + i * ny;

        dFdEta = (double **) malloc(nx * sizeof(double *));
        dFdEta[0] = malloc(nx * ny * sizeof(double));
        for(i = 1; i < nx; ++i)
                dFdEta[i] = dFdEta[0] + i * ny;

	flow_left = (double **) malloc(nx * sizeof(double *));
        flow_left[0] = malloc(nx * ny * sizeof(double));
        for(i = 1; i < nx; ++i)
                flow_left[i] = flow_left[0] + i * ny;


	flow_down = (double **) malloc(nx * sizeof(double *));
        flow_down[0] = malloc(nx * ny * sizeof(double));
        for(i = 1; i < nx; ++i)
                flow_down[i] = flow_down[0] + i * ny;

        nn_x_left = (int **) malloc(nx * sizeof(int *));
        nn_x_left[0] = malloc(nx * ny * sizeof(int));
        for(i = 1; i < nx; ++i)
                nn_x_left[i] = nn_x_left[0] + i * ny;

        nn_x_right = (int **) malloc(nx * sizeof(int *));
        nn_x_right[0] = malloc(nx * ny * sizeof(int));
        for(i = 1; i < nx; ++i)
                nn_x_right[i] = nn_x_right[0] + i * ny;

        nn_y_up = (int **) malloc(nx * sizeof(int *));
        nn_y_up[0] = malloc(nx * ny * sizeof(int));
        for(i = 1; i < nx; ++i)
                nn_y_up[i] = nn_y_up[0] + i * ny;
	
        nn_y_down = (int **) malloc(nx * sizeof(int *));
        nn_y_down[0] = malloc(nx * ny * sizeof(int));
        for(i = 1; i < nx; ++i)
                nn_y_down[i] = nn_y_down[0] + i * ny;
                

/*  Identifying charged surface points and remembering potential there and also calculating local lipid composition*/
/*  which I call order parameter */

        for (i=0;i<nx;++i) {
                for (j=0;j<ny;++j) {
                        for (k=0;k<nz;++k) {
                                if ((origin_z+hz*k) == surf_z) {
                                        SurfPotential[i][j] = dataPotential[i*ny*nz+j*nz+k];
                                        order_parameter[i][j] = fabs(dataCharge[i*ny*nz+j*nz+k]*area_per_lipid/valency);
                                }
                        }
                }
        }



	for (i=0;i<nx;++i) {		
		for (j=0;j<ny;++j) {
			if (i==0) 				/* This if-else structure identifies 4 nearest neighbors, */		
				nn_x_left[i][j]=nx-1;		/* nn_x_left, nn_x_right, nn_y_up, nn_y_down, */
			else					/* on 2-D square lattice by taking care of boundary conditions */
				nn_x_left[i][j]=i-1;
			if (i==nx-1) 
				nn_x_right[i][j]=0;
			else
				nn_x_right[i][j]=i+1;
			if (j==0) 
				nn_y_up[i][j]=ny-1;
			else
				nn_y_up[i][j]=j-1;
			if (j==ny-1) 
				nn_y_down[i][j]=0;
			else
				nn_y_down[i][j]=j+1;
 
//			niu1[i][j] = GRandom(0,1);		/* Gaussian random number*/
//			niu2[i][j] = GRandom(0,1);
			up = (order_parameter[i][j])*(1-average_op);
			down = (average_op)*(1-order_parameter[i][j]);
			dFdEta[i][j] = log(up/down) - valency*SurfPotential[i][j];  /* calculating gradients of F with respect to eta */
		}
	}


/* Dynamic update of the system */
        for (i=0;i<nx;++i) {
                for (j=0;j<ny;++j) {
			dFdEta_gradient_left = dFdEta[i][j] - dFdEta[nn_x_left[i][j]][j];
			dFdEta_gradient_down = dFdEta[i][j] - dFdEta[i][nn_y_down[i][j]];
			density_left = (order_parameter[nn_x_left[i][j]][j] + order_parameter[i][j]) *0.5;
			density_down = (order_parameter[i][nn_y_down[i][j]] + order_parameter[i][j]) *0.5;
			flow_left[i][j] = dFdEta_gradient_left * density_left;
			flow_down[i][j] = dFdEta_gradient_down * density_down;
		}
	}

	for (i=0;i<nx;++i) {
		for (j=0;j<ny;++j) {
//			noise = niu1[nn_x_right[i][j]][j] - niu1[i][j] + niu2[i][nn_y_down[i][j]]- niu2[i][j];
			laplacian_term = (- flow_left[i][j] + flow_left[nn_x_right[i][j]][j] - flow_down[i][j] + flow_down[i][nn_y_up[i][j]])/2;
			order_parameter[i][j] = order_parameter[i][j] + dt* laplacian_term + const_rand * noise;
		}					 	
	}



/* Creating updated charge map on the surface and then writing the entire charge map for APBS */

        for (i=0;i<nx;++i) {
                for (j=0;j<ny;++j) {
                        for (k=0;k<nz;++k) {
                                if ((origin_z+hz*k) == surf_z) {
					if (average_sigma < 0) 
						dataCharge[i*ny*nz+j*nz+k] = -order_parameter[i][j]*valency/area_per_lipid;
					else
						dataCharge[i*ny*nz+j*nz+k] = order_parameter[i][j]*valency/area_per_lipid;
                                }
                        }
                }
        }


        WriteMap(outptr, grid_points, dataCharge);

/* Releasing all the memory*/

	free(dataPotential);
	free(dataCharge);
	free(dataKappa);
	free(order_parameter);
	free(flow_left);
	free(flow_down);
	free(SurfPotential);
	free(dFdEta);
	free(niu1);
	free(niu2);
	free(nn_x_right);
	free(nn_x_left);
	free(nn_y_up);
	free(nn_y_down);
	fclose(paraptr);
	fclose(potptr);
	fclose(chargeptr);
	fclose(kappaptr);
	fclose(outptr);

}

/************************************************************************************************************/

double Random(double range)
{
  return (random()*range)/2147483647L;
}

/************************************************************************************************************/

void Randomize()
{
  srandom(time(NULL));
}

/************************************************************************************************************/

double GRandom(double mean, double sigma)
{
  double v1 = Random(2.0) - 1.0;
  double v2 = Random(2.0) - 1.0;
  double rsq = v1*v1 + v2*v2;

  while (rsq >= 1.0 || rsq == 0.0)
    {
      v1 = Random(2.0) - 1.0;
      v2 = Random(2.0) - 1.0;
      rsq = v1*v1 + v2*v2;
    }

  return mean+v1*sigma*sqrt(-2.0*log(rsq)/rsq);
}

/*** TO READ APBS MAPS ********************/


void ReadMap (FILE *inptr, double grid_points, double *data) {

        int i;

        for (i=0; i<header_lines; ++i) {
                fscanf (inptr, " %[^\n]", headers[i]);
        }


        i = 0;

        do {

              if (i==(grid_points-2))
                   fscanf (inptr, "%lf %lf\n", &data[i], &data[i+1]);
              else
                   fscanf (inptr, "%lf%lf%lf\n", &data[i], &data[i+1], &data[i+2]);


                i = i+3;

        } while (i < grid_points);


        for (i=0; i<trailor; ++i) {
                fscanf (inptr, " %[^\n]", tails[i]);
        }

}

/*** TO WRITE APBS MAPS ********************/


void WriteMap (FILE *outptr, double grid_points, double *data) {

        int i;

        for (i=0; i<header_lines; ++i) {
                fprintf (outptr, "%s\n", headers[i]);
        }


        i = 0;

        do {

              if (i==(grid_points-2))
                   fprintf (outptr, "%16.14f %16.14f\n", data[i], data[i+1]);
              else
                   fprintf (outptr, "%16.14f %16.14f %16.14f\n", data[i], data[i+1], data[i+2]);


                i = i+3;

        } while (i < grid_points);


        for (i=0; i<trailor; ++i) {
                fprintf (outptr, "%s\n", tails[i]);
        }

}

