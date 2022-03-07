/*	Integration of tabulated function in n dimensions using B-splines */

#include <stdio.h>
#include <math.h>
#include <malloc.h>
#define PI 3.1415
#define eps_die 78.54
#define R 8.314
#define kbol 1.38
#define e_charge 1.6
#define eps_zero 8.85
#define header_lines 11
#define trailor 5


int bspint(int n, double x[], double f[], int k, double *a, int la, double c[],
	double *xf, int *no, int *iflg, int inc[]);
int gaubnd(int n, int kb, int num, double *a, double *x, double *det,
		int *idet, int inc[], int lj, int *iflg);
int bsplin(double *x, int nx, int k, double xb, int nderiv, double b[],
		double db[], double ddb[], int *left);
double bspqd(int n, double *x, int k, double wt[], double xl, double xu, int *ier);
double bspqdn(int n, int nk[], double *x, int nxd, int k, double *wt,
	double xl[], double xu[], int *ier);
int bspintn(int n, int nk[], double *x, int nxd, double *f, int k, double *ax,
	double *c, double *xf, int mk[], int *intx);
void ReadMap(FILE *inptr, double grid_points, double *data);
void WriteMap(FILE *outptr, double grid_points, double *data);

char headers[11][150], tails[5][150];

main (int argc, char *argv[])
{
	int i,i1,j,k,ndim, np1,np2, n1,n2, iflg, if1, **intx;
	int m1,m2,nderiv,ier,n, *nk,*mk;
	double xb1,h,xb2, *xl,*xu, aeps, fb, hb, **ax, ***cof, gb, mb, average_sigma, pos_exp, neg_exp, surf_z, sigma_ion;
	double eta_local, valency, area_per_lipid, j1, j2, j3, number_of_lipids, integrand, fourth_term, averaging_eta;  
	double **xf, dfxx, dfyy, dfxy, **x2, ***Ip, ***Is, ***Iv, ***Im, sigma, const_term4, average_eta, *dataPotential, *dataCharge, *dataKappa;
	double *origin, hx, hy, hz, n0, debye, bjerrum, T, jx, jy, jz, jp, const_term1, const_term2, const_term3, debye_inv_sq, lambda_D, *dataKappaIon;
	int grid_points, dim, counter, kk;
	FILE *paraptr, *paraintptr, *potptr, *chargeptr, *kappaptr, *kappaionptr;

	if ((paraintptr = fopen(argv[1], "r")) == NULL)
                printf ("\nERROR - Cannnot open the trajectory file\n");

        if ((chargeptr = fopen(argv[2], "r")) == NULL)
                printf ("\nERROR - Cannnot open charge map file\n");

	fscanf (paraintptr, "%d\n", &dim); /* get the dimensionality of the system*/

	xl = (double *) malloc(dim*(sizeof(double)));
	xu = (double *) malloc(dim*(sizeof(double)));
	nk = (int *) malloc(dim*(sizeof(int)));
	mk = (int *) malloc(dim*(sizeof(int)));
	origin = (double *) malloc(dim*(sizeof(double)));
	
/**************  Reading the rest of the parameter file  ***************/	

	fscanf (paraintptr, "%d\n", &k);   		/*  order of B-spline */
	fscanf (paraintptr, "%d\n", &nk[0]);  		/*  grid points in x direction */		
	fscanf (paraintptr, "%d\n", &nk[1]);		/*  grid points in y direction */
	fscanf (paraintptr, "%d\n", &nk[2]);  		/*  grid points in z direction */
	fscanf (paraintptr, "%lf\n", &hx); 		/*  mesh size in x direction, in A */		
	fscanf (paraintptr, "%lf\n", &hy);			/*  mesh size in y direction, in A */
	fscanf (paraintptr, "%lf\n", &hz); 		/*  mesh size in z direction, in A */
	fscanf (paraintptr, "%lf\n", &lambda_D); 		/*  The Debye length, in A */
	fscanf (paraintptr, "%lf\n", &T); 		        /*  Temperature */
	fscanf (paraintptr, "%lf\n", &xl[0]);  		/*  lower integration limit in x direction */
	fscanf (paraintptr, "%lf\n", &xu[0]);  		/*  upper integration limit in x direction */	
	fscanf (paraintptr, "%lf\n", &xl[1]);  		/*  lower integration limit in y direction */
	fscanf (paraintptr, "%lf\n", &xu[1]);		/*  upper integration limit in y direction */	
	fscanf (paraintptr, "%lf\n", &xl[2]);  		/*  lower integration limit in z direction */
	fscanf (paraintptr, "%lf\n", &xu[2]);  		/*  upper integration limit in z direction */	
        fscanf (paraintptr, "%lf\n", &surf_z);             /*  z coordinate of the surface */
        fscanf (paraintptr, "%lf\n", &average_sigma);      /*  average charge density on the surface */
        fscanf (paraintptr, "%lf\n", &sigma_ion);          /*  ion surface charge density */
        fscanf (paraintptr, "%lf\n", &origin[0]);           /*  origin x coordinate */
        fscanf (paraintptr, "%lf\n", &origin[1]);           /*  origin y coordinate */
        fscanf (paraintptr, "%lf\n", &origin[2]);           /*  origin z coordinate */
        fscanf (paraintptr, "%lf\n", &valency);           /*  valency */
        fscanf (paraintptr, "%lf\n", &area_per_lipid);           /*  area_per_lipid */
		
	number_of_lipids = nk[0]*nk[1]/area_per_lipid;
        bjerrum = 1.386*pow(10,6)/(eps_die*R*T);	 //this is the Bjerrum length in A, it is around 7.08 A at 300 K
	debye_inv_sq = 1/(lambda_D*lambda_D);
	const_term1 = 0.5;
        const_term2 = debye_inv_sq/(4*PI*bjerrum);	 //constant factor up front of mixing entropy term
	const_term3 = 1.0/area_per_lipid;
	const_term3 = 1.0/number_of_lipids;
	const_term4 = 0.5;
        average_eta = fabs(area_per_lipid*average_sigma/valency); /* this is average fraction of chargd lipids */
	grid_points = nk[0]*nk[1]*nk[2];

	dataCharge = (double *) malloc(grid_points*(sizeof(double)));

	//printf ("%lf  %lf  %lf %lf  %lf  %lf\n", sqrt(1/debye_inv_sq), bjerrum, const_term1, const_term2, average_eta, valency);
	
	ndim = nk[0];
	for (i=1; i<dim; ++i)
		if (nk[i]>ndim) ndim = nk[i];
		
	nderiv=0;
	iflg=0;

/*************** Allocating different multidimensional matrices ****************/

        Im = (double ***) malloc( nk[2] * sizeof(double **) );
        Im[0] = (double **) malloc( nk[2] * nk[1] * sizeof( double *) );
        Im[0][0] = (double *) malloc( nk[2] * nk[1] * nk[0] * sizeof(double) );
        for ( j = 1  ;  j < nk[1]  ;  ++j ) {
            Im[0][j] = Im[0][j-1] + nk[0];
        }
        for ( i = 1  ;  i < nk[2]  ;  ++i ) {
            Im[i] = Im[i-1] + nk[1];
            Im[i][0] = Im[i-1][nk[1]-1] + nk[0];
            for ( j = 1  ;  j < nk[1]  ;  ++j ) {
                Im[i][j] = Im[i][j-1] + nk[0];
            }
        }


        cof = (double ***) malloc( nk[2] * sizeof(double **) );
        cof[0] = (double **) malloc( nk[2] * nk[1] * sizeof( double *) );
        cof[0][0] = (double *) malloc( nk[2] * nk[1] * nk[0] * sizeof(double) );
        for ( j = 1  ;  j < nk[1]  ;  ++j ) {
            cof[0][j] = cof[0][j-1] + nk[0];
        }
        for ( i = 1  ;  i < nk[2]  ;  ++i ) {
            cof[i] = cof[i-1] + nk[1];
            cof[i][0] = cof[i-1][nk[1]-1] + nk[0];
            for ( j = 1  ;  j < nk[1]  ;  ++j ) {
                cof[i][j] = cof[i][j-1] + nk[0];
            }
        }

        intx = (int **) malloc(ndim * sizeof(int *));
	intx[0] = malloc(ndim * ndim * sizeof(int));
	for(i = 1; i < ndim; ++i)
		intx[i] = intx[0] + i * ndim;

	xf = (double **) malloc(dim * sizeof(double *));
	xf[0] = malloc(dim * ndim * sizeof(double));
	for(i = 1; i < dim; ++i)
		xf[i] = xf[0] + i * ndim;

        x2 = (double **) malloc(dim * sizeof(double *));
        x2[0] = malloc(dim * ndim * sizeof(double));
        for(i = 1; i < dim; ++i)
                x2[i] = x2[0] + i * ndim;

        ax = (double **) malloc(dim * sizeof(double *));
        ax[0] = malloc(dim * (ndim*k*dim) * sizeof(double));
        for(i = 1; i < dim; ++i)
                ax[i] = ax[0] + i * (ndim*k*dim);


/*** READ APBS MAPS ****/

        ReadMap(chargeptr, grid_points, dataCharge);


/*** ********* INPUT ENDS HERE ***********************************************/

	integrand = 0;

	for(i=0; i<nk[0]; ++i) x2[0][i]=origin[0]+i*hx;
	for(i=0; i<nk[1]; ++i) x2[1][i]=origin[1]+i*hy;
	for(i=0; i<nk[2]; ++i) x2[2][i]=origin[2]+i*hz;

	averaging_eta = 0;

        for (i=0; i<nk[0]; ++i) {
                for (j=0; j<nk[1]; ++j) {
                       for (kk=0; kk<nk[2]; ++kk) {
                                Im[kk][j][i] = 0;
                                if ((origin[2]+hz*kk) == surf_z) { /* membrane surface integrand */
                                      eta_local = fabs(area_per_lipid*dataCharge[i*nk[1]*nk[2]+j*nk[2]+kk]/valency);
                                      Im[kk][j][i] = eta_local*log(eta_local/average_eta) + (1-eta_local)*log((1-eta_local)/(1-average_eta));
				      averaging_eta = averaging_eta + eta_local;
					if (i==nk[0]/2 && (origin[1]+j*hy)>=0) {
	//					printf ("%lf %lf %lf %lf\n", origin[0]+i*hx, origin[1]+j*hy, origin[2]+kk*hz, Im[kk][j][i]);
						integrand = integrand + (origin[1]+j*hy)*Im[kk][j][i];
					}
                                }
                        }
                }
        }


	average_eta = averaging_eta/(nk[0]*nk[1]);
        fourth_term = number_of_lipids*(average_eta*log(average_eta)+(1-average_eta)*log(1-average_eta));


/********************* PROBLEM SET UP ENDS HERE *******************************/


/******** (Lipid Mixing term) ******/

        i=bspintn(dim,nk,&x2[0][0],ndim,&Im[0][0][0],k,&ax[0][0],&cof[0][0][0],
                   &xf[0][0],mk,&intx[0][0]);
//        printf(" ier = %d   No. of knots = %d %d %d\n", i, mk[0],mk[1],mk[2]);

        mb=bspqdn(dim,mk,&xf[0][0],ndim,k,&cof[0][0][0],xl,xu,&ier);
 //       for(i=0; i<dim; ++i) printf(" limits = %e,%e \n", xl[i],xu[i]);
        printf(" ier = %d  INTEGRAL = %lf  LIPID MIXING TERM = %lf  ", ier,mb,const_term3*mb);

	printf(" TOTAL INTEGRAL = %lf  ", const_term3*mb);
	printf(" FOURTH_TERM = %lf  ", fourth_term);
	printf(" INTEGRAND = %lf  %lf\n", integrand, integrand * 2 * PI * const_term3);

	free(Im);
	free(cof);
	free(xl);
	free(xu);
	free(nk);
	free(mk);
	free(intx);
	free(xf);
	free(x2);
	free(ax);
	free(dataCharge);
	fclose(paraintptr);
	fclose(chargeptr);

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
                   fprintf (outptr, "%lf  %lf\n", data[i], data[i+1]);
              else
                   fprintf (outptr, "%lf  %lf  %lf\n", data[i], data[i+1], data[i+2]);


                i = i+3;

        } while (i < grid_points);


        for (i=0; i<trailor; ++i) {
                fprintf (outptr, "%s\n", tails[i]);
        }

}

/*	To calculate coefficients for B-spline interpolation

	N : (input) Number of entries in the table
	X : (input) Array of length N containing the abscissas
	F : (input) Array of length N containing the function values
		F[I] is the tabulated function value at X[I].
	K : (input) Order of B-spline required. K=4 gives cubic B-splines
	A : (input/output) Array of length LA*3K containing the
		triangular decomposition of equation matrix in band form
		For IFLG=2, this array must be supplied, for other values
		of IFLG it is calculated by the function
	LA : (input) The second dimension of A as specified in calling function
		LA>=N
	C : (output) Coefficients of expansion, which will be calculated
		provided IFLG!=1
	XF : (input/output) Array of size NO, containing
		the knots used for B-spline calculations.
		The knots must be distinct and in ascending order.
		For IFLG=2, this array must be supplied, for other values
		of IFLG it is calculated by the function
	NO : (input/output) Number of knots for B-splines
		For IFLG=2, this number must be supplied, for other values
		of IFLG it is calculated by the function
	IFLG : (input/output) Integer specifying the type of calculation required
		IFLG=0 The matrix will be calculated and solved for coefficients
		IFLG=1 The matrix will be calculated and triangular decomposition
			is obtained, but coefficients are not calculated
		IFLG=2 The triangular decomposition of matrix is assumed
			to be available in A and coefficients C are calculated
		IFLG=-1 same as 0, except that no pivoting will be used
	INC : (input/output) Integer array containing information about
		pivoting during solution of system of linear equations
		For IFLG=2, this array must be supplied, for other values
		of IFLG it is calculated by the function
		
	Error status is returned by the value of the function BSPINT.
		0 value implies successful execution
		204 implies N<K or K<2
		other values may be set by BSPLIN or GAUBND

	Required functions : BSPLIN, GAUBND
*/

#include <math.h>

int gaubnd(int n, int kb, int num, double *a, double *x, double *det,
		int *idet, int inc[], int lj, int *iflg);
int bsplin(double *x, int nx, int k, double xb, int nderiv, double b[],
		double db[], double ddb[], int *left);

int bspint(int n, double x[], double f[], int k, double *a, int la, double c[],
	double *xf, int *no, int *iflg, int inc[])

{
	int i,j,i1,i2,kl,ku,ndb,nderiv,kb,left,num,idet;
	double xb,det;
	double *wk;

	if(n<k || k<2) return 204;

	if(*iflg < 2) {
/*	set up the knots for B-splines by dropping points near the ends */
		xf[0] = x[0];
		kl=(k-2)/2;
		ku=(k-1)/2;
		for(i=1+kl; i<=n-2-ku; ++i) xf[i-kl] = x[i];
		xf[n-kl-ku-1] = x[n-1];
		*no=n-kl-ku;
		ndb=n+1;
		nderiv=0;
		wk=(double *) calloc(2*n+5, sizeof(double));

/*	Set up the equation matrix for calculating coefficients of expansion
	The matrix is in band form A_{i,j} is stored in A[J-I+K-1][I] */
		for(i=0; i<n; ++i) {
			xb=x[i];
			j=bsplin(xf,*no,k,xb,nderiv,wk,(wk+ndb),(wk+ndb+2),&left);
			if(j>100) {free(wk); return j;}
			i1=i-k+1; if(i1<0) i1=0;
			i2=i+k-1; if(i2>n-1) i2=n-1;
			for(j=i1; j<=i2; ++j) a[la*(j-i+k-1)+i] = wk[j];
		}
		free(wk);
	}
	
/*	Solve the system of equations for a band matrix */
	num=1;
	kb=k-1;
	for(i=0; i<n; ++i) c[i]=f[i];
	i=gaubnd(n,kb,num,a,c,&det,&idet,inc,la,iflg);
	return i;
}



/*	To calculate coefficients for B-spline interpolation in n-dimensions
 
	N : (input) Number of dimensions
	NK : (input) Integer array of length N giving the number of
		tabular points in each dimension
	X : (input) Array of length NXD*N containing the abscissas
		X[j][i] is the ith abscissa in jth dimension
	NXD : (input) Second dimension of arrays X, XF, INTX as specified
		in the calling function. NXD must be greater than or equal
		to the maximum of NK[I]
	F : (input) Array of length NK[0]*NK[1]* ... *NK[N-1] containing
		the table of values. The dimensions of F in the calling
		function must exactly match the size of table so that
		there are no gaps in memory allocation.
	K : (input) Order of B-spline required. K=4 gives cubic B-splines
	AX : (output) Array of length NXD*3*K*N containing the
		triangular decomposition of equation matrix in band form
		for each dimension.
	C : (output) Array of length NK[0]*NK[1]*... *NK[N-1]
		containing the coefficients of expansion.
	XF : (output) Array of size NXD*N, containing
		the knots used for B-spline calculations.
		XF[j][i] is the ith knot along jth dimension.
	MK : (output) Integer array of length N containing number of
		knots for B-splines in each dimension.
	INTX : (output) Integer array of length NXD*N containing information
		about pivoting during solution of system of linear equations
		for each dimension.
		
	Error status is returned by the value of the function BSPINTN.
		0 value implies successful execution
		Nonzero values may be set by BSPINT, BSPLIN or GAUBND

	Required functions : BSPINT, BSPLIN, GAUBND
*/

#include <math.h>

int bspint(int n, double x[], double f[], int k, double *a, int la, double c[],
	double *xf, int *no, int *iflg, int inc[]);
int gaubnd(int n, int kb, int num, double *a, double *x, double *det,
		int *idet, int inc[], int lj, int *iflg);

int bspintn(int n, int nk[], double *x, int nxd, double *f, int k, double *ax,
	double *c, double *xf, int mk[], int *intx)

{
	int i,i1,i2,j,lj,iflg1,kb,nx,ny,no,ju,nx1,ny1,lj1,j1,idet;
	double det;
	double *wk;



/*	Calculate the triangular decomposition of matrices for interpolation
	along each dimension */
	for(i=0; i<n; ++i) {
		iflg1=1;
		lj=nk[i];
		j=bspint(lj,(x+i*nxd),f,k,(ax+i*nxd*3*k),lj,c,(xf+i*nxd),&mk[i],
			&iflg1,(intx+i*nxd));
		if(j>100) return j;
	}

	kb=k-1;
	nx=1;
	for(i=0; i<n-1; ++i) nx=nx*nk[i];
/*	no=nx*nk[n-1]+1 */
	ny=1;
	ju=n-1;
	wk=(double *) calloc(nx*nk[n-1], sizeof(double));


 
/*	If N is odd interpolate along last dimension outside the loop */
	if((n - 2*(n/2)) == 1) {
		lj=nk[n-1];
/*	Set up the RHS */
		for(i=0; i<nx; ++i) {
			for(j=0; j<nk[n-1]; ++j) wk[j+i*lj]=f[i+j*nx];
		}

		iflg1=2;
		i=gaubnd(nk[n-1],kb,nx,(ax+(n-1)*nxd*3*k),wk,&det,&idet,(intx+(n-1)*nxd),
			lj,&iflg1);
		if(i>100) {free(wk); return i;}
/*	Set up the RHS for next interpolation along N-2 th dimension */
		nx1=nx/nk[n-2];
		ny=nk[n-1];
		lj1=nk[n-2];
		for(i1=0; i1<nk[n-2]; ++i1) {
			for(i=0; i<nx1; ++i) {
				for(j=0; j<ny; ++j)
					c[i1+i*lj1+j*nx1*lj1] = wk[j+i*lj+i1*nx1*lj];
			}
		}
		nx=nx1;
		ju=n-2;
	}
	else {
/*	Set up the RHS for interpolation along N-1 th dimension */
		lj=nk[n-1];
		for(i=0; i<nx; ++i) {
			for(j=0; j<nk[n-1]; ++j) c[j+i*lj]=f[i+j*nx];
		}
	}
/*	Loop for interpolation in each dimension, each pass
	interpolates along 2 dimensions */
	for(j1=ju; j1>=0; j1=j1-2) {
		iflg1=2;
		lj=nk[j1];
		i=gaubnd(nk[j1],kb,nx*ny,(ax+j1*nxd*3*k),c,&det,&idet,(intx+j1*nxd),
			lj,&iflg1);
		if(i>100) {free(wk); return i;}

		
/*	Set up the RHS for interpolation along the next dimension */
		nx1=nx/nk[j1-1];
		ny1=ny*nk[j1];
		lj1=nk[j1-1];
		for(i1=0; i1<ny; ++i1) {
			for(i2=0; i2<nk[j1]; ++i2) {
				for(i=0; i<nk[j1-1]; ++i) {
					for(j=0; j<nx1; ++j)
						wk[i+j*lj1+i2*nx+i1*nx*nk[j1]]=c[i2+j*lj+i*nx1*lj+i1*nx*lj];
				}
			}
		}
		nx=nx1; ny=ny1;
		iflg1=2;
		lj=nk[j1-1];
		i=gaubnd(nk[j1-1],kb,nx*ny,(ax+(j1-1)*nxd*3*k),wk,&det,&idet,
			(intx+(j1-1)*nxd),lj,&iflg1);
		if(i>100) {free(wk); return i;}


		if(j1==1) {
/*	Store the coefficients in array */
			for(i=0; i<nk[0]; ++i) {
				for(j=0; j<ny; ++j) c[i+j*nk[0]]=wk[i+j*lj];
			}
		}
		else {

/*	Set up the RHS for interpolation along the next dimension */
			lj1=nk[j1-2];
			nx1=nx/nk[j1-2];
			ny1=ny*nk[j1-1];
			for(i1=0; i1<nk[j1-1]; ++i1) {
				for(i2=0; i2<nk[j1-2]; ++i2) {
					for(i=0; i<nx1; ++i) {
						for(j=0; j<ny; ++j)
							c[i2+i*lj1+i1*nx+j*nx*nk[j1-1]]=wk[i1+i*lj+i2*lj*nx1+j*lj*nx];
					}
				}
			}
			nx=nx1; ny=ny1;
		}
	}
	free(wk);
	return 0;
}



/*	To calculate the B-spline basis functions at a specified point

	X : (input) Array of length NX containing the knots.
		The knots must be distinct and in ascending order.
	NX : (input) Number of knots
	K : (input) Order of B-spline, 0< K, K=4 gives cubic B-splines
	XB : (input) The point at which B-spline basis functions are to be evaluated
	NDERIV : (input) Number of derivatives required
		NDERIV<=0 only B-splines are calculated
		NDERIV=1 first derivative is also calculated
		NDERIV>1 first and second derivatives are also calculated
	B : (output) Array of length NX+K-2 containing the value of
		B-spline basis functions
	DB : (output) Array of length NX+K-2 containing the value of
		the first derivative of B-spline basis functions (if NDERIV>0)
	DDB : (output) Array of length NX+K-2 containing the value of
		the second derivative of B-spline basis functions (if NDERIV>1)
	LEFT : (output) XB is located between X[LEFT] and X[LEFT+1]
		
	Error status is returned by the value of the function BSPLIN.
		0 value implies successful execution
		26 implies XB > X[NX-1]
		27 implies XB < X[0]
		203 implies NX<2, K<1

	Required functions : None
*/

#include <math.h>

int bsplin(double *x, int nx, int k, double xb, int nderiv, double b[],
		double db[], double ddb[], int *left)

{
	int i,j,igh, mid,nigh,lx,ier;
	double t1,t2,t3,p1,p2;
	double *wk, *dr, *dl;
	static int low = -1;

	if(nx <= 1 || k<1 ) return 203;
	ier=0;

/*	If the previous value of LOW is inadmissible, set the range to (0,N-1) */
	if(low<0 || low>=nx-1) {low=0; igh=nx-1;}
	else igh=low+1;

	while((xb<x[low] && xb<x[igh]) || (xb>x[low] && xb>x[igh])) {
/*	Extend the range */
		if( xb>x[low] ) {
/*	Extend the range on higher side */
			if(igh >= nx-1) {ier=26; low=nx-2; break;}
			else {
				nigh=igh+2*(igh-low); if(nx-1 < nigh) nigh=nx-1;
				low=igh; igh=nigh;
			}
		}

		else {
/*	Extend the range on lower side */
			if(low <= 0) {ier=27; igh=low+1; break;}
			else {
				nigh=low;
				low=low-2*(igh-low); if(low<0) low=0;
				igh=nigh;
			}
		}
	}


/*	Once the point is bracketed between two tabular points locate it by bisection */
	while((igh-low > 1) && (xb != x[low])) {
		mid=(low+igh)/2;
		if((xb<= x[mid]) == (xb<= x[low])) low=mid;
		else igh=mid;
	}

 
/*	Evaluate the B-spline basis functions

	Define the extra knots on either side of table
	Note that the function assumes knots from -K+2 to NX+K-1
	and the B-splines B_{i,k}, i ranges from 0 to NX+K-3 
	The knots are stored in scratch array wk. */

	wk=(double *) calloc(nx+2*k+2, sizeof(double));
	dr=(double *) calloc(nx+2*k+2, sizeof(double));
	dl=(double *) calloc(nx+2*k+2, sizeof(double));

	for(i=0; i<nx; ++i) wk[i+k]=x[i];
	for(i=1; i<=k; ++i) {
		wk[k-i]=x[0];
		wk[nx+i-1+k]=x[nx-1];
	}

	for(i=0; i<nx+k-2; ++i) {b[i]=0.0; db[i]=0.0; ddb[i]=0.0;}
	*left=low;
	lx=low-1;
	b[lx+1]=1;

/*	The recurrence relation for B-splines */
	for(j=1; j<=k-1; ++j) {
		dr[j] = wk[low+j+k] - xb;
		dl[j] = xb - wk[low+1-j+k];
		t1=0.0;
		for(i=1; i<=j; ++i) {
			t2=b[lx+i]/(dr[i]+dl[j+1-i]);
			b[lx+i]=t1+t2*dr[i];
			t1=t2*dl[j+1-i];
		}
		b[lx+j+1]=t1;
			
/*	Calculate the first derivative using recurrence relations */
		if(j == k-2 && nderiv > 0) {
			t1=0.0;
			for(i=1; i<=j+1; ++i) {
				t2=b[lx+i]/(wk[low+i+k]-wk[low+i+1]);
				db[lx+i]=(k-1)*(t1-t2);
				t1=t2;
			}
			db[lx+j+2]=(k-1)*t1;
		}
 
/*	Calculate the second derivative using recurrence relations */
		if(j == k-3 && nderiv>1) {
			t2=0.0; p1=0.0;
			for(i=1; i<=j+1; ++i) {
				t3=b[lx+i]/(wk[low+i+k]-wk[low+i+2]);
				p2=(t2-t3)/(wk[low+i+k]-wk[low+i+1]);
				ddb[lx+i]=(k-2)*(k-1)*(p1-p2);
				t2=t3; p1=p2;
			}
			p2=t2/(wk[low+j+2+k]-wk[low+j+3]);
			ddb[lx+j+2]=(k-2)*(k-1)*(p1-p2);
			ddb[lx+j+3]=(k-2)*(k-1)*p2;
		}
	}

/*	For K=2 the first derivative has to be calculated outside the loop */
	if(k == 2 && nderiv > 0) {
		t2=1./(wk[low+1+k]-wk[low+2]);
		db[lx+1]=-t2;
		db[lx+2]=t2;
	}

/*	For K=3 the second derivative has to be calculated outside the loop */
	if(k == 3 && nderiv > 1) {
		t3=1./(wk[low+1+k]-wk[low+3]);
		p2=-t3/(wk[low+1+k]-wk[low+2]);
		ddb[lx+1]=-2.*p2;
		p1=p2;
		p2=t3/(wk[low+2+k]-wk[low+3]);
		ddb[lx+2]=2.*(p1-p2);
		ddb[lx+3]=2.*p2;
	}
	free(dl); free(dr); free(wk);
	return ier;
}



/*	To calculate the integral of a function defined by B-spline expansion

	N : (input) Number of knots for B-spline expansion
	X : (input) Array of length N containing the knots.
		The knots must be distinct and in ascending order.
	K : (input) Order of B-splines, K=4 for cubic B-splines
	WT : (input) Array of length N+K-2 containing the coefficients
		of B-spline expansion
	XL : (input) Lower limit of integration
	XU : (input) Upper limit of integration
	IER : (output) Error parameter, IER=0 implies successful execution
		IER=31 implies that XL is outside the range of table
		IER=32 implies that XU is outside the range of table

	BSPQD = Integral of \sum WT[I]\phi_I(x) over [XL,XU]

	Required functions : BSPLIN
*/

#include <math.h>

int bsplin(double *x, int nx, int k, double xb, int nderiv, double b[],
                double db[], double ddb[], int *left);

double bspqd(int n, double *x, int k, double wt[], double xl, double xu, int *ier)

{
	int i,j,k1,jk,lf1,lf2,nderiv;
	double f,f1,f2;
	double *wk1, *wk2, *wk;

	*ier=0;
	if(xl == xu) return 0.0;
 
/*	Calculate the B-spline basis of order K+1 at XL and XU */
	k1=k+1;
	nderiv=0;
	wk=(double *)calloc(n+3*k, sizeof(double));
	wk1=(double *)calloc(n+k, sizeof(double));
	wk2=(double *)calloc(n+k, sizeof(double));

	*ier=bsplin(x,n,k1,xl,nderiv,wk1,wk,wk,&lf1);
	if(*ier>100) {free(wk2); free(wk1); free(wk); return 0.0;}

	*ier=bsplin(x,n,k1,xu,nderiv,wk2,wk,wk,&lf2);
	if(*ier>100) {free(wk2); free(wk1); free(wk); return 0.0;}

	if(xl<x[0] || xl>x[n-1]) *ier=31;
	if(xu<x[0] || xu>x[n-1]) *ier=32;

	for(i=0; i<n; ++i) wk[k+i]=x[i];
	for(i=1; i<=k; ++i) {
		wk[k-i]=x[0];
		wk[n+i-1+k]=x[n-1];
	}
 
/*	The sum for x=XL */
	f1=0.0;
	for(i=lf1; i<=lf1+k; ++i) {
		f=0.0;
		for(j=1; j<=i; ++j) f=f+wt[j-1]*(wk[k+j]-wk[j]);
		f1=f1+f*wk1[i]/k;
	}

/*	The sum for x=XU */
	f2=0.0;
	for(i=lf2; i<=lf2+k; ++i) {
		f=0.0;
		for(j=1; j<=i; ++j) f=f+wt[j-1]*(wk[k+j]-wk[j]);
		f2=f2+f*wk2[i]/k;
	}

	free(wk2); free(wk1); free(wk);
	return f2-f1;
}



/*	To calculate the integral of a function defined by B-spline expansion
	over hyper-rectangular region in N dimensions

	N : (input) Number of dimensions
	NK : (input) Integer array of length N containing the number of knots
		for B-spline representation along each direction
	X : (input) Array of length NXD*N containing the knots along each direction
		X[j][i] is the ith knot along jth dimension
	NXD : (input) Second dimension of array X as declared in the calling function
		NXD must be greater than or equal to maximum of NK[I]
	K : (input) Order of B-splines, K=4 for cubic B-splines
	WT : (input) Array of length
		(NK[0]+K-2)(NK[1]+K-2) ... (NK[N-1]+K-2)
		containing coefficients of B-spline expansion. 
		The coefficients are assumed to be stored in natural
		Fortran order with no gaps in data. In the calling
		function the array should have dimension
		WT[NK[N-1]+K-2]...[NK[1]+K-2][NK(0)+K-2]
	XL : (input) Array of length N containing the lower limits of
		integration along each direction
	XU : (input) Array of length N containing the upper limits of
		integration along each direction
	IER : (output) Error parameter, IER=0 implies successful execution
		Nonzero values may be set by BSPQD which is called

	BSPQDN = Integral over[XL[0],XU[0]] x ... x [XL[N-1],XU[N-1]] of
	\sum WT[i_N-1]...[i_1][i_0]\phi_i0(x_0)\phi_i1(x_1) ... \phi_in-1(x_n-1) 

	Required functions : BSPLIN, BSPQD
*/

#include <math.h>

double bspqd(int n, double *x, int k, double wt[], double xl, double xu, int *ier);

double bspqdn(int n, int nk[], double *x, int nxd, int k, double *wt,
	double xl[], double xu[], int *ier)
 
{
	int i,j,n1,nk1,n0,na,nb,nt;
	double ri;
	double *wk;

	n1=1;
	for(i=1; i<n; ++i) n1=n1*(nk[i]+k-2);
	wk=(double *) calloc(2*n1+6, sizeof(double));
	n0=n1+3;
	
/*	Integration along the first dimension */
	nk1=nk[0]+k-2;
	for(j=0; j<n1; ++j) {
		wk[j]=bspqd(nk[0],&x[0],k,&wt[j*nk1],xl[0],xu[0],ier);
		if(*ier >100) {free(wk); return 0.0;}
	}
 
	na=0; nb=n0;
/*	Integrate along the remaining dimensions */
	for(i=1; i<n; ++i) {
		nk1=nk[i]+k-2;
		n1=n1/nk1;
		for(j=0; j<n1; ++j) {
			wk[nb+j]=bspqd(nk[i],&x[i*nxd],k,&wk[na+j*nk1],xl[i],xu[i],ier);
			if(*ier >100) {free(wk); return 0.0;}
		}
		nt=na; na=nb; nb=nt;
	}
	ri=wk[na];
	free(wk);
	return ri;
}




/*     Solution of a system of linear equations using Gaussian elimination
     	for a band matrix

	N : (input) Number of equations to be solved
	KB : (input) Bandwidth of matrix A[i][j]=0 if abs(I-J)>KB
	NUM : (input) Number of different sets (each with N equations) of
		equations to be solved
	A : (input/output) The matrix of coefficients of size LJ*(3*KB+1)
		A[J-I+KB][I] is the coefficient of x_J in Ith equation
		at output it will contain the triangular decomposition
	X : (input/output) The matrix containing right hand sides (size LJ*NUM)
		X[j][i] is the ith element of jth right hand side
		at output it will contain the solutions
	DET, IDET : (output) The determinant of the matrix = DET*2**IDET
	INC : (output) Integer array of length N containing information about
		interchanges performed during elimination
	LJ : (input) Second dimension of arrays A and X in calling function
	IFLG : (input) Integer variable used as a flag to specify the type
		of computation required
     	If IFLG=-1, both elimination and solution are calculated
     		without pivoting and IFLG is set to 2
	If IFLG=0, both elimination and solution are computed
     		with partial pivoting and IFLG is set to 2
     	If IFLG=1, only elimination is done with pivoting and IFLG is set to 2
     	If IFLG>=2 only solution is calculated, the triangular
     		decomposition should have been calculated earlier
		
	Error status is returned by the value of the function GAUBND.
		0 value implies successful execution
		104 implies N<=0 or N>LJ or KB>N
		124 implies some pivot turned out to be zero and hence
	     		matrix must be nearly singular

	Required functions : None
*/

#include <math.h>

int gaubnd(int n, int kb, int num, double *a, double *x, double *det,
		int *idet, int inc[], int lj, int *iflg)

{
	int i,j,k,km,l,kb1,m1,m2;
	double r1, t1;
	double *wk;

	if(n<=0 || n>lj || kb>n) return 104;

	kb1=kb+1;
	if(*iflg < 2) {
/*     Perform elimination */
		for(i=0; i<n ; ++i) {
			for(j=2*kb+1; j<=3*kb; ++j) a[j*lj+i] = 0.0;
		}

		*det=1.0; *idet=0;
		wk=(double *) calloc(3*kb+1, sizeof(double));
		for(k=0; k<n-1; ++k) {
/*     Find the maximum element in the Kth column  */
			r1=0.0; km=k;
			if(*iflg >= 0) {
				m1=k+kb; if(n-1 < m1) m1=n-1;
				for(l=k; l<=m1; ++l) {
					if(fabs(a[l+(k-l+kb)*lj]) > r1) {
						r1=fabs(a[l+(k-l+kb)*lj]);
						km=l;
					}
				}
			}

			inc[k]=km;
			if(km != k) {
/*     Interchange the rows if needed  */
				m1=2*kb+k; if(n-1 < m1) m1=n-1;
				for(l=k; l<=m1; ++l) wk[l-k]=a[k+(l-k+kb)*lj];
				for(l=k; l<=m1; ++l) {
					a[k+(l-k+kb)*lj]=a[km+(l-km+kb)*lj];
					a[km+(l-km+kb)*lj]=wk[l-k];
				}
				*det= -(*det);
			}

			*det = (*det)*a[kb*lj+k];
			if( a[kb*lj+k] == 0.0) {free(wk); return 124;}
			if(*det != 0.0) {
/*     Scale the value of the determinant   */
				while(fabs(*det) > 32.0) {
					*det = (*det)*0.03125e0; *idet = *idet + 5;
				}

				while(fabs(*det) < 0.03125e0) {
					*det = (*det)*32.0; *idet = *idet - 5;
				}
			}

			m1=k+kb; if(n-1 < m1) m1=n-1;
			for(l=k+1; l<= m1; ++l) {
				a[l+(k-l+kb)*lj] = a[l+(k-l+kb)*lj]/a[kb*lj+k];
				m2=k+2*kb; if(n-1 < m2) m2=n-1;
				for(i=k+1; i<=m2; ++i)
					a[l+(i-l+kb)*lj]=a[l+(i-l+kb)*lj]-a[l+(k-l+kb)*lj]*a[k+(i-k+kb)*lj];
			}
		}

		free(wk);
		*det = (*det)*a[(n-1)+kb*lj];
		inc[n-1]=n-1;
		if(a[(n-1)+kb*lj] == 0.0) return 124;

		if(*iflg==1) {*iflg=2; return 0;}
		*iflg=2;
	}
		
/*     Solution for the NUM different right-hand sides */
	for(j=0; j<num; ++j) {
/*     Forward substitution  */
		for(k=0; k<n-1; ++k) {
			if(k != inc[k]) {
				t1= x[j*lj+k];
				x[j*lj+k] = x[j*lj+inc[k]];
				x[j*lj+inc[k]] = t1;
			}
			m1=k+kb; if(n-1 < m1) m1=n-1;
			for(l=k+1; l<=m1; ++l)
				x[j*lj+l]=x[j*lj+l]-a[l+(k-l+kb)*lj]*x[j*lj+k];
		}

/*     back-substitution  */
		x[j*lj+n-1] = x[j*lj+n-1]/a[(n-1)+kb*lj];
		for(k=n-2; k>=0; --k) {
			m1=k+2*kb; if(n-1 < m1) m1=n-1;
			for(l=m1; l>=k+1; --l)
				x[j*lj+k]=x[j*lj+k]-x[j*lj+l]*a[k+(l-k+kb)*lj];
			x[j*lj+k] = x[j*lj+k]/a[kb*lj+k];
		}
	}
	return 0;
}

