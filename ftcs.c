#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>
#include<assert.h>
#define SWAP(a,b) {double temp=(a);(a)=(b);(b)=temp;} 

#define DEBUG 0 

void forwardTime(double ***matrix, double ***matrix_np1, int nx, int ny, int nz, double dx, double dy, double dz, double C, double q, int border, int NT);
void crank(double ***A, double ***B, int nx, int ny, int nz, double alpha, double dx, double C, int border, int NT);
void adi(double ***matrix, double ***matrix_np1, int nx, int ny, int nz, double alpha, double dx, double C, int border, int NT);
double ***initGaussians(int nx, int ny, int nz, double dx, double dy, double dz, double noise);
void borders(double ***matrix, double ***matrix_np1, double *zcoords, double *zcoords_np1, int nx, int ny, int nz, double dx, double dy, double dz, double noise, char *border);
void freeMats( double ***matrix, double ***matrix_np1, int nx, int ny);
void print3DMatrix(double ***matrix, int nx, int ny, int nz, double dx, double dy, double dz, int when);
void rref(double **matrix, int n, int m);
void rowSwap(double **A, int swap, int with, int n, int m);
void setMatrix(double **A, int n,  int m);
void onesPivot(double **A, int row, int col, int m);
void multAdd(double **A, int row, int row2, int col, int m);
double distanceFormula(double *mat_old, double *mat_new, int length);
void welcome(int *nx, int *ny, int *nz, double *alpha, double *q, double *dt, double *dx, int *algo, int *border );
void rref2(double **A, int m, int n);
double *dvector(long nl, long nh);


int main(int argc, char *argv[]){
    int i, j, l, m, n;
    /* Matrix Dimensions*/
    int nx=4, ny=4, nz=4;

    /*parameters to be changed by user*/
    int NT; /* Max Number of timesteps*/
    double alpha = .01;
    double dt = .001;
    double q = .005;
    double dx = .01; 
    double dz = .01;
    double dy = .01;
    double C = alpha*dt/pow(dx,2);
    int algo = 0;
    int border = 0;
    srand(time(NULL));/*make some noise*/
    double noise =  .09 * (double)(rand()/(double)RAND_MAX);

    welcome(&nx, &ny, &nz, &alpha, &q, &dt, &dx, &algo, &border );

    /*    if(argc >1){
	  FILE *fp;
	  fp  = fopen(argv[1], "r");
	  fscanf(fp, "%d %d %d %d %lf %lf %lf %d %d", &nx, &ny, &nz, &NT, &dt, &alpha, &q, &border, &algo);
	  fclose(fp);
	  }
	  */
    assert(alpha * dt <= (0.5 * dx * dx));

    /*2 3D matrices init under gaussian conditions + boundaries*/
    double ***matrix =  initGaussians(nx, ny, nz, dx, dy, dz, noise);
    double ***matrix_np1 =  initGaussians(nx, ny, nz, dx, dy, dz, noise);
    int when = 0;
    print3DMatrix(matrix, nx, ny, nz, dx, dx, dx, when);

    double t_start, t_end;
    if( algo ==0 ){
	NT = 10000;
	t_start = clock();
	forwardTime(matrix, matrix_np1, nx, ny, nz, dx, dy, dz, C, q,border, NT);
	t_end = clock();
    }	
    else if( algo == 1 ){
	NT = 5;
	t_start = clock();
	crank(matrix, matrix_np1, nx, ny, nz, alpha, dx, C, border, NT);
	t_end = clock();
    }
    else if( algo == 2){
	NT = 5;
	t_start = clock();
	adi(matrix, matrix_np1, nx, ny, nz, alpha, dx, C, border, NT);
	t_end = clock();

    }
    else{
	printf("error in parameters");
	return 0;
    }
    freeMats(matrix, matrix_np1, nx, ny);
    printf("Total Time Elapsed.....\n");
    printf("%f\n",(t_end - t_start)/CLOCKS_PER_SEC);
    return 0;
}

void crank(double ***matrix, double ***matrix_np1, int nx, int ny, int nz, double alpha, double dx, double C, int border, int NT){
    int i, j, k, l, f, g, h;
    double d;
    int b0;
    int bn;
    int when = 1;
    double ***temp;
    double *tempb;
    double *mat_b = (double*) malloc(sizeof(double)* nx*ny*nz);
    double *mat_newb = (double*) malloc(sizeof(double)* nx*ny*nz);
    double **mat_A = (double**) malloc(sizeof(double*) *nx*ny*nz);
    double **augmented = (double**) malloc(sizeof(double*) *nx*ny*nz);
    double *data = (double*) malloc(sizeof(double) * nx*nx*ny*ny*nz*nz);
    double *augdata = (double*) malloc(sizeof(double) * nx*nx*ny*ny*nz*nz+nx*ny*nz);/*care for last column - b*/
	    /*now need to make augmented matrix to pass into rref*/
	    for(h=0;h<nx*ny*nz;h++){
		*(augmented+h) = augdata + h*nx*ny*nz+1;/*care for augmented with b*/
	    }
    switch(border){
	case 0:
	    b0 = 0;
	    bn = 0;
	    break;
	case 1:
	    b0 = 3;
	    bn = 3;
	    break;
	case 2:
	    b0 = nx;
	    bn = 0;
	    break;
	default :
	    b0 = 0;
	    bn = 0;
    }
    print3DMatrix(matrix, nx, ny, nz, dx, dx, dx, when);
    when++;
    /*malloc space for A*/
    for(i=0;i<nx*ny*nz;i++){
	*(mat_A+i) = data + i*nx*ny*nz;
    }
    /*make matrix b which is 1-D from init 3-D matrix --flattened*/
    for(i = 0 ; i< nx ; i++){
	for(j = 0 ; j < ny ; j++){
	    for(k = 0 ; k < nz ; k++){
		*(mat_b + i*ny*nz + j*nz + k) = *(*(*(matrix + i) + j) +k);
	    }
	}
    }

    for(k=0;k<NT;k++){


	for(i = 0; i< nx*ny*nz ; i++){
	    for(j = 0;j< nx*ny*nz ; j++){
		if( i == j){
		    *(*(mat_A+i)+j) = /*mat_b[i]**/1+6.0*C;
		}	
		else if (i == j-1 ){
		    if( i != 0%nz){
			*(*(mat_A+i)+j) = /**(mat_b + i-1)**/(-1.0*C);
		    }
		    else if(i ==0%nz && border== 2){
			*(*(mat_A+i)+j) = /**(mat_b + (i-1)%nz)**/(-1.0*C);

		    }

		}
		else if (i == j+1){
		    if( i != (nz-1)%(nz-1)){
			*(*(mat_A+i)+j) =/* *(mat_b + i+1)**/(-1.0*C);
		    }
		    else if(i == (nz-1)%(nz-1) && border == 2){
			*(*(mat_A+i)+j) =/* *(mat_b + (i+1)%nz)**/(-1.0*C);

		    }

		}
		else if(i == j-ny*nz){
		    if( i > ny*nz){
			*(*(mat_A+i)+j) =/* *(mat_b + i-ny*nz)**/(-1.0*C);
		    }
		    else if(i <= ny*nz && border == 2){
			*(*(mat_A+i)+j) =/* *(mat_b + (i-ny*nz)%(nx*ny*nz))**/(-1.0*C);

		    }
		}
		else if(i == j+ny*nz){
		    if( i < ny*nz){
			*(*(mat_A+i)+j) = /**(mat_b + i+ny*nz)**/(-1.0*C);
		    }
		    else if(i >= ny*nz && border == 2){
			*(*(mat_A+i)+j) = /**(mat_b + (i+ny*nz)%(nx*ny*nz))**/(-1.0*C);

		    }

		}
		else if(i == j+ny){
		    if( i != (ny-1)%(ny-1)){
			*(*(mat_A+i)+j) =/* *(mat_b + i+ny)**/(-1.0*C);
		    }
		    else if(i ==(ny-1)%(ny-1) && border == 2){
			*(*(mat_A+i)+j) =/* *(mat_b + (i+ny)%ny)**/(-1.0*C);

		    }

		}
		else if(i == j-ny){
		    if( i != 0%ny){
			*(*(mat_A+i)+j) =/* *(mat_b + i-ny)**/(-1.0*C);
		    }
		    else if(i == 0%ny && border == 2){
			*(*(mat_A+i)+j) =/* *(mat_b + (i-ny)%ny)**/(-1.0*C);

		    }

		}
		else{
		    *(*(mat_A+i)+j) = 0.0;

		}
	    }

	    /*at this point A is created*/
	    for(g = 0; g< nx*ny*nz ; g++){
		for(h = 0;h< nx*ny*nz ; h++){
		    if(h == nx*ny*nz-1){
			*(*(augmented+g)+h) = *(mat_b+g);
		    }
		    else{
			*(*(augmented+g)+h) = *(*(mat_A+g)+h);
		    }
		}
	    }
	    /*at this point i have Ax=b in augmented matrix form to use rref on {augmented | b}  */
	    rref(augmented, nx*ny*nz,  nx*ny*nz+1);
	    /*now i need to extract the nx*ny*nz+1 column to cpy back to my old mat_b, mat_A*/
	    for(g = 0; g< nx*ny*nz ; g++){
		for(h = 0;h< nx*ny*nz ; h++){
		    if(h == nx*ny*nz-1){
			*(mat_newb+g)  = *(*(augmented+g)+h);
		    }
	/*	    else{
			*(*(mat_A+g)+h) = *(*(augmented+g)+h);
		    }*/
		}
	    }
	    d = distanceFormula(mat_newb, mat_b, nx*ny*nz);
	    //printf("%f \n",d);
	    tempb = mat_b;
	    mat_b = mat_newb;
	    mat_newb = tempb;

	    for(f = 0 ; f< nx ; f++){
		for(g = 0 ; g < ny ; g++){
		    for( h = 0 ; h < nz ; h++){
			*(*(*(matrix + f) + g) +h) =  *(mat_b + f*ny*nz + g*nz + h);
		    }
		}
	    }
	}
	if( k==0){
	    print3DMatrix(matrix, nx, ny, nz, dx, dx, dx, when);
	    when++;
	}
    }
}


void adi(double ***matrix, double ***matrix_np1, int nx, int ny, int nz, double alpha, double dx, double C, int border, int NT){
    int i, j, k, l, f, g, h;
    int dir;
    int when =1;
    double d;
    double ***temp;
    double *tempb;
    double *mat_b = (double*) malloc(sizeof(double)* nx*ny*nz);
    double *mat_newb = (double*) malloc(sizeof(double)* nx*ny*nz);
    double **mat_A = (double**) malloc(sizeof(double*) *nx*ny*nz);
    double **augmented = (double**) malloc(sizeof(double*) *nx*ny*nz);
    double *data = (double*) malloc(sizeof(double) * nx*nx*ny*ny*nz*nz);
    double *augdata = (double*) malloc(sizeof(double) * nx*nx*ny*ny*nz*nz+nx*ny*nz);/*care for last column - b*/
    print3DMatrix(matrix, nx, ny, nz, dx, dx, dx, when);
    when++;
	    /*now need to make augmented matrix to pass into rref*/
	    for(h=0;h<nx*ny*nz;h++){
		*(augmented+h) = augdata + h*nx*ny*nz+1;/*care for augmented with b*/
	    }

    /*malloc space for A*/
    for(i=0;i<nx*ny*nz;i++){
	*(mat_A+i) = data + i*nx*ny*nz;
    }
    /*make matrix b which is 1-D from init 3-D matrix --flattened*/
    for(i = 0 ; i< nx ; i++){
	for(j = 0 ; j < ny ; j++){
	    for(k = 0 ; k < nz ; k++){
		*(mat_b + i*ny*nz + j*nz + k) = *(*(*(matrix + i) + j) +k);
	    }
	}
    }
    /*run in each direction, without X, then without Y, then without Z, so A is made up of 3 different 2 D variations*/
    for(k=0;k<NT;k++){
	for(dir=0;dir<3;dir++){
	    switch(dir){
		case 0:
		    for(i = 0; i< nx*ny*nz ; i++){
			for(j = 0;j< nx*ny*nz ; j++){
			    if( i == j){
				*(*(mat_A+i)+j) = /*mat_b[i]**/1+4.0*C;
			    }	
			    else if (i == j-1){
				if( i != 0%nz){
				    *(*(mat_A+i)+j) = /**(mat_b + i-1)**/(-1.0*C);
				}
				else if(i ==0%nz && border== 2){
				    *(*(mat_A+i)+j) = /**(mat_b + (i-1)%nz)**/(-1.0*C);

				}

			    }
			    else if (i == j+1){
				if( i != (nz-1)%(nz-1)){
				    *(*(mat_A+i)+j) = /**(mat_b + i+1)**/(-1.0*C);
				}
				else if(i == (nz-1)%(nz-1) && border == 2){
				    *(*(mat_A+i)+j) =/* *(mat_b + (i+1)%nz)**/(-1.0*C);

				}

			    }
			    else if(i == j+ny){
				if( i != (ny-1)%(ny-1)){
				    *(*(mat_A+i)+j) = /**(mat_b + i+ny)**/(-1.0*C);
				}
				else if(i ==(ny-1)%(ny-1) && border == 2){
				    *(*(mat_A+i)+j) =/* *(mat_b + (i+ny)%ny)**/(-1.0*C);

				}

			    }
			    else if(i == j-ny){
				if( i != 0%ny){
				    *(*(mat_A+i)+j) = /**(mat_b + i-ny)**/(-1.0*C);
				}
				else if(i == 0%ny && border == 2){
				    *(*(mat_A+i)+j) = /**(mat_b + (i-ny)%ny)**/(-1.0*C);

				}

			    }
			    else{
				*(*(mat_A+i)+j) = 0.0;

			    }
			}
		    }
		    break;
		case 1:
		    for(i = 0; i< nx*ny*nz ; i++){
			for(j = 0;j< nx*ny*nz ; j++){
			    if( i == j){
				*(*(mat_A+i)+j) =/* mat_b[i]**/1+4.0*C;
			    }	
			    else if (i == j-1){
				if( i != 0%nz){
				    *(*(mat_A+i)+j) = /**(mat_b + i-1)**/(-1.0*C);
				}
				else if(i ==0%nz && border== 2){
				    *(*(mat_A+i)+j) = /**(mat_b + (i-1)%nz)**/(-1.0*C);

				}

			    }
			    else if (i == j+1){
				if( i != (nz-1)%(nz-1)){
				    *(*(mat_A+i)+j) = /**(mat_b + i+1)**/(-1.0*C);
				}
				else if(i == (nz-1)%(nz-1) && border == 2){
				    *(*(mat_A+i)+j) = /**(mat_b + (i+1)%nz)**/(-1.0*C);

				}

			    }
			    else if(i == j-ny*nz){
				if( i > ny*nz){
				    *(*(mat_A+i)+j) = /**(mat_b + i-ny*nz)**/(-1.0*C);
				}
				else if(i <= ny*nz && border == 2){
				    *(*(mat_A+i)+j) =/* *(mat_b + (i-ny*nz)%(nx*ny*nz))**/(-1.0*C);

				}

			    }
			    else if(i == j+ny*nz){
				if( i < ny*nz){
				    *(*(mat_A+i)+j) = /**(mat_b + i+ny*nz)**/(-1.0*C);
				}
				else if(i >= ny*nz && border == 2){
				    *(*(mat_A+i)+j) = /**(mat_b + (i+ny*nz)%(nx*ny*nz))**/(-1.0*C);

				}

			    }
			    else{
				*(*(mat_A+i)+j) = 0.0;

			    }
			}
		    }
		    break;
		case 2:
		    for(i = 0; i< nx*ny*nz ; i++){
			for(j = 0;j< nx*ny*nz ; j++){
			    if( i == j){
				*(*(mat_A+i)+j) = /*mat_b[i]**/1+4.0*C;
			    }	
			    else if(i == j-ny*nz){
				if( i > ny*nz){
				    *(*(mat_A+i)+j) = /**(mat_b + i-ny*nz)**/(-1.0*C);
				}
				else if(i <= ny*nz && border == 2){
				    *(*(mat_A+i)+j) = /**(mat_b + (i-ny*nz)%(nx*ny*nz))**/(-1.0*C);

				}

			    }
			    else if(i == j+ny*nz){
				if( i < ny*nz){
				    *(*(mat_A+i)+j) = /**(mat_b + i+ny*nz)**/(-1.0*C);
				}
				else if(i >= ny*nz && border == 2){
				    *(*(mat_A+i)+j) = /**(mat_b + (i+ny*nz)%(nx*ny*nz))**/(-1.0*C);

				}

			    }
			    else if(i == j+ny){
				if( i != (ny-1)%(ny-1)){
				    *(*(mat_A+i)+j) = /**(mat_b + i+ny)**/(-1.0*C);
				}
				else if(i ==(ny-1)%(ny-1) && border == 2){
				    *(*(mat_A+i)+j) = /**(mat_b + (i+ny)%ny)**/(-1.0*C);

				}

			    }
			    else if(i == j-ny){
				if( i != 0%ny){
				    *(*(mat_A+i)+j) = /**(mat_b + i-ny)**/(-1.0*C);
				}
				else if(i == 0%ny && border == 2){
				    *(*(mat_A+i)+j) = /**(mat_b + (i-ny)%ny)**/(-1.0*C);

				}

			    }
			    else{
				*(*(mat_A+i)+j) = 0.0;

			    }
			}
		    }
		    break;

	    }



	    /*at this point A is created*/
	    for(g = 0; g< nx*ny*nz ; g++){
		for(h = 0;h< nx*ny*nz ; h++){
		    if(h == nx*ny*nz-1){
			*(*(augmented+g)+h) = *(mat_b+g);
		    }
		    else{
			*(*(augmented+g)+h) = *(*(mat_A+g)+h);
		    }
		}
	    }
	    /*at this point i have Ax=b in augmented matrix form to use rref on {augmented | b}  */
	    rref(augmented, nx*ny*nz,  nx*ny*nz+1);
	    /*now i need to extract the nx*ny*nz+1 column to cpy back to my old mat_b, mat_A*/
	    for(g = 0; g< nx*ny*nz ; g++){
		for(h = 0;h< nx*ny*nz ; h++){
		    if(h == nx*ny*nz-1){
			*(mat_newb+g)  = *(*(augmented+g)+h);
		    }
		 /*   else{
			*(*(mat_A+g)+h) = *(*(augmented+g)+h);
		    }*/
		}
	    }
	    tempb = mat_b;
	    mat_b = mat_newb;
	    mat_newb = tempb;

	}

	for(f = 0 ; f< nx ; f++){
	    for(g = 0 ; g < ny ; g++){
		for( h = 0 ; h < nz ; h++){
		    *(*(*(matrix + f) + g) +h) =  *(mat_b + f*ny*nz + g*nz + h);
		}
	    }
	}
	if(k == 1){
	    print3DMatrix(matrix, nx, ny, nz, dx, dx, dx, when);
	    when++;
	}
    }

}

void forwardTime(double ***matrix, double ***matrix_np1, int nx, int ny, int nz, double dx, double dy, double dz, double C, double q, int border, int NT){
    int i=0, j=0, k=0, l=0;
    double ***temp;
    double b0 = 0;
    double bn = 0;
    int when = 1;
    switch(border){
	case 0:
	    b0 = 0;
	    bn = 0;
	    break;
	case 1:
	    b0 = 3;
	    bn = 3;
	    break;
	case 2:
	    b0 = nx;
	    bn = 0;
	    break;
	default :
	    b0 = 0;
	    bn = 0;
    }

    double up = 0.0, down = 0.0, left = 0.0, right = 0.0, front = 0.0, back = 0.0;

    for(l=0;l<NT;l++){
	temp = matrix;
	matrix = matrix_np1;
	matrix_np1 = temp;
	for(i = 0 ; i< nx ; i++){
	    for(j = 0 ; j < ny ; j++){
		for(k = 0 ; k < nz ; k++){
		    if(i==0){
			up = b0;
		    }
		    else{
			up = matrix[i-1][j][k];
		    }
		    if(i==nx-1){
			down = bn;
		    }
		    else{
			down = matrix[i+1][j][k];
		    }
		    if(j==0){
			front = b0;
		    }
		    else{
			front = matrix[i][j-1][k];
		    }
		    if(j == ny-1){
			back = bn;
		    }
		    else{
			back = matrix[i][j+1][k];
		    }
		    if(k==0){
			left = b0;
		    }
		    else{
			left = matrix[i][j][k-1];
		    }
		    if(k==nz-1){
			right = bn;
		    }
		    else{
			right = matrix[i][j][k+1];
		    }
		    matrix_np1[i][j][k] = matrix[i][j][k] + C * left + C*back +C* down -(C*6.0 * matrix[i][j][k]) + C*right + C*front + C*up + q;
		    if(l == 13 ||l==107 || l==777){
			print3DMatrix(matrix_np1, nx, ny, nz, dx, dy, dz, when);
			when++;
		    }
		}

	    }
	}
    }

}


double ***initGaussians(int nx, int ny, int nz, double dx, double dy, double dz, double noise){
    int i, j, l, m, n;
    double ***matrix = (double ***) malloc(sizeof(double**) * nx);
    double *zcoords =  (double*) malloc(sizeof(double) *nx * ny * nz);
    double cx2, cy2, cz2;
    double midx = nx*dx/2.0;/*phnyical size divided by 2*/
    double midy = ny*dy/2.0;
    double midz = nz*dz/2.0;
    double std2 = 0.005;/*std dev times 2*/
    double e = 0.0;
    /*malloc the space for the 3D array*/
    /* contiguous memory bitches */
    for(i=0;i<nx;i++){
	*(matrix+i) = malloc(ny*sizeof(double*));

	for(j=0;j<ny;j++){
	    *(*(matrix+i) + j) =/* (double*)malloc(nz*sizeof(double));*/zcoords+(i*nz*ny)+(j*nz);
	}
    }
    for(l = 0;l<nx;l++){
	cx2 = (l*dx-midx)*(l*dx-midx);
	for(m=0;m<ny;m++){
	    cy2 = (m*dy-midy)*(m*dy-midy);
	    for(n=0;n<nz;n++){
		cz2 = (n*dz-midz)*(n*dz-midz);
		*(*(*(matrix+l)+m)+n) = exp(-cx2/std2-cy2/std2-cz2/std2);/*gaussian*/
	    }
	}
    }
    return matrix;
}

void freeMats( double ***matrix, double ***matrix_np1, int nx, int ny){
    int i;
    /*free my mallocs forever*/
    for(i = 0; i < nx; i++)
    {
	free(matrix[i]);
	free(matrix_np1[i]);
    }
    free (matrix);
    free (matrix_np1);

}


double distanceFormula(double *mat_old, double *mat_new, int length){
    int i;
    double sum = 0.0;
    for(i=0;i<length;i++){
	sum += (*(mat_old+i) - *(mat_new+i)) * (*(mat_old+i) - *(mat_new+i));
    }	    
    sum = sqrt(sum);
    return sum;
}


void rref(double **matrix, int n, int m){
    int i, j;
    /*larger values of potential pivot on top*/
    setMatrix(matrix, n, m);
    int col = 0;
    int row = 0;
    for(row=0;row<n;row++){
	while(*(*(matrix+row)+col) == 0 && col<m){
	    col++;
	}
	onesPivot(matrix, row, col, m);
	for(j = row+1;j<n;j++){
	    if(*(*(matrix + j)+col) !=0){
		multAdd(matrix, row, j, col, m);
	    }
	}
	if(row != 0){
	    for(j = row;j>0;j--){
		multAdd(matrix, row, j-1, col, m);
	    }
	}
	col++;
    }

}

void setMatrix(double **A, int n, int m){
    int i, j;
    for(i = 0;i<n-1;i++){
	for(j=i+1;j<n;j++){
	    if(*(*(A+i)) < *(*(A+j))){
		rowSwap(A, i, j, n, m);
	    }
	}
    }

}

void onesPivot(double **A, int row, int col, int m){
    int i, j;
    double pivot = *(*(A+row)+col);
    for(i=col;i<m;i++){
	*(*(A+row)+i) = (double)(*(*(A+row)+i) / pivot);/*this will divide all by the pivot*/
    }
}

void multAdd(double **A, int row, int row2, int col,  int m){
    int i;
    double eliminate = *(*(A+row2)+col);

    for(i=col;i<m;i++){
	*(*(A+row2)+i) = (*(*(A+row)+i)*eliminate*-1) + *(*(A+row2)+i);
    }
}

void rowSwap(double **A, int swap, int with, int n, int m){
    int i, j;
    double *temp = (double *) malloc (sizeof(double)*m);
    temp = *(A+swap);
    *(A+swap) = *(A+with);
    *(A+with) = temp;
}


void print3DMatrix(double ***matrix, int nx, int ny, int nz, double dx, double dy, double dz, int when){
    int l, m, n;
    FILE * init;
    FILE * init1;
    FILE * init2;
    FILE * init3;
    for(l = 0;l<nx;l++){
	for(m=0;m<ny;m++){
	    for(n=0;n<nz;n++){
		if(n == nz/2){/*slice at arbitrary (middle) z coordinate*/
		    switch(when){
			case 0:
			    init = fopen ("ftcs0.txt", "a");
			    fprintf(init, "%lf %lf %lf %c",l*dx, m*dy, *(*(*(matrix+l)+m)+n),'\n' );
			    fclose(init);
			    break;
			case 1:
			    init1 = fopen ("ftcs1.txt", "a");
			    fprintf(init1, "%lf %lf %lf %c",l*dx, m*dy, *(*(*(matrix+l)+m)+n),'\n' );
			    fclose(init1);
			    break;

			case 2:
			    init2 = fopen ("ftcs2.txt", "a");
			    fprintf(init2, "%lf %lf %lf %c",l*dx, m*dy, *(*(*(matrix+l)+m)+n),'\n' );
			    fclose(init2);
			    break;
			case 3:
			    init3 = fopen ("ftcs3.txt", "a");
			    fprintf(init3, "%lf %lf %lf %c",l*dx, m*dy, *(*(*(matrix+l)+m)+n),'\n' );
			    fclose(init3);
			    break;

		    }
		    //		    printf("%f %f %f \n", l*dx, m*dy,*(*(*(matrix+l)+m)+n));
		}
	    }
	}
    }

}


void welcome(int *nx, int *ny, int *nz, double *alpha, double *q, double *dt, double *dx, int *algo, int *border ){
    printf("\n\n\n\n\n********************************************\n");
    printf("* FTCS - Crank-Nicholson - ADI Heat Solver *\n");
    printf("********************************************\n");
    printf("* The Initial Conditions Have Been Set To: *\n");
    printf("nx: %d\n",*nx);
    printf("ny: %d\n",*ny);
    printf("nz: %d\n",*nz);
    printf("alpha: %f\n", *alpha);
    printf("q: %f\n",*q);
    printf("dt: %f\n",*dt);
    printf("dx=dy=dx: %f\n\n",*dx);
    int *answer = malloc(sizeof(int));;

    printf("To change these values by passing in a file, type 1 (return), otherwise type 0 (return): \n");
    scanf("%d", answer);
    if(*answer ==1){
	printf("USAGE: int nx int ny int nz double dt double alpha double q int method int boundary\n");
	printf("Please enter the complete file name> \n");
	/*create file pointer, bring in file with inits via argv
	 *          *      */
	char *name = malloc(sizeof(char)*255);
	scanf("%s",name);
	printf("%s",name);
	FILE *fp;
	fp  = fopen(name, "w+");
	fscanf(fp, "%d %d %d %lf %lf %lf %d %d", nx, ny, nz, dt, alpha,q, border, algo);
	fclose(fp);

    }
    printf("Select Method\n0 FTCS\n1 Crank Nicholson\n2 ADI\n> ");
    scanf("%d", algo);
    printf("Select Boundary Conditions\n0 Zeroes\n1 Constant (3)\n2 Periodic\n > ");
    scanf("%d", border);
}


void rref2(double **A, int m, int n){
    double* tmp=dvector(1,n);
    double tol=0.000000000000001;
    double p;
    long i,j,k,l;
    i=1;
    j=1;
    while ((i <= m) && (j<= n)){
	/* Initialize pivot */
	p=0.; k=i;
	/* Extract best pivot, and it's row index 
	 *        (even non-zero entries.  Called partial pivoting -- improves numerical stability)
	 *            */
	for(l=i;l<=m;l++){
	    if ( p < fabs( A[l][j] ) ){
		p=A[l][j];
		k=l;
	    }
	}
	if (fabs(p)<=tol){
	    for(l=i;l<=m;l++){
		A[l][j]=0.;
	    }
	    j++;
	    continue;
	}
	if(k!=i) /* Swap ith and kth rows if needed */
	    for(l=j;l<=n;l++) SWAP( A[i][l] , A[k][l] );

	/* Divide pivot row by pivot element */
	p=A[i][j];
	for (l=j;l<=n;l++) A[i][l] = A[i][l]/p;

	/* Subtract multiples of the pivot row from all
	 *        other rows  */
	for (k=1;k<=m;k++){
	    if (k!=i){
		for (l=j;l<=n;l++)
		    tmp[l]=A[k][j]*A[i][l];
		for (l=j;l<=n;l++)
		    A[k][l]=A[k][l]-tmp[l];
	    }
	}
	if (DEBUG)
	    i++;
	j++;
    }
}

double *dvector(long nl, long nh)
    /* allocate a double vector with subscript range v[nl..nh] */
{
    double *v;

    v=(double *)malloc((size_t) ((nh-nl+1+1)*sizeof(double)));
    if (!v) printf("allocation failure in dvector()");
    return v-nl+1;
}
