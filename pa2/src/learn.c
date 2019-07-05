#include <stdio.h>
#include <stdlib.h>
#include <math.h>

	//function to malloc all double arrays created
double** Malloc2DDoubleArr(int arrSizeX, int arrSizeY) {
	double** ArrMal;
	int i;
	ArrMal = (double**) malloc(arrSizeX*sizeof(double*));
	for (i = 0; i < arrSizeX; i++){
		ArrMal[i] = (double*) malloc(arrSizeY*sizeof(double));
	}

	return ArrMal;
} 




int main (int argc, char** argv) {


	int K, N;
	FILE *inputFile = fopen(argv[1], "r");
	if (inputFile == NULL) {
		printf("error\n");
		return 0;
	}


	//first string in file
	char * commandStr = (char*) malloc( sizeof(char) * 10);

	fscanf(inputFile,"%s\n", commandStr);
	fscanf(inputFile, "%d", &K);
	fscanf(inputFile, "%d", &N);




	double **original = Malloc2DDoubleArr(N, K+1);
	double **matrix = Malloc2DDoubleArr(N, K+1);
	double **matrixT = Malloc2DDoubleArr(K+1, N);
	double **Y = Malloc2DDoubleArr(N, 1);
	double **W = Malloc2DDoubleArr(K+1, 1);

	double **mmT = Malloc2DDoubleArr(N, N);
	double **mmTAug = Malloc2DDoubleArr(N, (N*2)+1);
	double **inverse = Malloc2DDoubleArr(N, N);
	double **invmT = Malloc2DDoubleArr(K+1, N);







	//put numbers in matrix
	int row, col;
	for (row = 0; row < N; row++) {
		for (col = 0; col < K + 1; col++){
			fscanf(inputFile, "%lf%*[,]", &original[row][col]);     	
		}	
	}

	for (row =0;row<N;row++){
		matrix[row][0]=1;
		Y[row][0]=original[row][K];
	}
	for (row = 0; row < N; row++) {
		for (col = 0; col < K; col++){
			matrix[row][col+1]=original[row][col];
		}
	}

	//print y
	/* for (row = 0; row < N; row++) { */
	/*     printf("Y\n"); */
	/*     printf("-%lf", Y[row][0]); */
	/* 	   } */

	//print matrix
	/* for (row = 0; row < N; row++) { */
	/*   printf("\n"); */
	/*   for (col = 0; col < K + 1; col++){ */
	/*     printf("-%lf", matrix[row][col]); */
	/*   } */
	/* } */



	//transpose matrix to matrixT
	for (row = 0; row < N; row++) {
		for (col = 0; col < K+1; col++){
			matrixT[col][row]=matrix[row][col];
		}
	}



	//print matrix trans
	/* for (row = 0; row < K+1; row++) { */
	/*   printf("\n"); */
	/*   for (col = 0; col < N; col++){ */
	/*     printf("-%9.6f", matrixT[row][col]); */
	/*   } */
	/* } */





	//multiply matrix by matrixT to get mmT
	int i,j,p;

	for(i=0;i<K+1;i++){
		for(j=0;j<K+1;j++){
			mmT[i][j]=0;
			for(p=0;p<N;p++){
				mmT[i][j]+=matrixT[i][p]*matrix[p][j];
			}
		}
	}



	//print mmT
	/* printf("mmT\n"); */
	/* for (row = 0; row < K+1; row++) { */
	/*   printf("\n"); */
	/*   for (col = 0; col <K+1; col++){ */
	/*     printf("- %lf", mmT[row][col]); */
	/*   } */
	/* } */





	//augment mmT
	for (row = 0; row < K+1; row++) {
		for (col = 0; col < K+1; col++){
			mmTAug[row][col]=mmT[row][col];
		}
		for (col = K+2; col < ((K+1)*2)+1; col++){
			mmTAug[row][col]=0;
		}
	}
	col=K+1;
	for (row = 0; row < K+1; row++) {
		mmTAug[row][col]=1;
		col++;
	}




	//print augmented mmT
	/* printf("\n"); */
	/* printf("augmented matrix\n"); */
	/* for (row = 0; row < K+1; row++) { */
	/*   printf("\n"); */
	/*   for (col = 0; col <(K+1)*2; col++){ */
	/*     printf("-%9.6f", mmTAug[row][col]); */
	/*   } */
	/* } */





	//gauss jordan elimination to find inverse
	int x,y, tempR, tempC;
	double over=0;
	double temp;
	for(row=0;row<(K+1);row++){

		for(col=0;col<(K+1)*2;col++){


			if(row==col){
				if(mmTAug[row][col]==0){
					mmTAug[row][col]++;
				}
				temp = mmTAug[row][col];
				x=0;

				for(x=0;x<(K+1)*2;x++){
					mmTAug[row][x]/=temp;
				}

				tempR=row;
				tempC=0;

				for(tempR=row+1;tempR<(K+1);tempR++){
					over=mmTAug[tempR][col]/mmTAug[row][row];

					for(tempC=0;tempC<(K+1)*2;tempC++){
						mmTAug[tempR][tempC]=mmTAug[tempR][tempC]-(mmTAug[row][tempC]*over);
					}

				}

			}
		}
	}

	x=2;
	y=2;

	for(x=K;x>=0;x--){
		for(y=K;y>=0;y--){
			if(x!=y && mmTAug[x][y]!=0){
				double change = mmTAug[x][y];
				int column=0;

				for(column=0;column<(K+1)*2;column++){
					mmTAug[x][column]=mmTAug[x][column]-(change*mmTAug[y][column]);
				}
			}
		}
	}

	for(row=0;row<K+1;row++){
		for(col=0;col<K+1;col++){
			inverse[row][col]=mmTAug[row][col+K+1];
		}
	}

	//print  augmented mmT after gauss 
	/* printf("\n"); */
	/* printf("augmented after gauss\n"); */
	/* for (row = 0; row < K+1; row++) { */
	/*   printf("\n"); */
	/*   for (col = 0; col <(K+1)*2; col++){ */
	/*     printf("-%lf", mmTAug[row][col]); */
	/*   } */
	/* } */

	//print inverse
	/* 	printf("\n"); */
	/* printf("inverse\n"); */
	/* for (row = 0; row < K+1; row++) { */
	/*   printf("\n"); */
	/*   for (col = 0; col <K+1; col++){ */
	/*     printf("- %lf", inverse[row][col]); */
	/*   } */
	/* } */




	//multiply inverse by matrixT
	for(i=0;i<K+1;i++){
		for(j=0;j<N;j++){
			invmT[i][j]=0;
			for(p=0;p<K+1;p++){
				invmT[i][j]+=inverse[i][p]*matrixT[p][j];
			}
		}
	}



	//print invmT
	/* printf("\n"); */
	/* printf("invmt\n"); */
	/* for (row = 0; row < K+1; row++) { */
	/*   printf("\n"); */
	/*   for (col = 0; col <N; col++){ */
	/*     printf("-%lf", invmT[row][col]); */
	/*   } */
	/* } */




	//multiply previous result by prices to find weight
	for(i=0;i<K+1;i++){
		for(j=0;j<1;j++){
			W[i][j]=0;
			for(p=0;p<N;p++){
				W[i][j]+=invmT[i][p]*Y[p][j];
			}
		}
	}



	//print weight
	/* printf("\n"); */
	/* printf("weight\n"); */
	/* for (row = 0; row < K+1; row++) { */
	/*   printf("\n"); */
	/*     printf("-%lf", W[row][0]); */
	/*   } */



	//scan last file
	int examples, M;

	FILE *testFile = fopen(argv[2], "r");
	if (testFile == NULL) {
		printf("error\n");
		return 0;
	}
	char * sCommandStr = (char*) malloc( sizeof(char) * 10);

	fscanf(testFile,"%s\n", sCommandStr);
	fscanf(testFile, "%d \n", &M);
	fscanf(testFile, "%d \n", &examples);



	double** X= Malloc2DDoubleArr(examples, M+1);
	double** price= Malloc2DDoubleArr(examples, 1);



	for (row = 0; row < examples; row++) {
		X[row][0]=1;
		for (col = 1; col < M + 1; col++){
			fscanf(testFile, "%lf%*[,]", &X[row][col]);     	
		}	
	}




	//multiply X by W to find Y
	for(i=0;i<examples;i++){
		for(j=0;j<1;j++){
			price[i][j]=0;
			for(p=0;p<M+1;p++){
				price[i][j]+=X[i][p]*W[p][j];
			}
		}
	}

	//print result
	for(row=0;row<examples;row++){
		printf("%0.0lf\n", price[row][0]);
	}	


	//freeing allocated memory
	free(original);
	free(matrix);
	free(matrixT);
	free(Y);
	free(W);
	free(mmT);
	free(mmTAug);
	free(inverse);
	free(invmT);
	free(X);
	free(price);
	free(commandStr);
	free(sCommandStr);



	return 0;

}

