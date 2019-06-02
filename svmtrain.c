#include "global.h"

void svmtrain(double *X_new[], int y_new[], int num, int new_model) {
	struct parameter pm;
	int i, j;
	pm.max_reserve_vectors = 3000;
	pm.deps = 1e-3;
	pm.type = 5;

	if (new_model) {
		FILE * fp = fopen("result2.csv", "r");
		if (fp == NULL) {
			printf("open error!\n");
			return;
		}
		pm.initial_C = 10; //setting parameter C 
		pm.b = 0;
		pm.scale = 2;  // setting kernel scale
		pm.kernel_evals = 0;
		pm.perturbations = 0;
		pm.data_col = 70; //num of data
		pm.data_row = 12; //num of feature
		
		pm.num_m = 0;
		pm.num_e = 0;
		pm.num_r = 0;
		pm.num_u = pm.data_col;
		pm.num_uind = 0;

		pm.a = (double*)malloc(sizeof(double)*pm.data_col);
		pm.g = (double*)malloc(sizeof(double)*pm.data_col);
		pm.C = (double*)malloc(sizeof(double)*pm.data_col);
		pm.uind = (int*)malloc(sizeof(int)*pm.data_col);
		pm.y = (int*)malloc(sizeof(int)*pm.data_col);
		pm.Q = (double**)malloc(sizeof(double*) * 1);
		pm.Rs = (double**)malloc(sizeof(double*) * 2);
		pm.X = (double**)malloc(sizeof(double*)*pm.data_row);

		for (i = 1; i < 5; i++) {
			pm.ind[i] = (int*)malloc(sizeof(int)*pm.data_col);
		}
		for (i = 0; i < 2; i++) {
			pm.Rs[i] = (double*)malloc(sizeof(double) * 2);
		}
		pm.Q[0] = (double*)malloc(sizeof(double)*pm.data_col);

		for (i = 0; i < pm.data_row; i++) {
			pm.X[i] = (double*)malloc(sizeof(double)*pm.data_col);
		}

		for (i = 0; i < pm.data_col; i++) {
			pm.a[i] = 0;
			pm.g[i] = 0;
			pm.C[i] = pm.initial_C;

			for (j = 0; j < pm.data_row + 1; j++) {
				if (j == pm.data_row) {
					fscanf(fp, "%d\n", &pm.y[i]);
				}
				else {
					fscanf(fp, "%lf,", &pm.X[j][i]);
				}
			}
			pm.Q[0][i] = pm.y[i];
			pm.ind[4][i] = i;
		}

		pm.Rs[0][0] = INF;

		while (pm.num_u != 0) {
			learn(pm.ind[4][0], 1,&pm);
		}

		printf("\nMARGIN : %d\n", pm.num_m);
		printf("ERROR : %d\n",pm.num_e);
		printf("RESERVE : %d\n", pm.num_r);
		printf("Kernel evaluations: %d\n", pm.kernel_evals);

		//perturbk(1.5,&pm);
		//perturbc(30, &pm);
		fclose(fp);
		saveclass(pm);

		for (i = 0; i < pm.data_col; i++) {
			printf("%.4f\n", pm.a[i]);
		}
		free(pm.a);
		free(pm.C);
		free(pm.g);
		free(pm.uind);
		free(pm.y);
		
		for (i = 0; i < pm.num_m + 1; i++) {
			free(pm.Q[i]);
			free(pm.Rs[i]);
		}

		for (i = 0; i < pm.data_row; i++) {
			free(pm.X[i]);
		}
		for (i = 1; i < 5; i++) {
			free(pm.ind[i]);
		}
		free(pm.Q);
		free(pm.Rs);
		free(pm.X);
	}
	else {
		loadclass(&pm,num);

		double ** temp = (double**)malloc(sizeof(double*)*pm.num_m);
		int *indc_s = (int*)malloc(sizeof(int)*num);
		for (i = 0; i < pm.num_m; i++) {
			temp[i] = (double*)malloc(sizeof(double)*num);
		}
		pm.data_col += num;
		
		for (i = 0; i < num; i++) {
			pm.ind[UNLEARNED][pm.num_u++] = pm.data_col - num + i;
			indc_s[i] = pm.data_col - num + i;
		}

		for (i = 0; i < num; i++) {
			for (j = 0; j < pm.data_row; j++) {
				pm.X[j][pm.data_col - num + i] = X_new[j][i];
			}
			pm.y[pm.data_col - num + i] = y_new[i];
		}

		kernel(pm.ind[MARGIN], indc_s, pm.num_m, num, temp,&pm);

		for (i = 0; i < num; i++) {
			pm.Q[pm.data_col - num + i][0] = y_new[i];
			for (j = 1; j < pm.num_m+1; j++) {
				pm.Q[pm.data_col - num + i][j] = pm.y[pm.ind[MARGIN][j - 1]] * y_new[i] * temp[j - 1][i];
			}
		}

		for (i = 0; i < pm.num_m; i++) {
			free(temp[i]);
		}
		free(temp);
		free(indc_s);

		while (pm.num_u != 0) {
			//learn(pm.ind[4][0], 1);
			pm.num_u--;
		}
		saveclass(pm);

		free(pm.a);
		free(pm.C);
		free(pm.g);
		free(pm.uind);
		free(pm.y);

		for (i = 0; i < pm.num_m + 1; i++) {
			free(pm.Q[i]);
			free(pm.Rs[i]);
		}

		for (i = 0; i < pm.data_row; i++) {
			free(pm.X[i]);
		}
		for (i = 1; i < 5; i++) {
			free(pm.ind[i]);
		}
		free(pm.Q);
		free(pm.Rs);
		free(pm.X);
	}
}