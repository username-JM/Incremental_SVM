#include "global.h"

void loadclass(struct parameter * pm,int num) {
	int i, j;
	FILE * fp = fopen("svm.csv", "r");
	if (fp == NULL) {
		printf("load error!\n");
		return;
	}

	fscanf(fp, "%d,", &((*pm).data_col));
	fscanf(fp, "%d,", &((*pm).data_row));
	fscanf(fp, "%d,", &((*pm).num_m));
	fscanf(fp, "%d,", &((*pm).num_e));
	fscanf(fp, "%d,", &((*pm).num_r));
	fscanf(fp, "%d,", &((*pm).num_u));
	fscanf(fp, "%d,", &((*pm).num_uind));
	fscanf(fp, "%d,", &((*pm).kernel_evals));
	fscanf(fp, "%d,", &((*pm).perturbations));

	fscanf(fp, "%lf,", &((*pm).initial_C));
	fscanf(fp, "%lf\n", &((*pm).scale));
	
	(*pm).a = (double*)malloc(sizeof(double)*((*pm).data_col+num));
	(*pm).g = (double*)malloc(sizeof(double)*((*pm).data_col+num));
	(*pm).C = (double*)malloc(sizeof(double)*((*pm).data_col + num));
	(*pm).uind = (int*)malloc(sizeof(int)*((*pm).data_col + num));
	(*pm).y = (int*)malloc(sizeof(int)*((*pm).data_col + num));
	(*pm).Q = (double**)malloc(sizeof(double*) * ((*pm).num_m+1));
	(*pm).Rs = (double**)malloc(sizeof(double*) * ((*pm).num_m+1));
	(*pm).X = (double**)malloc(sizeof(double*)*(*pm).data_row);

	for (i = 1; i < 5; i++) {
		(*pm).ind[i] = (int*)malloc(sizeof(int)*((*pm).data_col + num));
	}

	for (i = 0; i < (*pm).data_row; i++) {
		(*pm).X[i] = (double*)malloc(sizeof(double)*((*pm).data_col + num));
	}

	for (i = 0; i < (*pm).num_m + 1; i++) {
		(*pm).Q[i] = (double*)malloc(sizeof(double)*((*pm).data_col + num));
		(*pm).Rs[i] = (double*)malloc(sizeof(double)*((*pm).num_m + 1));
	}

	for (i = 0; i < (*pm).data_col; i++) {
		if (i == (*pm).data_col-1) {
			fscanf(fp, "%lf\n", &((*pm).a[i]));
		}
		else {
			fscanf(fp, "%lf,", &((*pm).a[i]));
		}
	}
	fscanf(fp, "%lf\n", &((*pm).b));

	for (i = 0; i < (*pm).data_col; i++) {
		if (i == (*pm).data_col-1) {
			fscanf(fp, "%lf\n", &((*pm).g[i]));
		}
		else {
			fscanf(fp, "%lf,", &((*pm).g[i]));
		}
	}

	for (i = 0; i < (*pm).num_m; i++) {
		if (i == (*pm).num_m-1) {
			fscanf(fp, "%d\n", &((*pm).ind[1][i]));
		}
		else {
			fscanf(fp, "%d,", &((*pm).ind[1][i]));
		}
	}

	for (i = 0; i < (*pm).num_e; i++) {
		if (i == (*pm).num_e - 1) {
			fscanf(fp, "%d\n", &((*pm).ind[2][i]));
		}
		else {
			fscanf(fp, "%d,", &((*pm).ind[2][i]));
		}
	}

	for (i = 0; i < (*pm).num_r; i++) {
		if (i == (*pm).num_r - 1) {
			fscanf(fp, "%d\n", &((*pm).ind[3][i]));
		}
		else {
			fscanf(fp, "%d,", &((*pm).ind[3][i]));
		}
	}

	for (i = 0; i < (*pm).num_u; i++) {
		if (i == (*pm).num_u - 1) {
			fscanf(fp, "%d\n", &((*pm).ind[4][i]));
		}
		else {
			fscanf(fp, "%d,", &((*pm).ind[4][i]));
		}
	}

	for (i = 0; i < (*pm).num_uind; i++) {
		if (i == (*pm).num_uind - 1) {
			fscanf(fp, "%d\n", &((*pm).uind[i]));
		}
		else {
			fscanf(fp, "%d,", &((*pm).uind[i]));
		}
	}

	for (i = 0; i < (*pm).num_m + 1; i++) {
		for (j = 0; j < (*pm).data_col; j++) {
			if (j == (*pm).data_col-1) {
				fscanf(fp, "%lf\n", &((*pm).Q[i][j]));
			}
			else {
				fscanf(fp, "%lf,", &((*pm).Q[i][j]));
			}
		}
	}

	for (i = 0; i < (*pm).num_m + 1; i++) {
		for (j = 0; j < (*pm).num_m + 1; j++) {
			if (j == (*pm).num_m) {
				fscanf(fp, "%lf\n", &((*pm).Rs[i][j]));
			}
			else {
				fscanf(fp, "%lf,", &((*pm).Rs[i][j]));
			}
		}
	}

	for (i = 0; i < (*pm).data_col; i++) {
		for (j = 0; j < (*pm).data_row + 1; j++) {
			if (j == (*pm).data_row) {
				fscanf(fp, "%d\n", &((*pm).y[i]));
			}
			else {
				fscanf(fp, "%lf,", &((*pm).X[j][i]));
			}
		}
	}
	fclose(fp);
}