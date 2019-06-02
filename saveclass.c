#include "global.h"

void saveclass(struct parameter pm) {
	int i, j;
	FILE * fp = fopen("svm.csv", "w+");

	if (fp == NULL) {
		printf("load error!\n");
		return;
	}

	fprintf(fp, "%d,", pm.data_col);
	fprintf(fp, "%d,", pm.data_row);
	fprintf(fp, "%d,", pm.num_m);
	fprintf(fp, "%d,", pm.num_e);
	fprintf(fp, "%d,", pm.num_r);
	fprintf(fp, "%d,", pm.num_u);
	fprintf(fp, "%d,", pm.num_uind);

	fprintf(fp, "%d,", pm.kernel_evals);
	fprintf(fp, "%d,", pm.perturbations);
	fprintf(fp, "%lf,", pm.initial_C);
	fprintf(fp, "%lf\n",pm.scale);

	for (i = 0; i < pm.data_col-1; i++) {
		fprintf(fp, "%lf,", pm.a[i]);
	}
	fprintf(fp, "%lf\n", pm.a[pm.data_col-1]);
	fprintf(fp, "%lf\n", pm.b);

	for (i = 0; i < pm.data_col - 1; i++) {
		fprintf(fp, "%lf,", pm.g[i]);
	}
	fprintf(fp, "%lf\n", pm.g[pm.data_col - 1]);


	for (j = 0; j < pm.num_m; j++) {
		if (j == pm.num_m - 1) {
			fprintf(fp, "%d\n", pm.ind[1][j]);
		}
		else {
			fprintf(fp, "%d,", pm.ind[1][j]);
		}
	}
	for (j = 0; j < pm.num_e; j++) {
		if (j == pm.num_e - 1) {
			fprintf(fp, "%d\n", pm.ind[2][j]);
		}
		else {
			fprintf(fp, "%d,", pm.ind[2][j]);
		}
	}
	for (j = 0; j < pm.num_r; j++) {
		if (j == pm.num_r - 1) {
			fprintf(fp, "%d\n", pm.ind[3][j]);
		}
		else {
			fprintf(fp, "%d,", pm.ind[3][j]);
		}
	}
	for (j = 0; j < pm.num_u; j++) {
		if (j == pm.num_u - 1) {
			fprintf(fp, "%d\n", pm.ind[4][j]);
		}
		else {
			fprintf(fp, "%d,", pm.ind[4][j]);
		}
	}

	for (j = 0; j < pm.num_uind; j++) {
		if (j == pm.num_uind - 1) {
			fprintf(fp, "%d\n", pm.uind[j]);
		}
		else {
			fprintf(fp, "%d,", pm.uind[j]);
		}
	}

	for (i = 0; i < pm.num_m + 1; i++) {
		for (j = 0; j < pm.data_col; j++) {
			if (j == pm.data_col - 1) {
				fprintf(fp, "%lf\n", pm.Q[i][j]);
			}
			else {
				fprintf(fp, "%lf,", pm.Q[i][j]);
			}
		}
	}

	for (i = 0; i < pm.num_m + 1; i++) {
		for (j = 0; j < pm.num_m + 1; j++) {
			if (j == pm.num_m) {
				fprintf(fp, "%lf\n", pm.Rs[i][j]);
			}
			else {
				fprintf(fp, "%lf,", pm.Rs[i][j]);
			}
		}
	}
	for (i = 0; i < pm.data_col; i++) {
		for (j = 0; j < pm.data_row+1; j++) {
			if (j == pm.data_row) {
				fprintf(fp, "%d\n", pm.y[i]);
			}
			else {
				fprintf(fp, "%lf,", pm.X[j][i]);
			}
		}
	}
	fclose(fp);
}
