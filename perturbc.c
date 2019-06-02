#include "global.h"

void perturbc(double C_new,struct parameter * pm) {
	int kernel_evals_begin = (*pm).kernel_evals;
	double * lambda = (double*)malloc(sizeof(double)*(*pm).data_col);
	double * SQL = (double*)malloc(sizeof(double)*(*pm).data_col);
	double p = 1;
	double Syl = 0;
	double disp_p_delta = 0.2;
	double disp_p_count = 1;
	int i, j,num_MVs;
	char s[100];
	(*pm).num_u = 0;
	for (i = 0; i < (*pm).data_col; i++) {
		lambda[i] = C_new - (*pm).C[i];
	} // lambda 계산

	if ((*pm).num_e == 0) {
		int check = 0;
		int * inde = (int *)malloc(sizeof(int)*(*pm).data_col);
		double * delta_p;

		for (j = 0; j < (*pm).data_col; j++) {
			if (lambda[j] != 0) {
				inde[check++] = j;
			} // index번호는 inde배열에, check는 개수
		}
		delta_p = (double*)malloc(sizeof(double)*check);
		for (j = 0; j < check; j++) {
			delta_p[j] = ((*pm).a[inde[j]] - (*pm).C[inde[j]]) / lambda[inde[j]];

			if (delta_p[j] > 0 && delta_p[j] < p) {
				p = delta_p[j];
				i = j;
			}
		}
		for (j = 0; j < (*pm).data_col; j++) {
			(*pm).C[j] += lambda[j] * p;
		} // C값 결정하기

		if (p < 1) {
			int indco, removed_i;
			bookkeeping(inde[i], MARGIN, ERROR, &indco, &removed_i,pm);
			updateRQ(0,0,indco,0,pm);
		}

		free(delta_p);
		free(inde);
	}
	else {
		p = 0;
	}
	if (p < 1) {
		double ** temp = (double**)malloc(sizeof(double*)*(*pm).data_col);
		int * index = (int*)malloc(sizeof(int)*(*pm).data_col);
		for (i = 0; i < (*pm).data_col; i++) {
			temp[i] = (double*)malloc(sizeof(double)*(*pm).num_e);
			index[i] = i;
			SQL[i] = 0;
		}
		kernel(index, (*pm).ind[ERROR], (*pm).data_col, (*pm).num_e, temp,pm);
		for (i = 0; i < (*pm).data_col; i++) {
			for (j = 0; j < (*pm).num_e; j++) {
				temp[i][j] *= (*pm).y[i] * (*pm).y[(*pm).ind[ERROR][j]];
				SQL[i] += temp[i][j] * lambda[(*pm).ind[ERROR][j]];
			}
		}
		for (i = 0; i < (*pm).num_e; i++) {
			SQL[(*pm).ind[ERROR][i]] += (*pm).deps * lambda[(*pm).ind[ERROR][i]];
			Syl += (*pm).y[(*pm).ind[ERROR][i]] * lambda[(*pm).ind[ERROR][i]];
		}

		free(index);
		for (i = 0; i < (*pm).data_col; i++) {
			free(temp[i]);
		}
		free(temp);
	}

	sprintf(s, "p = %.2f", p);
	printf("%s\n", s);
	num_MVs = (*pm).num_m;
	(*pm).perturbations = 0;

	while (p < 1) {
		double * beta = (double*)malloc(sizeof(double)*(num_MVs + 1));
		double * gamma = (double*)malloc(sizeof(double)*((*pm).data_col));
		memset(beta, 0, sizeof(double)*(num_MVs + 1));
		memset(gamma, 0, sizeof(double)*((*pm).data_col));

		(*pm).perturbations += 1;

		if (num_MVs > 0) {
			double * v = (double*)malloc(sizeof(double)*(num_MVs + 1));
			memset(v, 0, sizeof(double)*(num_MVs + 1));
			int * ind_temp = (int*)malloc(sizeof(int)*((*pm).num_e + (*pm).num_u + (*pm).num_r));
			double sum = 0;

			if (p < 1 - 1.000e-03) {
				v[0] = -1 * Syl;
				for (i = 0; i < (*pm).data_col; i++) {
					sum += (*pm).y[i] * (*pm).a[i];
				}
				sum = sum / (1 - p);
				v[0] = v[0] - sum;
			}
			else {
				v[0] = -1 * Syl;
			}

			for (i = 1; i < num_MVs + 1; i++) {
				v[i] = -1 * SQL[(*pm).ind[MARGIN][i - 1]];
			}
			for (i = 0; i < num_MVs + 1; i++) {
				for (j = 0; j < num_MVs + 1; j++) {
					beta[i] += (*pm).Rs[i][j] * v[j];
				}
			}

			for (i = 0; i < (*pm).num_e + (*pm).num_r + (*pm).num_u; i++) {
				if (i < (*pm).num_e) {
					ind_temp[i] = (*pm).ind[ERROR][i];
				}
				else if ((i < (*pm).num_e + (*pm).num_r) && (i >= (*pm).num_e)) {
					ind_temp[i] = (*pm).ind[RESERVE][i - (*pm).num_e];
				}
				else {
					ind_temp[i] = (*pm).ind[UNLEARNED][i - (*pm).num_e - (*pm).num_r];
				}
			}
			if ((*pm).num_e + (*pm).num_r + (*pm).num_u > 0) {
				for (i = 0; i < (*pm).num_e + (*pm).num_r + (*pm).num_u; i++) {
					for (j = 0; j < num_MVs + 1; j++) {
						gamma[ind_temp[i]] += (*pm).Q[j][ind_temp[i]] * beta[j];
					}
					gamma[ind_temp[i]] += SQL[ind_temp[i]];
				}
			}
			free(v);
			free(ind_temp);
		}
		else {
			beta[0] = 0;
			for (i = 0; i < (*pm).data_col; i++) {
				gamma[i] = SQL[i];
			}
		}
		double min_delta_p;
		int indss, cstatus, nstatus, indco, removed_i;
		min_delta_p_c(p, gamma, beta, lambda, &min_delta_p, &indss, &cstatus, &nstatus,pm);

		if ((*pm).num_e > 0) {
			for (i = 0; i < (*pm).num_e; i++) {
				(*pm).a[(*pm).ind[ERROR][i]] += lambda[(*pm).ind[ERROR][i]] * min_delta_p;
			}
		}
		if (num_MVs > 0) {
			for (i = 0; i < num_MVs; i++) {
				(*pm).a[(*pm).ind[MARGIN][i]] += beta[i + 1] * min_delta_p;
			}
		}
		(*pm).b += beta[0] * min_delta_p;
		for (i = 0; i < (*pm).data_col; i++) {
			(*pm).g[i] += gamma[i] * min_delta_p;
			(*pm).C[i] += lambda[i] * min_delta_p;
		}
		p += min_delta_p;
		bookkeeping(indss, cstatus, nstatus, &indco, &removed_i,pm);
		
		if ((cstatus == MARGIN) && (nstatus == ERROR)) {
			for (i = 0; i < (*pm).data_col; i++) {
				SQL[i] += (*pm).Q[indco][i] * lambda[indss];
			}
			Syl += (*pm).y[indss] * lambda[indss];
		}
		for (i = 0; i < (*pm).num_m; i++) {
			(*pm).g[(*pm).ind[MARGIN][i]] = 0;
		}
		if (nstatus == MARGIN) {
			num_MVs += 1;
			if (num_MVs > 1) {
				int * indc_s = (int *)malloc(sizeof(int) * 1);
				double ** temps = (double**)malloc(sizeof(double*)*1);
				temps[0] = (double*)malloc(sizeof(double) * 1);
				indc_s[0] = indss;

				for (i = 0; i < num_MVs; i++) {
					beta[i] = 0;
					for (j = 0; j < num_MVs; j++) {
						beta[i] += (*pm).Rs[i][j] * (*pm).Q[j][indss];
					}
					beta[i] *= -1;
				}
				kernel(indc_s, indc_s, 1, 1, temps,pm);
				free(gamma);

				gamma = (double *)malloc(sizeof(double));
				gamma[0] = temps[0][0] + (*pm).deps;
				for (i = 0; i < num_MVs; i++) {
					gamma[0] += (*pm).Q[i][indss] * beta[i];
				}


				free(indc_s);
				free(temps[0]);
				free(temps);
			}
			updateRQ(beta, gamma, indss, 1,pm);
		}
		else if (cstatus == MARGIN) {
			num_MVs -= 1;
			updateRQ(beta, gamma, indco, 0,pm);
		}
		if ((cstatus == ERROR) && (nstatus == MARGIN)) {
			for (i = 0; i < (*pm).data_col; i++) {
				SQL[i] = SQL[i] - (*pm).Q[num_MVs][i] * lambda[indss];
			}
			Syl = Syl - (*pm).y[indss] * lambda[indss];
		}

		if (p >= disp_p_delta * disp_p_count) {
			disp_p_count = disp_p_count + 1;
			sprintf(s, "p = %.2f", p);
			printf("%s\n", s);
		}
		free(beta);
		free(gamma);
	}
	printf("Perturbations complete!\n");

	printf("Margin vectors: \t\t%d\n", (*pm).num_m);
	printf("Error vectors: \t\t%d\n", (*pm).num_e);
	printf("Reserve vectors: \t%d\n", (*pm).num_r);
	printf("Kernel evaluations: \t%d\n", (*pm).kernel_evals - kernel_evals_begin);
	free(lambda);
	free(SQL);
}