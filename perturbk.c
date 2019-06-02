#include "global.h"

void perturbk(double new_scale,struct parameter * pm) {
	int number = 0;
	int kernel_evals_begin = (*pm).kernel_evals;
	int i, j, check, num_MVs, num_learned;
	int K_size = (*pm).num_m + (*pm).num_e;
	char s[100];
	int inda_size = (*pm).num_m + (*pm).num_e + (*pm).num_r;
	double sum;
	double **y_mul; 
	double ** temp;
	int * inda = (int *)malloc(sizeof(int)*((*pm).num_m + (*pm).num_e + (*pm).num_r));
	int * flag;
	int * index_i = (int*)malloc(sizeof(int)*(*pm).data_col);
	int * kernel_second = (int*)malloc(sizeof(int)*(*pm).data_col);
	double * SQL = (double*)malloc(sizeof(double)*(*pm).data_col);
	double * f = (double*)malloc(sizeof(double)*((*pm).num_m + (*pm).num_e + (*pm).num_r));
	double ** K = (double**)malloc(sizeof(double)*((*pm).num_m + (*pm).num_e));
	double* lambda = (double*)malloc(sizeof(double)*(*pm).data_col);
	double Syl = 0;
	double p_s;

	for (i = 0; i < (*pm).num_m + 1; i++) {
		free((*pm).Q[i]);
		free((*pm).Rs[i]);
	}
	free((*pm).Q);
	free((*pm).Rs);

	for (j = 0; j < (*pm).num_m + (*pm).num_e; j++) {
		K[j] = (double*)malloc(sizeof(double)*((*pm).num_m + (*pm).num_e + (*pm).num_r));
	}
	memset(SQL, 0, sizeof(double)*(*pm).data_col);
	memset(lambda, 0, sizeof(double)*(*pm).data_col);

	(*pm).scale = new_scale;
	for (j = 0; j < (*pm).num_r + (*pm).num_m + (*pm).num_e; j++) {
		if (j < (*pm).num_m) {
			inda[j] = (*pm).ind[MARGIN][j];
		}else if ((j >= (*pm).num_m) && (j < (*pm).num_m + (*pm).num_e)) {
			inda[j] = (*pm).ind[ERROR][j - (*pm).num_m];
		}
		else {
			inda[j] = (*pm).ind[RESERVE][j - (*pm).num_m - (*pm).num_e];
		}
	}
	svmeval2(inda, f, K,pm);
	for (j = 0; j < (*pm).num_m + (*pm).num_e + (*pm).num_r; j++) {
		(*pm).g[inda[j]] = (*pm).y[inda[j]] * f[j] + (*pm).a[inda[j]] * (*pm).deps - 1;
	}
	flag = (int*)malloc(sizeof(int)*(*pm).data_col);
	check = 0;

	for (j = 0; j < (*pm).num_e; j++) {
		flag[j] = ((*pm).g[(*pm).ind[ERROR][j]] >= 0);
		if (flag[j] == 1) {
			index_i[check++] = j;
		}
	}
	if (check > 0) {
		for (j = 0; j < check; j++) {
			(*pm).ind[UNLEARNED][(*pm).num_u] = (*pm).ind[ERROR][index_i[j]];
			(*pm).num_u++;
			lambda[(*pm).ind[ERROR][index_i[j]]] = -1 * (*pm).a[(*pm).ind[ERROR][index_i[j]]];
		}
		y_mul = (double**)malloc(sizeof(double*)*check);
		for (j = 0; j < check; j++) {
			y_mul[j] = (double*)malloc(sizeof(double)*inda_size);
		}
		for (i = 0; i < check; i++) {
			for (j = 0; j < inda_size; j++) {
				y_mul[i][j] = (*pm).y[(*pm).ind[ERROR][index_i[i]]] * (*pm).y[inda[j]];
				y_mul[i][j] *= K[(*pm).num_m+index_i[i]][j];
			}
		}
		for (i = 0; i < inda_size; i++) {
			for (j = 0; j < check; j++) {
				SQL[inda[i]] += y_mul[j][i] * lambda[(*pm).ind[ERROR][index_i[j]]];
				if (i == 0) {
					Syl += (*pm).y[(*pm).ind[ERROR][index_i[j]]] * lambda[(*pm).ind[ERROR][index_i[j]]];
				}
			}
		}


		for (i = 0; i < check; i++) {
			(*pm).ind[ERROR][index_i[i]] = -1;
		}
		for (i = 0; i < check; i++) {
			for (j = 0; j < (*pm).num_e - 1; j++) {
				if ((*pm).ind[ERROR][j] == -1) {
					(*pm).ind[ERROR][j] = (*pm).ind[ERROR][j + 1];
					(*pm).ind[ERROR][j + 1] = -1;
				}
			}
		}
		for (i = 0; i < check; i++) {
			free(y_mul[i]);
		}
		free(y_mul);
		(*pm).num_e -= check;
	}

	check = 0;
	for (j = 0; j < (*pm).num_r; j++) {
		flag[j] = ((*pm).g[(*pm).ind[RESERVE][j]] <= 0);
		if (flag[j] == 1) {
			index_i[check++] = j;
		}
	}
	if (check > 0) {
		y_mul = (double**)malloc(sizeof(double*)*inda_size);
		temp = (double**)malloc(sizeof(double*)*inda_size);
		for (i = 0; i < inda_size; i++) {
			y_mul[i] = (double*)malloc(sizeof(double)*check);
			temp[i] = (double*)malloc(sizeof(double)*check);
		}

		for (j = 0; j < check; j++) {
			(*pm).ind[UNLEARNED][(*pm).num_u] = (*pm).ind[RESERVE][index_i[j]];
			(*pm).num_u++;
			lambda[(*pm).ind[RESERVE][index_i[j]]] = (*pm).C[(*pm).ind[RESERVE][index_i[j]]];
		}
		for (i = 0; i < check; i++) {
			kernel_second[i] = (*pm).ind[RESERVE][index_i[i]];
		}
		kernel(inda, kernel_second, inda_size, check, temp,pm);

		for (i = 0; i < inda_size; i++) {
			for (j = 0; j < check; j++) {
				y_mul[i][j] = (*pm).y[inda[i]] * (*pm).y[(*pm).ind[RESERVE][index_i[j]]] * temp[i][j];
			}
		}

		for (i = 0; i < inda_size; i++) {
			for (j = 0; j < check; j++) {
				SQL[inda[i]] += y_mul[i][j] * lambda[(*pm).ind[RESERVE][index_i[j]]];
				if (i == 0) {
					Syl += (*pm).y[(*pm).ind[RESERVE][index_i[j]]] * lambda[(*pm).ind[RESERVE][index_i[j]]];
				}
			}
		}
		for (i = 0; i < check; i++) {
			(*pm).ind[RESERVE][index_i[i]] = -1;
		}
		for (i = 0; i < check; i++) {
			for (j = 0; j < (*pm).num_r - 1; j++) {
				if ((*pm).ind[RESERVE][j] == -1) {
					(*pm).ind[RESERVE][j] = (*pm).ind[RESERVE][j + 1];
					(*pm).ind[RESERVE][j + 1] = -1;
				}
			}
		}
		for (j = 0; j < inda_size; j++) {
			free(y_mul[j]);
			free(temp[j]);
		}
		free(y_mul);
		free(temp);
		(*pm).num_r -= check;
	}
	check = 0;
	for (i = 0; i < (*pm).num_m; i++) {
		flag[i] = ((*pm).g[(*pm).ind[MARGIN][i]] > 0);
		if (flag[i] == 1) {
			index_i[check++] = i;
		}
	}
	if (check > 0) {
		y_mul = (double**)malloc(sizeof(double)*check);
		for (i = 0; i < check; i++) {
			lambda[(*pm).ind[MARGIN][index_i[i]]] = -1 * (*pm).a[(*pm).ind[MARGIN][index_i[i]]];
			y_mul[i] = (double*)malloc(sizeof(double)*inda_size);
		}

		for (i = 0; i < check; i++) {
			for (j = 0; j < inda_size; j++) {
				y_mul[i][j] = (*pm).y[(*pm).ind[MARGIN][index_i[i]]] * (*pm).y[inda[j]] * K[index_i[i]][j];
			}
		}

		for (i = 0; i < check; i++) {
			for (j = 0; j < inda_size; j++) {
				SQL[inda[j]] += y_mul[i][j] * lambda[(*pm).ind[MARGIN][index_i[i]]];
			}
			Syl += (*pm).y[(*pm).ind[MARGIN][index_i[i]]] * lambda[(*pm).ind[MARGIN][index_i[i]]];
		}
		for (i = 0; i < check; i++) {
			free(y_mul[i]);
		}
		free(y_mul);
	}
	check = 0;
	for (i = 0; i < (*pm).num_m; i++) {
		flag[i] += 1;
		if (flag[i] == 1) {
			index_i[check++] = i;
		}
	}
	if (check > 0) {
		y_mul = (double**)malloc(sizeof(double*)*check);
		for (i = 0; i < check; i++) {
			y_mul[i] = (double*)malloc(sizeof(double)*inda_size);
			lambda[(*pm).ind[MARGIN][index_i[i]]] = (*pm).C[(*pm).ind[MARGIN][index_i[i]]] - (*pm).a[(*pm).ind[MARGIN][index_i[i]]];
		}
		for (i = 0; i < check; i++) {
			for (j = 0; j < inda_size; j++) {
				y_mul[i][j] = (*pm).y[(*pm).ind[MARGIN][index_i[i]]] * (*pm).y[inda[j]];
				y_mul[i][j] *= K[index_i[i]][j];
			}
		}
		for (i = 0; i < inda_size; i++) {
			sum = 0;
			for (j = 0; j < check; j++) {
				sum += y_mul[j][i] * lambda[(*pm).ind[MARGIN][index_i[j]]];
				if (i == 0) {
					Syl += (*pm).y[(*pm).ind[MARGIN][index_i[j]]] * lambda[(*pm).ind[MARGIN][index_i[j]]];
				}
			}
			SQL[inda[i]] += sum;
		}
		for (i = 0; i < check; i++) {
			free(y_mul[i]);
		}
		free(y_mul);
	}
	for (i = 0; i < (*pm).num_m; i++) {
		(*pm).ind[UNLEARNED][(*pm).num_u] = (*pm).ind[MARGIN][i];
		(*pm).num_u++;
	}
	(*pm).num_m = 0;
	for (i = 0; i < (*pm).num_u; i++) {
		SQL[(*pm).ind[UNLEARNED][i]] += (*pm).deps * lambda[(*pm).ind[UNLEARNED][i]];
	}
	sprintf(s, "Number of unlearned vectors : %d\n", (*pm).num_u);
	printf("%s", s);

	p_s = 0;
	num_MVs = (*pm).num_m;
	num_learned = 0;
	(*pm).perturbations = 0;
	(*pm).Q = (double**)malloc(sizeof(double*) * 1);
	(*pm).Rs = (double**)malloc(sizeof(double*) * 2);
	(*pm).Q[0] = (double*)malloc(sizeof(double)*(*pm).data_col);
	(*pm).Rs[0] = (double*)malloc(sizeof(double)*2);
	(*pm).Rs[1] = (double*)malloc(sizeof(double)*2);
	for (i = 0; i < (*pm).data_col; i++) {
		(*pm).Q[0][i] = (*pm).y[i];
	}

	while (((*pm).num_u > 0) || ((p_s < 1) && ((*pm).num_u == 0))) {
		double * beta = (double*)malloc(sizeof(double)*((*pm).num_m+1));
		int * ind_temp = (int*)malloc(sizeof(int)*((*pm).num_e + (*pm).num_r + (*pm).num_u));
		double * gamma = (double*)malloc(sizeof(double)*(*pm).data_col);
		memset(beta, 0, sizeof(double)*((*pm).num_m + 1));
		memset(gamma, 0, sizeof(double)*((*pm).data_col));

		(*pm).perturbations += 1;

		if (num_MVs > 0) {
			double * v = (double*)malloc(sizeof(double)*(num_MVs+1));
			double v_temp = 0;

			if (p_s < 1 - 2.2204e-16) {
				v[0] = -1 * Syl;
				for (i = 0; i < (*pm).data_col; i++) {
					v_temp += (*pm).y[i] * (*pm).a[i];
				}
				v[0] = v[0] - v_temp / (1 - p_s);
			}
			else {
				v[0] = -1 * Syl;
			}

			for (i = 0; i < num_MVs + 1; i++) {
				if (i + 1 < num_MVs + 1) {
					v[i + 1] = -1 * SQL[(*pm).ind[MARGIN][i]];
				}
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
				else if ((i >= (*pm).num_e) && (i < (*pm).num_e + (*pm).num_r)) {
					ind_temp[i] = (*pm).ind[RESERVE][i - (*pm).num_e];
				}
				else {
					ind_temp[i] = (*pm).ind[UNLEARNED][i - (*pm).num_e - (*pm).num_r];
				}
			}

			for (i = 0; i < (*pm).num_e + (*pm).num_r + (*pm).num_u; i++) {
				for (j = 0; j < num_MVs + 1; j++) {
					gamma[ind_temp[i]] += (*pm).Q[j][ind_temp[i]] * beta[j];
				}
				gamma[ind_temp[i]] += SQL[ind_temp[i]];
			}
			free(v);
		}
		else {
			beta[0] = 0;
			for (i = 0; i < (*pm).data_col; i++) {
				gamma[i] = SQL[i];
			}
		}
		double min_dps;
		int indss, cstatus, nstatus;
		int indco, removed_i;

		min_delta_p_s(p_s,gamma,beta,lambda,&min_dps,&indss,&cstatus,&nstatus,pm);

		printf("%d   min_dps : %lf %d %d %d\n", number,min_dps,indss,cstatus,nstatus);
		number++;
		if ((*pm).num_u > 0) {
			for (i = 0; i < (*pm).num_u; i++) {
				(*pm).a[(*pm).ind[UNLEARNED][i]] += lambda[(*pm).ind[UNLEARNED][i]] * min_dps;
			}
		}
		if (num_MVs > 0) {
			for (i = 0; i < num_MVs; i++) {
				(*pm).a[(*pm).ind[MARGIN][i]] += beta[i + 1] * min_dps;
			}
		}
		(*pm).b += beta[0] * min_dps;
		for (i = 0; i < (*pm).data_col; i++) {
			(*pm).g[i] += gamma[i] * min_dps;
		}
		p_s += min_dps;
		bookkeeping(indss, cstatus, nstatus, &indco, &removed_i,pm);
		for (i = 0; i < num_MVs; i++) {
			(*pm).g[(*pm).ind[MARGIN][i]] = 0;
		}
		if (nstatus == MARGIN) {
			num_MVs += 1;
			if (num_MVs > 1) {
				int * indc_s1 = (int*)malloc(sizeof(int) * 1);
				temp = (double**)malloc(sizeof(double*)*1);
				temp[0] = (double*)malloc(sizeof(double)*1);
				indc_s1[0] = indss;

				for (i = 0; i < num_MVs; i++) {
					beta[i] = 0;
					for (j = 0; j < num_MVs; j++) {
						beta[i] += -1 * (*pm).Rs[i][j] * (*pm).Q[j][indss];
					}
				}
				kernel(indc_s1, indc_s1, 1, 1, temp,pm);
				gamma[0] = 0;
				for (i = 0; i < num_MVs; i++) {
					gamma[0] += (*pm).Q[i][indss] * beta[i];
				}
				gamma[0] += (*pm).deps + temp[0][0];
				free(indc_s1);
				free(temp[0]);
				free(temp);
			}
			updateRQ(beta, gamma, indss, 1,pm);
		}
		else if (cstatus == MARGIN) {
			num_MVs -= 1;
			updateRQ(beta, gamma, indco, 0,pm);
		}
		if (cstatus == UNLEARNED) {
			num_learned += 1;
			if (nstatus == MARGIN) {
				for (i = 0; i < (*pm).data_col; i++) {
					SQL[i] -= (*pm).Q[num_MVs][i] * lambda[indss];
				}
			}
			else {
				temp = (double**)malloc(sizeof(double*)*(*pm).data_col);
				int* indices = (int*)malloc(sizeof(int)*(*pm).data_col);
				for (i = 0; i < (*pm).data_col; i++) {
					temp[i] = (double*)malloc(sizeof(double)*1);
					indices[i] = i;
				}
				int indc_s[1];
				indc_s[0] = indss;
				kernel(indices, indc_s, (*pm).data_col, 1, temp,pm);
				for (i = 0; i < (*pm).data_col; i++) {
					SQL[i] -= (*pm).y[i] * (*pm).y[indss] * temp[i][0] * lambda[indss];
				}
				SQL[indss] -= (*pm).deps * lambda[indss];
				for (i = 0; i < (*pm).data_col; i++) {
					free(temp[i]);
				}
				free(temp);
				free(indices);
			}
			Syl -= (*pm).y[indss] * lambda[indss];
		}
		free(beta);
		free(gamma);
		free(ind_temp);
	}
	printf("Perturbation complete!\n");
	printf("MARGIN vectors : %d\n", (*pm).num_m);
	printf("ERROR vectors : %d\n", (*pm).num_e);
	printf("RESERVE vectors : %d\n", (*pm).num_r);
	printf("kernel evaluations : %d\n", (*pm).kernel_evals - kernel_evals_begin);

	free(f);
	for (i = 0; i < K_size; i++) {
		free(K[i]);
	}
	free(K);
	free(inda);
	free(SQL);
	free(lambda);
	free(flag);
	free(index_i);
	free(kernel_second);
}