#include "global.h"

void min_delta_p_c(double p_c, double * gamma, double * beta, double * lambda,
	double * min_dpc,int * indss, int * cstatus, int * nstatus, struct parameter * pm) {
	double delta_mr, delta_me, delta_rm, delta_em;
	int indss_s[5];
	int cstatus_s[5];
	int nstatus_s[5];

	memset(indss_s, 0, 20);
	memset(cstatus_s, 0, 20);
	memset(nstatus_s, 0, 20);

	double * beta_s = (double *)malloc(sizeof(double)*(*pm).num_m);
	double * delta_me_s = (double*)malloc(sizeof(double)*(*pm).num_m);
	int * flags = (int*)malloc(sizeof(int)*(*pm).num_m);
	int i, j;

	double delta_p_c = 1 - p_c;

	if ((*pm).num_m > 0) {
		double * a_s1 = (double*)malloc(sizeof(double)*(*pm).num_m);
		double * a_s2 = (double*)malloc(sizeof(double)*(*pm).num_m);

		for (j = 0; j < (*pm).num_m; j++) {
			beta_s[j] = beta[j + 1];
			flags[j] = (beta_s[j] < 0);
		}

		for (j = 0; j < (*pm).num_m; j++) {
			a_s1[j] = (*pm).a[(*pm).ind[MARGIN][j]];
			a_s2[j] = 0;
		}
		min_delta(flags, a_s1, a_s2, beta_s, (*pm).num_m, &delta_mr, &i);
		if (delta_mr < INF) {
			indss_s[1] = (*pm).ind[MARGIN][i];
			cstatus_s[1] = MARGIN;
			nstatus_s[1] = RESERVE;
		}
		free(a_s1);
		free(a_s2);
	} // MARGIN ¿¡¼­ RESERVE
	else { // no margin
		delta_mr = INF;
	}
	
	if ((*pm).num_m > 0) {
		double * lambda_s = (double*)malloc(sizeof(double)*(*pm).num_m);
		double * v = (double*)malloc(sizeof(double)*(*pm).num_m);
		int * not_z = (int*)malloc(sizeof(int)*(*pm).num_m);
		int check = 0;

		for (j = 0; j < (*pm).num_m; j++) {
			lambda_s[j] = lambda[(*pm).ind[MARGIN][j]];
			v[j] = beta_s[j] - lambda_s[j];
		}

		for (j = 0; j < (*pm).num_m; j++) {
			if (v[j] > 2.2204e-16) {
				flags[j] = 1;
				check = 1;
			}
			else {
				flags[j] = 0;
			}
		}
		if (check == 1) {
			check = 0;
			for (j = 0; j < (*pm).num_m; j++) {
				if (v[j] > 0) {
					not_z[check++] = j;
				}
			}
			for (j = 0; j < check; j++) {
				delta_me_s[j] = ((*pm).C[(*pm).ind[MARGIN][not_z[j]]] - (*pm).a[(*pm).ind[MARGIN][not_z[j]]]) / v[not_z[j]];
				if (j == 0) {
					delta_me = delta_me_s[j];
					i = not_z[j];
				}
				else if(delta_me > delta_me_s[j]) {
					delta_me = delta_me_s[j];
					i = not_z[j];
				}
			}
			if (delta_me < INF) {
				indss_s[2] = (*pm).ind[MARGIN][i];
				cstatus_s[2] = MARGIN;
				nstatus_s[2] = ERROR;
			}
		}
		else {
			delta_me = INF;
		}
		free(not_z);
		free(lambda_s); //error
		free(v);
	}
	else {
		delta_me = INF;
	}

	double * gamma_e = (double*)malloc(sizeof(double)*(*pm).num_e);
	int * flags_e = (int*)malloc(sizeof(int)*(*pm).num_e);
	double * g_e1 = (double*)malloc(sizeof(double)*(*pm).num_e);
	double * g_e2 = (double*)malloc(sizeof(double)*(*pm).num_e);

	for (j = 0; j < (*pm).num_e; j++) {
		gamma_e[j] = gamma[(*pm).ind[ERROR][j]];
		if (gamma_e[j] > 0) {
			flags_e[j] = 1;
		}
		else {
			flags_e[j] = 0;
		}
		g_e1[j] = (*pm).g[(*pm).ind[ERROR][j]];
		g_e2[j] = 0;
	}

	min_delta(flags_e, g_e1, g_e2, gamma_e, (*pm).num_e, &delta_em, &i);
	if (delta_em < INF) {
		indss_s[3] = (*pm).ind[ERROR][i];
		cstatus_s[3] = ERROR;
		nstatus_s[3] = MARGIN;
	}

	double * gamma_r = (double *)malloc(sizeof(double)*(*pm).num_r);
	int * flags_r = (int*)malloc(sizeof(int)*(*pm).num_r);
	double * g_r1 = (double*)malloc(sizeof(double)*(*pm).num_r);
	double * g_r2 = (double*)malloc(sizeof(double)*(*pm).num_r);

	for (j = 0; j < (*pm).num_r; j++) {
		gamma_r[j] = gamma[(*pm).ind[RESERVE][j]];
		if ((*pm).g[(*pm).ind[RESERVE][j]] >= 0 && gamma_r[j] < 0) {
			flags_r[j] = 1;
		}
		else {
			flags_r[j] = 0;
		}
		g_r1[j] = (*pm).g[(*pm).ind[RESERVE][j]];
		g_r2[j] = 0;
	}
	min_delta(flags_r, g_r1, g_r2, gamma_r, (*pm).num_r, &delta_rm, &i);
	if (delta_rm < INF) {
		indss_s[4] = (*pm).ind[RESERVE][i];
		cstatus_s[4] = RESERVE;
		nstatus_s[4] = MARGIN;
	}

	double arr[5];
	arr[0] = delta_p_c;
	arr[1] = delta_mr;
	arr[2] = delta_me;
	arr[3] = delta_em;
	arr[4] = delta_rm;

	for (j = 0; j < 5; j++) {
		if (j == 0) {
			*min_dpc = arr[0];
			i = j;
		}
		else if (*min_dpc > arr[j]) {
			*min_dpc = arr[j];
			i = j;
		}
	}
	*indss = indss_s[i];
	*cstatus = cstatus_s[i];
	*nstatus = nstatus_s[i];

	free(flags);
	free(beta_s);
	free(delta_me_s);
	free(gamma_e);
	free(flags_e);
	free(g_e1);
	free(g_e2);
	free(gamma_r);
	free(flags_r);
	free(g_r1);
	free(g_r2);
}