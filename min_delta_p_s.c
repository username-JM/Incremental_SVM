#include "global.h"

void min_delta_p_s(double p_s, double* gamma, double* beta, double* lambda,
	double* min_dps, int* indss, int* cstatus, int* nstatus,struct parameter* pm) {

	int indss_s[7];
	int cstatus_s[7];
	int nstatus_s[7];

	memset(indss_s, 0, sizeof(int) * 7);
	memset(cstatus_s, 0, sizeof(int) * 7);
	memset(nstatus_s, 0, sizeof(int) * 7);
	
	double delta_ue, delta_um, delta_ur, delta_mer, delta_em, delta_rm;
	double delta_p_s;
	int i, j;
	if ((*pm).num_u > 0) {
		double * lambda_u = (double*)malloc(sizeof(double)*(*pm).num_u);
		double * gamma_u = (double*)malloc(sizeof(double)*(*pm).num_u);
		double * g_u = (double*)malloc(sizeof(double)*(*pm).num_u);
		int * flags = (int*)malloc(sizeof(int)*(*pm).num_u);
		double * a_s1 = (double*)malloc(sizeof(double)*(*pm).num_u);
		double * a_s2 = (double*)malloc(sizeof(double)*(*pm).num_u);

		delta_p_s = INF;

		for (j = 0; j < (*pm).num_u; j++) {
			lambda_u[j] = lambda[(*pm).ind[UNLEARNED][j]];
			flags[j] = (lambda_u[j] > 0);

			a_s1[j] = (*pm).a[(*pm).ind[UNLEARNED][j]];
			a_s2[j] = (*pm).C[(*pm).ind[UNLEARNED][j]];
		}
		min_delta(flags, a_s1, a_s2, lambda_u, (*pm).num_u, &delta_ue, &i);
		if (delta_ue < INF) {
			indss_s[0] = (*pm).ind[UNLEARNED][i];
			cstatus_s[0] = UNLEARNED;
			nstatus_s[0] = ERROR;
		}
		//  UNLEARNED =======================> ERROR
		/////////////////////////////////////////////////////////////////////////////////
		for (j = 0; j < (*pm).num_u; j++) {
			gamma_u[j] = gamma[(*pm).ind[UNLEARNED][j]];
			g_u[j] = (*pm).g[(*pm).ind[UNLEARNED][j]];
			flags[j] = (gamma_u[j] * g_u[j] < 0);
			a_s2[j] = 0;
		}
		min_delta(flags, g_u, a_s2, gamma_u, (*pm).num_u, &delta_um, &i);
		if (delta_um < INF) {
			indss_s[1] = (*pm).ind[UNLEARNED][i];
			cstatus_s[1] = UNLEARNED;
			nstatus_s[1] = MARGIN;
		}
		//  UNLEARNED =======================> MARGIN
		//////////////////////////////////////////////////////////////////////////////////

		for (j = 0; j < (*pm).num_u; j++) {
			flags[j] = (lambda_u[j] < 0);
		}
		min_delta(flags, a_s1, a_s2, lambda_u, (*pm).num_u, &delta_ur, &i);
		if (delta_ur < INF) {
			indss_s[2] = (*pm).ind[UNLEARNED][i];
			cstatus_s[2] = UNLEARNED;
			nstatus_s[2] = RESERVE;
		}
		//  UNLEARNED =======================> RESERVE
		///////////////////////////////////////////////////////////////////////////////////
		free(g_u);
		free(gamma_u);
		free(lambda_u);
		free(flags);
		free(a_s1);
		free(a_s2);
	}
	else {
		delta_p_s = 1 - p_s;
		delta_ue = INF;
		delta_um = INF;
		delta_ur = INF;
	}
	if ((*pm).num_m > 0) {
		double * beta_s = (double*)malloc(sizeof(double)*(*pm).num_m);
		int * flags_m = (int*)malloc(sizeof(int)*(*pm).num_m);
		double * a_m1 = (double*)malloc(sizeof(double)*(*pm).num_m);
		double * a_m2 = (double*)malloc(sizeof(double)*(*pm).num_m);

		for (j = 0; j < (*pm).num_m; j++) {
			beta_s[j] = beta[j + 1];
			flags_m[j] = ((beta_s[j]) != 0);
			a_m1[j] = (*pm).a[(*pm).ind[MARGIN][j]];
			a_m2[j] = (*pm).C[(*pm).ind[MARGIN][j]] * (beta_s[j] > 0);
		}
		min_delta(flags_m, a_m1, a_m2, beta_s, (*pm).num_m, &delta_mer, &i);
		if (delta_mer < INF) {
			indss_s[4] = (*pm).ind[MARGIN][i];
			cstatus_s[4] = MARGIN;
			nstatus_s[4] = ERROR * (beta_s[i] > 0) + RESERVE * (beta_s[i] < 0);
		}
		free(beta_s);
		free(flags_m);
		free(a_m1);
		free(a_m2);
	}
	else {
		delta_mer = INF;
	}
	//  MARGIN ===========================================> ERROR,RESERVE
	/////////////////////////////////////////////////////////////////////////////////

	double * gamma_e = (double*)malloc(sizeof(double)*(*pm).num_e);
	int * flags_e = (int*)malloc(sizeof(int)*(*pm).num_e);
	double * a_e1 = (double*)malloc(sizeof(double)*(*pm).num_e);
	double * a_e2 = (double*)malloc(sizeof(double)*(*pm).num_e);

	double * gamma_r = (double*)malloc(sizeof(double)*(*pm).num_r);
	int * flags_r = (int*)malloc(sizeof(int)*(*pm).num_r);
	double * a_r1 = (double*)malloc(sizeof(double)*(*pm).num_r);
	double * a_r2 = (double*)malloc(sizeof(double)*(*pm).num_r);

	for (j = 0; j < (*pm).num_e; j++) {
		gamma_e[j] = gamma[(*pm).ind[ERROR][j]];
		flags_e[j] = (gamma_e[j] > 0);
		a_e1[j] = (*pm).g[(*pm).ind[ERROR][j]];
		a_e2[j] = 0;
	}
	min_delta(flags_e, a_e1, a_e2, gamma_e, (*pm).num_e, &delta_em, &i);
	if (delta_em < INF) {
		indss_s[5] = (*pm).ind[ERROR][i];
		cstatus_s[5] = ERROR;
		nstatus_s[5] = MARGIN;
	}
	//  ERROR ===============================================> MARGIN
	/////////////////////////////////////////////////////////////////////////////////

	for (j = 0; j < (*pm).num_r; j++) {
		gamma_r[j] = gamma[(*pm).ind[RESERVE][j]];
		flags_r[j] = (((*pm).g[(*pm).ind[RESERVE][j]] >= 0) && (gamma_r[j] < 0));
		a_r1[j] = (*pm).g[(*pm).ind[RESERVE][j]];
		a_r2[j] = 0;
	}
	min_delta(flags_r, a_r1, a_r2, gamma_r, (*pm).num_r, &delta_rm, &i);
	if (delta_rm < INF) {
		indss_s[6] = (*pm).ind[RESERVE][i];
		cstatus_s[6] = RESERVE;
		nstatus_s[6] = MARGIN;
	}
	//  RESERVE ==============================================> MARGIN
	///////////////////////////////////////////////////////////////////////////////////

	double arr[7];
	arr[0] = delta_ue;
	arr[1] = delta_um;
	arr[2] = delta_ur;
	arr[3] = delta_p_s;
	arr[4] = delta_mer;
	arr[5] = delta_em;
	arr[6] = delta_rm;

	*min_dps = arr[0];
	i = 0;
	for (j = 1; j < 7; j++) {
		if (*min_dps > arr[j]) {
			*min_dps = arr[j];
			i = j;
		}
	}

	*indss = indss_s[i];
	*cstatus = cstatus_s[i];
	*nstatus = nstatus_s[i];

	free(gamma_e);
	free(flags_e);
	free(a_e1);
	free(a_e2);
	free(gamma_r);
	free(flags_r);
	free(a_r1);
	free(a_r2);
}