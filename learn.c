#include "global.h"

void learn(int indc, int rflag,struct parameter * pm) {
	int i, j;
	double f_c;
	double ** k = (double **)malloc(sizeof(double*)*((*pm).num_m+ (*pm).num_e));
	for (i = 0; i < (*pm).num_m + (*pm).num_e; i++) {
		k[i] = (double*)malloc(sizeof(double) * 1);
	}

	int *indc_temp = (int*)malloc(sizeof(int)*1);
	indc_temp[0] = indc;

	svmeval(indc_temp, &f_c, k,pm);

	free(indc_temp);

	(*pm).g[indc] = (*pm).y[indc] * f_c - 1;
	// g값 계산 완료

	int nstatus, cstatus;
	if ((*pm).g[indc] >= 0) {
		int a, b;
		bookkeeping(indc, UNLEARNED, RESERVE, &a, &b,pm);
		nstatus = RESERVE;
		return;
	} // bookkeeping 

	int num_MVs = (*pm).num_m;
	double * Qc_M = (double*)malloc(sizeof(double)*MAX_NUM);
	double * Qc_E = (double*)malloc(sizeof(double)*MAX_NUM);
	double * Qc_R = (double*)malloc(sizeof(double)*MAX_NUM);
	double Qcc;

	if (num_MVs == 0) {
		if ((*pm).num_e > 0) {
			double ** temp;
			int *indc_s = (int*)malloc(sizeof(int) * 1);
			indc_s[0] = indc;
			temp = (double**)malloc(sizeof(double*)*(*pm).num_e);
			for (i = 0; i < (*pm).num_e; i++) {
				temp[i] = (double *)malloc(sizeof(double) * 1);
			}
			kernel((*pm).ind[ERROR], indc_s, (*pm).num_e, 1, temp,pm);

			for (i = 0; i < (*pm).num_e; i++) {
				Qc_E[i] = (*pm).y[(*pm).ind[ERROR][i]] * (*pm).y[indc] + temp[i][0];
			}

			free(indc_s);
			for (i = 0; i < (*pm).num_e; i++) {
				free(temp[i]);
			}
			free(temp);
		}
	}
	else {
		for (i = 0; i < num_MVs; i++) {
			Qc_M[i] = (*pm).y[(*pm).ind[MARGIN][i]] * (*pm).y[indc] * k[i][0];
		}
		if ((*pm).num_e > 0) {
			for (i = 0; i < (*pm).num_e; i++) {
				Qc_E[i] = (*pm).y[(*pm).ind[ERROR][i]] * (*pm).y[indc] * k[num_MVs + i][0];
			}
		}
	}
	for (i = 0; i < (*pm).num_m + (*pm).num_e; i++) {
		free(k[i]);
	}
	free(k);

	if ((*pm).num_r > 0) {
		double ** temp;
		int *indc_s = (int*)malloc(sizeof(int)*1);
		indc_s[0] = indc;
		temp = (double**)malloc(sizeof(double*)*(*pm).num_r);
		for (i = 0; i < (*pm).num_r; i++) {
			temp[i] = (double *)malloc(sizeof(double) * 1);
		}
		kernel((*pm).ind[RESERVE], indc_s, (*pm).num_r,1 ,temp,pm);

		for (i = 0; i < (*pm).num_r; i++) {
			Qc_R[i] = (*pm).y[(*pm).ind[RESERVE][i]] * (*pm).y[indc] * temp[i][0];
		}
			
		free(indc_s);
		for (i = 0; i < (*pm).num_r; i++) {
			free(temp[i]);
		}
		free(temp);
	}
		
	int *indc_s1 = (int*)malloc(sizeof(int)*1);
	double** Qcc_kernel = (double**)malloc(sizeof(double*)*1);
	Qcc_kernel[0] = (double*)malloc(sizeof(double) * 1);
	indc_s1[0] = indc;

	kernel(indc_s1, indc_s1,1, 1,Qcc_kernel,pm);
	Qcc = Qcc_kernel[0][0] + (*pm).deps;

	free(indc_s1);
	free(Qcc_kernel[0]);
	free(Qcc_kernel);
	 // Qc 및 Qcc를 구하는 코드

	int converged = 0;
	while (!converged) {
		(*pm).perturbations += 1;
		double * beta = (double *)malloc(sizeof(double)*((*pm).num_m + 1));
		double * gamma = (double *)malloc(sizeof(double)*((*pm).data_col));
		int * ind_temp = (int *)malloc(sizeof(int) * ((*pm).num_e + (*pm).num_r + 1));
		double * q_beta = (double *)malloc(sizeof(double)*((*pm).num_e+ (*pm).num_r+ 1));

		for (i = 0; i < (*pm).data_col; i++) {
			gamma[i] = 0;
		}

		if (num_MVs > 0) {
			for (i = 0; i < (*pm).num_m + 1; i++) {
				for (j = 0; j < (*pm).num_m + 1; j++) {
					if (j == 0) {
						beta[i] = 0;
						beta[i] += (*pm).Rs[i][j] * (*pm).y[indc];
					}
					else {
						beta[i] += (*pm).Rs[i][j] * Qc_M[j-1];
					}
				}
				beta[i] *= -1;
			}
			for (i = 0; i < (*pm).num_e; i++) {
				ind_temp[i] = (*pm).ind[ERROR][i];
			}
			for (i = 0; i < (*pm).num_r; i++) {
				ind_temp[i + (*pm).num_e] = (*pm).ind[RESERVE][i];
			}
			ind_temp[(*pm).num_r + (*pm).num_e] = indc;
			 // ind_temp 및 beta 값 설정하기

			for (i = 0; i < (*pm).num_e + (*pm).num_r + 1; i++) {
				q_beta[i] = 0;
				for (j = 0; j < (*pm).num_m + 1; j++) {
					q_beta[i] += (*pm).Q[j][ind_temp[i]] * beta[j];
				}

				if (i < (*pm).num_e) {
					gamma[ind_temp[i]] = Qc_E[i] + q_beta[i];
				}
				else if (i >= (*pm).num_e && i < (*pm).num_e + (*pm).num_r) {
					gamma[ind_temp[i]] = Qc_R[i - (*pm).num_e] + q_beta[i];
				}
				else {
					gamma[ind_temp[i]] = Qcc + q_beta[i];
				}
			}

		}// 마진이 하나라도 있는 경우
		else {
			beta[0] = (*pm).y[indc];
			for (i = 0; i < (*pm).data_col; i++) {
				gamma[i] =((*pm).y[indc] * (*pm).y[i]);
			}
		} // svm이 없는 경우

		double min_delta_param;
		int indss;
		min_delta_acb(indc, gamma, beta, 1, rflag, &min_delta_param, &indss, &cstatus, &nstatus,pm);

		printf("delta p, indss, cstatus, nstatus\n");
		printf("%f, %d, %d, %d\n", min_delta_param, indss, cstatus, nstatus);
		
		if (num_MVs > 0) {
			(*pm).a[indc] = (*pm).a[indc] + min_delta_param;
			for (i = 0; i < (*pm).num_m; i++) {
				(*pm).a[(*pm).ind[MARGIN][i]] += beta[i + 1] * min_delta_param;
			}
		}// alpha값 업데이트 
		(*pm).b += beta[0] * min_delta_param;
		for (i = 0; i < (*pm).data_col; i++) {
			(*pm).g[i] += gamma[i] * min_delta_param;
		}//b,g값 업데이트

		if (indss == indc) {
			converged = 1;
		}
		if (converged) {
			cstatus = UNLEARNED;
			if (nstatus == MARGIN) {
				Qc_M[(*pm).num_m] = Qcc;
			}
			else if (nstatus == ERROR) {
				Qc_E[(*pm).num_e] = Qcc;
			}
		} // unlearned를 svm에 추가시키는 경우
		else {
			int indss_temp;
			if (cstatus == MARGIN) {
				for (i = 0; i < (*pm).num_m; i++) {
					if ((*pm).ind[cstatus][i] == indss) {
						indss_temp = i;
					}
				}
				if (nstatus == RESERVE) {
					Qc_R[(*pm).num_r] = Qc_M[indss_temp];
					if (indss_temp == (*pm).num_m - 1) {
						Qc_M[indss_temp] = 0;
					}else {
						for (i = indss_temp; i < (*pm).num_m - 1; i++) {
							Qc_M[i] = Qc_M[i + 1];
							Qc_M[i + 1] = 0;
						}
					}
				}
				else if (nstatus == ERROR) {
					Qc_E[(*pm).num_e] = Qc_M[indss_temp];
					if (indss_temp == (*pm).num_m - 1) {
						Qc_M[indss_temp] = 0;
					}
					else {
						for (i = indss_temp; i < (*pm).num_m - 1; i++) {
							Qc_M[i] = Qc_M[i + 1];
							Qc_M[i + 1] = 0;
						}
					}
				}
			}
			else if (cstatus == ERROR) {
				for (i = 0; i < (*pm).num_e; i++) {
					if ((*pm).ind[cstatus][i] == indss) {
						indss_temp = i;
					}
				}
				Qc_M[(*pm).num_m] = Qc_E[indss_temp];

				if (indss_temp == (*pm).num_e - 1) {
					Qc_E[indss_temp] = 0;
				}
				else {
					for (i = indss_temp; i < (*pm).num_e-1; i++) {
						Qc_E[i] = Qc_E[i + 1];
						Qc_E[i + 1] = 0;
					}
				}
			}
			else if (cstatus == RESERVE) {
				for (i = 0; i < (*pm).num_r; i++) {
					if ((*pm).ind[cstatus][i] == indss) {
						indss_temp = i;
					}
				}
				Qc_M[(*pm).num_m] = Qc_R[indss_temp];

				if (indss_temp == (*pm).num_r - 1) {
					Qc_R[indss_temp] = 0;
				}
				else {
					for (i = indss_temp; i < (*pm).num_r - 1; i++) {
						Qc_R[i] = Qc_R[i + 1];
						Qc_R[i + 1] = 0;
					}
				}
			}
		} //기존에있던 svm에서 index를 이동하는 경우
		
		int indco, removed_i;
		bookkeeping(indss, cstatus, nstatus, &indco, &removed_i,pm);
		if ((nstatus == RESERVE) && (removed_i > 0)) {
			if (removed_i == (*pm).num_r - 1) {
				Qc_R[removed_i] = 0;
			} 
			else {
				for (i = 0; i < (*pm).num_r - 1; i++) {
					Qc_R[i] = Qc_R[i + 1];
					Qc_R[i + 1] = 0;
				}
			}
		} // reserve영역에서 삭제하는 경우
		for (i = 0; i < (*pm).num_m; i++) {
			(*pm).g[(*pm).ind[MARGIN][i]] = 0;
		} // margin g값 0

		if (nstatus == MARGIN) {
			num_MVs += 1;
			if (num_MVs > 1) {
				if (converged) {
					gamma[0] = gamma[indss];
				}
				else {
					for (i = 0; i < num_MVs; i++) {
						beta[i] = 0;
						for (j = 0; j < num_MVs; j++) {
							beta[i] += (*pm).Rs[i][j] * (*pm).Q[j][indss];
						}
						beta[i] = beta[i] * -1;
					}
					int *indc_s1 = (int*)malloc(sizeof(int)*1);
					double **temp = (double**)malloc(sizeof(double*)*1);
					temp[0] = (double*)malloc(sizeof(double) * 1);
					indc_s1[0] = indss;

					gamma[0] = 0;
					kernel(indc_s1, indc_s1,1,1, temp,pm);
					for (i = 0; i < (*pm).num_m; i++) {
						gamma[0] += (*pm).Q[i][indss] * beta[i];
					}
					gamma[0] = gamma[0] + (*pm).deps + temp[0][0];

					free(indc_s1);
					free(temp[0]);
					free(temp);
				}
			}
			updateRQ(beta, gamma, indss,1,pm);
		}
		else if (cstatus == MARGIN) {
			num_MVs -= 1;
			updateRQ(beta,gamma,indco,0,pm);
		}

		free(beta);
		free(gamma);
		free(q_beta);
		free(ind_temp);
	}
	free(Qc_M);
	free(Qc_E);
	free(Qc_R);
}