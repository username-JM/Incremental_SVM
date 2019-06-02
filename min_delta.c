#include "global.h"

void min_delta(int * flags, double * psi_initial, double * psi_final, double * psi_sens,
	int num ,double * min_d, int * k) {
	int i, j;
	int check = 0;
	int check2 = 0;
	int * index = (int*)malloc(sizeof(int)*num);
	int *min_k;
	double max_sens;
	double * deltas = (double *)malloc(sizeof(double)*num);

	for (j = 0; j < num; j++) {
		if (flags[j] == 1) {
			index[check] = j;
			check++;
		}
	} //flags�� 1�� index�� �����ϰ� �������ֱ�

	if (check >= 1) { // 1�� index�� ������ ���
		for (j = 0; j < check; j++) {
			deltas[j] = (psi_final[index[j]] - psi_initial[index[j]]) / psi_sens[index[j]];
			if (j == 0) {
				*min_d = deltas[j];
				i = j;
			}
			else if (*min_d > deltas[j]) {
				*min_d = deltas[j];
				i = j;
			}
		}
		*k = index[i];
		min_k = (int *)malloc(sizeof(int)*check);
		for (j = 0; j < check; j++) {
			if (deltas[j] == *min_d) {
				min_k[check2++] = j;
			} // �ּ� deltas�� ���� index�� ������ ���ֱ�
		}
		if (check2 > 1) { // �������� deltas�� �������ΰ��
			for (j = 0; j < check2; j++) {
				if (j == 0) {
					max_sens = psi_sens[index[min_k[j]]];
					if (max_sens < 0) {
						max_sens = max_sens * -1;
					}
					i = j;
					*k = index[min_k[i]];
				}
				if (psi_sens[index[min_k[j]]] < 0) {
					psi_sens[index[min_k[j]]] *= -1;
				}

				else if (max_sens < (psi_sens[index[min_k[j]]])) {
					max_sens = psi_sens[index[min_k[j]]];
					i = j;
					*k = index[min_k[i]];
				}
			}
		}

		free(min_k);
	}
	else {
		*min_d = INF;
		*k = -1;
	}
	free(deltas);
	free(index);
}