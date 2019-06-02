#include <stdio.h>

// RESERVE 영역에서 데이터를 지울 때에 사용하는 정렬 함수

void quicksort(double g_sorted[], int rem_idx[], int left, int right) {
	int L = left, R = right;
	double g_temp;
	int rem_temp;
	double pivot = g_sorted[(left + right) / 2];

	while (L <= R) {
		while (g_sorted[L] < pivot) {
			L++;
		}
		while (g_sorted[R] > pivot) {
			R--;
		}
	}

	if (L <= R) {
		if (L != R) {
			g_temp = g_sorted[L];
			rem_temp = rem_idx[L];
			g_sorted[L] = g_sorted[R];
			rem_idx[L] = rem_idx[R];
			g_sorted[R] = g_temp;
			rem_idx[R] = rem_temp;
		}
		L++;
		R--;
	}
	
	if (left < R) {
		quicksort(g_sorted, rem_idx, left, R);
	}
	if (L < right) {
		quicksort(g_sorted, rem_idx, L, right);
	}
}