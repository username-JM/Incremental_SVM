#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#pragma warning(disable:4996)
#define MAX_NUM 1000
#define MARGIN 1
#define ERROR 2
#define RESERVE 3
#define UNLEARNED 4
#define INF 9999

struct parameter {
	double initial_C;
	double b; double * C;
	double deps; double scale;
	int kernel_evals;
	int max_reserve_vectors;
	int type;
	int perturbations;
	int *ind[5];
	int num_m; int num_e; int num_r; int num_u;  int num_uind;
	int data_row; int data_col;
	double * a;
	double * g;
	double ** Q;
	double ** Rs;
	double ** X;
	int * uind;
	int * y;
};

void loadclass(struct parameter * pm, int num);
void saveclass(struct parameter pm);
void svmtrain(double *X_new[], int y_new[], int num, int new_model);
void learn(int indc, int rflag, struct parameter * pm);
void kernel(int ind_X[], int ind_y[], int len_ind_x, int len_ind_y, double ** K, struct parameter* pm);
void bookkeeping(int indss, int cstatus, int nstatus, int *indco, int *i, struct parameter * pm);
void min_delta_acb(int indc, double* gamma, double* beta, int polc, int rflag, double* min_dacb,
	int*indss_s, int* cstatus_s, int* nstatus_s, struct parameter * pm);
void svmeval(int * indc, double * f_c, double ** K, struct parameter * pm);
void updateRQ(double * beta, double * gamma, int indc, int flag, struct parameter * pm);
void min_delta(int * flags, double * psi_initial, double * psi_final, double * psi_sens,
	int num, double * min_d, int * k);
void svmeval2(int * indc, double * f_c, double ** K, struct parameter * pm);
void min_delta_p_s(double p_s, double* gamma, double* beta, double* lambda,
	double* min_dps, int* indss, int* cstatus, int* nstatus,struct parameter * pm);
void perturbk(double new_scale,struct parameter * pm);
void perturbc(double C_new, struct parameter * pm);
void min_delta_p_c(double p_c, double * gamma, double * beta, double * lambda,
	double * min_dpc, int * indss, int * cstatus, int * nstatus, struct parameter * pm);