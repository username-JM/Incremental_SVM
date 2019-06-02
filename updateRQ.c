#include"global.h"

void updateRQ(double * beta, double * gamma, int indc, int flag,struct parameter * pm)
{
	int expand, i, j;
	int size_gamma = 1;
	double * temp_beta;
	double ** temp_Q, ** temp_Rs, ** temp;
	int * temp_indc, *temp_ind, * stripped;
	int mat_size = (*pm).num_m;

	if (flag == 1)
		expand = 1;
	else if (flag == 0)
		expand = 0;
	else
		//error
		expand = 0;

	if (expand == 1)
	{
		if (gamma[0] < (*pm).deps && (*pm).num_m > 1)
			gamma[0] = (*pm).deps;

		temp_indc = (int *)malloc(sizeof(int));
		temp_indc[0] = indc;
		if (mat_size > 1)
		{
			(*pm).Rs = (double **)realloc((*pm).Rs, sizeof(double*)*(mat_size + 1));
			for (i = 0; i < mat_size; i++) {
				(*pm).Rs[i] = (double *)realloc((*pm).Rs[i], sizeof(double)*(mat_size + 1));  // 늘린부분 초기화
			}
			(*pm).Rs[mat_size] = (double*)malloc(sizeof(double)*(mat_size + 1));
			for (i = 0; i < mat_size+1; i++) {
				(*pm).Rs[i][mat_size] = 0;
				(*pm).Rs[mat_size][i] = 0;
			}
			
			temp_beta = (double *)malloc(sizeof(double)*(mat_size + 1));
			for (i = 0; i < mat_size; i++)
				temp_beta[i] = beta[i];
			temp_beta[mat_size] = 1;

			for (i = 0; i < mat_size + 1; i++)
			{
				for (j = 0; j < mat_size + 1; j++) {
					(*pm).Rs[i][j] += (temp_beta[i] * temp_beta[j]) / gamma[0];
				}
			}
	
			free(temp_beta);
		}
		else
		{
			kernel(temp_indc, temp_indc, 1, 1, (*pm).Rs,pm);
			(*pm).Rs[0][0] += (*pm).deps;
			(*pm).Rs[0][0] *= -1;
			(*pm).Rs = (double**)realloc((*pm).Rs, sizeof(double *)*(mat_size + 1));
			(*pm).Rs[0] = (double*)realloc((*pm).Rs[0], sizeof(double) * 2);
			(*pm).Rs[0][1] = (*pm).y[indc];
			(*pm).Rs[1] = (double *)malloc(sizeof(double)*(mat_size + 1));
			(*pm).Rs[1][1] = 0;
			(*pm).Rs[1][0] = (*pm).y[indc];
		}

		temp_ind = (int *)malloc(sizeof(int)*(*pm).data_col);
		for (i = 0; i < (*pm).data_col; i++)
			temp_ind[i] = i;
		temp_Q = (double**)malloc(sizeof(double*));
		temp_Q[0] = (double*)malloc(sizeof(double)*(*pm).data_col);
		kernel(temp_indc, temp_ind, 1, (*pm).data_col, temp_Q,pm);
		free(temp_ind);

		for (i = 0; i < (*pm).data_col; i++) {
			if (temp_Q[0][i] > 2.2204e-16 || temp_Q[0][i] < -2.2204e-16) {
				temp_Q[0][i] *= ((*pm).y[indc] * (*pm).y[i]);
			}
			
		}

		(*pm).Q = (double **)realloc((*pm).Q, sizeof(double*)*(mat_size + 1));
		(*pm).Q[mat_size] = temp_Q[0];
		(*pm).Q[mat_size][indc] += (*pm).deps;
		mat_size++;

	}


	else
	{
		if (mat_size+2 > 2)
		{
			stripped = (int *)malloc(sizeof(int) * (mat_size + 1));
			for (i = 0; i < indc; i++)
				stripped[i] = i;
			i++;
			for (; i < mat_size+2; i++)
				stripped[i - 1] = i;
			temp_Rs = (double **)malloc(sizeof(double*)*(mat_size + 1));
			for (i = 0; i < mat_size+1; i++)
				temp_Rs[i] = (double *)malloc(sizeof(double)*(mat_size + 1));

			for (i = 0; i < mat_size + 1; i++)
			{
				for (j = 0; j < mat_size + 1; j++)
				{
					temp_Rs[i][j] = (*pm).Rs[stripped[i]][stripped[j]];
				}
			}
			for (i = 0; i < mat_size +1; i++)
			{
				for (j = 0; j < mat_size + 1; j++) {
					temp_Rs[i][j] -= ((*pm).Rs[stripped[i]][indc] * (*pm).Rs[stripped[j]][indc]) / (*pm).Rs[indc][indc];
				}
			}
			for (i = 0; i < mat_size+2; i++)
				free((*pm).Rs[i]);
			free((*pm).Rs);
			(*pm).Rs = temp_Rs;
		}

		else
			(*pm).Rs[0][0] = INF;

		temp = (double **)malloc(sizeof(double*)*(mat_size + 1));
		int check = 0;
		for (i = 0; i < mat_size + 2; i++)
		{
			if (i != indc) {
				temp[check++] = (*pm).Q[i];
			}
		}

		free((*pm).Q[indc]);
		free((*pm).Q);
		(*pm).Q = temp;
	}
}
