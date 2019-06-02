#include "global.h"

void kernel(int ind_X[], int ind_y[], int len_ind_x, int len_ind_y, double ** K,struct parameter * pm)
{
	int num_feature = (*pm).data_row;
	int i,j,k;
	double sum = 0;
	for(i=0; i<len_ind_x; i++)
	{
		for(j=0; j<len_ind_y; j++)
		{
			for(k=0; k<num_feature; k++)
			{
				sum += (*pm).X[k][ind_X[i]] * (*pm).X[k][ind_y[j]];
			}
			K[i][j] = sum;
			sum = 0;
		}
	}
	
	double * tmp_x = (double *)malloc(sizeof(double)*len_ind_x);
	double * tmp_y = (double *)malloc(sizeof(double)*len_ind_y);
	for(i=0; i<len_ind_x; i++) tmp_x[i] = 0;
	for(i=0; i<len_ind_y; i++) tmp_y[i] = 0;
	
	if(((*pm).type > 1) && ((*pm).type < 5))
	{
		//not yet
	}
	else if((*pm).type == 5)
	{
		for(i=0; i<len_ind_x; i++)
		{
			for(j=0; j<len_ind_y; j++)
			{
				K[i][j] *= 2;
			}
		}
		
		for(i=0; i<len_ind_x; i++)
		{
			for(j=0; j<num_feature; j++)
			{
				tmp_x[i] += (*pm).X[j][ind_X[i]] * (*pm).X[j][ind_X[i]];
			}
		}
		
		for(i=0; i<len_ind_x; i++)
		{
			for(j=0; j<len_ind_y; j++)
			{
				K[i][j] -= tmp_x[i];
			}
		}
		
		for(i=0; i<len_ind_y; i++)
		{
			for(j=0; j<num_feature; j++)
			{
				tmp_y[i] += (*pm).X[j][ind_y[i]] * (*pm).X[j][ind_y[i]];
			}
		}
		
		for(i=0; i<len_ind_x; i++)
		{
			for(j=0; j<len_ind_y; j++)
			{
				K[i][j] -= tmp_y[j];
			}
		}
		
		for(i=0; i<len_ind_x; i++)
		{
			for(j=0; j<len_ind_y; j++)
			{
				K[i][j] = exp(K[i][j] / (2* (*pm).scale));
			}
		}
	}
	
	else
	{
		//type error
	}
	(*pm).kernel_evals += len_ind_x * len_ind_y;
}
