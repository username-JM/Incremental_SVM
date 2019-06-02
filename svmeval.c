#include"global.h"

void svmeval(int * indc, double * f_c, double ** K,struct parameter * pm)
{
	int * indu = (int *)malloc(sizeof(int)*(*pm).num_u);
	int * indr = (int *)malloc(sizeof(int)*(*pm).num_r);
	int * indme = (int *)malloc(sizeof(int)*((*pm).num_e + (*pm).num_m));

	int num_indu = 0;
	int num_indr = 0;
	int i, j;
	double ** temp_K;

	for (i = 0; i < (*pm).data_col; i++)
	{
		for (j = 0; j < (*pm).num_u; j++)
		{
			if ((*pm).a[(*pm).ind[UNLEARNED][j]] > 0)
			{
				indu[num_indu] = (*pm).ind[UNLEARNED][j];
				num_indu++;
			}
		}
	}
	for (i = 0; i < (*pm).data_col; i++)
	{
		for (j = 0; j < (*pm).num_r; j++)
		{
			if ((*pm).a[(*pm).ind[RESERVE][j]] > 0)
			{
				indu[num_indr] = (*pm).ind[RESERVE][j];
				num_indr++;
			}
		}
	}

	for (i = 0; i < (*pm).num_m; i++)
	{
		indme[i] = (*pm).ind[MARGIN][i];
	}
	for (; i < (*pm).num_m + (*pm).num_e; i++)
	{
		indme[i] = (*pm).ind[ERROR][i - (*pm).num_m];
	}
	
	*f_c = (*pm).b;
	if ((*pm).num_m + (*pm).num_e > 0)
	{
		kernel(indme, indc, (*pm).num_e + (*pm).num_m, 1, K,pm);
		for (i = 0; i < (*pm).num_m + (*pm).num_e; i++)
			*f_c += K[i][0] * (*pm).y[indme[i]] * (*pm).a[indme[i]];
	}

	if (num_indu > 0)
	{
		temp_K = (double **)malloc(sizeof(double*)*num_indu);
		for (i = 0; i < num_indu; i++)
			temp_K[i] = (double*)malloc(sizeof(double));

		kernel(indu, indc, num_indu, 1, temp_K,pm);
		for (i = 0; i < num_indu; i++)
			*f_c += temp_K[i][0] * (*pm).y[indu[i]] * (*pm).a[indu[i]];

		for (i = 0; i < num_indu; i++)
			free(temp_K[i]);
		free(temp_K);
	}

	if (num_indr > 0)
	{
		temp_K = (double **)malloc(sizeof(double*)*num_indr);
		for (i = 0; i < num_indr; i++)
			temp_K[i] = (double*)malloc(sizeof(double));

		kernel(indr, indc, num_indr, 1, temp_K,pm);
		for (i = 0; i < num_indr; i++)
			*f_c += temp_K[i][0] * (*pm).y[indr[i]] * (*pm).a[indr[i]];

		for (i = 0; i < num_indr; i++)
			free(temp_K[i]);
		free(temp_K);
	}
	free(indr);
	free(indu);
	free(indme);
}
