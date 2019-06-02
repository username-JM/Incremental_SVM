#include "global.h"

void min_delta_acb(int indc, double* gamma, double* beta, int polc, int rflag, double* min_dacb,
	int*indss_s, int* cstatus_s, int* nstatus_s,struct parameter * pm)
{
	int indss[6];
	int cstatus[6];
	int nstatus[6];

	int i;
	for (i = 0; i < 6; i++)
	{
		indss[i] = 0;
		cstatus[i] = 0;
		nstatus[i] = 0;
	}

	for (i = 0; i < 3; i++)
	{
		indss[i] = indc;
	}

	double delta_m, delta_e, delta_r;

	if (polc == 1)
	{//incremental
		delta_m = (-1)*(*pm).g[indc] / gamma[indc];
		nstatus[0] = MARGIN;
		if ((*pm).num_m + 1 > 1)
		{
			delta_e = (*pm).C[indc] - (*pm).a[indc];
			nstatus[1] = ERROR;
		}
		else
		{
			delta_e = INF;
		}
		delta_r = INF;
	}
	else
	{//decremental
		if ((*pm).g[indc] > 0)
		{
			delta_m = (*pm).g[indc] / gamma[indc];
			nstatus[0] = MARGIN;
		}
		else
			delta_m = INF;
		if ((*pm).a[indc] <= (*pm).C[indc])
		{
			delta_e = INF;
			delta_r = (*pm).a[indc];
			if ((*pm).g[indc] > 0)
				nstatus[2] = RESERVE;
			else
				nstatus[2] = UNLEARNED;
		}
		else
		{
			delta_e = (*pm).a[indc] - (*pm).C[indc];
			delta_r = INF;
			nstatus[1] = ERROR;
		}
	}

	double beta_s[MAX_NUM];
	int flags[MAX_NUM] ;
	
	memset(beta_s, 0, sizeof(beta_s));
	memset(flags, 0, sizeof(flags));

	double delta_mer;
	int i_num;
	if ((*pm).num_m + 1 > 1)
	{//length(beta) of num is num_MV+1
		for (i = 1; i < (*pm).num_m + 1; i++)
		{
			beta_s[i-1] = polc*beta[i];
			if (beta_s[i-1] != 0)
				flags[i-1] = 1;
			else
				flags[i-1] = 0;
		}
		double send_C[MAX_NUM];
		memset(send_C, 0, sizeof(send_C));

		for (i = 0; i < (*pm).num_m; i++)
		{
			if (beta_s[i] > 0)
			{
				send_C[i] = (*pm).C[(*pm).ind[MARGIN][i]];
			}
		}

		double a_s[MAX_NUM];
		memset(a_s, 0, sizeof(a_s));

		for (i = 0; i < (*pm).num_m; i++)
		{
			a_s[i] = (*pm).a[(*pm).ind[MARGIN][i]];
		}

		min_delta(flags, a_s, send_C, beta_s, (*pm).num_m, &delta_mer, &i_num);

		if (delta_mer < INF)
		{
			indss[3] = (*pm).ind[MARGIN][i_num];
			cstatus[3] = MARGIN;

			if (beta_s[i_num] > 0)
			{
				nstatus[3] = ERROR;
			}
			else if (beta_s[i_num] < 0)
				nstatus[3] = RESERVE;
		}
	}
	else
		delta_mer = INF;

	double gamma_e[MAX_NUM];
	for (i = 0; i < (*pm).num_e; i++)
	{
		gamma_e[i] = polc*gamma[(*pm).ind[ERROR][i]];
		if (gamma_e[i] > 0)
		{
			flags[i] = 1;
		}
		else
		{
			flags[i] = 0;
		}
	}

	double delta_em;
	double send_z[MAX_NUM];
	memset(send_z, 0, sizeof(send_z));

	double tmp1[MAX_NUM];
	memset(tmp1, 0, sizeof(tmp1));

	for (i = 0; i < (*pm).num_e; i++)
	{
		tmp1[i] = (*pm).g[(*pm).ind[ERROR][i]];
	}

	min_delta(flags, tmp1, send_z, gamma_e, (*pm).num_e, &delta_em, &i_num);

	if (delta_em < INF) 
	{
		indss[4] = (*pm).ind[ERROR][i_num];
		cstatus[4] = ERROR;
		nstatus[4] = MARGIN;
	}

	double gamma_r[MAX_NUM];
	double delta_rm;
	if (rflag)
	{
		for (i = 0; i < (*pm).num_r; i++)
		{
			gamma_r[i] = polc*gamma[(*pm).ind[RESERVE][i]];
			if ((*pm).g[(*pm).ind[RESERVE][i]] >= 0 && gamma_r[i] < 0)
			{
				flags[i] = 1;
			}
			else
			{
				flags[i] = 0;
			}
		}
		double send_r[MAX_NUM];
		memset(send_r, 0, sizeof(send_r));

		double tmp2[MAX_NUM];
		memset(tmp2, 0, sizeof(tmp2));
		int tmp2_idx = 0;
		for (i = 0; i < (*pm).num_r; i++)
		{
			tmp2[i] = (*pm).g[(*pm).ind[RESERVE][i]];
		}

		min_delta(flags, tmp2, send_r, gamma_r, (*pm).num_r, &delta_rm, &i_num);

		if (delta_rm < INF)
		{
			indss[5] = (*pm).ind[RESERVE][i_num];
			cstatus[5] = RESERVE;
			nstatus[5] = MARGIN;
		}
	}
	else
		delta_rm = INF;

	int min_ind;
	double temp1, temp2;

	//앞에 세 개중 가장 작은
	if (delta_m < delta_e&&delta_m < delta_r)
		temp1 = delta_m;
	else if (delta_e < delta_m&&delta_e < delta_r)
		temp1 = delta_e;
	else
		temp1 = delta_r;

	//그 뒤 세 개중 가장 작은
	if (delta_mer < delta_em&&delta_mer < delta_rm)
		temp2 = delta_mer;
	else if (delta_em < delta_mer&&delta_em < delta_rm)
		temp2 = delta_em;
	else
		temp2 = delta_rm;

	if (temp1 > temp2)
		*min_dacb = temp2;
	else
		*min_dacb = temp1;

	//min_ind찾기
	if (*min_dacb == temp1)
	{
		if (temp1 == delta_m)
			min_ind = 0;
		else if (temp1 == delta_e)
			min_ind = 1;
		else
			min_ind = 2;
	}
	else
	{
		if (temp2 == delta_mer)
			min_ind = 3;
		else if (temp2 == delta_em)
			min_ind = 4;
		else
			min_ind = 5;
	}

	*indss_s = indss[min_ind];
	*cstatus_s = cstatus[min_ind];
	*nstatus_s = nstatus[min_ind];

	if (polc != 1)
		*min_dacb = (double)polc*(*min_dacb);

	return;
}