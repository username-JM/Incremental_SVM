#include "global.h"

void move_ind(int* inda, int* indb, int indc, int cstatus, int nstatus, struct parameter * pm);
int move_indr(int* inda, int indss, int cstatus,struct parameter * pm);

void bookkeeping(int indss, int cstatus, int nstatus, int *indco, int *i,struct parameter * pm)
{
	*indco = -1;
	*i = 0;

	if (cstatus != nstatus)
	{
		switch (nstatus) {
		case RESERVE:
			(*pm).a[indss] = 0; break;
		case ERROR:
			(*pm).a[indss] = (*pm).C[indss]; break;
		default:
			 break;
		}

		if (cstatus == MARGIN)
		{//change the matlab library function "find".
			int j = 0;
			for (j = 0; j< (*pm).num_m; j++)
			{//ind[MARGIN][j]!=NULL
				if (indss == (*pm).ind[MARGIN][j])
				{
					*indco = j + 1;
				}
			}
		}

		switch (nstatus)
		{
		case RESERVE:
			*i = move_indr((*pm).ind[cstatus], indss, cstatus,pm);
			(*pm).num_r = (*pm).num_r + *i; break;
		case ERROR:
			move_ind((*pm).ind[cstatus], (*pm).ind[nstatus], indss, cstatus, nstatus,pm); break;
		case MARGIN:
			move_ind((*pm).ind[cstatus], (*pm).ind[nstatus], indss, cstatus, nstatus,pm); break;
		case UNLEARNED:
			move_ind((*pm).ind[cstatus], (*pm).ind[nstatus], indss, cstatus, nstatus,pm); break;
		default:
			printf("error"); break;
		}
	}

	return;
}

int move_indr(int* inda, int indss, int cstatus, struct parameter *pm)
{
	int removed_i = 0;
	int num_RVs_orig = (*pm).num_r;

	move_ind(inda, (*pm).ind[RESERVE], indss, cstatus, RESERVE,pm);

	if ((*pm).num_r > (*pm).max_reserve_vectors) {
		int i, k;
		k = -1;
		int rem_idx[MAX_NUM] = { 0, };
		double * g_sorted = (double *)malloc(sizeof(double)*(*pm).num_r);

		for (i = 0; i < (*pm).num_r; i++) {
			g_sorted[i] = (*pm).g[(*pm).ind[RESERVE][i]];
			rem_idx[i] = i;
		}
		quicksort(g_sorted, rem_idx, 0, (*pm).num_r - 1);
		int removed = rem_idx[(*pm).num_r - 1];
		if (removed <= num_RVs_orig)
			k = 0;
		if (k == 0) 
			removed_i = removed	;
		if (removed == (*pm).num_r - 1) {
			(*pm).ind[RESERVE][removed] = -1;
			(*pm).num_r--;
		}
		else {
			for (i = removed; i < (*pm).num_r - 1; i++) {
				(*pm).ind[RESERVE][i] = (*pm).ind[RESERVE][i + 1];
			}
			(*pm).num_r--;
			free(g_sorted);
		}
	}
	return removed_i;
}

void move_ind(int* inda, int* indb, int indc, int cstatus, int nstatus, struct parameter * pm)
{
	//move_ind(ind[cstatus], ind[nstatus], indss, cstatus, nstatus);
	//indss is example changing status.
	int new_inds[MAX_NUM] = { 0, };
	int new_idx;
	int i;
	if (nstatus == RESERVE)
	{
		indb[(*pm).num_r] = indc;
		(*pm).num_r++;
		new_idx = 0;
		if (cstatus == RESERVE)
		{
			for (i = 0; i < (*pm).num_r; i++)
			{
				if (inda[i] != indc)
				{
					//add inda[i] to new_inds.
					new_inds[new_idx] = inda[i];
					new_idx++;
				}
			}
			(*pm).num_r--;
		}
		else if (cstatus == ERROR)
		{
			for (i = 0; i < (*pm).num_e; i++)
			{
				if (inda[i] != indc)
				{
					//add inda[i] to new_inds.
					new_inds[new_idx] = inda[i];
					new_idx++;
				}
			}
			(*pm).num_e--;
		}
		else if (cstatus == MARGIN)
		{
			for (i = 0; i < (*pm).num_m; i++)
			{
				if (inda[i] != indc)
				{
					//add inda[i] to new_inds.
					new_inds[new_idx] = inda[i];
					new_idx++;
				}
			}
			(*pm).num_m--;
		}
		else
		{
			for (i = 0; i < (*pm).num_u; i++)
			{
				if (inda[i] != indc)
				{
					//add inda[i] to new_inds.
					new_inds[new_idx] = inda[i];
					new_idx++;
				}
			}
			(*pm).num_u--;
		}
	}
	else if (nstatus == MARGIN)
	{
		indb[(*pm).num_m] = indc;
		(*pm).num_m++;

		int i;
		new_idx = 0;
		if (cstatus == RESERVE)
		{
			for (i = 0; i < (*pm).num_r; i++)
			{
				if (inda[i] != indc)
				{
					//add inda[i] to new_inds.
					new_inds[new_idx] = inda[i];
					new_idx++;
				}
			}
			(*pm).num_r--;
		}
		else if (cstatus == ERROR)
		{
			for (i = 0; i < (*pm).num_e; i++)
			{
				if (inda[i] != indc)
				{
					//add inda[i] to new_inds.
					new_inds[new_idx] = inda[i];
					new_idx++;
				}
			}
			(*pm).num_e--;
		}
		else if (cstatus == MARGIN)
		{
			for (i = 0; i < (*pm).num_m; i++)
			{
				if (inda[i] != indc)
				{
					//add inda[i] to new_inds.
					new_inds[new_idx] = inda[i];
					new_idx++;
				}
			}
			(*pm).num_m--;
		}
		else
		{
			for (i = 0; i < (*pm).num_u; i++)
			{
				if (inda[i] != indc)
				{
					//add inda[i] to new_inds.
					new_inds[new_idx] = inda[i];
					new_idx++;
				}
			}
			(*pm).num_u--;
		}
	}
	else if (nstatus == ERROR)
	{
		indb[(*pm).num_e] = indc;
		(*pm).num_e++;

		int i;
		new_idx = 0;
		if (cstatus == RESERVE)
		{
			for (i = 0; i < (*pm).num_r; i++)
			{
				if (inda[i] != indc)
				{
					//add inda[i] to new_inds.
					new_inds[new_idx] = inda[i];
					new_idx++;
				}
			}
			(*pm).num_r--;
		}
		else if (cstatus == ERROR)
		{
			for (i = 0; i < (*pm).num_e; i++)
			{
				if (inda[i] != indc)
				{
					//add inda[i] to new_inds.
					new_inds[new_idx] = inda[i];
					new_idx++;
				}
			}
			(*pm).num_e--;
		}
		else if (cstatus == MARGIN)
		{
			for (i = 0; i < (*pm).num_m; i++)
			{
				if (inda[i] != indc)
				{
					//add inda[i] to new_inds.
					new_inds[new_idx] = inda[i];
					new_idx++;
				}
			}
			(*pm).num_m--;
		}
		else
		{
			for (i = 0; i < (*pm).num_u; i++)
			{
				if (inda[i] != indc)
				{
					//add inda[i] to new_inds.
					new_inds[new_idx] = inda[i];
					new_idx++;
				}
			}
			(*pm).num_u--;
		}
	}
	else
	{
		indb[(*pm).num_u] = indc;
		(*pm).num_u++;

		int i;
		new_idx = 0;
		if (cstatus == RESERVE)
		{
			for (i = 0; i < (*pm).num_r; i++)
			{
				if (inda[i] != indc)
				{
					//add inda[i] to new_inds.
					new_inds[new_idx] = inda[i];
					new_idx++;
				}
			}
			(*pm).num_r--;
		}
		else if (cstatus == ERROR)
		{
			for (i = 0; i < (*pm).num_e; i++)
			{
				if (inda[i] != indc)
				{
					//add inda[i] to new_inds.
					new_inds[new_idx] = inda[i];
					new_idx++;
				}
			}
			(*pm).num_e--;
		}
		else if (cstatus == MARGIN)
		{
			for (i = 0; i < (*pm).num_m; i++)
			{
				if (inda[i] != indc)
				{
					//add inda[i] to new_inds.
					new_inds[new_idx] = inda[i];
					new_idx++;
				}
			}
			(*pm).num_m--;
		}
		else
		{
			for (i = 0; i < (*pm).num_u; i++)
			{
				if (inda[i] != indc)
				{
					//add inda[i] to new_inds.
					new_inds[new_idx] = inda[i];
					new_idx++;
				}
			}
			(*pm).num_u--;
		}
	}

	for (i = 0; i < new_idx; i++) {
		inda[i] = new_inds[i];
	}
}