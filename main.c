#include "global.h"
int main() {
	int num = 0;
	double *new_X[2];
	int new_y[2];
	int new_model = 1;
	
	svmtrain(new_X,new_y,num,new_model);
}

