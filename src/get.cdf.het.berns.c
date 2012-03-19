#include <stdio.h>
#include <R.h>
#include <Rmath.h> 

void cdf(int *n, double *j, double *theta_n, double *alpha_n, double *result);
double subr(int n, double j, double *theta_n, double *alpha_n);

/* have subroutine to carry out recursion*/
double subr(int n, double j, double *theta_n, double *alpha_n)
{
    double sum_alphas = 0;
    for(int i = 1; i <= n; i++)
    {
	sum_alphas = sum_alphas + alpha_n[i-1];
    }
    if(n == 1)
    {
	if(j < 0)
	{
	    return 0;
	}
	else 
	    {
		if(j >= sum_alphas)
		{
		    return 1;
		}
		else
		{
		    return (1-theta_n[0]);
		}
	    }
    }
    else
    {
	/* create arrays which remove the last element */
	double theta_n_1[n-1];
	double alpha_n_1[n-1];
	for(int k = 1; k <= n-1; k++)
	{
	    theta_n_1[k-1] = theta_n[k-1];
	    alpha_n_1[k-1] = alpha_n[k-1];
	}
	int n_1 = n-1;
	double new_j = j-alpha_n[n-1];
	return (1-theta_n[n-1])*subr(n_1,j,theta_n_1,alpha_n_1) +
	    theta_n[n-1]*subr(n_1,new_j,theta_n_1,alpha_n_1);
    }
}

void cdf(int *n, double *j, double *theta_n, double *alpha_n, double *result)
{
    *result = subr(*n, *j, theta_n, alpha_n);
}
