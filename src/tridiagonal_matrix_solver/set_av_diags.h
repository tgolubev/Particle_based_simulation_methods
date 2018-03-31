#ifndef SET_AV_DIAGS_H
#define SET_AV_DIAGS_H


double* set_main_diag(double *epsilon, double *a, int num_elements);
double* set_upper_diag(double *epsilon, double *b, int num_elements);
double* set_lower_diag(double *epsilon, double *c, int num_elements);

#endif // SET_AV_DIAGS_H
