/* This file was automatically generated.  Do not edit! */
#if 0
int main();
#else
struct tensor {    
  double v[3][3];    
  int magic;    
tensor():magic(0){};    
};
#endif

double S3QR_model(double vol,double QS,double QG,double RG,double V2,double nu);
double S3PR_model(double vol,double QS,double QG,double RG,double V2,double nu);
double S3PQ_model(double vol,double QS,double QG,double RG,double V2,double nu);
double S3PQR_model(double Cs3pqr,double vol,double QS,double QG,double RG,double V2,double p,double q,double r,double nu);
double Sigma_model(double vol,double QS,double QG,double RG,double V2,double nu);
double WALE_model(double vol,double QS,double QG,double V2,double nu);
double Vreman_model(double vol,double QS,double QG,double V2,double nu);
double Verstappen_model(double d[3],double RS,double QS,double nu);
double calc_lambda_Delta(double d[3]);
double Smag_model(double Cs,double vol,double Q,double nu);
double Smag_model(double vol,double QS,double nu);
void solve_cubic_eq(double a,double b,double c,double d,double x[3]);
double calc_RG(struct tensor *G);
double calc_RS(struct tensor *G);
double calc_QO(double QS,double QG);
double calc_V2(struct tensor *G);
double calc_QG(struct tensor *G);
double calc_QS(struct tensor *G);
double calc_third_invariant(struct tensor *T);
double calc_second_invariant(struct tensor *T);
double calc_first_invariant(struct tensor *T);
double calc_determinant(struct tensor *t);
void calc_prod_tensors(struct tensor *a,struct tensor *b,struct tensor *ab);
double calc_trace(struct tensor *t);
void calc_Omega(struct tensor *G,struct tensor *O);
void calc_S(struct tensor *G,struct tensor *S);
void calc_force_symmetry(struct tensor *A,struct tensor *symA);
void calc_transpose_tensor(struct tensor *A,struct tensor *At);
void cl_3tensors(double coeffA,struct tensor *A,double coeffB,struct tensor *B,double coeffC,struct tensor *C,struct tensor *ABC);
void cl_2tensors(double coeffA,struct tensor *A,double coeffB,struct tensor *B,struct tensor *AB);
void print_tensor(struct tensor *t);
void check_tensor(struct tensor *t);
void force_traceless(struct tensor *A);
void alloc_fill_rand_tensor(struct tensor *t,int traceless);
double generate_a_random_number(double lim0,double lim1);
double random_number(double lim0,double lim1);
void crash(const char *fmt,...);
int is_zero(double a);
