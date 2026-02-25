#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h> 
#include <stdarg.h>    

#define PI  3.14159265358979323846
#define PI2 6.28318530717958647692

#define magic_tensor   239843

#define MYZERO 1e-14

#include "tensor_algebra.hpp"

int is_zero(double a) {
  if(fabs(a)<MYZERO) return(1); else return(0);
}

#if 0
void crash(const char *fmt,...){   
    va_list ap;   

    va_start(ap,fmt);   
    vprintf(fmt,ap);   
    va_end(ap);   

    exit(0);   
}
#endif

double random_number(double lim0,double lim1) {
  return(lim0+(lim1-lim0)*((double)rand())/((double)RAND_MAX));
}

double generate_a_random_number(double lim0,double lim1) {
  double r;
  double THRESHOLD=2e-7;
  
  r=random_number(lim0,lim1);
  while(2.0*fabs(r)/(fabs(lim0)+fabs(lim1))<THRESHOLD) r=random_number(lim0,lim1);

  return(r);
}

/************************************************************************************************/    
/*************************** TENSOR ALGEBRA *****************************************************/    
/************************************************************************************************/    

#if 0
struct tensor {    
  double v[3][3];    
  int magic;    
tensor():magic(0){};    
};    
#endif

void fill_rand_tensor(struct tensor *t,int traceless) {
  check_tensor(t);

  for(int i=0;i<=2;i++) {
    for(int j=0;j<=2;j++) {
      t->v[i][j]=generate_a_random_number(-1.0,1.0);
    }	
  }
  
  t->magic=magic_tensor;

  if(traceless) force_traceless(t);  
}

void alloc_fill_rand_tensor(struct tensor *t,int traceless) {
  int i,j;

  if(t->magic==magic_tensor) crash("fill_rand_tensor\n");

  for(i=0;i<=2;i++) {
    for(j=0;j<=2;j++) {
      t->v[i][j]=generate_a_random_number(-1.0,1.0);
    }	
  }
  
  t->magic=magic_tensor;

  if(traceless) force_traceless(t);  
}

void check_tensor(struct tensor * /*t*/) {
  //int i;
  //double trace=0.0;
  //if(t->magic!=magic_tensor) crash("check_tensor\n");

  //for(i=0;i<=2;i++) trace+=t->v[i][i];

  //if(!is_zero(trace)) crash("check_tensor\t it is not traceless!!!\n");
}

void print_tensor(struct tensor *t) {
  int i;

  check_tensor(t);
  printf("\n");
  for(i=0;i<=2;i++) printf("%e %e %e\n",t->v[i][0],t->v[i][1],t->v[i][2]);
}

void cl_2tensors(double coeffA,struct tensor *A,
		 double coeffB,struct tensor *B,
		 struct tensor *AB) {
  int i,j;

  check_tensor(A);
  check_tensor(B);

  for(i=0;i<=2;i++) {
    for(j=0;j<=2;j++) {
      AB->v[i][j]=coeffA*(A->v[i][j])+coeffB*(B->v[i][j]);
    }
  }

  AB->magic=magic_tensor;
}

void cl_3tensors(double coeffA,struct tensor *A,
		 double coeffB,struct tensor *B,
		 double coeffC,struct tensor *C,
		 struct tensor *ABC) {
  int i,j;

  check_tensor(A);
  check_tensor(B);
  check_tensor(C);

  for(i=0;i<=2;i++) {
    for(j=0;j<=2;j++) {
      ABC->v[i][j]=coeffA*(A->v[i][j])+coeffB*(B->v[i][j])+coeffC*(C->v[i][j]);
    }
  }

  ABC->magic=magic_tensor;
}

void calc_transpose_tensor(struct tensor *A,struct tensor *At) {
  int i,j;

  check_tensor(A);

  for(i=0;i<=2;i++) {
    for(j=0;j<=2;j++) {
      At->v[i][j]=A->v[j][i];
    }
  }

  At->magic=magic_tensor;
}

void calc_force_symmetry(struct tensor *A,struct tensor *symA) {
  struct tensor At;
  check_tensor(A);

  calc_transpose_tensor(A,&At);

  cl_2tensors(0.5,A,0.5,&At,symA); 
}

void calc_S(struct tensor *G,struct tensor *S) {
  check_tensor(G);

  calc_force_symmetry(G,S);
}

void calc_Omega(struct tensor *G,struct tensor *O) {
  struct tensor S;
  check_tensor(G);

  calc_S(G,&S);

  cl_2tensors(1.0,G,-1.0,&S,O);
}

double calc_trace(struct tensor *t) {
  int i;
  double trace=0.0;  //inicialitzem
  check_tensor(t);
		
  for(i=0;i<=2;i++) trace+=t->v[i][i];

  return(trace);
}

void calc_prod_tensors(struct tensor *a,struct tensor *b,struct tensor *ab) {
  int i,j,k;
	
  check_tensor(a);
  check_tensor(b);

  for(i=0;i<=2;i++) {
    for(j=0;j<=2;j++) {
      ab->v[i][j]=0.0;  //inicialitzem
      for(k=0;k<=2;k++) ab->v[i][j]+=a->v[i][k]*b->v[k][j];	
    }
  }

  ab->magic=magic_tensor;
}

void force_traceless(struct tensor *A) {
  int i;
  double trace;
  check_tensor(A);

  trace=calc_trace(A);

  for(i=0;i<=2;i++) A->v[i][i]-=(trace/3.0);
}

double calc_determinant(struct tensor *t) {
  double det=0.0;

  det+=t->v[0][0]*t->v[1][1]*t->v[2][2];
  det+=t->v[0][1]*t->v[1][2]*t->v[2][0];
  det+=t->v[0][2]*t->v[1][0]*t->v[2][1];
	
  det-=t->v[0][2]*t->v[1][1]*t->v[2][0];
  det-=t->v[0][1]*t->v[1][0]*t->v[2][2];
  det-=t->v[0][0]*t->v[1][2]*t->v[2][1];	
	
  return(det);
}

double calc_first_invariant(struct tensor *T) {
  return(calc_trace(T));
}

double calc_second_invariant(struct tensor *T) {
  double trT,trT2;
  struct tensor T2;

  calc_prod_tensors(T,T,&T2);

  trT=calc_trace(T);
  trT2=calc_trace(&T2);

  return(0.5*(trT*trT-trT2));
}

double calc_third_invariant(struct tensor *T) {
  return(calc_determinant(T));
}

double calc_QS(struct tensor *G) {
  struct tensor S;

  check_tensor(G);
  calc_S(G,&S);

  return(calc_second_invariant(&S));
}

double calc_QG(struct tensor *G) {

  check_tensor(G);

  return(calc_second_invariant(G));
}

double calc_V2(struct tensor *G) {   /*XAVI95.x2*/
  struct tensor S,O,S2,O2,S2O2;   /*XAVI95.x2*/
  double trS2O2,QS,QG,QO;   /*XAVI95.x2*/

  check_tensor(G);   /*XAVI95.x2*/

  calc_S(G,&S);   /*XAVI95.x2*/
  calc_Omega(G,&O);   /*XAVI95.x2*/

  calc_prod_tensors(&S,&S,&S2);     /*XAVI95.x2*/
  calc_prod_tensors(&O,&O,&O2);     /*XAVI95.x2*/

  calc_prod_tensors(&S2,&O2,&S2O2);     /*XAVI95.x2*/

  trS2O2=calc_trace(&S2O2);   /*XAVI95.x2*/

  QS=calc_QS(G);   /*XAVI95.x2*/
  QG=calc_QG(G);   /*XAVI95.x2*/
  QO=calc_QO(QS,QG);   /*XAVI95.x2*/

  return(4.0*(trS2O2-2.0*QS*QO));   /*XAVI95.x2*/
}   /*XAVI95.x2*/

double calc_RS(struct tensor *G) {
  struct tensor S;

  check_tensor(G);
  calc_S(G,&S);
 
  return(calc_third_invariant(&S));
}

double calc_RG(struct tensor *G) {

  check_tensor(G);
 
  return(calc_third_invariant(G));
}

/*a x^3 + b x^2 + c x + d = 0*/
/*The roots are already ordered: x[0]<=x[1]<=x[2]*/   
void solve_cubic_eq(double a,double b,double c,double d,double x[3]){
    double A,B,C;
    double R,Q;
    double sqrtQ;   
    double theta;

    /*normalitzem*/
    A=b/a;
    B=c/a;
    C=d/a;

    Q=(A*A-3.0*B)/9.0;
    R=(2.0*A*A*A - 9.0*A*B + 27*C)/54.0;

    /*In this case, we assume that A=0.0 and instead of solving the cubic equation
      x^3+Ax^2+Bx+C=0
      we solve the quadratic equation
      x^2+Ax+B=0*/   
    if(fabs(C/A)<1e-11) {     /*before the threshold was 1e-13*/ /*XAVI95.x3*/  
      double sol1=0.5*(-A+sqrt(A*A-4.0*B));   
      double sol0=0.5*(-A-sqrt(A*A-4.0*B));   
      if(sol0>0.0) {   
	x[0]=sol1; x[1]=sol0; x[2]=0.0;     
	return;   
      }   
      
      if(sol1<0.0) {   
	x[0]=0.0; x[1]=sol1; x[2]=sol0;     
	return;   
      }   

      x[0]=sol1; x[1]=0.0; x[2]=sol0;   

      return;   
    }   

    if(R*R>=Q*Q*Q) crash("solve_cubic_eq\t does not have real roots  R=%e Q=%e  R^2-Q^3=%e a=%e b=%e c=%e d=%e theta=%e\n",   
			 R,Q,R*R-Q*Q*Q,a,b,c,d,acos(R/pow(Q,1.5)));   

    sqrtQ=sqrt(Q);   

    theta=acos(R/pow(sqrtQ,3.0));   

    //In this way the roots are already ordered x[0]>=x[1]>=x[2]   
    x[2]=-2.0*sqrtQ*cos(theta/3.0      )-A/3.0;   
    x[0]=-2.0*sqrtQ*cos((theta+PI2)/3.0)-A/3.0;   
    x[1]=-2.0*sqrtQ*cos((theta-PI2)/3.0)-A/3.0;   
}

double calc_QO(double QS,double QG) {    
  return(QG-QS);    
}    

/************************************************************************************************/    
/*************************** EDDY-VISCOSITY *****************************************************/    
/************************************************************************************************/    

double Smag_model(double Cs,double vol,double Q,double nu) {  
  return(nu+Cs*Cs*pow(vol,2.0/3.0)*sqrt(-4.0*Q));  
} 
 
double calc_lambda_Delta(double d[3]) {   
  return(-(4.0)*(1.0/(d[0]*d[0])+1.0/(d[1]*d[1])+1.0/(d[2]*d[2])));  
}  

double Verstappen_model(double d[3],double RS,double QS,double nu) {  
  double lambda_Delta=calc_lambda_Delta(d);  
  double nue=1.5*(1.0/lambda_Delta)*fabs(RS)/QS;

  if(fabs(QS)>1e-15 && fabs(RS)>1e-15) return(nu+nue);  
  else                                 return(nu);  
}  

double Vreman_model(double vol,double QS,double QG,double V2,double nu) {      
  double QO=calc_QO(QS,QG);    
  double CVr=0.26457513110645905905;  /*sqrt(0.07)=0.26457513110645905905*/    
  double nue;    

  if(fabs(QG)<1e-15) return(nu);  
  
  nue=(nu+CVr*CVr*pow(vol,2.0/3.0)*sqrt((QG*QG+V2)/(2.0*(QO-QS))));       /*XAVI95.x2*/
  
  return(nue);    
}      

double WALE_model(double vol,double QS,double QG,double V2,double nu) {    
  double CW=0.5; 
  double aux,nue;    

  aux=0.5*V2+2.0/3.0*QG*QG;       /*XAVI95.x2*/

  if(fabs(aux)<1e-15) return(0.0);  /*XAVI95.x3*/

  nue=(nu+CW*CW*pow(vol,2.0/3.0)*(pow(aux,1.5)/(pow(-2.0*QS,2.5)+pow(aux,1.25))));    

  return(nue);    
}    

double Smag_model(double vol,double QS,double nu) {  
  double Cs=0.17;

  return(nu+Cs*Cs*pow(vol,2.0/3.0)*sqrt(-4.0*QS));  
}  

double Sigma_model(double vol,double QS,double QG,double RG,double V2,double nu) {    
  double CSig=1.5;    
  double P,Q,R;  //the three invariants of the GG^t tensor     
  double Pscaled,Qscaled,Rscaled;  //the three invariants of the the scaling*GG^t tensor    
  double eigenval[3];  //eigenvalues of the GG^t tensor    
  double eigenval_scaled[3];  //eigenvalues of the scaling*GG^t tensor    
  double sigma[3];    
  double nue;    
  double scaling;    
  int checking=1;   
  int aa;    
  
  P=2.0*(calc_QO(QS,QG)-QS);    
  Q=V2+QG*QG;    
  R=RG*RG;    
  
  if(Q<1e-15 && R<1e-15) return(nu);    /*XAVI95.x3*/
  
  scaling=1.0/P;   
  
  Pscaled=scaling*P;   
  Qscaled=scaling*scaling*Q;   
  Rscaled=scaling*scaling*scaling*R;   

  /*solving the cubic equation with the invariants of scaling*GG^t;
    the results are the ordered eigenvalues of the scaling*GG^t tensor*/   
  solve_cubic_eq(1.0,-Pscaled,Qscaled,-Rscaled,eigenval_scaled);    

  /*computing the eigenvalues of GG^t*/   
  for(aa=0;aa<=2;aa++) eigenval[aa]=eigenval_scaled[aa]/scaling;   

  /*computing the singular eigenvalues of G*/    
  for(aa=0;aa<=2;aa++) sigma[aa]=sqrt(eigenval[aa]);    

  /*checkings*/    
  if(checking) {   
    for(aa=0;aa<=1;aa++) {    
      if(sigma[aa]<sigma[aa+1]) crash("Alguna cosa a fallat amb les sigmes!!!\n");    
    }   

    double errorP=fabs(P-(eigenval[0]+eigenval[1]+eigenval[2]));   
    double errorQ=fabs(Q-(eigenval[0]*eigenval[1]+eigenval[0]*eigenval[2]+eigenval[1]*eigenval[2]));   
    double errorR=fabs(R-(eigenval[0]*eigenval[1]*eigenval[2]));

    if(fabs(errorP/P)>1e-7 && errorP>1e-9) crash("Something wrong with P!!! P=%e Q=%e R=%e\n",P,Q,R);    
    if(fabs(errorQ/P)>1e-7 && errorQ>1e-7) crash("Something wrong with Q!!! P=%e Q=%e R=%e l1*l2+l1*l3+l2*l3=%e\n",P,Q,R,(eigenval[0]*eigenval[1]+eigenval[0]*eigenval[2]+eigenval[1]*eigenval[2]));    
    if(fabs(errorR/P)>1e-7 && errorR>1e-5) printf("Something wrong with R!!! P=%e Q=%e R=%e l1*l2*l3=%e\n",P,Q,R,(eigenval[0]*eigenval[1]*eigenval[2]));    
  }   

  nue=(nu+CSig*CSig*pow(vol,2.0/3.0)*((sigma[2]*(sigma[0]-sigma[1])*(sigma[1]-sigma[2]))/(sigma[0]*sigma[0])));    

  if(fabs(nue-nue)<=MYZERO) crash("nue=%e nu=%e   eigenval=%e %e %e   P=%e Q=%e R=%e\n",nue,nu,eigenval[0],eigenval[1],eigenval[2],P,Q,R);    

  return(nue);    
}    

double S3PQR_model(double Cs3pqr,double vol,   
		   double QS,double QG,double RG,double V2,   
		   double p,double q,double r,   
		   double nu) {       
  double QO=calc_QO(QS,QG);    
  double P=2.0*(QO-QS);   
  double Q=V2+QG*QG;      /*XAVI95.x2*/
  double R=RG*RG;   
  double nue;   

  if(fabs(QG)<1e-15) return(nu);  

  nue=(nu+Cs3pqr*Cs3pqr*pow(vol,2.0/3.0)*pow(P,p)*pow(Q,q)*pow(R,r));   

  return(nue);   
}   

double S3PQ_model(double vol,double QS,double QG,double RG,double V2,double nu) {   
  double Cs3pq=0.45825756949558400066;  /*sqrt(3)*sqrt(0.07)=0.45825756949558400066, where sqrt(0.07) is Vreman's constant*/ 
  double p,q,r;   
  double nue;   

  //S3PQ model   
  p=-2.5; q=1.5; r=0.0;   
  
  nue=S3PQR_model(Cs3pq,vol,QS,QG,RG,V2,p,q,r,nu);   

  return(nue);   
}   

double S3PR_model(double vol,double QS,double QG,double RG,double V2,double nu) {   
  double Cs3pr=0.45825756949558400066;  /*sqrt(3)*sqrt(0.07)=0.45825756949558400066, where sqrt(0.07) is Vreman's constant*/ 
  double p,q,r;   
  double nue;   

  //S3PR model   
  p=-1.0; q=0.0; r=0.5;   
  
  nue=S3PQR_model(Cs3pr,vol,QS,QG,RG,V2,p,q,r,nu);   

  return(nue);   
}   

double S3QR_model(double vol,double QS,double QG,double RG,double V2,double nu) {   
  double Cs3qr=0.45825756949558400066;  /*sqrt(3)*sqrt(0.07)=0.45825756949558400066, where sqrt(0.07) is Vreman's constant*/ 
  double p,q,r;   
  double nue;   

  //S3QR model   
  p=0.0; q=-1.0; r=0.83333333333333333333333333333; /*5/6*/   
  
  //Vreman's model    
  //p=-0.5; q=0.5; r=0.0;   

  nue=S3PQR_model(Cs3qr,vol,QS,QG,RG,V2,p,q,r,nu);   

  return(nue);   
}   

/***********************************************************************************************/
/***********************************************************************************************/
/***********************************************************************************************/

void test_models(int ntests) {
  struct tensor G;
  double QG,QS,RG,RS,V2;
  double nueSmag,nueWALE,nueVre,nueSig,nueVer,nueS3PQ,nueS3PR,nueS3QR;
  double dissSmag,dissWALE,dissVre,dissSig,dissVer,dissS3PQ,dissS3PR,dissS3QR;
  double dx=0.01;
  double nu=1e-5;  
  double vol=dx*dx*dx;
  double d[3]={dx,dx,dx};

  printf("Testing models...\n");

  dissSmag=dissWALE=dissVre=dissSig=dissVer=dissS3PQ=dissS3PR=dissS3QR=0.0;  /*initialization*/
  
  alloc_fill_rand_tensor(&G,1/*traceless*/);

  
  for(int n=1;n<=ntests;n++) {
    fill_rand_tensor(&G,1/*traceless*/);
    
    QG=calc_QG(&G);
    QS=calc_QS(&G);
    RG=calc_RG(&G);
    RS=calc_RS(&G);
    V2=calc_V2(&G);
    
    nueSmag=Smag_model(vol,QS,nu);
    nueWALE=WALE_model(vol,QS,QG,V2,nu);
    nueVre=Vreman_model(vol,QS,QG,V2,nu);
    nueSig=Sigma_model(vol,QS,QG,RG,V2,nu);
    nueVer=Verstappen_model(d,RS,QS,nu);
    nueS3PQ=S3PQ_model(vol,QS,QG,RG,V2,nu);
    nueS3PR=S3PR_model(vol,QS,QG,RG,V2,nu);
    nueS3QR=S3QR_model(vol,QS,QG,RG,V2,nu);

    dissSmag-=nueSmag*QS;
    dissWALE-=nueWALE*QS;
    dissVre -=nueVre *QS;
    dissSig -=nueSig *QS;
    dissVer -=nueVer *QS;
    dissS3PQ-=nueS3PQ*QS;
    dissS3PR-=nueS3PR*QS;
    dissS3QR-=nueS3QR*QS;

    if(nueVer>nueSmag) crash("nueVer=%e nueSmag=%e\n",nueVer,nueSmag);
    if(nueS3PQ>nueVre) crash("nueS3PQ=%e nueVre=%e\n",nueS3PQ,nueVre);
    if(nueS3PR>nueVre) crash("nueS3PR=%e nueVre=%e\n",nueS3PR,nueVre);
    if(nueS3QR>nueVre) crash("nueS3QR=%e nueVre=%e\n",nueS3QR,nueVre);
  }
  
  dissSmag/=(double)ntests;
  dissWALE/=((double)ntests);
  dissVre /=((double)ntests);
  dissSig /=((double)ntests);
  dissVer /=((double)ntests);
  dissS3PQ/=((double)ntests);
  dissS3PR/=((double)ntests);
  dissS3QR/=((double)ntests);
  
  printf("Models tested! :-)\n");
  
  printf("dissSmag=%e dissWALE=%e\n",dissSmag,dissWALE);
  printf("dissWALE/dissSmag=%e\n",dissWALE/dissSmag);
  printf("dissVre /dissSmag=%e\n",dissVre /dissSmag);
  printf("dissSig /dissSmag=%e\n",dissSig /dissSmag);
  printf("dissVer /dissSmag=%e\n",dissVer /dissSmag);
  printf("dissS3PQ/dissSmag=%e\n",dissS3PQ/dissSmag);
  printf("dissS3PR/dissSmag=%e\n",dissS3PR/dissSmag);
  printf("dissS3QR/dissSmag=%e\n",dissS3QR/dissSmag);
}

#if 0
int main() {
  struct tensor G;
  double QG,QS,RG,RS,V2;
  double nueSmag,nueWALE,nueVre,nueSig,nueVer,nueS3PQ,nueS3PR,nueS3QR;
  double dx=0.01;
  double nu=1e-5;  
  double vol=dx*dx*dx;
  double d[3]={dx,dx,dx};

  srand (time(NULL));		

  test_models(10000000);

  alloc_fill_rand_tensor(&G,1/*traceless*/);
  
  QG=calc_QG(&G);
  QS=calc_QS(&G);
  RG=calc_RG(&G);
  RS=calc_RS(&G);
  V2=calc_V2(&G);

  nueSmag=Smag_model(vol,QS,nu);
  nueWALE=WALE_model(vol,QS,QG,V2,nu);
  nueVre=Vreman_model(vol,QS,QG,V2,nu);
  nueSig=Sigma_model(vol,QS,QG,RG,V2,nu);
  nueVer=Verstappen_model(d,RS,QS,nu);
  nueS3PQ=S3PQ_model(vol,QS,QG,RG,V2,nu);
  nueS3PR=S3PR_model(vol,QS,QG,RG,V2,nu);
  nueS3QR=S3QR_model(vol,QS,QG,RG,V2,nu);

  print_tensor(&G);
  printf("QG=%e QS=%e RG=%e RS=%e V2=%e\n",QG,QS,RG,RS,V2);
  printf("nueSmag=%e nueWALE=%e nueVre=%e nueSig=%e nueVer=%e nueS3PQ=%e nueS3PR=%e nueS3QR=%e\n",nueSmag,nueWALE,nueVre,nueSig,nueVer,nueS3PQ,nueS3PR,nueS3QR);

  return(0);
}
#endif
