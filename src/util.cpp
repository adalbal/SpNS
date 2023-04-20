#include "config.h"

#include <stdio.h>
#include <stdlib.h> 
#include <stdarg.h>
#include <mpi.h>
#include <fftw3-mpi.h>

#include "util.h"

int threads_ok;
#if QA
    static int shutitup = 0; //Avoid output of logfile
#else
    static int shutitup = 1; //Force output of logfile
#endif
static int info = 0; //if 1, then additional info is written in log file
static int init = 0; //initially 0 and, immediately after mpi has been initiated, it is set to 1 (and pprintf starts writing in each log files)
static int numprocs, myid;
static FILE *afile;

static int NoParallel(void){ return 0; } //unused yet 

//It returns least common multiple: lcm = (n1*n2)/hcf
int lcm (const int& n1, const int& n2) {
	int hcf = n1;
	int temp = n2;

	while (hcf != temp) {
		if (hcf > temp) hcf -= temp;
		else temp -= hcf;
	}
	return (n1 * n2) / hcf; //hcf: highest common factor
}

int MyID(){return myid;}
int NumProc(){return numprocs;}

void open_logfile (void) {
	if (!shutitup) {
		char fname[1000];
		sprintf(fname,"./%s.%d","stdout",myid);
		afile=fopen(fname,"w"); //text
		if(afile==NULL) {
			fprintf(stderr,"4\n");
			exit(-1);
		}
	}
}
static void close_logfile(void){
	if (!shutitup) {
		fclose(afile);
	}
}
static void abort_paralel(void){
	if(NoParallel()) return;
	MPI_Abort(MPI_COMM_WORLD,-1);
}
int pprintf0(const char *fmt,...){
	int r;
	va_list ap;
	if(shutitup||!init) return 0;

	va_start(ap,fmt);
	r=vfprintf(afile,fmt,ap);
	va_end(ap);

	if(myid==0){
		va_start(ap,fmt);
		r=vfprintf(stdout,fmt,ap);
		va_end(ap);
	}
 
	if(myid==0)fflush(stdout); 
	fflush(afile);
 
	return(r);
}
int pprintf(const char *fmt,...){
	int r;
	va_list ap;
	if(shutitup||!init) return 0;

	va_start(ap,fmt);

	r=vfprintf(afile,fmt,ap);
	va_end(ap);

	fflush(afile);
	return(r);
}
 
void crash(const char *fmt,...){
	FILE *f=NULL;
	va_list ap;
	int i;

	fflush(NULL);
	for(i=0; i<=2; i++){ /* treiem el missatge per la sortida d'errors i per la sortida estandar */
		if(i==2 && afile==stdout) continue; /* si utilitzem stdout, unicament per la sortida d'errors */
		switch(i){
		case 0: f=stdout; break;
		case 1: f=stderr; break;
		case 2: f=afile;  break;
		}
		fprintf(f,"crash:quisoc=%d: ",myid);
		va_start(ap,fmt);
		vfprintf(f,fmt,ap);
		va_end(ap);
		fprintf(f,"\n");
		fflush(f);
	}

	abort_paralel();
	exit(0);
}
 
void barrier(MPI_Comm c){
	int ok=MPI_Barrier(c);
	if(ok!=MPI_SUCCESS) crash("barrier/1");
}
 
void init_parallel(const bool& isOpenMP, const int& OMP_Threads, int *argc,char ***argv){
	if (NoParallel()) {
		myid = 0;
		numprocs = 1;
	} else if(isOpenMP) {
		//Multi-threaded MPI FFTW initialization
		int provided;
		MPI_Init_thread(argc, argv, MPI_THREAD_FUNNELED, &provided);
		threads_ok = provided >= MPI_THREAD_FUNNELED;
		if (threads_ok) threads_ok = init_threads();
		mpi_init();
		if (threads_ok) plan_with_nthreads(OMP_Threads); //From here, all FFTW3 plans will consider (at most) OMP_Threads threads
		//Rest of MPI related variables
		if(MPI_Comm_rank(MPI_COMM_WORLD,&myid) != MPI_SUCCESS){fprintf(stderr,"2\n"); exit(-1); };
		if(MPI_Comm_size(MPI_COMM_WORLD,&numprocs) != MPI_SUCCESS){fprintf(stderr,"3\n"); exit(-1); };
	} else {
		//MPI-only FFTW initialization
		if(MPI_Init(argc,argv) != MPI_SUCCESS){fprintf(stderr,"1\n"); exit(-1); }
		mpi_init();
		//Rest of MPI related variables
		if(MPI_Comm_rank(MPI_COMM_WORLD,&myid) != MPI_SUCCESS){fprintf(stderr,"2\n"); exit(-1); };
		if(MPI_Comm_size(MPI_COMM_WORLD,&numprocs) != MPI_SUCCESS){fprintf(stderr,"3\n"); exit(-1); };
	}

	if(NoParallel()) return;

	init=1;

	open_logfile(); //Opens log file where pprintf will write

	if(info) pprintf("Output initiated. Process %d \n", myid);

	barrier(MPI_COMM_WORLD);
	if(info) pprintf("Barrier and initial checkint done\n");

	{
		int q;
		if(*argc) {
            if(info) pprintf("Arguments:\n");
			for(q=0; q<=*argc-1; q++) {
                if(info) pprintf("q=%d argv[q]=<%s> \n", q, (*argv)[q]);
            }
		} else{
            if(info) pprintf("No arguments\n");
        }
	}
}

void end_parallel(void){
	if(info) pprintf("Entered end_parallel()\n");
	barrier(MPI_COMM_WORLD);
	if(info) pprintf("Barrier and final checkint done\n");
	close_logfile();
	MPI_Finalize();
}
