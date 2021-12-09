#ifndef UTIL_HEADER
#define UTIL_HEADER

//=====================================================================================================================
// Preprocessor definitions
//=====================================================================================================================
#ifndef BYTES //BYTES can be defined via Makefile doing: "make float" or "make double"
	#define BYTES 8 //By default, double-precision
#endif
#if BYTES == 4
	#define REAL float
	#define COMPLEX fftwf_complex
	#define REAL_MPI MPI_FLOAT
	#define alloc_real fftwf_alloc_real
	#define alloc_complex fftwf_alloc_complex
	#define plan fftwf_plan
	#define execute fftwf_execute
	#define mpi_execute_dft_r2c fftwf_mpi_execute_dft_r2c
	#define mpi_execute_dft_c2r fftwf_mpi_execute_dft_c2r
	#define mpi_plan_dft_r2c_3d fftwf_mpi_plan_dft_r2c_3d
	#define mpi_plan_dft_c2r_3d fftwf_mpi_plan_dft_c2r_3d
	#define plan_dft_r2c_3d fftwf_plan_dft_r2c_3d
	#define destroy_plan fftwf_destroy_plan
	#define init_threads fftwf_init_threads
	#define mpi_init fftwf_mpi_init
	#define plan_with_nthreads fftwf_plan_with_nthreads
	#define mpi_local_size_3d_transposed fftwf_mpi_local_size_3d_transposed
	#define MALLOC fftwf_malloc
	#define FREE fftwf_free
#elif BYTES == 8
	#define REAL double
	#define COMPLEX fftw_complex
	#define REAL_MPI MPI_DOUBLE
	#define alloc_real fftw_alloc_real
	#define alloc_complex fftw_alloc_complex
	#define plan fftw_plan
	#define execute fftw_execute
	#define mpi_execute_dft_r2c fftw_mpi_execute_dft_r2c
	#define mpi_execute_dft_c2r fftw_mpi_execute_dft_c2r
	#define mpi_plan_dft_r2c_3d fftw_mpi_plan_dft_r2c_3d
	#define mpi_plan_dft_c2r_3d fftw_mpi_plan_dft_c2r_3d
	#define plan_dft_r2c_3d fftw_plan_dft_r2c_3d
	#define destroy_plan fftw_destroy_plan
	#define init_threads fftw_init_threads
	#define mpi_init fftw_mpi_init
	#define plan_with_nthreads fftw_plan_with_nthreads
	#define mpi_local_size_3d_transposed fftw_mpi_local_size_3d_transposed
	#define MALLOC fftw_malloc
	#define FREE fftw_free
#endif

#define MCW MPI_COMM_WORLD
#define PI 3.141592653589793238462643

// Variable Argument Macro (VA_MACRO) upto 3 arguments
#define NUM_ARGS_(_1, _2, _3, TOTAL, ...) TOTAL
#define NUM_ARGS(...) NUM_ARGS_(__VA_ARGS__, 3, 2, 1)
#define CONCATE_(X, Y) X##Y  // Fixed the double '_' from previous code
#define CONCATE(MACRO, NUMBER) CONCATE_(MACRO, NUMBER)
#define VA_MACRO(MACRO, ...) CONCATE(MACRO, NUM_ARGS(__VA_ARGS__))(__VA_ARGS__)
// Definition of FOURIER_FREQ
#define FOURIER_FREQ(...) VA_MACRO(FOURIER_FREQ, __VA_ARGS__)
#define FOURIER_FREQ2(a,M) ((a<=M/2)?(a):(a-M)) //Dimension non-split over processes
#define FOURIER_FREQ3(a,start,M) (((start+a)<=M/2)?(start+a):((start+a)-M)) //Dimension split over processes
// General LOOPS
#define LOOP_2D(i0,n00,n01, i1,n10,n11) for(int i0=n00;i0<n01;i0++) for(int i1=n10;i1<n11;i1++)
#define LOOP_3D(i0,n00,n01, i1,n10,n11, i2,n20,n21) for(int i0=n00;i0<n01;i0++) for(int i1=n10;i1<n11;i1++) for(int i2=n20;i2<n21;i2++)
#define LOOP_FOURIER_k1k2k3 for(int b=0, k2=FOURIER_FREQ(b,local_1_start,My); b<local_n1; b++, k2=FOURIER_FREQ(b,local_1_start,My)) \
							for(int a=0, k1=FOURIER_FREQ(a,Mx); a<Mx; a++, k1=FOURIER_FREQ(a,Mx)) \
							for(int k3=0, ind=(b*Mx+a)*(Mz_2+1)+k3; k3<(Mz_2+1); k3++, ind=(b*Mx+a)*(Mz_2+1)+k3)
#define LOOP_FOURIER for(int b=0; b<local_n1; b++) \
					 for(int a=0; a<Mx; a++) \
					 for(int k3=0, ind=(b*Mx+a)*(Mz_2+1)+k3; k3<(Mz_2+1); k3++, ind=(b*Mx+a)*(Mz_2+1)+k3)
#define LOOP_REAL for (int i=0; i<local_n0; i++) \
				  for (int j=0; j<My; j++) \
				  for (int k=0, ind=(i*My+j)*2*(Mz_2+1)+k; k<Mz; k++, ind=(i*My+j)*2*(Mz_2+1)+k)

using namespace std;

//=====================================================================================================================
// General functions
//=====================================================================================================================

int lcm (const int& n1, const int& n2); // returns least common multiple

void crash(const char *fmt,...); // abort execution with error
int pprintf(const char *fmt,...); //prints to process log
int pprintf0(const char *fmt,...); // prints to log and zero procs to stdout

int MyID(); // reports MPI id
int NumProc(); // reports MPI size

void open_logfile(void); //Opens output logfile
void init_parallel(const bool& isOpenMP, const int& OMP_Threads, int *argc, char ***argv); // init of MPI+OMP
void end_parallel(void); // finish MPI and IO 

#endif
