#ifndef UTIL_HEADER
#define UTIL_HEADER

//=====================================================================================================================
// Preprocessor definitions
//=====================================================================================================================
#ifndef QA //QA can be set to true via Makefile doing: "make QA"
	#define QA 0 //By default, QA=false
#endif

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
// Worksharing constructs ("omp for"), not "omp parallel for": they assume an
// enclosing "#pragma omp parallel" is already open (opened once per iteration in
// New_Fractional_Step_Method) and just divide iterations among the existing team --
// no per-loop fork/join. With no enclosing parallel region (one-off init/I-O code),
// "omp for" just runs serially, which is correct, if not parallel, there.
#define OMP_SCHEDULE static
#define LOOP_2D(i0,n00,n01, i1,n10,n11) _Pragma("omp for collapse(2) schedule(OMP_SCHEDULE)") \
                                        for(int i0=n00;i0<n01;i0++) for(int i1=n10;i1<n11;i1++)
#define LOOP_3D(i0,n00,n01, i1,n10,n11, i2,n20,n21) _Pragma("omp for collapse(3) schedule(OMP_SCHEDULE)") \
                                                    for(int i0=n00;i0<n01;i0++) for(int i1=n10;i1<n11;i1++) for(int i2=n20;i2<n21;i2++)
// NOTE: OpenMP's worksharing constructs require canonical for-loop form, which
// rules out the old multi-variable comma-declaration style -- every call site must
// compute ind (and k1/k2 for the _k1k2k3 variant) itself as the first statement(s)
// of the loop body; see existing call sites for the pattern.
// collapse(1) (spreading work over local_n1/local_n0 only) beat collapse(2)/(3) at
// every scale tested -- more cache/prefetch-friendly than interleaved flat-index
// chunks. CAVEAT: if local_n1/local_n0 is ever smaller than the thread count, some
// threads get zero iterations; collapse(3) doesn't have this failure mode.
#define LOOP_FOURIER_k1k2k3 _Pragma("omp for collapse(1) schedule(OMP_SCHEDULE)") \
							for(int b=0; b<local_n1; b++) \
							for(int a=0; a<Mx; a++) \
							for(int k3=0; k3<(Mz_2+1); k3++)
#define LOOP_FOURIER _Pragma("omp for collapse(1) schedule(OMP_SCHEDULE)") \
					 for(int b=0; b<local_n1; b++) \
					 for(int a=0; a<Mx; a++) \
					 for(int k3=0; k3<(Mz_2+1); k3++)
// nowait variant: ONLY safe when the very next construct is another LOOP_FOURIER or
// LOOP_FOURIER_k1k2k3 (same trip count/schedule -> thread t provably gets the same
// index range both times) with nothing else -- no master/FFTW/MPI/differently-shaped
// loop -- in between.
#define LOOP_FOURIER_NOWAIT _Pragma("omp for collapse(1) schedule(OMP_SCHEDULE) nowait") \
					 for(int b=0; b<local_n1; b++) \
					 for(int a=0; a<Mx; a++) \
					 for(int k3=0; k3<(Mz_2+1); k3++)
#define LOOP_REAL _Pragma("omp for collapse(1) schedule(OMP_SCHEDULE)") \
				  for (int i=0; i<local_n0; i++) \
				  for (int j=0; j<My; j++) \
				  for (int k=0; k<Mz; k++)
// Reduction variant of LOOP_REAL: for loops that accumulate into a scalar declared
// outside the loop (e.g. a max-timestep). "clause" is whatever goes inside
// reduction(...), e.g. LOOP_REAL_REDUCE(max:local_MaxVel).
#define LOOP_REAL_REDUCE(clause) _Pragma(STR(omp for collapse(1) schedule(OMP_SCHEDULE) reduction(clause))) \
				  for (int i=0; i<local_n0; i++) \
				  for (int j=0; j<My; j++) \
				  for (int k=0; k<Mz; k++)
#define STR_HELPER(x) #x
#define STR(x) STR_HELPER(x)

// Copies 3 same-length COMPLEX arrays (the u/v/w velocity-triplet pattern used all
// over the hot path). Master-only memcpy + barrier, not split across the team via
// omp-for: measured consistently worse split (0.2-16% slower at every thread count
// on MN5) -- single-core memory bandwidth isn't the bottleneck here.
#define COPY3_COMPLEX(du,dv,dw, su,sv,sw, n) \
	_Pragma("omp master") \
	{ \
		memcpy((du), (su), sizeof(COMPLEX) * (n)); \
		memcpy((dv), (sv), sizeof(COMPLEX) * (n)); \
		memcpy((dw), (sw), sizeof(COMPLEX) * (n)); \
	} \
	_Pragma("omp barrier")
// nowait variant: like LOOP_FOURIER_NOWAIT, safe only when nothing reads/writes
// du/dv/dw before some later barrier already provides ordering (any omp barrier is
// a full team-wide flush, not just an ordering point for the data that triggered it).
#define COPY3_COMPLEX_NOWAIT(du,dv,dw, su,sv,sw, n) \
	_Pragma("omp master") \
	{ \
		memcpy((du), (su), sizeof(COMPLEX) * (n)); \
		memcpy((dv), (sv), sizeof(COMPLEX) * (n)); \
		memcpy((dw), (sw), sizeof(COMPLEX) * (n)); \
	}
// Debug macros
#define FUCK printf("rank %d ------> %s,\t %s(),\t line: %d\n", myrank, __FILE__, __FUNCTION__, __LINE__);

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
