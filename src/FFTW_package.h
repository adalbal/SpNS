#ifndef __FFTW_PACKAGE__
#define __FFTW_PACKAGE__

#include <fftw3-mpi.h>

#include "util.h"

using namespace std;

//=====================================================================================================================
// Preprocessor definitions (FFTW related. The rest are in util.h)
//=====================================================================================================================
#define FE FFTW_ESTIMATE
#define FM FFTW_MEASURE
#define FP FFTW_PATIENT
#define FDI FFTW_DESTROY_INPUT
#define FMTI FFTW_MPI_TRANSPOSED_IN
#define FMTO FFTW_MPI_TRANSPOSED_OUT

//=====================================================================================================================
// FFT transforms package
//=====================================================================================================================
class FFT_package {
	public :
		FFT_package(const ptrdiff_t Mx_, const ptrdiff_t My_, const ptrdiff_t Mz_):
			Mx(Mx_), //Number of points to transform in each dimension (De-aliasing mesh size: M = 3*N / 2)
			My(My_),
			Mz(Mz_),
			Mz_2(Mz_/2),
			norm_Fourier(Mx_*My_*Mz_), //Total number of values transformed
			numprocs(NumProc()), //Total number of processes
			myrank(MyID()) //Rank of the current process
			{
				FFT_init();
			}; //def. constructor

		~FFT_package() {
			FFT_destroy();
		}; //def. destructor

		void FFT_init(); //initialization for all plans
		void FFT_destroy(); //destruction of all plans

		inline const ptrdiff_t& getalloc_local() const;
		inline const ptrdiff_t& getlocal_n0() const;
		inline const ptrdiff_t& getlocal_0_start() const;
		inline const ptrdiff_t& getlocal_n1() const;
		inline const ptrdiff_t& getlocal_1_start() const;

		// Both plans transpose their input (X-dim <-> Y-dim)
		void FFTWr2c(REAL *rin, COMPLEX *cout); //r2c 3D FFT (destroying input)
		void IFFTWc2r(COMPLEX *cin, REAL *rout); //c2r 3D IFFT (destroying input)
		//Normalization of transformed vectors (to be used wisely)
		void Fourier_Normalization(REAL* phi, bool const * const relevant); //Normalizes "relevant" components of real vectors
		void Fourier_Normalization(COMPLEX* phi, bool const * const relevant); //Normalizes "relevant" components of complex vectors
		void Fourier_Normalization(REAL* phi); //Normalizes all components of real vectors

	private :
		const int Mx, My, Mz, Mz_2;
		const REAL norm_Fourier;
		const int numprocs, myrank;
		ptrdiff_t alloc_local; //Total number of real values to be allocated by current process (2*alloc_local complex values)
		ptrdiff_t local_n0; //Number of different X-frequencies stored by current process - PHYSICAL
		ptrdiff_t local_0_start; //Value of smallest X-frequency stored by current process - PHYSICAL
		ptrdiff_t local_n1; //Number of different Y-frequencies stored by current process - FOURIER (transposed)
		ptrdiff_t local_1_start; //Value of smallest X-frequency stored by current process - FOURIER (transposed)

		plan r2c_plan, c2r_plan; //FFTW plans

		// Buffers for plan arrays
		COMPLEX *complex_3D;
		REAL *real_3D;
};

const ptrdiff_t& FFT_package::getalloc_local() const { return alloc_local;};
const ptrdiff_t& FFT_package::getlocal_n0() const { return local_n0;};
const ptrdiff_t& FFT_package::getlocal_0_start() const { return local_0_start;};
const ptrdiff_t& FFT_package::getlocal_n1() const { return local_n1;};
const ptrdiff_t& FFT_package::getlocal_1_start() const { return local_1_start;};


#endif //__FFTW_PACKAGE__

