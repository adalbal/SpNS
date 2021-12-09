#include "FFTW_package.h"

//=====================================================================================================================
// Fast Fourier Transforms
//=====================================================================================================================
//r2c 3D FFT (input becomes garbage!)
void FFT_package::FFTWr2c(REAL *rin, COMPLEX *cout){
	mpi_execute_dft_r2c(r2c_plan, rin, cout);
};

//c2r 3D IFFT (input becomes garbage!)
void FFT_package::IFFTWc2r(COMPLEX *cin, REAL *rout){
	mpi_execute_dft_c2r(c2r_plan, cin, rout);
};

//=====================================================================================================================
// Normalization of "relevant" (due to dealiasing or non-cubic domains) components of vectors
//=====================================================================================================================
void FFT_package::Fourier_Normalization(REAL* phi) {
	LOOP_3D(i,0,local_n0, j,0,My, k,0,Mz) {
		int ind=(i*My+j)*2*(Mz_2+1)+k; //2*(Mz/2+1) instead of Mz because of mpi-FFTW padding!
		phi[ind] /= norm_Fourier;
	}
};
void FFT_package::Fourier_Normalization(REAL* phi, bool const * const relevant) {
	LOOP_3D(i,0,local_n0, j,0,My, k,0,Mz) {
		int ind=(i*My+j)*2*(Mz_2+1)+k; //2*(Mz/2+1) instead of Mz because of mpi-FFTW padding!
		if (relevant[ind]) {
			phi[ind] /= norm_Fourier;
		}
	}
};
void FFT_package::Fourier_Normalization(COMPLEX* phi, bool const * const relevant) {
	LOOP_3D(b,0,local_n1, a,0,Mx, k3,0,(Mz_2+1)) {
		int ind = (b*Mx+a)*(Mz_2+1)+k3;
		if (relevant[ind]) {
			for (int ic=0; ic<=1; ic++){ //ic=0 => Real part, ic=1 => Imaginary part
				phi[ind][ic] /= norm_Fourier;
			}
		}
	}
};

//=====================================================================================================================
// Initialization and destruction of FFTW_package class
//=====================================================================================================================
void FFT_package::FFT_init(){ // initialization of sizes and allocation of plans for FFT transforms
	//Get all mpi-related values
	alloc_local = mpi_local_size_3d_transposed(Mx, My, Mz_2+1, MCW, &local_n0, &local_0_start, &local_n1, &local_1_start);
	//Allocate memory
	complex_3D = alloc_complex(alloc_local);
	real_3D  = alloc_real(2*alloc_local);
	//Create plans
	r2c_plan = mpi_plan_dft_r2c_3d(Mx, My, Mz, real_3D, complex_3D, MCW, FMTO|FM|FDI); //Input preserved (by default r2c FFTW_DESTROY_RESERVED)
	c2r_plan = mpi_plan_dft_c2r_3d(Mx, My, Mz, complex_3D, real_3D, MCW, FMTI|FM); //Input destroyed (By default FFTW_DESTROY_INPUT)
};

void FFT_package::FFT_destroy(){
	FREE(complex_3D);
	FREE(real_3D);
};
