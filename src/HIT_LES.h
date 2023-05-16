#ifndef __HIT_LES__
#define __HIT_LES__

#include <iostream>
#include <cmath>
#include <cstring> //memcpy
#include <fstream> //to read input file
#include<functional> //function

#include "fftw3-mpi.h"

#include "util.h"
#include "FFTW_package.h"

using namespace std;

//=====================================================================================================================
// Homogeneous Isotropic Turbulence class (implementing LES)
//=====================================================================================================================
class HIT {
	public:
		HIT(const ptrdiff_t& Nx_, const ptrdiff_t& Ny_, const ptrdiff_t& Nz_,
			const ptrdiff_t& Mx_, const ptrdiff_t& My_, const ptrdiff_t& Mz_,
			const REAL& nu_, const REAL& C_Smag_, const REAL& C_At_,
			const int& Lx_factor_, const int& Ly_factor_, const int& Lz_factor_,
			const REAL& omega_x_, const REAL& omega_y_, const REAL& omega_z_,
			const char *ASCII_Input_Filename_, const char *Binary_Input_Filename_, const ptrdiff_t& Nx_file_, const ptrdiff_t& Ny_file_, const ptrdiff_t& Nz_file_,
			const bool& isReLambda_, const bool& isInitialFieldFromBinaryFile_, const bool& isInitialFieldFromASCIIFile_, const bool& isNotCubic_, const bool& isRotating_, const bool& isLES_, const bool& isComplexConjugateCorrected_, const bool& isSelfAdaptive_):
			fftw(Mx_,My_,Mz_), //Class to compute FFT and IFFT
			ActualNx((Nx_%2==0) ? Nx_/Lx_factor_ : (Nx_-1)/Lx_factor_+1), //Number of dealiased modes to be resolved considering length factors
			ActualNy((Ny_%2==0) ? Ny_/Ly_factor_ : (Ny_-1)/Ly_factor_+1),
			ActualNz((Nz_%2==0) ? Nz_/Lz_factor_ : (Nz_-1)/Lz_factor_+1),
			ActualMx((Mx_%2==0) ? Mx_/Lx_factor_ : (Mx_-1)/Lx_factor_+1), //Number of aliased modes to be resolved considering length factors
			ActualMy((My_%2==0) ? My_/Ly_factor_ : (My_-1)/Ly_factor_+1),
			ActualMz((Mz_%2==0) ? Mz_/Lz_factor_ : (Mz_-1)/Lz_factor_+1),
			Nx(Nx_), //Number of dealiased modes to be calculated
			Ny(Ny_),
			Nz(Nz_),
			Nx_2(Nx_/2),
			Ny_2(Ny_/2),
			Nz_2(Nz_/2),
			Mx(Mx_), //De-aliasing mesh size = 3*N/2
			My(My_),
			Mz(Mz_),
			Mx_2(Mx_/2), //Half de-aliasing mesh size = N/2
			Mz_2(Mz_/2),
			last_rad(min((Nx_-1)/2, min((Ny_-1)/2, (Nz_-1)/2))), //Radium of the biggest circumference incribed in a Nx*Ny*Nz parallelogram
			numprocs(NumProc()), //Total number of MPI processes
			myrank(MyID()), //Current MPI process ID
			nu(nu_), //Kinetic viscosity
			delta(2.0*PI/cbrt(Lx_factor_*Ly_factor_*Lz_factor_*ActualMx*ActualMy*ActualMz)), //LES-model characteristic length
			C_Smag(C_Smag_), //LES-model constant
			C_At(C_At_), //Scaling factor to manually tune timestep
			Ax(2.0*PI/ActualMx/Lx_factor_),
			Ay(2.0*PI/ActualMy/Ly_factor_),
			Az(2.0*PI/ActualMz/Lz_factor_),
			Lx_factor(Lx_factor_), //Length factors
			Ly_factor(Ly_factor_),
			Lz_factor(Lz_factor_),
			omega_x(omega_x_), //Angular rotation velocity of the system
			omega_y(omega_y_),
			omega_z(omega_z_),
			ASCII_Input_Filename(ASCII_Input_Filename_), //ASCII velocity input filename
			Binary_Input_Filename(Binary_Input_Filename_), //Binary velocity input filename
			ActualNx_file((Nx_file_%2==0) ? Nx_file_/Lx_factor_ : (Nx_file_-1)/Lx_factor_+1), //Number of dealiased modes to be resolved considering length factors
			ActualNy_file((Ny_file_%2==0) ? Ny_file_/Ly_factor_ : (Ny_file_-1)/Ly_factor_+1),
			ActualNz_file((Nz_file_%2==0) ? Nz_file_/Lz_factor_ : (Nz_file_-1)/Lz_factor_+1),
			Nx_file(Nx_file_), //Augmented mesh size of the velocity input file (field of velocityes repeated as many times as )
			Ny_file(Ny_file_),
			Nz_file(Nz_file_),
			isReLambda(isReLambda_), //if true, then nu_ is undefined and setnu() needs to be called (in MAIN.cpp)
			isInitialFieldFromBinaryFile(isInitialFieldFromBinaryFile_),
			isInitialFieldFromASCIIFile(isInitialFieldFromASCIIFile_),
			isNotCubic(isNotCubic_),
			isRotating(isRotating_),
			isLES(isLES_),
			isComplexConjugateCorrected(isComplexConjugateCorrected_),
			isSelfAdaptive(isSelfAdaptive_),
			time(0.0) //Simulated time
		{
			HIT_init();
		}; //def. constructor

		~HIT() {
			HIT_destroy();
		}; //def. destructor

		void HIT_init(); //initialization of class
		void HIT_destroy(); //destruction of class

		//===================================================
		// INPUTTING INITIAL VELOCITY FIELD
		//===================================================
		//=========================
		// A) Initialization of velocity field in Fourier space using K41 distibution |k|^{-5/3}
		void Input_K41_Field();
		//=========================
		// B) Initialization of velocity field in Fourier space from an input file containing velocities in Fourier space
		void Input_Real_Field();
			// B.1) Load from file (a velocity field in real space)
			void Read_Input_File_Vxyz (REAL* u_file, REAL* v_file, REAL* w_file, const int& dumb_columns);
			// B.2) Truncate/Extend input velocity field in Fourier space (accordingly to Nx, Ny, Nz)
			void Truncate_Input_File_Fourier_Coefficients (COMPLEX const * const ukx_file, COMPLEX const * const uky_file, COMPLEX const * const ukz_file);

		//===================================================
		// FRACTIONAL STEP METHOD
		//===================================================
		void New_Fractional_Step_Method();
		//=========================
		// 1.) R VECTOR:
		//Convective term in Fourier space
		void Recalculate_Velocity_Gradients_Fourier_Coefficients();
		void Recalculate_Velocity_Antitransform();
		void Recalculate_Convective_Fourier_Coefficients();
		//Divergence of SGS tensor in Fourier space (LES modelling)
		void Recalculate_Eddy_Viscosity();
		void Recalculate_SGS_Tensor_Components();
		void Recalculate_SGS_Tensor_Transform();
		void Recalculate_SGS_Tensor_Divergence_Fourier_Coefficients();
		//R vector
		void Recalculate_R_vector_Fourier_Coefficients();
		//=========================
		// 2.) PREDICTOR VELOCITY:
		void Recalculate_Predictor_Velocity_Fourier_Coefficients();
		//=========================
		// 3.) DIVERGENCE-FREE PROJECTION OF PREDICTOR VELOCITY:
		void Recalculate_Divergence_Free_Projection();

		//===================================================
		// CALCULATE NEW TIMESTEP:
		//===================================================
		//Obtention of new timestep
		void Recalculate_TimeStep();
		// A) Calculation of new self-adaptive timestep
		void Recalculate_SelfAdaptive_TimeStep();
			//Auxiliary related functions
			inline REAL FuL(REAL x, REAL x0, REAL x1, REAL f0, REAL f1);
			inline REAL FuQ(REAL x, REAL x0, REAL x1);
			inline REAL FuG(REAL x, REAL a, REAL b, REAL c, REAL x0, REAL x1, REAL f0, REAL f1);
			inline REAL Topt(REAL phi);
			inline REAL Kopt(REAL phi);
		// B) Calculation of new CFL timestep
		void Recalculate_CFL_TimeStep();

		//===================================================
		// FORCING TERM
		//===================================================
		// Imposition of a given Energy distribution ("normalizing" current Fourier coefficients)
		void Forcing_Energy_Cascade(REAL const * const Forced_Ek, const ptrdiff_t& last_input_rad, const bool& isNullifyMissingEk);
			// Imposition of correct complex conjugation for k3=0
			void Complex_Conjugate_Correction();
		// Calculation of Reynolds lambda
		REAL Recalculate_Reynolds_Lambda();
		// Imposition of a given Reynolds lambda
		void Forcing_Reynolds_Lambda(const REAL ReLambda_);
			// Calculation of Enstrophy
			REAL Recalculate_Enstrophy(const char* filename = NULL);

		//===================================================
		// POST-PROCESS
		//===================================================
		// Integrate complex fields over spherical shells from k=1 to k=last_rad
		REAL Integrate_Field(const char* filename, std::function<REAL(int a, int b, int k3)>& funcFieldNorm, REAL* inAcumField);
		inline REAL Integrate_Field(const char* filename, std::function<REAL(int a, int b, int k3)>& funcFieldNorm);
		// Calculation of kinetic energy cascade
		REAL Recalculate_Energy_Cascade(const char* filename = NULL);
		void Calculate_Ek_init_file(COMPLEX const * const uk_file, COMPLEX const * const vk_file, COMPLEX const * const wk_file);
		// I/O functions (implemented in HIT_LES_output.cpp)
		void DealiasedComplex3DimField_to_BinaryFile(COMPLEX const * const field_x, COMPLEX const * const field_y, COMPLEX const * const field_z, char const * const filename) const;
		void Real3DimField_to_BinaryFile(REAL const * const field_x, REAL const * const field_y, REAL const * const field_z, char const * const filename) const;
		void FourierVelocity_to_BinaryFile(char const * const filename) const;
		void RealVelocity_to_BinaryFile(char const * const filename) const;
		void Field2File(REAL const * const field, char const * const filename) const;
		void Field2File(COMPLEX const * const field, char const * const filename) const;

		//===================================================
		// GETS, SETS, EXTRACTS
		//===================================================
		inline const ptrdiff_t& getNx() const;
		inline const ptrdiff_t& getNy() const;
		inline const ptrdiff_t& getNz() const;
		inline const ptrdiff_t& getlast_rad() const;
		inline const REAL& getnu() const;
		inline void setnu(const REAL& nu_);
		inline const REAL& getReLambda() const;
		inline const REAL& getAt() const;
		inline const REAL& gettime() const;
		inline const REAL& getEk_Tot();
		inline const REAL& getEk_init() const;
		inline const REAL& getEk_init_file() const;
		inline const REAL& getEnstrophy_init_file() const; //Actual dissipation: epsilon = nu*Enstrophy
		inline const REAL& getEk(ptrdiff_t K) const;

	private:
		FFT_package fftw;
		const ptrdiff_t ActualNx, ActualNy, ActualNz, ActualMx, ActualMy, ActualMz, Nx, Ny, Nz, Nx_2, Ny_2, Nz_2, Mx, My, Mz, Mx_2, Mz_2, last_rad;
		const int numprocs, myrank;
		REAL nu, ReLambda;
		const REAL delta, C_Smag;
		const REAL C_At, Ax, Ay, Az; // Any timestep
		const int Lx_factor, Ly_factor, Lz_factor;
		const REAL omega_x, omega_y, omega_z;
		const char *ASCII_Input_Filename, *Binary_Input_Filename;
		const ptrdiff_t ActualNx_file, ActualNy_file, ActualNz_file, Nx_file, Ny_file, Nz_file;
		const bool isReLambda,isInitialFieldFromBinaryFile,isInitialFieldFromASCIIFile,isNotCubic,isRotating,isLES,isComplexConjugateCorrected,isSelfAdaptive;
		REAL time;
		REAL MaxVel, At; //Any timestep
		REAL C_visc, C_conv; //CFL timestep
		REAL Ah, t1, kappa, *k, *phi, *c, phi0, _kappa05; //Self-adaptive timestep
		//General case
		ptrdiff_t alloc_local, local_n0, local_0_start, local_n1, local_1_start;
		REAL *u, *v, *w;
		REAL *gradu[3], *gradv[3], *gradw[3];
		REAL *convx, *convy, *convz;
		COMPLEX *uk_aux, *vk_aux, *wk_aux, *uk_0, *vk_0, *wk_0, *uk_1, *vk_1, *wk_1;
		COMPLEX *gradu_k[3], *gradv_k[3], *gradw_k[3];
		COMPLEX *convx_k, *convy_k, *convz_k;
		COMPLEX *Rx0_k, *Ry0_k, *Rz0_k, *Rx1_k, *Ry1_k, *Rz1_k;
		bool *dealiased, *conjugate;
		REAL *rad2, *local_acumField, *global_acumField, *Ek;
		REAL Ek_Tot, Ek_init, Ek_init_file, Enstrophy_init_file;
		int num_dealiased, num_dealiased_real;
		int *local_n0_, *local_n1_;
		MPI_Offset offset_Fourier, offset_Real;
		//LES
		REAL *nu_eddy, *SGS_11, *SGS_21, *SGS_22, *SGS_31, *SGS_32, *SGS_33;
		COMPLEX *SGS_11_k, *SGS_21_k, *SGS_22_k, *SGS_31_k, *SGS_32_k, *SGS_33_k, *divSGSx_k, *divSGSy_k, *divSGSz_k;
		//Complex conjugate correction
		COMPLEX *local_uk_k3_0, *local_vk_k3_0, *local_wk_k3_0, *uk_k3_0, *vk_k3_0, *wk_k3_0;
		int *displs_k3_0, *recvcounts_k3_0;
		// Complex field integrals
		std::function<REAL(int a, int b, int k3)> KineticEnergy;
		std::function<REAL(int a, int b, int k3)> Enstrophy;
	};

//Auxiliary functions related to the calculation of self-adaptive timestep
inline REAL HIT::FuL(REAL x, REAL x0, REAL x1, REAL f0, REAL f1){
	return f0+(x-x0)*(f1-f0)/(x1-x0);
}
inline REAL HIT::FuQ(REAL x, REAL x0, REAL x1){
	return (x-x0)*(x-x1);
}
inline REAL HIT::FuG(REAL x, REAL a, REAL b, REAL c, REAL x0, REAL x1, REAL f0, REAL f1){
	return (a*x*x+b*x+c)*FuQ(x,x0,x1)+FuL(x,x0,x1,f0,f1);
}
inline REAL HIT::Topt(REAL phi0){
	if(phi0>=0.0&&phi0<phi[0]) return FuG(phi0,0.0,c[0],c[1],0.0,phi[0],4.0/3.0,t1);
	if(phi0>=phi[0]&&phi0<=PI/2.0) return FuG(phi0,c[2],c[3],c[4],phi[0],PI/2.0,t1,1.0);
	crash("Topt fail!"); return 0;
}
inline REAL HIT::Kopt(REAL phi0){
	if(phi0>=0.0&&phi0<=phi[0]) return 1.0;
	if(phi0>phi[0]&&phi0<=phi[1]) return FuG(phi0,c[5],c[6],c[7],phi[0],phi[1],1.0,k[0]);
	if(phi0>phi[1]&&phi0<=phi[2]) return FuG(phi0,c[8],c[9],c[10],phi[1],phi[2],k[0],k[1]);
	if(phi0>phi[2]&&phi0<=PI/2.0) return FuG(phi0,c[11],c[12],c[13],phi[2],PI/2.0,k[1],0.0);
	crash("Kopt fail!"); return 0;
}
//Post-process
REAL HIT::Integrate_Field(const char* filename, std::function<REAL(int a, int b, int k3)>& funcFieldNorm) {
	return Integrate_Field(filename, funcFieldNorm, global_acumField);
};
//Gets, sets and extracts
const ptrdiff_t& HIT::getNx() const { return Nx;};
const ptrdiff_t& HIT::getNy() const { return Ny;};
const ptrdiff_t& HIT::getNz() const { return Nz;};
const ptrdiff_t& HIT::getlast_rad() const { return last_rad;};
const REAL& HIT::getnu() const { return nu;};
void HIT::setnu(const REAL& nu_) { nu=nu_;};
const REAL& HIT::getReLambda() const { return ReLambda;};
const REAL& HIT::getAt() const { return At;};
const REAL& HIT::gettime() const { return time;};
const REAL& HIT::getEk_Tot() { return Ek_Tot;};
const REAL& HIT::getEk_init() const { return Ek_init;};
const REAL& HIT::getEk_init_file() const { return Ek_init_file;};
const REAL& HIT::getEnstrophy_init_file() const { return Enstrophy_init_file;};
const REAL& HIT::getEk(ptrdiff_t K) const { return Ek[K];};

#endif // __HIT_LES__
