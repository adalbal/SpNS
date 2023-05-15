//spectral box Time step self adaptative parallel LES
#include "config.h"

#include <iostream>
#include <cmath>

#include <time.h>
#include <stdlib.h>
#include <stdarg.h>
#include <fftw3-mpi.h>

#include "util.h"
#include "parser.h"
#include "HIT_LES.h"

using namespace std;

//=====================================================================================================================
// Input functions
//=====================================================================================================================
void Read_Ek_Input_File (const int& Last_K, REAL* Ek_input, const int& dumb_columns, const char* EkFilename);

//=====================================================================================================================
// main program
//=====================================================================================================================
int main (int argc, char **argv){
	//===================================================
	// PARAMETERS
	//===================================================
	clock_t timebeg = clock(); //Initial time
	// -----------------------------------------------------------
	// File parameters -------------------------------------------
	//OBLIGATORY:
	//Mesh size
	int Nx, Ny, Nz;
	//Kinematic viscosity
	double nu;
	//Reynolds lambda
	double ReLambda;
	int ReLambda_Freq;
	//Last instant of time to be simulated
	double Final_Time;
	//------------
	//OPTIONAL:
	//Outputting obtained results
	char VelPhysfolder[256], VelFourfolder[256], Ekfolder[256];
	int Velocity_Physical, Velocity_Fourier, Energy_Cascade;  // (1:Output field / elsewhere:No output of the field)
	int Velocity_Phys_Freq, Velocity_Four_Freq, Energy_Cascade_Freq; // Output frequency (in terms of iteration number)
	//Inputting velocity in phycical space from file
	char ASCII_Input_Filename[256];
	char Binary_Input_Filename[256];
	int Nx_file, Ny_file, Nz_file;
	//Inputting energy cascade to be forced from file
	char EkFilename[256];
	int last_input_rad;
	int nullify_missingEk;
	double Final_Forced_Time;
	//Extra constants
	int Lx_factor, Ly_factor, Lz_factor; //Length factors
	double omega_x, omega_y, omega_z; //Angular rotation velocity of the system
	double C_Smag; //Smagorinsky model constant
	double C_At; //Timestep scaling factor
	int SelfAdaptiveTimestep; //Temporal scheme (1:self-adaptive / elsewhere:CFL+Adams-Bashforth)
	int ComplexConjugateCorrection; //Correction/Omission of complex conjugates at plane k3=0 (z=0)
	int OMP_Threads; //Number of OpenMP threads assigned to each MPI-process
	// -----------------------------------------------------------
	// Derived parameters ----------------------------------------
	bool isReLambda, isReLambda_iter, isInitialFieldFromBinaryFile, isInitialFieldFromASCIIFile, isForcedEk, isNullifyMissingEk, isNotCubic, isRotating, isLES, isSelfAdaptive, isComplexConjugateCorrected, isOpenMP;
	bool isVelPhysOut, isVelPhysOut_iter, isVelFourOut, isVelFourOut_iter, isEkOut, isEkOut_iter;
	int Mx, My, Mz; //De-aliasing mesh size = 3*(N/2) != (3*N)/2 !!!!
	// -----------------------------------------------------------
	// Output parameters -----------------------------------------
	//Final filenames (result of concatenating strings)
	char VelPhysfilename[512];
	char VelFourfilename[512];
	char Ekfilename[512];
	//Output files can be easily stored in folders modifying base_* variables
	const string base_VelPhysfilename = "Velocity";
	const string base_VelFourfilename = "Velocity_Fourier";
	const string base_Ekfilename = "Energy_Cascade";
	//File extensions
	const string ASCII_fileformat = "dat"; //ASCII
	const string Binary_fileformat = "bin"; //Binary
	int VelPhys_file_iter = 0;
	int VelFour_file_iter = 0;
	int Ek_file_iter = 0;
	// -----------------------------------------------------------

	{ //read of number of OMP threads per process in order to initiate parallelization (and parallel output of log file)
		char fparam[128];
		if(argc>1) sprintf(fparam, "%s", argv[1]);
		else sprintf(fparam, "%s","./parameters");

		FILE *File; //All-purpose file descriptor
		File=fopen(fparam,"r");
		if(File==NULL){ crash("Parameters file is missing!!! Include \"%s\", please\n", fparam); return 0;}
		fclose(File);

		CParamInitManager PM;
		FileBuffer FB(fparam, IO_CRASH);
		PM.RequestParameter(OMP_Threads,"OMP_Threads",TYPE_INT,IO_DONTCRASH,GT,0); //Number of OMP threads assigned to each MPI process
		PM.ReadParamsFromBuffer(FB);

		if (PM["OMP_Threads"].GetIsSet() && (OMP_Threads > 1)) {
			isOpenMP = true;
		} else {
			isOpenMP = false;
		}
	}

	//Beginning of parallelization
	init_parallel(isOpenMP, OMP_Threads, &argc, &argv);
	const int myrank = MyID();
	const int numprocs = NumProc();
	{ //Parameter read from file "parameters" and initialization of derived variables
		char fparam[128];
		if(argc>1) sprintf(fparam, "%s", argv[1]);
		else sprintf(fparam, "%s","./parameters");

		FILE *File; //All-purpose file descriptor
		File=fopen(fparam,"r");
		if(File==NULL){ crash("Parameters file is missing!!! Include \"%s\", please\n", fparam); return 0;}
		fclose(File);

		CParamInitManager PM;
		FileBuffer FB(fparam, IO_CRASH);

		PM.RequestParameter(Nx,"Nx",TYPE_INT,IO_CRASH,GT,0);
		PM.RequestParameter(Ny,"Ny",TYPE_INT,IO_CRASH,GT,0);
		PM.RequestParameter(Nz,"Nz",TYPE_INT,IO_CRASH,GT,0);
		PM.RequestParameter(nu,"nu",TYPE_DOUBLE,IO_DONTCRASH,GT,0.0);
		PM.RequestParameter(ReLambda,"Reynolds_lambda",TYPE_DOUBLE,IO_DONTCRASH,GT,0.0);
		PM.RequestParameter(ReLambda_Freq,"Reynolds_lambda_Freq",TYPE_INT,IO_DONTCRASH,GE,0);
		PM.RequestParameter(Final_Time,"Final_Time",TYPE_DOUBLE,IO_CRASH,GT,0.0);

		PM.RequestParameter(Velocity_Physical,"Velocity_Physical",TYPE_INT,IO_DONTCRASH); //(1=true/elsewhere=false)
		PM.RequestParameter(Velocity_Fourier,"Velocity_Fourier",TYPE_INT,IO_DONTCRASH); //(1=true/elsewhere=false)
		PM.RequestParameter(Energy_Cascade,"Energy_Cascade",TYPE_INT,IO_DONTCRASH); //(1=true/elsewhere=false)
		PM.RequestParameter(Velocity_Phys_Freq,"Velocity_Phys_Freq",TYPE_INT,IO_DONTCRASH,GE,0);
		PM.RequestParameter(Velocity_Four_Freq,"Velocity_Four_Freq",TYPE_INT,IO_DONTCRASH,GE,0);
		PM.RequestParameter(Energy_Cascade_Freq,"Energy_Cascade_Freq",TYPE_INT,IO_DONTCRASH,GE,0);
		PM.RequestParameter(VelPhysfolder,"Velocity_Phys_Folder",TYPE_WORD,IO_DONTCRASH);
		PM.RequestParameter(VelFourfolder,"Velocity_Four_Folder",TYPE_WORD,IO_DONTCRASH);
		PM.RequestParameter(Ekfolder,"Energy_Cascade_Folder",TYPE_WORD,IO_DONTCRASH);

		PM.RequestParameter(ASCII_Input_Filename,"ASCII_Input_Filename",TYPE_WORD,IO_DONTCRASH);
		PM.RequestParameter(Binary_Input_Filename,"Binary_Input_Filename",TYPE_WORD,IO_DONTCRASH);
		PM.RequestParameter(Nx_file,"Nx_file",TYPE_INT,IO_DONTCRASH,GT,0);
		PM.RequestParameter(Ny_file,"Ny_file",TYPE_INT,IO_DONTCRASH,GT,0);
		PM.RequestParameter(Nz_file,"Nz_file",TYPE_INT,IO_DONTCRASH,GT,0);

		PM.RequestParameter(EkFilename,"InputEkFilename",TYPE_WORD,IO_DONTCRASH);
		PM.RequestParameter(last_input_rad,"Last_K",TYPE_INT,IO_DONTCRASH,GT,0);
		PM.RequestParameter(nullify_missingEk,"Nullify_MissingK",TYPE_INT,IO_DONTCRASH,GE,0);
		PM.RequestParameter(Final_Forced_Time,"Final_Forced_Time",TYPE_DOUBLE,IO_DONTCRASH,GE,0.0);

		PM.RequestParameter(Lx_factor,"Lx",TYPE_INT,IO_DONTCRASH,GT,0);
		PM.RequestParameter(Ly_factor,"Ly",TYPE_INT,IO_DONTCRASH,GT,0);
		PM.RequestParameter(Lz_factor,"Lz",TYPE_INT,IO_DONTCRASH,GT,0);
		PM.RequestParameter(omega_x,"omega_x",TYPE_DOUBLE,IO_DONTCRASH);
		PM.RequestParameter(omega_y,"omega_y",TYPE_DOUBLE,IO_DONTCRASH);
		PM.RequestParameter(omega_z,"omega_z",TYPE_DOUBLE,IO_DONTCRASH);
		PM.RequestParameter(C_Smag,"C_Smag",TYPE_DOUBLE,IO_DONTCRASH,GE,0.0);
		PM.RequestParameter(C_At,"C_At",TYPE_DOUBLE,IO_DONTCRASH,GT,0.0);
		PM.RequestParameter(SelfAdaptiveTimestep,"SelfAdaptiveTimestep",TYPE_INT,IO_DONTCRASH); //(1=true/elsewhere=false)
		PM.RequestParameter(ComplexConjugateCorrection,"ComplexConjugateCorrection",TYPE_INT,IO_DONTCRASH); //(1=true/elsewhere=false)

		PM.ReadParamsFromBuffer(FB);

		if (PM["nu"].GetIsSet() && PM["Reynolds_lambda"].GetIsSet()) crash("Kinematic viscosity and Reynolds lambda cannot be defined simoultaniously!!!\n");
		if (!PM["nu"].GetIsSet() && !PM["Reynolds_lambda"].GetIsSet()) crash("Either kinematic viscosity or Reynolds lambda has to be defined!!!\n");

		Mx = (3*Nx)/2; My = (3*Ny)/2; Mz = (3*Nz)/2;

		if (PM["Reynolds_lambda"].GetIsSet()) {
			isReLambda = true;
			isReLambda_iter = ((PM["Reynolds_lambda_Freq"].GetIsSet() && (ReLambda_Freq > 0)) ? true : false);
		} else {
			isReLambda = false;
			isReLambda_iter = false;
		}

		if (PM["Nx_file"].GetIsSet() && PM["Ny_file"].GetIsSet() && PM["Nz_file"].GetIsSet()) {
			if (PM["Binary_Input_Filename"].GetIsSet()) {
				isInitialFieldFromBinaryFile = true;
				isInitialFieldFromASCIIFile = false;
				if (PM["ASCII_Input_Filename"].GetIsSet()) {
					crash("\"%s\" and \"%s\" have been defined as input files!!\n", Binary_Input_Filename, ASCII_Input_Filename);
				}
			} else if (PM["ASCII_Input_Filename"].GetIsSet()) {
				isInitialFieldFromASCIIFile = true;
				isInitialFieldFromBinaryFile = false;
			} else {
				isInitialFieldFromASCIIFile = false;
				isInitialFieldFromBinaryFile = false;
			}
		} else {
			isInitialFieldFromASCIIFile = false;
			isInitialFieldFromBinaryFile = false;
			Nx_file = -1; //unnecessary, but to avoid valgrind complains in HIT's ctor
			Ny_file = -1; //unnecessary, but to avoid valgrind complains in HIT's ctor
			Nz_file = -1; //unnecessary, but to avoid valgrind complains in HIT's ctor
		}
		if (PM["InputEkFilename"].GetIsSet() && PM["Last_K"].GetIsSet() && PM["Final_Forced_Time"].GetIsSet()) {
			isForcedEk = true;
			isNullifyMissingEk = (PM["Nullify_MissingK"].GetIsSet() && nullify_missingEk ? true : false);
			Final_Forced_Time = ((Final_Forced_Time > 0.0) ? Final_Forced_Time : Final_Time);
		} else {
			isForcedEk = false;
			isNullifyMissingEk = false;
		}
		if ((PM["Lx"].GetIsSet() && (Lx_factor > 1)) ||
			(PM["Ly"].GetIsSet() && (Ly_factor > 1)) ||
			(PM["Lz"].GetIsSet() && (Lz_factor > 1))) {
			// In case only significant length factor has been set
			Lx_factor = ((PM["Lx"].GetIsSet()) ? Lx_factor : 1);
			Ly_factor = ((PM["Ly"].GetIsSet()) ? Ly_factor : 1);
			Lz_factor = ((PM["Lz"].GetIsSet()) ? Lz_factor : 1);
			int lcm_Lxyz = lcm(Lx_factor, lcm(Ly_factor, Lz_factor)); //Least common multiple of Lx, Ly and Lz
			// The following change is made on the length factors: Let's assume that (Lx,Ly,Lz) = (L,3L,12L),
			// then we assume that new length is L'=12L (the biggest factor). Thus, (Lx,Ly,Lz) = (L'/12,L'/4,L').
			// This allows us to implement different length factors by setting to zero all coefficients such
			// that k1%12!=0 and k2%4!=0 (bearing that Kx_k1=2PIk1/Ax and Ky_k2=2PIk2/Ay).
			Lx_factor = lcm_Lxyz / Lx_factor;
			Ly_factor = lcm_Lxyz / Ly_factor;
			Lz_factor = lcm_Lxyz / Lz_factor;
			if ((Lx_factor > 1) || (Ly_factor > 1) || (Lz_factor > 1)) {
				isNotCubic = true;
				//Redefine "ActualM" (see HIT_LES.h) to avoid odd values of M. Then, ActualM=2*k => M=k*L_factor.
				//Consequently, the physical domain will be coherently divided into L_factor slices of k nodes, and
				//the velocity field will be repeated (k times) in that direction (only relevant for L_factor>1)
				Mx = ((Lx_factor>1 && Mx%2!=0) ? Mx+1 : Mx);
				My = ((Ly_factor>1 && My%2!=0) ? My+1 : My);
				Mz = ((Lz_factor>1 && Mz%2!=0) ? Mz+1 : Mz);
				//Update accordingly Nx,Ny,Nz and Mx,My,Mz (and, if needed Nx_file,Ny_file,Nz_file)
				Nx = ((Nx%2==0) ? Nx*Lx_factor : (Nx-1)*Lx_factor+1);
				Ny = ((Ny%2==0) ? Ny*Ly_factor : (Ny-1)*Ly_factor+1);
				Nz = ((Nz%2==0) ? Nz*Lz_factor : (Nz-1)*Lz_factor+1);
				Mx = ((Mx%2==0) ? Mx*Lx_factor : (Mx-1)*Lx_factor+1);
				My = ((My%2==0) ? My*Ly_factor : (My-1)*Ly_factor+1);
				Mz = ((Mz%2==0) ? Mz*Lz_factor : (Mz-1)*Lz_factor+1);
				if (isInitialFieldFromASCIIFile || isInitialFieldFromBinaryFile) {
					Nx_file = ((Nx_file%2==0) ? Nx_file*Lx_factor : (Nx_file-1)*Lx_factor+1);
					Ny_file = ((Ny_file%2==0) ? Ny_file*Ly_factor : (Ny_file-1)*Ly_factor+1);
					Nz_file = ((Nz_file%2==0) ? Nz_file*Lz_factor : (Nz_file-1)*Lz_factor+1);
				}
			} else { //In case input was cubic, e.g. Lx=Ly=Lz=2L
				isNotCubic = false;
			}
		} else {
			isNotCubic = false;
			// Must be initialized to its default value for hit.Real3DimField_to_BinaryFile(...)
			Lx_factor = 1;
			Ly_factor = 1;
			Lz_factor = 1;
		}
		if ((PM["omega_x"].GetIsSet() && (fabs(omega_x) > 0.0)) ||
			(PM["omega_y"].GetIsSet() && (fabs(omega_y) > 0.0)) ||
			(PM["omega_z"].GetIsSet() && (fabs(omega_z) > 0.0))) {
			isRotating = true;
			// In case only significant rotation velocity component has been set
			omega_x = ((PM["omega_x"].GetIsSet()) ? omega_x : 0.0);
			omega_y = ((PM["omega_y"].GetIsSet()) ? omega_y : 0.0);
			omega_z = ((PM["omega_z"].GetIsSet()) ? omega_z : 0.0);
		} else {
			isRotating = false;
		}
		if (PM["C_Smag"].GetIsSet() && (C_Smag > 0)) {
			isLES = true;
		} else {
			isLES = false;
		}
		if (PM["SelfAdaptiveTimestep"].GetIsSet() && (SelfAdaptiveTimestep == 1)) {
			isSelfAdaptive = true;
		} else {
			isSelfAdaptive = false;
		}
		if (PM["ComplexConjugateCorrection"].GetIsSet() && (ComplexConjugateCorrection == 1)) {
			isComplexConjugateCorrected = true;
		} else {
			isComplexConjugateCorrected = false;
		}
		if (PM["Velocity_Physical"].GetIsSet() && (Velocity_Physical == 1)) {
			isVelPhysOut = true;
			isVelPhysOut_iter = ((PM["Velocity_Phys_Freq"].GetIsSet() && (Velocity_Phys_Freq > 0)) ? true : false);
			if (!PM["Velocity_Phys_Folder"].GetIsSet()) snprintf(VelPhysfolder, 256, ".");
		} else {
			isVelPhysOut = false;
			isVelPhysOut_iter = false;
		}
		if (PM["Velocity_Fourier"].GetIsSet() && (Velocity_Fourier == 1)) {
			isVelFourOut = true;
			isVelFourOut_iter = ((PM["Velocity_Four_Freq"].GetIsSet() && (Velocity_Four_Freq > 0)) ? true : false);
			if (!PM["Velocity_Four_Folder"].GetIsSet()) snprintf(VelFourfolder, 256, ".");
		} else {
			isVelFourOut = false;
			isVelFourOut_iter = false;
		}
		if (PM["Energy_Cascade"].GetIsSet() && (Energy_Cascade == 1)) {
			isEkOut = true;
			isEkOut_iter = ((PM["Energy_Cascade_Freq"].GetIsSet() && (Energy_Cascade_Freq > 0)) ? true : false);
			if (!PM["Energy_Cascade_Folder"].GetIsSet()) snprintf(Ekfolder , 256, ".");
		} else {
			isEkOut = false;
			isEkOut_iter = false;
		}
	}

	//===================================================
	// INITIALIZATION OF OBJECTS
	//===================================================
	//Homogeneous Isotropic turbulence solver
	HIT hit(Nx,Ny,Nz,Mx,My,Mz,nu,C_Smag,C_At,Lx_factor,Ly_factor,Lz_factor,omega_x,omega_y,omega_z,ASCII_Input_Filename,Binary_Input_Filename,Nx_file,Ny_file,Nz_file,isReLambda,isInitialFieldFromBinaryFile,isInitialFieldFromASCIIFile,isNotCubic,isRotating,isLES,isComplexConjugateCorrected,isSelfAdaptive);

	//Energy cascade to be forced (if so)
	REAL *Forced_Ek = NULL;
	if (isForcedEk) {
		//Memory allocation
		Forced_Ek = alloc_real(last_input_rad+1);
		//Reading input file
		Read_Ek_Input_File (last_input_rad, Forced_Ek, 1, EkFilename);
	}

	//If Reynolds lambda is passed, then nu has to be set accordingly
	if (isReLambda) {
		REAL Forced_Ek_Tot, Forced_pseudoEpsilon_Tot;
		if (isForcedEk) {
			Forced_Ek_Tot = 0.0;
			Forced_pseudoEpsilon_Tot = 0.0;
			for (int K=1; K<=last_input_rad; K++) { //Forced_Ek[0] assumed to be 0.0
				Forced_Ek_Tot += Forced_Ek[K];
				Forced_pseudoEpsilon_Tot += (K * K * Forced_Ek[K]);
			}
		} else {
			Forced_Ek_Tot = hit.getEk_init_file();
			Forced_pseudoEpsilon_Tot = hit.getpseudoEpsilon_init_file();
		}
		nu = Forced_Ek_Tot / ReLambda * sqrt(10.0 / 3.0 / Forced_pseudoEpsilon_Tot);
		hit.setnu(nu);
		//As it depends on nu, the initial time-step was undefined and needs to be recomputed
		hit.Recalculate_TimeStep();
	}

	if (myrank == 0){
		printf("\n################################################\n");
		printf("### BRUTAL FFTW SPECTRAL CODE (BFSC) 2019-V2 ###\n");
		printf("###------------------------------------------###\n");
		printf("### MPI tasks: %4d processes                ###\n", numprocs);
		if (isOpenMP) {
			printf("###------------------------------------------###\n");
			printf("### OMP threads: %2d threads/process          ###\n", OMP_Threads);
		}
		if (isLES) {
			printf("###------------------------------------------###\n");
			printf("### Smagorinsky constant: %.4f             ###\n", C_Smag);
		}
		if (isReLambda) {
			printf("###------------------------------------------###\n");
			printf("### Kinematic viscosity: %.2e            ###\n", nu);
		}
		if (isRotating) {
			printf("###------------------------------------------###\n");
			printf("### Rotation velocity:                       ###\n");
			printf("### w_x:%.1e,  w_y:%.1e,  w_z:%.1e ###\n", omega_x, omega_y, omega_z);
		}
		printf("###------------------------------------------###\n");
		printf("### Total number of modes:                   ###\n");
		printf("### Nx:%4d,  Ny:%4d,  Nz:%4d              ###\n", Nx, Ny, Nz);
		printf("### Mx:%4d,  My:%4d,  Mz:%4d              ###\n", Mx, My, Mz); //Values used in FFT
		if (isNotCubic) {
			int ActualNx = ((Nx%2==0) ? Nx/Lx_factor : (Nx-1)/Lx_factor+1);
			int ActualNy = ((Ny%2==0) ? Ny/Ly_factor : (Ny-1)/Ly_factor+1);
			int ActualNz = ((Nz%2==0) ? Nz/Lz_factor : (Nz-1)/Lz_factor+1);
			int ActualMx = ((Mx%2==0) ? Mx/Lx_factor : (Mx-1)/Lx_factor+1);
			int ActualMy = ((My%2==0) ? My/Ly_factor : (My-1)/Ly_factor+1);
			int ActualMz = ((Mz%2==0) ? Mz/Lz_factor : (Mz-1)/Lz_factor+1);
			printf("###------------------------------------------###\n");
			printf("### Number of relevant modes:                ###\n");
			printf("### Nx:%4d,  Ny:%4d,  Nz:%4d              ###\n", ActualNx, ActualNy, ActualNz); //Actual relevant values (when the domain is not cubic)
			printf("### Mx:%4d,  My:%4d,  Mz:%4d              ###\n", ActualMx, ActualMy, ActualMz);
		}
		printf("################################################\n\n");
	}

	//===================================================
	// MAIN ITERATIVE LOOP
	//===================================================
	//Auxiliary variables
	int iter = 0; //Number of iteration
	clock_t timebegloop = clock(); //Time before main loop
#if(QA) //USED FOR QA TESTS
	pprintf("\n\n");
#endif

	while (hit.gettime() < Final_Time) {
		//Periodic NS resolution
		hit.New_Fractional_Step_Method();
		//Force energy cascade
		if (isForcedEk && (hit.gettime() < Final_Forced_Time)) {
			hit.Forcing_Energy_Cascade(Forced_Ek, last_input_rad, isNullifyMissingEk);
		}
		//Force Reynolds lambda
		if (isReLambda_iter) {
			if (iter % ReLambda_Freq == 0) {
				if(!myrank) printf("\nnu: %f -> ", hit.getnu());
				hit.Recalculate_Kinematic_Viscosity(ReLambda);
				//Update time-step to new nu
				hit.Recalculate_TimeStep();
				if(!myrank) printf("%f\n", hit.getnu());
			}
		}
		//Update of iteration counter
		iter++;
#if(QA) //USED FOR QA TESTS
		if (iter % 10 == 0) {
			hit.Recalculate_Energy_Cascade();
			pprintf("Time: %f,    iter: %5d,    Ek: %f\n", hit.gettime(), iter, hit.getEk_Tot());
		}
#endif
		if (iter % 100 == 0) {
			hit.Recalculate_Energy_Cascade();
			if (myrank == 0) {
				printf("Time: %f,    iter: %5d,    Ek: %f\n", hit.gettime(), iter, hit.getEk_Tot());
			}
		}
		//Data output (Energy cascade and Fourier velocity use iter+1 instead of iter because new real velocity is not computed
		//until New_Fractional_Step_Method() is called (IFFT of Fourier velocity done in Recalculate_R_vector_Fourier_Coefficients())
		//Thus, in order that all output fields correspond to the same iteration, this modification has to be made)
		if (isEkOut_iter) {
			if ((iter+1) % Energy_Cascade_Freq == 0) {
				snprintf(Ekfilename, 512, "%s/%s_%d.%s", Ekfolder, base_Ekfilename.c_str(), Ek_file_iter, ASCII_fileformat.c_str());
				hit.Recalculate_Energy_Cascade(Ekfilename); //must be called from all ranks
				Ek_file_iter++;
			}
		}
		if (isVelFourOut_iter) {
			if ((iter+1) % Velocity_Four_Freq == 0) {
				snprintf(VelFourfilename, 512, "%s/%s_%d.%s", VelFourfolder, base_VelFourfilename.c_str(), VelFour_file_iter, Binary_fileformat.c_str());
				hit.FourierVelocity_to_BinaryFile(VelFourfilename); //must be called from all ranks
				VelFour_file_iter++;
			}
		}
		if (isVelPhysOut_iter) {
			if (iter % Velocity_Phys_Freq == 0) {
				snprintf(VelPhysfilename, 512, "%s/%s_%d.%s", VelPhysfolder, base_VelPhysfilename.c_str(), VelPhys_file_iter, Binary_fileformat.c_str());
				hit.RealVelocity_to_BinaryFile(VelPhysfilename); //must be called from all ranks
				VelPhys_file_iter++;
			}
		}
	}
	clock_t timeend = clock(); //Time after main loop

	//===================================================
	// FINAL OUTPUT OF RESULTS
	//===================================================
	//Output of energy cascade (ASCII)
	if (isEkOut) {
		snprintf(Ekfilename, 512, "%s/%s.%s", Ekfolder, base_Ekfilename.c_str(), ASCII_fileformat.c_str());
		hit.Recalculate_Energy_Cascade(Ekfilename);
	}
	//Output of velocity Fourier coefficients (binary)
	if (isVelFourOut) {
		snprintf(VelFourfilename, 512, "%s/%s.%s", VelFourfolder, base_VelFourfilename.c_str(), Binary_fileformat.c_str());
		hit.FourierVelocity_to_BinaryFile(VelFourfilename);
	}
	//Output of velocity field (binary)
	if (isVelPhysOut) {
		hit.Recalculate_R_vector_Fourier_Coefficients(); //In order to compute IFFT of velocity (otherwise Fourier velocity-iteration n, and real veocity-iteration n-1)
		snprintf(VelPhysfilename, 512, "%s/%s.%s", VelPhysfolder, base_VelPhysfilename.c_str(), Binary_fileformat.c_str());
		hit.RealVelocity_to_BinaryFile(VelPhysfilename);
	}
	//Terminal print of execution times
	hit.Recalculate_Energy_Cascade();
	double local_timeTot = ((double) (timeend - timebeg)) / CLOCKS_PER_SEC;
	double local_timeLoop = ((double) (timeend - timebegloop)) / CLOCKS_PER_SEC;
	double timeTot = 0.0;
	double timeLoop = 0.0;
	MPI_Allreduce(&local_timeTot, &timeTot, 1, REAL_MPI, MPI_MAX, MCW);
	MPI_Allreduce(&local_timeLoop, &timeLoop, 1, REAL_MPI, MPI_MAX, MCW);
	if (myrank == 0) {
		printf("\nFinal time: %f\t Final kinetic energy: %f\n", hit.gettime(), hit.getEk_Tot());
		printf("Number of iterations: %d\t Total execution time: %f\t Main loop time: %f\n", iter, timeTot, timeLoop);
		printf("TIME/ITERATION: %f\n\n", (timeLoop/iter));
	}
#if(QA) //USED FOR QA TESTS
	pprintf("\nFinal time: %f\t Final kinetic energy: %f\n", hit.gettime(), hit.getEk_Tot());
	pprintf("Number of iterations: %d\t Total execution time: %f\t Main loop time: %f\n", iter, timeTot, timeLoop);
	pprintf("TIME/ITERATION: %f\n\n", (timeLoop/iter));
#endif

	//Free memory
	FREE(Forced_Ek);

	//Ending of parallelization
	end_parallel();

	return 0;
}


//=====================================================================================================================
// Auxiliary functions
//=====================================================================================================================
//Read and load of energy cascade from input file of known size
void Read_Ek_Input_File (const int& Last_K, REAL* Ek_input, const int& dumb_columns, const char* EkFilename) {
	const int myrank = MyID();

	if (!myrank) {
		ifstream infile(EkFilename);
		if(infile.fail()) {
			crash("Read_Ek_Input_File(): Forced energy cascade input file is missing!!! Include \"%s\", please\n", EkFilename);
		}

		int row = 0;
		int col = 0;
		REAL shit;
		Ek_input[0] = 0.0;
		while((!infile.eof()) && (row < Last_K)) {
			if (col < dumb_columns) {
				infile >> shit;
				if (static_cast<int>(shit) != row+1) crash("Read_Ek_Input_File(): invalid entry in (row,col)=(%d,%d)\n", row, col);
				col++;
			} else {
				infile >> Ek_input[row+1]; //Ek for |k|=i is stored in Ek_input[i] and Ek[0] is unused (and must not be set in input file)
				row++;
				col=0;
			}
		}
		infile.close();
	}

	//Broadcast energy cascade to all ranks
	MPI_Bcast(Ek_input, Last_K+1, MPI_REAL, 0, MPI_COMM_WORLD);
};
