#include "HIT_LES.h"

//=====================================================================================================================
// MPI I/O
//=====================================================================================================================

// Output velocity Fourier coefficients (into a binary file)
void HIT::FourierVelocity_to_BinaryFile(char const * const filename) const {
	DealiasedComplex3DimField_to_BinaryFile(uk_1, vk_1, wk_1, filename);
};
// Output velocity field (into a binary file)
void HIT::RealVelocity_to_BinaryFile(char const * const filename) const {
	Real3DimField_to_BinaryFile(u, v, w, filename);
};

//It only outputs significant dealiased terms, i.e.,
//	-> k1<=(Nx_2+1) instead of (Mx_2+1) (idem for y and z-axis)
//	-> k1%Lx_factor==0 (idem for y and z-axis)
//Field is stored in parallel in a shared binary file
void HIT::DealiasedComplex3DimField_to_BinaryFile(COMPLEX const * const field_x, COMPLEX const * const field_y, COMPLEX const * const field_z, char const * const filename) const {
	//Opening MPI file
	MPI_File fh;
	MPI_Status status;
	int err = MPI_File_open(MCW, filename, MPI_MODE_CREATE | MPI_MODE_EXCL | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
	if (err) {
		MPI_File_close(&fh);
		MPI_File_open(MCW, filename, MPI_MODE_DELETE_ON_CLOSE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
		MPI_File_close(&fh);
		MPI_File_open(MCW, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
	}
	MPI_Offset bytes_per_component = sizeof(COMPLEX) * ActualNx * ActualNy * (ActualNz/2+1);
	//Allocate memory for current's process dealiased field[ind] coefficients
	COMPLEX *field_dealiased = alloc_complex(num_dealiased);
	//Parallel output of each component of the complex field
	for(int coord=0; coord<3; coord++){
		//Temporary alias
		COMPLEX const * field = NULL;
		//Alias assignation
		switch(coord){
			case 0: field = field_x; break;
			case 1: field = field_y; break;
			case 2: field = field_z; break;
		}
		//Initialize field_dealiased
		int ind_dealiased = 0; //Note that ind_dealiased<=(num_dealiased-1)
		LOOP_FOURIER {
			if (dealiased[ind]) {
				for (int ic=0; ic<=1; ic++) { //ic=0 => Real part, ic=1 => Imaginary part
					field_dealiased[ind_dealiased][ic] = field[ind][ic];
				}
				ind_dealiased++;
			}
		}
		if (ind_dealiased!=num_dealiased) crash("The number of dealiased complex terms does not match!\n");
		//Output of field_dealiased
		MPI_File_seek(fh, offset_Fourier+coord*bytes_per_component, MPI_SEEK_SET);
		MPI_File_write(fh, field_dealiased, num_dealiased*2, REAL_MPI, &status); //1 COMPLEX = 2 doubles
	}
	//Memory freeing
	FREE(field_dealiased);
	//Closing MPI file
	MPI_File_close(&fh);
};
//It outputs all the terms of a real field (the consideration of dealiased terms does not make sense here)
//Field is stored in parallel in a shared binary file
void HIT::Real3DimField_to_BinaryFile(REAL const * const field_x, REAL const * const field_y, REAL const * const field_z, char const * const filename) const {
	//Opening MPI file
	MPI_File fh;
	MPI_Status status;
	int err = MPI_File_open(MCW, filename, MPI_MODE_CREATE | MPI_MODE_EXCL | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
	if (err) {
		MPI_File_close(&fh);
		MPI_File_open(MCW, filename, MPI_MODE_DELETE_ON_CLOSE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
		MPI_File_close(&fh);
		MPI_File_open(MCW, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
	}
	//Output of 3D real field
	MPI_Offset bytes_per_component = sizeof(REAL) * ActualMx * ActualMy * ActualMz;
	//Allocate memory for current's process dealiased field[ind] coefficients
	REAL *field_to_output = alloc_real(num_dealiased_real);
	//Parallel output of each component of the real field
	for (int coord=0; coord<3; coord++) {
		//Temporary alias
		REAL const * field = NULL;
		//Alias assignation
		switch(coord){
			case 0: field = field_x; break;
			case 1: field = field_y; break;
			case 2: field = field_z; break;
		}
		//Initialize field_to_output
		int ind_dealiased = 0; //Note that ind_dealiased<=(num_dealiased_real-1)
		if (isNotCubic) { //Removing fftw padding and irrelevant terms with respect to the length factors
			LOOP_REAL {
				if ((local_0_start+i)<ActualMx && j<ActualMy && k<ActualMz) {
					field_to_output[ind_dealiased] = field[ind];
					ind_dealiased++;
				}
			}
		} else { //Removing fftw padding
			LOOP_REAL {
				field_to_output[ind_dealiased] = field[ind];
				ind_dealiased++;
			}
		}
		if (ind_dealiased!=num_dealiased_real) crash("The number of real terms does not match!\n");
		//Output of field_to_output
		MPI_File_seek(fh, offset_Real+coord*bytes_per_component, MPI_SEEK_SET);
		MPI_File_write(fh, field_to_output, num_dealiased_real, REAL_MPI, &status);
	}
	//Memory freeing
	FREE(field_to_output);
	//Closing MPI file
	MPI_File_close(&fh);
};

//=====================================================================================================================
// DEBUG I/O FUNCTIONS
//=====================================================================================================================

//Output of a real field into an ASCII file by rank 0 (not efficient)
void HIT::Field2File (REAL const * const field, char const * const filename) const {
	//Allocate memory for auxiliar pointers
	int *displs_RealField = (int*) MALLOC(sizeof(int) * numprocs); //Position of the first byte gathered
	int *recvcounts_RealField = (int*) MALLOC(sizeof(int) * numprocs); //Number of MPI_BYTEs (for k3=0)
	//Initialization of recvcounts_RealField and displs_RealField (used in MPI_Allgatherv clause)
	//MPI_Gather is necessary because local_n1 differs from different processes
	if(myrank==0){
		for(int proc=0; proc<numprocs; proc++) {
			recvcounts_RealField[proc] = sizeof(REAL) * local_n0_[proc] * My * 2*(Mz_2+1); //as MPI_BYTEs are sent in block, extra padding terms (from Mz to 2*(Mz_2+1)) have to be sent as well (they are not contiguous in memory)
			displs_RealField[proc] = ((proc == 0) ? 0 : (displs_RealField[proc-1]+recvcounts_RealField[proc-1]));
		}
	}
	MPI_Bcast(recvcounts_RealField, numprocs, MPI_INT, 0, MCW);
	MPI_Bcast(displs_RealField,  numprocs, MPI_INT, 0, MCW);
	//Gathering of data by process 0
	REAL *field_Gathered=NULL;
	if (myrank == 0) {
		field_Gathered = alloc_real(Mx * My * 2*(Mz_2+1));
	}
	MPI_Gatherv(field, recvcounts_RealField[myrank], MPI_BYTE, field_Gathered, recvcounts_RealField, displs_RealField, MPI_BYTE, 0, MCW);
	//Freeing auxiliary memory
	FREE(displs_RealField);
	FREE(recvcounts_RealField);
	//Data output
	if (myrank == 0) {
		ofstream OutputField;
		OutputField.open(filename);
		LOOP_3D(i,0,Mx, j,0,My, k,0,Mz) {
//			if ((local_0_start+i)<ActualMx && j<ActualMy && k<ActualMz) { //If isNonCubic, avoids dummy repeated velocities (ActualN instead of N)
				int l1=FOURIER_FREQ(i,Mx);
				int l2=FOURIER_FREQ(j,My);
				int l3=FOURIER_FREQ(k,Mz);
				int ind = (i*My+j)*2*(Mz_2+1)+k;
				OutputField << l1 << " " << l2 << " " << l3 << "\t" << i << " " << j << " " << k << "\t" << field_Gathered[ind] << endl;
//			}
		}
		OutputField.close();
		//Memory freeing
		FREE(field_Gathered);
	}
};
//Output of a complex field into an ASCII file by rank 0 (not efficient)
void HIT::Field2File (COMPLEX const * const field, char const * const filename) const {
	//Allocate memory for auxiliar pointers
	int *displs_ComplexField = (int*) MALLOC(sizeof(int) * numprocs); //Position of the first byte gathered
	int *recvcounts_ComplexField = (int*) MALLOC(sizeof(int) * numprocs); //Number of MPI_BYTEs (for k3=0)
	//Initialization of recvcounts_ComplexField and displs_ComplexField (used in MPI_Allgatherv clause)
	//MPI_Gather is necessary because local_n1 differs from different processes
	if(myrank==0){
		for(int proc=0; proc<numprocs; proc++) {
			recvcounts_ComplexField[proc] = sizeof(COMPLEX) * Mx * local_n1_[proc] * (Mz_2+1);
			displs_ComplexField[proc] = ((proc == 0) ? 0 : (displs_ComplexField[proc-1]+recvcounts_ComplexField[proc-1]));
		}
	}
	MPI_Bcast(recvcounts_ComplexField, numprocs, MPI_INT, 0, MCW);
	MPI_Bcast(displs_ComplexField,  numprocs, MPI_INT, 0, MCW);
	//Gathering of data by process 0
	COMPLEX *field_Gathered=NULL;
	if (myrank == 0) {
		field_Gathered = alloc_complex(Mx * My * (Mz_2+1));
	}
	MPI_Gatherv(field, recvcounts_ComplexField[myrank], MPI_BYTE, field_Gathered, recvcounts_ComplexField, displs_ComplexField, MPI_BYTE, 0, MCW);
	//Freeing auxiliary memory
	FREE(displs_ComplexField);
	FREE(recvcounts_ComplexField);
	//Data output by rank 0. dealiased[ind] cannot be used as it is defined locally at each processor
	if (myrank == 0) {
		ofstream OutputField;
		OutputField.open(filename);
		LOOP_3D(b,0,My, a,0,Mx, k3,0,(Mz_2+1)) {
			int ind=(b*Mx+a)*(Mz_2+1)+k3;
			int k1=FOURIER_FREQ(a,Mx);
			int k2=FOURIER_FREQ(b,My);
			if ((abs(k1)<=Nx_2 && abs(k2)<=Ny_2 && k3<=Nz_2) && !(Nx%2==0 && k1==-Nx_2) && !(Ny%2==0 && k2==-Ny_2)
				&& (k1%Lx_factor == 0) && (k2%Ly_factor == 0) && (k3%Lz_factor == 0)) { //Avoids aliased terms
				OutputField << k1 << " " << k2 << " " << k3 << "\t" << field_Gathered[ind][0] << " " << field_Gathered[ind][1] << endl;
			}
		}
		OutputField.close();
		//Memory freeing
		FREE(field_Gathered);
	}
};
