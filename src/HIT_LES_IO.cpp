#include "HIT_LES.h"

//=====================================================================================================================
// MPI I/O
//=====================================================================================================================

// Output velocity Fourier coefficients (into a binary file) - collective routine!
void HIT::FourierVelocity_to_BinaryFile(char const * const filename) const {
	DealiasedComplex3DimField_to_BinaryFile(uk_1, vk_1, wk_1, filename);
};
// Output velocity field (into a binary file) - collective routine!
void HIT::RealVelocity_to_BinaryFile(char const * const filename) const {
	Real3DimField_to_BinaryFile(u, v, w, filename);
};

//It only outputs significant dealiased terms, i.e.,
//	-> k1<=(Nx_2+1) instead of (Mx_2+1) (idem for y and z-axis)
//	-> k1%Lx_factor==0 (idem for y and z-axis)
//Field is stored in parallel in a shared binary file - collective routine!
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
		//Not called from inside any parallel region -- this relies on serial execution:
		//ind_dealiased is a running compaction index, not independent per iteration.
		LOOP_FOURIER {
			const int ind = (b*Mx + a) * (Mz_2 + 1) + k3;
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
//Field is stored in parallel in a shared binary file - collective routine!
void HIT::Real3DimField_to_BinaryFile(REAL const * const field_x, REAL const * const field_y, REAL const * const field_z, char const * const filename) const {
    //Objective: get a 3D slice from a local_n0*My*(Mz+pad) array of size local_n0*My*Mz to remove padding
    MPI_Datatype slice_1D, slice_1D_pad, slice_2D, slice_3D;
    int padding_z = 2*(Mz/2 + 1) - Mz;
    //First slice in z-axis
    MPI_Type_vector(Mz, 1, 1, REAL_MPI, &slice_1D);
    //Second slice in y-axis
    MPI_Type_create_resized(slice_1D, 0, (Mz + padding_z) * BYTES, &slice_1D_pad);
    MPI_Type_contiguous(My, slice_1D_pad, &slice_2D);
    //Third slice in x-axis
    MPI_Type_contiguous(local_n0, slice_2D, &slice_3D);
    MPI_Type_commit(&slice_3D);
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
	MPI_Offset components = Mx * My * Mz;
    MPI_Offset offset = local_0_start * My * Mz;
	//Parallel output of each component of the real field
    MPI_File_set_view(fh, 0, REAL_MPI, REAL_MPI, "native", MPI_INFO_NULL);
	for (int coord=0; coord<3; coord++) {
		//Temporary alias
		REAL const * field = NULL;
		//Alias assignation
		switch(coord){
			case 0: field = field_x; break;
			case 1: field = field_y; break;
			case 2: field = field_z; break;
		}
		MPI_File_seek(fh, offset + components*coord, MPI_SEEK_SET);
		MPI_File_write(fh, field, 1, slice_3D, &status);
	}
	//Closing MPI file
	MPI_File_close(&fh);
    MPI_Type_free(&slice_3D);
};

void HIT::BinaryFile_to_Real3DimField (REAL * const field_x, REAL * const field_y, REAL * const field_z) {
    unsigned int total_file_num = Nx_file * Ny_file * Nz_file;
    unsigned int expected_filesize        = sizeof(REAL)   * 3 * total_file_num;
    unsigned int expected_single_filesize = sizeof(float)  * 3 * total_file_num;
    unsigned int expected_double_filesize = sizeof(double) * 3 * total_file_num;
    MPI_File fh;
    MPI_Status status;
    MPI_File_open(MCW, Binary_Input_Filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    //Check the size of the input file
    MPI_Offset filesize;
    MPI_File_get_size(fh, &filesize);
    if ((filesize != expected_single_filesize) && (filesize != expected_double_filesize)) crash("Incorrect size of \"%s\"!! %lld bytes instead of %d bytes (for single-precision) or %d bytes (for double-precision).\n", Binary_Input_Filename, filesize, expected_single_filesize, expected_double_filesize);
    MPI_Offset offset_Real_x  = local_0_start_file * Ny_file * Nz_file;
    MPI_Offset offset_Real_y  = offset_Real_x + total_file_num;
    MPI_Offset offset_Real_z  = offset_Real_y + total_file_num;
    int padding_z = 2*(Nz_file/2 + 1) - Nz_file;

    int num_to_allocate = local_n0_file * Ny_file * (Nz_file + padding_z);
    //Binary read of input velocity field
    MPI_Datatype Filetype;
    if (filesize != expected_filesize) { //calculation precision != input precision
        switch (BYTES) {
            case 4: //single-precision simulation but double-precision input
                {
                    Filetype = MPI_DOUBLE;
                    //Create padded memory buffer
                    MPI_Datatype slice_1D, slice_2D, slice_3D;

                    //First slice in z-axis
                    MPI_Type_vector(Nz_file, 1, 1, Filetype, &slice_1D);
                    //Second slice in y-axis
                    MPI_Type_create_hvector(Ny_file, 1, (Nz_file + padding_z) * 8, slice_1D, &slice_2D);
                    //Third slice in x-axis
                    MPI_Type_create_hvector(local_n0_file, 1, Ny_file * (Nz_file + padding_z) * 8, slice_2D, &slice_3D);
                    MPI_Type_commit(&slice_3D);
                    
                    double *x_input = NULL, *y_input = NULL, *z_input = NULL;
                    x_input = fftw_alloc_real(num_to_allocate);
                    y_input = fftw_alloc_real(num_to_allocate);
                    z_input = fftw_alloc_real(num_to_allocate);
                    MPI_File_set_view(fh, 0, Filetype, Filetype, "native", MPI_INFO_NULL);
                    MPI_File_seek(fh, offset_Real_x, MPI_SEEK_SET);
                    MPI_File_seek(fh, offset_Real_y, MPI_SEEK_SET);
                    MPI_File_read(fh, y_input, 1, slice_3D, &status);
                    MPI_File_seek(fh, offset_Real_z, MPI_SEEK_SET);
                    MPI_File_read(fh, z_input, 1, slice_3D, &status);
                    for (int ind=0; ind<num_to_allocate; ind++) {
                        field_x[ind] = x_input[ind];
                        field_y[ind] = y_input[ind];
                        field_z[ind] = z_input[ind];
                    }
                    fftw_free(x_input); fftw_free(y_input); fftw_free(z_input);
                    MPI_Type_free(&slice_3D);
                    break;
                }
            case 8: //double-precision simulation but single-precision input
                {
                    Filetype = MPI_FLOAT;
                    //Create padded memory buffer
                    MPI_Datatype slice_1D, slice_2D, slice_3D;
                    //First slice in z-axis
                    MPI_Type_vector(Nz_file, 1, 1, Filetype, &slice_1D);
                    //Second slice in y-axis
                    MPI_Type_create_hvector(Ny_file, 1, (Nz_file + padding_z) * 4, slice_1D, &slice_2D);
                    //Third slice in x-axis
                    MPI_Type_create_hvector(local_n0_file, 1, Ny_file * (Nz_file + padding_z) * 4, slice_2D, &slice_3D);
                    MPI_Type_commit(&slice_3D);
                    
                    double *x_input = NULL, *y_input = NULL, *z_input = NULL;
                    x_input = fftw_alloc_real(num_to_allocate);
                    y_input = fftw_alloc_real(num_to_allocate);
                    z_input = fftw_alloc_real(num_to_allocate);
                    MPI_File_set_view(fh, 0, Filetype, Filetype, "native", MPI_INFO_NULL);
                    MPI_File_seek(fh, offset_Real_x, MPI_SEEK_SET);
                    MPI_File_read(fh, x_input, 1, slice_3D, &status);
                    MPI_File_seek(fh, offset_Real_y, MPI_SEEK_SET);
                    MPI_File_read(fh, y_input, 1, slice_3D, &status);
                    MPI_File_seek(fh, offset_Real_z, MPI_SEEK_SET);
                    MPI_File_read(fh, z_input, 1, slice_3D, &status);
                    for (int ind=0; ind<num_to_allocate; ind++) {
                        field_x[ind] = x_input[ind];
                        field_y[ind] = y_input[ind];
                        field_z[ind] = z_input[ind];
                    }
                    fftw_free(x_input); fftw_free(y_input); fftw_free(z_input);
                    MPI_Type_free(&slice_3D);
                    break;

                }
        }
    } else { //calculation precision == input precision
        Filetype = REAL_MPI;
        //Create padded memory buffer
        MPI_Datatype slice_1D, slice_2D, slice_3D;
        //First slice in z-axis
        MPI_Type_vector(Nz_file, 1, 1, Filetype, &slice_1D);
        //Second slice in y-axis
        MPI_Type_create_hvector(Ny_file, 1, (Nz_file + padding_z) * BYTES, slice_1D, &slice_2D);
        //Third slice in x-axis
        MPI_Type_create_hvector(local_n0_file, 1, Ny_file * (Nz_file + padding_z) * BYTES, slice_2D, &slice_3D);
        MPI_Type_commit(&slice_3D);

        MPI_File_set_view(fh, 0, Filetype, Filetype, "native", MPI_INFO_NULL);
        MPI_File_seek(fh, offset_Real_x, MPI_SEEK_SET);
        MPI_File_read(fh, field_x, 1, slice_3D, &status);
        MPI_File_seek(fh, offset_Real_y, MPI_SEEK_SET);
        MPI_File_read(fh, field_y, 1, slice_3D, &status);
        MPI_File_seek(fh, offset_Real_z, MPI_SEEK_SET);
        MPI_File_read(fh, field_z, 1, slice_3D, &status);
        MPI_Type_free(&slice_3D);
    }
    MPI_File_close(&fh);
}

void HIT::ASCIIFile_to_Real3DimField (REAL * const field_x, REAL * const field_y, REAL * const field_z) {
    unsigned int total_file_read = local_n0_file * ActualNy_file * ActualNz_file;
    //File read by only one MPI process (proc 0)
    unsigned int row_offset = local_0_start_file * ActualNy_file * ActualNz_file;
    int padding_z = 2 * (ActualNz_file/2 + 1) - ActualNz_file;

    REAL shit;
    ifstream infile;
    infile.open(ASCII_Input_Filename);
    if(infile.fail()) {
        crash("Error opening \"%s\"!!\n", ASCII_Input_Filename);
    }
    unsigned int row = 0;
    while((!infile.eof()) && (row < row_offset)) {
        infile >> shit >> shit >> shit;
        row++;
    }
    unsigned int index = 0;
    unsigned int num_read = 0;
    while((!infile.eof()) && (num_read < total_file_read)) {
        if(num_read%ActualNz_file == 0 && num_read != 0){
            for(int i=0;i<padding_z;++i){
                field_x[index] = 0;
                field_y[index] = 0;
                field_z[index] = 0;
                index++;
            }
            infile >> field_x[index] >> field_y[index] >> field_z[index];
            num_read++;
            index++;
        }
        else{
            infile >> field_x[index] >> field_y[index] >> field_z[index];
            num_read++;
            index++;
        }
    }
    infile.close();
}

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
