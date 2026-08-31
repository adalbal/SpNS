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
//Field is stored in parallel in a shared binary file - collective routine!
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
	MPI_Offset bytes_per_component = sizeof(REAL) * ActualMx * ActualMy * 2*(ActualMz/2 + 1);
    MPI_Offset offset_Real_test = sizeof(REAL) * local_0_start * ActualMy * 2*(ActualMz/2 + 1);
    MPI_Offset num_Real_test = local_n0 * ActualMy * 2*(ActualMz/2 + 1);
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
		MPI_File_seek(fh, offset_Real_test+coord*bytes_per_component, MPI_SEEK_SET);
		MPI_File_write(fh, field, num_Real_test, REAL_MPI, &status);
	}
	//Closing MPI file
	MPI_File_close(&fh);
};

void HIT::BinaryFile_to_Real3DimField (REAL* field_x, REAL* field_y, REAL* field_z) {
    unsigned int total_file_num = ActualNx_file * ActualNy_file * 2 * (ActualNz_file/2 + 1);
    //if (myrank == 0) std::cout << total_file_num << std::endl;
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
    //Read of input file (only by rank 0)
    MPI_Offset offset_Real_x  = sizeof(REAL) * local_0_start_file * ActualNy_file * 2 * (ActualNz_file/2 + 1);
    MPI_Offset offset_Real_y  = offset_Real_x + sizeof(REAL) * total_file_num;
    MPI_Offset offset_Real_z  = offset_Real_y + sizeof(REAL) * total_file_num;

    int num_to_read = local_n0_file * ActualNy_file * 2 * (ActualNz_file/2 + 1);
    //Binary read of input velocity field
    if (filesize != expected_filesize) { //calculation precision != input precision
        switch (BYTES) {
            case 4: //single-precision simulation but double-precision input
                {
                    double *x_input = NULL, *y_input = NULL, *z_input = NULL;
                    x_input = fftw_alloc_real(num_to_read);
                    y_input = fftw_alloc_real(num_to_read);
                    z_input = fftw_alloc_real(num_to_read);
                    MPI_File_seek(fh, offset_Real_x * 2, MPI_SEEK_SET);
                    MPI_File_read(fh, x_input, num_to_read, MPI_DOUBLE, &status);
                    MPI_File_seek(fh, offset_Real_y * 2, MPI_SEEK_SET);
                    MPI_File_read(fh, y_input, num_to_read, MPI_DOUBLE, &status);
                    MPI_File_seek(fh, offset_Real_z * 2, MPI_SEEK_SET);
                    MPI_File_read(fh, z_input, num_to_read, MPI_DOUBLE, &status);
                    for (int ind=0; ind<num_to_read; ind++) {
                        field_x[ind] = x_input[ind];
                        field_y[ind] = y_input[ind];
                        field_z[ind] = z_input[ind];
                    }
                    fftw_free(x_input); fftw_free(y_input); fftw_free(z_input);
                    break;
                }
            case 8: //double-precision simulation but single-precision input
                {

                    float *x_input = NULL, *y_input = NULL, *z_input = NULL;
                    x_input = fftwf_alloc_real(num_to_read);
                    y_input = fftwf_alloc_real(num_to_read);
                    z_input = fftwf_alloc_real(num_to_read);
                    MPI_File_seek(fh, offset_Real_x / 2, MPI_SEEK_SET);
                    MPI_File_read(fh, x_input, num_to_read, MPI_FLOAT, &status);
                    MPI_File_seek(fh, offset_Real_y / 2, MPI_SEEK_SET);
                    MPI_File_read(fh, y_input, num_to_read, MPI_FLOAT, &status);
                    MPI_File_seek(fh, offset_Real_z / 2, MPI_SEEK_SET);
                    MPI_File_read(fh, z_input, num_to_read, MPI_FLOAT, &status);
                    for (int ind=0; ind<num_to_read; ind++) {
                        field_x[ind] = x_input[ind];
                        field_y[ind] = y_input[ind];
                        field_z[ind] = z_input[ind];
                    }
                    fftwf_free(x_input); fftwf_free(y_input); fftwf_free(z_input);
                    break;

                }
        }
    } else { //calculation precision == input precision
        MPI_File_seek(fh, offset_Real_x, MPI_SEEK_SET);
        MPI_File_read(fh, field_x, num_to_read, REAL_MPI, &status);
        MPI_File_seek(fh, offset_Real_y, MPI_SEEK_SET);
        MPI_File_read(fh, field_y, num_to_read, REAL_MPI, &status);
        MPI_File_seek(fh, offset_Real_z, MPI_SEEK_SET);
        MPI_File_read(fh, field_z, num_to_read, REAL_MPI, &status);
    }
    MPI_File_close(&fh);
}

void HIT::ASCIIFile_to_Real3DimField (REAL* field_x, REAL* field_y, REAL* field_z) {
    unsigned int total_file_num = ActualNx_file * ActualNy_file * ActualNz_file;
    //File read by only one MPI process (proc 0)
    REAL* x_read;
    REAL* y_read;
    REAL* z_read;
    if (myrank == 0) {
		x_read = alloc_real(total_file_num);
        y_read = alloc_real(total_file_num);
        z_read = alloc_real(total_file_num);
        unsigned int row = 0;
        ifstream infile;
        infile.open(ASCII_Input_Filename);
        if(infile.fail()) {
            crash("Error opening \"%s\"!!\n", ASCII_Input_Filename);
        }
        while((!infile.eof()) && (row < total_file_num)) {
            infile >> x_read[row] >> y_read[row] >> z_read[row];
            row++;
        }
        infile.close();
    }
    //Data distributed to the rest
    ptrdiff_t message_size = ActualNz_file;
    int proc_send = 0;
    MPI_Status status;

    MPI_Request * x_send_requests;
    MPI_Request * y_send_requests;
    MPI_Request * z_send_requests;
    if (myrank == 0) {
        x_send_requests = new MPI_Request[ActualNx_file * ActualNy_file];
        y_send_requests = new MPI_Request[ActualNx_file * ActualNy_file];
        z_send_requests = new MPI_Request[ActualNx_file * ActualNy_file];
    }
    MPI_Request * x_recv_requests = new MPI_Request[local_n0_file_[myrank] * ActualNy_file];
    MPI_Request * y_recv_requests = new MPI_Request[local_n0_file_[myrank] * ActualNy_file];
    MPI_Request * z_recv_requests = new MPI_Request[local_n0_file_[myrank] * ActualNy_file];
    int num_request_send = 0;
    int num_request_recv = 0;
    //Process 0 executes send requests
    if (myrank == 0) {
        REAL* x_read_send = x_read;
        REAL* y_read_send = y_read;
        REAL* z_read_send = z_read;
        for (int i = 0; i < ActualNx_file; i++) {
            int proc_recv = -1;
            int proc_max_freq = 0;
            for(int k = 0; k < numprocs; k++){
                proc_max_freq += local_n0_file_[k];
                if(i < proc_max_freq){
                    proc_recv = k;
                    break;
                }
            }
            for (int j = 0; j < ActualNy_file; j++) {
                MPI_Isend(x_read_send, message_size, REAL_MPI, proc_recv, 0, MCW, &x_send_requests[num_request_send]);
                MPI_Isend(y_read_send, message_size, REAL_MPI, proc_recv, 0, MCW, &y_send_requests[num_request_send]);
                MPI_Isend(z_read_send, message_size, REAL_MPI, proc_recv, 0, MCW, &z_send_requests[num_request_send]);
                num_request_send++;
                x_read_send += ActualNz_file;
                y_read_send += ActualNz_file;
                z_read_send += ActualNz_file;
            }
        }
    }
    //Executing receive requests
    REAL* field_x_recv = field_x;
    REAL* field_y_recv = field_y;
    REAL* field_z_recv = field_z;
    for (int i = 0; i < local_n0_file; i++) {
        for (int j = 0; j < ActualNy_file; j++) {
            MPI_Irecv(field_x_recv, message_size, COMPLEX_MPI, proc_send, 0, MCW, &x_recv_requests[num_request_recv]);
            MPI_Irecv(field_y_recv, message_size, COMPLEX_MPI, proc_send, 0, MCW, &y_recv_requests[num_request_recv]);
            MPI_Irecv(field_z_recv, message_size, COMPLEX_MPI, proc_send, 0, MCW, &z_recv_requests[num_request_recv]);
            num_request_recv++;
            field_x_recv += 2*(ActualNz_file/2+ 1);
            field_y_recv += 2*(ActualNz_file/2+ 1);
            field_z_recv += 2*(ActualNz_file/2+ 1);
        }
    }
    //Wait for comunications to complete
    for (int i = 0; i < num_request_send-1; i++) {
        MPI_Wait(&x_send_requests[i], &status);
        MPI_Wait(&y_send_requests[i], &status);
        MPI_Wait(&z_send_requests[i], &status);
    }
    for (int i = 0; i < num_request_recv-1; i++) {
        MPI_Wait(&x_recv_requests[i], &status);
        MPI_Wait(&y_recv_requests[i], &status);
        MPI_Wait(&z_recv_requests[i], &status);
    }
    MPI_Barrier(MCW);
    //std::cout << myrank << " Sends completed" << std::endl;
    if (myrank == 0) {
        delete [] x_read;
        delete [] y_read;
        delete [] z_read;
        delete [] x_send_requests;
        delete [] y_send_requests;
        delete [] z_send_requests;
    }
    delete [] x_recv_requests;
    delete [] y_recv_requests;
    delete [] z_recv_requests;
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
