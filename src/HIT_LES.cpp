#include "HIT_LES.h"

int iter = 0;

//=====================================================================================================================
// FRACTIONAL STEP METHOD
//=====================================================================================================================

void HIT::New_Fractional_Step_Method() {
	//Calculation of R vector
	Recalculate_R_vector_Fourier_Coefficients();
	//Calculation of the predictor velocity at the following timestep (from R vector)
	Recalculate_Predictor_Velocity_Fourier_Coefficients();
	//Calculation of the new timesteps velocity distribution (from predictor velocity)
	Recalculate_Divergence_Free_Projection();
	//Obtention of new timestep
	if (isSelfAdaptive) {
		Recalculate_SelfAdaptive_TimeStep();
	} else {
		Recalculate_CFL_TimeStep();
	}
	//Update of total simulated time
	time=time+At;
	//Update of intern iteration counter
	iter++;
};

//===================================================
// R VECTOR
//===================================================

//=========================
// 1.) CONVECTIVE TERM IN FOURIER SPACE
// Comment Recalculate_Velocity_Gradients_Fourier_Coefficients 
void HIT::Recalculate_Velocity_Gradients_Fourier_Coefficients() {
	LOOP_FOURIER_k1k2k3{ //Velocity spatial derivatives in Fourier space
		if (dealiased[ind]) {
			gradu_k[0][ind][0] = -k1 * uk_1[ind][1];  gradu_k[0][ind][1] = k1 * uk_1[ind][0]; //dudx_k
			gradu_k[1][ind][0] = -k2 * uk_1[ind][1];  gradu_k[1][ind][1] = k2 * uk_1[ind][0]; //dudy_k
			gradu_k[2][ind][0] = -k3 * uk_1[ind][1];  gradu_k[2][ind][1] = k3 * uk_1[ind][0]; //dudz_k

			gradv_k[0][ind][0] = -k1 * vk_1[ind][1];  gradv_k[0][ind][1] = k1 * vk_1[ind][0]; //dvdx_k
			gradv_k[1][ind][0] = -k2 * vk_1[ind][1];  gradv_k[1][ind][1] = k2 * vk_1[ind][0]; //dvdy_k
			gradv_k[2][ind][0] = -k3 * vk_1[ind][1];  gradv_k[2][ind][1] = k3 * vk_1[ind][0]; //dvdz_k

			gradw_k[0][ind][0] = -k1 * wk_1[ind][1];  gradw_k[0][ind][1] = k1 * wk_1[ind][0]; //dwdx_k
			gradw_k[1][ind][0] = -k2 * wk_1[ind][1];  gradw_k[1][ind][1] = k2 * wk_1[ind][0]; //dwdy_k
			gradw_k[2][ind][0] = -k3 * wk_1[ind][1];  gradw_k[2][ind][1] = k3 * wk_1[ind][0]; //dwdz_k
		} else { //For de-aliasing purposes
			for (int coord=0; coord<3; coord++) {
				for (int ic=0; ic<=1; ic++){ //ic=0 => Real part, ic=1 => Imaginary part
					gradu_k[coord][ind][ic] = 0.0;
					gradv_k[coord][ind][ic] = 0.0;
					gradw_k[coord][ind][ic] = 0.0;
				}
			}
		}
	}
};
//Calculation of physical velocity using Fourier series expansion
void HIT::Recalculate_Velocity_Antitransform() {
	//Auxiliary copy needed because c2r leaves garbage at input
	memcpy(uk_aux, uk_1, sizeof(COMPLEX) * alloc_local);
	memcpy(vk_aux, vk_1, sizeof(COMPLEX) * alloc_local);
	memcpy(wk_aux, wk_1, sizeof(COMPLEX) * alloc_local);
	//IFFT: Vxy
	fftw.IFFTWc2r(uk_aux,u);
	fftw.IFFTWc2r(vk_aux,v);
	fftw.IFFTWc2r(wk_aux,w);
};
//Recalculate_Convective_Fourier_Coefficients():
//	1.) Calculation of gradu, gradv and gradw
//	2.) Calculation convective term in Physical space
//	3.) Calculation convective term in Fourier space
void HIT::Recalculate_Convective_Fourier_Coefficients() {
	//Recalculate gradients of velocity in Fourier space
	Recalculate_Velocity_Gradients_Fourier_Coefficients();
	//Obtention of the velocity field in the Physical space (via the IFFT)
	Recalculate_Velocity_Antitransform();
	//Calculation of the convective term in Fourier space. Loop over [coord] (1=x, 2=y, 3=z)
	for(int coord=0; coord<3; coord++){
		//Temporary aliases
		REAL *convxyz = NULL;
		REAL *graduvw[3] = {NULL,NULL,NULL};
		COMPLEX *convxyz_k = NULL;
		COMPLEX *graduvw_k[3] = {NULL,NULL,NULL};
		//Alias assignation
		switch(coord){
			case 0: 
				for (int sub_coord=0; sub_coord<3; sub_coord++) { //x-velocity gradient
					graduvw[sub_coord] = gradu[sub_coord];
					graduvw_k[sub_coord] = gradu_k[sub_coord];
				}
				convxyz = convx;
				convxyz_k = convx_k;
				break;
			case 1: 
				for (int sub_coord=0; sub_coord<3; sub_coord++) { //y-velocity gradient
					graduvw[sub_coord] = gradv[sub_coord];
					graduvw_k[sub_coord] = gradv_k[sub_coord];
				}
				convxyz = convy;
				convxyz_k = convy_k;
				break;
			case 2: 
				for (int sub_coord=0; sub_coord<3; sub_coord++) { //z-velocity gradient
					graduvw[sub_coord] = gradw[sub_coord];
					graduvw_k[sub_coord] = gradw_k[sub_coord];
				}
				convxyz = convz;
				convxyz_k = convz_k;
				break;
		}
		//IFFT of each component of the gradient of velocity (c2r)
		for(int sub_coord=0; sub_coord<3; sub_coord++){
			fftw.IFFTWc2r(graduvw_k[sub_coord],graduvw[sub_coord]);
		}
		//Calculation of the [coord]-convective term in Physical space
		LOOP_REAL {
			convxyz[ind] = u[ind] * graduvw[0][ind] + v[ind] * graduvw[1][ind] + w[ind] * graduvw[2][ind];
		}
		//Calculation of the [coord]-convective term in Fourier space
		fftw.FFTWr2c(convxyz,convxyz_k); //FFT (r2c)
		fftw.Fourier_Normalization(convxyz_k,dealiased); //Normalization of relevant coefficients
	}
};

//=========================
// 2.) DIVERGENCE OF SGS TENSOR IN FOURIER SPACE (LES MODELLING)
// Calculation of eddy viscosity
void HIT::Recalculate_Eddy_Viscosity() {
	// Smagorinsky model
	const REAL ct = C_Smag * C_Smag * delta * delta;
	LOOP_REAL {
		nu_eddy[ind] = 2.0 * gradu[0][ind] * gradu[0][ind];
		nu_eddy[ind] += 2.0 * gradv[1][ind] * gradv[1][ind];
		nu_eddy[ind] += 2.0 * gradw[2][ind] * gradw[2][ind];
		nu_eddy[ind] += (gradu[1][ind] + gradv[0][ind]) * (gradu[1][ind] + gradv[0][ind]);
		nu_eddy[ind] += (gradu[2][ind] + gradw[0][ind]) * (gradu[2][ind] + gradw[0][ind]);
		nu_eddy[ind] += (gradv[2][ind] + gradw[1][ind]) * (gradv[2][ind] + gradw[1][ind]);
		nu_eddy[ind] = ct * sqrt(nu_eddy[ind]);
	}
};
// Calculation of SGS tensor components in Physical space
void HIT::Recalculate_SGS_Tensor_Components() {
	Recalculate_Eddy_Viscosity();
	LOOP_REAL {
		SGS_11[ind] = (-2.0) * nu_eddy[ind] * gradu[0][ind];
		SGS_21[ind] = (-1.0) * nu_eddy[ind] * (gradu[1][ind] + gradv[0][ind]);
		SGS_22[ind] = (-2.0) * nu_eddy[ind] * gradv[1][ind];
		SGS_31[ind] = (-1.0) * nu_eddy[ind] * (gradu[2][ind] + gradw[0][ind]);
		SGS_32[ind] = (-1.0) * nu_eddy[ind] * (gradv[2][ind] + gradw[1][ind]);
		SGS_33[ind] = (-2.0) * nu_eddy[ind] * gradw[2][ind];
	}
};
// Calculation of SGS tensor components in Fourier space
void HIT::Recalculate_SGS_Tensor_Transform() {
	//FFT: SGS_Tensor
	fftw.FFTWr2c(SGS_11,SGS_11_k);
	fftw.FFTWr2c(SGS_21,SGS_21_k);
	fftw.FFTWr2c(SGS_22,SGS_22_k);
	fftw.FFTWr2c(SGS_31,SGS_31_k);
	fftw.FFTWr2c(SGS_32,SGS_32_k);
	fftw.FFTWr2c(SGS_33,SGS_33_k);
	//Normalization
	fftw.Fourier_Normalization(SGS_11_k,dealiased);
	fftw.Fourier_Normalization(SGS_21_k,dealiased);
	fftw.Fourier_Normalization(SGS_22_k,dealiased);
	fftw.Fourier_Normalization(SGS_31_k,dealiased);
	fftw.Fourier_Normalization(SGS_32_k,dealiased);
	fftw.Fourier_Normalization(SGS_33_k,dealiased);
};
// Calculation of SGS tensor divergence in Fourier space
void HIT::Recalculate_SGS_Tensor_Divergence_Fourier_Coefficients() {
	//SGS tensor components in Physical space
	Recalculate_SGS_Tensor_Components();
	//SGS tensor components in Fourier space
	Recalculate_SGS_Tensor_Transform();
	//SGS tensor divergence components in Fourier space
	LOOP_FOURIER_k1k2k3 {
		if (dealiased[ind]) {
			// divSGSx_k
			divSGSx_k[ind][0] = - k1 * SGS_11_k[ind][1] - k2 * SGS_21_k[ind][1] - k3 * SGS_31_k[ind][1];
			divSGSx_k[ind][1] =   k1 * SGS_11_k[ind][0] + k2 * SGS_21_k[ind][0] + k3 * SGS_31_k[ind][0];
			// divSGSy_k // SGS_12 = SGS_21
			divSGSy_k[ind][0] = - k1 * SGS_21_k[ind][1] - k2 * SGS_22_k[ind][1] - k3 * SGS_32_k[ind][1];
			divSGSy_k[ind][1] =   k1 * SGS_21_k[ind][0] + k2 * SGS_22_k[ind][0] + k3 * SGS_32_k[ind][0];
			// divSGSz_k // SGS_13 = SGS_31 // SGS_23 = SGS_32
			divSGSz_k[ind][0] = - k1 * SGS_31_k[ind][1] - k2 * SGS_32_k[ind][1] - k3 * SGS_33_k[ind][1];
			divSGSz_k[ind][1] =   k1 * SGS_31_k[ind][0] + k2 * SGS_32_k[ind][0] + k3 * SGS_33_k[ind][0];
		}
	}
};

//=========================
// 3.) R VECTOR
void HIT::Recalculate_R_vector_Fourier_Coefficients() {
	//Update last timestep value of R vector
	memcpy(Rx0_k, Rx1_k, sizeof(COMPLEX) * alloc_local);
	memcpy(Ry0_k, Ry1_k, sizeof(COMPLEX) * alloc_local);
	memcpy(Rz0_k, Rz1_k, sizeof(COMPLEX) * alloc_local);
	//Calculate new R vector in Fourier space
	Recalculate_Convective_Fourier_Coefficients();
	REAL ct;
	LOOP_FOURIER {
		if (dealiased[ind]) {
			ct = nu * rad2[ind];
			for (int ic=0; ic<=1; ic++){ //ic=0 => Real part, ic=1 => Imaginary part
				Rx1_k[ind][ic] = - convx_k[ind][ic] - ct * uk_1[ind][ic];
				Ry1_k[ind][ic] = - convy_k[ind][ic] - ct * vk_1[ind][ic];
				Rz1_k[ind][ic] = - convz_k[ind][ic] - ct * wk_1[ind][ic];
			}
		}
	}
	if (isLES) {
		Recalculate_SGS_Tensor_Divergence_Fourier_Coefficients();
		LOOP_FOURIER {
			if (dealiased[ind]) {
				for (int ic=0; ic<=1; ic++){ //ic=0 => Real part, ic=1 => Imaginary part
					Rx1_k[ind][ic] -= divSGSx_k[ind][ic];
					Ry1_k[ind][ic] -= divSGSy_k[ind][ic];
					Rz1_k[ind][ic] -= divSGSz_k[ind][ic];
				}
			}
		}
	}
	if (isRotating) {
		//Considering the system is rotating with angular rotation velocity (omega_x,omega_y,omega_z)
		LOOP_FOURIER {
			if (dealiased[ind]) {
				for (int ic=0; ic<=1; ic++){ //ic=0 => Real part, ic=1 => Imaginary part
					Rx1_k[ind][ic] += 2 * (omega_z * vk_1[ind][ic] - omega_y * wk_1[ind][ic]);
					Ry1_k[ind][ic] += 2 * (omega_x * wk_1[ind][ic] - omega_z * uk_1[ind][ic]);
					Rz1_k[ind][ic] += 2 * (omega_y * uk_1[ind][ic] - omega_x * vk_1[ind][ic]);
				}
			}
		}
	}
};

//===================================================
// PREDICTOR VELOCITY
//===================================================

//Recalculation Predictor Velocity Fourier Coefficients, for a determined (k1,k2,k3)
// INPUT:
//		ukxyz_aux ->	Garbage after c2r fftw_plan
//		ukxyz_0 ->		Previous distribution of velocity (in Fourier space)
//		ukxyz_1 ->		Current distribution of velocity (in Fourier space)
// OUTPUT:
//		ukxyz_aux ->	Previous distribution of velocity (in Fourier space)
//		ukxyz_0 ->		Current distribution of velocity (in Fourier space)
//		ukxyz_1 ->		New Predictor velocity (in Fourier space)
void HIT::Recalculate_Predictor_Velocity_Fourier_Coefficients() {
	if (isSelfAdaptive) {
		//Update old values
		if (iter > 0) { //Initial old velocity field equal to 0.0 unstabilities
			memcpy(uk_aux, uk_0, sizeof(COMPLEX) * alloc_local);
			memcpy(vk_aux, vk_0, sizeof(COMPLEX) * alloc_local);
			memcpy(wk_aux, wk_0, sizeof(COMPLEX) * alloc_local);
		}
		memcpy(uk_0, uk_1, sizeof(COMPLEX) * alloc_local);
		memcpy(vk_0, vk_1, sizeof(COMPLEX) * alloc_local);
		memcpy(wk_0, wk_1, sizeof(COMPLEX) * alloc_local);
		//Related to self-adaptive timestep
		_kappa05 = (1.0/(kappa+0.5));
		//Calculation of new coefficients
		LOOP_FOURIER {
			if (dealiased[ind]) {
				for (int ic=0; ic<=1; ic++){ //ic=0 => Real part, ic=1 => Imaginary part
					//Adaptive timestep
					_kappa05 = (1.0/(kappa+0.5));
					uk_1[ind][ic]=_kappa05*(2.0*kappa*uk_0[ind][ic]+At*((1.0+kappa)*Rx1_k[ind][ic]
											   -kappa*Rx0_k[ind][ic])-(kappa-0.5)*uk_aux[ind][ic]);
					vk_1[ind][ic]=_kappa05*(2.0*kappa*vk_0[ind][ic]+At*((1.0+kappa)*Ry1_k[ind][ic]
											   -kappa*Ry0_k[ind][ic])-(kappa-0.5)*vk_aux[ind][ic]);
					wk_1[ind][ic]=_kappa05*(2.0*kappa*wk_0[ind][ic]+At*((1.0+kappa)*Rz1_k[ind][ic]
											   -kappa*Rz0_k[ind][ic])-(kappa-0.5)*wk_aux[ind][ic]);
				}
			}
		}
	} else {
		//Update old values
		memcpy(uk_0, uk_1, sizeof(COMPLEX) * alloc_local);
		memcpy(vk_0, vk_1, sizeof(COMPLEX) * alloc_local);
		memcpy(wk_0, wk_1, sizeof(COMPLEX) * alloc_local);
		//Calculation of new coefficients
		LOOP_FOURIER {
			if (dealiased[ind]) {
				for (int ic=0; ic<=1; ic++){ //ic=0 => Real part, ic=1 => Imaginary part
					//Adams-Bashforth method
					uk_1[ind][ic] = uk_0[ind][ic] + At * (1.5 * Rx1_k[ind][ic] - 0.5 * Rx0_k[ind][ic]);
					vk_1[ind][ic] = vk_0[ind][ic] + At * (1.5 * Ry1_k[ind][ic] - 0.5 * Ry0_k[ind][ic]);
					wk_1[ind][ic] = wk_0[ind][ic] + At * (1.5 * Rz1_k[ind][ic] - 0.5 * Rz0_k[ind][ic]);
				}
			}
		}
	}
};

//===================================================
// DIVERGENCE-FREE PROJECTION OF PREDICTOR VELOCITY
//===================================================

//Divergence-free projection of predictor velocity (result of substracting the gradient of pseudo-pressure -> operator (k^t·k)/(k·k))
void HIT::Recalculate_Divergence_Free_Projection() {
	REAL aux;
	//Projection of Predictor velocity and calculation of uk_n+1 (stored in ukxy_1 replacing predictor velocity)
	LOOP_FOURIER_k1k2k3 {
		if (dealiased[ind]) {
			for (int ic=0; ic<=1; ic++){ //ic=0 => Real part, ic=1 => Imaginary part
				aux = ((rad2[ind]>0) ? ((k1*uk_1[ind][ic]+k2*vk_1[ind][ic]+k3*wk_1[ind][ic])/rad2[ind]) : 0.0);
				uk_1[ind][ic] -= aux * k1;
				vk_1[ind][ic] -= aux * k2;
				wk_1[ind][ic] -= aux * k3;
			}
		}
	}
};

//===================================================
// CALCULATE NEW TIMESTEP:
//===================================================
//Calculation of new self-adaptive timestep
void HIT::Recalculate_SelfAdaptive_TimeStep() {
	//Compute local maximums at each process
	REAL local_MaxVel = 0.0;
	LOOP_REAL {
		local_MaxVel = max(local_MaxVel, (fabs(u[ind])/Ax + fabs(v[ind])/Ay + fabs(w[ind])/Az));
	}
	//Compute maximum accros processes
	MPI_Allreduce(&local_MaxVel, &MaxVel, 1, REAL_MPI, MPI_MAX, MCW);
	//Recalculate self-adaptive timestep
	phi0 = atan(MaxVel*Ah*Ah/(0.35*nu*12.0)); //0.35 per la constant Cc
	At = Topt(phi0)/sqrt((MaxVel*MaxVel)/(0.35*0.35)+(144.0*nu*nu/(Ah*Ah*Ah*Ah)));
	kappa = Kopt(phi0);
	At = C_At*At;
};

//Calculation of new CFL timestep
void HIT::Recalculate_CFL_TimeStep() {
	//Compute local maximums at each process
	REAL local_MaxVel = 0.0;
	LOOP_REAL {
		local_MaxVel = max(local_MaxVel, max(fabs(u[ind])/Ax, max(fabs(v[ind])/Ay, fabs(w[ind])/Az)));
	}
	//Compute maximum accros processes
	MPI_Allreduce(&local_MaxVel, &MaxVel, 1, REAL_MPI, MPI_MAX, MCW);
	//Recalculate CFL timestep
	C_conv = 0.35 / MaxVel;
	At = min(C_visc, C_conv);
	At = C_At * At;
};

//=====================================================================================================================
// FORCING TERM
//=====================================================================================================================

// Imposition of a given Energy distribution ("normalizing" current Fourier coefficients)
void HIT::Forcing_Energy_Cascade(REAL const * const Forced_Ek, const ptrdiff_t& last_input_rad) {
	//Update Ek in order to calculate forcing term scaling factors FT[k]
	Recalculate_Energy_Cascade();
	//Internal temporary variables
	int last_eff_rad = min(last_input_rad, last_rad);
	REAL *FT = alloc_real(last_eff_rad+1);
	//Calculation of FT corresponding to the current energy cascade (Ek)
	FT[0] = 0.0; //Dummy case
	for (int K=1; K<=last_eff_rad; K++) {
		FT[K] = ((Ek[K]>0.0) ? sqrt(Forced_Ek[K] / Ek[K]) : 0.0); //Beware of Ek[K]=0!!!
	}
	//Scaling of velocity to force the energy cascade to be as Forced_Ek
	LOOP_FOURIER {
		if (dealiased[ind]) {
			int K = static_cast<int>(round(sqrt(rad2[ind])));
			if (K <= last_eff_rad) {
				for (int ic=0; ic<=1; ic++){ //ic=0 => Real part, ic=1 => Imaginary part
					uk_1[ind][ic] *= FT[K];
					vk_1[ind][ic] *= FT[K];
					wk_1[ind][ic] *= FT[K];
				}
			} else { //Extra modes (with respect to Forced_Ek)
				for (int ic=0; ic<=1; ic++){ //ic=0 => Real part, ic=1 => Imaginary part
					uk_1[ind][ic] = 0.0;
					vk_1[ind][ic] = 0.0;
					wk_1[ind][ic] = 0.0;
				}
			}
		}
	}
	FREE(FT);
	if (isComplexConjugateCorrected) {
		//Correction of explicit complex conjugates (due to numerical instabilities)
		Complex_Conjugate_Correction();
	}
};
// Imposition of correct complex conjugation for k3=0
void HIT::Complex_Conjugate_Correction() {
	//Load of velocity Fourier coefficients at k3=0 into local_uvwk_k3_0
	LOOP_2D(b,0,local_n1, a,0,Mx) { //For plane k3=0
		int ind = (b*Mx+a)*(Mz_2+1)+0;
		int ind_k3_0 = b*Mx+a;
		for (int ic=0; ic<=1; ic++){ //ic=0 => Real part, ic=1 => Imaginary part
			local_uk_k3_0[ind_k3_0][ic] = uk_1[ind][ic];
			local_vk_k3_0[ind_k3_0][ic] = vk_1[ind][ic];
			local_wk_k3_0[ind_k3_0][ic] = wk_1[ind][ic];
		}
	}
	//Gathering of velocity Fourier coefficients at plane k3=0 and distribution to all processes
	MPI_Allgatherv(local_uk_k3_0, recvcounts_k3_0[myrank], MPI_BYTE, uk_k3_0, recvcounts_k3_0, displs_k3_0, MPI_BYTE, MCW);
	MPI_Allgatherv(local_vk_k3_0, recvcounts_k3_0[myrank], MPI_BYTE, vk_k3_0, recvcounts_k3_0, displs_k3_0, MPI_BYTE, MCW);
	MPI_Allgatherv(local_wk_k3_0, recvcounts_k3_0[myrank], MPI_BYTE, wk_k3_0, recvcounts_k3_0, displs_k3_0, MPI_BYTE, MCW);
	//Complex conjugate correction for k3=0
	LOOP_2D(b,0,local_n1, a,0,Mx) { //For plane k3=0
		int k1 = FOURIER_FREQ(a,Mx);
		int k2 = FOURIER_FREQ(b,local_1_start,My);
		if (k1<0 || (k1==0 && k2<0)) {
			int k1_conj = -k1, k2_conj = -k2;
			int ind = (b*Mx+a)*(Mz_2+1)+0;
			int a_conj = k1_conj; //k1 <= 0 => k1_conj >= 0
			int b_conj = ((k2_conj >= 0) ? k2_conj : (My+k2_conj));
			int ind_k3_0 = b_conj*Mx+a_conj;
			uk_1[ind][0] = uk_k3_0[ind_k3_0][0];  uk_1[ind][1] = -uk_k3_0[ind_k3_0][1];
			vk_1[ind][0] = vk_k3_0[ind_k3_0][0];  vk_1[ind][1] = -vk_k3_0[ind_k3_0][1];
			wk_1[ind][0] = wk_k3_0[ind_k3_0][0];  wk_1[ind][1] = -wk_k3_0[ind_k3_0][1];
		}
	}
};


//=====================================================================================================================
// INPUTTING INITIAL VELOCITY FIELD
//=====================================================================================================================

//Initialization of velocity field in Fourier space using dummy distibution 1/|k|
void HIT::Input_Dummy_Field() {
	LOOP_FOURIER {
		if (dealiased[ind]) {
			uk_1[ind][0] = ((rad2[ind]>0.0) ? (1.0/(rad2[ind])) : 0.0);  uk_1[ind][1] = 0.0;
			vk_1[ind][0] = ((rad2[ind]>0.0) ? (1.0/(rad2[ind])) : 0.0);  vk_1[ind][1] = 0.0;
			wk_1[ind][0] = ((rad2[ind]>0.0) ? (1.0/(rad2[ind])) : 0.0);  wk_1[ind][1] = 0.0;
		} else { //Aliased terms
			uk_1[ind][0] = 0.0;  uk_1[ind][1] = 0.0;
			vk_1[ind][0] = 0.0;  vk_1[ind][1] = 0.0;
			wk_1[ind][0] = 0.0;  wk_1[ind][1] = 0.0;
		}
	}
};

//Initialization of velocity field in Fourier space from an input file containing velocities in Fourier space
void HIT::Input_Real_Field() {
	//Create pointers for Fourier-space velocity
	COMPLEX *uk_file = alloc_complex(Nx_file * Ny_file * (Nz_file/2+1));
	COMPLEX *vk_file = alloc_complex(Nx_file * Ny_file * (Nz_file/2+1));
	COMPLEX *wk_file = alloc_complex(Nx_file * Ny_file * (Nz_file/2+1));
	//Create pointers for physical-space velocity
	REAL *u_file = NULL, *v_file = NULL, *w_file = NULL;
	if (myrank == 0) {
		u_file = alloc_real(Nx_file * Ny_file * Nz_file);
		v_file = alloc_real(Nx_file * Ny_file * Nz_file);
		w_file = alloc_real(Nx_file * Ny_file * Nz_file);
	}
	//Create fftw3 plans (FFT of physical velocities)
	plan r2c_u_file = NULL, r2c_v_file = NULL, r2c_w_file = NULL;
	if (myrank == 0) {
		r2c_u_file = plan_dft_r2c_3d(Nx_file, Ny_file, Nz_file, u_file, uk_file, FE|FDI);
		r2c_v_file = plan_dft_r2c_3d(Nx_file, Ny_file, Nz_file, v_file, vk_file, FE|FDI);
		r2c_w_file = plan_dft_r2c_3d(Nx_file, Ny_file, Nz_file, w_file, wk_file, FE|FDI);
	}
	//Load of input file (velocity in physical space). Only process 0 reads the (whole) file
	int dumb_columns = 0;
	Read_Input_File_Vxyz(u_file, v_file, w_file, dumb_columns);
	//Execute FFT
	if (myrank == 0) {
		execute(r2c_u_file);
		execute(r2c_v_file);
		execute(r2c_w_file);
		//Freeing of mallocated memory
		destroy_plan(r2c_u_file);
		destroy_plan(r2c_v_file);
		destroy_plan(r2c_w_file);
		FREE(u_file);
		FREE(v_file);
		FREE(w_file);
		//Normalization following FFT
		for (int ind=0; ind<Nx_file*Ny_file*(Nz_file/2+1); ind++) {
			for (int ic=0; ic<=1; ic++){ //ic=0 => Real part, ic=1 => Imaginary part
				uk_file[ind][ic] /= (Nx_file*Ny_file*Nz_file);
				vk_file[ind][ic] /= (Nx_file*Ny_file*Nz_file);
				wk_file[ind][ic] /= (Nx_file*Ny_file*Nz_file);
			}
		}
		Calculate_Ek_init_file(uk_file, vk_file, wk_file);
	}
	//Share Ek_init_file and pseudoEpsilon_init_file from process 0 to the rest of processes
	MPI_Bcast(&Ek_init_file, 1, REAL_MPI, 0, MCW);
	MPI_Bcast(&pseudoEpsilon_init_file, 1, REAL_MPI, 0, MCW);
	//Share vectors uk_file, vk_file and wk_file from process 0 to the rest of processes
	MPI_Bcast(uk_file, sizeof(COMPLEX) *  Nx_file * Ny_file * (Nz_file/2+1), MPI_BYTE, 0, MCW);
	MPI_Bcast(vk_file, sizeof(COMPLEX) *  Nx_file * Ny_file * (Nz_file/2+1), MPI_BYTE, 0, MCW);
	MPI_Bcast(wk_file, sizeof(COMPLEX) *  Nx_file * Ny_file * (Nz_file/2+1), MPI_BYTE, 0, MCW);
	//Truncation/Extension of input field (v in Fourier space)
	Truncate_Input_File_Fourier_Coefficients(uk_file, vk_file, wk_file);
	//Freeing of mallocated memory
	FREE(uk_file);
	FREE(vk_file);
	FREE(wk_file);
};

//=========================
// 1.) LOAD PHYSICAL VELOCITY FROM FILE (All by read by one rank!!)
// Read and load of velocity fields from input files of known size
void HIT::Read_Input_File_Vxyz (REAL* u_file, REAL* v_file, REAL* w_file, const int& dumb_columns) {
	REAL *u_read = NULL, *v_read = NULL, *w_read = NULL;
	unsigned int total_file_num = ActualNx_file * ActualNy_file * ActualNz_file;
	if (myrank == 0) {
		if (isNotCubic) { //isNotCubic => ActualN <= N
			u_read = alloc_real(ActualNx_file * ActualNy_file * ActualNz_file);
			v_read = alloc_real(ActualNx_file * ActualNy_file * ActualNz_file);
			w_read = alloc_real(ActualNx_file * ActualNy_file * ActualNz_file);
		} else { //Cubic => ActualN = N
			u_read = u_file;
			v_read = v_file;
			w_read = w_file;
		}
	}
	if (isInitialFieldFromBinaryFile) {
		unsigned int expected_filesize = sizeof(REAL) * 3 * total_file_num;
		unsigned int expected_single_filesize = sizeof(float) * 3 * total_file_num;
		unsigned int expected_double_filesize = sizeof(double) * 3 * total_file_num;
		MPI_File fh;
		MPI_Status status;
		MPI_File_open(MCW, Binary_Input_Filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
		//Check the size of the input file
		MPI_Offset filesize;
		MPI_File_get_size(fh, &filesize);
		if ((filesize != expected_single_filesize) && (filesize != expected_double_filesize)) crash("Incorrect size of \"%s\"!! %lld bytes instead of %d bytes (for single-precision) or %d bytes (for double-precision).\n", Binary_Input_Filename, filesize, expected_single_filesize, expected_double_filesize);
		//Read of input file (only by rank 0)
		MPI_Offset offset_Real_read = 0;
		MPI_File_seek(fh, offset_Real_read, MPI_SEEK_SET);
		int num_to_read;
		if (myrank == 0) num_to_read = total_file_num; //everything read by rank 0
		else num_to_read = 0;
		//Binary read of input velocity field
		if (filesize != expected_filesize) { //calculation precision != input precision
			switch (BYTES) {
				case 4: //single-precision simulation but double-precision input
				{
					double *u_input = NULL, *v_input = NULL, *w_input = NULL;
					if (myrank == 0) {
						u_input = fftw_alloc_real(total_file_num);
						v_input = fftw_alloc_real(total_file_num);
						w_input = fftw_alloc_real(total_file_num);
					}
					MPI_File_read_all(fh, u_input, num_to_read, MPI_DOUBLE, &status);
					MPI_File_read_all(fh, v_input, num_to_read, MPI_DOUBLE, &status);
					MPI_File_read_all(fh, w_input, num_to_read, MPI_DOUBLE, &status);
					if (myrank == 0) {
						for (int ind=0; ind<ActualNx_file*ActualNy_file*ActualNz_file; ind++) {
							u_read[ind] = u_input[ind];
							v_read[ind] = v_input[ind];
							w_read[ind] = w_input[ind];
						}
						fftw_free(u_input); fftw_free(v_input); fftw_free(w_input);
					}
					break;
				}
				case 8: //double-precision simulation but single-precision input
				{
					float *u_input = NULL, *v_input = NULL, *w_input = NULL;
					if (myrank == 0) {
						u_input = fftwf_alloc_real(total_file_num);
						v_input = fftwf_alloc_real(total_file_num);
						w_input = fftwf_alloc_real(total_file_num);
					}
					MPI_File_read_all(fh, u_input, num_to_read, MPI_FLOAT, &status);
					MPI_File_read_all(fh, v_input, num_to_read, MPI_FLOAT, &status);
					MPI_File_read_all(fh, w_input, num_to_read, MPI_FLOAT, &status);
					if (myrank == 0) {
						for (int ind=0; ind<ActualNx_file*ActualNy_file*ActualNz_file; ind++) {
							u_read[ind] = u_input[ind];
							v_read[ind] = v_input[ind];
							w_read[ind] = w_input[ind];
						}
						fftwf_free(u_input); fftwf_free(v_input); fftwf_free(w_input);
					}
					break;
				}
			}
		} else { //calculation precision == input precision
			MPI_File_read_all(fh, u_read, num_to_read, REAL_MPI, &status);
			MPI_File_read_all(fh, v_read, num_to_read, REAL_MPI, &status);
			MPI_File_read_all(fh, w_read, num_to_read, REAL_MPI, &status);
		}
		MPI_File_close(&fh);
	} else if (isInitialFieldFromASCIIFile) {
		if (myrank == 0) {
			int row = 0;
			int col = 0;
			REAL shit;
			ifstream infile;
			infile.open(ASCII_Input_Filename);
			if(infile.fail()) {
				crash("Error opening \"%s\"!!\n", ASCII_Input_Filename);
			}
			while((!infile.eof()) && (row < ActualNx_file*ActualNy_file*ActualNz_file)) {
				if (col < dumb_columns) {
					infile >> shit;
					col++;
				} else {
					infile >> u_read[row] >> v_read[row] >> w_read[row];
					row++;
					col=0;
				}
			}
			infile.close();
		}
	}
	if (isNotCubic) {
		if(myrank == 0) {
			LOOP_3D (i,0,Nx_file, j,0,Ny_file, k,0,Nz_file) {
				int ind = (i*Ny_file+j)*Nz_file+k;
				if (i<ActualNx_file && j<ActualNy_file && k<ActualNz_file) {
					int ind_read = (i*ActualNy_file+j)*ActualNz_file+k;
					u_file[ind] = u_read[ind_read];
					v_file[ind] = v_read[ind_read];
					w_file[ind] = w_read[ind_read];
				} else {
					int i_ = ((i < ActualNx_file) ? i : i-ActualNx_file);
					int j_ = ((j < ActualNy_file) ? j : j-ActualNy_file);
					int k_ = ((k < ActualNz_file) ? k : k-ActualNz_file);
					int ind_ = (i_*Ny_file+j_)*Nz_file+k_;
					u_file[ind] = u_file[ind_];
					v_file[ind] = v_file[ind_];
					w_file[ind] = w_file[ind_];
				}
			}
			FREE(u_read); FREE(v_read); FREE(w_read);
		}
	}
};

//=========================
// 2.) TRUNCATION/EXTENSION IN FOURIER SPACE OF INPUT VELOCITY FIELD
// Truncate/Extend input velocity field in Fourier space (accordingly to Nx, Ny, Nz)
// Afterwards, uk_1, vk_1 and wk_1 are completely initialized
void HIT::Truncate_Input_File_Fourier_Coefficients (COMPLEX const * const uk_file, COMPLEX const * const vk_file, COMPLEX const * const wk_file) {
	// TRUNCATION:
	LOOP_FOURIER_k1k2k3 {
		for (int ic=0; ic<=1; ic++) { //ic=0 => Real part, ic=1 => Imaginary part
			if (dealiased[ind]) {
				if (((abs(k1) <= (Nx_file-1)/2) || (k1 == Nx_file/2)) && // In case Nx_file is even, or clause: (k1 == Nx_file/2)
					((abs(k2) <= (Ny_file-1)/2) || (k2 == Ny_file/2)) && // In case Ny_file is even, or clause: (k2 == Ny_file/2)
					((abs(k3) <= (Nz_file-1)/2) || (k3 == Nz_file/2))) { // In case Nz_file is even, or clause: (k3 == Nz_file/2)
					int a_file = ((k1 >= 0) ? a : (Nx_file+k1)); // k1<0 => Nx_file+k1 = Nx_file-|k1|
					int b_file = ((k2 >= 0) ? (local_1_start+b) : (Ny_file+k2)); // k2<0 => Ny_file+k2 = Ny_file-|k2|
					int c_file = ((k3 >= 0) ? k3 : (Nz_file+k3)); // k3<0 => Nz_file+k3 = Nz_file-|k3|
					int ind_file = (a_file*Ny_file+b_file)*(Nz_file/2+1)+c_file;
					uk_1[ind][ic] = uk_file[ind_file][ic];
					vk_1[ind][ic] = vk_file[ind_file][ic];
					wk_1[ind][ic] = wk_file[ind_file][ic];
				} else { //Extension of input data (when N* > N*_file)
					uk_1[ind][ic] = 0.0;
					vk_1[ind][ic] = 0.0;
					wk_1[ind][ic] = 0.0;	
				}
			} else { //Aliased terms (abs(k*) = N*_2+1, ..., M*_2)
				uk_1[ind][ic] = 0.0;
				vk_1[ind][ic] = 0.0;
				wk_1[ind][ic] = 0.0;					
			}
		}
	}
};

//=====================================================================================================================
// AUXILIAR FUNCTIONS
//=====================================================================================================================

//Initialization of the clas HIT (called uniquely in the constructor)
void HIT::HIT_init() {
	//===================================================
	// 1. GENERAL MEMORY ALLOCATION
	//===================================================
	// fftw3 parameters
	alloc_local = fftw.getalloc_local();
	local_n0 = fftw.getlocal_n0();
	local_0_start = fftw.getlocal_0_start();
	local_n1 = fftw.getlocal_n1();
	local_1_start = fftw.getlocal_1_start();
	// double*
	u = alloc_real(2 * alloc_local);
	v = alloc_real(2 * alloc_local);
	w = alloc_real(2 * alloc_local);
	for(int coord=0; coord<3; coord++) {
		gradu[coord] = alloc_real(2*alloc_local);
		gradv[coord] = alloc_real(2*alloc_local);
		gradw[coord] = alloc_real(2*alloc_local);
	}
	convx = alloc_real(2 * alloc_local);
	convy = alloc_real(2 * alloc_local);
	convz = alloc_real(2 * alloc_local);
	rad2 = alloc_real(2 * alloc_local);
	local_Ek = alloc_real(last_rad+1);
	Ek = alloc_real(last_rad+1);
	// COMPLEX*
	uk_aux = alloc_complex(alloc_local);
	vk_aux = alloc_complex(alloc_local);
	wk_aux = alloc_complex(alloc_local);
	uk_0 = alloc_complex(alloc_local);
	vk_0 = alloc_complex(alloc_local);
	wk_0 = alloc_complex(alloc_local);
	uk_1 = alloc_complex(alloc_local);
	vk_1 = alloc_complex(alloc_local);
	wk_1 = alloc_complex(alloc_local);
	for(int coord=0; coord<3; coord++) {
		gradu_k[coord] = alloc_complex(alloc_local);
		gradv_k[coord] = alloc_complex(alloc_local);
		gradw_k[coord] = alloc_complex(alloc_local);
	}
	convx_k = alloc_complex(alloc_local);
	convy_k = alloc_complex(alloc_local);
	convz_k = alloc_complex(alloc_local);
	Rx0_k = alloc_complex(alloc_local);
	Ry0_k = alloc_complex(alloc_local);
	Rz0_k = alloc_complex(alloc_local);
	Rx1_k = alloc_complex(alloc_local);
	Ry1_k = alloc_complex(alloc_local);
	Rz1_k = alloc_complex(alloc_local);
	//bool*
	dealiased = (bool*) MALLOC(sizeof(bool) * alloc_local);
	conjugate = (bool*) MALLOC(sizeof(bool) * alloc_local);

	//===================================================
	// 2. INITIALIZATIONS OF AUXILIAR INTERNAL VARIABLES
	//===================================================
	//rad2: Initialization of vector radium (to avoid repeated calculations)
	LOOP_FOURIER_k1k2k3 {
		rad2[ind] = static_cast<REAL>(k1*k1+k2*k2+k3*k3);
	}
	//num_dealiased_real: number of relevant terms (considering the length factors) per process
	if (isNotCubic) {
		num_dealiased_real = 0;
		LOOP_3D (i,0,local_n0, j,0,My, k,0,Mz) {
			if ((local_0_start+i)<ActualMx && j<ActualMy && k<ActualMz) {
				num_dealiased_real++;
			}
		}
	} else {
		num_dealiased_real = local_n0 * My * Mz;
	}

	//dealiased: Initialization of a bool vector to know wether a mode is significant. Two cases:
	//A) Exclusion of aliased terms:
	//	-If N is odd, all k in [-N_2, N_2] are included, and abs(k) in [N_2, M_2] are excluded
	//	-If N is even, k=-N_2 is also excluded (complex z-dimension is only computed from [0,N_2])
	//(*If Nz is even, no exclusion has to be made as only k3>0 are considered and k3 == -Nz_2 cannot happen)
	//B) Exclusion of undesired coefficients (due to non-cubic domains, i.e., length factors >1):
	//	-FFTW assumes the same length for x, y and z-directions (in the exponential term of DFTs)
	//	-This is overcome by setting to 0 all modes irrelevant to each direction (dealiased[ind]=false)
	//C) num_dealiased: Number of dealiased terms present in current process
	num_dealiased = 0;
	LOOP_FOURIER_k1k2k3 {
		//Exclusion of aliased terms (3/2 rule implemented)
		if (abs(k1)<=Nx_2 && abs(k2)<=Ny_2 && k3<=Nz_2) {
			dealiased[ind] = true;
			if (Nx%2==0 && k1==-Nx_2) {
				dealiased[ind] =false; //Exclusion in case Nx is even
			}
			if (Ny%2==0 && k2==-Ny_2) {
				dealiased[ind] =false; //Exclusion in case Ny is even
			}
		} else {
			dealiased[ind] = false;
		}
		//Until here, (although it has been separated to make it clearer) dealiased[ind] is equivalent
		//to: "((abs(k1)<=(Nx-1)/2)||(k1==Nx_2)) && ((abs(k2)<=(Ny-1)/2)||(k2==Ny_2)) && (k3<=Nz_2)"
		if (isNotCubic) {
			if ((k1%Lx_factor != 0) || (k2%Ly_factor != 0) || (k3%Lz_factor != 0)){
				dealiased[ind] = false; //Exclusion of undesired coefficients due to length factors
			}
		}
		if (dealiased[ind]) num_dealiased++; //It needs to be at the end!!!
	}

	//conjugate: Initialization of a bool vector to know wether a mode also has an associated conjugate (and 
	//it has to be counted twice when computing the kinetic energy) or not
	LOOP_FOURIER_k1k2k3 {
		if (dealiased[ind]) {
			if ((k3!=0) && (abs(k1)<=(Nx-1)/2) && (abs(k2)<=(Ny-1)/2) && (abs(k3)<=(Nz-1)/2)) {
				conjugate[ind] = true;
			} else {
				conjugate[ind] = false; //Exclusion in case Nx, Ny or Nz are even
			}
		} else {
			conjugate[ind] = false; //Dummy case
		}
	}

	//local_n0_: pointer containing in local_n0_[i] the value of local_n0 corresponding to process i
	local_n0_ = (int*) MALLOC(sizeof(int) * numprocs);
	MPI_Gather(&local_n0, 1, MPI_INT, local_n0_, 1, MPI_INT, 0, MCW);
	MPI_Bcast(local_n0_, numprocs, MPI_INT, 0, MCW);

	//local_n1_: pointer containing in local_n1_[i] the value of local_n1 corresponding to process i
	local_n1_ = (int*) MALLOC(sizeof(int) * numprocs);
	MPI_Gather(&local_n1, 1, MPI_INT, local_n1_, 1, MPI_INT, 0, MCW);
	MPI_Bcast(local_n1_, numprocs, MPI_INT, 0, MCW);

	//offset_Fourier: Position of the first byte to be written by the current proces. Used to MPI_File_write uk_1, vk_1
	//and wk_1. Only dealiased terms are written!!!
	MPI_Offset *offset_Fourier_ = (MPI_Offset*) MALLOC(sizeof(MPI_Offset) * numprocs);
	int *num_dealiased_ = NULL;
	if (myrank == 0) {
		num_dealiased_ = (int*) MALLOC(sizeof(int) * numprocs);
	}
	//Rank 0 gathers each num_dealiased into the temporary pointer *num_dealiased_
	MPI_Gather(&num_dealiased, 1, MPI_INT, num_dealiased_, 1, MPI_INT, 0, MCW);
	if(myrank == 0){
		for(int proc=0; proc<numprocs; proc++) {
			offset_Fourier_[proc] = ((proc == 0) ? 0 : (offset_Fourier_[proc-1] + sizeof(COMPLEX) * static_cast<MPI_Offset>(num_dealiased_[proc-1])));
		}
	FREE(num_dealiased_);
	}
	MPI_Bcast(offset_Fourier_, sizeof(MPI_Offset)*numprocs, MPI_BYTE, 0, MCW);
	//Final initialization of the variables
	offset_Fourier = offset_Fourier_[myrank];
	//Memory freeing of auxiliar pointers
	FREE(offset_Fourier_);

	//offset_Real: Exactly as offset_Fourier but for real fields (where distinction between 
	//aliased/dealiased terms does not make sense, except for non-cubic domains, where terms
	//such that k%L_factor!=0 have to be treated as the complex aliased ones)
	//Allocate memory
	MPI_Offset *offset_Real_ = (MPI_Offset*) MALLOC(sizeof(MPI_Offset) * numprocs);
	int *num_dealiased_real_ = NULL;
	if (myrank == 0) {
		num_dealiased_real_ = (int*) MALLOC(sizeof(int) * numprocs);
	}
	//Rank 0 gathers each num_dealiased_real into the temporary pointer *num_dealiased_real_
	MPI_Gather(&num_dealiased_real, 1, MPI_INT, num_dealiased_real_, 1, MPI_INT, 0, MCW);
	if(myrank == 0){
		for(int proc=0; proc<numprocs; proc++) {
			offset_Real_[proc] = ((proc == 0) ? 0 : (offset_Real_[proc-1] + sizeof(REAL) * static_cast<MPI_Offset>(num_dealiased_real_[proc-1])));
		}
	FREE(num_dealiased_real_);
	}
	MPI_Bcast(offset_Real_, sizeof(MPI_Offset)*numprocs, MPI_BYTE, 0, MCW);
	//Final initialization of the variables
	offset_Real = offset_Real_[myrank];
	//Memory freeing of auxiliar pointers
	FREE(offset_Real_);

	//===================================================
	// 3. Initial velocity field
	//===================================================
	if(isInitialFieldFromBinaryFile || isInitialFieldFromASCIIFile) {
		//Input velocity (it needs to be after initialization of *dealiased)
		Input_Real_Field();
	} else {
		Input_Dummy_Field();
	}

	//Trivial initial distribution of R vector
	for (int ind=0; ind<local_n1*Mx*(Mz_2+1); ind++) {
		for (int ic=0; ic<=1; ic++){ //ic=0 => Real part, ic=1 => Imaginary part
			Rx1_k[ind][ic] = 0.0;
			Ry1_k[ind][ic] = 0.0;
			Rz1_k[ind][ic] = 0.0;
//			uk_0[ind][ic] = 0.0;
//			vk_0[ind][ic] = 0.0;
//			wk_0[ind][ic] = 0.0;
		}
	}

	//===================================================
	// 4.1 LES MEMORY ALLOCATION
	//===================================================
	if (isLES) {
		// double*
		nu_eddy = alloc_real(2 * alloc_local);
		SGS_11 = alloc_real(2 * alloc_local);
		SGS_21 = alloc_real(2 * alloc_local);
		SGS_22 = alloc_real(2 * alloc_local);
		SGS_31 = alloc_real(2 * alloc_local);
		SGS_32 = alloc_real(2 * alloc_local);
		SGS_33 = alloc_real(2 * alloc_local);
		// COMPLEX*
		SGS_11_k = alloc_complex(alloc_local);
		SGS_21_k = alloc_complex(alloc_local);
		SGS_22_k = alloc_complex(alloc_local);
		SGS_31_k = alloc_complex(alloc_local);
		SGS_32_k = alloc_complex(alloc_local);
		SGS_33_k = alloc_complex(alloc_local);
		divSGSx_k = alloc_complex(alloc_local);
		divSGSy_k = alloc_complex(alloc_local);
		divSGSz_k = alloc_complex(alloc_local);
	}

	//===================================================
	// 4.2 COMPLEX CONJUGATE CORRECTION
	//===================================================
	if (isComplexConjugateCorrected) {
		// COMPLEX*
		local_uk_k3_0 = alloc_complex(Mx * local_n1 * 1); //Only for k3=0
		local_vk_k3_0 = alloc_complex(Mx * local_n1 * 1); //Only for k3=0
		local_wk_k3_0 = alloc_complex(Mx * local_n1 * 1); //Only for k3=0
		uk_k3_0 = alloc_complex(Mx * My * 1); //Only for k3=0
		vk_k3_0 = alloc_complex(Mx * My * 1); //Only for k3=0
		wk_k3_0 = alloc_complex(Mx * My * 1); //Only for k3=0
		//int*
		displs_k3_0 = (int*) MALLOC(sizeof(int) * numprocs); //Position of the first byte gathered
		recvcounts_k3_0 = (int*) MALLOC(sizeof(int) * numprocs); //Number of MPI_BYTEs (for k3=0)
		//Initialization of recvcounts_k3_0 and displs_k3_0 (used in MPI_Allgatherv clause)
		//MPI_Gather is necessary because local_n1 differs from different processes
		if(myrank==0){
			for(int proc=0; proc<numprocs; proc++) {
				recvcounts_k3_0[proc] = sizeof(COMPLEX) * Mx * local_n1_[proc]; //Equal to: sizeof(COMPLEX) * Mx * local_n1 * 1 (*1 because only for k3=0)
				displs_k3_0[proc] = ((proc == 0) ? 0 : (displs_k3_0[proc-1]+recvcounts_k3_0[proc-1]));
			}
		}
		MPI_Bcast(recvcounts_k3_0, numprocs, MPI_INT, 0, MCW);
		MPI_Bcast(displs_k3_0,  numprocs, MPI_INT, 0, MCW);
	}

	//===================================================
	// 4.3 TIMESTEP
	//===================================================
	if (isSelfAdaptive) {
		//Related constants:
		Ah = (Ax + Ay + Az) / 3.0;
		t1 = 0.9302468;
		kappa = 0.5; //Initial value of kappa
		// double*
		k = alloc_real(2);
		phi = alloc_real(3);
		c = alloc_real(14);
		// Initializaiton
		k[0] = 0.73782212;
		k[1] = 0.44660387;
		phi[0] = atan(164.0/99.0);
		phi[1] = PI/3.0;
		phi[2] = PI*3.0*3.0/(5.0*5.0);
		c[0] = 0.0647998;
		c[1] = -0.386022;
		c[2] = 3.72945;
		c[3] = -9.38143;
		c[4] = 7.06574;
		c[5] = 2403400.0;
		c[6] = -5018490.0;
		c[7] = 2620140.0;
		c[8] = 2945.0;
		c[9] = -6665.76;
		c[10] = 3790.54;
		c[11] = 4.80513;
		c[12] = -16.9473;
		c[13] = 15.0155;
		//Initial self-adaptive timestep
		Recalculate_Velocity_Antitransform(); //As only uk_1,vk_1,wk_1 were already available
		Recalculate_SelfAdaptive_TimeStep();
	} else {
		//Related constants:
		C_visc = 0.2 * min(Ax*Ax, min (Ay*Ay, Az*Az)) / nu;
		//Initial CFL timestep
		Recalculate_Velocity_Antitransform(); //As only uk_1,vk_1,wk_1 were already available
		Recalculate_CFL_TimeStep();
	}

	Calculate_Ek_init();
#if(0)
	if (myrank == 0) {
		cout << "Initial kinetic energy: " << Ek_init << endl;
	}
#endif
};
//Deletion of the clas HIT (called uniquely in the destructor)
void HIT::HIT_destroy() {
	//Freeing memory
	FREE(u); FREE(v); FREE(w);
	FREE(uk_aux); FREE(vk_aux); FREE(wk_aux);
	FREE(uk_0); FREE(vk_0); FREE(wk_0);
	FREE(uk_1); FREE(vk_1); FREE(wk_1);
	FREE(convx); FREE(convy); FREE(convz);
	FREE(convx_k); FREE(convy_k); FREE(convz_k);
	FREE(Rx0_k); FREE(Ry0_k); FREE(Rz0_k);
	FREE(Rx1_k); FREE(Ry1_k); FREE(Rz1_k);
	FREE(dealiased); FREE(rad2); FREE(local_Ek); FREE(Ek);
	for (int coord=0; coord<3; coord++) {
		FREE(gradu[coord]); FREE(gradv[coord]); FREE(gradw[coord]);
		FREE(gradu_k[coord]); FREE(gradv_k[coord]); FREE(gradw_k[coord]);
	}
	if (isLES) {
		FREE(nu_eddy);
		FREE(SGS_11); FREE(SGS_21); FREE(SGS_22);
		FREE(SGS_31); FREE(SGS_32); FREE(SGS_33);
		FREE(SGS_11_k); FREE(SGS_21_k); FREE(SGS_22_k);
		FREE(SGS_31_k); FREE(SGS_32_k); FREE(SGS_33_k);
		FREE(divSGSx_k); FREE(divSGSy_k); FREE(divSGSz_k);
	}
	if (isComplexConjugateCorrected) {
		FREE(local_uk_k3_0); FREE(local_vk_k3_0); FREE(local_wk_k3_0);
		FREE(uk_k3_0); FREE(vk_k3_0); FREE(wk_k3_0);
		FREE(displs_k3_0); FREE(recvcounts_k3_0);
	}
	if (isSelfAdaptive) {
		FREE(k); FREE(phi); FREE(c);
	}
};

//Calculation of kinetic energy (total, Ek_Tot, and modal, Ek[K])
void HIT::Recalculate_Energy_Cascade() {
	//Reset of vector local_Ek[K]
	for (int K=0; K<=last_rad; K++) {
		local_Ek[K] = 0.0;
	}
	//Recalculation of vector local_Ek[K]
	LOOP_FOURIER {
		if (dealiased[ind]) {
			int K = static_cast<int>(round(sqrt(rad2[ind])));
			if (K<=last_rad) {
				REAL aux =  (uk_1[ind][0] * uk_1[ind][0]) + (uk_1[ind][1] * uk_1[ind][1]);
				aux += (vk_1[ind][0] * vk_1[ind][0]) + (vk_1[ind][1] * vk_1[ind][1]);
				aux += (wk_1[ind][0] * wk_1[ind][0]) + (wk_1[ind][1] * wk_1[ind][1]);
				if (conjugate[ind]) {
					local_Ek[K] += aux;
				} else {
					local_Ek[K] += (0.5 * aux);
				}
			}
		}
	}
	//Sum accross processes
	for (int K=0; K<=last_rad; K++) {
		MPI_Allreduce(&local_Ek[K], &Ek[K], 1, REAL_MPI, MPI_SUM, MCW);
	}
};
//Calculation of total kinetic energy from Ek[K]
void HIT::Recalculate_Total_Kinetic_Energy() {
	Ek_Tot = 0.0;
	for (int K=0; K<=last_rad; K++) {
		Ek_Tot += Ek[K];
	}
};
// Calculation of the initial kinetic energy of the system
void HIT::Calculate_Ek_init () {
	//Calculates Ek_Tot given the distribution of velocity in Fourier space
	Recalculate_Energy_Cascade();
	//Calculates the total kinetic energy summing Ek accros |k|
	Recalculate_Total_Kinetic_Energy();
	//Updates the value of Ek_init with the current total kinetic energy
	Ek_init = Ek_Tot;
};

// Calculation of the total kinetic energy of initial velocity field loaded from file (all done by process 0 in Input_Real_Field())
void HIT::Calculate_Ek_init_file (COMPLEX const * const uk_file, COMPLEX const * const vk_file, COMPLEX const * const wk_file) {
	Ek_init_file = 0.0;
	pseudoEpsilon_init_file = 0.0;
	//Sum accross Fourier modes
	LOOP_3D (a,0,Nx_file, b,0,Ny_file, k3,0,Nz_file/2+1) {
		int ind = (a*Ny_file+b)*(Nz_file/2+1)+k3;
		int k1 = FOURIER_FREQ(a,Nx_file);
		int k2 = FOURIER_FREQ(b,Ny_file);
		int rad2_file = k1*k1+k2*k2+k3*k3;
		bool relevant_file = ((k1%Lx_factor == 0) || (k2%Ly_factor == 0) || (k3%Lz_factor == 0));
		if (relevant_file) {
			//bool conjugate_file = ((k3!=0) && (abs(k1)<=(Nx_file-1)/2) && (abs(k2)<=(Ny_file-1)/2) && (abs(k3)<=(Nz_file-1)/2));
			bool conjugate_file = ((k3!=0) && (k1!=16) && (k2!=16) && (k3!=16));
			REAL Ek_file =  (uk_file[ind][0] * uk_file[ind][0]) + (uk_file[ind][1] * uk_file[ind][1]);
			Ek_file += (vk_file[ind][0] * vk_file[ind][0]) + (vk_file[ind][1] * vk_file[ind][1]);
			Ek_file += (wk_file[ind][0] * wk_file[ind][0]) + (wk_file[ind][1] * wk_file[ind][1]);
			if (conjugate_file) { //Ek = 0.5 * |uk|^2 (conjugate=true => implicit avoided mode and Ek counted twice)
				Ek_init_file += Ek_file;
				pseudoEpsilon_init_file += (rad2_file * Ek_file);
			} else {
				Ek_init_file += (0.5 * Ek_file);
				pseudoEpsilon_init_file += (rad2_file * 0.5 * Ek_file);
			}
		}
	}
};

