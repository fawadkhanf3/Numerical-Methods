#include "mex.h"
#include "conio.h"
#include "math.h"
#include "Data/atm76.h"
#include "Data/constants.h"
#include "Utilities/utility.h"

// Function Declaration
double **RK4_Fixed(double *tspan, double *xinit, int &length);
double *eoms(double t, double y[6]);

double ECI2INV[3][3];
const double w[3] = { Gwx, 0, 0 }, h = 0.5;
const int nvar = 6;
double Alt = 12e3;


// Mex-Function (Gateway to MATLAB)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *Tspan, *y0; // Pointers to Input Arrays
	double *ty;         // Pointer to Output Arrays

	Tspan = mxGetPr(prhs[0]); // Input Array 1
	y0 = mxGetPr(prhs[1]);    // Input Array 2

	int length;
	double **yp = RK4_Fixed(Tspan, y0, length);

	plhs[0] = mxCreateDoubleMatrix(length, 7, mxREAL); // Output Array
	ty = mxGetPr(plhs[0]);

	for (int i = 0; i < length; i++)
	{
		ty[0 * length + i] = yp[i][0];
		ty[1 * length + i] = yp[i][1];
		ty[2 * length + i] = yp[i][2];
		ty[3 * length + i] = yp[i][3];
		ty[4 * length + i] = yp[i][4];
		ty[5 * length + i] = yp[i][5];
		ty[6 * length + i] = yp[i][6];
	}

}

// Function Definition
double **RK4_Fixed(double *tspan, double *xinit, int &length) {

	double temp[6], K[4][nvar];
	double tsave[9000], ysave[6][9000], t, y[nvar];
	int icount = 0;

	ECI2INV[0][0] = cos(azim0)*cos(latd0);
	ECI2INV[0][1] = sin(azim0);
	ECI2INV[0][2] = cos(azim0)*sin(latd0);
	ECI2INV[1][0] = -sin(azim0)*cos(latd0);
	ECI2INV[1][1] = cos(azim0);
	ECI2INV[1][2] = -sin(azim0)*sin(latd0);
	ECI2INV[2][0] = -sin(latd0);
	ECI2INV[2][1] = 0.0;
	ECI2INV[2][2] = cos(latd0);

	t = tsave[icount] = tspan[0];
	y[0] = ysave[0][icount] = xinit[0];
	y[1] = ysave[1][icount] = xinit[1];
	y[2] = ysave[2][icount] = xinit[2];
	y[3] = ysave[3][icount] = xinit[3];
	y[4] = ysave[4][icount] = xinit[4];
	y[5] = ysave[5][icount] = xinit[5];

	while (Alt >= ALT_T) {


		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < nvar; j++) {
				K[i][j] = 0.0;
			}
		}

		// k1 = f(t, x)
		double *dots = eoms(t, y);
		for (int k = 0; k < nvar; k++) {K[0][k] = dots[k];}

		// k2 = f(t+1/2*h,x+1/2*k1*h)
		for (int k = 0; k < nvar; k++) {temp[k] = y[k] + 0.5*h*K[0][k];}

		dots = eoms(t + 0.5*h, temp);

		for (int k = 0; k < nvar; k++) {K[1][k] = dots[k];}

		// k3 = f(t+1/2*h,x+1/2*k2*h)

		for (int k = 0; k < nvar; k++) {temp[k] = y[k] + 0.5*h*K[1][k];}

		dots = eoms(t + 0.5*h, temp);

		for (int k = 0; k < nvar; k++) {K[2][k] = dots[k];}

		// k4 = f(t+h,x+k3*h)
		for (int k = 0; k < nvar; k++) {temp[k] = y[k] + h*K[2][k];}

		dots = eoms(t + h, temp);

		for (int k = 0; k < nvar; k++) {K[3][k] = dots[k];}

		// Correction
		for (int k = 0; k < nvar; k++) {
			temp[k] = y[k] + 1.0/6.0*h*(K[0][k] + 2.0*K[1][k] + 2.0*K[2][k] + K[3][k]);
		}

		icount++;
		tsave[icount] = t = t + h;
		for (int k = 0; k < nvar; k++) {
			ysave[k][icount] = y[k] = temp[k];
		}

		double latc, lngc, dumy, latd, dum, Pmag;
		Pmag = sqrt(SQR(y[0]) + SQR(y[1]) + SQR(y[2]));
		latc = atan2(y[0], sqrt(SQR(y[1]) + SQR(y[2])));
		lngc = atan2(y[1], -y[2]);
		dumy = (RE / Pmag)*((Gf*sin(2.0*latc)) + (SQR(Gf) * sin(4.0*latc))*(RE / Pmag - 0.2500));
		latd = latc + asin(dumy);
		dum = 1.0 - (Gf*SQR(sin(latc))) - (0.5*SQR(Gf) * SQR(sin(2.0*latc)))*(RE / Pmag - 0.2500);
		Alt = (Pmag - RE*dum);
		if (Alt <= ALT_T) break;

	}

	icount++;
	int counts = (icount)*(nvar + 1);
	double **t_y = new double*[counts];
	for (int i = 0; i < icount; i++)
	{
		t_y[i] = new double[nvar + 1];
		t_y[i][0] = tsave[i];
		t_y[i][1] = ysave[0][i];
		t_y[i][2] = ysave[1][i];
		t_y[i][3] = ysave[2][i];
		t_y[i][4] = ysave[3][i];
		t_y[i][5] = ysave[4][i];
		t_y[i][6] = ysave[5][i];

	}
	// Outputs
	length = icount;

	return t_y;
}
double *eoms(double t, double y[6]) {

	//
	double Px, Py, Pz, Vx, Vy, Vz, Pmag;
	Px = y[0];	Py = y[1];	Pz = y[2];
	Vx = y[3];	Vy = y[4];	Vz = y[5];
	Pmag = sqrt(SQR(Px) + SQR(Py) + SQR(Pz));

	double P[3] = { Px,Py,Pz }, V[3] = { Vx,Vy,Vz };

	double gx, gy, gz, dum1, dum2;
	dum1 = 1.5*J2*SQR(RE / Pmag);
	dum2 = 5 * SQR(Px / Pmag);
	gx = -GM*(Px / CUB(Pmag))*(1 + dum1*(3 - dum2));
	gy = -GM*(Py / CUB(Pmag))*(1 + dum1*(1 - dum2));
	gz = -GM*(Pz / CUB(Pmag))*(1 + dum1*(1 - dum2));

	//
	double latc, lngc, dumy, latd, dum, Altm;
	latc = atan2(Px, sqrt(SQR(Py) + SQR(Pz)));
	lngc = atan2(Py, -Pz);
	dumy = (RE / Pmag)*((Gf*sin(2.0*latc)) + (SQR(Gf) * sin(4.0*latc))*(RE / Pmag - 0.2500));
	latd = latc + asin(dumy);
	dum = 1.0 - (Gf*SQR(sin(latc))) - (0.5*SQR(Gf) * SQR(sin(2.0*latc)))*(RE / Pmag - 0.2500);
	Altm = (Pmag - RE*dum) / 1000.0;

	//
	double Ve[3], Ve_inv[3];
	Ve[0] = V[0] - (w[1] * P[2] - w[2] * P[1]);
	Ve[1] = V[1] - (w[2] * P[0] - w[0] * P[2]);
	Ve[2] = V[2] - (w[0] * P[1] - w[1] * P[0]);
	MatMult(3, 1, 3, *ECI2INV, Ve, Ve_inv);

	//
	double psi_inv, gama_inv;
	psi_inv = atan2(Ve_inv[1], Ve_inv[0]);
	gama_inv = asin(Ve_inv[2] / (sqrt(SQR(Ve_inv[0]) + SQR(Ve_inv[1]) + SQR(Ve_inv[2]))));

	//
	double aa[3][3], bb[3][3], Binv2vel[3][3],
		Tvel2eci[3][3], Vm, INV2ECI[3][3], vel2Binv[3][3];

	aa[0][0] = cos(psi_inv); aa[0][1] = sin(psi_inv); aa[0][2] = 0.0;
	aa[1][0] = -sin(psi_inv); aa[1][1] = cos(psi_inv); aa[1][2] = 0.0;
	aa[2][0] = 0.0;	aa[2][1] = 0.0; aa[2][2] = 1.0;

	bb[0][0] = cos(gama_inv); bb[0][1] = 0.0; bb[0][2] = sin(gama_inv);
	bb[1][0] = 0.0;	bb[1][1] = 1.0;	bb[1][2] = 0.0;
	bb[2][0] = -sin(gama_inv); bb[2][1] = 0.0; bb[2][2] = cos(gama_inv);

	MatMult(3, 3, 3, *bb, *aa, *Binv2vel);
	MatMult(3, 3, 3, *bb, *aa, *Binv2vel);

	Transpose(3, 3, ECI2INV, INV2ECI);
	Transpose(3, 3, Binv2vel, vel2Binv);

	MatMult(3, 3, 3, *INV2ECI, *vel2Binv, *Tvel2eci);

	Vm = sqrt(SQR(Ve[0]) + SQR(Ve[1]) + SQR(Ve[2]));

	// Atmospheric Data
	CAtm us76;
	double Temp, sos, rho, machc, atm_out[3];
	us76.atm(Altm, atm_out);
	Temp = atm_out[0];
	rho = atm_out[2];
	sos = 20.04680275919057*sqrt(Temp);
	machc = Vm / sos;

	//
	double Fxa, Fza, Faero_tb[3];
	Fxa = 0.5*rho*SQR(Vm) * Sref * Cx / mass_rv;
	Fza = 0.5*rho*SQR(Vm) * Sref * Cz / mass_rv;
	Faero_tb[0] = Fxa;	Faero_tb[1] = 0; Faero_tb[2] = Fza;

	//
	double body2vel[3][3], bodyrol[3][3], dummy[3], Der[3], DD[3];
	body2vel[0][0] = cos(alp_t); body2vel[0][1] = 0.0; body2vel[0][2] = sin(alp_t);
	body2vel[1][0] = 0.0; body2vel[1][1] = 1.0; body2vel[1][2] = 0.0;
	body2vel[2][0] = -sin(alp_t); body2vel[2][1] = 0.0;	body2vel[2][2] = cos(alp_t);

	bodyrol[0][0] = 1.0; bodyrol[0][1] = 0.0; bodyrol[0][2] = 0.0;
	bodyrol[1][0] = 0.0; bodyrol[1][1] = cos(phi_r); bodyrol[1][2] = -sin(phi_r);
	bodyrol[2][0] = 0.0; bodyrol[2][1] = sin(phi_r); bodyrol[2][2] = cos(phi_r);

	MatMult(3, 1, 3, *body2vel, Faero_tb, dummy);
	MatMult(3, 1, 3, *bodyrol, dummy, Der);
	MatMult(3, 1, 3, *Tvel2eci, Der, DD);

	//
	double *dydx = new double[6];
	dydx[0] = Vx; dydx[1] = Vy;	dydx[2] = Vz;
	dydx[3] = gx + DD[0]; dydx[4] = gy + DD[1]; dydx[5] = gz + DD[2];

	return dydx;
}