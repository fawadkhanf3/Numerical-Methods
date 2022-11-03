#include "../Data/constants.h"
#include "../Data/atm76.h"
#include "../Utilities/utility.h"

struct eoms {

	void operator() (const Doub x, VecDoub_I &y, VecDoub_O &dydx) {

		//
		double ECI2INV[3][3];
		ECI2INV[0][0] = cos(azim0)*cos(latd0);
		ECI2INV[0][1] = sin(azim0);
		ECI2INV[0][2] = cos(azim0)*sin(latd0);
		ECI2INV[1][0] = -sin(azim0)*cos(latd0);
		ECI2INV[1][1] = cos(azim0);
		ECI2INV[1][2] = -sin(azim0)*sin(latd0);
		ECI2INV[2][0] = -sin(latd0);
		ECI2INV[2][1] = 0.0;
		ECI2INV[2][2] = cos(latd0);

		//
		double Px, Py, Pz, Vx, Vy, Vz, Pmag;
		Px = y[0];
		Py = y[1];
		Pz = y[2];
		Vx = y[3];
		Vy = y[4];
		Vz = y[5];
		Pmag = sqrt(SQR(Px) + SQR(Py) + SQR(Pz));

		double P[3] = { Px,Py,Pz }, V[3] = { Vx,Vy,Vz };

		//
		double gx, gy, gz, dum1, dum2;
		dum1 = 1.5*J2*SQR(RE / Pmag);
		dum2 = 5 * SQR(Px / Pmag);
		gx = -GM*(Px / CUB(Pmag))*(1 + dum1*(3 - dum2));
		gy = -GM*(Py / CUB(Pmag))*(1 + dum1*(1 - dum2));
		gz = -GM*(Pz / CUB(Pmag))*(1 + dum1*(1 - dum2));

		//
		double latc, lngc, dumy, latd, dum, Alt;
		latc = atan2(Px, sqrt(SQR(Py) + SQR(Pz)));
		lngc = atan2(Py, -Pz);
		dumy = (RE / Pmag)*((Gf*sin(2.0*latc)) + (SQR(Gf) * sin(4.0*latc))*(RE / Pmag - 0.2500));
		latd = latc + asin(dumy);
		dum = 1.0 - (Gf*SQR(sin(latc))) - (0.5*SQR(Gf) * SQR(sin(2.0*latc)))*(RE / Pmag - 0.2500);
		Alt = (Pmag - RE*dum) / 1000.0;
        
		//
		const double w[3] = { Gwx, 0, 0 };
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
		
		aa[0][0] = cos(psi_inv);
		aa[0][1] = sin(psi_inv);
		aa[0][2] = 0.0;
		aa[1][0] = -sin(psi_inv);
		aa[1][1] = cos(psi_inv);
		aa[1][2] = 0.0;
		aa[2][0] = 0.0;
		aa[2][1] = 0.0;
		aa[2][2] = 1.0;

		bb[0][0] = cos(gama_inv);
		bb[0][1] = 0.0;
		bb[0][2] = sin(gama_inv);
		bb[1][0] = 0.0;
		bb[1][1] = 1.0; 
		bb[1][2] = 0.0;
		bb[2][0] = -sin(gama_inv);
		bb[2][1] = 0.0;
		bb[2][2] = cos(gama_inv);
		
		MatMult(3, 3, 3, *bb, *aa, *Binv2vel);
		MatMult(3, 3, 3, *bb, *aa, *Binv2vel);

		Transpose(3, 3, ECI2INV, INV2ECI);
		Transpose(3, 3, Binv2vel, vel2Binv);

		MatMult(3, 3, 3, *INV2ECI, *vel2Binv, *Tvel2eci);

		Vm = sqrt(SQR(Ve[0]) + SQR(Ve[1]) + SQR(Ve[2]));

		// Atmospheric Data
		CAtm us76;
		double Temp, sos, rho, machc, atm_out[3];
		us76.atm(Alt, atm_out);
		Temp  = atm_out[0];
		rho   = atm_out[2];
		sos   = 20.04680275919057*sqrt(Temp);
		machc = Vm / sos;

		//
		double Fxa, Fza, Faero_tb[3];
		Fxa = 0.5*rho*SQR(Vm) * Sref * Cx / mass_rv;
		Fza = 0.5*rho*SQR(Vm) * Sref * Cz / mass_rv;
		Faero_tb[0] = Fxa;
		Faero_tb[1] = 0;
		Faero_tb[2] = Fza;

		//
		double body2vel[3][3], bodyrol[3][3], dummy[3], Der[3], DD[3];
		body2vel[0][0] = cos(alp_t);
		body2vel[0][1] = 0.0;
		body2vel[0][2] = sin(alp_t);
		body2vel[1][0] = 0.0; 
		body2vel[1][1] = 1.0; 
		body2vel[1][2] = 0.0;
		body2vel[2][0] = -sin(alp_t);
		body2vel[2][1] = 0.0;
		body2vel[2][2] = cos(alp_t);

		bodyrol[0][0] = 1.0; 
		bodyrol[0][1] = 0.0;
		bodyrol[0][2] = 0.0; 
		bodyrol[1][0] = 0.0;
		bodyrol[1][1] = cos(phi_r);
		bodyrol[1][2] = -sin(phi_r);
		bodyrol[2][0] = 0.0; 
		bodyrol[2][1] = sin(phi_r);
		bodyrol[2][2] = cos(phi_r);

		MatMult(3, 1, 3, *body2vel, Faero_tb, dummy);
		MatMult(3, 1, 3, *bodyrol, dummy, Der);
		MatMult(3, 1, 3, *Tvel2eci, Der, DD);

		//
		dydx[0] = Vx;
		dydx[1] = Vy;
		dydx[2] = Vz;
		dydx[3] = gx + DD[0];
		dydx[4] = gy + DD[1];
		dydx[5] = gz + DD[2];

	}
};
