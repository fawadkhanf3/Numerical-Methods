#include "mex.h"
#include "conio.h"
#include "Utilities/nr3.h"
#include "Numerical Methods/odeint.h"
#include "Numerical Methods/stepper.h"
#include "Numerical Methods/stepperdopr5.h"
#include "Numerical Methods/stepperbs.h"
#include "Numerical Methods/stepperdopr853.h"
#include "Model/eoms3.h"

// Function Declaration
double **ode45cpp(double *tspan, double *xinit, int &length, int &successful, int &failed);

// Mex-Function (Gateway to MATLAB)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *Tspan, *y0; // Pointers to Input Arrays
	double *ty;         // Pointer to Output Arrays

	Tspan = mxGetPr(prhs[0]); // Input Array 1
	y0 = mxGetPr(prhs[1]);    // Input Array 2

	int length, successful_steps, failed_steps;
	double **yp = ode45cpp(Tspan, y0, length, successful_steps, failed_steps);

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

	printf("\n Successful Steps     = %d \n", successful_steps);
	printf("\n Failed Steps         = %d \n", failed_steps);
	printf("\n Function Evaluations = %d \n", (successful_steps + failed_steps) * 6);

}

// Function Definition
double **ode45cpp(double *tspan, double *xinit, int &length, int &successful, int &failed)
{
	//** Initializations **\\

	const int nvar = 6; // Number of States

	const double atol = 1.0e-5, rtol = atol, h1 = 0.01, hmin = 0.0,
		x1 = tspan[0], x2 = tspan[1]; // Time Span [x1 x2];

    // Initial Conditions
	VecDoub ystart(nvar);
	ystart[0] = xinit[0];
	ystart[1] = xinit[1];
	ystart[2] = xinit[2];
	ystart[3] = xinit[3];
	ystart[4] = xinit[4];
	ystart[5] = xinit[5];

	//** Sim Setup **\\

	Output out(-1); //OR// Output out(20); //(For Dense Output)

	// Equations of Motion (Differential Equations)
	eoms d;

	// Solvers \\

	//Odeint<StepperDopr5<eoms> > ode(ystart, x1, x2, atol, rtol, h1, hmin, out, d);
	//Odeint<StepperDopr853<eoms> > ode(ystart, x1, x2, atol, rtol, h1, hmin, out, d);
    Odeint<StepperBS<eoms> > ode(ystart, x1, x2, atol, rtol, h1, hmin, out, d);

	// Integration
	ode.integrate();

	// Outputs
	successful = ode.nok;
	failed     = ode.nbad;
	length     = out.count;

	int counts     = (out.count)*(nvar + 1);

	double **t_y = new double*[counts];

	for (int i = 0; i < out.count; i++)
	{
		t_y[i] = new double[7];
		t_y[i][0] = out.xsave[i];
		t_y[i][1] = out.ysave[0][i];
		t_y[i][2] = out.ysave[1][i];
		t_y[i][3] = out.ysave[2][i];
		t_y[i][4] = out.ysave[3][i];
		t_y[i][5] = out.ysave[4][i];
		t_y[i][6] = out.ysave[5][i];
	}

	return t_y;
}