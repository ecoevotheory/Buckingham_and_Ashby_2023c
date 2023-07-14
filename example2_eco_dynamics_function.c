/***********************************************************************************************************
 * [t,S,I,EQFLAG] = example2_eco_dynamics_function(t_max,a,c,alphabar,sigma,b,K,A,B,cb1,cb2,res,alpha,alphapower,eqtol,init_pop,strain_totalH,strain_totalP)
 ***********************************************************************************************************/

/* Compile in Matlab using mex example2_eco_dynamics_function.c */

#include <mex.h>
#include <math.h>


/***********************************
 * Constant parameter values
 ***********************************/
#define MAXSTEPS 1e6 /* Maximum number of steps for ODE solver */
#define INTERVAL 1e2 /* Check if the system is close to equilibrium */
#define EPS 1e-6 /* ODE solver tolerance */
#define TINY 1e-30 /* Constant value for solver */
#define TINY2 1e-30 /* Constant value for solver to avoid tolerance issues */
/* RK solver parameters */
#define b21 0.2
#define b31 3.0/40.0
#define b32 9.0/40.0
#define b41 0.3
#define b42 -0.9
#define b43 1.2
#define b51 -11.0/54.0
#define b52 2.5
#define b53 -70.0/27.0
#define b54 35.0/27.0
#define b61 1631.0/55296
#define b62 175.0/512.0
#define b63 575.0/13824.0
#define b64 44275.0/110592
#define b65 253.0/4096.0
#define c1 37.0/378.0
#define c3 250.0/621.0
#define c4 125.0/594.0
#define c6 512.0/1771.0
#define dc5 -277.00/14336

/*************************************
 * Define structure for model parameters
 *************************************/
struct PARAM{
    double t_max;
    double a;
    double c;
    double alphabar;
    double sigma; 
    double b;
    double K;
    double A;
    double B;
    double cb1;
    double cb2;
    double eqtol;
    int strain_totalH;
    int strain_totalP;
};

/*************************************
 * Function prototypes
 *************************************/
int my_rungkut (double *T, double *S_out, double *I_out, double *EQFLAG, double *init_pop, double *res, double *alpha, double *alphapower, struct PARAM *p);
void rkqs(double *S, double *I, double *DSDT, double *DIDT, double *h, double *hnext, double *S_SCALE, double *I_SCALE, double *res,double *alpha, double *alphapower, struct PARAM *p);
void rkck(double *S, double *I, double *DSDT, double *DIDT, double *Sout, double *Iout, double *Serr, double *Ierr, double h, double *res, double *alpha, double *alphapower, struct PARAM *p);
void dynamic(double *S, double *I, double *DSDT, double *DIDT, double *res, double *alpha, double *alphapower, struct PARAM *p);
double FMAX(double, double);
double FMIN(double, double);

/*************************************
 * Main function
 *************************************/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    double *T, *S, *I, *EQFLAG, *init_pop, *res, *alpha, *alphapower, *parameter;
    double *t_temp, *S_temp, *I_temp;
    int i, j, k, colLen, maxsteps;
    struct PARAM p;


   /* Allocate inputs */
    if(nrhs!=18){
       mexErrMsgTxt("Incorrect number of input arguments!\n");
    }

    else{
        parameter= mxGetPr(prhs[0]);
        p.t_max= *parameter;
        parameter= mxGetPr(prhs[1]);
        p.a= *parameter;
        parameter= mxGetPr(prhs[2]);
        p.c= *parameter;
        parameter= mxGetPr(prhs[3]);
        p.alphabar= *parameter;
        parameter= mxGetPr(prhs[4]);        
        p.sigma= *parameter;
        parameter= mxGetPr(prhs[5]);
        p.b= *parameter;
        parameter= mxGetPr(prhs[6]);
        p.K= *parameter;
	parameter= mxGetPr(prhs[7]);
        p.A= *parameter;
	parameter= mxGetPr(prhs[8]);
        p.B= *parameter;
	parameter= mxGetPr(prhs[9]);
        p.cb1= *parameter;
	parameter= mxGetPr(prhs[10]);
        p.cb2= *parameter;
        res= mxGetPr(prhs[11]);
        alpha= mxGetPr(prhs[12]);
	alphapower= mxGetPr(prhs[13]);
        parameter= mxGetPr(prhs[14]);
        p.eqtol= *parameter;
        init_pop= mxGetPr(prhs[15]);    
        parameter= mxGetPr(prhs[16]);
        p.strain_totalH= (int)*parameter;
        parameter= mxGetPr(prhs[17]);
        p.strain_totalP= (int)*parameter;
    }

    maxsteps = (int)MAXSTEPS;
    
    /* Allocate memory */
    t_temp = malloc(maxsteps*sizeof(double));
    S_temp = malloc(maxsteps*(p.strain_totalH)*(p.strain_totalP)*sizeof(double));
    I_temp = malloc(maxsteps*(p.strain_totalH)*(p.strain_totalP)*sizeof(double));
    
    /* Initialise this output */           
    plhs[3] = mxCreateDoubleMatrix(1, 1, mxREAL);
    EQFLAG = mxGetPr(plhs[3]);
    EQFLAG[0] = 0;


    /* Call ODE solver */
    colLen = my_rungkut(t_temp, S_temp, I_temp, EQFLAG, init_pop, res, alpha, alphapower, &p); 

    /* Create outputs */
    plhs[0] = mxCreateDoubleMatrix(colLen, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(colLen, p.strain_totalH*p.strain_totalP, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(colLen, p.strain_totalH*p.strain_totalP, mxREAL);
    T = mxGetPr(plhs[0]);
    S = mxGetPr(plhs[1]);
    I = mxGetPr(plhs[2]);

    /* Copy data to outputs */
    for (i=0;i<colLen;i++){
        T[i] = t_temp[i];
        for (j=0;j<p.strain_totalH;j++) {
	    for (k=0;k<p.strain_totalP;k++) {
		S[i+j*colLen+k*colLen*p.strain_totalH]=S_temp[i+j*maxsteps+k*maxsteps*p.strain_totalH];
		I[i+j*colLen+k*colLen*p.strain_totalH]=I_temp[i+j*maxsteps+k*maxsteps*p.strain_totalH];
	    }
        }
    }
    
    /* Free memory */
    free(t_temp);
    free(S_temp);
    free(I_temp); 
    
    return;

}

/*****************************************
 * ODE solver
 ****************************************/
int my_rungkut (double *T, double *S_out, double *I_out, double *EQFLAG, double *init_pop, double *res, double *alpha, double *alphapower, struct PARAM *p){
    
    double *S, *I, *DSDT, *DIDT, *S_SCALE, *I_SCALE;
    double *SMIN, *SMAX, *IMIN, *IMAX, hnext[1], h[1];
    double t, nextcheck;
    int i, j, exitflag, count, maxsteps;
    
    /* Allocate memory */
    S = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    I = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    DSDT = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    DIDT = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    S_SCALE = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    I_SCALE = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    SMIN = malloc(p->strain_totalH*p->strain_totalP*sizeof(double)); 
    SMAX = malloc(p->strain_totalH*p->strain_totalP*sizeof(double)); 
    IMIN = malloc(p->strain_totalH*p->strain_totalP*sizeof(double)); 
    IMAX = malloc(p->strain_totalH*p->strain_totalP*sizeof(double)); 
    
    /* Other parameters */
    exitflag = 1;
    count=0;
    h[0] = 1e-3;
    hnext[0] = 1e-3;
    t=0;
    nextcheck = INTERVAL;
    maxsteps = (int)MAXSTEPS;
    
    /* Initialise populations */
    for(i=0;i<p->strain_totalH;i++){
	for (j=0;j<p->strain_totalP;j++){
            S[i+j*p->strain_totalH] = init_pop[2*(i+j*p->strain_totalH)];
            I[i+j*p->strain_totalH] = init_pop[2*(i+j*p->strain_totalH)+1];
	}
    }
    
    /* Initialise equilibrium arrays */
    for(i=0;i<p->strain_totalH;i++){
	for (j=0;j<p->strain_totalP;j++){
            SMIN[i+j*p->strain_totalH] = S[i+j*p->strain_totalH];
            SMAX[i+j*p->strain_totalH] = S[i+j*p->strain_totalH];
            IMIN[i+j*p->strain_totalH] = I[i+j*p->strain_totalH];
            IMAX[i+j*p->strain_totalH] = I[i+j*p->strain_totalH];
	}
    }
    
    /* Update output */
    T[0]=t;
    for (i=0; i<p->strain_totalH; i++) {
	for (j=0; j<p->strain_totalP;j++){
            S_out[(i+j*p->strain_totalH)*maxsteps] = S[i+j*p->strain_totalH];
            I_out[(i+j*p->strain_totalH)*maxsteps] = I[i+j*p->strain_totalH];
	}
    }

    
    /* Main loop: */
    do{
        /* This ensures the final step lands us on the final time point */
        if(1.1*hnext[0]>(p->t_max-t)){
            hnext[0] = p->t_max-t;
            h[0] = p->t_max-t;
            t=p->t_max;
            exitflag=0;
        }
        else{
            h[0] = hnext[0];
            t+=h[0];
        }
        if(t>=p->t_max) {
            t=p->t_max;
            exitflag=0;
        }
        /* This is where the equations are first solved */

        dynamic(S, I, DSDT, DIDT, res, alpha, alphapower, p);

        
        /* Adjust the step size to maintain accuracy */
        for (i=0; i<p->strain_totalH; i++){
	    for (j=0; j<p->strain_totalP; j++){
                S[i+j*p->strain_totalH] = FMAX(S[i+j*p->strain_totalH],0);
                S_SCALE[i+j*p->strain_totalH]=fabs(S[i+j*p->strain_totalH])+fabs(DSDT[i+j*p->strain_totalH]*(*h))+TINY;
		I[i+j*p->strain_totalH] = FMAX(I[i+j*p->strain_totalH],0);
                I_SCALE[i+j*p->strain_totalH]=fabs(I[i+j*p->strain_totalH])+fabs(DIDT[i+j*p->strain_totalH]*(*h))+TINY;
	    }
        }
       

        /* RK solver & adaptive step-size */

        rkqs(S, I, DSDT, DIDT, h, hnext, S_SCALE, I_SCALE, res, alpha, alphapower, p); 

        /* Make sure nothing has gone negative */
        for (i=0; i<p->strain_totalH; i++){
	    for (j=0; j<p->strain_totalP; j++){
                S[i+j*p->strain_totalH] = FMAX(S[i+j*p->strain_totalH],0);
                I[i+j*p->strain_totalH] = FMAX(I[i+j*p->strain_totalH],0);
	    }
        }
        
        /* Update output */
        count++;
        T[count] = t;
        for (i=0; i<p->strain_totalH; i++) {
	    for (j=0; j<p->strain_totalP; j++) {
                S_out[count + (i+j*p->strain_totalH)*maxsteps] = S[i+j*p->strain_totalH];
                I_out[count + (i+j*p->strain_totalH)*maxsteps] = I[i+j*p->strain_totalH];
	    }
        }


        /* For equilibrium check */
        for (i=0; i<p->strain_totalH; i++){
	    for (j=0; j<p->strain_totalP; j++) {
                SMIN[i+j*p->strain_totalH] = FMIN(SMIN[i+j*p->strain_totalH],S[i+j*p->strain_totalH]);
                SMAX[i+j*p->strain_totalH] = FMAX(SMAX[i+j*p->strain_totalH],S[i+j*p->strain_totalH]);
                IMIN[i+j*p->strain_totalH] = FMIN(IMIN[i+j*p->strain_totalH],I[i+j*p->strain_totalH]);
                IMAX[i+j*p->strain_totalH] = FMAX(IMAX[i+j*p->strain_totalH],I[i+j*p->strain_totalH]);
	    }

        }


        /* Check if we're close to equilibrium */
        if(t>nextcheck){
            exitflag = 0;
            for (i=0; i<p->strain_totalH; i++){
		for (j=0; j<p->strain_totalP; j++){
                    if(fabs(SMAX[i+j*p->strain_totalH]-SMIN[i+j*p->strain_totalH])>p->eqtol || fabs(IMAX[i+j*p->strain_totalH]-IMIN[i+j*p->strain_totalH])>p->eqtol){
                        exitflag = 1;
                        break;
		    }
                }
            }
            /* If close to equilibrium, then break */
            if(exitflag==0){
                t=p->t_max; 
                T[count] = t; 
                EQFLAG[0] = 1; 
                break; 
            } 
            
            /* If not, then reset min/max values for each class */
            nextcheck+=INTERVAL;
            for (i=0; i<p->strain_totalH; i++){
		for (j=0; j<p->strain_totalP; j++){
                    SMIN[i+j*p->strain_totalH] = S[i+j*p->strain_totalH];
                    SMAX[i+j*p->strain_totalH] = S[i+j*p->strain_totalH];
                    IMIN[i+j*p->strain_totalH] = I[i+j*p->strain_totalH];
                    IMAX[i+j*p->strain_totalH] = I[i+j*p->strain_totalH];
		}
            }

        }
    }while(count<(maxsteps-1) && t<=p->t_max && exitflag);
    count++;

    
    
    /* Free memory */
    free(S);
    free(I);
    free(DSDT);
    free(DIDT);
    free(S_SCALE);
    free(I_SCALE);
    free(SMIN);
    free(IMIN);
    free(SMAX);
    free(IMAX);
    
    return count;
}

/***************************************
 * This generates the adaptive step-size
 **************************************/
void rkqs(double *S, double *I, double *DSDT, double *DIDT, double *h, double *hnext, double *S_SCALE, double *I_SCALE, double *res, double *alpha, double *alphapower, struct PARAM *p)
{
    double *S_temp, *I_temp, *S_err, *I_err;
    double htemp, errmax;
    int i, j, count;
    
    /* Allocate memory */
    S_temp = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    I_temp = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    S_err = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    I_err = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    
    count = 0;
    for(;;)
    {
        rkck(S, I, DSDT, DIDT, S_temp, I_temp, S_err, I_err, *h, res, alpha, alphapower, p); 
        
        errmax= 0.0;
        for(i=0;i<p->strain_totalH;i++){
	    for (j=0; j<p->strain_totalP; j++){

                errmax= FMAX(errmax, fabs(S_err[i+j*p->strain_totalH]/(S_SCALE[i+j*p->strain_totalH])));   
                errmax= FMAX(errmax, fabs(I_err[i+j*p->strain_totalH]/(I_SCALE[i+j*p->strain_totalH])));      
	    }
        }

        errmax/= EPS;

        if(errmax<=1.0) break;
        htemp= 0.9*(*h)*pow(errmax, -0.25);
        *h= (*h>=0.0 ? FMAX(htemp, 0.1*(*h)) : FMIN(htemp, 0.1*(*h)));
        count++;
            

        if(count>1e4){
            mexErrMsgTxt("stuck in loop!\n");
            break;
        }
    }    
    if(errmax > 1.89E-4) {
        *hnext= 0.9*(*h)*pow(errmax, -0.2);
    }
    else {
        *hnext= 5.0*(*h);
    }    
    *hnext = FMAX(*hnext, p->t_max/MAXSTEPS);
    
    for(i=0;i<p->strain_totalH;i++){
	for (j=0; j<p->strain_totalP; j++){
            S[i+j*p->strain_totalH] = S_temp[i+j*p->strain_totalH];
            I[i+j*p->strain_totalH] = I_temp[i+j*p->strain_totalH];
	}
    }
    
    /* Free memory */
    free(S_temp);
    free(I_temp);
    free(S_err);
    free(I_err);
}

/**************************************
 * Standard RK solver
 **************************************/
void rkck(double *S, double *I, double *DSDT, double *DIDT, double *Sout, double *Iout, double *Serr, double *Ierr, double h, double *res, double *alpha, double *alphapower, struct PARAM *p){
    int i, j;
    double *Sk1, *Sk2, *Sk3, *Sk4, *Sk5, *Sk6, *Stemp;
    double *Ik1, *Ik2, *Ik3, *Ik4, *Ik5, *Ik6, *Itemp;
    double dc1=c1-2825.0/27648.0, dc3=c3-18575.0/48384.0, dc4=c4-13525.0/55296.0, 
            dc6=c6-0.25;
    
    /* Allocate memory */    
    Sk1 = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    Sk2 = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    Sk3 = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    Sk4 = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    Sk5 = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    Sk6 = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    Stemp = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));

    Ik1 = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    Ik2 = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    Ik3 = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    Ik4 = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    Ik5 = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    Ik6 = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    Itemp = malloc(p->strain_totalH*p->strain_totalP*sizeof(double));
    
    for(i=0;i<p->strain_totalH;i++){
	for (j=0; j<p->strain_totalP; j++){
            Stemp[i+j*p->strain_totalH] = S[i+j*p->strain_totalH] + b21*h*DSDT[i+j*p->strain_totalH];
            Itemp[i+j*p->strain_totalH] = I[i+j*p->strain_totalH] + b21*h*DIDT[i+j*p->strain_totalH];
	}
    }

    dynamic(Stemp, Itemp, Sk2, Ik2, res, alpha, alphapower, p);
    
    for(i=0;i<p->strain_totalH;i++){
	for (j=0; j<p->strain_totalP; j++){
            Stemp[i+j*p->strain_totalH] = S[i+j*p->strain_totalH]+h*(b31*DSDT[i+j*p->strain_totalH]+b32*Sk2[i+j*p->strain_totalH]);
            Itemp[i+j*p->strain_totalH] = I[i+j*p->strain_totalH]+h*(b31*DIDT[i+j*p->strain_totalH]+b32*Ik2[i+j*p->strain_totalH]);
	}
    }    

    dynamic(Stemp, Itemp, Sk3, Ik3, res, alpha, alphapower, p);
    
    for(i=0;i<p->strain_totalH;i++){
	for (j=0; j<p->strain_totalP; j++){
            Stemp[i+j*p->strain_totalH] = S[i+j*p->strain_totalH]+h*(b41*DSDT[i+j*p->strain_totalH]+b42*Sk2[i+j*p->strain_totalH]+b43*Sk3[i+j*p->strain_totalH]);
            Itemp[i+j*p->strain_totalH] = I[i+j*p->strain_totalH]+h*(b41*DIDT[i+j*p->strain_totalH]+b42*Ik2[i+j*p->strain_totalH]+b43*Ik3[i+j*p->strain_totalH]);
	}
    }

    dynamic(Stemp, Itemp, Sk4, Ik4, res, alpha, alphapower, p);
    
    for(i=0;i<p->strain_totalH;i++){
	for (j=0; j<p->strain_totalP; j++){
            Stemp[i+j*p->strain_totalH] = S[i+j*p->strain_totalH]+h*(b51*DSDT[i+j*p->strain_totalH]+b52*Sk2[i+j*p->strain_totalH]+b53*Sk3[i+j*p->strain_totalH]+b54*Sk4[i+j*p->strain_totalH]);
            Itemp[i+j*p->strain_totalH] = I[i+j*p->strain_totalH]+h*(b51*DIDT[i+j*p->strain_totalH]+b52*Ik2[i+j*p->strain_totalH]+b53*Ik3[i+j*p->strain_totalH]+b54*Ik4[i+j*p->strain_totalH]);
	}
    }

    dynamic(Stemp, Itemp, Sk5, Ik5, res, alpha, alphapower, p); 
    
    for(i=0;i<p->strain_totalH;i++){
	for (j=0; j<p->strain_totalP; j++){ 
            Stemp[i+j*p->strain_totalH] = S[i+j*p->strain_totalH]+h*(b61*DSDT[i+j*p->strain_totalH]+b62*Sk2[i+j*p->strain_totalH]+b63*Sk3[i+j*p->strain_totalH]+b64*Sk4[i+j*p->strain_totalH]+b65*Sk5[i+j*p->strain_totalH]);
            Itemp[i+j*p->strain_totalH] = I[i+j*p->strain_totalH]+h*(b61*DIDT[i+j*p->strain_totalH]+b62*Ik2[i+j*p->strain_totalH]+b63*Ik3[i+j*p->strain_totalH]+b64*Ik4[i+j*p->strain_totalH]+b65*Ik5[i+j*p->strain_totalH]);
	}
    }

    dynamic(Stemp, Itemp, Sk6, Ik6, res, alpha, alphapower, p);
    
    for(i=0;i<p->strain_totalH;i++){
	for (j=0; j<p->strain_totalP; j++){
            Sout[i+j*p->strain_totalH]= S[i+j*p->strain_totalH]+h*(c1*DSDT[i+j*p->strain_totalH]+c3*Sk3[i+j*p->strain_totalH]+c4*Sk4[i+j*p->strain_totalH]+c6*Sk6[i+j*p->strain_totalH]);
            Serr[i+j*p->strain_totalH]= h*(dc1*DSDT[i+j*p->strain_totalH]+dc3*Sk3[i+j*p->strain_totalH]+dc4*Sk4[i+j*p->strain_totalH]+dc5*Sk5[i+j*p->strain_totalH]+dc6*Sk6[i+j*p->strain_totalH]);
	    Iout[i+j*p->strain_totalH]= I[i+j*p->strain_totalH]+h*(c1*DIDT[i+j*p->strain_totalH]+c3*Ik3[i+j*p->strain_totalH]+c4*Ik4[i+j*p->strain_totalH]+c6*Ik6[i+j*p->strain_totalH]);
            Ierr[i+j*p->strain_totalH]= h*(dc1*DIDT[i+j*p->strain_totalH]+dc3*Ik3[i+j*p->strain_totalH]+dc4*Ik4[i+j*p->strain_totalH]+dc5*Ik5[i+j*p->strain_totalH]+dc6*Ik6[i+j*p->strain_totalH]);
	}
    }

    /* Free memory */
    free(Sk1);
    free(Sk2);
    free(Sk3);
    free(Sk4);
    free(Sk5);
    free(Sk6);
    free(Stemp);
    free(Ik1);
    free(Ik2);
    free(Ik3);
    free(Ik4);
    free(Ik5);
    free(Ik6);
    free(Itemp);
}

/**************************************
 * Population and evolutionary dynamics
 **************************************/
void dynamic(double *S, double *I, double *DSDT, double *DIDT, double *res, double *alpha, double *alphapower, struct PARAM *p){
    

    int i, j;

    double *Iallparasites, *totalparasites, N, betaallinfecteds;

    Iallparasites = malloc(p->strain_totalH*sizeof(double));
    totalparasites = malloc(p->strain_totalP*sizeof(double));

    /* Population sums */
    N = 0;
    for(i=0;i<p->strain_totalH;i++){ 
        for (j=0; j<p->strain_totalP; j++){
	    N = N + I[i+j*p->strain_totalH] + S[i+j*p->strain_totalH];
	}
    }

    betaallinfecteds = 0;
    for(i=0;i<p->strain_totalH;i++){ 
        for (j=0; j<p->strain_totalP; j++){
	    betaallinfecteds = betaallinfecteds + ((p->c*alpha[j])/(p->a+alpha[j]))*(1-p->K*alphapower[j])*I[i+j*p->strain_totalH];
	}
    }

    for(i=0;i<p->strain_totalH;i++){ 
	    Iallparasites[i]=0;
    }
    for(i=0;i<p->strain_totalH;i++){ 
        for (j=0; j<p->strain_totalP; j++){
	    Iallparasites[i]=Iallparasites[i]+I[i+j*p->strain_totalH];
	}
    }

    for(j=0;j<p->strain_totalP;j++){ 
	    totalparasites[j]=0;
    }
    for(j=0;j<p->strain_totalP;j++){ 
        for (i=0; i<p->strain_totalH; i++){
	    totalparasites[j]=totalparasites[j]+I[i+j*p->strain_totalH];
	}
    }

    for(i=0;i<p->strain_totalH;i++){
	for (j=0; j<p->strain_totalP; j++){
	    DSDT[i+j*p->strain_totalH] = 0;
	    DIDT[i+j*p->strain_totalH] = (1-res[i])*((p->c*alpha[j])/(p->a+alpha[j]))*(1-p->K*alphapower[j])*totalparasites[j]*S[i]-(p->A+p->B*N)*I[i+j*p->strain_totalH]-alpha[j]*I[i+j*p->strain_totalH];
	}

 	DSDT[i] = (1-((p->cb1*(1-exp(p->cb2*res[i])))/(1-exp(p->cb2))))*p->b*(S[i]+Iallparasites[i])-(1-res[i])*betaallinfecteds*S[i]-(p->A+p->B*N)*S[i];
	   
    }


}


/***************************************
 * Return maximum of two inputs
 ***************************************/
double FMAX(double l, double r)
{
    if(l>r)return l;
    else   return r;
}

/***************************************
 * Return minimum of two inputs
 ***************************************/
double FMIN(double l, double r)
{
    if(l<r)return l;
    else   return r;
}
