/***********************************************************************************************************
 * [t,SJ,SA,IJ,IA,EQFLAG] = example1_eco_dynamics_function(t_max,beta0max,beta0min,p0,p1,h0,h1,a0,a1,alpha,b,q,epsilon,resJ_current,resA_current,eqtol,init_pop,strain_totalJ,strain_totalA)
 ***********************************************************************************************************/

/* Compile in Matlab using mex example1_eco_dynamics_function.c */

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
    double beta0max;
    double beta0min;
    double p0;
    double p1; 
    double h0;
    double h1;
    double a0;
    double a1;
    double alpha;
    double b;
    double q;
    double epsilon;
    double eqtol;
    int strain_totalJ;
    int strain_totalA;
};

/*************************************
 * Function prototypes
 *************************************/
int my_rungkut (double *T, double *SJ_out, double *SA_out, double *IJ_out, double *IA_out, double *EQFLAG, double *init_pop, double *betaJ, double *betaA, struct PARAM *p);
void rkqs(double *SJ, double *SA, double *IJ, double *IA,  double *DSJDT, double *DSADT, double *DIJDT, double *DIADT,  double *h, double *hnext, double *SJ_SCALE, double *SA_SCALE, double *IJ_SCALE, double *IA_SCALE, double *betaJ,double *betaA, struct PARAM *p);
void rkck(double *SJ, double *SA,  double *IJ, double *IA,  double *DSJDT, double *DSADT,  double *DIJDT, double *DIADT, double *SJout, double *SAout,  double *IJout, double *IAout, double *SJerr, double *SAerr,  double *IJerr, double *IAerr,  double h, double *betaJ, double *betaA, struct PARAM *p);
void dynamic(double *SJ, double *SA, double *IJ, double *IA, double *DSJDT, double *DSADT, double *DIJDT, double *DIADT, double *betaJ, double *betaA, struct PARAM *p);
double FMAX(double, double);
double FMIN(double, double);

/*************************************
 * Main function
 *************************************/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    double *T, *SJ, *SA, *IJ, *IA, *EQFLAG, *init_pop, *betaJ, *betaA, *parameter;
    double *t_temp, *SJ_temp, *SA_temp, *IJ_temp, *IA_temp;
    int i, j, k, colLen, maxsteps;
    struct PARAM p;


   /* Allocate inputs */
    if(nrhs!=19){
       mexErrMsgTxt("Incorrect number of input arguments!\n");
    }

    else{
        parameter= mxGetPr(prhs[0]);
        p.t_max= *parameter;
        parameter= mxGetPr(prhs[1]);
        p.beta0max= *parameter;
        parameter= mxGetPr(prhs[2]);
        p.beta0min= *parameter;
        parameter= mxGetPr(prhs[3]);
        p.p0= *parameter;
        parameter= mxGetPr(prhs[4]);        
        p.p1= *parameter;
        parameter= mxGetPr(prhs[5]);
        p.h0= *parameter;
        parameter= mxGetPr(prhs[6]);
        p.h1= *parameter;
        parameter= mxGetPr(prhs[7]);
        p.a0= *parameter;
	parameter= mxGetPr(prhs[8]);
        p.a1= *parameter;
	parameter= mxGetPr(prhs[9]);
        p.alpha= *parameter;
	parameter= mxGetPr(prhs[10]);
        p.b= *parameter;
	parameter= mxGetPr(prhs[11]);
        p.q= *parameter;
	parameter= mxGetPr(prhs[12]);
        p.epsilon= *parameter;
        betaJ= mxGetPr(prhs[13]);
        betaA= mxGetPr(prhs[14]);
        parameter= mxGetPr(prhs[15]);
        p.eqtol= *parameter;
        init_pop= mxGetPr(prhs[16]);    
        parameter= mxGetPr(prhs[17]);
        p.strain_totalJ= (int)*parameter;
        parameter= mxGetPr(prhs[18]);
        p.strain_totalA= (int)*parameter;
    }

    maxsteps = (int)MAXSTEPS;
    
    /* Allocate memory */
    t_temp = malloc(maxsteps*sizeof(double));
    SJ_temp = malloc(maxsteps*(p.strain_totalJ)*(p.strain_totalA)*sizeof(double));
    IJ_temp = malloc(maxsteps*(p.strain_totalJ)*(p.strain_totalA)*sizeof(double));
    SA_temp = malloc(maxsteps*(p.strain_totalJ)*(p.strain_totalA)*sizeof(double));
    IA_temp = malloc(maxsteps*(p.strain_totalJ)*(p.strain_totalA)*sizeof(double));
    
    /* Initialise this output */           
    plhs[5] = mxCreateDoubleMatrix(1, 1, mxREAL);
    EQFLAG = mxGetPr(plhs[5]);
    EQFLAG[0] = 0;

    /* Call ODE solver */
    colLen = my_rungkut(t_temp, SJ_temp, SA_temp, IJ_temp, IA_temp,EQFLAG, init_pop, betaJ,betaA, &p); 

    /* Create outputs */
    plhs[0] = mxCreateDoubleMatrix(colLen, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(colLen, p.strain_totalJ*p.strain_totalA, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(colLen, p.strain_totalJ*p.strain_totalA, mxREAL);
    plhs[3] = mxCreateDoubleMatrix(colLen, p.strain_totalJ*p.strain_totalA, mxREAL);
    plhs[4] = mxCreateDoubleMatrix(colLen, p.strain_totalJ*p.strain_totalA, mxREAL);
    T = mxGetPr(plhs[0]);
    SJ = mxGetPr(plhs[1]);
    SA = mxGetPr(plhs[2]);
    IJ = mxGetPr(plhs[3]);
    IA = mxGetPr(plhs[4]);

    /* Copy data to outputs */
    for (i=0;i<colLen;i++){
        T[i] = t_temp[i];
        for (j=0;j<p.strain_totalJ;j++) {
	    for (k=0;k<p.strain_totalA;k++) {
		SJ[i+j*colLen+k*colLen*p.strain_totalJ]=SJ_temp[i+j*maxsteps+k*maxsteps*p.strain_totalJ];
		SA[i+j*colLen+k*colLen*p.strain_totalJ]=SA_temp[i+j*maxsteps+k*maxsteps*p.strain_totalJ];
		IJ[i+j*colLen+k*colLen*p.strain_totalJ]=IJ_temp[i+j*maxsteps+k*maxsteps*p.strain_totalJ];
		IA[i+j*colLen+k*colLen*p.strain_totalJ]=IA_temp[i+j*maxsteps+k*maxsteps*p.strain_totalJ];
	    }
        }
    }
    
    /* Free memory */
    free(t_temp);
    free(SJ_temp);
    free(IJ_temp); 
    free(SA_temp);
    free(IA_temp); 
    
    return;

}

/*****************************************
 * ODE solver
 ****************************************/
int my_rungkut (double *T, double *SJ_out, double *SA_out, double *IJ_out, double *IA_out,  double *EQFLAG, double *init_pop, double *betaJ, double *betaA, struct PARAM *p){
    
    double *SJ, *SA,  *IJ, *IA, *DSJDT, *DSADT, *DIJDT, *DIADT,  *SJ_SCALE, *SA_SCALE, *IJ_SCALE, *IA_SCALE;
    double *SJMIN, *SJMAX, *SAMIN, *SAMAX, *IJMIN, *IJMAX, *IAMIN, *IAMAX, hnext[1], h[1];
    double t, nextcheck;
    int i, j, exitflag, count, maxsteps;
    
    /* Allocate memory */
    SJ = malloc(p->strain_totalJ*p->strain_totalA*sizeof(double));
    SA = malloc(p->strain_totalJ*p->strain_totalA*sizeof(double));
    IJ = malloc(p->strain_totalJ*p->strain_totalA*sizeof(double));
    IA = malloc(p->strain_totalJ*p->strain_totalA*sizeof(double));
    DSJDT = malloc(p->strain_totalJ*p->strain_totalA*sizeof(double));
    DSADT = malloc(p->strain_totalJ*p->strain_totalA*sizeof(double));
    DIJDT = malloc(p->strain_totalJ*p->strain_totalA*sizeof(double));
    DIADT = malloc(p->strain_totalJ*p->strain_totalA*sizeof(double));
    SJ_SCALE = malloc(p->strain_totalJ*p->strain_totalA*sizeof(double));
    SA_SCALE = malloc(p->strain_totalJ*p->strain_totalA*sizeof(double));
    IJ_SCALE = malloc(p->strain_totalJ*p->strain_totalA*sizeof(double));
    IA_SCALE = malloc(p->strain_totalJ*p->strain_totalA*sizeof(double));
    SJMIN = malloc(p->strain_totalJ*p->strain_totalA*sizeof(double)); 
    SJMAX = malloc(p->strain_totalJ*p->strain_totalA*sizeof(double)); 
    SAMIN = malloc(p->strain_totalJ*p->strain_totalA*sizeof(double)); 
    SAMAX = malloc(p->strain_totalJ*p->strain_totalA*sizeof(double)); 
    IJMIN = malloc(p->strain_totalJ*p->strain_totalA*sizeof(double)); 
    IJMAX = malloc(p->strain_totalJ*p->strain_totalA*sizeof(double)); 
    IAMIN = malloc(p->strain_totalJ*p->strain_totalA*sizeof(double)); 
    IAMAX = malloc(p->strain_totalJ*p->strain_totalA*sizeof(double)); 
    
    /* Other parameters */
    exitflag = 1;
    count=0;
    /* k=1;*/
    h[0] = 1e-3;
    hnext[0] = 1e-3;
    t=0;
    nextcheck = INTERVAL;
    maxsteps = (int)MAXSTEPS;
    
    /* Initialise populations */
    for(i=0;i<p->strain_totalJ;i++){
	for (j=0;j<p->strain_totalA;j++){
            SJ[i+j*p->strain_totalJ] = init_pop[4*(i+j*p->strain_totalJ)];
            SA[i+j*p->strain_totalJ] = init_pop[4*(i+j*p->strain_totalJ)+1];
            IJ[i+j*p->strain_totalJ] = init_pop[4*(i+j*p->strain_totalJ)+2];
            IA[i+j*p->strain_totalJ] = init_pop[4*(i+j*p->strain_totalJ)+3];
	}
    }
    
    /* Initialise equilibrium arrays */
    for(i=0;i<p->strain_totalJ;i++){
	for (j=0;j<p->strain_totalA;j++){
            SJMIN[i+j*p->strain_totalJ] = SJ[i+j*p->strain_totalJ];
            SJMAX[i+j*p->strain_totalJ] = SJ[i+j*p->strain_totalJ];
            SAMIN[i+j*p->strain_totalJ] = SA[i+j*p->strain_totalJ];
            SAMAX[i+j*p->strain_totalJ] = SA[i+j*p->strain_totalJ];
            IJMIN[i+j*p->strain_totalJ] = IJ[i+j*p->strain_totalJ];
            IJMAX[i+j*p->strain_totalJ] = IJ[i+j*p->strain_totalJ];
            IAMIN[i+j*p->strain_totalJ] = IA[i+j*p->strain_totalJ];
            IAMAX[i+j*p->strain_totalJ] = IA[i+j*p->strain_totalJ];
	}
    }
    
    /* Update output */
    T[0]=t;
    for (i=0; i<p->strain_totalJ; i++) {
	for (j=0; j<p->strain_totalA;j++){
            SJ_out[(i+j*p->strain_totalJ)*maxsteps] = SJ[i+j*p->strain_totalJ];
            SA_out[(i+j*p->strain_totalJ)*maxsteps] = SA[i+j*p->strain_totalJ];
            IJ_out[(i+j*p->strain_totalJ)*maxsteps] = IJ[i+j*p->strain_totalJ];
            IA_out[(i+j*p->strain_totalJ)*maxsteps] = IA[i+j*p->strain_totalJ];
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

        dynamic(SJ, SA,  IJ, IA,  DSJDT, DSADT,  DIJDT, DIADT, betaJ,betaA, p);

        
        /* Adjust the step size to maintain accuracy */
        for (i=0; i<p->strain_totalJ; i++){
	    for (j=0; j<p->strain_totalA; j++){
                SJ[i+j*p->strain_totalJ] = FMAX(SJ[i+j*p->strain_totalJ],0);
                SJ_SCALE[i+j*p->strain_totalJ]=fabs(SJ[i+j*p->strain_totalJ])+fabs(DSJDT[i+j*p->strain_totalJ]*(*h))+TINY;
                SA[i+j*p->strain_totalJ] = FMAX(SA[i+j*p->strain_totalJ],0);
                SA_SCALE[i+j*p->strain_totalJ]=fabs(SA[i+j*p->strain_totalJ])+fabs(DSADT[i+j*p->strain_totalJ]*(*h))+TINY;
                IJ[i+j*p->strain_totalJ] = FMAX(IJ[i+j*p->strain_totalJ],0);
                IJ_SCALE[i+j*p->strain_totalJ]=fabs(IJ[i+j*p->strain_totalJ])+fabs(DIJDT[i+j*p->strain_totalJ]*(*h))+TINY;
                IA[i+j*p->strain_totalJ] = FMAX(IA[i+j*p->strain_totalJ],0);
                IA_SCALE[i+j*p->strain_totalJ]=fabs(IA[i+j*p->strain_totalJ])+fabs(DIADT[i+j*p->strain_totalJ]*(*h))+TINY;
	    }
        }
        
        /* RK solver & adaptive step-size */

        rkqs(SJ, SA, IJ, IA, DSJDT, DSADT, DIJDT, DIADT, h, hnext, SJ_SCALE, SA_SCALE, IJ_SCALE, IA_SCALE, betaJ,betaA, p);

        /* Make sure nothin has gone negative */
        for (i=0; i<p->strain_totalJ; i++){
	    for (j=0; j<p->strain_totalA; j++){
                SJ[i+j*p->strain_totalJ] = FMAX(SJ[i+j*p->strain_totalJ],0);
                SA[i+j*p->strain_totalJ] = FMAX(SA[i+j*p->strain_totalJ],0);
                IJ[i+j*p->strain_totalJ] = FMAX(IJ[i+j*p->strain_totalJ],0);
                IA[i+j*p->strain_totalJ] = FMAX(IA[i+j*p->strain_totalJ],0);
	    }
        }
        
        /* Update output */
        count++;
        T[count] = t;
        for (i=0; i<p->strain_totalJ; i++) {
	    for (j=0; j<p->strain_totalA; j++) {
                SJ_out[count + (i+j*p->strain_totalJ)*maxsteps] = SJ[i+j*p->strain_totalJ];
                SA_out[count + (i+j*p->strain_totalJ)*maxsteps] = SA[i+j*p->strain_totalJ];
                IJ_out[count + (i+j*p->strain_totalJ)*maxsteps] = IJ[i+j*p->strain_totalJ];
                IA_out[count + (i+j*p->strain_totalJ)*maxsteps] = IA[i+j*p->strain_totalJ];
	    }
        }


        /* For equilibrium check */
        for (i=0; i<p->strain_totalJ; i++){
	    for (j=0; j<p->strain_totalA; j++) {
                SJMIN[i+j*p->strain_totalJ] = FMIN(SJMIN[i+j*p->strain_totalJ],SJ[i+j*p->strain_totalJ]);
                SJMAX[i+j*p->strain_totalJ] = FMAX(SJMAX[i+j*p->strain_totalJ],SJ[i+j*p->strain_totalJ]);
                SAMIN[i+j*p->strain_totalJ] = FMIN(SAMIN[i+j*p->strain_totalJ],SA[i+j*p->strain_totalJ]);
                SAMAX[i+j*p->strain_totalJ] = FMAX(SAMAX[i+j*p->strain_totalJ],SA[i+j*p->strain_totalJ]);
                IJMIN[i+j*p->strain_totalJ] = FMIN(IJMIN[i+j*p->strain_totalJ],IJ[i+j*p->strain_totalJ]);
                IJMAX[i+j*p->strain_totalJ] = FMAX(IJMAX[i+j*p->strain_totalJ],IJ[i+j*p->strain_totalJ]);
                IAMIN[i+j*p->strain_totalJ] = FMIN(IAMIN[i+j*p->strain_totalJ],IA[i+j*p->strain_totalJ]);
                IAMAX[i+j*p->strain_totalJ] = FMAX(IAMAX[i+j*p->strain_totalJ],IA[i+j*p->strain_totalJ]);
	    }

        }


        /* Check if we're close to equilibrium */
        if(t>nextcheck){
            exitflag = 0;
            for (i=0; i<p->strain_totalJ; i++){
		for (j=0; j<p->strain_totalA; j++){
                    if(fabs(SJMAX[i+j*p->strain_totalJ]-SJMIN[i+j*p->strain_totalJ])>p->eqtol || fabs(SAMAX[i+j*p->strain_totalJ]-SAMIN[i+j*p->strain_totalJ])>p->eqtol || fabs(IJMAX[i+j*p->strain_totalJ]-IJMIN[i+j*p->strain_totalJ])>p->eqtol || fabs(IAMAX[i+j*p->strain_totalJ]-IAMIN[i+j*p->strain_totalJ])>p->eqtol ){
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
            for (i=0; i<p->strain_totalJ; i++){
		for (j=0; j<p->strain_totalA; j++){
                    SJMIN[i+j*p->strain_totalJ] = SJ[i+j*p->strain_totalJ];
                    SJMAX[i+j*p->strain_totalJ] = SJ[i+j*p->strain_totalJ];
                    SAMIN[i+j*p->strain_totalJ] = SA[i+j*p->strain_totalJ];
                    SAMAX[i+j*p->strain_totalJ] = SA[i+j*p->strain_totalJ];
                    IJMIN[i+j*p->strain_totalJ] = IJ[i+j*p->strain_totalJ];
                    IJMAX[i+j*p->strain_totalJ] = IJ[i+j*p->strain_totalJ];
                    IAMIN[i+j*p->strain_totalJ] = IA[i+j*p->strain_totalJ];
                    IAMAX[i+j*p->strain_totalJ] = IA[i+j*p->strain_totalJ];
		}
            }

        }
    }while(count<(maxsteps-1) && t<=p->t_max && exitflag);
    count++;

    
    
    /* Free memory */
    free(SJ);
    free(SA);
    free(IJ);
    free(IA);
    free(DSJDT);
    free(DSADT);
    free(DIJDT);
    free(DIADT);
    free(SJ_SCALE);
    free(SA_SCALE);
    free(IJ_SCALE);
    free(IA_SCALE);
    free(SJMIN);
    free(SAMIN);
    free(IJMIN);
    free(IAMIN);
    free(SJMAX);
    free(SAMAX);
    free(IJMAX);
    free(IAMAX);
    
    return count;
}

/***************************************
 * This generates the adaptive step-size
 **************************************/
void rkqs(double *SJ, double *SA, double *IJ, double *IA,   double *DSJDT, double *DSADT,  double *DIJDT, double *DIADT,  double *h, double *hnext, double *SJ_SCALE, double *SA_SCALE, double *IJ_SCALE, double *IA_SCALE, double *betaJ, double *betaA, struct PARAM *p)
{
    double *SJ_temp, *SA_temp,  *IJ_temp, *IA_temp, *SJ_err, *SA_err, *IJ_err, *IA_err;
    double htemp, errmax;
    int i, j, count;
    
    /* Allocate memory */
    SJ_temp = malloc(p->strain_totalJ*p->strain_totalA*sizeof(double));
    SA_temp = malloc(p->strain_totalJ*p->strain_totalA*sizeof(double));
    IJ_temp = malloc(p->strain_totalJ*p->strain_totalA*sizeof(double));
    IA_temp = malloc(p->strain_totalJ*p->strain_totalA*sizeof(double));
    SJ_err = malloc(p->strain_totalJ*p->strain_totalA*sizeof(double));
    SA_err = malloc(p->strain_totalJ*p->strain_totalA*sizeof(double));
    IJ_err = malloc(p->strain_totalJ*p->strain_totalA*sizeof(double));
    IA_err = malloc(p->strain_totalJ*p->strain_totalA*sizeof(double));
    
    count = 0;
    for(;;)
    {
        rkck(SJ, SA, IJ, IA, DSJDT, DSADT, DIJDT, DIADT, SJ_temp, SA_temp,  IJ_temp, IA_temp, SJ_err, SA_err,  IJ_err, IA_err,*h, betaJ,betaA, p); 
        
        errmax= 0.0;
        for(i=0;i<p->strain_totalJ;i++){
	    for (j=0; j<p->strain_totalA; j++){

                errmax= FMAX(errmax, fabs(SJ_err[i+j*p->strain_totalJ]/(SJ_SCALE[i+j*p->strain_totalJ]))); 
                errmax= FMAX(errmax, fabs(SA_err[i+j*p->strain_totalJ]/(SA_SCALE[i+j*p->strain_totalJ])));   
                errmax= FMAX(errmax, fabs(IJ_err[i+j*p->strain_totalJ]/(IJ_SCALE[i+j*p->strain_totalJ])));   
                errmax= FMAX(errmax, fabs(IA_err[i+j*p->strain_totalJ]/(IA_SCALE[i+j*p->strain_totalJ])));     
	    }
        }

        errmax/= EPS;

        if(errmax<=1.0) break;
        htemp= 0.9*(*h)*pow(errmax, -0.25);
        *h= (*h>=0.0 ? FMAX(htemp, 0.1*(*h)) : FMIN(htemp, 0.1*(*h)));
        count++;
            

        if(count>1e4){
            printf("%f\n",errmax);
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
    
    for(i=0;i<p->strain_totalJ;i++){
	for (j=0; j<p->strain_totalA; j++){
            SJ[i+j*p->strain_totalJ] = SJ_temp[i+j*p->strain_totalJ];
            SA[i+j*p->strain_totalJ] = SA_temp[i+j*p->strain_totalJ];
            IJ[i+j*p->strain_totalJ] = IJ_temp[i+j*p->strain_totalJ];
            IA[i+j*p->strain_totalJ] = IA_temp[i+j*p->strain_totalJ];
	}
    }
    
    /* Free memory */
    free(SJ_temp);
    free(IJ_temp);
    free(SJ_err);
    free(IJ_err);
    free(SA_temp);
    free(IA_temp);
    free(SA_err);
    free(IA_err);
}

/**************************************
 * Standard RK solver
 **************************************/
void rkck(double *SJ, double *SA, double *IJ, double *IA, double *DSJDT, double *DSADT, double *DIJDT, double *DIADT, double *SJout, double *SAout, double *IJout, double *IAout, double *SJerr, double *SAerr, double *IJerr, double *IAerr, double h, double *betaJ, double *betaA, struct PARAM *p){
    int i, j;
    double *SJk1, *SJk2, *SJk3, *SJk4, *SJk5, *SJk6, *SJtemp;
    double *IJk1, *IJk2, *IJk3, *IJk4, *IJk5, *IJk6, *IJtemp;
    double *SAk1, *SAk2, *SAk3, *SAk4, *SAk5, *SAk6, *SAtemp;
    double *IAk1, *IAk2, *IAk3, *IAk4, *IAk5, *IAk6, *IAtemp;
    double dc1=c1-2825.0/27648.0, dc3=c3-18575.0/48384.0, dc4=c4-13525.0/55296.0, 
            dc6=c6-0.25;
    
    /* Allocate memory */    
    SJk1 = malloc(p->strain_totalJ*p->strain_totalA*sizeof(double));
    SJk2 = malloc(p->strain_totalJ*p->strain_totalA*sizeof(double));
    SJk3 = malloc(p->strain_totalJ*p->strain_totalA*sizeof(double));
    SJk4 = malloc(p->strain_totalJ*p->strain_totalA*sizeof(double));
    SJk5 = malloc(p->strain_totalJ*p->strain_totalA*sizeof(double));
    SJk6 = malloc(p->strain_totalJ*p->strain_totalA*sizeof(double));
    SJtemp = malloc(p->strain_totalJ*p->strain_totalA*sizeof(double));
    
    SAk1 = malloc(p->strain_totalJ*p->strain_totalA*sizeof(double));
    SAk2 = malloc(p->strain_totalJ*p->strain_totalA*sizeof(double));
    SAk3 = malloc(p->strain_totalJ*p->strain_totalA*sizeof(double));
    SAk4 = malloc(p->strain_totalJ*p->strain_totalA*sizeof(double));
    SAk5 = malloc(p->strain_totalJ*p->strain_totalA*sizeof(double));
    SAk6 = malloc(p->strain_totalJ*p->strain_totalA*sizeof(double));
    SAtemp = malloc(p->strain_totalJ*p->strain_totalA*sizeof(double));

    IJk1 = malloc(p->strain_totalJ*p->strain_totalA*sizeof(double));
    IJk2 = malloc(p->strain_totalJ*p->strain_totalA*sizeof(double));
    IJk3 = malloc(p->strain_totalJ*p->strain_totalA*sizeof(double));
    IJk4 = malloc(p->strain_totalJ*p->strain_totalA*sizeof(double));
    IJk5 = malloc(p->strain_totalJ*p->strain_totalA*sizeof(double));
    IJk6 = malloc(p->strain_totalJ*p->strain_totalA*sizeof(double));
    IJtemp = malloc(p->strain_totalJ*p->strain_totalA*sizeof(double));

    IAk1 = malloc(p->strain_totalJ*p->strain_totalA*sizeof(double));
    IAk2 = malloc(p->strain_totalJ*p->strain_totalA*sizeof(double));
    IAk3 = malloc(p->strain_totalJ*p->strain_totalA*sizeof(double));
    IAk4 = malloc(p->strain_totalJ*p->strain_totalA*sizeof(double));
    IAk5 = malloc(p->strain_totalJ*p->strain_totalA*sizeof(double));
    IAk6 = malloc(p->strain_totalJ*p->strain_totalA*sizeof(double));
    IAtemp = malloc(p->strain_totalJ*p->strain_totalA*sizeof(double));
    
    for(i=0;i<p->strain_totalJ;i++){
	for (j=0; j<p->strain_totalA; j++){
            SJtemp[i+j*p->strain_totalJ] = SJ[i+j*p->strain_totalJ] + b21*h*DSJDT[i+j*p->strain_totalJ];
            SAtemp[i+j*p->strain_totalJ] = SA[i+j*p->strain_totalJ] + b21*h*DSADT[i+j*p->strain_totalJ];
            IJtemp[i+j*p->strain_totalJ] = IJ[i+j*p->strain_totalJ] + b21*h*DIJDT[i+j*p->strain_totalJ];
            IAtemp[i+j*p->strain_totalJ] = IA[i+j*p->strain_totalJ] + b21*h*DIADT[i+j*p->strain_totalJ];
	}
    }
    dynamic(SJtemp, SAtemp, IJtemp, IAtemp, SJk2, SAk2, IJk2, IAk2,  betaJ, betaA, p);


    
    for(i=0;i<p->strain_totalJ;i++){
	for (j=0; j<p->strain_totalA; j++){
            SJtemp[i+j*p->strain_totalJ] = SJ[i+j*p->strain_totalJ]+h*(b31*DSJDT[i+j*p->strain_totalJ]+b32*SJk2[i+j*p->strain_totalJ]);
	    SAtemp[i+j*p->strain_totalJ] = SA[i+j*p->strain_totalJ]+h*(b31*DSADT[i+j*p->strain_totalJ]+b32*SAk2[i+j*p->strain_totalJ]);
            IJtemp[i+j*p->strain_totalJ] = IJ[i+j*p->strain_totalJ]+h*(b31*DIJDT[i+j*p->strain_totalJ]+b32*IJk2[i+j*p->strain_totalJ]);
            IAtemp[i+j*p->strain_totalJ] = IA[i+j*p->strain_totalJ]+h*(b31*DIADT[i+j*p->strain_totalJ]+b32*IAk2[i+j*p->strain_totalJ]);
	}
    }    
    dynamic(SJtemp, SAtemp, IJtemp, IAtemp, SJk3, SAk3, IJk3, IAk3,  betaJ, betaA, p);
    
    for(i=0;i<p->strain_totalJ;i++){
	for (j=0; j<p->strain_totalA; j++){
            SJtemp[i+j*p->strain_totalJ] = SJ[i+j*p->strain_totalJ]+h*(b41*DSJDT[i+j*p->strain_totalJ]+b42*SJk2[i+j*p->strain_totalJ]+b43*SJk3[i+j*p->strain_totalJ]);
	    SAtemp[i+j*p->strain_totalJ] = SA[i+j*p->strain_totalJ]+h*(b41*DSADT[i+j*p->strain_totalJ]+b42*SAk2[i+j*p->strain_totalJ]+b43*SAk3[i+j*p->strain_totalJ]);
            IJtemp[i+j*p->strain_totalJ] = IJ[i+j*p->strain_totalJ]+h*(b41*DIJDT[i+j*p->strain_totalJ]+b42*IJk2[i+j*p->strain_totalJ]+b43*IJk3[i+j*p->strain_totalJ]);
            IAtemp[i+j*p->strain_totalJ] = IA[i+j*p->strain_totalJ]+h*(b41*DIADT[i+j*p->strain_totalJ]+b42*IAk2[i+j*p->strain_totalJ]+b43*IAk3[i+j*p->strain_totalJ]);
	}
    }


    dynamic(SJtemp, SAtemp, IJtemp, IAtemp, SJk4, SAk4, IJk4, IAk4,  betaJ, betaA, p);
    
    for(i=0;i<p->strain_totalJ;i++){
	for (j=0; j<p->strain_totalA; j++){
            SJtemp[i+j*p->strain_totalJ] = SJ[i+j*p->strain_totalJ]+h*(b51*DSJDT[i+j*p->strain_totalJ]+b52*SJk2[i+j*p->strain_totalJ]+b53*SJk3[i+j*p->strain_totalJ]+b54*SJk4[i+j*p->strain_totalJ]);
	    SAtemp[i+j*p->strain_totalJ] = SA[i+j*p->strain_totalJ]+h*(b51*DSADT[i+j*p->strain_totalJ]+b52*SAk2[i+j*p->strain_totalJ]+b53*SAk3[i+j*p->strain_totalJ]+b54*SAk4[i+j*p->strain_totalJ]);
            IJtemp[i+j*p->strain_totalJ] = IJ[i+j*p->strain_totalJ]+h*(b51*DIJDT[i+j*p->strain_totalJ]+b52*IJk2[i+j*p->strain_totalJ]+b53*IJk3[i+j*p->strain_totalJ]+b54*IJk4[i+j*p->strain_totalJ]);
            IAtemp[i+j*p->strain_totalJ] = IA[i+j*p->strain_totalJ]+h*(b51*DIADT[i+j*p->strain_totalJ]+b52*IAk2[i+j*p->strain_totalJ]+b53*IAk3[i+j*p->strain_totalJ]+b54*IAk4[i+j*p->strain_totalJ]);
	}
    }



    dynamic(SJtemp, SAtemp, IJtemp, IAtemp, SJk5, SAk5, IJk5, IAk5,  betaJ, betaA, p); 
    
    for(i=0;i<p->strain_totalJ;i++){
	for (j=0; j<p->strain_totalA; j++){ 
            SJtemp[i+j*p->strain_totalJ] = SJ[i+j*p->strain_totalJ]+h*(b61*DSJDT[i+j*p->strain_totalJ]+b62*SJk2[i+j*p->strain_totalJ]+b63*SJk3[i+j*p->strain_totalJ]+b64*SJk4[i+j*p->strain_totalJ]+b65*SJk5[i+j*p->strain_totalJ]);
	    SAtemp[i+j*p->strain_totalJ] = SA[i+j*p->strain_totalJ]+h*(b61*DSADT[i+j*p->strain_totalJ]+b62*SAk2[i+j*p->strain_totalJ]+b63*SAk3[i+j*p->strain_totalJ]+b64*SAk4[i+j*p->strain_totalJ]+b65*SAk5[i+j*p->strain_totalJ]);
            IJtemp[i+j*p->strain_totalJ] = IJ[i+j*p->strain_totalJ]+h*(b61*DIJDT[i+j*p->strain_totalJ]+b62*IJk2[i+j*p->strain_totalJ]+b63*IJk3[i+j*p->strain_totalJ]+b64*IJk4[i+j*p->strain_totalJ]+b65*IJk5[i+j*p->strain_totalJ]);
            IAtemp[i+j*p->strain_totalJ] = IA[i+j*p->strain_totalJ]+h*(b61*DIADT[i+j*p->strain_totalJ]+b62*IAk2[i+j*p->strain_totalJ]+b63*IAk3[i+j*p->strain_totalJ]+b64*IAk4[i+j*p->strain_totalJ]+b65*IAk5[i+j*p->strain_totalJ]);
	}
    }

    dynamic(SJtemp, SAtemp, IJtemp, IAtemp, SJk6, SAk6, IJk6, IAk6,  betaJ, betaA, p);
    
    for(i=0;i<p->strain_totalJ;i++){
	for (j=0; j<p->strain_totalA; j++){
            SJout[i+j*p->strain_totalJ]= SJ[i+j*p->strain_totalJ]+h*(c1*DSJDT[i+j*p->strain_totalJ]+c3*SJk3[i+j*p->strain_totalJ]+c4*SJk4[i+j*p->strain_totalJ]+c6*SJk6[i+j*p->strain_totalJ]);
            SJerr[i+j*p->strain_totalJ]= h*(dc1*DSJDT[i+j*p->strain_totalJ]+dc3*SJk3[i+j*p->strain_totalJ]+dc4*SJk4[i+j*p->strain_totalJ]+dc5*SJk5[i+j*p->strain_totalJ]+dc6*SJk6[i+j*p->strain_totalJ]);
            SAout[i+j*p->strain_totalJ]= SA[i+j*p->strain_totalJ]+h*(c1*DSADT[i+j*p->strain_totalJ]+c3*SAk3[i+j*p->strain_totalJ]+c4*SAk4[i+j*p->strain_totalJ]+c6*SAk6[i+j*p->strain_totalJ]);
            SAerr[i+j*p->strain_totalJ]= h*(dc1*DSADT[i+j*p->strain_totalJ]+dc3*SAk3[i+j*p->strain_totalJ]+dc4*SAk4[i+j*p->strain_totalJ]+dc5*SAk5[i+j*p->strain_totalJ]+dc6*SAk6[i+j*p->strain_totalJ]);
            IJout[i+j*p->strain_totalJ]= IJ[i+j*p->strain_totalJ]+h*(c1*DIJDT[i+j*p->strain_totalJ]+c3*IJk3[i+j*p->strain_totalJ]+c4*IJk4[i+j*p->strain_totalJ]+c6*IJk6[i+j*p->strain_totalJ]);
            IJerr[i+j*p->strain_totalJ]= h*(dc1*DIJDT[i+j*p->strain_totalJ]+dc3*IJk3[i+j*p->strain_totalJ]+dc4*IJk4[i+j*p->strain_totalJ]+dc5*IJk5[i+j*p->strain_totalJ]+dc6*IJk6[i+j*p->strain_totalJ]);
            IAout[i+j*p->strain_totalJ]= IA[i+j*p->strain_totalJ]+h*(c1*DIADT[i+j*p->strain_totalJ]+c3*IAk3[i+j*p->strain_totalJ]+c4*IAk4[i+j*p->strain_totalJ]+c6*IAk6[i+j*p->strain_totalJ]);
            IAerr[i+j*p->strain_totalJ]= h*(dc1*DIADT[i+j*p->strain_totalJ]+dc3*IAk3[i+j*p->strain_totalJ]+dc4*IAk4[i+j*p->strain_totalJ]+dc5*IAk5[i+j*p->strain_totalJ]+dc6*IAk6[i+j*p->strain_totalJ]);
	}
    }

    /* Free memory */
    free(SJk1);
    free(SJk2);
    free(SJk3);
    free(SJk4);
    free(SJk5);
    free(SJk6);
    free(SJtemp);
    free(IJk1);
    free(IJk2);
    free(IJk3);
    free(IJk4);
    free(IJk5);
    free(IJk6);
    free(IJtemp);
    free(SAk1);
    free(SAk2);
    free(SAk3);
    free(SAk4);
    free(SAk5);
    free(SAk6);
    free(SAtemp);
    free(IAk1);
    free(IAk2);
    free(IAk3);
    free(IAk4);
    free(IAk5);
    free(IAk6);
    free(IAtemp);
}

/**************************************
 * Population and evolutionary dynamics
 **************************************/
void dynamic(double *SJ, double *SA, double *IJ, double *IA,  double *DSJDT, double *DSADT, double *DIJDT, double *DIADT, double *betaJ, double *betaA, struct PARAM *p){
    
    int i, j;
    double N, *sumbetaI, *allIJ;

    sumbetaI = malloc(p->strain_totalJ*sizeof(double));
    allIJ = malloc(p->strain_totalA*sizeof(double));

    /* Population sums */
    N = 0;
    for(i=0;i<p->strain_totalJ;i++){ 
        for (j=0; j<p->strain_totalA; j++){
	    N = N + IJ[i+j*p->strain_totalJ] + SJ[i+j*p->strain_totalJ];
	}
    }

    for(j=0;j<p->strain_totalA;j++){ 
	    allIJ[j]=0;
    }
    for(j=0;j<p->strain_totalA;j++){ 
        for (i=0; i<p->strain_totalJ; i++){
	    allIJ[j]=allIJ[j]+IJ[i+j*p->strain_totalJ];
	}
    }


    for(i=0;i<p->strain_totalJ;i++){ 
	    sumbetaI[i]=0;
    }
    for(i=0;i<p->strain_totalJ;i++){ 
        for (j=0; j<p->strain_totalA; j++){
	    sumbetaI[i]=sumbetaI[i]+allIJ[j]*0.25*(p->beta0min+((p->beta0max-p->beta0min)/(p->p0-p->p1))*(betaA[j]-p->p1))*(betaA[j]-betaJ[i]+4);
	}
    }


    /* ODEs */
    for(i=0;i<p->strain_totalJ;i++){
	for (j=0; j<p->strain_totalA; j++){
	    DSJDT[i+j*p->strain_totalJ] = 0;
            DSADT[i+j*p->strain_totalJ] = 0;
	    DIJDT[i+j*p->strain_totalJ] = 0.25*(p->beta0min+((p->beta0max-p->beta0min)/(p->p0-p->p1))*(betaA[j]-p->p1))*(betaA[j]-betaJ[i]+4)*allIJ[j]*SJ[i] -(p->alpha+p->b)*IJ[i+j*p->strain_totalJ];
	    DIADT[i+j*p->strain_totalJ] = 0;
	}
	DSJDT[i] = (p->a0+((p->a1-p->a0)/(p->h0-p->h1))*(betaJ[i]-p->h1)+p->epsilon*(betaJ[i]-p->h1)*(betaJ[i]-p->h1)-p->q*N)*SJ[i]-p->b*SJ[i]-sumbetaI[i]*SJ[i];
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
