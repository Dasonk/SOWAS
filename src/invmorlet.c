#include <R.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#ifndef Macintosh
#include <sys/file.h>
#include <sys/types.h>
#endif

#include <time.h>

#ifndef Macintosh
#include <sys/time.h>
#endif

#define NR_END 1
#define FREE_ARG char*
#define float double

double *dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) Rprintf("allocation failure in dvector()");
	return v-nl+NR_END;
}

double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
	if (!m) Rprintf("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) Rprintf("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

int **imatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	int **m;

	/* allocate pointers to rows */
	m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
	if (!m) Rprintf("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;


	/* allocate rows and set pointers to them */
	m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
	if (!m[nrl]) Rprintf("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}


void free_dvector(double *v, long nl, long nh)
/* free a double vector allocated with dvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)
/* free an int matrix allocated by imatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}


#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

void four1(double data[], unsigned long nn, int isign)
{
	unsigned long n,mmax,m,j,istep,i;
	double wtemp,wr,wpr,wpi,wi,theta;
	double tempr,tempi;

	n=nn << 1;
	j=1;
	for (i=1;i<n;i+=2) {
		if (j > i) {
			SWAP(data[j],data[i]);
			SWAP(data[j+1],data[i+1]);
		}
		m=n >> 1;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	mmax=2;
	while (n > mmax) {
		istep=mmax << 1;
		theta=isign*(6.28318530717959/mmax);
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1;m<mmax;m+=2) {
			for (i=m;i<=n;i+=istep) {
				j=i+mmax;
				tempr=wr*data[j]-wi*data[j+1];
				tempi=wr*data[j+1]+wi*data[j];
				data[j]=data[i]-tempr;
				data[j+1]=data[i+1]-tempi;
				data[i] += tempr;
				data[i+1] += tempi;
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
		}
		mmax=istep;
	}
}
#undef SWAP


//Morlet Wavelet in Fourier Domain
double morletFD(double x, int c, double w0)
{

  double r;

  if (x>0)
    r=pow(PI,-0.25)*exp(-(x-w0)*(x-w0)/2);
  else
    r=0;
  
  return r;

}

//Scaling of Wavelet
double scaling(double s, int scal, int N, double dt)
{

  double factor;

  // scal = 0 Scaling preserving peak hight
  // scal = 1 Scaling preserving flat white noise

  factor = 1./(1.*N) * sqrt(2.*PI*pow(s,(double)scal)/dt);

  return factor;

}


// Continuous Wavelet Transform
void cwt(double *xre, int N2, double dt, double w0, double *s, int S, double **cwc, int scal)
{

  int i=0,j=0,k=0;
  double *w,*x,wk=0,fs;
  
  int  N=pow(2,1+(int)(log((double)N2)/log(2)-0.00000000001));

  w=dvector(1,2*N);
  x=dvector(1,2*N);

  //Complex signal with IM=0
  for (i=1;i<=N2;i++){
    x[2*i-1]=xre[i];
    x[2*i]=0;
  }

  //Zeropadding
  if (N>N2){
    for (i=N2+1;i<=N;i++){
      x[2*i-1]=0;
      x[2*i]=0;
    }
  }

  // Fouriertrafo
  four1(x,N,-1);

  for (i=1;i<=S;i++){ //Scale s[i]

    fs= w0/(2*PI);
   
    //Sampling Wavelet and multiplying with FT of time series
    for (k=1;k<=N;k++){

      //respecting order of frequencies in FT-routine
      if (k<=N/2)
	wk=2.*PI*(k-1)/(1.*N*dt);
      else
	wk=-2.*PI*(double)(N+1-k)/(1.*N*dt);

      //real part 
      w[2*k-1] = scaling(fs*s[i],scal,N,dt) * x[2*k-1] * morletFD(wk*fs*s[i],0,w0);
      //imaginary part
      w[2*k]   = scaling(fs*s[i],scal,N,dt) * x[2*k] * morletFD(wk*fs*s[i],1,w0);
      
    }
    
    // Inverse Trafo (not really inverse, but because of Normalization)
    four1(w,N,1);
    
    for (k=1;k<=N2;k++){
      cwc[2*k-1][i]=w[2*k-1];//real part
      cwc[2*k][i]=w[2*k];    //imaginary part
    }
      
  }

  free_dvector(w,1,2*N);
  free_dvector(x,1,2*N);

}

double scalinv(double s,int scal){

  return 1./sqrt(pow(s,(double)scal));

}


// Inverse Wavelet Transformation
void invcwt(double **cwc, int N, double dt, double w0, double *s, int S, double *cts, int scal)
{

  int i=0,j=0,k=0;
  double *rewt,*imwt,*w,*x,*sv,**wtr,**wti,wk=0,fs;

  wtr=dmatrix(1,2*N,1,1);
  wti=dmatrix(1,2*N,1,1);

  sv=dvector(1,1);
  rewt=dvector(1,N);
  imwt=dvector(1,N);

  for (k=1;k<=2*N;k++)
    cts[k]=0;

  for (i=1;i<=S;i++){

    sv[1]=s[i];

    for (k=1;k<=N;k++){
      rewt[k]=cwc[2*k-1][i]; //real part
      imwt[k]=cwc[2*k][i];   //imaginary part

    }
    
    cwt(rewt,N,dt,w0,sv,1,wtr,0);
    cwt(imwt,N,dt,w0,sv,1,wti,0);

    for (k=1;k<=N;k++){

      cts[2*k-1] += scalinv(sv[1],scal)*(wtr[2*k-1][1]-wti[2*k][1]);//real part
      cts[2*k]   += scalinv(sv[1],scal)*(wti[2*k][1]+wti[2*k-1][1]);//imaginary part

    }

  }

  free_dmatrix(wtr,1,2*N,1,1);
  free_dmatrix(wti,1,2*N,1,1);
  free_dvector(sv,1,1);
  free_dvector(rewt,1,N);
  free_dvector(imwt,1,N);

}


void invmorlet(double *wtReInput, double *wtImInput, int *pn, 
double *pdt, double *ps0, int *pnoctave, int *pnvoice, 
	       double *pw0, double *invReOut, double *invImOut){

  int i,j,m,n=*pn,noctave=*pnoctave,nvoice=*pnvoice;
  double **wt,*ts,dt=*pdt,s0=*ps0,*s,w0=*pw0;

  m=noctave*nvoice+1;

  wt=dmatrix(1,2*n,1,m);
  ts=dvector(1,2*n);
  s=dvector(1,m);

  for (i=1;i<=2*n;i++)
    ts[i]=0;

  //Get Input
  for (i=1;i<=n;i++){
    for (j=1;j<=m;j++){
      wt[2*i-1][j]=wtReInput[(i-1)*m+j-1];
      wt[2*i][j]=wtImInput[(i-1)*m+j-1];
    }
  }

  for (i=1;i<=m;i++)
    s[i]=s0*exp(((i-1.)/(1.*nvoice))*log(2.));

  invcwt(wt,n,dt,w0,s,m,ts,1);

  //Output
  for (i=1;i<=n;i++){
    invReOut[i-1]=ts[2*i-1];
    invImOut[i-1]=ts[2*i];
  }
  
  free_dmatrix(wt,1,2*n,1,m);
  free_dvector(ts,1,2*n);
  free_dvector(s,1,m);

}
