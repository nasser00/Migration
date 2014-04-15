/* Post stack migration (Split Step and PSPI)
*/
/*
  Author Nasser Kazemi
  Copyright (C) 2014 University of Alberta
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 
 Contact information:
 (electronic):
 Nasser Kazemi
 kazemino@ualberta.ca
 
 (paper):
 Nasser Kazemi
 Department of Physics
 CCIS,4-183
 University of Alberta
 Edmonton Alberta CANADA
 T6G 2E1
*/

#include <rsf.h>
#include <fftw3.h>
#include "myfree.h"

#ifndef PI
#define PI (3.141592653589793)
#endif


void  my_forward_fft(sf_complex **m,float **d,int nt,float dt,int nx, int padt);

void my_phase_shift(sf_complex **m,sf_complex czero,int iw,int iz,float *omega,float *kx,int nk,int nx,float v_ave,sf_complex *in2a, sf_complex *in2b,fftwf_plan p2a,fftwf_plan p2b,float dz);

void my_ref_vel (float *vref,float **vmig, int iz,int nx,int nref);

void  my_interpolation (sf_complex **m,sf_complex **tempx,float **vmig,float *vref,int nref,int iz,int iw,int nx);

void  my_v_ave (float v_ave,float **vmig,int iz,int nx);

void  my_split_step_correction (sf_complex **m,sf_complex *c,float **vmig,float v_ave,int iz,float dz,int iw,float dw,int nx);
int main(int argc, char* argv[])
{
  /*define variables*/
  int nx,nt;
  int n1,n2;
  float d1,o1,d2,o2;
  int padt,padx;
  int ntfft,*n,nw,nk;
  float **d,**vel,**vmig,**M,*vref,v_ave;
  float *kx,*omega,dkx,dw;
  sf_complex **m,*in2a,*in2b,**tempx,*c,czero;
  float fmin,fmax,f_low,f_high;
  int if_low,if_high;
  int ix,iw,ik;
  float dt,dx,ox,dz,zmax;
  fftwf_plan p2a,p2b;
  sf_file in,out,velfile;
  int iz,nz,nref,i;
  int option;
  
  /*define sf input output*/
  sf_init (argc,argv);
  in = sf_input("in");
  out = sf_output("out");
  velfile = sf_input("velfile");
  if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
  if (!sf_histfloat(in,"d1",&d1)) sf_error("No d1= in input");
  if (!sf_histfloat(in,"o1",&o1)) o1=0.;
  if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
  if (!sf_histfloat(in,"d2",&d2)) sf_error("No d2= in input");
  if (!sf_histfloat(in,"o2",&o2)) o2=0.;

  if (!sf_histint(velfile,"n1",&nz)) sf_error("No n1= in vel");
  if (!sf_histfloat(velfile,"d1",&dz)) sf_error("No n1= in vel");

  dt = d1;
  dx = d2;
  ox = o2;
  nx = n2;
  nt = n1;
  padt = 2;
  padx = 2;
  ntfft = padt*nt;
  nw=ntfft/2+1;
  nk = padx*nx;
  dw = 2*PI/ntfft/dt;
  dkx = 2*PI/nk/dx;

  sf_putint(out,"n1",nz);
  sf_putfloat(out,"d1",dz);
  sf_putstring(out,"label1","z");
  sf_putstring(out,"unit1","m");
  sf_putstring(out,"title","migrated");

  if (!sf_getfloat("fmax",&fmax)) fmax = 0.5/d1; /* max frequency to process */
  if (fmax > 0.5/d1) fmax = 0.5/d1;
  if (!sf_getfloat("fmin",&fmin)) fmin = 0.1; /* min frequency to process */
  if (!sf_getfloat("Zmax",&zmax)) zmax = (nz-1)*dz; /* max Depth to migrate */
  if (!sf_getint("Number_of_reference_vel",&nref)) nref = 2; /* Number of reference velocities */
  if (!sf_getint("option",&option)) option = 2; /* option=1 for PSPI and option=2 for split step */
  /*define axis variables*/
 
  dkx=(float) 2*PI/nk/dx;
  dw=(float) 2*PI/ntfft/dt;
  /*allocate memory to dynamic arrays*/
  d = sf_floatalloc2(nt,nx);
  vel = sf_floatalloc2(nz,nx);
  vmig = sf_floatalloc2(nz,nx);
  m = sf_complexalloc2(nw,nx);
  kx= sf_floatalloc (nk);
  omega= sf_floatalloc (nw);
  vref=sf_floatalloc (nref);
  in2a = sf_complexalloc(nk);
  in2b = sf_complexalloc(nk);
  n = sf_intalloc(1);
  M= sf_floatalloc2(nz,nx);
  tempx=sf_complexalloc2(nx,nref);
  c = sf_complexalloc(nx);
  /*read input files*/ 
  sf_floatread(d[0],nx*nt,in);
  sf_floatread(vel[0],nx*nz,velfile);
    
  /* apply fft in time direction*/
  my_forward_fft(m,d,nt,dt,nx,padt);

  /* This part is important: we need to define the horizontal wavenumber and frequency axes right.*/

dw = 2*PI/ntfft/dt;

 dkx = 2*PI/nk/dx;

for (iw=0;iw<nw;iw++){

    omega[iw] = dw*iw;

}
for (ik=0;ik<nk;ik++){ 

  if (ik<nk/2) kx[ik] = dkx*ik;

  else         kx[ik] = -(dkx*nk - dkx*ik);

}

  /* Define minimum and maximum frequency index to process*/   
  
  f_low = fmin;   /* min frequency to process */
  f_high = fmax;  /* max frequency to process */
  
  if(f_low>0){
  	if_low = trunc(f_low*dt*ntfft);
  }	
  else{ 
  	if_low = 0;
  }
  if(f_high*dt*ntfft+1<nw){
  	if_high = trunc(f_high*dt*ntfft)+1;
  }
  else{
 	 if_high = nw;
  }
  
  
  __real__ czero = 0;
  __imag__ czero = 0;
  n[0] = nk;
  p2a = fftwf_plan_dft(1, n, (fftwf_complex*)in2a, (fftwf_complex*)in2a, FFTW_FORWARD, FFTW_ESTIMATE);
  p2b = fftwf_plan_dft(1, n, (fftwf_complex*)in2b, (fftwf_complex*)in2b, FFTW_BACKWARD, FFTW_ESTIMATE);

fftwf_execute(p2a); /* FFT x to k */
fftwf_execute(p2b); /* FFT x to k */

  /* Define initial migrated model as zeros*/
 
  for (iz=0; iz<nz; iz++) {	
    for (ix=0; ix<nx; ix++) M[ix][iz] = 0.0;
  }
  
  for (iz=0; iz<nz;iz++){
    for (ix=0;ix<nx;ix++){
      vmig[ix][iz]=vel[ix][iz]/2;
    }
  }
  

/* option==1 is PSPI and option==2 is Split Step */

if (option==1){
   for (iw=if_low;iw<if_high;iw++){
     for (iz=0; iz<nz;iz++){	
	        my_ref_vel (vref,vmig,iz,nx,nref);
		 for (i=0;i<nref;i++){
			v_ave=vref[i];
			 my_phase_shift(m,czero,iw,iz,omega,kx,nk,nx,v_ave,in2a,in2b,p2a,p2b,dz);
			 for (ix=0;ix<nx;ix++) tempx[i][ix]= in2b[ix];
		}
         	my_interpolation (m,tempx,vmig,vref,nref,iz,iw,nx);

   			/* Apply imaging condition*/	
		 for (ix=0;ix<nx;ix++) M[ix][iz]=M[ix][iz]+2*crealf(m[ix][iw]);
     }
    fprintf(stderr,"\r progress = %6.2f%%",(float) 100*(iw-if_low)/(if_high-if_low));
   }
  fprintf(stderr,"\r                            \n");
  }

if (option==2){
  for (iw=if_low;iw<if_high;iw++){
     for (iz=0; iz<nz;iz++){
		v_ave=vmig[0][iz];	
	        my_v_ave (v_ave,vmig,iz,nx);
	        my_phase_shift(m,czero,iw,iz,omega,kx,nk,nx,v_ave,in2a,in2b,p2a,p2b,dz);
			for (ix=0;ix<nx;ix++) {
					c[ix]= in2b[ix];
					}

         	my_split_step_correction (m,c,vmig,v_ave,iz,dz,iw,dw,nx);

   			/* Apply imaging condition*/	
		for (ix=0;ix<nx;ix++) M[ix][iz]=M[ix][iz]+2*crealf(m[ix][iw]);
      }
    fprintf(stderr,"\r progress = %6.2f%%",(float) 100*(iw-if_low)/(if_high-if_low));
   }
  fprintf(stderr,"\r                            \n");
 }


  sf_floatwrite(M[0],nz*nx,out);

  
  fftwf_destroy_plan(p2a);
  fftwf_free(in2a); 
  fftwf_destroy_plan(p2b);
  fftwf_free(in2b);	
  
  exit (0);
}



void my_forward_fft(sf_complex **m,float **d,int nt,float dt,int nx, int padt)
{
  int nw;
  sf_complex *out1a;
  float *in1a;
  int ntfft,ix,it,iw;
  fftwf_plan p1a;
  
  ntfft = nt*padt;
  nw = ntfft/2 + 1;
  out1a = sf_complexalloc(nw);
  in1a = sf_floatalloc(ntfft);
  p1a = fftwf_plan_dft_r2c_1d(ntfft, in1a, (fftwf_complex*)out1a, FFTW_ESTIMATE);

  for (ix=0;ix<nx;ix++){
  	for(it=0;it<nt;it++) {
		in1a[it] = d[ix][it];
		}
    	for(it=nt;it<ntfft;it++){
		in1a[it] = 0;
		}
    fftwf_execute(p1a); 
    for(iw=0;iw<nw;iw++) {
	m[ix][iw] = out1a[iw]; 
	}
  }
  fftwf_destroy_plan(p1a);
  fftwf_free(in1a); fftwf_free(out1a);
  return;

}
void my_phase_shift(sf_complex **m,sf_complex czero,int iw,int iz,float *omega,float *kx,int nk,int nx,float v_ave,sf_complex *in2a, sf_complex *in2b,fftwf_plan p2a,fftwf_plan p2b,float dz)
{
  int ix;
  float shift;
  sf_complex L;

  for(ix=0;ix<nx;ix++) in2a[ix] = m[ix][iw];
  for(ix=nx;ix<nk;ix++) in2a[ix] = czero;
  fftwf_execute(p2a); /* FFT x to k */
  for (ix=0;ix<nk;ix++){
    shift= ((omega[iw]*omega[iw])/(v_ave*v_ave)) - (kx[ix]*kx[ix]);
    if (shift>=0){
      __real__ L= cos(sqrtf(shift)*dz);
      __imag__ L= sin(sqrtf(shift)*dz);
    } 
    else{
      L=czero;
    }
    in2b[ix]=in2a[ix]*L/((float) nk);
  }
  fftwf_execute(p2b);  

  return;
}

void my_ref_vel (float *vref,float **vmig, int iz,int nx,int nref)
{
float min_vz,max_vz,b,dv;
int i,ix;
min_vz=vmig[0][iz];
 max_vz=vmig[0][iz];
for (ix=0;ix<nx;ix++){
b=vmig[ix][iz];
if (b <= min_vz )min_vz= b;
if (b >= max_vz) max_vz= b;
}
dv=(max_vz-min_vz)/ (float) (nref-1);
for (i=0;i<nref;i++) vref[i]=min_vz+ (float) i*dv;
return;
}

void  my_interpolation (sf_complex **m,sf_complex **tempx,float **vmig,float *vref,int nref,int iz,int iw,int nx)
{
int i1,i2,ix;
float dv,a1,a2;
i1=0;
i2=1;
dv=(vref[nref-1]-vref[0])/(float) (nref-1);
      /* apply interpolation to the downward continuated data*/
for (ix=0;ix<nx;ix++){
 	if (dv>= (vmig[ix][iz]/(float) 1000) ){ 
	i1=trunc((vmig[ix][iz]-vref[0])/dv);
		if (i1<0)i1=0;
		if (i1>nref-2)i1=nref-2;
	i2=i1+1;
	a1= (vref[i2]-vmig[ix][iz])/(vref[i2]-vref[i1]);
	a2=(vmig[ix][iz]-vref[i1])/(vref[i2]-vref[i1]);
	__real__ m[ix][iw]= a1* crealf(tempx[i1][ix])+a2*crealf(tempx[i2][ix]);
	__imag__ m[ix][iw]= a1*cimagf(tempx[i1][ix])+a2*cimagf(tempx[i2][ix]);
	}
	if (dv< (vmig[ix][iz]/(float) 1000) ) {
	__real__ m[ix][iw]= (crealf(tempx[0][ix])+crealf(tempx[nref-1][ix]))/2;
	__imag__ m[ix][iw]= (cimagf(tempx[0][ix])+cimagf(tempx[nref-1][ix]))/2;
	}
}
return;
}

void  my_v_ave (float v_ave,float **vmig,int iz,int nx)
{
int ix;
	for (ix=1;ix<nx;ix++) {
		v_ave=v_ave+vmig[ix][iz];
		}
v_ave=v_ave/ ((float) nx);
return;
}


void  my_split_step_correction (sf_complex **m,sf_complex *c,float **vmig,float v_ave,int iz,float dz,int iw,float dw,int nx)
{
int ix;
float dp;
sf_complex L;

      /* apply split step correction to the downward continuated data*/
for (ix=0;ix<nx;ix++){
 	dp=(v_ave-vmig[ix][iz])/(vmig[ix][iz]*v_ave);/* dp=1/v-1/v_ave;*/
        __real__ L= cos(iw*dw*dp*dz);
        __imag__ L= sin(iw*dw*dp*dz);
	 m[ix][iw]=c[ix]*L;
	}
return;
}
