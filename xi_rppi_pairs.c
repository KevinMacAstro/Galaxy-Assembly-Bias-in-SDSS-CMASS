#include <malloc.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PI 3.1415926535897932384626433832795029

void FO_error(char *filename)
{
  printf("Cannot open file %s\n",filename);
  exit(0);
}

double *vector(int nl,int nh)
{
        double *v;

        v=(double *)malloc((unsigned) (nh-nl+1)*sizeof(double));
        if (!v) {printf("allocation failure in vector()\n"); exit(1);}
        return v-nl;
}

void free_vector(double *v,int nl,int nh)
{
        free((char*) (v+nl));
}

int *ivector(int nl,int nh)
{
        int *v;

        v=(int *)malloc((unsigned) (nh-nl+1)*sizeof(int));
        if (!v) {printf("allocation failure in ivector()\n"); exit(1);}
        return v-nl;
}

void free_ivector(int *v,int nl,int nh)
{
        free((char*) (v+nl));
}

double **matrix(nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch;
{
        int i;
        double **m;

        m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
        if (!m) {printf("allocation failure 1 in matrix()\n"); exit(1);}
        m -= nrl;

        for(i=nrl;i<=nrh;i++) {
                m[i]=(double *) malloc((unsigned) (nch-ncl+1)*sizeof(double));
                if (!m[i]) {printf("allocation failure 2 in matrix()\n"); exit(1);}
                m[i] -= ncl;
        }
        return m;
}

void free_matrix(m,nrl,nrh,ncl,nch)
double **m;
int nrl,nrh,ncl,nch;
{
        int i;

        for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
        free((char*) (m+nrl));
}

void helpmessage(char *com)
{
  printf("\nGiven the simulation output, calculate the correlation function.\n");
  printf("Format of the input file: No. #Particle x y z?\n\n"
         "Usage: %s inputfile boxsize N rmin rmax [#bin]\n"
         "  boxsize - boxsize of the simulation (h^-1Mpc)\n"
         "  N       - the first N halos in the input file (0 means all halos)\n"
         "  rmin    - outer radius of the first bin\n"
         "  rmax    - outer radius of the last bin\n"
         "  #bin    - number of bins (logarithmically uniform spaced) [defalt=10]\n\n",com);
  printf("Output Format: r_outer r_ave xi xi_err Npair Npspair\n"
         "  r_outer, r_ave  - outer and volumn-weighted radius of each bin\n"
         "  xi, xi_err      - correlation function with errorbars\n"
         "  Npair, Npspair  - numbers of pairs and poisson pairs in each bin\n\n"); 
  exit(0);
}

void PairCount(double *rp_cen, double *pi_cen, double **xi,double **Npair, double *wp,double *wp_error,
               double rp_min,  double rp_max,  int Nbin_rp,
               double   pi_min,  double   pi_max,  int Nbin_pi,
               double *x1, double *y1, double *z1, int N1, 
               double *x2, double *y2, double *z2, int N2, 
               int auto_opt,
               double Lx,  double Ly,  double Lz) 
{
 int    i,ic,j,jc,j0;
 double drp,mu,dpi,rp,s_pi,d1x,d2x,d1y,d2y,d1z,d2z,l,s,sx,sy,sz,lx,ly,lz;
 double rp_lo,rp_hi,navg1,navg2,*wpNpair;
 double S,V,factor,Nranpair;

 if(auto_opt!=0 && N1!=N2) { 
   printf("PairCount> Error, N1!=N2 for option auto\n");
   exit(0);
 }


 wpNpair      = vector(1,Nbin_rp);


 // rp_min=pow(10.0,lgrp_min);
 //rp_max=pow(10.0,lgrp_max);
 drp=log10(pow(10,rp_max)/pow(10,rp_min))/Nbin_rp;/*Bin width perpendicular to los*/
 dpi  =(  pi_max-  pi_min)/Nbin_pi;/*Bin width line of sight*/
 for(i=1;i<=Nbin_rp;i++) rp_cen[i]=pow(10,(rp_min+(i-0.5)*drp));
 for(i=1;i<=Nbin_pi;i++) pi_cen[i]=pi_min+(i-0.5)*dpi;
 for(i=1;i<=Nbin_rp;i++){
    wp_error[i]=0.0;
    for(j=1;j<=Nbin_pi;j++){ 
       xi[i][j]=0.0;
       Npair[i][j]=0.0;
    }
 }
 for(i=1;i<=N1;i++) {
    if(auto_opt==0) j0=1;
    else j0=i+1;          /** for auto correlation **/

    for(j=j0;j<=N2;j++) {
      // dz=z1[i]-z2[j];
       /** apply periodic boundary condition for all directions **/
       //if(dx> 0.5*Lx) dx=dx-Lx;
       //if(dx<-0.5*Lx) dx=dx+Lx;
       //if(dy> 0.5*Ly) dy=dy-Ly;
       //if(dy<-0.5*Ly) dy=dy+Ly;
       //if(dz> 0.5*Lz) dz=dz-Lz;
       //if(dz<-0.5*Lz) dz=dz+Lz;
       
       //s_pi=fabs(dz);	 
       d1x=z1[i]*sin(y1[i])*cos(x1[i]);
       d1y=z1[i]*sin(y1[i])*sin(x1[i]);
       d1z=z1[i]*cos(y1[i]);
       d2x=z2[j]*sin(y2[j])*cos(x2[j]);
       d2y=z2[j]*sin(y2[j])*sin(x2[j]);
       d2z=z2[j]*cos(y2[j]);
       sx=d2x-d1x;
       sy=d2y-d1y;
       sz=d2z-d1z;
       lx=d1x+sx/2.;
       ly=d1y+sy/2.;
       lz=d1z+sz/2.;
       s=sqrt(sx*sx+sy*sy+sz*sz);
       l=sqrt(lx*lx+ly*ly+lz*lz);
       mu=fabs((lx*sx+ly*sy+lz*sz)/(s*l));    
       s_pi=s*mu;
       rp=s*sqrt(1-mu*mu);


//rp=2*l*tan(dphi/2);



	//if(z1[i]<z2[j]) l=z1[i]+s_pi/2.;
      //else l=z2[j]+s_pi/2.;
       //rp=2*l*tan(dphi/2.);

       if(log10(rp)>=rp_min && log10(rp)<=rp_max && s_pi>=pi_min && s_pi<=pi_max) {
         ic=floor(log10(rp/pow(10,rp_min))/drp  +1.0);
         jc=floor((s_pi  -pi_min)/dpi  +1.0);
         if (ic<1) ic=1;
         else if(ic>Nbin_rp) ic=Nbin_rp;
         if (jc<1) jc=1;
         else if(jc>Nbin_pi) jc=Nbin_pi;
         //if(1<=ic && ic<=Nbin_rp && 1<=jc && jc<=Nbin_pi)
         Npair[ic][jc]=Npair[ic][jc]+1.0;
       } 
    }
 }

 /** compute the random pairs, obtain the correlation function **/
// V=Lx*Ly*Lz;
// navg1=N1/V;
// navg2=N2/V;
// if(auto_opt) factor=0.5;
// else         factor=1.0;
// for(i=1;i<=Nbin_rp;i++) {
// 	rp_lo=rp_min+(i-1.0)*drp;
// 	rp_hi=rp_lo+drp;
 	//rp_lo=pow(10.0,rp_lo);
 	//rp_hi=pow(10.0,rp_hi);
// 	S=M_PI*(pow(pow(10,rp_hi),2)-pow(pow(10,rp_lo),2));/*Area of a washer*/
// 	wp[i]=0.0;
//	wpNpair[i]=0.0;
// 	for(j=1;j<=Nbin_pi;j++) {
/// 		Nranpair = factor*N1*navg2*S*dpi*2.0; 
//        	xi[i][j] = Npair[i][j]/Nranpair - 1.0; 
//       		wp[i]    = wp[i] + xi[i][j] ;
//		wpNpair[i]=wpNpair[i]+Npair[i][j];
//	}
//	wp_error[i]=sqrt(wpNpair[i])/Nranpair;
//	 wp[i]=wp[i]*2.0*dpi;
// }


}

int main(int argc, char *argv[])
{
 int    i,j,auto_opt;
 double rp_min,rp_max,pi_min,pi_max;
 double Lx,Ly,Lz;

 double *wp_error,*rp_cen,*pi_cen,**Npair,**xi,*wp;
 int    Nbin_rp,Nbin_pi;

 double xtmp,ytmp,ztmp;
 double *x1,*y1,*z1,*x2,*y2,*z2;
 int    N1,N2;

 double Npair_tot;

 char   file1[128],file2[128],filename[128];
 FILE   *fp1,*fp2,*fp;

 int    Nx,Ny;
 float  ftmp;

 if(argc!=13) { 
   printf("Count pairs ...\n"
          "Usage: %s catalog1 catalog2 ln(rp_min) ln(rp_max) pi_min pi_max "
          "Nbin_rp Nbin_pi Lx Ly Lz Fileout4xi2d.dat\n",argv[0]);
   exit(0);
 }

 strcpy(file1,argv[1]);
 strcpy(file2,argv[2]); 
 rp_min=atof(argv[3]);
 rp_max=atof(argv[4]);
 pi_min  =atof(argv[5]);
 pi_max  =atof(argv[6]);
 Nbin_rp=atoi(argv[7]);
 Nbin_pi=atoi(argv[8]);
 Lx     =atof(argv[9]);
 Ly     =atof(argv[10]);
 Lz     =atof(argv[11]);

 if(rp_min>=rp_max) {
    printf("Well, rp_min < rp_max, please.\n");
    exit(0);
 }
 if(pi_min>=pi_max) {
    printf("Well, pi_min < pi_max, please.\n");
    exit(0);
 }

 if(Nbin_rp<=1 || Nbin_pi<=1) {
    printf("Well, Nbin_rp and Nbin_pi > 1, please.\n");
    exit(0);
 }

/*** read files and define arrays **/
 if(strcmp(file1,file2)==0) { 
    auto_opt=1;
    if(!(fp1=fopen(file1,"rt"))) {
       printf("cannot find %s\n",file1);
       exit(0);
    }
    fscanf(fp1,"%lf %lf %lf",&xtmp,&ytmp,&ztmp);
    N1=0;
    while(!feof(fp1)) {
          N1=N1+1;
          fscanf(fp1,"%lf %lf %lf",&xtmp,&ytmp,&ztmp);
    }
    rewind(fp1);
    N2=N1;
 }
 else {
    auto_opt=0;
    if(!(fp1=fopen(file1,"rt"))) {
         printf("cannot find %s\n",file1);
         exit(0);
    }
    if(!(fp2=fopen(file2,"rt"))) {
         printf("cannot find %s\n",file2);
         exit(0);
    }

    fscanf(fp1,"%lf %lf %lf",&xtmp,&ytmp,&ztmp);
    N1=0;
    while(!feof(fp1)) { 
            N1=N1+1;
            fscanf(fp1,"%lf %lf %lf",&xtmp,&ytmp,&ztmp);
    }
    rewind(fp1);

    fscanf(fp2,"%lf %lf %lf",&xtmp,&ytmp,&ztmp);
    N2=0;
    while(!feof(fp2)) {
            N2=N2+1;
            fscanf(fp2,"%lf %lf %lf",&xtmp,&ytmp,&ztmp);
    }
    rewind(fp2);
 }

 rp_cen= vector(1,Nbin_rp);
 pi_cen  = vector(1,Nbin_pi);
 wp      = vector(1,Nbin_rp);
 wp_error= vector(1,Nbin_rp);
 xi      = matrix(1,Nbin_rp,1,Nbin_pi);
 x1      = vector(1,N1);
 y1      = vector(1,N1);
 z1      = vector(1,N1);
 x2      = vector(1,N2);
 y2      = vector(1,N2);
 z2      = vector(1,N2);
 Npair = matrix(1,Nbin_rp,1,Nbin_pi);


 if(auto_opt==1) { 
    for(i=1;i<=N1;i++) { 
        fscanf(fp1,"%lf %lf %lf",&(x1[i]),&(y1[i]),&(z1[i]));
        x2[i]=x1[i];
        y2[i]=y1[i];
        z2[i]=z1[i];
    }
    fclose(fp1);
 }
 else {
    for(i=1;i<=N1;i++)
         fscanf(fp1,"%lf %lf %lf",&(x1[i]),&(y1[i]),&(z1[i]));
    for(i=1;i<=N2;i++)
         fscanf(fp2,"%lf %lf %lf",&(x2[i]),&(y2[i]),&(z2[i]));
    fclose(fp1);
    fclose(fp2);
 }
/*
 if(auto_opt==1) Npair_tot=0.5*N1*(N1-1.0);
 else            Npair_tot=(N1+0.0)*(N2+0.0);
 printf("N1=%d N2=%d Npair_tot=%e\n",N1,N2,Npair_tot);
 exit(0);
*/

 printf("# %s %s ln[r_p]=(%6.2f,%6.2f) r_pi=(%6.2f,%6.2f) %d r_p bins %d r_pi bins L_x=%5.2fh^-1Mpc L_y=%5.2fh^-1Mpc L_z=%5.2fh^-1Mpc %d Halos\n",file1,file2,rp_min,rp_max,pi_min,pi_max,Nbin_rp,Nbin_pi,Lx,Ly,Lz,N2);

 PairCount(rp_cen, pi_cen, xi, Npair, wp, wp_error,
           rp_min, rp_max, Nbin_rp,
             pi_min,   pi_max, Nbin_pi,
           x1, y1, z1, N1,
           x2, y2, z2, N2, auto_opt,
           Lx, Ly, Lz);

 if(auto_opt==1) Npair_tot=0.5*N1*(N1-1.0);
 else            Npair_tot=(N1+0.0)*(N2+0.0);

 /** output the pair count: lgr Npair_normalized Npair **/ 
 /*
 printf("%d ",Nbin_rp);
 for(i=1;i<=Nbin_rp;i++) printf("%f ",lgrp_cen[i]);
 printf("\n");

 printf("%d ",Nbin_pi);
 for(i=1;i<=Nbin_pi;i++) printf("%f ",pi_cen[i]);
 printf("\n");
 for(i=1;i<=Nbin_rp;i++) 
    for(j=1;j<=Nbin_pi;j++) 
       printf("%12.6e\n",xi[i][j]);
 */




 /** xi(rp,pi) **/
 //Nx=Nbin_rp;
 //Ny=Nbin_pi;
 //strcpy(filename,"xi_rp_pi.bin");
 //if(!(fp=fopen(filename,"wb"))) FO_error(filename);    /* write BINARY file*/
 //fwrite(&Nx,4,1,fp);                                   /* SM header */
 //fwrite(&Ny,4,1,fp);
 //for(j=1;j<=Ny;j++)
 //  for(i=1;i<=Nx;i++) {
 //     ftmp=xi[i][j];
 //     fwrite(&ftmp,4,1,fp);
 //  }
 //fclose(fp);

 char fileout[128];
 strcpy(fileout,argv[12]);
 if(!(fp=fopen(fileout,"wt"))) FO_error(fileout);

 for (j=1;j<=Nbin_pi;j++)
   for(i=1;i<=Nbin_rp;i++)
     fprintf(fp,"%f %f %f\n",rp_cen[i],pi_cen[j],Npair[i][j]);
   
 fclose(fp);






 


 free_vector(rp_cen, 1,Nbin_rp);
 free_vector(pi_cen,   1,Nbin_pi);
 free_vector(wp,       1,Nbin_rp);
 free_vector(wp_error, 1, Nbin_rp);
 free_matrix(xi,     1,Nbin_rp,1,Nbin_pi);
 free_matrix(Npair,     1,Nbin_rp,1,Nbin_pi);
 free_vector(x1,     1,N1);
 free_vector(y1,     1,N1);
 free_vector(z1,     1,N1);
 free_vector(x2,     1,N2);
 free_vector(y2,     1,N2);
 free_vector(z2,     1,N2);
} 
