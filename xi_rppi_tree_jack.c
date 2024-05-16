#include <malloc.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NR_END 1
#define FREE_ARG char*
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
int **imatrix(nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch;
{
        int i;
        int **m;

        m=(int **) malloc((unsigned) (nrh-nrl+1)*sizeof(int*));
        if (!m) {printf("allocation failure 1 in matrix()\n"); exit(1);}
        m -= nrl;

        for(i=nrl;i<=nrh;i++) {
                m[i]=(int *) malloc((unsigned) (nch-ncl+1)*sizeof(int));
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

void free_imatrix(m,nrl,nrh,ncl,nch)
int **m;
int nrl,nrh,ncl,nch;
{
        int i;

        for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
        free((char*) (m+nrl));
}

double ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a double 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
        long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
        double ***t;

        /* allocate pointers to pointers to rows */
        t=(double ***) malloc((size_t)((nrow+NR_END)*sizeof(double**)));
        if (!t) {printf("allocation failure 1 in f3tensor()\n"); exit(1);}
        t += NR_END;
        t -= nrl;

        /* allocate pointers to rows and set pointers to them */
        t[nrl]=(double **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double*)));
        if (!t[nrl]) {printf("allocation failure 2 in f3tensor()\n");exit(1);}
        t[nrl] += NR_END;
        t[nrl] -= ncl;

        /* allocate rows and set pointers to them */
        t[nrl][ncl]=(double *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(double)));
        if (!t[nrl][ncl]) {printf("allocation failure 3 in f3tensor()\n");exit(1);}
        t[nrl][ncl] += NR_END;
        t[nrl][ncl] -= ndl;

        for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
        for(i=nrl+1;i<=nrh;i++) {
                t[i]=t[i-1]+ncol;
                t[i][ncl]=t[i-1][ncl]+ncol*ndep;
                for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
        }

        /* return pointer to array of pointers to rows */
        return t;
}

void free_f3tensor(double ***t, long nrl, long nrh, long ncl, long nch,
        long ndl, long ndh)
/* free a float f3tensor allocated by f3tensor() */
{
        free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
        free((FREE_ARG) (t[nrl]+ncl-NR_END));
        free((FREE_ARG) (t+nrl-NR_END));
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

void PairCount(double *rp_cen, double *pi_cen,
	       double ***Npair,
               double rp_min,  double rp_max,  int Nbin_rp,
               double pi_min,  double pi_max,  int Nbin_pi,
               double *x1, double *y1, double *z1,int *jack1, int N1, 
               double *x2, double *y2, double *z2,int *jack2, int N2, 
               int auto_opt,
               double Lx,  double Ly,  double Lz, int jack_samp, int **cellNEI, int **galIND_1,int **galIND_2, int N3, int N4, int N5) 
{
 int    i,ic,j,jc,j0,ii,k,kc;
 double drp,mu,dpi,rp,s_pi,d1x,d1y,d1z,d2x,d2y,d2z,s,l,lx,ly,lz,sx,sy,sz;
 double rp_lo,rp_hi,navg1,navg2;
 double S,V,factor;

 
 if(auto_opt!=0 && N1!=N2) { 
   printf("PairCount> Error, N1!=N2 for option auto\n");
   exit(0);
 }

 drp=log10(pow(10,rp_max)/pow(10,rp_min))/Nbin_rp;/*Bin width perpendicular to los*/
 dpi  =(  pi_max-  pi_min)/Nbin_pi;/*Bin width line of sight*/
 for(i=1;i<=Nbin_rp;i++) rp_cen[i]=pow(10,(rp_min+(i-0.5)*drp));
 for(i=1;i<=Nbin_pi;i++) pi_cen[i]=pi_min+(i-0.5)*dpi;
 for(i=1;i<=Nbin_rp;i++){ 
    for(j=1;j<=Nbin_pi;j++){
	for(k=1;k<=jack_samp;k++){
       		Npair[i][j][k]=0.0;
	}
    }
 }

 for(ii=1;ii<=N3;ii++) {
	if(galIND_1[ii][1]==0) continue;
	for(k=1;k<=27;k++){
		if(cellNEI[ii][k]==0) break;
		if(galIND_2[cellNEI[ii][k]][1]==0) continue;
		if(auto_opt==1 && cellNEI[ii][k]<ii) continue;
		if(auto_opt==1 && cellNEI[ii][k]==ii){
			if(galIND_2[ii][1]==galIND_2[ii][2]) continue;
			for(i=galIND_1[ii][1];i<=galIND_1[ii][2];i++){
				for(j=(i+1);j<=galIND_2[cellNEI[ii][k]][2];j++){
       			
			d1x=x1[i];
       			d1y=y1[i];
       			d1z=z1[i];
       			d2x=x2[j];
       			d2y=y2[j];
       			d2z=z2[j];
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


       			if(log10(rp)>=rp_min && log10(rp)<=rp_max && s_pi>=pi_min && s_pi<=pi_max) {
         		ic=floor(log10(rp/pow(10,rp_min))/drp  +1.0);
         		jc=floor((s_pi  -pi_min)/dpi  +1.0);
         		if (ic<1) ic=1;
         		else if(ic>Nbin_rp) ic=Nbin_rp;
         		if (jc<1) jc=1;
         		else if(jc>Nbin_pi) jc=Nbin_pi;
                        for(kc=1;kc<=jack_samp;kc++) {
                        	if(jack1[i] != kc && jack2[j] != kc) Npair[ic][jc][kc]=Npair[ic][jc][kc]+1;
                                } 
			}
				}
			}
		}
       		else {
			for(i=galIND_1[ii][1];i<=galIND_1[ii][2];i++){
                        	for(j=galIND_2[cellNEI[ii][k]][1];j<=galIND_2[cellNEI[ii][k]][2];j++){	
			d1x=x1[i];
                        d1y=y1[i];
                        d1z=z1[i];
                        d2x=x2[j];
                        d2y=y2[j];
                        d2z=z2[j];
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

                        
                        if(log10(rp)>=rp_min && log10(rp)<=rp_max && s_pi>=pi_min && s_pi<=pi_max) {
                        ic=floor(log10(rp/pow(10,rp_min))/drp  +1.0);
                        jc=floor((s_pi  -pi_min)/dpi  +1.0);
                        if (ic<1) ic=1;
                        else if(ic>Nbin_rp) ic=Nbin_rp;
                        if (jc<1) jc=1;
                        else if(jc>Nbin_pi) jc=Nbin_pi;
                        for(kc=1;kc<=jack_samp;kc++) {
                        	if(jack1[i] != kc && jack2[j] != kc) Npair[ic][jc][kc]=Npair[ic][jc][kc]+1;
                                }                     
			}
				}
			}
		}
      }
     }
}


int main(int argc, char *argv[])
{
 int    i,j,k,auto_opt;
 double rp_min,rp_max,pi_min,pi_max;
 double Lx,Ly,Lz;
 double *rp_cen,*pi_cen;
 int    Nbin_rp,Nbin_pi,jack_samp;

 int *jack1,*jack2,**galIND_1, **cellNEI,**galIND_2,*begin_1,*end_1,*begin_2,*end_2;
 int jacktmp,begintmp,endtmp,cn1tmp,cn2tmp,cn3tmp,cn4tmp,cn5tmp,cn6tmp,cn7tmp,cn8tmp,cn9tmp,cn10tmp,cn11tmp,cn12tmp,cn13tmp,cn14tmp,cn15tmp,cn16tmp,cn17tmp,cn18tmp,cn19tmp,cn20tmp,cn21tmp,cn22tmp,cn23tmp,cn24tmp,cn25tmp,cn26tmp,cn27tmp,*cn1,*cn2,*cn3,*cn4,*cn5,*cn6,*cn7,*cn8,*cn9,*cn10,*cn11,*cn12,*cn13,*cn14,*cn15,*cn16,*cn17,*cn18,*cn19,*cn20,*cn21,*cn22,*cn23,*cn24,*cn25,*cn26,*cn27;
 double xtmp,ytmp,ztmp;
 double *x1,*y1,*z1,*x2,*y2,*z2;
 int    N1,N2,N3,N4,N5;

 double Npair_tot;
 double ***Npair;

 char   file1[128],file2[128],file3[128],file4[128],file5[128],filename[128];
 FILE   *fp1,*fp2,*fp3,*fp4,*fp5,*fp;

 int    Nx,Ny;
 float  ftmp;

 if(argc!=17) { 
   printf("Count pairs ...\n"
          "Usage: %s catalog1 catalog2 galIND_file1 galIND_file2 cellNEI_file1 s_min(ln) s_max(ln) u_min (linear) u_max(linear) "
          "Nbin_s Nbin_u Lx Ly Lz Fileout.dat #jack samples\n",argv[0]);
   exit(0);
 }

 strcpy(file1,argv[1]);
 strcpy(file2,argv[2]); 
 strcpy(file3,argv[3]);
 strcpy(file4,argv[4]);
 strcpy(file5,argv[5]);
 rp_min  =atof(argv[6]);
 rp_max  =atof(argv[7]);
 pi_min  =atof(argv[8]);
 pi_max  =atof(argv[9]);
 Nbin_rp =atoi(argv[10]);
 Nbin_pi =atoi(argv[11]);
 Lx     =atof(argv[12]);
 Ly     =atof(argv[13]);
 Lz     =atof(argv[14]);
 jack_samp =atof(argv[16]);

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
    fscanf(fp1,"%lf %lf %lf %d",&xtmp,&ytmp,&ztmp,&jacktmp);
    N1=0;
    while(!feof(fp1)) {
          N1=N1+1;
          fscanf(fp1,"%lf %lf %lf %d",&xtmp,&ytmp,&ztmp,&jacktmp);
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

    fscanf(fp1,"%lf %lf %lf %d",&xtmp,&ytmp,&ztmp,&jacktmp);
    N1=0;
    while(!feof(fp1)) { 
            N1=N1+1;
            fscanf(fp1,"%lf %lf %lf %d",&xtmp,&ytmp,&ztmp,&jacktmp);
    }
    rewind(fp1);

    fscanf(fp2,"%lf %lf %lf %d",&xtmp,&ytmp,&ztmp,&jacktmp);
    N2=0;
    while(!feof(fp2)) {
            N2=N2+1;
            fscanf(fp2,"%lf %lf %lf %d",&xtmp,&ytmp,&ztmp,&jacktmp);
    }
    rewind(fp2);
 }
	//**Read in tree cell files**//
    
    if(strcmp(file3,file4)==0) {
    if(!(fp3=fopen(file3,"rt"))) {
       printf("cannot find %s\n",file3);
       exit(0);
    }
      fscanf(fp3,"%d %d",&begintmp,&endtmp);
      N3=0;
      while(!feof(fp3)) {
            N3=N3+1;
            fscanf(fp3,"%d %d",&begintmp,&endtmp);
      }
      rewind(fp3);
      N4=N3;
   }
 else {
    if(!(fp3=fopen(file3,"rt"))) {
         printf("cannot find %s\n",file3);
         exit(0);
    }
    if(!(fp4=fopen(file4,"rt"))) {
         printf("cannot find %s\n",file4);
         exit(0);
    }

      fscanf(fp3,"%d %d",&begintmp,&endtmp);
      N3=0;
      while(!feof(fp3)) {
            N3=N3+1;
            fscanf(fp3,"%d %d",&begintmp,&endtmp);
      }
      rewind(fp3);

      fscanf(fp4,"%d %d",&begintmp,&endtmp);
      N4=0;
      while(!feof(fp4)) {
            N4=N4+1;
            fscanf(fp4,"%d %d",&begintmp,&endtmp);
      }
      rewind(fp4);
}
    if(!(fp5=fopen(file5,"rt"))) {
       printf("cannot find %s\n",file5);
       exit(0);
    }
        fscanf(fp5,"%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d", &cn1tmp,&cn2tmp,&cn3tmp,&cn4tmp,&cn5tmp,&cn6tmp,&cn7tmp,&cn8tmp,&cn9tmp,&cn10tmp,&cn11tmp,&cn12tmp,&cn13tmp,&cn14tmp,&cn15tmp,&cn16tmp,&cn17tmp,&cn18tmp,&cn19tmp,&cn20tmp,&cn21tmp,&cn22tmp,&cn23tmp,&cn24tmp,&cn25tmp,&cn26tmp,&cn27tmp);
      N5=0;
      while(!feof(fp5)) {
            N5=N5+1;
        fscanf(fp5,"%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d", &cn1tmp,&cn2tmp,&cn3tmp,&cn4tmp,&cn5tmp,&cn6tmp,&cn7tmp,&cn8tmp,&cn9tmp,&cn10tmp,&cn11tmp,&cn12tmp,&cn13tmp,&cn14tmp,&cn15tmp,&cn16tmp,&cn17tmp,&cn18tmp,&cn19tmp,&cn20tmp,&cn21tmp,&cn22tmp,&cn23tmp,&cn24tmp,&cn25tmp,&cn26tmp,&cn27tmp);
      }
      rewind(fp5);

if (N3!=N5) printf("Catalog 1 does not have same cell number as cell neighbor table.\nN3= %d\nN5= %d\n",N3,N5);
if (N4!=N5) printf("Catalog 2 does not have same cell number as cell neighbor table.\nN4= %d\nN5= %d\n",N4,N5);
    
 rp_cen= vector(1,Nbin_rp);
 pi_cen  = vector(1,Nbin_pi);
 Npair   =f3tensor(1,Nbin_rp,1,Nbin_pi,1,jack_samp);
 x1      = vector(1,N1);
 y1      = vector(1,N1);
 z1      = vector(1,N1);
 x2      = vector(1,N2);
 y2      = vector(1,N2);
 z2      = vector(1,N2);
 jack1    = ivector(1,N1);
 jack2    = ivector(1,N2);
 begin_1   = ivector(1,N3);
 end_1     = ivector(1,N3);
 begin_2   = ivector(1,N4);
 end_2     = ivector(1,N4);
 cn1     = ivector(1,N5);
 cn2 = ivector(1,N5);
 cn3 = ivector(1,N5);
 cn4 = ivector(1,N5);
 cn5 = ivector(1,N5);
 cn6 = ivector(1,N5);
 cn7  = ivector(1,N5);
 cn8 = ivector(1,N5);
 cn9 = ivector(1,N5);
 cn10 = ivector(1,N5);
 cn11 = ivector(1,N5);
 cn12 = ivector(1,N5);
 cn13 = ivector(1,N5);
 cn14 = ivector(1,N5);
 cn15 = ivector(1,N5);
 cn16 = ivector(1,N5);
 cn17 = ivector(1,N5);
 cn18 = ivector(1,N5);
 cn19 = ivector(1,N5);
 cn20 = ivector(1,N5);
 cn21 = ivector(1,N5);
 cn22 = ivector(1,N5);
 cn23 = ivector(1,N5);
 cn24 = ivector(1,N5);
 cn25 = ivector(1,N5);
 cn26 = ivector(1,N5);
 cn27 = ivector(1,N5);
 galIND_1 = imatrix(1,N3,1,2);
 cellNEI = imatrix(1,N5,1,27);
 galIND_2 = imatrix(1,N4,1,2);


 if(auto_opt==1) { 
    for(i=1;i<=N1;i++) { 
        fscanf(fp1,"%lf %lf %lf %d",&(x1[i]),&(y1[i]),&(z1[i]),&(jack1[i]));
        x2[i]=x1[i];
        y2[i]=y1[i];
        z2[i]=z1[i];
	jack2[i]=jack1[i];
    }
    fclose(fp1);

    for(i=1;i<=N3;i++) {
        fscanf(fp3,"%d %d",&(begin_1[i]),&(end_1[i]));
        galIND_1[i][1]=begin_1[i];
        galIND_1[i][2]=end_1[i];
	galIND_2[i][1]=begin_1[i];
	galIND_2[i][2]=end_1[i];
        }
        fclose(fp3); 
	}
 else {
    for(i=1;i<=N1;i++)
         fscanf(fp1,"%lf %lf %lf %d",&(x1[i]),&(y1[i]),&(z1[i]),&(jack1[i]));
    for(i=1;i<=N2;i++)
         fscanf(fp2,"%lf %lf %lf %d",&(x2[i]),&(y2[i]),&(z2[i]),&(jack2[i]));
    fclose(fp1);
    fclose(fp2);
 
	for(i=1;i<=N3;i++) {
		fscanf(fp3,"%d %d",&(begin_1[i]),&(end_1[i]));
		galIND_1[i][1]=begin_1[i];
		galIND_1[i][2]=end_1[i];

     	}
   	fclose(fp3);

        for(i=1;i<=N4;i++) {
                fscanf(fp4,"%d %d",&(begin_2[i]),&(end_2[i]));
                galIND_2[i][1]=begin_2[i];
                galIND_2[i][2]=end_2[i];

        }
        fclose(fp4);
   }
   
     for(i=1;i<=N5;i++) {
                fscanf(fp5,"%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d", &(cn1[i]),&(cn2[i]),&(cn3[i]),&(cn4[i]),&(cn5[i]),&(cn6[i]),&(cn7[i]),&(cn8[i]),&(cn9[i]),&(cn10[i]),&(cn11[i]),&(cn12[i]),&(cn13[i]),&(cn14[i]),&(cn15[i]),&(cn16[i]),&(cn17[i]),&(cn18[i]),&(cn19[i]),&(cn20[i]),&(cn21[i]),&(cn22[i]),&(cn23[i]),&(cn24[i]),&(cn25[i]),&(cn26[i]),&(cn27[i]));
                cellNEI[i][1]=cn1[i];
                cellNEI[i][2]=cn2[i];
                cellNEI[i][3]=cn3[i];
                cellNEI[i][4]=cn4[i];
                cellNEI[i][5]=cn5[i];
                cellNEI[i][6]=cn6[i];
                cellNEI[i][7]=cn7[i];
                cellNEI[i][8]=cn8[i];
                cellNEI[i][9]=cn9[i];
                cellNEI[i][10]=cn10[i];
                cellNEI[i][11]=cn11[i];
                cellNEI[i][12]=cn12[i];
                cellNEI[i][13]=cn13[i];
                cellNEI[i][14]=cn14[i];
                cellNEI[i][15]=cn15[i];
                cellNEI[i][16]=cn16[i];
                cellNEI[i][17]=cn17[i];
                cellNEI[i][18]=cn18[i];
                cellNEI[i][19]=cn19[i];
                cellNEI[i][20]=cn20[i];
                cellNEI[i][21]=cn21[i];
                cellNEI[i][22]=cn22[i];
                cellNEI[i][23]=cn23[i];
                cellNEI[i][24]=cn24[i];
                cellNEI[i][25]=cn25[i];
                cellNEI[i][26]=cn26[i];
                cellNEI[i][27]=cn27[i];
        }
        fclose(fp5);
 

printf("# %s %s ln(r_p)=(%6.2f,%6.2f) r_pi=(%6.2f,%6.2f) %d r_p bins %d r_pi bins L_x=%5.2fh^-1Mpc L_y=%5.2fh^-1Mpc L_z=%5.2fh^-1Mpc %d Halos1 %d Halos2 %d #jack samples\n",file1,file2,rp_min,rp_max,pi_min,pi_max,Nbin_rp,Nbin_pi,Lx,Ly,Lz,N1,N2,jack_samp);

 PairCount(rp_cen, pi_cen,Npair,
           rp_min, rp_max, Nbin_rp,
             pi_min,   pi_max, Nbin_pi,
           x1, y1, z1,jack1, N1,
           x2, y2, z2,jack2, N2, auto_opt,
           Lx, Ly, Lz,jack_samp,cellNEI,galIND_1,galIND_2,N3,N4,N5);

 if(auto_opt==1) Npair_tot=0.5*N1*(N1-1.0);
 else            Npair_tot=(N1+0.0)*(N2+0.0);

 char fileout[128];
 strcpy(fileout,argv[15]);
 if(!(fp=fopen(fileout,"wt"))) FO_error(fileout);

 for (j=1;j<=Nbin_pi;j++) {
   for(i=1;i<=Nbin_rp;i++) {
	for(k=1;k<=jack_samp;k++) {
     fprintf(fp,"%10.6f %10.6f %d %e\n",rp_cen[i],pi_cen[j],k, Npair[i][j][k]);
	}
   }
 }
 fclose(fp);

 free_vector(rp_cen, 1,Nbin_rp);
 free_vector(pi_cen,   1,Nbin_pi);
 free_f3tensor(Npair,1,Nbin_rp,1,Nbin_pi,1,jack_samp);
 free_vector(x1,     1,N1);
 free_vector(y1,     1,N1);
 free_vector(z1,     1,N1);
 free_vector(x2,     1,N2);
 free_vector(y2,     1,N2);
 free_vector(z2,     1,N2);
 free_ivector(jack1,     1,N1);
 free_ivector(jack2,     1,N2);
 free_ivector(begin_1,     1,N3);
 free_ivector(end_1,     1,N3);
 free_ivector(begin_2,     1,N4);
 free_ivector(end_2,     1,N4);
 free_ivector(cn1,     1,N5);
 free_ivector(cn2,     1,N5);
 free_ivector(cn3,     1,N5);
 free_ivector(cn4,     1,N5);
 free_ivector(cn5,     1,N5);
 free_ivector(cn6,     1,N5);
 free_ivector(cn7,     1,N5);
 free_ivector(cn8,     1,N5);
 free_ivector(cn9,     1,N5);
 free_ivector(cn10,     1,N5);
 free_ivector(cn11,     1,N5);
 free_ivector(cn12,     1,N5);
 free_ivector(cn13,     1,N5);
 free_ivector(cn14,     1,N5);
 free_ivector(cn15,     1,N5);
 free_ivector(cn16,     1,N5);
 free_ivector(cn17,     1,N5);
 free_ivector(cn18,     1,N5);
 free_ivector(cn19,     1,N5);
 free_ivector(cn20,     1,N5);
 free_ivector(cn21,     1,N5);
 free_ivector(cn22,     1,N5);
 free_ivector(cn23,     1,N5);
 free_ivector(cn24,     1,N5); 
 free_ivector(cn25,     1,N5);
 free_ivector(cn26,     1,N5);
 free_ivector(cn27,     1,N5);
 free_imatrix(galIND_1,1,N3,1,2);
 free_imatrix(cellNEI,1,N5,1,27);
 free_imatrix(galIND_2,1,N4,1,2);

	printf("Program Finished.");
} 
