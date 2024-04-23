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

void PairCount(double *s_cen, double *u_cen,
	       double **Npair,
               double s_min,  double s_max,  int Nbin_s,
               double   u_min,  double   u_max,  int Nbin_u,
               double *x1, double *y1, double *z1, int N1, 
               double *x2, double *y2, double *z2, int N2, 
               int auto_opt,
               double Lx,  double Ly,  double Lz, int **cellNEI, int **galIND_1,int **galIND_2, int N3, int N4, int N5) 
{
 int    i,ic,j,jc,j0,ii,k;
 double ds,du,rpi,dx,dy,dz,d1x,d1y,d1z,d2x,d2y,d2z,s,u,l,lx,ly,lz,sx,sy,sz;
 double slnbinwidth,s_lo,*s_hi,dv,navg1,navg2;
 double S,V,factor;

 
 if(auto_opt!=0 && N1!=N2) { 
   printf("PairCount> Error, N1!=N2 for option auto\n");
   exit(0);
 }

 ds=log10(pow(10,s_max)/pow(10,s_min))/Nbin_s;
 du=(u_max-u_min)/Nbin_u;
 s_hi=vector(1,Nbin_s);

 for(i=1;i<=Nbin_s;i++) s_cen[i]=pow(10,(s_min+(i-0.5)*ds));
 for(i=1;i<=Nbin_u;i++) u_cen[i]=u_min+(i-0.5)*du;
 for(i=1;i<=Nbin_s;i++)
    for(j=1;j<=Nbin_u;j++){ 
       Npair[i][j]=0.0; 
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
       			u=fabs((lx*sx+ly*sy+lz*sz)/(s*l));    

			if(log10(s)>=s_min && log10(s)<=s_max && u>=u_min && u<=u_max) {
				ic=floor(log10(s/pow(10,s_min))/ds  +1.0);
				jc=floor((u-u_min)/du+1.0);
				if (ic<1) ic=1;
				else if(ic>Nbin_s) ic=Nbin_s;
				if (jc<1) jc=1;
				else if (jc>Nbin_u) jc=Nbin_u;
        			Npair[ic][jc]=Npair[ic][jc]+1;
				}
			}
		}}
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
                        u=fabs((lx*sx+ly*sy+lz*sz)/(s*l));

                        if(log10(s)>=s_min && log10(s)<=s_max && u>=u_min && u<=u_max) {
                                ic=floor(log10(s/pow(10,s_min))/ds  +1.0);
                                jc=floor((u-u_min)/du+1.0);
                                if (ic<1) ic=1;
                                else if(ic>Nbin_s) ic=Nbin_s;
                                if (jc<1) jc=1;
                                else if (jc>Nbin_u) jc=Nbin_u;
                                Npair[ic][jc]=Npair[ic][jc]+1;
				}
			}
		}}
      }
     }
   }
   //	if(auto_opt==1) {
//		for(ic=1;ic<=Nbin_s;ic++) {
//			for(jc=1;jc<=Nbin_u;jc++){
  //                         Npair[ic][jc]=0.5*Npair[ic][jc];
//		}
//		}
//		}


int main(int argc, char *argv[])
{
 int    i,j,auto_opt;
 double s_min,s_max,u_min,u_max;
 double Lx,Ly,Lz;
 double *s_cen,*u_cen;
 int    Nbin_s,Nbin_u;

 int **galIND_1, **cellNEI,**galIND_2,*begin_1,*end_1,*begin_2,*end_2;
 int begintmp,endtmp,cn1tmp,cn2tmp,cn3tmp,cn4tmp,cn5tmp,cn6tmp,cn7tmp,cn8tmp,cn9tmp,cn10tmp,cn11tmp,cn12tmp,cn13tmp,cn14tmp,cn15tmp,cn16tmp,cn17tmp,cn18tmp,cn19tmp,cn20tmp,cn21tmp,cn22tmp,cn23tmp,cn24tmp,cn25tmp,cn26tmp,cn27tmp,*cn1,*cn2,*cn3,*cn4,*cn5,*cn6,*cn7,*cn8,*cn9,*cn10,*cn11,*cn12,*cn13,*cn14,*cn15,*cn16,*cn17,*cn18,*cn19,*cn20,*cn21,*cn22,*cn23,*cn24,*cn25,*cn26,*cn27;
 double xtmp,ytmp,ztmp;
 double *x1,*y1,*z1,*x2,*y2,*z2;
 int    N1,N2,N3,N4,N5;

 double Npair_tot,**Npair;

 char   file1[128],file2[128],file3[128],file4[128],file5[128],filename[128];
 FILE   *fp1,*fp2,*fp3,*fp4,*fp5,*fp;

 int    Nx,Ny;
 float  ftmp;

 if(argc!=16) { 
   printf("Count pairs ...\n"
          "Usage: %s catalog1 catalog2 galIND_file1 galIND_file2 cellNEI_file1 s_min(ln) s_max(ln) u_min (linear) u_max(linear) "
          "Nbin_s Nbin_u Lx Ly Lz Fileout.dat\n",argv[0]);
   exit(0);
 }

 strcpy(file1,argv[1]);
 strcpy(file2,argv[2]); 
 strcpy(file3,argv[3]);
 strcpy(file4,argv[4]);
 strcpy(file5,argv[5]);
 s_min  =atof(argv[6]);
 s_max  =atof(argv[7]);
 u_min  =atof(argv[8]);
 u_max  =atof(argv[9]);
 Nbin_s =atoi(argv[10]);
 Nbin_u =atoi(argv[11]);
 Lx     =atof(argv[12]);
 Ly     =atof(argv[13]);
 Lz     =atof(argv[14]);

 if(s_min>=s_max) {
    printf("Well, s_min < s_max, please.\n");
    exit(0);
 }
 if(u_min>=u_max) {
    printf("Well, u_min < u_max, please.\n");
    exit(0);
 }

 if(Nbin_s<=1 || Nbin_u<=1) {
    printf("Well, Nbin_s and Nbin_u > 1, please.\n");
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
    
 s_cen= vector(1,Nbin_s);
 u_cen  = vector(1,Nbin_u);
 Npair   =matrix(1,Nbin_s,1,Nbin_u);
 x1      = vector(1,N1);
 y1      = vector(1,N1);
 z1      = vector(1,N1);
 x2      = vector(1,N2);
 y2      = vector(1,N2);
 z2      = vector(1,N2);
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
        fscanf(fp1,"%lf %lf %lf",&(x1[i]),&(y1[i]),&(z1[i]));
        x2[i]=x1[i];
        y2[i]=y1[i];
        z2[i]=z1[i];
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
         fscanf(fp1,"%lf %lf %lf",&(x1[i]),&(y1[i]),&(z1[i]));
    for(i=1;i<=N2;i++)
         fscanf(fp2,"%lf %lf %lf",&(x2[i]),&(y2[i]),&(z2[i]));
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
 

printf("# %s %s ln(s)=(%6.2f,%6.2f) u=(%6.2f,%6.2f) %d s bins %d u bins L_x=%5.2fh^-1Mpc L_y=%5.2fh^-1Mpc L_z=%5.2fh^-1Mpc %d Halos\n",file1,file2,s_min,s_max,u_min,u_max,Nbin_s,Nbin_u,Lx,Ly,Lz,N2);

 PairCount(s_cen, u_cen,Npair,
           s_min, s_max, Nbin_s,
             u_min,   u_max, Nbin_u,
           x1, y1, z1, N1,
           x2, y2, z2, N2, auto_opt,
           Lx, Ly, Lz,cellNEI,galIND_1,galIND_2,N3,N4,N5);

 if(auto_opt==1) Npair_tot=0.5*N1*(N1-1.0);
 else            Npair_tot=(N1+0.0)*(N2+0.0);

 char fileout[128];
 strcpy(fileout,argv[15]);
 if(!(fp=fopen(fileout,"wt"))) FO_error(fileout);

 for (j=1;j<=Nbin_u;j++)
   for(i=1;i<=Nbin_s;i++)
     fprintf(fp,"%10.6f %10.6f %e\n",s_cen[i],u_cen[j],Npair[i][j]);
 fclose(fp);

 free_vector(s_cen, 1,Nbin_s);
 free_vector(u_cen,   1,Nbin_u);
 free_matrix(Npair,1,Nbin_s,1,Nbin_u);
 free_vector(x1,     1,N1);
 free_vector(y1,     1,N1);
 free_vector(z1,     1,N1);
 free_vector(x2,     1,N2);
 free_vector(y2,     1,N2);
 free_vector(z2,     1,N2);
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
