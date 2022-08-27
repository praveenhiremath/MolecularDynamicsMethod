/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   On February 9, 2015, added 'CrackAniso' style for body centred cubic materials by Praveenkumar Hiremath (TU Bergakademie Freiberg | Computational Material Science (Years 2013-15) | www.tu-freiberg.de |)
------------------------------------------------------------------------- */
#include <iostream>
#include "lmptype.h"
#include "mpi.h"
#include "stdlib.h"
#include "string.h"
#include "displace_atoms.h"
#include "atom.h"
#include "modify.h"
#include "domain.h"
#include "lattice.h"
#include "comm.h"
#include "irregular.h"
#include "group.h"
#include "math_const.h"
#include "random_park.h"
#include "error.h"
#include "force.h"
#include <cmath>


using namespace std;
using namespace LAMMPS_NS;
using namespace MathConst;

enum{MOVE,RAMP,RANDOM,ROTATE,CRACKANISO};  //CRACKANISO style is new here

/* ---------------------------------------------------------------------- */

DisplaceAtoms::DisplaceAtoms(LAMMPS *lmp) : Pointers(lmp) {}

/* ---------------------------------------------------------------------- */

void DisplaceAtoms::command(int narg, char **arg)
{
  int i;

  if (domain->box_exist == 0)
    error->all(FLERR,"Displace_atoms command before simulation box is defined");
  if (narg < 2) error->all(FLERR,"Illegal displace_atoms command");
  if (modify->nfix_restart_peratom)
    error->all(FLERR,"Cannot displace_atoms after "
               "reading restart file with per-atom info");

  if (comm->me == 0 && screen) fprintf(screen,"Displacing atoms ...\n");

  // group and style

  int igroup = group->find(arg[0]);
  if (igroup == -1) error->all(FLERR,"Could not find displace_atoms group ID");
  int groupbit = group->bitmask[igroup];

  int style=-1;
  if (strcmp(arg[1],"move") == 0) style = MOVE;
  else if (strcmp(arg[1],"ramp") == 0) style = RAMP;
  else if (strcmp(arg[1],"random") == 0) style = RANDOM;
  else if (strcmp(arg[1],"rotate") == 0) style = ROTATE;
  else if (strcmp(arg[1],"CrackAniso") == 0) style = CRACKANISO;   // Added by Praveenkumar Hiremath
  else error->all(FLERR,"Illegal displace_atoms command");

  // set option defaults

  scaleflag = 1;

  // read options from end of input line

  if (style == MOVE) options(narg-5,&arg[5]);
  else if (style == RAMP) options(narg-8,&arg[8]);
  else if (style == RANDOM) options(narg-6,&arg[6]);
  else if (style == ROTATE) options(narg-9,&arg[9]);
  else if (style == CRACKANISO) options(narg-15,&arg[15]);     // Added by Praveenkumar Hiremath

  // setup scaling

  double xscale,yscale,zscale;
  if (scaleflag) {
    xscale = domain->lattice->xlattice;
    yscale = domain->lattice->ylattice;
    zscale = domain->lattice->zlattice;
  }
  else xscale = yscale = zscale = 1.0;

  // move atoms by 3-vector

  if (style == MOVE) {

    double delx = xscale*force->numeric(FLERR,arg[2]);
    double dely = yscale*force->numeric(FLERR,arg[3]);
    double delz = zscale*force->numeric(FLERR,arg[4]);

    double **x = atom->x;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        x[i][0] += delx;
        x[i][1] += dely;
        x[i][2] += delz;
      }
    }
  }


// displace atoms by the amounts ux,uy in X and Y directions respectively. This block of code for 'CrackAniso' style is added by Praveenkumar Hiremath.
  if (style == CRACKANISO) {
double K,r,t,lamda,rho,d,n,m,s11,s12,s44,o,p,q,u,v,w,s11_New,s12_New,s22_New,s16_New,s26_New,s66_New;
double A[6];
int j;

//inputs:Algorithm:Step 13
s11=force->numeric(FLERR,arg[5]);  //Compliance constant s11
s12=force->numeric(FLERR,arg[6]);  //Compliance constant s12
s44=force->numeric(FLERR,arg[7]);  //Compliance constant s44

K=force->numeric(FLERR,arg[2]);    //Stress intensity factor

o=force->numeric(FLERR,arg[12]);   //x component of Crack ropagation direction
p=force->numeric(FLERR,arg[13]);   //y component of Crack ropagation direction
q=force->numeric(FLERR,arg[14]);   //z component of Crack ropagation direction

u=force->numeric(FLERR,arg[9]);    //x component of Crack plane direction
v=force->numeric(FLERR,arg[10]);   //y component of Crack plane direction
w=force->numeric(FLERR,arg[11]);   //z component of Crack plane direction

//Function for calculating compliance constants in rotated system with plane strain condition
rotation(o,p,q,u,v,w,s11,s12,s44,s11_New,s12_New,s22_New,s16_New,s26_New,s66_New,A);

cout << "s11=" << A[0] << "\n" << "s12=" << A[1] << "\n" << "s22=" << A[2] << "\n" <<  "s16=" << A[3] << "\n" << "s26=" << A[4] << "\n" << "s66=" << A[5] <<"\n";

//Re-assignment of Compliance constants stored in Array 
s11_New=A[0];     
s12_New=A[1];
s22_New=A[2];
s16_New=A[3];
s26_New=A[4];
s66_New=A[5];

double e=2*sqrt(s11_New*s22_New);

double f=(2*s12_New)+s66_New;
//Calculating λ and ρ 
lamda=s11_New/s22_New;   rho=f/(e);

d=pow(lamda,(-0.25));   n=sqrt((1+rho)*0.5);    m=sqrt((abs(1-rho)*0.5));

j=1;

//Calculating μ1 and μ2
double mu1[2],mu2[2];
//1<rho<infinity    [0]->Indicates real and [1]->Indicates imaginary
if(rho>1)                  //(a)->Report
{
mu1[1]=(j*d)*(n+m);
mu2[1]=j*d*(n-m);
mu1[0]=0;
mu2[0]=0;
cout<<"region of rho rho>1" << "\n";
}

//-1<rho<1
if(rho>-1 && rho<1)        //(b)->Report
{
mu1[0]=m*d;
mu1[1]=j*n*d;
mu2[0]=(-m)*d;
mu2[1]=j*n*d;
cout<<"region of rho -1<rho<1" << "\n";
}

//Isotropic limit
if(rho==1)                 //(c)->Report
{
mu1[1]=j*d;
mu1[0]=0;

mu2[1]=j*d;
mu2[0]=0;
cout<<"region of rho rho=1" << "\n";
}

//Calculating p1,p2,q1,q2
double temp1[2],temp2[2],temp3[2],temp4[2];
complex_square(temp1,mu1);

double p1[2];
p1[0]=(s11_New*temp1[0])+s12_New-(s16_New*mu1[0]);
p1[1]=(s11_New*temp1[1])-(s16_New*mu1[1]);


complex_square(temp2,mu2);
double p2[2];
p2[0]=(s11_New*temp2[0])+s12_New-(s16_New*mu2[0]);
p2[1]=(s11_New*temp2[1])-(s16_New*mu2[1]);


//q1=(s12*mu1)+(s22/mu1)-s26
complex_mul_inverse(temp3,mu1);
double q1[2],tempq1[2],Multi1[2];
tempq1[0]=(s12_New*temp1[0])+s22_New-(s26_New*mu1[0]);
tempq1[1]=(s12_New*temp1[1])-(s26_New*mu1[1]);
complex_multi(Multi1,tempq1,temp3);
q1[0]=Multi1[0];
q1[1]=Multi1[1];


//q2=(s12*mu2)+(s22/mu2)-s26                                
complex_mul_inverse(temp4,mu2);
double q2[2],tempq2[2],Multi2[2];
tempq2[0]=(s12_New*temp2[0])+s22_New-(s26_New*mu2[0]);
tempq2[1]=(s12_New*temp2[1])-(s26_New*mu2[1]);
complex_multi(Multi2,tempq2,temp4);
q2[0]=Multi2[0];
q2[1]=Multi2[1];


    double **x = atom->x;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;


for (i = 0; i < nlocal; i++)  {

r=sqrt(pow(x[i][0],2)+pow(x[i][1],2));          //r=square root of(x^2+y^2)
t=atan2(x[i][1],x[i][0]);                       //t=tan inverse(y/x)

// 1/(mu1-mu2)
double temp5[2],temp6[2];
complex_subtract(temp5,mu1,mu2);
complex_mul_inverse(temp6,temp5);

//mu1*p2
double temp7[2];
complex_multi(temp7,mu1,p2);
double cosine_t=x[i][0]/r;
double sine_t=x[i][1]/r;
double trigno1[2],temp8[2];   
trigno1[0]=cos(t)+(mu2[0]*sin(t));
trigno1[1]=mu2[1]*sin(t);
complex_root(temp8,trigno1);

//(1/(mu1-mu2))*(mu1*p2)*(cos(t)+mu2*sin(t))^0.5
double Multi3[2],Multi4[2];
complex_multi(Multi3,temp6,temp7);
complex_multi(Multi4,Multi3,temp8);

// mu2*p1
double temp9[2];
complex_multi(temp9,mu2,p1);

//(cos(t)+mu1*sin(t))^0.5
double trigno2[2],temp10[2];
trigno2[0]=cos(t)+(mu1[0]*sin(t));
trigno2[1]=mu1[1]*sin(t);
complex_root(temp10,trigno2);

//(1/(mu1-mu2))*(mu2*p1)*(cos(t)+mu1*sin(t))^0.5
double Multi5[2],Multi6[2];
complex_multi(Multi5,temp6,temp9);
complex_multi(Multi6,Multi5,temp10);

//Real component in X displacement 
double UX[2];
complex_subtract(UX,Multi4,Multi6);

//..................UY...................

// mu1*q2 
double temp11[2];
complex_multi(temp11,mu1,q2);

//(1/(mu1-mu2))*(mu1*q2)*(cos(t)+mu1*sin(t))^0.5
double Multi7[2],Multi8[2];
complex_multi(Multi7,temp6,temp11);
complex_multi(Multi8,Multi7,temp8);

// mu2*q1
double temp12[2];
complex_multi(temp12,mu2,q1);

//(1/(mu1-m2))*(mu2*q1)*(cos(t)+mu1*sin(t))^0.5
double Multi9[2],Multi10[2];
complex_multi(Multi9,temp6,temp12);
complex_multi(Multi10,Multi9,temp10);

//Real component in Y displacement 
double UY[2];
complex_subtract(UY,Multi8,Multi10);
//Displacement field solutions
      double delx = K*sqrt(2*r/3.142)*(UX[0]);           //(Eq:3)->Report
 
      double dely = K*sqrt(2*r/3.142)*(UY[0]);           //(Eq:4)->Report
      double delz = 0; 

     if (mask[i] & groupbit) {
        x[i][0] += delx;
        x[i][1] += dely;
        x[i][2] += delz;

      }
    }
  }     // Until here added by Praveenkumar Hiremath

  // move atoms in ramped fashion

  if (style == RAMP) {

    int d_dim;
    if (strcmp(arg[2],"x") == 0) d_dim = 0;
    else if (strcmp(arg[2],"y") == 0) d_dim = 1;
    else if (strcmp(arg[2],"z") == 0) d_dim = 2;
    else error->all(FLERR,"Illegal displace_atoms ramp command");

    double d_lo,d_hi;
    if (d_dim == 0) {
      d_lo = xscale*force->numeric(FLERR,arg[3]);
      d_hi = xscale*force->numeric(FLERR,arg[4]);
    } else if (d_dim == 1) {
      d_lo = yscale*force->numeric(FLERR,arg[3]);
      d_hi = yscale*force->numeric(FLERR,arg[4]);
    } else if (d_dim == 2) {
      d_lo = zscale*force->numeric(FLERR,arg[3]);
      d_hi = zscale*force->numeric(FLERR,arg[4]);
    }

    int coord_dim;
    if (strcmp(arg[5],"x") == 0) coord_dim = 0;
    else if (strcmp(arg[5],"y") == 0) coord_dim = 1;
    else if (strcmp(arg[5],"z") == 0) coord_dim = 2;
    else error->all(FLERR,"Illegal displace_atoms ramp command");

    double coord_lo,coord_hi;
    if (coord_dim == 0) {
      coord_lo = xscale*force->numeric(FLERR,arg[6]);
      coord_hi = xscale*force->numeric(FLERR,arg[7]);
    } else if (coord_dim == 1) {
      coord_lo = yscale*force->numeric(FLERR,arg[6]);
      coord_hi = yscale*force->numeric(FLERR,arg[7]);
    } else if (coord_dim == 2) {
      coord_lo = zscale*force->numeric(FLERR,arg[6]);
      coord_hi = zscale*force->numeric(FLERR,arg[7]);
    }

    double **x = atom->x;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    double fraction,dramp;

    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        fraction = (x[i][coord_dim] - coord_lo) / (coord_hi - coord_lo);
        fraction = MAX(fraction,0.0);
        fraction = MIN(fraction,1.0);
        dramp = d_lo + fraction*(d_hi - d_lo);
        x[i][d_dim] += dramp;
      }
    }
  }

  // move atoms randomly
  // makes atom result independent of what proc owns it via random->reset()

  if (style == RANDOM) {
    RanPark *random = new RanPark(lmp,1);

    double dx = xscale*force->numeric(FLERR,arg[2]);
    double dy = yscale*force->numeric(FLERR,arg[3]);
    double dz = zscale*force->numeric(FLERR,arg[4]);
    int seed = force->inumeric(FLERR,arg[5]);
    if (seed <= 0) error->all(FLERR,"Illegal displace_atoms random command");

    double **x = atom->x;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        random->reset(seed,x[i]);
        x[i][0] += dx * 2.0*(random->uniform()-0.5);
        x[i][1] += dy * 2.0*(random->uniform()-0.5);
        x[i][2] += dz * 2.0*(random->uniform()-0.5);
      }
    }

    delete random;
  }
 
  // rotate atoms by right-hand rule by theta around R
  // P = point = vector = point of rotation
  // R = vector = axis of rotation
  // R0 = runit = unit vector for R
  // D = X - P = vector from P to X
  // C = (D dot R0) R0 = projection of atom coord onto R line
  // A = D - C = vector from R line to X
  // B = R0 cross A = vector perp to A in plane of rotation
  // A,B define plane of circular rotation around R line
  // X = P + C + A cos(theta) + B sin(theta)

  if (style == ROTATE) {
    double axis[3],point[3];
    double a[3],b[3],c[3],d[3],disp[3],runit[3];
    
    int dim = domain->dimension;
    point[0] = xscale*force->numeric(FLERR,arg[2]);
    point[1] = yscale*force->numeric(FLERR,arg[3]);
    point[2] = zscale*force->numeric(FLERR,arg[4]);
    axis[0] = force->numeric(FLERR,arg[5]);
    axis[1] = force->numeric(FLERR,arg[6]);
    axis[2] = force->numeric(FLERR,arg[7]);
    double theta = force->numeric(FLERR,arg[8]);
    if (dim == 2 && (axis[0] != 0.0 || axis[1] != 0.0))
      error->all(FLERR,"Invalid displace_atoms rotate axis for 2d");

    double len = sqrt(axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2]);
    if (len == 0.0)
      error->all(FLERR,"Zero length rotation vector with displace_atoms");
    runit[0] = axis[0]/len;
    runit[1] = axis[1]/len;
    runit[2] = axis[2]/len;

    double sine = sin(MY_PI*theta/180.0);
    double cosine = cos(MY_PI*theta/180.0);
    double ddotr;

    double **x = atom->x;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        d[0] = x[i][0] - point[0];
        d[1] = x[i][1] - point[1];
        d[2] = x[i][2] - point[2];
        ddotr = d[0]*runit[0] + d[1]*runit[1] + d[2]*runit[2];
        c[0] = ddotr*runit[0];
        c[1] = ddotr*runit[1];
        c[2] = ddotr*runit[2];
        a[0] = d[0] - c[0];
        a[1] = d[1] - c[1];
        a[2] = d[2] - c[2];
        b[0] = runit[1]*a[2] - runit[2]*a[1];
        b[1] = runit[2]*a[0] - runit[0]*a[2];
        b[2] = runit[0]*a[1] - runit[1]*a[0];
        disp[0] = a[0]*cosine  + b[0]*sine;
        disp[1] = a[1]*cosine  + b[1]*sine;
        disp[2] = a[2]*cosine  + b[2]*sine;
        x[i][0] = point[0] + c[0] + disp[0];
        x[i][1] = point[1] + c[1] + disp[1];
        if (dim == 3) x[i][2] = point[2] + c[2] + disp[2];
      }
    }
  }

  // move atoms back inside simulation box and to new processors
  // use remap() instead of pbc() in case atoms moved a long distance
  // use irregular() in case atoms moved a long distance

  double **x = atom->x;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;
  for (i = 0; i < nlocal; i++) domain->remap(x[i],image[i]);

  if (domain->triclinic) domain->x2lamda(atom->nlocal);
  domain->reset_box();
  Irregular *irregular = new Irregular(lmp);
  irregular->migrate_atoms(1);
  delete irregular;
  if (domain->triclinic) domain->lamda2x(atom->nlocal);

  // check if any atoms were lost

  bigint natoms;
  bigint nblocal = atom->nlocal;
  MPI_Allreduce(&nblocal,&natoms,1,MPI_LMP_BIGINT,MPI_SUM,world);
  if (natoms != atom->natoms && comm->me == 0) {
    char str[128];
    sprintf(str,"Lost atoms via displace_atoms: original " BIGINT_FORMAT
            " current " BIGINT_FORMAT,atom->natoms,natoms);
    error->warning(FLERR,str);
  }
}

/* ----------------------------------------------------------------------
   parse optional parameters at end of displace_atoms input line
------------------------------------------------------------------------- */

void DisplaceAtoms::options(int narg, char **arg)
{
  if (narg < 0) error->all(FLERR,"Illegal displace_atoms command");

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"units") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal displace_atoms command");
      if (strcmp(arg[iarg+1],"box") == 0) scaleflag = 0;
      else if (strcmp(arg[iarg+1],"lattice") == 0) scaleflag = 1;
      else error->all(FLERR,"Illegal displace_atoms command");
      iarg += 2;
    } else error->all(FLERR,"Illegal displace_atoms command");
  }
}


//Function to calculate compliance constants in rotated system with plane strain condition. This block of code is added by Praveenkumar Hiremath
void DisplaceAtoms::rotation(double o,double p,double q,double u,double v,double w,double s11,double s12,double s44,double s11_New,double s12_New,double s22_New,double s16_New,double s26_New,double s66_New,double A[])
{
double S_New[3][3][3][3]={0},T[3][3]={0},S_Old[3][3][3][3]={0};
double a,b,c;
int g,h,m,n,s,t,k,l;

//compliance constants in original system
double s22=s11,s33=s11;
double s13=s12,s21=s12,s23=s12,s31=s12,s32=s12;
double s55=s44,s66=s44;
double s14=0,s15=0,s16=0,s24=0,s25=0,s26=0,s34=0,s35=0,s36=0,s41=0,s42=0,s43=0,s45=0,s46=0;
double s51=0,s52=0,s53=0,s54=0,s56=0,s61=0,s62=0,s63=0,s64=0,s65=0;

//Building 4th order compliance tensor in original system --> (*) condition in Report.
S_Old[0][0][0][0]=s11;
S_Old[0][0][1][1]=s12;
S_Old[0][0][2][2]=s13;
S_Old[0][0][0][1]=s16/2;     S_Old[0][0][1][0]=s16/2;
S_Old[0][1][0][1]=s66/4;     S_Old[0][1][1][0]=s66/4;    S_Old[1][0][0][1]=s66/4;    S_Old[1][0][1][0]=s66/4;
S_Old[0][1][1][1]=s62/2;     S_Old[1][0][1][1]=s62/2;
S_Old[0][1][2][2]=s63/2;     S_Old[1][0][2][2]=s63/2;
S_Old[1][1][1][1]=s22;
S_Old[1][1][2][2]=s23;
S_Old[1][1][0][1]=s26/2;     S_Old[1][1][1][0]=s26/2;
S_Old[2][2][2][2]=s33;
S_Old[1][1][0][0]=s21;
S_Old[2][2][0][0]=s31;
S_Old[2][2][1][1]=s32;
S_Old[1][2][1][2]=s44/4;     S_Old[2][1][2][1]=s44/4;    S_Old[1][2][2][1]=s44/4;    S_Old[2][1][1][2]=s44/4;
S_Old[0][2][0][2]=s55/4;     S_Old[2][0][2][0]=s55/4;    S_Old[0][2][2][0]=s55/4;    S_Old[2][0][0][2]=s55/4;



//to find the third axis(crack front) of the new coordinate system-> Cross Product:
//a=(p*w)-(q*v);  b=(u*q)-(w*o);  c=(o*v)-(p*u);
a=(p*w)-(q*v);  b=(u*q)-(w*o);  c=(o*v)-(p*u);
cout<<"a="<<a<<"\t"<<"b="<<b<<"\t"<<"c="<<c<<"\n";

//Normalization of the above vectors to form Orthonormal basis set:
double X1=o/sqrt(pow(o,2)+pow(p,2)+pow(q,2)); double Y1=p/sqrt(pow(o,2)+pow(p,2)+pow(q,2)); double Z1=q/sqrt(pow(o,2)+pow(p,2)+pow(q,2));
double X3=a/sqrt(pow(a,2)+pow(b,2)+pow(c,2)); double Y3=b/sqrt(pow(a,2)+pow(b,2)+pow(c,2)); double Z3=c/sqrt(pow(a,2)+pow(b,2)+pow(c,2));
double X2=u/sqrt(pow(u,2)+pow(v,2)+pow(w,2)); double Y2=v/sqrt(pow(u,2)+pow(v,2)+pow(w,2)); double Z2=w/sqrt(pow(u,2)+pow(v,2)+pow(w,2));

//Rotation matrix: T(ij)=x'(i).x(j) Here i is associated with the rotated system axes and j with original.
T[0][0]=X1;    T[0][1]=Y1;     T[0][2]=Z1;
T[1][0]=X2;    T[1][1]=Y2;     T[1][2]=Z2;
T[2][0]=X3;    T[2][1]=Y3;     T[2][2]=Z3;

/*
cout << "Rotation matrix" << "\n" ;
for(s=0;s<3;s++)
{ 
for(t=0;t<3;t++)
{ 
cout<< T[s][t] <<"\t" ;
}
cout << "\n";
}*/


s=0; t=0; k=0; l=0;
//compliance constants in rotated coordinate system
for(s=0;s<3;s++)
{ 
 for(t=0;t<3;t++)
 {
  for(k=0;k<3;k++)
  {
   for(l=0;l<3;l++)
   {
for(g=0;g<3;g++)
{ 
 for(h=0;h<3;h++)
 {
  for(m=0;m<3;m++)
  {
   for(n=0;n<3;n++)
   {

  
S_New[s][t][k][l]+= T[s][g]*T[t][h]*S_Old[g][h][m][n]*T[k][m]*T[l][n];  //Rotation


   }
  }
 }
}
   }
  }
 }
}


double sN11,sN12,sN13,sN14,sN15,sN16;
double sN21,sN22,sN23,sN24,sN25,sN26;
double sN31,sN32,sN33,sN34,sN35,sN36;
double sN41,sN42,sN43,sN44,sN45,sN46;
double sN51,sN52,sN53,sN54,sN55,sN56;
double sN61,sN62,sN63,sN64,sN65,sN66;

//No plane stress and No plane strain   (*) condition in Report.
sN11=S_New[0][0][0][0]; sN12=S_New[0][0][1][1]; sN13=S_New[0][0][2][2]; sN14=2*S_New[0][0][1][2]; sN15=2*S_New[0][0][0][2]; sN16=2*S_New[0][0][0][1];
sN21=S_New[1][1][0][0]; sN22=S_New[1][1][1][1]; sN23=S_New[1][1][2][2]; sN24=2*S_New[1][1][1][2]; sN25=2*S_New[1][1][0][2]; sN26=2*S_New[1][1][0][1];
sN31=S_New[2][2][0][0]; sN32=S_New[2][2][1][1]; sN33=S_New[2][2][2][2]; sN34=2*S_New[2][2][1][2]; sN35=2*S_New[2][2][0][2]; sN36=2*S_New[2][2][0][1];
sN41=2*S_New[1][2][0][0]; sN42=2*S_New[1][2][1][1]; sN43=2*S_New[1][2][2][2]; sN44=4*S_New[1][2][1][2]; sN45=4*S_New[1][2][0][2]; sN46=4*S_New[1][2][0][1];
sN51=2*S_New[0][2][0][0]; sN52=2*S_New[0][2][1][1]; sN53=2*S_New[0][2][2][2]; sN54=4*S_New[0][2][1][2]; sN55=4*S_New[0][2][0][2]; sN56=4*S_New[0][2][0][1];
sN61=2*S_New[0][1][0][0]; sN62=2*S_New[0][1][1][1]; sN63=2*S_New[0][1][2][2]; sN64=4*S_New[0][1][1][2]; sN65=4*S_New[0][1][0][2]; sN66=4*S_New[0][1][0][1];

//For plane strain condition

double A1=(1-((sN45*sN45)/(sN44*sN55)));            //(xi)-> Report
double A2=(sN45*sN14)/(sN44*sN55)-(sN15/sN55);      //(xii)-> Report
double A3=(sN45*sN24)/(sN44*sN55)-(sN25/sN55);      //(xiii)-> Report
double A4=(sN45*sN34)/(sN44*sN55)-(sN35/sN55);      //(xiv)-> Report
double A5=(sN45*sN46)/(sN44*sN55)-(sN56/sN55);      //(xv)-> Report


double B1=(1-((sN34*sN34)/(sN44*sN33)));            //(xvi)-> Report 
double B2=(sN43*sN14)/(sN44*sN33)-(sN13/sN33);      //(xvii)-> Report
double B3=(sN43*sN24)/(sN44*sN33)-(sN23/sN33);      //(xviii)-> Report
double B4=(sN43*sN45)/(sN44*sN33)-(sN53/sN33);      //(xix)-> Report
double B5=(sN43*sN46)/(sN44*sN33)-(sN63/sN33);      //(xx)-> Report


double del1=(B1-((B4*A4)/A1));             //(i)-> Report
double del2=(B2+((B4*A2)/A1));             //(ii)-> Report
double del3=(B3+((B4*A3)/A1));             //(iii)-> Report 
double del4=(B5+((B4*A5)/A1));             //(iv)-> Report


double gamma1=(-1/sN44)*(sN14+((A2*sN45)/A1)+((A4*sN45*del2)/(A1*del1))+((sN34*del2)/del1));   //(v)-> Report
double gamma2=(-1/sN44)*(sN24+((A3*sN45)/A1)+((A4*sN45*del3)/(A1*del1))+((sN34*del3)/del1));   //(vi)-> Report
double gamma3=(-1/sN44)*(sN46+((A5*sN45)/A1)+((A4*sN45*del4)/(A1*del1))+((sN34*del4)/del1));   //(vii)-> Report


double K1=((A2/A1)+((A4*del2)/del1));             //(viii)-> Report
double K2=((A3/A1)+((A4*del3)/del1));             //(ix)-> Report
double K3=((A5/A1)+((A4*del4)/del1));             //(x)-> Report


s11_New=sN11+((sN13*del2)/del1)+(sN14*gamma1)+(sN15*K1);             //(1)-> Report
s12_New=sN12+((sN13*del3)/del1)+(sN14*gamma2)+(sN15*K2);             //(2)-> Report
s16_New=sN16+((sN13*del4)/del1)+(sN14*gamma3)+(sN15*K3);             //(3)-> Report
double s21_New=sN21+((sN23*del2)/del1)+(sN24*gamma1)+(sN25*K1);      //(4)-> Report
s22_New=sN22+((sN23*del3)/del1)+(sN24*gamma2)+(sN25*K2);             //(5)-> Report  
s26_New=sN26+((sN23*del4)/del1)+(sN24*gamma3)+(sN25*K3);             //(6)-> Report
s66_New=sN66+((sN63*del4)/del1)+(sN64*gamma3)+(sN65*K3);             //(7)-> Report

//assigning compliance constants to array
A[0]=s11_New;
A[1]=s12_New;
A[2]=s22_New;
A[3]=s16_New;
A[4]=s26_New;
A[5]=s66_New;

return;
}   //until here by Praveenkumar Hiremath

//Function for adding complex numbers. This block of code is added by Praveenkumar Hiremath.
void DisplaceAtoms::complex_add(double sum[],double Cmplx1[],double Cmplx2[])
{
sum[0]=Cmplx1[0]+Cmplx2[0];       
sum[1]=Cmplx1[1]+Cmplx2[1];       

return;
}

//Function for subtracting complex numbers
void DisplaceAtoms::complex_subtract(double sub[],double Cmplx1[],double Cmplx2[])
{
sub[0]=Cmplx1[0]-Cmplx2[0];       
sub[1]=Cmplx1[1]-Cmplx2[1];       

return;
}

//Function for multiplying complex numbers
void DisplaceAtoms::complex_multi(double mul[],double Cmplx1[],double Cmplx2[])
{
mul[0]=(Cmplx1[0]*Cmplx2[0])-(Cmplx1[1]*Cmplx2[1]);       
mul[1]=(Cmplx1[0]*Cmplx2[1])+(Cmplx1[1]*Cmplx2[0]);       

return;
}

//Function for multiplicative inverse of a complex number
void DisplaceAtoms::complex_mul_inverse(double inv[],double Cmplx[])
{
inv[0]=Cmplx[0]/(pow(Cmplx[0],2)+pow(Cmplx[1],2));       
inv[1]=(-Cmplx[1])/(pow(Cmplx[0],2)+pow(Cmplx[1],2));       

return;
}

//Function for squaring a complex number
void DisplaceAtoms::complex_square(double square[],double Cmplx1[])
{
square[0]=(Cmplx1[0]*Cmplx1[0])-(Cmplx1[1]*Cmplx1[1]);       
square[1]=(2*Cmplx1[0]*Cmplx1[1]);                                 

return;
}

//Function for square root of a complex number
void DisplaceAtoms::complex_root(double sq_root[],double Cmplx1[])
{
if(Cmplx1[1]<0)
{
sq_root[0]=sqrt((sqrt(pow(Cmplx1[0],2)+pow(Cmplx1[1],2))+Cmplx1[0])*0.5);
sq_root[1]=-sqrt((sqrt(pow(Cmplx1[0],2)+pow(Cmplx1[1],2))-Cmplx1[0])*0.5);
}
else
{
sq_root[0]=sqrt((sqrt(pow(Cmplx1[0],2)+pow(Cmplx1[1],2))+Cmplx1[0])*0.5);
sq_root[1]=sqrt((sqrt(pow(Cmplx1[0],2)+pow(Cmplx1[1],2))-Cmplx1[0])*0.5);
}

return;
}   // Until here added by Praveenkumar Hiremath































