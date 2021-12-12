#ifndef FIBER_UTILS_H
#define FIBER_UTILS_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(FIBER_ENABLE_CUDA)
  #include <cufft.h>
  #define fiber_copy_cpu2gpu(h_, g_, size_)  cudaMemcpy((g_), (h_), (size_), cudaMemcpyHostToDevice)
  #define fiber_copy_gpu2cpu(g_, h_, size_)  cudaMemcpy((h_), (g_), (size_), cudaMemcpyDeviceToHost)
#else 
  #define fiber_copy_cpu2gpu(h_, g_, size_)
  #define fiber_copy_gpu2cpu(h_, g_, size_)
#endif    

typedef struct{
    double r;
    double i;
} fiber_complex;


#define BIG 1.0e20
int factors[60];
int n_iterations = 10;

double surfarea(int i, int j, int k, int nx, int ny, int nz)
{
  if (nx == 1 && i != 1) return BIG;
  if (ny == 1 && j != 1) return BIG;
  if (nz == 1 && k != 1) return BIG;

  double dx = 1.0*nx/i;
  double dy = 1.0*ny/j;
  double dz = 1.0*nz/k;
  return dx*dy + dy*dz + dx*dz;
}


void proc3d(int nx, int ny, int nz, int *pgrid, int nprocs)
{

  int ipx,ipy,ipz,nremain;
  double boxx,boxy,boxz,surf;
  double xprd = nx;
  double yprd = ny;
  double zprd = nz;

  double bestsurf = 2.0 * (xprd*yprd + yprd*zprd + zprd*xprd);

  ipx = 1;
  int count = 0;
  while (ipx <= nprocs) {
    count++;
    if (nprocs % ipx == 0) {
      nremain = nprocs/ipx;
      ipy = 1;
      while (ipy <= nremain) {
        if (nremain % ipy == 0) {
          ipz = nremain/ipy;
          boxx = xprd/ipx;
          boxy = yprd/ipy;
          boxz = zprd/ipz;

          surf = boxx*boxy + boxy*boxz + boxz*boxx;
          if (surf < bestsurf) {
            bestsurf = surf;
            pgrid[0] = ipx;
            pgrid[1] = ipy;
            pgrid[2] = ipz;
          }
        }
        ipy++;
      }
    }
    ipx++;
  }

  if (pgrid[0]*pgrid[1]*pgrid[2] != nprocs)
    printf("Computed proc grid does not match nprocs \n");
}







void procfactors(int nx, int ny, int nz, int* np, int nprocs, int nfactor)
{
  int i, j, jk, ifac, jfac, kfac;
  double newarea;

  int sqroot = (int) sqrt(nprocs) + 1;
  if (sqroot*sqroot > nprocs) sqroot--;

  double minarea = 2.0*nx*ny + 2.0*ny*nz + 2.0*nx*nz;

  for (i = 0; i < nfactor; i++) {
    ifac = factors[i];
    jk = nprocs/ifac;
    for (j = i; j < nfactor; j++) {
      jfac = factors[j];
      kfac = jk/jfac;
      if (ifac*jfac*kfac != nprocs) continue;
      if (ifac > jfac || jfac > kfac) continue;

      newarea = surfarea(ifac, jfac, kfac, nx, ny, nz);
      if (newarea < minarea) {
        minarea = newarea;
        np[0] = ifac;
        np[1] = jfac;
        np[2] = kfac;
      }

      newarea = surfarea(ifac, kfac, jfac, nx, ny, nz);
      if (newarea < minarea) {
        minarea = newarea;
        np[0] = ifac;
        np[1] = kfac;
        np[2] = jfac;
      }

      newarea = surfarea(jfac,ifac, kfac, nx, ny, nz);
      if (newarea < minarea) {
        minarea = newarea;
        np[0] = jfac;
        np[1] = ifac;
        np[2] = kfac;
      }

      newarea = surfarea(jfac, kfac,ifac, nx, ny, nz);
      if (newarea < minarea) {
        minarea = newarea;
        np[0] = jfac;
        np[1] = kfac;
        np[2] = ifac;
      }

      newarea = surfarea(kfac,ifac, jfac, nx, ny, nz);
      if (newarea < minarea) {
        minarea = newarea;
        np[0] = kfac;
        np[1] = ifac;
        np[2] = jfac;
      }

      newarea = surfarea(kfac, jfac,ifac, nx, ny, nz);
      if (newarea < minarea) {
        minarea = newarea;
        np[0] = kfac;
        np[1] = jfac;
        np[2] = ifac;
      }
    }
  }

}



int factorize(int n)
{
  int sqroot = (int) sqrt(n) + 1;
  if (sqroot*sqroot > n) sqroot--;

  int nfactor = 0;
  for (int i = 1; i <= sqroot; i++) {
    if (n % i) continue;
    if (nfactor < 60) factors[nfactor++] = i;
  }

  return(nfactor);
}

int fiber_get_backend(char * backend)
{
      if (strcmp(backend,"heffte") == 0)
         return 0;
      else if (strcmp(backend,"fftmpi") == 0)
         return 1;
      else if (strcmp(backend,"accfft") == 0)
         return 2;
      else if (strcmp(backend,"p3dfft") == 0)
         return 3;
      else if (strcmp(backend,"ffte") == 0)
         return 4;
      else if (strcmp(backend,"swfft") == 0)
         return 5;
      else if (strcmp(backend,"decomp2d") == 0)
         return 6;
      else if (strcmp(backend,"nb3dfft") == 0)
         return 7;
      else if (strcmp(backend,"fftw") == 0)
         return 8;
      else if (strcmp(backend,"fftadvmpi") == 0)
         return 9;         
      else
         printf("Invalid Backend \n" );
}


// MPI error handling
void error_all(const char *str)
{
  MPI_Barrier(MPI_COMM_WORLD);
  int me;
  MPI_Comm_rank(MPI_COMM_WORLD,&me);
  if (me == 0) printf("ERROR: %s\n",str);
  MPI_Finalize();
  exit(1);
}

void error_one(const char *str)
{
  int me;
  MPI_Comm_rank(MPI_COMM_WORLD,&me);
  printf("ERROR on proc %d: %s\n",me,str);
  MPI_Abort(MPI_COMM_WORLD,1);
}


#endif//! FIBER_UTILS_H
