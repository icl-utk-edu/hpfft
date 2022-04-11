#ifndef FIBER_UTILS_H
#define FIBER_UTILS_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct{
    double r;
    double i;
} fiber_complex;

// constants
  const char *syntax =
    "Example of correct syntax: test3D_C2C -lib heffte   -backend fftw    -size nx ny nz      -pgrid P Q    -iter Niter \n"
    "                                      -comm a2a/a2av/a2aw/p2p    -reshape pencil/slabs/bricks   -trace   \n";

enum {
 option_fft_op,         // 0
 option_nx,             // 1
 option_ny,             // 2
 option_nz,             // 3
 option_backend,        // 4
 option_grid_p,         // 5
 option_grid_q,         // 6
 option_grid_r,         // 7
 option_niter,          // 8
 option_A,              // 9,  unasigned
 option_B,              // 10, unasigned
 option_C,              // 11, unasigned
 option_D,              // 12, unasigned
 option_physical,       // 13
};

const int n_options = 14;

#define BIG 1.0e20
int factors[60];
int n_iterations = 10;

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


// Command line parsing
void fiber_parse_options(int argc, char** argv, int * backend_options, char * lib_name, char * lib_1d_backend){
  int iarg = 1;
  while (iarg < argc) {
    if (strcmp(argv[iarg],"-h") == 0) {
      error_all(syntax);
    } else if (strcmp(argv[iarg],"-lib") == 0) {
      if (iarg+2 > argc) error_all(syntax);
      strcpy (lib_name, argv[iarg+1]);
      iarg += 2;
    } else if (strcmp(argv[iarg],"-backend") == 0) {
      if (iarg+2 > argc) error_all(syntax);
      strcpy (lib_1d_backend, argv[iarg+1]);
      iarg += 2;
    } else if (strcmp(argv[iarg],"-size") == 0) {
      if (iarg+4 > argc) error_all(syntax);
      backend_options[1] = atoi(argv[iarg+1]);
      backend_options[2] = atoi(argv[iarg+2]);
      backend_options[3] = atoi(argv[iarg+3]);
      iarg += 4;
    } else if (strcmp(argv[iarg],"-pgrid") == 0) {
      if (iarg+3 > argc) error_all(syntax);
      backend_options[5] = atoi(argv[iarg+1]);
      backend_options[6] = atoi(argv[iarg+2]);
      iarg += 3;
    } else if (strcmp(argv[iarg],"-iter") == 0) {
      if (iarg+2 > argc) error_all(syntax);
      backend_options[8] = atoi(argv[iarg+1]);
      iarg += 2;
    } else if (strcmp(argv[iarg],"-comm") == 0) {
      if (iarg+2 > argc) error_all(syntax);
      if (strcmp(argv[iarg+1],"a2a") == 0) backend_options[10] = 0; 
      else if (strcmp(argv[iarg+1],"a2av") == 0) backend_options[10] = 1;
      else if (strcmp(argv[iarg+1],"a2aw") == 0) backend_options[10] = 2;
      else if (strcmp(argv[iarg+1],"p2p")  == 0) backend_options[10] = 3;
      else if (strcmp(argv[iarg+1],"p2p_pl")  == 0) backend_options[10] = 4;
      else error_all(syntax);
      iarg += 2;
    } else if (strcmp(argv[iarg],"-reshape") == 0) {
      if (iarg+2 > argc) error_all(syntax);
      if (strcmp(argv[iarg+1],"slabs") == 0) backend_options[11] = 1; 
      else if (strcmp(argv[iarg+1],"pencils") == 0) backend_options[11] = 2; 
      else if (strcmp(argv[iarg+1],"bricks") == 0) backend_options[11] = 3; 
      else error_all(syntax);
      iarg += 2;
    } else if (strcmp(argv[iarg],"-trace") == 0) {
      backend_options[15] = 1;
      iarg += 1;
    } else error_all(syntax);
  }
}


#if defined(FIBER_ENABLE_CUDA)
  #include <cufft.h>
  #define fiber_copy_cpu2gpu(h_, g_, size_)  cudaMemcpy((g_), (h_), (size_), cudaMemcpyHostToDevice)
  #define fiber_copy_gpu2cpu(g_, h_, size_)  cudaMemcpy((h_), (g_), (size_), cudaMemcpyDeviceToHost)
#else 
  #define fiber_copy_cpu2gpu(h_, g_, size_)
  #define fiber_copy_gpu2cpu(h_, g_, size_)
#endif    

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
      else if (strcmp(backend,"fftwpp") == 0)
         return 10;                  
      else
         printf("Invalid Backend \n" );
}


int fiber_get_1d_backend(char * backend)
{
      if (strcmp(backend,"stock") == 0)
         return 0;
      else if (strcmp(backend,"fftw") == 0)
         return 1;
      else if (strcmp(backend,"mkl") == 0)
         return 2;
      else if (strcmp(backend,"cufft") == 0)
         return 3;
      else if (strcmp(backend,"rocfft") == 0)
         return 4;
      else if (strcmp(backend,"onemkl") == 0)
         return 5;
      else if (strcmp(backend,"kiss") == 0)
         return 6;         
      else
         printf("Invalid 1-D Backend \n" );
}



#endif//! FIBER_UTILS_H
