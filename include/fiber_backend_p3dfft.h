/*
* ---------------
* P3DFFT backend
* ---------------
*/
#ifndef FIBER_BACKEND_P3DFFT_H
#define FIBER_BACKEND_P3DFFT_H

#include <stdio.h>
#include <p3dfft.h>

//=====================  Complex-to-Complex transform =========================

void compute_z2z_p3dfft( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
			 void const *in, void *out, double *timer, int const pgrid_in[3], int const pgrid_out[3])
{

  int i,j;
  int dmap1[3],dmap2[3],mo1[3],mo2[3],type_ids1[3],type_ids2[3],type_forward,type_backward;
  int gdims[3],gdims2[3],ldims1[3];

  for(i=0;i<3;i++) {
    dmap1[i] = mo1[i] = i;
    for(j=0;j<3;j++)
      if(pgrid_in[i] == pgrid_out[j])
	dmap2[j] = mo2[j] = i;
  }

  p3dfft_setup();

  type_ids1[0] = P3DFFT_CFFT_FORWARD_D;
  type_ids1[1] = P3DFFT_CFFT_FORWARD_D;
  type_ids1[2] = P3DFFT_CFFT_FORWARD_D;

  type_forward = p3dfft_init_3Dtype(type_ids1);

  ldims1[0] = inbox_high[0]-inbox_low[0]+1;
  ldims1[1] = inbox_high[1]-inbox_low[1]+1;
  ldims1[2] = inbox_high[2]-inbox_low[2]+1;

  for(i=0;i<3;i++)
    gdims[i] = ldims1[i] * pgrid_in[i];

  int Pgrid = p3dfft_init_proc_grid(pgrid_in,comm);

  Grid *Xpencil = p3dfft_init_data_grid(gdims,-1,Pgrid,dmap1,mo1);

  Grid *Zpencil = p3dfft_init_data_grid(gdims,-1,Pgrid,dmap2,mo2);


    // Plan creation ...
    // void *plan;

    MPI_Barrier(comm);
    timer[0] = -MPI_Wtime();

    // plan create
    Plan3D trans_f = p3dfft_plan_3Dtrans(Xpencil,Zpencil,type_forward);

    MPI_Barrier(comm);
    timer[0] += MPI_Wtime();



    // FFT execution ...
    MPI_Barrier(comm);
    timer[1] = -MPI_Wtime();

    // compute FFT 
    p3dfft_exec_3Dtrans_double(trans_f,in,out,0);

    MPI_Barrier(comm);
    timer[1] = +MPI_Wtime();

    // Delete plan ...

    //    printf("P3DFFT: Get Complex-to-Complex using sample/C/test3d_c2c.c  \n");
    //MPI_Abort(comm, 1);
}

//=====================  Real-to-Complex transform =========================


void compute_d2z_p3dfft( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
			 double const *in, void *out, double *timer,int pgrid_in[3],int pgrid_out[3])
{

  int i,j,d;
  int dmap1[3],dmap2[3],mo1[3],mo2[3],type_ids[3],type_forward;
  int gdims[3],gdims2[3],ldims1[3],trans_b,Pgrid;
  

  for(i=0;i<3;i++) {
    dmap1[i] = mo1[i] = i;
    for(j=0;j<3;j++)
      if(pgrid_in[i] == pgrid_out[j])
	dmap2[j] = mo2[j] = i;
  }

  p3dfft_setup();

  for(i=0;i<3;i++)
    if(pgrid_in[i] == 1) {
      d = i;
      break;
    }

  type_ids[0] = P3DFFT_CFFT_FORWARD_D;
  type_ids[1] = P3DFFT_CFFT_FORWARD_D;
  type_ids[2] = P3DFFT_CFFT_FORWARD_D;
  type_ids[d] = P3DFFT_R2CFFT_D;

  type_forward = p3dfft_init_3Dtype(type_ids);

  ldims1[0] = inbox_high[0]-inbox_low[0]+1;
  ldims1[1] = inbox_high[1]-inbox_low[1]+1;
  ldims1[2] = inbox_high[2]-inbox_low[2]+1;

  for(i=0;i<3;i++)
    gdims[i] = gdims2[i] = ldims1[i] * pgrid_in[i];

  Pgrid = p3dfft_init_proc_grid(pgrid_in,comm);

  Grid *Xpencil = p3dfft_init_data_grid(gdims,-1,Pgrid,dmap1,mo1);
  
  gdims2[d] = gdims2[d]/2+1;
  Grid *Zpencil = p3dfft_init_data_grid(gdims2,d,Pgrid,dmap2,mo2);

    // Plan creation ...
    // void *plan;

  MPI_Barrier(comm);
  timer[0] = -MPI_Wtime();
  
  // plan create
  Plan3D trans_f = p3dfft_plan_3Dtrans(Xpencil,Zpencil,type_forward);
  
  MPI_Barrier(comm);
  timer[0] += MPI_Wtime();
  
  
  
  // FFT execution ...
  MPI_Barrier(comm);
  timer[1] = -MPI_Wtime();
  
  // compute FFT 
  p3dfft_exec_3Dtrans_double(trans_f,in,out,0);
  
  MPI_Barrier(comm);
  timer[1] = +MPI_Wtime();

    // Delete plan ...

    //    printf("P3DFFT: Get Real-to-Complex using sample/C/test3d_r2c.c  \n");
    //MPI_Abort(comm, 1);
}

//=====================  Complex-to-Real transform =========================


void compute_z2d_p3dfft( int const inbox_low[3], int const inbox_high[3],
                  int const outbox_low[3], int const outbox_high[3], 
                  MPI_Comm const comm,
			 double const *in, void *out, double *timer,int pgrid_in[3],int pgrid_out[3])
{

  int i,j,d;
  int dmap1[3],dmap2[3],mo1[3],mo2[3],type_ids[3],type_backward;
  int gdims[3],gdims2[3],ldims1[3],ldims2[3],trans_b,Pgrid;


  for(i=0;i<3;i++) {
    dmap1[i] = mo1[i] = i;
    for(j=0;j<3;j++)
      if(pgrid_in[i] == pgrid_out[j])
	dmap2[j] = mo2[j] = i;
  }

  p3dfft_setup();

  for(i=0;i<3;i++)
    if(pgrid_out[i] == 1) {
      d = i;
      break;
    }

  type_ids[0] = P3DFFT_CFFT_FORWARD_D;
  type_ids[1] = P3DFFT_CFFT_FORWARD_D;
  type_ids[2] = P3DFFT_CFFT_FORWARD_D;
  type_ids[d] = P3DFFT_C2RFFT_D;

  type_backward = p3dfft_init_3Dtype(type_ids);

  ldims2[0] = outbox_high[0]-outbox_low[0]+1;
  ldims2[1] = outbox_high[1]-outbox_low[1]+1;
  ldims2[2] = outbox_high[2]-outbox_low[2]+1;

  for(i=0;i<3;i++)
    gdims[i] = gdims2[i] = ldims1[i] * pgrid_out[i];
  gdims[d] = gdims2[d]/2+1;
  Pgrid = p3dfft_init_proc_grid(pgrid_in,comm);

  Grid *Zpencil = p3dfft_init_data_grid(gdims,-1,Pgrid,dmap1,mo1);
  
  Grid *Xpencil = p3dfft_init_data_grid(gdims2,d,Pgrid,dmap2,mo2);

    // Plan creation ...
    // void *plan;

  MPI_Barrier(comm);
  timer[0] = -MPI_Wtime();
  
    // plan create
  trans_b = p3dfft_plan_3Dtrans(Zpencil,Xpencil,type_backward);
  
  MPI_Barrier(comm);
  timer[0] += MPI_Wtime();
  
  
  
  // FFT execution ...
  MPI_Barrier(comm);
  timer[1] = -MPI_Wtime();
  
  // compute FFT 
  p3dfft_exec_3Dtrans_double(trans_b,in,out,0);
  
  MPI_Barrier(comm);
  timer[1] = +MPI_Wtime();

  // Delete plan ...

    //    printf("P3DFFT: Get Complex-to-Real using sample/C/test3d_r2c.c  \n");
    //MPI_Abort(comm, 1);
}




#endif  //! FIBER_BACKEND_P3DFFT_H
