/*
* ---------------
* FFTE backend
* ---------------
*/
#ifndef FIBER_BACKEND_FFTE_H
#define FIBER_BACKEND_FFTE_H

#include <stdio.h>

#if defined(FIBER_ENABLE_FFTE)

extern void pzfft3d_(void const *in, void *out, int *nx, int *ny, int *nz,
    MPI_Fint fcomm, int *comm_size, int *opt);
extern void pdzfft3d_(double const *in, void *out, int *nx, int *ny, int *nz,
    MPI_Fint fcomm, int *comm_rank, int *comm_size, int *opt);

/* =====================  Complex-to-Complex transform ========================= */

void compute_z2z_ffte(int const inbox_low[3], int const inbox_high[3],
                    int const outbox_low[3], int const outbox_high[3], MPI_Comm const comm,
                    void const *in, void *out, int *ffte_options, double *timer) {

    int comm_size, opt, nd[3];
    MPI_Fint fcomm = MPI_Comm_c2f(comm);

    /*
    Plan creation ...
    void *plan;
    */

    MPI_Barrier(comm);
    timer[0] = -MPI_Wtime();

    /*
    plan create
    */

    for (int i = 0; i < 3; ++i)
        nd[i] = inbox_high[i] - inbox_low[i] + 1;

    MPI_Comm_size(comm, &comm_size);
    opt = 0;
    pzfft3d_(in, out, nd+0, nd+1, nd+2, fcomm, &comm_size, &opt);

    MPI_Barrier(comm);
    timer[0] += MPI_Wtime();

    /*
    FFT execution ...
    */
    MPI_Barrier(comm);
    timer[1] = -MPI_Wtime();

    /*
    compute FFT
    */
    opt = 1;
    pzfft3d_(in, out, nd+0, nd+1, nd+2, fcomm, &comm_size, &opt);

    MPI_Barrier(comm);
    timer[1] = +MPI_Wtime();

    /*
    Delete plan ...
    */

    printf("Benchmarking FFTE: Get Complex-to-Complex backend using ptest3d.f\n");
}


/* =====================  Real-to-Complex transform ========================= */

void compute_d2z_ffte(int const inbox_low[3], int const inbox_high[3],
    int const outbox_low[3], int const outbox_high[3], MPI_Comm const comm,
    double const *in, void *out, int *ffte_options, double *timer) {

    int comm_rank, comm_size, opt, nd[3];
    MPI_Fint fcomm = MPI_Comm_c2f(comm);

    /*
    Plan creation ...
    void *plan;
    */

    MPI_Barrier(comm);
    timer[0] = -MPI_Wtime();

    /*
    plan create
    */

    for (int i = 0; i < 3; ++i)
        nd[i] = inbox_high[i] - inbox_low[i] + 1;

    MPI_Comm_rank(comm, &comm_rank);
    MPI_Comm_size(comm, &comm_size);
    opt = 0;
    pdzfft3d_(in, out, nd+0, nd+1, nd+2, fcomm, &comm_rank, &comm_size, &opt);

    MPI_Barrier(comm);
    timer[0] += MPI_Wtime();

    /*
    FFT execution ...
    */
    MPI_Barrier(comm);
    timer[1] = -MPI_Wtime();

    /*
    compute FFT
    */
    opt = 1;
    pdzfft3d_(in, out, nd+0, nd+1, nd+2, fcomm, &comm_rank, &comm_size, &opt);

    MPI_Barrier(comm);
    timer[1] = +MPI_Wtime();

    /*
    Delete plan ...
    */

    printf("Benchmarking FFTE: Get Real-to-Complex backend using prtest3d.f\n");
}

/* =====================  Complex-to-Real transform ========================= */

void compute_z2d_ffte(int const inbox_low[3], int const inbox_high[3],
                    int const outbox_low[3], int const outbox_high[3], MPI_Comm const comm,
                    void const *in, double *out, int *ffte_options, double *timer) 
{
// Missing!


}

#else 

void compute_z2z_ffte(int const inbox_low[3], int const inbox_high[3],
                    int const outbox_low[3], int const outbox_high[3], MPI_Comm const comm,
                    void const *in, void *out, int *ffte_options, double *timer) 
{}

void compute_d2z_ffte(int const inbox_low[3], int const inbox_high[3],
                    int const outbox_low[3], int const outbox_high[3], MPI_Comm const comm,
                    double const *in, void *out, int *ffte_options, double *timer) 
{}

void compute_z2d_ffte(int const inbox_low[3], int const inbox_high[3],
                    int const outbox_low[3], int const outbox_high[3], MPI_Comm const comm,
                    void const *in, double *out, int *ffte_options, double *timer) 
{}

#endif

#endif  /* ! FIBER_BACKEND_FFTE_H */
