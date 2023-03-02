#ifndef FIBER_FORTRAN_HEADER_INCLUDED
#define FIBER_FORTRAN_HEADER_INCLUDED

/* Mangling for Fortran global symbols without underscores. */
#define FIBER_FORTRAN_GLOBAL(name,NAME) name##_

/* Mangling for Fortran global symbols with underscores. */
#define FIBER_FORTRAN_GLOBAL_(name,NAME) name##_

/* Mangling for Fortran module symbols without underscores. */
#define FIBER_FORTRAN_MODULE(mod_name,name, mod_NAME,NAME) __##mod_name##_MOD_##name

/* Mangling for Fortran module symbols with underscores. */
#define FIBER_FORTRAN_MODULE_(mod_name,name, mod_NAME,NAME) __##mod_name##_MOD_##name

/*--------------------------------------------------------------------------*/
/* Mangle some symbols automatically.                                       */
#define fiber_pzfft3d FIBER_FORTRAN_GLOBAL(pzfft3d, PZFFT3D)
#define fiber_pdzfft3d FIBER_FORTRAN_GLOBAL(pdzfft3d, PDZFFT3D)

#endif
