!module decomp_2d_ciface
!  use iso_c_binding
!  implicit none
!
!  contains

subroutine decomp_2d_c_init(nx, ny, nz, p_row, p_col) &
  bind(C, name="decomp_2d_init")
  use decomp_2d
  use iso_c_binding
  implicit none

  integer(C_INT),   value :: nx
  integer(C_INT),   value :: ny
  integer(C_INT),   value :: nz
  integer(C_INT),   value :: p_row
  integer(C_INT),   value :: p_col

  call decomp_2d_init(nx, ny, nz, p_row, p_col)
end subroutine decomp_2d_c_init

subroutine decomp_2d_c_fft_init(cpencil) &
  bind(C, name="decomp_2d_fft_init")
  use decomp_2d_fft
  use iso_c_binding
  implicit none

  integer(C_INT),   value :: cpencil

  integer                 :: fpencil

  ! This test is based on the values of the enum created in the header
  ! file
  if(cpencil.eq.1) then
    fpencil = PHYSICAL_IN_Z
  else
    fpencil = PHYSICAL_IN_X
  end if

  call decomp_2d_fft_init( fpencil )
end subroutine decomp_2d_c_fft_init

subroutine decomp_2d_c_get_local_sizes(cpencil, csize_0, csize_1, csize_2) &
  bind(C, name="decomp_2d_get_local_sizes")
  use decomp_2d_fft, only: PHYSICAL_IN_X, PHYSICAL_IN_Z
  use decomp_2d
  use iso_c_binding
  implicit none

  integer(C_INT),   value       :: cpencil
  integer(C_INT),   intent(out) :: csize_0
  integer(C_INT),   intent(out) :: csize_1
  integer(C_INT),   intent(out) :: csize_2

  integer                       :: fpencil
  ! This test is based on the values of the enum created in the header
  ! file
  if(cpencil.eq.1) then
    fpencil = PHYSICAL_IN_Z
  else
    fpencil = PHYSICAL_IN_X
  end if

  if(fpencil.eq.PHYSICAL_IN_X) then
    csize_0 = xsize(1)
    csize_1 = xsize(2)
    csize_2 = xsize(3)
  else
    csize_0 = zsize(1)
    csize_1 = zsize(2)
    csize_2 = zsize(3)
  end if

 !cxsize = xsize(1) ! For now, we try the first entry of this 3-entry array
 !cysize = ysize(2) ! For now, we try the first entry of this 3-entry array
 !czsize = zsize(3) ! For now, we try the first entry of this 3-entry array
end subroutine decomp_2d_c_get_local_sizes

!subroutine decomp_2d_c_fft_3d_c2c(cdata_in, cdata_out, direction) &
subroutine decomp_2d_c_fft_3d_c2c(nx_in, ny_in, nz_in, &
  cdata_in, nx_out, ny_out, nz_out, cdata_out, direction) &
  bind(C, name="decomp_2d_fft_3d_c2c")
  use decomp_2d
  use decomp_2d_fft
  use iso_c_binding
  use ISO_Fortran_env, only: stderr => ERROR_UNIT
  implicit none

  integer(C_INT),   value :: nx_in
  integer(C_INT),   value :: ny_in
  integer(C_INT),   value :: nz_in
  type(C_PTR),      value :: cdata_in
  integer(C_INT),   value :: nx_out
  integer(C_INT),   value :: ny_out
  integer(C_INT),   value :: nz_out
  type(C_PTR),      value :: cdata_out
  integer(C_INT),   value :: direction

  complex(mytype),  pointer :: fdata_in(:,:,:)
  complex(mytype),  pointer :: fdata_out(:,:,:)
  integer                   :: fdirection

  if(direction.eq.1) then
    fdirection = DECOMP_2D_FFT_FORWARD
  else
    fdirection = DECOMP_2D_FFT_BACKWARD
  end if

  if(.not. C_ASSOCIATED(cdata_in)) then
    write (stderr,*) "Error, cdata_in provided by the use is empty"
  endif
  call C_F_POINTER(cdata_in, fdata_in, &
    shape=(/ nx_in, ny_in, nz_in /)) ! Must provide size

  if(.not. C_ASSOCIATED(cdata_out)) then
    write (stderr,*) "Error, cdata_out provided by the use is empty"
  endif
  call C_F_POINTER(cdata_out, fdata_out, &
    shape=(/ nx_out, ny_out, nz_out /)) ! Must provide size

  call decomp_2d_fft_3d(fdata_in, fdata_out, fdirection)
 !call fft_3d_c2c(fdata_in, fdata_out, fdirection)

end subroutine decomp_2d_c_fft_3d_c2c

subroutine decomp_2d_c_fft_finalize() &
  bind(C, name="decomp_2d_fft_finalize")
  use decomp_2d_fft
  use iso_c_binding
  implicit none

  call decomp_2d_fft_finalize()
end subroutine decomp_2d_c_fft_finalize

subroutine decomp_2d_c_finalize() &
  bind(C, name="decomp_2d_finalize")
  use decomp_2d
  use decomp_2d
  use iso_c_binding
  implicit none

  call decomp_2d_finalize()
end subroutine decomp_2d_c_finalize

!end module decomp_2d_ciface
