!
! Fortran 77 solver of eigenvalue probem
!
      subroutine solve_eigenvalue_problem(
     +           nz, ny, nx, ng, nmix,
     +           flux, imap, sigt, sigp,
     +           nsigs, sigs,
     +           pitch, dz
     +           )
       
      implicit none
      integer nz, ny, nx, ng, nmix
      ! neutron flux array: flux(nz,ny,nx,ng)
      real*8 flux(:,:,:,:)
      ! map of material indexes: map(nz,ny,nx)
      integer imap(:,:,:)
      ! total cross section: sigt(nmix,ng)
      real*8 sigt(:,:)
      ! production cross section: sigp(nmix,ng)
      real*8 sigp(:,:)
      ! number of entries in scattering cross section matrix: sigs(nmix)
      integer nsigs(:)
      ! scattering cross section matrix: sigs(nmix,max(nsigs]))
      real*8 sigs(:,:)
      ! subassembly pitch
      real*8 pitch
      ! axial nodalization: dz(nz-2)
      real*8 dz(:)

      integer i

      end
