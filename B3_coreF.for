      subroutine solve_eigenvalue_problem(
     +           nz, ny, nx, ng, nmix,
     +           flux, imap, sigt, sigp,
     +           pitch, dz
     +           )
       
      implicit none
      integer nz, ny, nx, ng, nmix, nsigs
      ! neutron flux array: flux(nz,ny,nx,ng)
      real*8 flux(:,:,:,:)
      ! map of material indexes: map(nz,ny,nx)
      integer imap(:,:,:)
      ! total cross section: sigt(nmix,ng)
      real*8 sigt(:,:)
      ! production cross section: sigp(nmix,ng)
      real*8 sigp(:,:)
      ! subassembly pitch
      real*8 pitch
      ! axial nodalization: dz(nz-2)
      real*8 dz(:)

      integer i

      end
