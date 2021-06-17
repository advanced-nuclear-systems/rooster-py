!
! Fortran 95 solver of eigenvalue probem
!
subroutine solve_eigenvalue_problem(geom, nz, ny, nx, nt, ng, nmix, flux, imap, &
                                  & sigt, sigtra, sigp, nsigs, fsigs, tsigs, sigs, &
                                  & nsign2n, fsign2n, tsign2n, sign2n, &
                                  & chi, pitch, dz)

use omp_lib
 
implicit none

! geometry flag ('squar', 'hex01', 'hex06', 'hex24')
character*5 geom
! number of nodes in z, y, x dimensions and number of triangles per hexagon, number of groups and number of mixes
integer nz, ny, nx, nt, ng, nmix
! neutron flux array: flux(nz,ny,nx,nt,ng)
real*8 flux(:,:,:,:,:)
! map of material indexes: imap(nz,ny,nx)
integer imap(:,:,:)
! total cross section: sigt(nmix,ng)
real*8 sigt(:,:)
! transport cross section: sigtra(nmix,ng)
real*8 sigtra(:,:)
! production cross section: sigp(nmix,ng)
real*8 sigp(:,:)
! number of entries in scattering cross section matrix: sigs(nmix)
integer nsigs(:)
! index of energy group from which scattering occurs: fsigs(nmix,max(nsigs))
integer fsigs(:,:)
! index of energy group to which scattering occurs: fsigs(nmix,max(nsigs))
integer tsigs(:,:)
! scattering cross section matrix: sigs(nmix,max(nsigs)))
real*8 sigs(:,:)
! number of entries in n2n cross section matrix: sign2n(nmix)
integer nsign2n(:)
! index of energy group from which n2n occurs: fsign2n(nmix,max(nsign2n))
integer fsign2n(:,:)
! index of energy group to which n2n occurs: fsigs(nmix,max(nsign2n))
integer tsign2n(:,:)
! n2n cross section matrix: sign2n(nmix,max(nsign2n)))
real*8 sign2n(:,:)
! fission spectrum: chi(nmix,ng)
real*8 chi(:,:)
! subassembly pitch
real*8 pitch
! axial nodalization: dz(nz-2)
real*8 dz(:)

! convergence flags
logical converge_k, converge_flux
! from and to indices for scattering matrix
integer f, t
! for do loop
integer ix, iy, iz, it, ig, indx
! index of mix in the node
integer imix
! counter of inner (flux) iterations
integer niteri
! counter of outer (fission source) iterations
integer nitero
! number of threads in case OMP is used
integer :: nthreads = 0
! absolute tolerance
real*8 atol
! node axial area-to-volume ratio
real*8 az_over_v
! node side area-to-volume ratio
real*8 aside_over_v
! distance between nodes in xy plane
real*8 db
! sum of contributions from diffusion terms from neighbouring nodes
real*8 dif
! flux for evaluating the difference from the previous iteration
real*8 fluxnew
! multiplication factor: keff
real*8 keff
! multiplication factor for evaluating the difference from the previous iteration
real*8 knew
! multiplier at flux
real*8 mlt
! fission source
real*8 qf(nz,ny,nx,nt)
! total fission source
real*8 tfs
! fission source
real*8 qfis
! n2n source
real*8 qn2n
! scattering source
real*8 qs
! relative tolerance
real*8 rtol
! removal cross section
real*8 sigr

!$ call OMP_set_dynamic(.true.)

! relative tolerance
rtol = 1.0e-8
! absolute tolerance
atol = 1.0e-8

! initialize fission source
qf = 1.

! change from python style to fortran style
imap = imap + 1

! side area to volume ratio of control volume 
if(geom == 'squar')then
   aside_over_v = 1./pitch
else if(geom == 'hex01')then
   aside_over_v = 2./(3.*pitch)
else if(geom == 'hex06')then
   aside_over_v = 4./pitch
else ! geom == 'hex24'
   aside_over_v = 8./pitch
end if

! eigenvalue keff equal to ratio of total fission source at two iterations. 
! flux is normalise to total fission source = 1 at previous iteration 
keff = 1.

! convergence flag
converge_k = .false.
! outer iteration counter
nitero = 0
do while(.not. converge_k .and. nitero < 1000)
   nitero = nitero + 1
   ! initialize flux convergence flag
   converge_flux = .false.
   ! inner iteration counter
   niteri = 0
   do while(.not. converge_flux .and. niteri < 5)
      niteri = niteri + 1
      converge_flux = .true.
      !$omp parallel do default(shared) &
      !$omp private(imix,az_over_v,mlt,dif,db,sigr,qs,f,t,qfis,fluxnew)
      do iz = 1,nz
         do iy = 1,ny
            do ix = 1,nx
               !$ if(omp_get_thread_num() == 0 .and. iz == 1 .and. iy == 1 .and. ix == 1) nthreads = OMP_get_num_threads()
            
               ! if (ix, iy, iz) is not a boundary condition node, i.e. not -1 (vac) and not -2 (ref)
               imix = imap(iz,iy,ix)
               if(imix > 0)then
                  ! node axial area-to-volume ratio
                  az_over_v = 1./dz(iz-1)
                  do it = 1,nt
                     do ig = 1,ng
                        mlt = 0.
                        dif = 0.
            
                        ! diffusion terms (mlt and dif) in different nodalizations and different directions
                        call mltdif(mlt, dif, iz, iy, ix, it, ig, imix, imap, nz, ny, nx, nt, ng, nmix, &
                                    pitch, dz, sigtra, flux, aside_over_v, az_over_v, geom)
            
                        ! removal xs
                        sigr = sigt(imix,ig)
                        do indx = 1,nsigs(imix)
                           f = fsigs(imix,indx)+1
                           t = tsigs(imix,indx)+1
                           if(f == ig .and. t == ig)then
                              sigr = sigr - sigs(imix,indx)
                           end if
                        end do
                           
                        mlt = mlt + sigr
            
                        ! scattering source
                        qs = 0.
                        do indx = 1,nsigs(imix)
                           f = fsigs(imix,indx)+1
                           t = tsigs(imix,indx)+1
                           if(f .ne. ig .and. t == ig)then
                              qs = qs + sigs(imix,indx) * flux(iz,iy,ix,it,f)
                           end if
                        end do
                        ! n2n source
                        qn2n = 0.
                        do indx = 1,nsign2n(imix)
                           f = fsign2n(imix,indx)+1
                           t = tsign2n(imix,indx)+1
                           if(f .ne. ig .and. t == ig)then
                              qn2n = qn2n + 2.*sign2n(imix,indx)*flux(iz,iy,ix,it,f)
                           end if
                        end do
                     
                        ! fission source
                        qfis = chi(imix,ig)*qf(iz,iy,ix,it)/keff
                     
                        ! neutron flux
                        fluxnew = (dif + qs + qn2n + qfis)/mlt
                        if(converge_flux)then
                           converge_flux = abs(fluxnew - flux(iz,iy,ix,it,ig)) < rtol*abs(fluxnew) + atol
                        end if
                        flux(iz,iy,ix,it,ig) = fluxnew
                     end do
                  end do
               end if
            end do
         end do
      end do
      !$omp end parallel do
   end do

   ! calculate node-wise fission source qf and total fission source tfs
   !$omp parallel do default(shared) private(imix)
   do iz = 1,nz
      do iy = 1,ny
         do ix = 1,nx
            ! if (ix, iy, iz) is not a boundary condition node, i.e. not -1 (vac) and not -2 (ref)
            imix = imap(iz,iy,ix)
            if(imix > 0)then
               do it = 1,nt
                  qf(iz,iy,ix,it) = 0.
                  do ig = 1,ng
                     qf(iz,iy,ix,it) = qf(iz,iy,ix,it) + sigp(imix,ig)*flux(iz,iy,ix,it,ig)
                  end do
               end do
            end if
         end do
      end do
   end do
   !$omp end parallel do

   ! calculate total fission source tfs
   tfs = 0.
   do iz = 1,nz
      do iy = 1,ny
         do ix = 1,nx
            ! if (ix, iy, iz) is not a boundary condition node, i.e. not -1 (vac) and not -2 (ref)
            imix = imap(iz,iy,ix)
            if(imix > 0)then
               do it = 1,nt
                  tfs = tfs + qf(iz,iy,ix,it)
               end do
            end if
         end do
      end do
   end do

   ! new k-effective is the ratio of total fission sources at the current (tfs) and previous (1.0) iterations
   knew = tfs
   converge_k = abs(knew - keff) < rtol*abs(knew) + atol
   keff = knew

   if(nthreads == 0)then
      write(*,'("keff: ",f13.6, " | niteri: ", i3, " | nitero: ", i3, " | ")') &
            keff,niteri,nitero
   else
      write(*,'("keff: ",f13.6, " | niteri: ", i3, " | nitero: ", i3, " | OMPthreads: ", i3, " | ")') &
            keff,niteri,nitero,nthreads
   end if
   if(isnan(keff)) stop

end do
  
! change from fortran style to python style
imap = imap - 1

end subroutine

!--------------------------------------------------------------------------------------------------
! Diffusion terms (mlt and dif) in different nodalizations and different directions
!
subroutine mltdif(mlt, dif, iz, iy, ix, it, ig, imix, imap, nz, ny, nx, nt, ng, nmix, &
                  pitch, dz, sigtra, flux, aside_over_v, az_over_v, geom)

implicit none
integer iz, iy, ix, it, ig, imix, nz, ny, nx, nt, ng, nmix, imap(nz,ny,nx)
real*8 mlt, dif, pitch, dz(nz-2), sigtra(nmix,ng), flux(nz,ny,nx,nt,ng), aside_over_v, az_over_v
character*5 geom

real*8 db

! diffusion terms in z direction for all geometries: mlt and dif
call difz(mlt,dif,iz,iz-1,iy,ix,it,ig,imix,imap,nz,ny,nx,nt,ng,nmix,dz,sigtra,flux,az_over_v)
call difz(mlt,dif,iz,iz+1,iy,ix,it,ig,imix,imap,nz,ny,nx,nt,ng,nmix,dz,sigtra,flux,az_over_v)

! square geometry
if(geom == 'squar')then
   ! diffusion terms in xy direction: mlt and dif
   call difxy(mlt,dif,iz,iy,ix-1,it,ig,imix,imap,nz,ny,nx,nt,ng,nmix,pitch,sigtra,flux,aside_over_v)
   call difxy(mlt,dif,iz,iy,ix+1,it,ig,imix,imap,nz,ny,nx,nt,ng,nmix,pitch,sigtra,flux,aside_over_v)
   call difxy(mlt,dif,iz,iy-1,ix,it,ig,imix,imap,nz,ny,nx,nt,ng,nmix,pitch,sigtra,flux,aside_over_v)
   call difxy(mlt,dif,iz,iy+1,ix,it,ig,imix,imap,nz,ny,nx,nt,ng,nmix,pitch,sigtra,flux,aside_over_v)

! hexagonal geometry with 1 node per hexagon
else if(geom == 'hex01')then
   ! diffusion terms in xy direction: mlt and dif
   if(mod(iy,2) == 0)then ! even
      call difxy(mlt,dif,iz,iy-1,ix,it,ig,imix,imap,nz,ny,nx,nt,ng,nmix,pitch,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,iy-1,ix+1,it,ig,imix,imap,nz,ny,nx,nt,ng,nmix,pitch,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,iy,ix-1,it,ig,imix,imap,nz,ny,nx,nt,ng,nmix,pitch,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,iy,ix+1,it,ig,imix,imap,nz,ny,nx,nt,ng,nmix,pitch,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,iy+1,ix,it,ig,imix,imap,nz,ny,nx,nt,ng,nmix,pitch,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,iy+1,ix+1,it,ig,imix,imap,nz,ny,nx,nt,ng,nmix,pitch,sigtra,flux,aside_over_v)
   else ! odd
      call difxy(mlt,dif,iz,iy-1,ix-1,it,ig,imix,imap,nz,ny,nx,nt,ng,nmix,pitch,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,iy-1,ix,it,ig,imix,imap,nz,ny,nx,nt,ng,nmix,pitch,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,iy,ix-1,it,ig,imix,imap,nz,ny,nx,nt,ng,nmix,pitch,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,iy,ix+1,it,ig,imix,imap,nz,ny,nx,nt,ng,nmix,pitch,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,iy+1,ix-1,it,ig,imix,imap,nz,ny,nx,nt,ng,nmix,pitch,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,iy+1,ix,it,ig,imix,imap,nz,ny,nx,nt,ng,nmix,pitch,sigtra,flux,aside_over_v)
   end if

! hexagonal geometry with 6 nodes per hexagon
else if(geom == 'hex06')then
   db = pitch/3.
   ! diffusion terms in xy direction: mlt and dif
   if(it == 1)then
      ! north-east
      call difxy(mlt,dif,iz,iy,ix,6,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      if(mod(iy,2) == 0)then ! even
         call difxy(mlt,dif,iz,iy-1,ix+1,4,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      else ! odd
         call difxy(mlt,dif,iz,iy-1,ix,4,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      end if
      call difxy(mlt,dif,iz,iy,ix,2,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
   else if(it == 2)then
      ! east
      call difxy(mlt,dif,iz,iy,ix,1,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,iy,ix+1,5,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,iy,ix,3,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
   else if(it == 3)then
      ! south-east
      call difxy(mlt,dif,iz,iy,ix,2,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      if(mod(iy,2) == 0)then ! even
         call difxy(mlt,dif,iz,iy+1,ix+1,6,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      else ! odd
         call difxy(mlt,dif,iz,iy+1,ix,6,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      end if
      call difxy(mlt,dif,iz,iy,ix,4,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
   else if(it == 4)then
      ! south-west
      call difxy(mlt,dif,iz,iy,ix,3,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      if(mod(iy,2) == 0)then ! even
         call difxy(mlt,dif,iz,iy+1,ix,1,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      else ! odd
         call difxy(mlt,dif,iz,iy+1,ix-1,1,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      end if
      call difxy(mlt,dif,iz,iy,ix,5,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
   else if(it == 5)then
      ! west
      call difxy(mlt,dif,iz,iy,ix,4,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,iy,ix-1,2,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,iy,ix,6,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
   else if(it == 6)then
      ! north-west
      call difxy(mlt,dif,iz,iy,ix,5,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      if(mod(iy,2) == 0)then ! even
         call difxy(mlt,dif,iz,iy-1,ix,3,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      else ! odd
         call difxy(mlt,dif,iz,iy-1,ix-1,3,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      end if
      call difxy(mlt,dif,iz,iy,ix,1,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
   end if

! hexagonal geometry with 24 nodes per hexagon
else if(geom == 'hex24')then
   db = pitch/6.
   ! diffusion terms in xy direction: mlt and dif
   if(it == 1)then
      if(mod(iy,2) == 0)then ! even
         call difxy(mlt,dif,iz,iy-1,ix,19,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      else ! odd
         call difxy(mlt,dif,iz,iy-1,ix-1,19,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      end if
      call difxy(mlt,dif,iz,iy,ix,2,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,iy,ix,7,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
   else if(it == 2)then
      call difxy(mlt,dif,iz,iy,ix,1,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      if(mod(iy,2) == 0)then ! even
         call difxy(mlt,dif,iz,iy-1,ix+1,21,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      else ! odd
         call difxy(mlt,dif,iz,iy-1,ix,21,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      end if
      call difxy(mlt,dif,iz,iy,ix,3,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
   else if(it == 3)then
      call difxy(mlt,dif,iz,iy,ix,2,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,iy,ix,4,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,iy,ix,9,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
   else if(it == 4)then
      call difxy(mlt,dif,iz,iy,ix,3,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      if(mod(iy,2) == 0)then ! even
         call difxy(mlt,dif,iz,iy-1,ix+1,23,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      else ! odd
         call difxy(mlt,dif,iz,iy-1,ix,23,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      end if
      call difxy(mlt,dif,iz,iy,ix,5,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
   else if(it == 5)then
      call difxy(mlt,dif,iz,iy,ix,4,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,iy,ix+1,13,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,iy,ix,11,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
   else if(it == 6)then
      if(mod(iy,2) == 0)then ! even
         call difxy(mlt,dif,iz,iy-1,ix,24,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      else ! odd
         call difxy(mlt,dif,iz,iy-1,ix-1,24,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      end if
      call difxy(mlt,dif,iz,iy,ix,7,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,iy,ix,13,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
   else if(it == 7)then
      call difxy(mlt,dif,iz,iy,ix,6,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,iy,ix,1,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,iy,ix,8,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
   else if(it == 8)then
      call difxy(mlt,dif,iz,iy,ix,7,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,iy,ix,9,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,iy,ix,15,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
   else if(it == 9)then
      call difxy(mlt,dif,iz,iy,ix,8,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,iy,ix,3,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,iy,ix,10,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
   else if(it == 10)then
      call difxy(mlt,dif,iz,iy,ix,9,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,iy,ix,11,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,iy,ix,17,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
   else if(it == 11)then
      call difxy(mlt,dif,iz,iy,ix,10,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,iy,ix,5,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,iy,ix,12,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
   else if(it == 12)then
      call difxy(mlt,dif,iz,iy,ix,11,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,iy,ix+1,20,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,iy,ix,19,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
   else if(it == 13)then
      call difxy(mlt,dif,iz,iy,ix-1,5,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,iy,ix,6,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,iy,ix,14,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
   else if(it == 14)then
      call difxy(mlt,dif,iz,iy,ix,13,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,iy,ix,15,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,iy,ix,20,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
   else if(it == 15)then
      call difxy(mlt,dif,iz,iy,ix,14,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,iy,ix,8,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,iy,ix,16,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
   else if(it == 16)then
      call difxy(mlt,dif,iz,iy,ix,15,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,iy,ix,17,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,iy,ix,22,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
   else if(it == 17)then
      call difxy(mlt,dif,iz,iy,ix,16,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,iy,ix,10,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,iy,ix,18,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
   else if(it == 18)then
      call difxy(mlt,dif,iz,iy,ix,17,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,iy,ix,19,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,iy,ix,24,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
   else if(it == 19)then
      call difxy(mlt,dif,iz,iy,ix,18,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,iy,ix,12,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      if(mod(iy,2) == 0)then ! even
         call difxy(mlt,dif,iz,iy+1,ix+1,1,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      else ! odd
         call difxy(mlt,dif,iz,iy+1,ix,1,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      end if
   else if(it == 20)then
      call difxy(mlt,dif,iz,iy,ix-1,12,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,iy,ix,14,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,iy,ix,21,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
   else if(it == 21)then
      call difxy(mlt,dif,iz,iy,ix,20,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,iy,ix,22,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      if(mod(iy,2) == 0)then ! even
         call difxy(mlt,dif,iz,iy+1,ix,2,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      else ! odd
         call difxy(mlt,dif,iz,iy+1,ix-1,2,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      end if
   else if(it == 22)then
      call difxy(mlt,dif,iz,iy,ix,21,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,iy,ix,16,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,iy,ix,23,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
   else if(it == 23)then
      call difxy(mlt,dif,iz,iy,ix,22,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,iy,ix,24,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      if(mod(iy,2) == 0)then ! even
         call difxy(mlt,dif,iz,iy+1,ix,4,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      else ! odd
         call difxy(mlt,dif,iz,iy+1,ix-1,4,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      end if
   else if(it == 24)then
      call difxy(mlt,dif,iz,iy,ix,23,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,iy,ix,18,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      if(mod(iy,2) == 0)then ! even
         call difxy(mlt,dif,iz,iy+1,ix+1,6,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      else ! odd
         call difxy(mlt,dif,iz,iy+1,ix,6,ig,imix,imap,nz,ny,nx,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      end if
   end if
else
   write(*,*)'***ERROR: unknown core geometry in B3_coreF: ', geom
   stop
end if
                  
end subroutine

!--------------------------------------------------------------------------------------------------
! Diffusion terms (mlt and dif) in xy directions
!
subroutine difxy(mlt, dif, jz, jy, jx, jt, ig, imix, imap, nz, ny, nx, nt, ng, nmix, p, sigtra, flux, a_over_v)

implicit none

integer jz, jy, jx, jt, ig, imix, nz, ny, nx, nt, ng, nmix, imap(nz,ny,nx)
real*8 mlt, dif, p, sigtra(nmix,ng), flux(nz,ny,nx,nt,ng), a_over_v

integer imix_n
real*8 db, D

imix_n = imap(jz,jy,jx)
if(imix_n == 0)then
   ! neighbour is vacuum
   db = 0.5*p + 0.71/sigtra(imix, ig)
   D = 1./(3.*sigtra(imix, ig))
   mlt = mlt + D * a_over_v / db
else if(imix_n .ne. -1)then
   ! neighbour is normal node (neither vacuum nor reflective)
   D = 2./(3.*sigtra(imix, ig) + 3.*sigtra(imix_n, ig))
   mlt = mlt + D * a_over_v / p
   dif = dif + D * flux(jz,jy,jx,jt,ig) * a_over_v / p
end if

end subroutine

!--------------------------------------------------------------------------------------------------
! Diffusion terms (mlt and dif) in z directions
!
subroutine difz(mlt, dif, iz, jz, jy, jx, jt, ig, imix, imap, nz, ny, nx, nt, ng, nmix, dz, sigtra, flux, a_over_v)

implicit none

integer iz, jz, jy, jx, jt, ig, imix, nz, ny, nx, nt, ng, nmix, imap(nz,ny,nx)
real*8 mlt, dif, dz(nz-2), sigtra(nmix,ng), flux(nz,ny,nx,nt,ng), a_over_v

! index of mix in the neighbouring node
integer imix_n
real*8 db, D

imix_n = imap(jz,jy,jx)
if(imix_n == 0)then
   ! neighbour is vacuum
   db = 0.5*dz(iz-1) + 0.71/sigtra(imix, ig)
   D = 1./(3.*sigtra(imix, ig))
   mlt = mlt + D * a_over_v / db
else if(imix_n .ne. -1)then
   ! neighbour is normal node (neither vacuum nor reflective)
   db = 0.5*(dz(iz-1) + dz(jz-1))
   D = 2.*db/(3.*sigtra(imix, ig)*dz(iz-1) + 3.*sigtra(imix_n, ig)*dz(jz-1))
   mlt = mlt + D * a_over_v / db
   dif = dif + D * flux(jz,jy,jx,jt,ig) * a_over_v / db
end if

end subroutine
