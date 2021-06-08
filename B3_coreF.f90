!
! Fortran 95 solver of eigenvalue probem
!
subroutine solve_eigenvalue_problem(geom, nz, ny, nx, nt, ng, nmix, flux, imap, &
                                  & sigt, sigtra, sigp, nsigs, fsigs, tsigs, sigs, &
                                  & nsign2n, fsign2n, tsign2n, sign2n, &
                                  & chi, pitch, dz)

implicit none

! geometry flag ('squ' and 'hex')
character*3 geom
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
! index of mix in the neighbouring node
integer imix_n
! counter of inner (flux) iterations
integer niteri
! counter of outer (fission source) iterations
integer nitero
! absolute tolerance
real*8 atol
! node axial area-to-volume ratio
real*8 az_over_v
! node side area-to-volume ratio
real*8 aside_over_v
! sum of contributions from diffusion terms from neighbouring nodes
real*8 dif
! diffusion coefficient
real*8 D
! distance between nodes in z or xy direction
real*8 db
! working variable for flux
real*8 flx
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

! relative tolerance
rtol = 1.0e-6
! absolute tolerance
atol = 1.0e-6

! initialize fission source and total fission source
qf = 1.
tfs = sum(qf)

! change from python style to fortran style
imap = imap + 1

! side area to volume ratio of control volume 
if(geom == 'squ')then
   aside_over_v = 1./pitch
else ! geom == 'hex'
   aside_over_v = 2./(3.*pitch)
end if

! eigenvalue keff equal to ratio of total fission source at two iterations. 
! flux is normalise to total fission cource = 1 at previous iteration 
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
   do while(.not. converge_flux .and. niteri < 10)
      niteri = niteri + 1
      converge_flux = .true.
      do iz = 1, nz
      do iy = 1, ny
      do ix = 1, nx
         ! if (ix, iy, iz) is not a boundary condition node, i.e. not -1 (vac) and not -2 (ref)
         imix = imap(iz,iy,ix)
         if(imix > 0)then
            ! node axial area-to-volume ratio
            az_over_v = 1./dz(iz-1)
            do it = 1, nt
               do ig = 1, ng
                  mlt = 0.
                  dif = 0.
               
                  ! diffusion term: bottom
                  imix_n = imap(iz-1,iy,ix)
                  if(imix_n == 0)then
                     ! neighbour is vacuum
                     db = 0.5*dz(iz-1) + 0.71/sigtra(imix, ig)
                     D = 1./(3.*sigtra(imix, ig))
                     flx = 0.
                  else if(imix_n .eq. -1)then
                     ! neighbour is reflective
                     db = 1.
                     D = 0.
                     flx = 0.
                  else
                     ! neighbour is normal node (neither vacuum nor reflective)
                     db = 0.5*(dz(iz-2) + dz(iz-1))
                     D = 2.*db/(3.*sigtra(imix, ig)*dz(iz-1) + 3.*sigtra(imix_n, ig)*dz(iz-2))
                     flx = flux(iz-1,iy,ix,it,ig)
                  end if
                  mlt = mlt + D * az_over_v / db
                  dif = dif + D * flx * az_over_v / db
               
                  ! diffusion term: top
                  imix_n = imap(iz+1,iy,ix)
                  if(imix_n == 0)then
                     ! neighbour is vacuum
                     db = 0.5*dz(iz) + 0.71/sigtra(imix, ig)
                     D = 1./(3.*sigtra(imix, ig))
                     flx = 0.
                  else if(imix_n .eq. -1)then
                     ! neighbour is reflective
                     db = 1.
                     D = 0.
                     flx = 0.
                  else
                     ! neighbour is normal node (neither vacuum nor reflective)
                     db = 0.5*(dz(iz-1) + dz(iz))
                     D = 2.*db/(3.*sigtra(imix, ig)*dz(iz-1) + 3.*sigtra(imix_n, ig)*dz(iz))
                     flx = flux(iz+1,iy,ix,it,ig)
                  end if
                  mlt = mlt + D * az_over_v / db
                  dif = dif + D * flx * az_over_v / db
               
                  if(geom == 'squ' .or. geom == 'hex')then
                     ! diffusion term: west
                     call difxy(mlt,dif,iz,iy,ix-1,it,ig,imix,imap,nz,ny,nx,nt,ng,nmix,pitch,sigtra,flux,aside_over_v)
                     !imix_n = imap(iz,iy,ix-1)
                     !if(imix_n == 0)then
                     !   ! neighbour is vacuum
                     !   db = 0.5*pitch + 0.71/sigtra(imix, ig)
                     !   D = 1./(3.*sigtra(imix, ig))
                     !   flx = 0.
                     !else if(imix_n .eq. -1)then
                     !   ! neighbour is reflective
                     !   db = 1.
                     !   D = 0.
                     !   flx = 0.
                     !else
                     !   ! neighbour is normal node (neither vacuum nor reflective)
                     !   db = pitch
                     !   D = 2./(3.*sigtra(imix, ig) + 3.*sigtra(imix_n, ig))
                     !   flx = flux(iz,iy,ix-1,it,ig)
                     !end if
                     !mlt = mlt + D * aside_over_v / db
                     !dif = dif + D * flx * aside_over_v / db
                     
                     ! diffusion term: east
                     call difxy(mlt,dif,iz,iy,ix+1,it,ig,imix,imap,nz,ny,nx,nt,ng,nmix,pitch,sigtra,flux,aside_over_v)
                     !imix_n = imap(iz,iy,ix+1)
                     !if(imix_n == 0)then
                     !   ! neighbour is vacuum
                     !   db = 0.5*pitch + 0.71/sigtra(imix, ig)
                     !   D = 1./(3.*sigtra(imix, ig))
                     !   flx = 0.
                     !else if(imix_n .eq. -1)then
                     !   ! neighbour is reflective
                     !   db = 1.
                     !   D = 0.
                     !   flx = 0.
                     !else
                     !   ! neighbour is normal node (neither vacuum nor reflective)
                     !   db = pitch
                     !   D = 2./(3.*sigtra(imix, ig) + 3.*sigtra(imix_n, ig))
                     !   flx = flux(iz,iy,ix+1,it,ig)
                     !end if
                     !mlt = mlt + D * aside_over_v / db
                     !dif = dif + D * flx * aside_over_v / db
                  end if
               
                  if(geom == 'squ')then
                     ! diffusion term: north (square geometry)
                     call difxy(mlt,dif,iz,iy-1,ix,it,ig,imix,imap,nz,ny,nx,nt,ng,nmix,pitch,sigtra,flux,aside_over_v)
                     !imix_n = imap(iz,iy-1,ix)
                     !if(imix_n == 0)then
                     !   ! neighbour is vacuum
                     !   db = 0.5*pitch + 0.71/sigtra(imix, ig)
                     !   D = 1./(3.*sigtra(imix, ig))
                     !   flx = 0.
                     !else if(imix_n .eq. -1)then
                     !   ! neighbour is reflective
                     !   db = 1.
                     !   D = 0.
                     !   flx = 0.
                     !else
                     !   ! neighbour is normal node (neither vacuum nor reflective)
                     !   db = pitch
                     !   D = 2./(3.*sigtra(imix, ig) + 3.*sigtra(imix_n, ig))
                     !   flx = flux(iz,iy-1,ix,it,ig)
                     !end if
                     !mlt = mlt + D * aside_over_v / db
                     !dif = dif + D * flx * aside_over_v / db
               
                     ! diffusion term: south
                     call difxy(mlt,dif,iz,iy+1,ix,it,ig,imix,imap,nz,ny,nx,nt,ng,nmix,pitch,sigtra,flux,aside_over_v)
                     !imix_n = imap(iz,iy+1,ix)
                     !if(imix_n == 0)then
                     !   ! neighbour is vacuum
                     !   db = 0.5*pitch + 0.71/sigtra(imix, ig)
                     !   D = 1./(3.*sigtra(imix, ig))
                     !   flx = 0.
                     !else if(imix_n .eq. -1)then
                     !   ! neighbour is reflective
                     !   db = 1.
                     !   D = 0.
                     !   flx = 0.
                     !else
                     !   ! neighbour is normal node (neither vacuum nor reflective)
                     !   db = pitch
                     !   D = 2./(3.*sigtra(imix, ig) + 3.*sigtra(imix_n, ig))
                     !   flx = flux(iz,iy+1,ix,it,ig)
                     !end if
                     !mlt = mlt + D * aside_over_v / db
                     !dif = dif + D * flx * aside_over_v / db
               
                  else if(geom == 'hex')then
               
                     ! diffusion term: north-west (hexagonal geometry)
                     if(mod(iy,2) == 0)then ! even
                     call difxy(mlt,dif,iz,iy-1,ix,it,ig,imix,imap,nz,ny,nx,nt,ng,nmix,pitch,sigtra,flux,aside_over_v)
                        !imix_n = imap(iz,iy-1,ix)
                     else ! odd
                     call difxy(mlt,dif,iz,iy-1,ix-1,it,ig,imix,imap,nz,ny,nx,nt,ng,nmix,pitch,sigtra,flux,aside_over_v)
                        !imix_n = imap(iz,iy-1,ix-1)
                     end if
                     !if(imix_n == 0)then
                     !   ! neighbour is vacuum
                     !   db = 0.5*pitch + 0.71/sigtra(imix, ig)
                     !   D = 1./(3.*sigtra(imix, ig))
                     !   flx = 0.
                     !else if(imix_n .eq. -1)then
                     !   ! neighbour is reflective
                     !   db = 1.
                     !   D = 0.
                     !   flx = 0.
                     !else
                     !   ! neighbour is normal node (neither vacuum nor reflective)
                     !   db = pitch
                     !   D = 2./(3.*sigtra(imix, ig) + 3.*sigtra(imix_n, ig))
                     !   if(mod(iy,2) == 0)then ! even
                     !      flx = flux(iz,iy-1,ix,it,ig)
                     !   else ! odd
                     !      flx = flux(iz,iy-1,ix-1,it,ig)
                     !   end if
                     !end if
                     !mlt = mlt + D * aside_over_v / db
                     !dif = dif + D * flx * aside_over_v / db
               
                     ! diffusion term: north-east (hexagonal geometry)
                     if(mod(iy,2) == 0)then ! even
                     call difxy(mlt,dif,iz,iy-1,ix+1,it,ig,imix,imap,nz,ny,nx,nt,ng,nmix,pitch,sigtra,flux,aside_over_v)
                        !imix_n = imap(iz,iy-1,ix+1)
                     else ! odd
                     call difxy(mlt,dif,iz,iy-1,ix,it,ig,imix,imap,nz,ny,nx,nt,ng,nmix,pitch,sigtra,flux,aside_over_v)
                        !imix_n = imap(iz,iy-1,ix)
                     end if
                     !if(imix_n == 0)then
                     !   ! neighbour is vacuum
                     !   db = 0.5*pitch + 0.71/sigtra(imix, ig)
                     !   D = 1./(3.*sigtra(imix, ig))
                     !   flx = 0.
                     !else if(imix_n .eq. -1)then
                     !   ! neighbour is reflective
                     !   db = 1.
                     !   D = 0.
                     !   flx = 0.
                     !else
                     !   ! neighbour is normal node (neither vacuum nor reflective)
                     !   db = pitch
                     !   D = 2./(3.*sigtra(imix, ig) + 3.*sigtra(imix_n, ig))
                     !   if(mod(iy,2) == 0)then ! even
                     !      flx = flux(iz,iy-1,ix+1,it,ig)
                     !   else ! odd
                     !      flx = flux(iz,iy-1,ix,it,ig)
                     !   end if
                     !end if
                     !mlt = mlt + D * aside_over_v / db
                     !dif = dif + D * flx * aside_over_v / db
               
                     ! diffusion term: south-west (hexagonal geometry)
                     if(mod(iy,2) == 0)then ! even
                     call difxy(mlt,dif,iz,iy+1,ix,it,ig,imix,imap,nz,ny,nx,nt,ng,nmix,pitch,sigtra,flux,aside_over_v)
                        !imix_n = imap(iz,iy+1,ix)
                     else ! odd
                     call difxy(mlt,dif,iz,iy+1,ix-1,it,ig,imix,imap,nz,ny,nx,nt,ng,nmix,pitch,sigtra,flux,aside_over_v)
                        !imix_n = imap(iz,iy+1,ix-1)
                     end if
                     !if(imix_n == 0)then
                     !   ! neighbour is vacuum
                     !   db = 0.5*pitch + 0.71/sigtra(imix, ig)
                     !   D = 1./(3.*sigtra(imix, ig))
                     !   flx = 0.
                     !else if(imix_n .eq. -1)then
                     !   ! neighbour is reflective
                     !   db = 1.
                     !   D = 0.
                     !   flx = 0.
                     !else
                     !   ! neighbour is normal node (neither vacuum nor reflective)
                     !   db = pitch
                     !   D = 2./(3.*sigtra(imix, ig) +3.*sigtra(imix_n, ig))
                     !   if(mod(iy,2) == 0)then ! even
                     !      flx = flux(iz,iy+1,ix,it,ig)
                     !   else ! odd
                     !      flx = flux(iz,iy+1,ix-1,it,ig)
                     !   end if
                     !end if
                     !mlt = mlt + D * aside_over_v / db
                     !dif = dif + D * flx * aside_over_v / db
               
                     ! diffusion term: to south-east (hexagonal geometry)
                     if(mod(iy,2) == 0)then ! even
                     call difxy(mlt,dif,iz,iy+1,ix+1,it,ig,imix,imap,nz,ny,nx,nt,ng,nmix,pitch,sigtra,flux,aside_over_v)
                        !imix_n = imap(iz,iy+1,ix+1)
                     else ! odd
                     call difxy(mlt,dif,iz,iy+1,ix,it,ig,imix,imap,nz,ny,nx,nt,ng,nmix,pitch,sigtra,flux,aside_over_v)
                        !imix_n = imap(iz,iy+1,ix)
                     end if
                     !if(imix_n == 0)then
                     !   ! neighbour is vacuum
                     !   db = 0.5*pitch + 0.71/sigtra(imix, ig)
                     !   D = 1./(3.*sigtra(imix, ig))
                     !   flx = 0.
                     !else if(imix_n .eq. -1)then
                     !   ! neighbour is reflective
                     !   db = 1.
                     !   D = 0.
                     !   flx = 0.
                     !else
                     !   ! neighbour is normal node (neither vacuum nor reflective)
                     !   db = pitch
                     !   D = 2./(3.*sigtra(imix, ig) + 3.*sigtra(imix_n, ig))
                     !   if(mod(iy,2) == 0)then ! even
                     !      flx = flux(iz,iy+1,ix+1,it,ig)
                     !   else ! odd
                     !      flx = flux(iz,iy+1,ix,it,ig)
                     !   end if
                     !end if
                     !mlt = mlt + D * aside_over_v / db
                     !dif = dif + D * flx * aside_over_v / db
                        
                     ! removal xs
                     sigr = sigt(imix,ig)
                     ! scattering source
                     qs = 0.
                     do indx = 1, nsigs(imix)
                        f = fsigs(imix,indx)+1
                        t = tsigs(imix,indx)+1
                        if(f .ne. ig .and. t == ig)then
                           qs = qs + sigs(imix,indx) * flux(iz,iy,ix,it,f)
                        end if
                        if(f == ig .and. t == ig)then
                           sigr = sigr - sigs(imix,indx)
                        end if
                     end do
                     ! n2n source
                     qn2n = 0.
                     do indx = 1, nsign2n(imix)
                        f = fsign2n(imix,indx)+1
                        t = tsign2n(imix,indx)+1
                        if(f .ne. ig .and. t == ig)then
                           qn2n = qn2n + 2.*sign2n(imix,indx)*flux(iz,iy,ix,it,f)
                        end if
                     end do
                        
                     mlt = mlt + sigr
               
                     ! fission source
                     qfis = chi(imix,ig)*qf(iz,iy,ix,it)/keff
               
                     ! neutron flux
                     fluxnew = (dif + qs + qn2n + qfis)/mlt
                     if(converge_flux)then
                        converge_flux = abs(fluxnew - flux(iz,iy,ix,it,ig)) < rtol*abs(fluxnew) + atol
                     end if
                     flux(iz,iy,ix,it,ig) = fluxnew
               
                  else
                     write(*,*)'***ERROR: unknown core geometry in B3_coreF: ', geom
                     stop
                  end if
               end do
            end do
         end if
      end do
      end do
      end do
   end do

   ! calculate node-wise fission source qf and total fission source tfs
   tfs = 0.
   do iz = 1, nz
   do iy = 1, ny
   do ix = 1, nx
      ! if (ix, iy, iz) is not a boundary condition node, i.e. not -1 (vac) and not -2 (ref)
      imix = imap(iz,iy,ix)
      if(imix > 0)then
         do it = 1, nt
            qf(iz,iy,ix,it) = 0.
            do ig = 1, ng
               qf(iz,iy,ix,it) = qf(iz,iy,ix,it) + sigp(imix,ig)*flux(iz,iy,ix,it,ig)
            end do
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

   write(*,*)'k-effective: ', keff, 'nitero = ', nitero
   !if(isnan(keff)) stop

end do

! change from fortran style to python style
imap = imap - 1

end subroutine

!--------------------------------------------------------------------------------------------------
subroutine difxy(mlt, dif, jz, jy, jx, jt, ig, imix, imap, nz, ny, nx, nt, ng, nmix, pitch, sigtra, flux, aside_over_v)

implicit none

integer jz, jy, jx, jt, ig, imix, nz, ny, nx, nt, ng, nmix, imap(nz,ny,nx)
real*8 mlt, dif, pitch, sigtra(nmix,ng), flux(nz,ny,nx,nt,ng), aside_over_v

integer imix_n
real*8 db, D, flx

imix_n = imap(jz,jy,jx)
if(imix_n == 0)then
   ! neighbour is vacuum
   db = 0.5*pitch + 0.71/sigtra(imix, ig)
   D = 1./(3.*sigtra(imix, ig))
   flx = 0.
else if(imix_n .eq. -1)then
   ! neighbour is reflective
   db = 1.
   D = 0.
   flx = 0.
else
   ! neighbour is normal node (neither vacuum nor reflective)
   db = pitch
   D = 2./(3.*sigtra(imix, ig) + 3.*sigtra(imix_n, ig))
   flx = flux(jz,jy,jx,jt,ig)
end if
mlt = mlt + D * aside_over_v / db
dif = dif + D * flx * aside_over_v / db

end subroutine
