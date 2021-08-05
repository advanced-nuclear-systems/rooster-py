!
! Fortran 95 solver of eigenvalue probem
!
subroutine solve_eigenvalue_problem(meth, geom, nz, ny, nx, nt, ng, nmix, flux, imap, &
                                  & sigt, sigtra, sigp, nsigs, fsigs, tsigs, sigs, &
                                  & nsign2n, fsign2n, tsign2n, sign2n, &
                                  & chi, sigf, pitch, dz)

use omp_lib
 
implicit none

! method flag ('MC', 'df')
character*2 meth
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
! fission cross section: sigf(nmix,ng)
real*8 sigf(:,:)
! number of entries in scattering cross section matrix: sigs(nmix)
integer nsigs(:)
! index of energy group from which scattering occurs: fsigs(nmix,max(nsigs))
integer fsigs(:,:)
! index of energy group to which scattering occurs: tsigs(nmix,max(nsigs))
integer tsigs(:,:)
! scattering cross section matrix: sigs(nmix,max(nsigs)))
real*8 sigs(:,:)
! number of entries in n2n cross section matrix: sign2n(nmix)
integer nsign2n(:)
! index of energy group from which n2n occurs: fsign2n(nmix,max(nsign2n))
integer fsign2n(:,:)
! index of energy group to which n2n occurs: tsign2n(nmix,max(nsign2n))
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
integer ix, iy, iz, it, ig, indx, i, j
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

! Monte Carlo specific variables

! number of neutrons at the end of cycle and neutrons born at the beginning of cycle, number of inactive source cycles to skip before starting k-eff accumulation and number of active cycles for k-eff accumulation (to be done input parameters)
integer num_neutrons, num_neutrons_born, num_cycles_inactive, num_cycles_active
parameter(num_neutrons_born = 50000, num_cycles_inactive = 100, num_cycles_active = 100)

! cycle index
integer icycle
! group number of neutron
integer igroup(num_neutrons_born*2)
! neutron free path
real*8 free_path
! k-effective in a cycle
real*8 :: keff_active_cycle(num_cycles_active) = 1.0
! averaged k-effective
real*8 :: keff_expected(num_cycles_active) = 1.0
! pi number
real*8 :: pi = 3.1415926535897932d0
! k-effective standard deviation
real*8 :: sigma_keff(num_cycles_active) = 0.0d0
! total n2n cross section
real*8 sign2n_sum
! scattering cross section
real*8 sigs1(ng)
! total scattering cross section
real*8 sigs_sum
! majorant
real*8 sigtmax(ng)
! virtual cross section
real*8 sigv
! neutron terminate probability
real*8 terminate_p
! angles and directions]
real*8 teta, phi, dir_x, dir_y, dir_z
! weight of neutrons
real*8 :: weight(num_neutrons_born*2) = 1.0d0
real*8 :: weight0(num_neutrons_born*2) = 1.0d0
! coordinates of hex centres
real*8 xnode(nx,ny), ynode(ny), zb(nz-1)
! coordinates of neutrons
real*8 x(num_neutrons_born*2), y(num_neutrons_born*2), z(num_neutrons_born*2)

logical :: virtual_collision = .false.
logical absorbed

! temporal variables
integer itmp(num_neutrons_born*2), num_new, N, iactive, ixp
real*8 r, tmp, xtmp(num_neutrons_born*2), ytmp(num_neutrons_born*2), &
       ztmp(num_neutrons_born*2), wtmp(num_neutrons_born*2),keff_cycle, siga

!$ call OMP_set_dynamic(.true.)

! change from python style to fortran style
imap = imap + 1

if(meth == 'MC')then

   ! MONTE CARLO SOLVER

   ! prepare flux array
   do iz = 1,nz
      do iy = 1,ny
         do ix = 1,nx
            do ig = 1,ng
               flux(iz,iy,ix,1,ig) = 0.0d0
            end do
         end do
      end do
   end do

   ! find coordinate of the hexagon centres
   do iy = 1,ny
      ynode(iy) = dble(iy-1)*dsqrt(3.0d0)*pitch/2.0d0
      do ix = 1,nx
         if(mod(iy,2) == 0)then ! even row
            xnode(ix,iy) = dble(ix-1)*pitch + pitch/2.0d0
         else ! odd row 
            xnode(ix,iy) = dble(ix-1)*pitch
         end if
      end do
   end do
   zb(1) = 0.0d0
   do iz = 2,nz-1
      zb(iz) = zb(iz-1) + dz(iz-1)
   end do
   ! define the majorant: the maximum total cross section vector 
   do ig = 1,ng
      sigtmax(ig) = 0.0d0
      do imix = 1,nmix
         sigtmax(ig) = dmax1(sigtmax(ig), sigt(imix,ig))
      end do
   end do

   ! neutrons are assumed born randomly distributed in the core with weight 1 with sampled fission energy spectrum
   num_neutrons = num_neutrons_born
   do i = 1,num_neutrons_born
      do
         ! generate random y coordinate
         call random_number(r)
         y(i) = (ynode(ny)+pitch)*r - 0.5d0*pitch
         iy = minloc(dabs(ynode - y(i)),1)
         
         ! generate random x coordinate
         call random_number(r)
         x(i) = (xnode(nx,iy)+pitch)*r - 0.5d0*pitch
         ix = minloc(dabs(xnode(:,iy) - x(i)),1)
         
         ! generate random z coordinate
         call random_number(r)
         z(i) = zb(nz-1)*r
         iz = minloc(dabs(zb - z(i)),1)
         if(z(i) > zb(iz)) iz = iz + 1

         imix = imap(iz,iy,ix)
         if(imix > 0 .and. sum(chi(imix,:)) > 0.0d0)then
            weight(i) = 1.0d0
            ! sample group by comparing cumulative sum of fission spectrum with random number
            call random_number(r)
            tmp = 0.0d0
            ig = ng + 1
            do while(tmp <= r)
               ig = ig - 1
               tmp = tmp + chi(imix,ig)
            end do
            igroup(i) = ig
            exit
         end if
      end do
   end do

   ! main batch cycle
   do icycle = 1,(num_cycles_inactive + num_cycles_active)

      ! normalize the weights of the neutrons to make the total weight equal to num_neutrons_born (equivalent to division by keff_cycle)
      tmp = 0.0d0
      do i = 1,num_neutrons
         tmp = tmp + weight(i)
      end do
      do i = 1,num_neutrons
         weight(i) = weight(i) / tmp * num_neutrons_born
         weight0(i) = weight(i)
      end do

      !$omp parallel do default(shared) &
      !$omp private(absorbed,r,free_path,virtual_collision,teta,phi, &
      !$omp dir_x,dir_y,dir_z,iy,ix,ixp,iz,imix,sigv,f,t,sigs1,sigs_sum, &
      !$omp tmp,ig,sign2n_sum,siga)
      ! loop over neutrons
      do i = 1,num_neutrons
          absorbed = .false.
          ! neutron random walk cycle: from emission to absorption
          do while(.not. absorbed)
             ! sample free path length according to the Woodcock method
             call random_number(r)
             free_path = -dlog(r)/sigtmax(igroup(i))
             if(.not. virtual_collision)then
                ! sample the direction of neutron flight assuming both fission and scattering are isotropic in the lab
                call random_number(r)
                teta = pi*r
                call random_number(r)
                phi = 2.0d0*pi*r
                dir_x = dsin(teta)*dcos(phi)
                dir_y = dsin(teta)*dsin(phi)
                dir_z = dcos(teta)
             end if
             ! fly
             x(i) = x(i) + free_path * dir_x
             y(i) = y(i) + free_path * dir_y
             z(i) = z(i) + free_path * dir_z

             ! find neutron position in structured mesh (ix, iy, iz)
             iy = minloc(dabs(ynode - y(i)),1)
             ix = minloc(dabs(xnode(:,iy) - x(i)),1)
             if(ynode(iy) - y(i) > 0.5d0*pitch/dsqrt(3.0d0))then
                ixp = minloc(dabs(xnode(:,iy-1) - x(i)),1)
                if(dabs(xnode(ixp,iy-1) - x(i)) < dabs(xnode(ix,iy) - x(i)))then
                   iy = iy - 1
                   ix = ixp
                end if
             else if(y(i) - ynode(iy) > 0.5d0*pitch/dsqrt(3.0d0))then
                ixp = minloc(dabs(xnode(:,iy+1) - x(i)),1)
                if(dabs(xnode(ixp,iy+1) - x(i)) < dabs(xnode(ix,iy) - x(i)))then
                   iy = iy + 1
                   ix = ixp
                end if
             end if
             iz = minloc(dabs(zb - z(i)),1)
             if(z(i) > zb(iz)) iz = iz + 1
             ! mixture index
             imix = imap(iz,iy,ix)

             if(imix == 0)then
                ! neutron is out of core
                weight(i) = 0.0d0
                absorbed = .true.
             else
                ! virtual cross section
                sigv = sigtmax(igroup(i)) - sigt(imix,igroup(i))
                ! sample the type of the collision: virtual (do nothing) or real
                call random_number(r)
                if(sigv/sigtmax(igroup(i)) >= r)then 
                   ! virtual collision
                   virtual_collision = .true.
                else
                   ! real collision
                   virtual_collision = .false.
                   ! score reactions with account for weight divided by the total cross section to estimate flux
                   if(icycle > num_cycles_inactive)then
                      flux(iz,iy,ix,1,igroup(i)) = flux(iz,iy,ix,1,igroup(i)) + weight(i)
                   end if
                   ! scattering xs for group igroup(i)
                   do ig = 1,ng
                      sigs1(ig) = 0.0d0
                   end do
                   do indx = 1,nsigs(imix)
                      f = fsigs(imix,indx)+1
                      t = tsigs(imix,indx)+1
                      if(f == igroup(i))sigs1(t) = sigs1(t) + sigs(imix,indx)
                   end do
                   ! total scattering xs for group igroup(i)
                   sigs_sum = 0.0d0
                   do ig = 1,ng
                      sigs_sum = sigs_sum + sigs1(ig)
                   end do
                   ! sample type of the collision: scattering or absorption]
                   call random_number(r)
                   if(sigs_sum/sigt(imix,igroup(i)) >= r)then ! isotropic scattering
                      ! sample group for the secondary neutron by comparing cumulative sum of partial scattering xs with random number
                      call random_number(r)
                      tmp = 0.0d0
                      ig = ng + 1
                      do while(tmp <= r)
                         ig = ig - 1
                         tmp = tmp + sigs1(ig)/sigs_sum
                      end do
                      igroup(i) = ig
                   else ! absorption
                      absorbed = .true.
                      ! total n2n XS
                      sign2n_sum = 0.0d0
                      do indx = 1,nsign2n(imix)
                         f = fsign2n(imix,indx)+1
                         if(f .eq. igroup(i))sign2n_sum = sign2n_sum + sign2n(imix,indx)
                      end do
                      siga = sigt(imix,igroup(i)) - sigs_sum
                      ! neutron is converted to the new fission neutron with the weight increased by eta
                      weight(i) = weight(i) * (sigp(imix,igroup(i)) + 2.0*sign2n_sum)/siga
                      if(weight(i) > 0.0)then
                         ! sample group for the new-born neutron by comparing cumulative sum of fission spectrum with random number
                         call random_number(r)
                         tmp = 0.0d0
                         ig = ng + 1
                         do while(tmp <= r)
                            ig = ig - 1
                            tmp = tmp + chi(imix,ig)
                         end do
                         igroup(i) = ig
                      end if
                   end if ! scattering or absorption
                end if ! virtual or real
             end if ! imix <= 0 or imix > 0
          end do ! random walk cycle
      end do ! neutron cycle
      !$omp end parallel do

      ! Russian roulette
      do i = 1,num_neutrons
         if(weight(i) > 0.0d0 .and. weight(i) < 1.0d0 .and. weight(i) < weight0(i))then
            terminate_p = 1.0d0 - weight(i)/weight0(i)
            call random_number(r)
            if(terminate_p >= r)then
               ! kill
               weight(i) = 0.0d0
            else
               ! restore the weight
               weight(i) = weight0(i)
            end if
         end if
      end do

      ! clean up absorbed or killed neutrons
      j = 0
      do i = 1,num_neutrons
         if(weight(i) > 0.0d0)then
            j = j + 1
            xtmp(j) = x(i)
            ytmp(j) = y(i)
            ztmp(j) = z(i)
            itmp(j) = igroup(i)
            wtmp(j) = weight(i)
         end if
      end do
      num_neutrons = j
      x = 0.0
      y = 0.0
      z = 0.0
      igroup = 0
      weight = 0.0
      do i = 1,num_neutrons
         x(i) = xtmp(i)
         y(i) = ytmp(i)
         z(i) = ztmp(i)
         igroup(i) = itmp(i)
         weight(i) = wtmp(i)
      end do
      
      ! split too "heavy" neutrons
      num_new = 0
      do i = 1,num_neutrons
          if(weight(i) > 1.0d0)then
             ! truncated integer value of the neutron weight
             N = floor(weight(i))
             ! sample the number of split neutrons
             call random_number(r)
             if(weight(i)-dble(N) >= r) N = N + 1
             if(N > 1)then
                ! change the weight of the split neutron
                weight(i) = weight(i)/dble(N)
                ! introduce new neutrons
                do j = 1,N-1
                    num_new = num_new + 1
                    x(num_neutrons + num_new) = x(i)
                    y(num_neutrons + num_new) = y(i)
                    z(num_neutrons + num_new) = z(i)
                    weight(num_neutrons + num_new) = weight(i)
                    igroup(num_neutrons + num_new) = igroup(i)
                end do
             end if
          end if
      end do
      ! increase the number of neutrons
      num_neutrons = num_neutrons + num_new

      ! k-eff in a cycle equals the total weight of the new generation over the total weight of the old generation (the old generation weight = num_neutrons_born)
      keff_cycle = 0.0d0
      do i = 1,num_neutrons
         keff_cycle = keff_cycle + weight(i)
      end do
      keff_cycle = keff_cycle / num_neutrons_born

      iactive = icycle - num_cycles_inactive
      if(iactive <= 0)then
         write(*,"('Inactive cycle = ',i4,'/',i4,' | k-eff cycle = ',f8.5,' | num_neutrons = ',i6)")&
               icycle, num_cycles_inactive, keff_cycle, num_neutrons
      else
         ! k-effective of the cycle
         keff_active_cycle(iactive) = keff_cycle
         ! k-effective of the problem
         keff_expected(iactive) = 0.0d0
         do i = 1, iactive
            keff_expected(iactive) = keff_expected(iactive) + keff_active_cycle(i)
         end do
         keff_expected(iactive) = keff_expected(iactive)/iactive
         ! standard deviation of k-effective
         sigma_keff(iactive) = 0.0d0
         do i = 1, iactive
            sigma_keff(iactive) = sigma_keff(iactive) + ( keff_active_cycle(i) - keff_expected(iactive) )**2
         end do
         sigma_keff(iactive) = dsqrt(sigma_keff(iactive) / dble(max(iactive-1,1)) / dble(iactive) )
         write(*,"('Active cycle = ',i6,'/',i6,' | k-eff cycle = ',f8.5,' | num_neutrons = ',i6,&
                 & ' | k-eff expected = ',f9.5,' | sigma = ',f9.5)")&
              & icycle-num_cycles_inactive, num_cycles_active, keff_cycle, num_neutrons,&
              & keff_expected(iactive), sigma_keff(iactive)
      end if
   end do ! batch cycle
   do iz = 1,nz
      do iy = 1,ny
         do ix = 1,nx
            imix = imap(iz,iy,ix)
            if(imix > 0)then
               do ig = 1,ng
                  flux(iz,iy,ix,1,ig) = flux(iz,iy,ix,1,ig) / sigt(imix,ig) / dz(iz-1)
               end do
            end if
         end do
      end do
   end do

else

   ! DIFFUSION SOLVER

   ! relative tolerance
   rtol = 1.0e-8
   ! absolute tolerance
   atol = 1.0e-8
   
   ! initialize fission source
   qf = 1.
   
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
         !$omp private(imix,az_over_v,mlt,dif,sigr,qs,f,t,qn2n,qfis,fluxnew)
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
end if

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