!
! Fortran 95 solver of eigenvalue probem
!
subroutine solve_eigenvalue_problem(meth, geom, nz, nx, ny, nt, ng, nmix, &
                                  & flux, flux_a, imap, &
                                  & sigt, sigtra, sigp, &
                                  & nsigsn, fsigsn, tsigsn, sigsn, &
                                  & nsign2n, fsign2n, tsign2n, sign2n, &
                                  & chi, pitch, dz, keff, keff_a)

use omp_lib
 
implicit none

! method flag ('MC', 'DI')
character*2 meth
! geometry flag ('squar', 'hex01', 'hex06', 'hex24')
character*5 geom
! number of nodes in z, y, x dimensions and number of triangles per hexagon, number of groups and number of mixes
integer nz, nx, ny, nt, ng, nmix
! neutron flux array: flux(nz,nx,ny,nt,ng)
real*8 flux(:,:,:,:,:)
! adjoint flux array: flux_a(nz,nx,ny,nt,ng)
real*8 flux_a(:,:,:,:,:)
! map of material indexes: imap(nz,nx,ny)
integer imap(:,:,:)
! total cross section: sigt(nmix,ng)
real*8 sigt(:,:)
! transport cross section: sigtra(nmix,ng)
real*8 sigtra(:,:)
! production cross section: sigp(nmix,ng)
real*8 sigp(:,:)
! number of entries in full scattering cross section matrix: nsigsn(nmix)
integer nsigsn(:)
! index of energy group from which scattering occurs: fsigsn(nlgndr,nmix,max(nsigsn))
integer fsigsn(:,:,:)
! index of energy group to which scattering occurs: tsigsn(nlgndr,nmix,max(nsigsn))
integer tsigsn(:,:,:)
! full scattering cross section matrix: sigsn(nlgndr,nmix,max(nsigsn)))
real*8 sigsn(:,:,:)
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
integer iy, ix, iz, it, ig, indx, i, j
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
real*8 keff(:)
! adjoint multiplication factor: keff_a
real*8 keff_a(:)
! multiplication factor for evaluating the difference from the previous iteration
real*8 knew
! multiplier at flux
real*8 mlt
! fission source
real*8 qf(nz,nx,ny,nt)
! fission source for adjoint problem
real*8 qff(nz,nx,ny,nt)
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
! sum of sigp for keff_a
real*8 sum_sigp

! Monte Carlo specific variables

! number of neutrons at the end of cycle and neutrons born at the beginning of cycle, number of inactive source cycles to skip before starting k-eff accumulation and number of active cycles for k-eff accumulation (to be done input parameters)
integer num_neutrons, num_neutrons_born, num_cycles_inactive, num_cycles_active
parameter(num_neutrons_born = 50000, num_cycles_inactive = 100, num_cycles_active = 1000)

! cycle index
integer icycle
! group number of neutron
integer igroup(num_neutrons_born*2)
! number of Legendre coefficients
integer :: nlgndr = 8
! number of scattering cosine bins
integer num_mubin
parameter(num_mubin = 200)
! neutron free path
real*8 free_path
! k-effective in a cycle
real*8 :: keff_active_cycle(num_cycles_active) = 1.0
! averaged k-effective
real*8 :: keff_expected(num_cycles_active) = 1.0
! scattering cosine
real*8 mu
! mu dependent pdf of scattering cross section from -> to
real*8 pdf_mu(nmix,num_mubin,ng,ng)
! Legendre coefficients
real*8 pl(8)
! pi number
real*8 :: pi = dacos(-1.0d0)
! k-effective standard deviation
real*8 :: sigma_keff(num_cycles_active) = 0.0d0
! total n2n cross section
real*8 sign2n_f(nmix,ng)
! n2n cross section for group ig
real*8 sign2n_ft(nmix,ng,ng)
! total scattering cross section for group f
real*8 sigs_f(nmix,ng)
! scattering cross section from -> to
real*8 sigs_ft(nmix,ng,ng),sigs1_ft(nmix,ng,ng)
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
real*8 ynode(ny,nx), xnode(nx), zb(nz-1)
! coordinates of neutrons
real*8 y(num_neutrons_born*2), x(num_neutrons_born*2), z(num_neutrons_born*2)

logical absorbed
logical isotropic
logical virtual_collision

! temporal variables
integer itmp(num_neutrons_born*2), num_new, N, iactive, iyp
real*8 r, tmp, xtmp(num_neutrons_born*2), ytmp(num_neutrons_born*2), &
       ztmp(num_neutrons_born*2), wtmp(num_neutrons_born*2), keff_cycle, siga, &
       phi_prime, teta_prime, dmu, dcos_phi, dsin_phi, r1, &
       dir_x_tmp, dir_y_tmp, dir_z_tmp

!$ call OMP_set_dynamic(.true.)

! change from python style to fortran style
imap = imap + 1

! scattering xs from group f to group t for mix imix
pdf_mu = 0.0d0
sigs_ft = 0.0d0
sigs1_ft = 0.0d0
sigs_f = 0.0d0
sign2n_ft = 0.0d0
sign2n_f = 0.0d0
do imix = 1,nmix
   ! scattering matrix
   do indx = 1,nsigsn(imix)
      f = fsigsn(2,imix,indx)+1
      t = tsigsn(2,imix,indx)+1
      sigs1_ft(imix,f,t) = sigs1_ft(imix,f,t) + sigsn(2,imix,indx)
      f = fsigsn(1,imix,indx)+1
      t = tsigsn(1,imix,indx)+1
      sigs_ft(imix,f,t) = sigs_ft(imix,f,t) + sigsn(1,imix,indx)
      ! angular distribution of scattering pdf
      mu = -1.0d0
      dmu = 2.0d0/(num_mubin-1)
      tmp = 0.0d0
      do i = 1,num_mubin
         do j = 1,nlgndr
             if(j == 1)then
                pl(1) = 1.0d0
             else if(j == 2)then
                pl(2) = mu
             else
                pl(j) = ( (2.0d0*dble(j)-1.0d0)*mu*pl(j-1) - (dble(j)-1.0d0)*pl(j-2) )/ dble(j)
             end if
             pdf_mu(imix,i,f,t) = pdf_mu(imix,i,f,t) + (dble(j) - 0.5d0)*sigsn(j,imix,indx)*pl(j)
         end do
         if(pdf_mu(imix,i,f,t) < 0.0d0)pdf_mu(imix,i,f,t) = 0.0d0
         tmp = tmp + pdf_mu(imix,i,f,t)
         mu = mu + dmu
      end do
      do i = 1,num_mubin
         pdf_mu(imix,i,f,t) = pdf_mu(imix,i,f,t)/tmp
      end do
   end do
   do ig = 1,ng
      do j = 1,ng
         sigs_f(imix,ig) = sigs_f(imix,ig) + sigs_ft(imix,ig,j)
      end do
   end do

   ! n2n matrix
   do indx = 1,nsign2n(imix)
      f = fsign2n(imix,indx)+1
      t = tsign2n(imix,indx)+1
      sign2n_ft(imix,f,t) = sign2n_ft(imix,f,t) + sign2n(imix,indx)
   end do
   do ig = 1,ng
      do j = 1,ng
         sign2n_f(imix,ig) = sign2n_f(imix,ig) + sign2n_ft(imix,ig,j)
      end do
   end do
end do
!open(1967,file='angular.txt')
!do imix = 1,nmix
!   do indx = 1,nsigsn(imix)
!      f = fsigsn(1,imix,indx)+1
!      t = tsigsn(1,imix,indx)+1
!      do i=1,num_mubin
!         if(pdf_mu(imix,i,f,t) .ne. 0.0d0)write(1967,*)imix,f,t,pdf_mu(imix,i,f,t)
!      end do
!      write(1967,*)
!   end do
!end do
!close(1967)

if(meth == 'MC')then

   ! MONTE CARLO SOLVER

   ! prepare flux array
   do iz = 1,nz
      do ix = 1,nx
         do iy = 1,ny
            do ig = 1,ng
               flux(iz,ix,iy,1,ig) = 0.0d0
            end do
         end do
      end do
   end do

   ! find coordinate of the hexagon centres
   do ix = 1,nx
      xnode(ix) = dble(ix-1)*dsqrt(3.0d0)*pitch/2.0d0
      do iy = 1,ny
         if(mod(ix,2) == 0)then ! even row
            ynode(iy,ix) = dble(iy-1)*pitch + pitch/2.0d0
         else ! odd row 
            ynode(iy,ix) = dble(iy-1)*pitch
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
         ! generate random x coordinate
         call random_number(r)
         x(i) = (xnode(nx)+pitch)*r - 0.5d0*pitch
         ix = minloc(dabs(xnode - x(i)),1)
         
         ! generate random y coordinate
         call random_number(r)
         y(i) = (ynode(ny,ix)+pitch)*r - 0.5d0*pitch
         iy = minloc(dabs(ynode(:,ix) - y(i)),1)
         
         ! generate random z coordinate
         call random_number(r)
         z(i) = zb(nz-1)*r
         iz = minloc(dabs(zb - z(i)),1)
         if(z(i) > zb(iz)) iz = iz + 1

         imix = imap(iz,ix,iy)
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
      !$omp private(absorbed,isotropic,virtual_collision,r,free_path,teta,phi,mu, &
      !$omp dir_x,dir_y,dir_z,dir_x_tmp,dir_y_tmp,dir_z_tmp, &
      !$omp phi_prime,teta_prime,dcos_phi,dsin_phi, &
      !$omp ix,iy,iyp,iz,imix,sigv,tmp,ig,siga,j,r1)
      ! loop over neutrons
      do i = 1,num_neutrons
          absorbed = .false.
          isotropic = .true.
          virtual_collision = .false.
          ! neutron random walk cycle: from emission to absorption
          do while(.not. absorbed)
             ! sample free path length according to the Woodcock method
             call random_number(r)
             free_path = -dlog(r)/sigtmax(igroup(i))
             if(.not. virtual_collision)then
                if(isotropic)then
                   ! sample the direction of neutron flight assuming isotropic fission
                   !call random_number(r)
                   !phi = 2.0d0*pi*r
                   !call random_number(r)
                   !teta = pi*r
                   !dir_x = dsin(teta)*dcos(phi)
                   !dir_y = dsin(teta)*dsin(phi)
                   !dir_z = dcos(teta)
                   ! sample the direction of neutron flight assuming anisotropic scattering
                   r = 1.0d0
                   r1 = 1.0d0
                   do while(r**2+r1**2 > 1.0d0)
                      call random_number(r)
                      r = 2.d0*r - 1.0d0
                      call random_number(r1)
                      r1 = 2.d0*r1 - 1.0d0
                   end do
                   dir_x = 2.d0*r**2 + 2.d0*r1**2 - 1.0d0
                   dir_y = r*dsqrt((1.0d0-dir_x**2)/(r**2+r1**2))
                   dir_z = r1*dsqrt((1.0d0-dir_x**2)/(r**2+r1**2))

                   !dcos_phi = dcos(phi)
                   !dsin_phi = dsin(phi)
                else
                   !! sample the direction of neutron flight assuming anisotropic scattering
                   !call random_number(r)
                   !phi_prime = 2.0d0*pi*r
                   !call random_number(r)
                   !teta_prime = pi*r
                   !!teta_prime = dacos(mu)
                   !
                   !tmp = dcos(teta_prime)*(1.0d0 - dcos(phi_prime))
                   !! Rodrigues formula (this is important for convergence not to have explicitly dcos(phi) and dsin(phi) in calculations)
                   !dir_x = dsin(teta+teta_prime)*dcos_phi*dcos(phi_prime) - &
                   !      & dsin_phi*dsin(teta_prime)*dsin(phi_prime) + &
                   !      & dir_x*tmp
                   !dir_y = dsin(teta+teta_prime)*dsin_phi*dcos(phi_prime) + &
                   !      & dcos_phi*dsin(teta_prime)*dsin(phi_prime) + &
                   !      & dir_y*tmp
                   !dir_z = dcos(teta+teta_prime)*dcos(phi_prime) + &
                   !      & dir_z*tmp
                   !if(i == 100)write(*,*)dir_x**2 + dir_y**2 + dir_z**2
                   !! for next calculation of dir_x, dir_y and dir_z by formulas above
                   !teta = dacos(dir_z)
                   !dcos_phi = dir_x/dsin(teta)
                   !dsin_phi = dir_y/dsin(teta)
                   ! sample the direction of neutron flight assuming anisotropic scattering
                   r = 1.0d0
                   r1 = 1.0d0
                   do while(r**2+r1**2 > 1.0d0)
                      call random_number(r)
                      r = 2.d0*r - 1.0d0
                      call random_number(r1)
                      r1 = 2.d0*r1 - 1.0d0
                   end do
                   
                   dir_x_tmp = dir_x*mu + dsqrt(1.0d0 - mu**2)*(r*dir_x*dir_z - r1*dir_y)/dsqrt(r**2+r1**2)/dsqrt(1.d0-dir_z**2)
                   dir_y_tmp = dir_y*mu + dsqrt(1.0d0 - mu**2)*(r*dir_y*dir_z + r1*dir_x)/dsqrt(r**2+r1**2)/dsqrt(1.d0-dir_z**2)
                   dir_z_tmp = dir_z*mu - dsqrt(1.0d0 - mu**2)*r*dsqrt(1.d0-dir_z**2)/dsqrt(r**2+r1**2)
                   dir_x = dir_x_tmp
                   dir_y = dir_y_tmp
                   dir_z = dir_z_tmp
                   !if(i == 100)write(*,*)dir_x**2 + dir_y**2 + dir_z**2
                   
                end if
             end if
             ! fly
             x(i) = x(i) + free_path * dir_x
             y(i) = y(i) + free_path * dir_y
             z(i) = z(i) + free_path * dir_z

             ! find neutron position in structured mesh (ix, iy, iz)
             ix = minloc(dabs(xnode - x(i)),1)
             iy = minloc(dabs(ynode(:,ix) - y(i)),1)
             if(xnode(ix) - x(i) > 0.5d0*pitch/dsqrt(3.0d0))then
                iyp = minloc(dabs(ynode(:,ix-1) - y(i)),1)
                if(dabs(ynode(iyp,ix-1) - y(i)) < dabs(ynode(iy,ix) - y(i)))then
                   ix = ix - 1
                   iy = iyp
                end if
             else if(x(i) - xnode(ix) > 0.5d0*pitch/dsqrt(3.0d0))then
                iyp = minloc(dabs(ynode(:,ix+1) - y(i)),1)
                if(dabs(ynode(iyp,ix+1) - y(i)) < dabs(ynode(iy,ix) - y(i)))then
                   ix = ix + 1
                   iy = iyp
                end if
             end if
             iz = minloc(dabs(zb - z(i)),1)
             if(z(i) > zb(iz)) iz = iz + 1
             
             ! mixture index
             imix = imap(iz,ix,iy)
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
                      flux(iz,ix,iy,1,igroup(i)) = flux(iz,ix,iy,1,igroup(i)) + weight(i)
                   end if
                   ! sample type of the collision: scattering or absorption]
                   call random_number(r)
                   if(sigs_f(imix,igroup(i))/sigt(imix,igroup(i)) >= r)then 
                      ! anisotropic scattering
                      isotropic = .false.

                      ! sample group for the secondary neutron by comparing cumulative sum of partial scattering xs with random number
                      call random_number(r)
                      tmp = 0.0d0
                      ig = ng + 1
                      do while(tmp <= r)
                         ig = ig - 1
                         tmp = tmp + sigs_ft(imix,igroup(i),ig)/sigs_f(imix,igroup(i))
                      end do

                      ! sample scattering angle cosine mu
                      !call random_number(r)
                      !tmp = 0.0d0
                      !j = 0
                      !do while(tmp <= r)
                      !   j = j + 1
                      !   tmp = tmp + pdf_mu(imix,j,igroup(i),ig)
                      !end do
                      !mu = dble(j - 1)*2.0d0/dble(num_mubin - 1) - 1.0d0

                      ! average scattering cosine
                      mu = sigs1_ft(imix,igroup(i),ig)/sigs_ft(imix,igroup(i),ig)
                      
                      ! sample scattering cosine mu isotropically
                      !call random_number(r)
                      !mu = 2.0d0*r - 1.0d0

                      igroup(i) = ig

                   else ! absorption
                      absorbed = .true.
                      siga = sigt(imix,igroup(i)) - sigs_f(imix,igroup(i))
                      ! sample type of the collision: n2n or non-n2n absorption
                      call random_number(r)
                      if(sign2n_f(imix,igroup(i))/siga >= r)then ! isotropic n2n
                         ! neutron is converted to the new neutron with the weight increased by 2
                         weight(i) = weight(i) * 2.0d0
                         ! sample group for the secondary neutrons by comparing cumulative sum of partial n2n xs with random number
                         call random_number(r)
                         tmp = 0.0d0
                         ig = ng + 1
                         do while(tmp <= r)
                            ig = ig - 1
                            tmp = tmp + sign2n_ft(imix,igroup(i),ig)/sign2n_f(imix,igroup(i))
                         end do
                         igroup(i) = ig
                      else
                         ! neutron is converted to the new fission neutron with the weight increased by eta
                         weight(i) = weight(i) * sigp(imix,igroup(i))/(siga - sign2n_f(imix,igroup(i)))
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
                      end if ! n2n or non-n2n absorption
                   end if ! scattering or absorption
                end if ! virtual or real
             end if ! imix <= 0 or imix > 0
          end do ! random walk cycle
      end do ! neutron cycle
      !$omp end parallel do

      ! Russian roulette
      do i = 1,num_neutrons
         if(weight(i) > 0.0d0 .and. weight(i) < 0.5d0 .and. weight(i) < weight0(i))then
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
      do ix = 1,nx
         do iy = 1,ny
            imix = imap(iz,ix,iy)
            if(imix > 0)then
               do ig = 1,ng
                  flux(iz,ix,iy,1,ig) = flux(iz,ix,iy,1,ig) / sigt(imix,ig) / dz(iz-1)
               end do
            end if
         end do
      end do
   end do

else

   ! DIFFUSION SOLVER

   ! verification test homogeneous cube
   ng = 2
   sigtra(1,1) = 0.2468
   sigtra(1,2) = 0.3084
   nsigsn(1) = 1
   fsigsn(1,1,1) = 0
   tsigsn(1,1,1) = 1
   sigsn(1,1,1) = 2.3e-3
   nsign2n(1) = 0
   sigp(1,1) = 2.41*2.42e-4
   sigp(1,2) = 2.41*4.08e-3
   sigt(1,1) = 1.382e-3 + sigsn(1,1,1)
   sigt(1,2) = 5.4869e-3
   chi(1,1) = 1
   chi(1,2) = 0 

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
   
   keff(1) = 1.

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
      do while(.not. converge_flux)
         niteri = niteri + 1
         converge_flux = .true.
         !$omp parallel do default(shared) &
         !$omp private(imix,az_over_v,mlt,dif,sigr,qs,f,t,qn2n,qfis,fluxnew)
         do iz = 1,nz
            do ix = 1,nx
               do iy = 1,ny
                  !$ if(omp_get_thread_num() == 0 .and. iz == 1 .and. ix == 1 .and. iy == 1) nthreads = OMP_get_num_threads()
               
                  ! if (iz, ix, iy) is not a boundary condition node, i.e. not 0 (vac) and not -1 (ref)
                  imix = imap(iz,ix,iy)
                  if(imix > 0)then
                     ! node axial area-to-volume ratio
                     az_over_v = 1./dz(iz-1)
                     do it = 1,nt
                        do ig = 1,ng
                           mlt = 0.
                           dif = 0.
               
                           ! diffusion terms (mlt and dif) in different nodalizations and different directions
                           call mltdif(mlt, dif, iz, ix, iy, it, ig, imix, imap, nz, nx, ny, nt, ng, nmix, &
                                       pitch, dz, sigtra, flux, aside_over_v, az_over_v, geom)
               
                           ! removal xs
                           sigr = sigt(imix,ig)
                           do indx = 1,nsigsn(imix)
                              f = fsigsn(1,imix,indx)+1
                              t = tsigsn(1,imix,indx)+1
                              if(f == ig .and. t == ig)then
                                 sigr = sigr - sigsn(1,imix,indx)
                              end if
                           end do
                              
                           mlt = mlt + sigr
               
                           ! scattering source
                           qs = 0.
                           do indx = 1,nsigsn(imix)
                              f = fsigsn(1,imix,indx)+1
                              t = tsigsn(1,imix,indx)+1
                              if(f .ne. ig .and. t == ig)then
                                 qs = qs + sigsn(1,imix,indx) * flux(iz,ix,iy,it,f)
                              end if
                           end do
                           ! n2n source
                           qn2n = 0.
                           do indx = 1,nsign2n(imix)
                              f = fsign2n(imix,indx)+1
                              t = tsign2n(imix,indx)+1
                              if(f .ne. ig .and. t == ig)then
                                 qn2n = qn2n + 2.*sign2n(imix,indx)*flux(iz,ix,iy,it,f)
                              end if
                           end do
                        
                           ! fission source
                           qfis = chi(imix,ig)*qf(iz,ix,iy,it)/keff(1)
                        
                           ! neutron flux
                           fluxnew = (dif + qs + qn2n + qfis)/mlt
                           if(converge_flux)then
                              converge_flux = abs(fluxnew - flux(iz,ix,iy,it,ig)) < rtol*abs(fluxnew) + atol
                           end if
                           flux(iz,ix,iy,it,ig) = fluxnew
                        end do
                     end do
                  end if
               end do
            end do
         end do
         !$omp end parallel do
      end do
   
      ! calculate node-wise fission source qf
      !$omp parallel do default(shared) private(imix)
      do iz = 1,nz
         do ix = 1,nx
            do iy = 1,ny
               ! if (iz, ix, iy) is not a boundary condition node, i.e. not 0 (vac) and not -1 (ref)
               imix = imap(iz,ix,iy)
               if(imix > 0)then
                  do it = 1,nt
                     qf(iz,ix,iy,it) = 0.
                     do ig = 1,ng
                        qf(iz,ix,iy,it) = qf(iz,ix,iy,it) + sigp(imix,ig)*flux(iz,ix,iy,it,ig)
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
         do ix = 1,nx
            do iy = 1,ny
               ! if (iz, ix, iy) is not a boundary condition node, i.e. not 0 (vac) and not -1 (ref)
               imix = imap(iz,ix,iy)
               if(imix > 0)then
                  do it = 1,nt
                     tfs = tfs + qf(iz,ix,iy,it)
                  end do
               end if
            end do
         end do
      end do
   
      ! new k-effective is the ratio of total fission sources at the current (tfs) and previous (1.0) iterations
      knew = tfs
      converge_k = abs(knew - keff(1)) < rtol*abs(knew) + atol
      keff(1) = knew
   
      if(nthreads == 0)then
         write(*,'("keff: ",f13.6, " | niteri: ", i3, " | nitero: ", i3, " | ")') &
               keff(1),niteri,nitero
      else
         write(*,'("keff: ",f13.6, " | niteri: ", i3, " | nitero: ", i3, " | OMPthreads: ", i3, " | ")') &
               keff(1),niteri,nitero,nthreads
      end if
      if(isnan(keff(1))) stop
   
   end do

   ! ADJOINT DIFFUSION SOLVER

   ! eigenvalue keff_a equal to ratio of total fission source at two iterations.
   ! flux is normalise to total fission source = 1 at previous iteration
   keff_a(1) = 1.

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
      do while(.not. converge_flux)
         niteri = niteri + 1
         converge_flux = .true.
         !$omp parallel do default(shared) &
         !$omp private(imix,az_over_v,mlt,dif,sigr,qs,f,t,qn2n,qfis,fluxnew)
         do iz = 1,nz
            do ix = 1,nx
               do iy = 1,ny
                  !$ if(omp_get_thread_num() == 0 .and. iz == 1 .and. ix == 1 .and. iy == 1) nthreads = OMP_get_num_threads()

                  ! if (iz, ix, iy) is not a boundary condition node, i.e. not 0 (vac) and not -1 (ref)
                  imix = imap(iz,ix,iy)
                  if(imix > 0)then
                     ! node axial area-to-volume ratio
                     az_over_v = 1./dz(iz-1)
                     do it = 1,nt
                        do ig = 1,ng
                           mlt = 0.
                           dif = 0.

                           ! diffusion terms (mlt and dif) in different nodalizations and different directions
                           call mltdif(mlt, dif, iz, ix, iy, it, ig, imix, imap, nz, nx, ny, nt, ng, nmix, &
                                       pitch, dz, sigtra, flux_a, aside_over_v, az_over_v, geom)

                           ! removal xs
                           sigr = sigt(imix,ig)
                           do indx = 1,nsigsn(imix)
                              f = fsigsn(1,imix,indx)+1
                              t = tsigsn(1,imix,indx)+1
                              if(f == ig .and. t == ig)then
                                 sigr = sigr - sigsn(1,imix,indx)
                              end if
                           end do

                           mlt = mlt + sigr

                           ! scattering source
                           qs = 0.
                           do indx = 1,nsigsn(imix)
                              f = fsigsn(1,imix,indx)+1
                              t = tsigsn(1,imix,indx)+1
                              if(t .ne. ig .and. f == ig)then
                                 qs = qs + sigsn(1,imix,indx) * flux_a(iz,ix,iy,it,t)
                              end if
                           end do
                           ! n2n source
                           qn2n = 0.
                           do indx = 1,nsign2n(imix)
                              f = fsign2n(imix,indx)+1
                              t = tsign2n(imix,indx)+1
                              if(t .ne. ig .and. f == ig)then
                                 qn2n = qn2n + 2.*sign2n(imix,indx)*flux_a(iz,ix,iy,it,t)
                              end if
                           end do

                           ! fission source
                           qfis = sigp(imix,ig)*qf(iz,ix,iy,it)/keff_a(1)

                           ! adjoint flux
                           fluxnew = (dif + qs + qn2n + qfis)/mlt
                           if(converge_flux)then
                              converge_flux = abs(fluxnew - flux_a(iz,ix,iy,it,ig)) < rtol*abs(fluxnew) + atol
                           end if
                           flux_a(iz,ix,iy,it,ig) = fluxnew
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
         do ix = 1,nx
            do iy = 1,ny
               ! if (iz, ix, iy) is not a boundary condition node, i.e. not 0 (vac) and not -1 (ref)
               imix = imap(iz,ix,iy)
               if(imix > 0)then
                  do it = 1,nt
                     qf(iz,ix,iy,it) = 0.
                     sum_sigp = 0.
                     do ig = 1,ng
                        qf(iz,ix,iy,it) = qf(iz,ix,iy,it) + chi(imix,ig)*flux_a(iz,ix,iy,it,ig)
                        sum_sigp = sum_sigp + sigp(imix,ig)
                     end do
                     qff(iz,ix,iy,it) = qf(iz,ix,iy,it)*sum_sigp
                  end do
               end if
            end do
         end do
      end do
      !$omp end parallel do

      ! calculate total fission source tfs
      tfs = 0.
      do iz = 1,nz
         do ix = 1,nx
            do iy = 1,ny
               ! if (iz, ix, iy) is not a boundary condition node, i.e. not 0 (vac) and not -1 (ref)
               imix = imap(iz,ix,iy)
               if(imix > 0)then
                  do it = 1,nt
                     tfs = tfs + qff(iz,ix,iy,it)
                  end do
               end if
            end do
         end do
      end do

      ! new k-effective is the ratio of total fission sources at the current (tfs) and previous (1.0) iterations
      knew = tfs
      converge_k = abs(knew - keff_a(1)) < rtol*abs(knew) + atol
      keff_a(1) = knew

      if(nthreads == 0)then
         write(*,'("keff_a: ",f13.6, " | niteri: ", i3, " | nitero: ", i3, " | ")') &
               keff_a(1),niteri,nitero
      else
         write(*,'("keff_a: ",f13.6, " | niteri: ", i3, " | nitero: ", i3, " | OMPthreads: ", i3, " | ")') &
               keff_a(1),niteri,nitero,nthreads
      end if
      if(isnan(keff_a(1))) stop

   end do
end if

! change from fortran style to python style
imap = imap - 1

end subroutine

!--------------------------------------------------------------------------------------------------
! Fortran 95 solver of reactor kinetics problem

subroutine solve_kinetic_problem(keff, geom, nz, nx, ny, nt, ng, nmix, flux, dfidt, imap, &
                               & sigt, sigtra, sigp, &
                               & nsigsn, fsigsn, tsigsn, sigsn, &
                               & nsign2n, fsign2n, tsign2n, sign2n, &
                               & chi, pitch, dz, cdnp, dcdnpdt, betaeff, dnplmb, ndnp)

use omp_lib

implicit none

! geometry flag ('squar', 'hex01', 'hex06', 'hex24')
character*5 geom
! number of nodes in z, y, x dimensions and number of triangles per hexagon, number of groups and number of mixes
integer nz, nx, ny, nt, ng, nmix
! map of material indexes: imap(nz,nx,ny)
integer imap(:,:,:)
! subassembly pitch
real*8 pitch
! axial nodalization: dz(nz-2)
real*8 dz(:)
! for do loop
integer iy, ix, iz, it, ig, indx, i, j, im
! index of mix in the node
integer imix
! node axial area-to-volume ratio
real*8 az_over_v
! node side area-to-volume ratio
real*8 aside_over_v
! neutron flux- array: flux(nz,nx,ny,nt,ng)
real*8 flux(:,:,:,:,:)
! neutron flux derivative array: dfidt(nz,nx,ny,nt,ng)
real*8 dfidt(:,:,:,:,:)
! delayed neutron precursor concentration map (nz,nx,ny,nt,ndnp)
real*8 cdnp(:,:,:,:,:)
! delayed neutron precursor derivative concentration map (nz,nx,ny,nt,ndnp)
real*8 dcdnpdt(:,:,:,:,:)
! cross-sections
! total cross section: sigt(nmix,ng)
real*8 sigt(:,:)
! transport cross section: sigtra(nmix,ng)
real*8 sigtra(:,:)
! removal cross section
real*8 sigr
! production cross section: sigp(nmix,ng)
real*8 sigp(:,:)
! number of entries in full scattering cross section matrix: nsigsn(nmix)
integer nsigsn(:)
! index of energy group from which scattering occurs: fsigsn(nlgndr,nmix,max(nsigsn))
integer fsigsn(:,:,:)
! index of energy group to which scattering occurs: tsigsn(nlgndr,nmix,max(nsigsn))
integer tsigsn(:,:,:)
! full scattering cross section matrix: sigsn(nlgndr,nmix,max(nsigsn)))
real*8 sigsn(:,:,:)
! number of entries in n2n cross section matrix: sign2n(nmix)
integer nsign2n(:)
! index of energy group from which n2n occurs: fsign2n(nmix,max(nsign2n))
integer fsign2n(:,:)
! index of energy group to which n2n occurs: tsign2n(nmix,max(nsign2n))
integer tsign2n(:,:)
! n2n cross section matrix: sign2n(nmix,max(nsign2n)))
real*8 sign2n(:,:)
! from and to indices for scattering matrix
integer f, t
! fission spectrum: chi(nmix,ng)
real*8 chi(:,:)
! sum of contributions from diffusion terms from neighbouring nodes
real*8 dif
! Beta effective for delayed neutron precursors
real*8 betaeff(:)
! Lambda values for precursors (1/s)
real*8 dnplmb(:)
! number of precursors
integer ndnp
! Term for precursor derivative calculation
real*8 prec_term
! multiplication factor: keff
real*8 keff(:)
! multiplier at flux
real*8 mlt
! fission source
real*8 qf(nz,nx,ny,nt)
! total fission source
real*8 tfs
! fission source
real*8 qfis
!diffusion term
real*8 qdif
! Delayed neutron source
real*8 q_delay
! n2n source
real*8 qn2n
! scattering source
real*8 qs
! removal rate
real*8 qr
! neutron velocities
real*8 vel(2)

! verification test homogeneous cube
ng = 2
sigtra(1,1) = 0.2468
sigtra(1,2) = 0.3084
nsigsn(1) = 1
fsigsn(1,1,1) = 0
tsigsn(1,1,1) = 1
sigsn(1,1,1) = 2.3e-3
nsign2n(1) = 0
sigp(1,1) = 2.41*2.42e-4
sigp(1,2) = 2.41*4.08e-3
sigt(1,1) = 1.382e-3 + sigsn(1,1,1)
sigt(1,2) = 5.4869e-3 - 0.369e-4 ! perturbation
chi(1,1) = 1
chi(1,2) = 0
vel(1) = 3.0e7
vel(2) = 2.2e5

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

! change from python style to fortran style
imap = imap + 1

do iz = 1,nz
   do ix = 1,nx
      do iy = 1,ny
         ! if (iz, ix, iy) is not a boundary condition node, i.e. not 0 (vac) and not -1 (ref)
         imix = imap(iz,ix,iy)
         if(imix > 0)then
            ! node axial area-to-volume ratio
            az_over_v = 1./dz(iz)
            do it = 1,nt
               ! fission source
               qf(iz,ix,iy,it) = 0.
               do ig = 1,ng
                  qf(iz,ix,iy,it) = qf(iz,ix,iy,it) + sigp(imix,ig)*flux(iz,ix,iy,it,ig)
               end do
               do ig = 1,ng
                  mlt = 0.
                  dif = 0.
                  ! diffusion terms (mlt and dif) in different nodalizations and different directions
                  call mltdif(mlt, dif, iz, ix, iy, it, ig, imix, imap, nz, nx, ny, nt, ng, nmix, &
                              pitch, dz, sigtra, flux, aside_over_v, az_over_v, geom)
                  qdif = dif - mlt * flux(iz,ix,iy,it,ig)

                  ! removal xs
                  sigr = sigt(imix,ig)
                  do indx = 1,nsigsn(imix)
                     f = fsigsn(1,imix,indx)+1
                     t = tsigsn(1,imix,indx)+1
                     if(f == ig .and. t == ig)then
                        sigr = sigr - sigsn(1,imix,indx)
                     end if
                  end do

                  ! multiply by flux to obtain qr for removal rate
                  qr = sigr * flux(iz,ix,iy,it,ig)

                  ! scattering source
                  qs = 0.
                  do indx = 1,nsigsn(imix)
                     f = fsigsn(1,imix,indx)+1
                     t = tsigsn(1,imix,indx)+1
                     if(f .ne. ig .and. t == ig)then
                        qs = qs + sigsn(1,imix,indx) * flux(iz,ix,iy,it,f)
                     end if
                  end do

                  ! n2n source
                  qn2n = 0.
                  do indx = 1,nsign2n(imix)
                     f = fsign2n(imix,indx)+1
                     t = tsign2n(imix,indx)+1
                     if(f .ne. ig .and. t == ig)then
                        qn2n = qn2n + 2.*sign2n(imix,indx)*flux(iz,ix,iy,it,f)
                     end if
                  end do

                  ! fission source
                  qfis = (1 - sum(betaeff))*chi(imix,ig)*qf(iz,ix,iy,it)/keff(1)

                  ! Delayed neutron source
                  q_delay = 0.
                  do im = 1,ndnp
                     q_delay = q_delay + dnplmb(im)*cdnp(iz,ix,iy,it,im)
                  end do
                  q_delay = q_delay*chi(imix,ig)

                  dfidt(iz,ix,iy,it,ig) = vel(ig)*(qdif + qs + qn2n - qr + qfis + q_delay)
               end do

               do im = 1,ndnp
                  dcdnpdt(iz,ix,iy,it,im) = betaeff(im)*qf(iz,ix,iy,it)/keff(1) - dnplmb(im)*cdnp(iz,ix,iy,it,im)
               end do

            end do
         end if
      end do
   end do
end do

! change from fortran style to python style
imap = imap - 1

end subroutine

!--------------------------------------------------------------------------------------------------
! Diffusion terms (mlt and dif) in different nodalizations and different directions
!
subroutine mltdif(mlt, dif, iz, ix, iy, it, ig, imix, imap, nz, nx, ny, nt, ng, nmix, &
                  pitch, dz, sigtra, flux, aside_over_v, az_over_v, geom)

implicit none
integer iz, ix, iy, it, ig, imix, nz, nx, ny, nt, ng, nmix, imap(nz,nx,ny)
real*8 mlt, dif, pitch, dz(nz), sigtra(nmix,ng), flux(nz,nx,ny,nt,ng), aside_over_v, az_over_v
character*5 geom

real*8 db

! diffusion terms in z direction for all geometries: mlt and dif
call difz(mlt,dif,iz,iz-1,ix,iy,it,ig,imix,imap,nz,nx,ny,nt,ng,nmix,dz,sigtra,flux,az_over_v)
call difz(mlt,dif,iz,iz+1,ix,iy,it,ig,imix,imap,nz,nx,ny,nt,ng,nmix,dz,sigtra,flux,az_over_v)

! square geometry
if(geom == 'squar')then
   ! diffusion terms in xy direction: mlt and dif
   call difxy(mlt,dif,iz,ix,iy-1,it,ig,imix,imap,nz,nx,ny,nt,ng,nmix,pitch,sigtra,flux,aside_over_v)
   call difxy(mlt,dif,iz,ix,iy+1,it,ig,imix,imap,nz,nx,ny,nt,ng,nmix,pitch,sigtra,flux,aside_over_v)
   call difxy(mlt,dif,iz,ix-1,iy,it,ig,imix,imap,nz,nx,ny,nt,ng,nmix,pitch,sigtra,flux,aside_over_v)
   call difxy(mlt,dif,iz,ix+1,iy,it,ig,imix,imap,nz,nx,ny,nt,ng,nmix,pitch,sigtra,flux,aside_over_v)

! hexagonal geometry with 1 node per hexagon
else if(geom == 'hex01')then
   ! diffusion terms in xy direction: mlt and dif
   if(mod(ix,2) == 0)then ! even
      call difxy(mlt,dif,iz,ix-1,iy,it,ig,imix,imap,nz,nx,ny,nt,ng,nmix,pitch,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,ix-1,iy+1,it,ig,imix,imap,nz,nx,ny,nt,ng,nmix,pitch,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,ix,iy-1,it,ig,imix,imap,nz,nx,ny,nt,ng,nmix,pitch,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,ix,iy+1,it,ig,imix,imap,nz,nx,ny,nt,ng,nmix,pitch,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,ix+1,iy,it,ig,imix,imap,nz,nx,ny,nt,ng,nmix,pitch,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,ix+1,iy+1,it,ig,imix,imap,nz,nx,ny,nt,ng,nmix,pitch,sigtra,flux,aside_over_v)
   else ! odd
      call difxy(mlt,dif,iz,ix-1,iy-1,it,ig,imix,imap,nz,nx,ny,nt,ng,nmix,pitch,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,ix-1,iy,it,ig,imix,imap,nz,nx,ny,nt,ng,nmix,pitch,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,ix,iy-1,it,ig,imix,imap,nz,nx,ny,nt,ng,nmix,pitch,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,ix,iy+1,it,ig,imix,imap,nz,nx,ny,nt,ng,nmix,pitch,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,ix+1,iy-1,it,ig,imix,imap,nz,nx,ny,nt,ng,nmix,pitch,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,ix+1,iy,it,ig,imix,imap,nz,nx,ny,nt,ng,nmix,pitch,sigtra,flux,aside_over_v)
   end if

! hexagonal geometry with 6 nodes per hexagon
else if(geom == 'hex06')then
   db = pitch/3.
   ! diffusion terms in xy direction: mlt and dif
   if(it == 1)then
      ! north-east
      call difxy(mlt,dif,iz,ix,iy,6,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      if(mod(ix,2) == 0)then ! even
         call difxy(mlt,dif,iz,ix-1,iy+1,4,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      else ! odd
         call difxy(mlt,dif,iz,ix-1,iy,4,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      end if
      call difxy(mlt,dif,iz,ix,iy,2,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
   else if(it == 2)then
      ! east
      call difxy(mlt,dif,iz,ix,iy,1,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,ix,iy+1,5,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,ix,iy,3,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
   else if(it == 3)then
      ! south-east
      call difxy(mlt,dif,iz,ix,iy,2,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      if(mod(ix,2) == 0)then ! even
         call difxy(mlt,dif,iz,ix+1,iy+1,6,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      else ! odd
         call difxy(mlt,dif,iz,ix+1,iy,6,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      end if
      call difxy(mlt,dif,iz,ix,iy,4,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
   else if(it == 4)then
      ! south-west
      call difxy(mlt,dif,iz,ix,iy,3,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      if(mod(ix,2) == 0)then ! even
         call difxy(mlt,dif,iz,ix+1,iy,1,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      else ! odd
         call difxy(mlt,dif,iz,ix+1,iy-1,1,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      end if
      call difxy(mlt,dif,iz,ix,iy,5,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
   else if(it == 5)then
      ! west
      call difxy(mlt,dif,iz,ix,iy,4,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,ix,iy-1,2,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,ix,iy,6,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
   else if(it == 6)then
      ! north-west
      call difxy(mlt,dif,iz,ix,iy,5,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      if(mod(ix,2) == 0)then ! even
         call difxy(mlt,dif,iz,ix-1,iy,3,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      else ! odd
         call difxy(mlt,dif,iz,ix-1,iy-1,3,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      end if
      call difxy(mlt,dif,iz,ix,iy,1,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
   end if

! hexagonal geometry with 24 nodes per hexagon
else if(geom == 'hex24')then
   db = pitch/6.
   ! diffusion terms in xy direction: mlt and dif
   if(it == 1)then
      if(mod(ix,2) == 0)then ! even
         call difxy(mlt,dif,iz,ix-1,iy,19,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      else ! odd
         call difxy(mlt,dif,iz,ix-1,iy-1,19,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      end if
      call difxy(mlt,dif,iz,ix,iy,2,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,ix,iy,7,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
   else if(it == 2)then
      call difxy(mlt,dif,iz,ix,iy,1,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      if(mod(ix,2) == 0)then ! even
         call difxy(mlt,dif,iz,ix-1,iy+1,21,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      else ! odd
         call difxy(mlt,dif,iz,ix-1,iy,21,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      end if
      call difxy(mlt,dif,iz,ix,iy,3,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
   else if(it == 3)then
      call difxy(mlt,dif,iz,ix,iy,2,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,ix,iy,4,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,ix,iy,9,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
   else if(it == 4)then
      call difxy(mlt,dif,iz,ix,iy,3,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      if(mod(ix,2) == 0)then ! even
         call difxy(mlt,dif,iz,ix-1,iy+1,23,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      else ! odd
         call difxy(mlt,dif,iz,ix-1,iy,23,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      end if
      call difxy(mlt,dif,iz,ix,iy,5,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
   else if(it == 5)then
      call difxy(mlt,dif,iz,ix,iy,4,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,ix,iy+1,13,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,ix,iy,11,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
   else if(it == 6)then
      if(mod(ix,2) == 0)then ! even
         call difxy(mlt,dif,iz,ix-1,iy,24,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      else ! odd
         call difxy(mlt,dif,iz,ix-1,iy-1,24,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      end if
      call difxy(mlt,dif,iz,ix,iy,7,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,ix,iy,13,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
   else if(it == 7)then
      call difxy(mlt,dif,iz,ix,iy,6,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,ix,iy,1,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,ix,iy,8,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
   else if(it == 8)then
      call difxy(mlt,dif,iz,ix,iy,7,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,ix,iy,9,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,ix,iy,15,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
   else if(it == 9)then
      call difxy(mlt,dif,iz,ix,iy,8,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,ix,iy,3,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,ix,iy,10,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
   else if(it == 10)then
      call difxy(mlt,dif,iz,ix,iy,9,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,ix,iy,11,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,ix,iy,17,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
   else if(it == 11)then
      call difxy(mlt,dif,iz,ix,iy,10,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,ix,iy,5,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,ix,iy,12,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
   else if(it == 12)then
      call difxy(mlt,dif,iz,ix,iy,11,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,ix,iy+1,20,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,ix,iy,19,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
   else if(it == 13)then
      call difxy(mlt,dif,iz,ix,iy-1,5,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,ix,iy,6,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,ix,iy,14,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
   else if(it == 14)then
      call difxy(mlt,dif,iz,ix,iy,13,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,ix,iy,15,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,ix,iy,20,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
   else if(it == 15)then
      call difxy(mlt,dif,iz,ix,iy,14,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,ix,iy,8,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,ix,iy,16,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
   else if(it == 16)then
      call difxy(mlt,dif,iz,ix,iy,15,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,ix,iy,17,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,ix,iy,22,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
   else if(it == 17)then
      call difxy(mlt,dif,iz,ix,iy,16,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,ix,iy,10,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,ix,iy,18,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
   else if(it == 18)then
      call difxy(mlt,dif,iz,ix,iy,17,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,ix,iy,19,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,ix,iy,24,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
   else if(it == 19)then
      call difxy(mlt,dif,iz,ix,iy,18,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,ix,iy,12,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      if(mod(ix,2) == 0)then ! even
         call difxy(mlt,dif,iz,ix+1,iy+1,1,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      else ! odd
         call difxy(mlt,dif,iz,ix+1,iy,1,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      end if
   else if(it == 20)then
      call difxy(mlt,dif,iz,ix,iy-1,12,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,ix,iy,14,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,ix,iy,21,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
   else if(it == 21)then
      call difxy(mlt,dif,iz,ix,iy,20,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,ix,iy,22,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      if(mod(ix,2) == 0)then ! even
         call difxy(mlt,dif,iz,ix+1,iy,2,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      else ! odd
         call difxy(mlt,dif,iz,ix+1,iy-1,2,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      end if
   else if(it == 22)then
      call difxy(mlt,dif,iz,ix,iy,21,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,ix,iy,16,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,ix,iy,23,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
   else if(it == 23)then
      call difxy(mlt,dif,iz,ix,iy,22,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,ix,iy,24,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      if(mod(ix,2) == 0)then ! even
         call difxy(mlt,dif,iz,ix+1,iy,4,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      else ! odd
         call difxy(mlt,dif,iz,ix+1,iy-1,4,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      end if
   else if(it == 24)then
      call difxy(mlt,dif,iz,ix,iy,23,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      call difxy(mlt,dif,iz,ix,iy,18,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      if(mod(ix,2) == 0)then ! even
         call difxy(mlt,dif,iz,ix+1,iy+1,6,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
      else ! odd
         call difxy(mlt,dif,iz,ix+1,iy,6,ig,imix,imap,nz,nx,ny,nt,ng,nmix,db,sigtra,flux,aside_over_v)
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
subroutine difxy(mlt, dif, jz, jx, jy, jt, ig, imix, imap, nz, nx, ny, nt, ng, nmix, p, sigtra, flux, a_over_v)

implicit none

integer jz, jy, jx, jt, ig, imix, nz, nx, ny, nt, ng, nmix, imap(nz,nx,ny)
real*8 mlt, dif, p, sigtra(nmix,ng), flux(nz,nx,ny,nt,ng), a_over_v

integer imix_n
real*8 db, D

imix_n = imap(jz,jx,jy)
if(imix_n == 0)then
   ! neighbour is vacuum
   db = 0.5*p !+ 0.71/sigtra(imix, ig)
   D = 1./(3.*sigtra(imix, ig))
   mlt = mlt + D * a_over_v / db
else if(imix_n .ne. -1)then
   ! neighbour is normal node (neither vacuum nor reflective)
   D = 2./(3.*sigtra(imix, ig) + 3.*sigtra(imix_n, ig))
   mlt = mlt + D * a_over_v / p
   dif = dif + D * flux(jz,jx,jy,jt,ig) * a_over_v / p
end if

end subroutine

!--------------------------------------------------------------------------------------------------
! Diffusion terms (mlt and dif) in z directions
!
subroutine difz(mlt, dif, iz, jz, jx, jy, jt, ig, imix, imap, nz, nx, ny, nt, ng, nmix, dz, sigtra, flux, a_over_v)

implicit none

integer iz, jz, jy, jx, jt, ig, imix, nz, nx, ny, nt, ng, nmix, imap(nz,nx,ny)
real*8 mlt, dif, dz(nz), sigtra(nmix,ng), flux(nz,nx,ny,nt,ng), a_over_v

! index of mix in the neighbouring node
integer imix_n
real*8 db, D

imix_n = imap(jz,jx,jy)
if(imix_n == 0)then
   ! neighbour is vacuum
   db = 0.5*dz(iz) !+ 0.71/sigtra(imix, ig)
   D = 1./(3.*sigtra(imix, ig))
   mlt = mlt + D * a_over_v / db
else if(imix_n .ne. -1)then
   ! neighbour is normal node (neither vacuum nor reflective)
   db = 0.5*(dz(iz) + dz(jz))
   D = 2.*db/(3.*sigtra(imix, ig)*dz(iz) + 3.*sigtra(imix_n, ig)*dz(jz))
   mlt = mlt + D * a_over_v / db
   dif = dif + D * flux(jz,jx,jy,jt,ig) * a_over_v / db
end if

end subroutine