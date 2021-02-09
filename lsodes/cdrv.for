      subroutine cdrv(n, r, c, ic, ia, ja, a, b, z, nsp, isp, rsp, esp, path, flag)

      integer r(1), c(1), ic(1), ia(1), ja(1), isp(1), esp, path, flag, d, u, q, row, tmp, ar, umax
      double precision a(1), b(1), z(1), rsp(1)

      if(path .lt. 1 .or. path .gt. 5)then
         flag = 10*n + 1
         RETURN
      end if
!     initialize and divide up temporary storage
      il = 1
      ijl = il + (n+1)
      iu = ijl + n
      iju = iu + (n+1)
      irl = iju + n
      jrl = irl + n
      jl = jrl + n
      if((path-1) * (path-5) .eq. 0)then
         max = (2*nsp + 1 - jl) - (n+1) - 5*n
         jlmax = max/2
         q = jl + jlmax
         ira = q + (n+1)
         jra = ira + n
         irac = jra + n
         iru = irac + n
         jru = iru + n
         jutmp = jru + n
         jumax = 2*nsp + 1 - jutmp
         esp = max/2
         if(jlmax.le.0 .or. jumax.le.0)then
!           error.. insufficient storage
            flag = 10*n + 1
            RETURN
         end if
         do i=1,n
            if(c(i).ne.i)then
               ar = nsp + 1 - n
               call nroc(n, ic, ia, ja, a, isp(il), rsp(ar), isp(iu), flag)
               if(flag .ne. 0)then
!                 error.. in nroc, nsfc, nnfc, or nnsc
                  RETURN 
               end if
               exit
            end if
         end do
         call  nsfc(n, r, ic, ia, ja, jlmax, isp(il), isp(jl), isp(ijl), jumax, isp(iu), isp(jutmp), isp(iju), isp(q), isp(ira), isp(jra), isp(irac), isp(irl), isp(jrl), isp(iru), isp(jru),  flag)
         if(flag .ne. 0)then
!           error.. in nroc, nsfc, nnfc, or nnsc
            RETURN 
         end if
         jlmax = isp(ijl+n-1)
         ju = jl + jlmax
         jumax = isp(iju+n-1)
         if(jumax .gt. 0)then
            do j=1,jumax
               isp(ju+j-1) = isp(jutmp+j-1)
            end do
         end if
      end if
      jlmax = isp(ijl+n-1)
      ju = jl  + jlmax
      jumax = isp(iju+n-1)
      l = (ju + jumax - 2 + 2)/2 + 1
      lmax = isp(il+n) - 1
      d = l + lmax
      u = d + n
      row = nsp + 1 - n
      tmp = row - n
      umax = tmp - u
      esp = umax - (isp(iu+n) - 1)

      if((path-1) * (path-2) .eq. 0)then
         if(umax.lt.0)then
!           error.. insufficient storage
            flag = 10*n + 1
            RETURN
         end if
         call nnfc(n, r, c, ic,  ia, ja, a, z, b, lmax, isp(il), isp(jl), isp(ijl), rsp(l), rsp(d), umax, isp(iu), isp(ju), isp(iju), rsp(u), rsp(row), rsp(tmp),  isp(irl), isp(jrl),  flag)
         if(flag .ne. 0)then
!           error.. in nroc, nsfc, nnfc, or nnsc
            RETURN 
         end if
      end if
      if ((path-3) .eq. 0) call nnsc(n,  r, c, isp(il), isp(jl), isp(ijl), rsp(l), rsp(d), isp(iu), isp(ju), isp(iju), rsp(u), z, b, rsp(tmp))
      if ((path-4) .eq. 0) call nntc(n,  r, c, isp(il), isp(jl), isp(ijl), rsp(l), rsp(d), isp(iu), isp(ju), isp(iju), rsp(u), z, b, rsp(tmp))
      RETURN
      end
