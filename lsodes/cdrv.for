! driver for subroutines for solving sparse nonsymmetric systems of linear equations (compressed pointer storage)
!
! parameters
! class abbreviations are
!    n - integer variable
!    f - real variable
!    v - supplies a value to the driver
!    r - returns a result from the driver
!    i - used internally by the driver
!    a - array
!
! the nonzero entries of the coefficient matrix m are stored row-by-row in the array a.  to identify the individual nonzero entries in each row, 
! we need to know in which column each entry lies.
! the column indices which correspond to the nonzero entries of m are stored in the array ja.  i.e., if  a(k) = m(i,j),  then ja(k) = j.  
! in addition, we need to know where each row starts and how long it is.
! the index positions in ja and a where the rows of m begin are stored in the array ia.  
! i.e., if m(i,j) is the first nonzero entry (stored) in the i-th row and a(k) = m(i,j), then! ia(i) = k.
! moreover, the index in ja and a of the first location following the last element in the last row is stored in ia(n+1).
! thus, the number of entries in the i-th row is given by ia(i+1) - ia(i), the nonzero entries of the i-th row are stored consecutively in
! a(ia(i)),  a(ia(i)+1),  ..., a(ia(i+1)-1), and the corresponding column indices are stored consecutively in ja(ia(i)), ja(ia(i)+1), ..., ja(ia(i+1)-1).
! for example, the 5 by 5 matrix
!             ( 1. 0. 2. 0. 0.)
!             ( 0. 3. 0. 0. 0.)
!         m = ( 0. 4. 5. 6. 0.)
!             ( 0. 0. 0. 7. 0.)
!             ( 0. 0. 0. 8. 9.)
! would be stored as
!         ia - 1  3  4  7  8 10
!         ja - 1  3  2  2  3  4  4  4  5
!          a - 1. 2. 3. 4. 5. 6. 7. 8. 9.         .
!
! nv    - n     - number of variables/equations.
! fva   - a     - nonzero entries of the coefficient matrix m, stored by rows. size = number of nonzero entries in m.
! nva   - ia    - pointers to delimit the rows in a. size = n+1.
! nva   - ja    - column numbers corresponding to the elements of a. size = size of a.
! fva   - b     - right-hand side b. b and z can the same array. size = n.
! fra   - z     - solution x.  b and z can be the same array. size = n.
!
! the rows and columns of the original matrix m can be reordered (e.g., to reduce fillin or ensure numerical stability) before calling the driver.
! if no reordering is done, then set r(i) = c(i) = ic(i) = i for i=1,...,n.  the solution z is returned in the original order.
! if the columns have been reordered (i.e.,  c(i).ne.i  for some i), then the driver will call a subroutine (nroc) which rearranges
! each row of ja and a, leaving the rows in the original order, but placing the elements of each row in increasing order with respect to the new ordering.
! if  path.ne.1,  then nroc is assumed to have been called already.
!
! nva   - r     - ordering of the rows of m. size = n.
! nva   - c     - ordering of the columns of m. size = n.
! nva   - ic    - inverse of the ordering of the columns of m.  i.e., ic(c(i)) = i  for i=1,...,n. size = n.
!
! the solution of the system of linear equations is divided into three stages:
! nsfc -- the matrix m is processed symbolically to determine where filling will occur during the numeric factorization.
! nnfc -- the matrix m is factored numerically into the product ldu of a unit lower triangular matrix l, a diagonal matrix d, and 
!         a unit upper triangular matrix u, and the system mx = b  is solved.
! nnsc -- the linear system  mx = b  is solved using the ldu or factorization from nnfc.
! nntc -- the transposed linear system  mt x = b  is solved using the ldu factorization from nnf.
!
! for several systems whose coefficient matrices have the same nonzero structure, nsfc need be done only once (for the first system).
! then nnfc is done once for each additional system. for several systems with the same coefficient matrix, nsfc and nnfc need be done only once (for the first system).
! then nnsc or nntc is done once for each additional right-hand side.
!
! nv    - path  - path specification.  values and their meanings are:
!    1  perform nroc, nsfc, and nnfc.
!    2  perform nnfc only  (nsfc is assumed to have been done in a manner compatible with the storage allocation used in the driver).
!    3  perform nnsc only  (nsfc and nnfc are assumed to have been done in a manner compatible with the storage allocation used in the driver).
!    4  perform nntc only  (nsfc and nnfc are assumed to have been done in a manner compatible with the storage allocation used in the driver).
!    5  perform nroc and nsfc.
!
! various errors are detected by the driver and the individual subroutines.
!
! nr  error flag.  values and their meanings are:
!         0     no errors detected
!         n+k   null row in a  --  row = k
!        2n+k   duplicate entry in a  --  row = k
!        3n+k   insufficient storage in nsfc  --  row = k
!        4n+1   insufficient storage in nnfc
!        5n+k   null pivot  --  row = k
!        6n+k   insufficient storage in nsfc  --  row = k
!        7n+1   insufficient storage in nnfc
!        8n+k   zero pivot  --  row = k
!       10n+1   insufficient storage in cdrv
!       11n+1   illegal path specification
!
! working storage is needed for the factored form of the matrix m plus various temporary vectors.
! the arrays isp and rsp should be equivalenced.  integer storage is allocated from the beginning of isp and real storage from the end of rsp.
!
! nv    - nsp   - declared dimension of rsp.  nsp generally must be larger than  8n+2 + 2k  (where  k = (number of nonzero entries in m)).
! nvira - isp   - integer working storage divided up into various arrays needed by the subroutines.  isp and rsp should be equivalenced. size = lratio*nsp.
! fvira - rsp   - real working storage divided up into various arrays needed by the subroutines.  isp and rsp should be equivalenced. size = nsp.
! nr    - esp   - if sufficient storage was available to perform the symbolic factorization (nsfc), then esp is set to the amount of excess storage provided 
!                 (negative if insufficient storage was available to perform the numeric factorization (nnfc)).
!
!
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
