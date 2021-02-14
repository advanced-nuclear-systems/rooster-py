
! symbolic ldu-factorization of nonsymmetric sparse matrix (compressed pointer storage)

      subroutine nsfc(n, r, ic, ia, ja, jlmax, il, jl, ijl, jumax, iu, ju, iju, q, ira, jra, irac, irl, jrl, iru, jru, flag)
      integer cend, qm, rend, rk, vj
      integer ia(1), ja(1), ira(1), jra(1), il(1), jl(1), ijl(1)
      integer iu(1), ju(1), iju(1), irl(1), jrl(1), iru(1), jru(1)
      integer r(1), ic(1), q(1), irac(1), flag
      logical exit17, exit 34

      np1 = n + 1
      jlmin = 1
      jlptr = 0
      il(1) = 1
      jumin = 1
      juptr = 0
      iu(1) = 1
      do k=1,n
         irac(k) = 0
         jra(k) = 0
         jrl(k) = 0
         jru(k) = 0
      end do

      do k=1,n
         rk = r(k)
         iak = ia(rk)
         if(iak .ge. ia(rk+1))then
            flag = n + rk
            return
         end if
         jaiak = ic(ja(iak))
         if(jaiak .gt. k)then
            flag = 5*n + k
            return
         end if
         jra(k) = irac(jaiak)
         irac(jaiak) = k
         ira(k) = iak
      end do

      do k=1,n
         q(np1) = np1
         luk = -1
         vj = irac(k)
         if(vj .ne. 0)then
            do while(.true.)
               qm = np1
               do while(.true.)
                  m = qm
                  qm =  q(m)
                  if(qm .ge. vj) exit
               end do
               if(qm .eq. vj)then
                  flag = 2*n + rk
                  return
               end if
               luk = luk + 1
               q(m) = vj
               q(vj) = qm
               vj = jra(vj)
               if(vj .eq. 0) exit
            end do
         end if
         lastid = 0
         lasti = 0
         ijl(k) = jlptr
         i = k
         do while(.true.)
            i = jru(i)
            if(i .eq. 0) exit
            qm = np1
            jmin = irl(i)
            jmax = ijl(i) + il(i+1) - il(i) - 1
            long = jmax - jmin
            if(long .ge. 0)then
               jtmp = jl(jmin)
               if(jtmp .ne. k)  long = long + 1
               if(jtmp .eq. k)  r(i) = -r(i)
               if(lastid .lt. long)then
                  lasti = i
                  lastid = long
               end if
               do j=jmin,jmax
                  vj = jl(j)
                  do while(.true.)
                     m = qm
                     qm = q(m)
                     if(qm .ge. vj) exit
                  end do
                  if(qm .ne. vj)then
                     luk = luk + 1
                     q(m) = vj
                     q(vj) = qm
                     qm = vj
                  end if
               end do
            end if
         end do

         qm = q(np1)
         if(qm .ne. k)then
            flag = 5*n + k
            return
         end if

         exit17 = .false.
         if(luk .ne. 0)then
            if(lastid .eq. luk)then
               irll = irl(lasti)
               ijl(k) = irll + 1
               if(jl(irll) .ne. k)  ijl(k) = ijl(k) - 1
            else
               if(jlmin .le. jlptr)then
                  qm = q(qm)
                  do j=jlmin,jlptr
                     if(jl(j) .ge. qm)then
                        if(jl(j) .eq. qm)then
                           ijl(k) = j
                           do i=j,jlptr
                              if(jl(i) .eq. qm)then
                                 qm = q(qm)
                                 if(qm .gt. n)then
                                    exit17 = .true.
                                    exit
                                 end if
                              else
                                 jlmin = jlptr + 1
                                 ijl(k) = jlmin
                                 if(luk .ne. 0)then
                                    jlptr = jlptr + luk
                                    if(jlptr .gt. jlmax)then
                                       flag = 3*n + k
                                       return
                                    end if
                                    qm = q(np1)
                                    do jj=jlmin,jlptr
                                       qm = q(qm)
                                       jl(jj) = qm
                                    end do
                                 end if
                                 exit17 = .true.
                                 exit
                              end if
                           end do
                           if(.not. exit17) jlptr = j - 1
                        end if
                        if(.not. exit17)then
                           jlmin = jlptr + 1
                           ijl(k) = jlmin
                           if(luk .ne. 0)then
                              jlptr = jlptr + luk
                              if(jlptr .gt. jlmax)then
                                 flag = 3*n + k
                                 return
                              end if
                              qm = q(np1)
                              do jj=jlmin,jlptr
                                 qm = q(qm)
                                 jl(jj) = qm
                              end do
                           end if
                        end if
                        exit17 = .true.
                        exit
                     end if
                  end do
               end if

               if(.not. exit17)then
                  jlmin = jlptr + 1
                  ijl(k) = jlmin
                  if(luk .ne. 0)then
                     jlptr = jlptr + luk
                     if(jlptr .gt. jlmax)then
                        flag = 3*n + k
                        return
                     end if
                     qm = q(np1)
                     do jj=jlmin,jlptr
                        qm = q(qm)
                        jl(jj) = qm
                     end do
                  end if
               end if
            end if
         end if

         irl(k) = ijl(k)
         il(k+1) = il(k) + luk     
         q(np1) = np1
         luk = -1
         rk = r(k)
         jmin = ira(k)
         jmax = ia(rk+1) - 1
         if(jmin .le. jmax)then
            do j=jmin,jmax
               vj = ic(ja(j))
               qm = np1
               do while(.true.)
                  m = qm
                  qm = q(m)
                  if(qm .ge. vj) exit
               end do
               if(qm .eq. vj)then
                  flag = 2*n + rk
                  return
               end if
               luk = luk + 1
               q(m) = vj
               q(vj) = qm
            end do
         end if
         lastid = 0
         lastid = 0
         lasti = 0
         iju(k) = juptr
         i = k
         i1 = jrl(k)

         do while(.true.)
            i = i1
            if(i .eq. 0) exit
            i1 = jrl(i)
            qm = np1
            jmin = iru(i)
            jmax = iju(i) + iu(i+1) - iu(i) - 1
            long = jmax - jmin
            if(long .ge. 0)then
               jtmp = ju(jmin)
               if(jtmp .ne. k)then
                  long = long + 1
                  cend = ijl(i) + il(i+1) - il(i)
                  irl(i) = irl(i) + 1
                  if(irl(i) .lt. cend)then
                     j = jl(irl(i))
                     jrl(i) = jrl(j)
                     jrl(j) = i
                  end if
               end if
               if(lastid .lt. long)then
                  lasti = i
                  lastid = long
               end if
               do j=jmin,jmax
                  vj = ju(j)
                  do while(.true.)
                     m = qm
                     qm = q(m)
                     if(qm .ge. vj) exit
                  end do
                  if(qm .ne. vj)then
                     luk = luk + 1
                     q(m) = vj
                     q(vj) = qm
                     qm = vj
                  end if
               end do
            end if
         end do

         if(il(k+1) .gt. il(k))then
            j = jl(irl(k))
            jrl(k) = jrl(j)
            jrl(j) = k
         end if
         qm = q(np1)
         if(qm .ne. k)then
            flag = 5*n + k
            return
         end if

         exit34 = .false.
         if(luk .ne. 0)then
            if(lastid .eq. luk)then
               irul = iru(lasti)
               iju(k) = irul + 1
               if(ju(irul) .ne. k)  iju(k) = iju(k) - 1
            else
               if(jumin .le. juptr)then
                  qm = q(qm)
                  do j=jumin,juptr
                     if(ju(j) .ge. qm)then
                        if(ju(j) .eq. qm)then
                           iju(k) = j
                           do i=j,juptr
                              if(ju(i) .eq. qm)then
                                 qm = q(qm)
                                 if(qm .gt. n)then
                                    exit34 = .true.
                                    exit
                                 end if
                              else
                                 jumin = juptr + 1
                                 iju(k) = jumin
                                 if(luk .ne. 0)then
                                    juptr = juptr + luk
                                    if(juptr .gt. jumax)then
                                       flag = 3*n + k
                                       return
                                    end if
                                    qm = q(np1)
                                    do jj=jumin,juptr
                                       qm = q(qm)
                                       ju(jj) = qm
                                    end do
                                 end if
                                 exit34 = .true.
                                 exit
                              end if
                           end do
                           if(.not. exit34) juptr = j - 1
                        end if
                        if(.not. exit34)then
                           jumin = juptr + 1
                           iju(k) = jumin
                           if(luk .ne. 0)then
                              juptr = juptr + luk
                              if(juptr .gt. jumax)then
                                 flag = 3*n + k
                                 return
                              end if
                              qm = q(np1)
                              do jj=jumin,juptr
                                 qm = q(qm)
                                 ju(jj) = qm
                              end do
                           end if
                        end if
                        exit34 = .true.
                        exit
                     end if
                  end do
               end if

               if(.not. exit34)then
                  jumin = juptr + 1
                  iju(k) = jumin
                  if(luk .ne. 0)then
                     juptr = juptr + luk
                     if(juptr .gt. jumax)then
                        flag = 3*n + k
                        return
                     end if
                     qm = q(np1)
                     do jj=jumin,juptr
                        qm = q(qm)
                        ju(jj) = qm
                     end do
                  end if
               end if
            end if
         end if

         iru(k) = iju(k)
         iu(k+1) = iu(k) + luk
         i = k
         do while(.true.)
            i1 = jru(i)
            if(r(i) .ge. 0)then
               rend = iju(i) + iu(i+1) - iu(i)
               if(iru(i) .lt. rend)then
                  j = ju(iru(i))
                  jru(i) = jru(j)
                  jru(j) = i
               end if
            else
               r(i) = -r(i)
            end if
            i = i1
            if(i .eq. 0) exit
            iru(i) = iru(i) + 1
         end do
         i = irac(k)
         if(i .ne. 0)then
            do while(.true.)
               i1 = jra(i)
               ira(i) = ira(i) + 1
               if(ira(i) .lt. ia(r(i)+1))then
                  irai = ira(i)
                  jairai = ic(ja(irai))
                  if(jairai .le. i)then
                     jra(i) = irac(jairai)
                     irac(jairai) = i
                  end if
               end if
               i = i1
               if(i .eq. 0) exit
            end do
         end if
      end do
      ijl(n) = jlptr
      iju(n) = juptr
      flag = 0
      return
      end
