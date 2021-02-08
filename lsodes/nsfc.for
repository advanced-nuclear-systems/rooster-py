      subroutine nsfc(n, r, ic, ia, ja, jlmax, il, jl, ijl, jumax, iu, ju, iju, q, ira, jra, irac, irl, jrl, iru, jru, flag)
      integer cend, qm, rend, rk, vj
      integer ia(1), ja(1), ira(1), jra(1), il(1), jl(1), ijl(1)
      integer iu(1), ju(1), iju(1), irl(1), jrl(1), iru(1), jru(1)
      integer r(1), ic(1), q(1), irac(1), flag

      np1 = n + 1
      jlmin = 1
      jlptr = 0
      il(1) = 1
      jumin = 1
      juptr = 0
      iu(1) = 1
      do 1 k=1,n
        irac(k) = 0
        jra(k) = 0
        jrl(k) = 0
   1    jru(k) = 0
      do 2 k=1,n
        rk = r(k)
        iak = ia(rk)
        if (iak .ge. ia(rk+1))  go to 101
        jaiak = ic(ja(iak))
        if (jaiak .gt. k)  go to 105
        jra(k) = irac(jaiak)
        irac(jaiak) = k
   2    ira(k) = iak

      do 41 k=1,n

        q(np1) = np1
        luk = -1
        vj = irac(k)
        if (vj .eq. 0)  go to 5
   3      qm = np1
   4      m = qm
          qm =  q(m)
          if (qm .lt. vj)  go to 4
          if (qm .eq. vj)  go to 102
            luk = luk + 1
            q(m) = vj
            q(vj) = qm
            vj = jra(vj)
            if (vj .ne. 0)  go to 3
   5    lastid = 0
        lasti = 0
        ijl(k) = jlptr
        i = k
   6      i = jru(i)
          if (i .eq. 0)  go to 10
          qm = np1
          jmin = irl(i)
          jmax = ijl(i) + il(i+1) - il(i) - 1
          long = jmax - jmin
          if (long .lt. 0)  go to 6
          jtmp = jl(jmin)
          if (jtmp .ne. k)  long = long + 1
          if (jtmp .eq. k)  r(i) = -r(i)
          if (lastid .ge. long)  go to 7
            lasti = i
            lastid = long
   7      do 9 j=jmin,jmax
            vj = jl(j)
   8        m = qm
            qm = q(m)
            if (qm .lt. vj)  go to 8
            if (qm .eq. vj)  go to 9
              luk = luk + 1
              q(m) = vj
              q(vj) = qm
              qm = vj
   9        continue
            go to 6
  10    qm = q(np1)
        if (qm .ne. k)  go to 105
        if (luk .eq. 0)  go to 17
        if (lastid .ne. luk)  go to 11
        irll = irl(lasti)
        ijl(k) = irll + 1
        if (jl(irll) .ne. k)  ijl(k) = ijl(k) - 1
        go to 17
  11    if (jlmin .gt. jlptr)  go to 15
        qm = q(qm)
        do 12 j=jlmin,jlptr
          if (jl(j) - qm)  12, 13, 15
  12      continue
        go to 15
  13    ijl(k) = j
        do 14 i=j,jlptr
          if (jl(i) .ne. qm)  go to 15
          qm = q(qm)
          if (qm .gt. n)  go to 17
  14      continue
        jlptr = j - 1
  15    jlmin = jlptr + 1
        ijl(k) = jlmin
        if (luk .eq. 0)  go to 17
        jlptr = jlptr + luk
        if (jlptr .gt. jlmax)  go to 103
          qm = q(np1)
          do 16 j=jlmin,jlptr
            qm = q(qm)
  16        jl(j) = qm
  17    irl(k) = ijl(k)
        il(k+1) = il(k) + luk

        q(np1) = np1
        luk = -1
        rk = r(k)
        jmin = ira(k)
        jmax = ia(rk+1) - 1
        if (jmin .gt. jmax)  go to 20
        do 19 j=jmin,jmax
          vj = ic(ja(j))
          qm = np1
  18      m = qm
          qm = q(m)
          if (qm .lt. vj)  go to 18
          if (qm .eq. vj)  go to 102
            luk = luk + 1
            q(m) = vj
            q(vj) = qm
  19      continue
  20    lastid = 0
        lasti = 0
        iju(k) = juptr
        i = k
        i1 = jrl(k)
  21      i = i1
          if (i .eq. 0)  go to 26
          i1 = jrl(i)
          qm = np1
          jmin = iru(i)
          jmax = iju(i) + iu(i+1) - iu(i) - 1
          long = jmax - jmin
          if (long .lt. 0)  go to 21
          jtmp = ju(jmin)
          if (jtmp .eq. k)  go to 22
            long = long + 1
            cend = ijl(i) + il(i+1) - il(i)
            irl(i) = irl(i) + 1
            if (irl(i) .ge. cend)  go to 22
              j = jl(irl(i))
              jrl(i) = jrl(j)
              jrl(j) = i
  22      if (lastid .ge. long)  go to 23
            lasti = i
            lastid = long
  23      do 25 j=jmin,jmax
            vj = ju(j)
  24        m = qm
            qm = q(m)
            if (qm .lt. vj)  go to 24
            if (qm .eq. vj)  go to 25
              luk = luk + 1
              q(m) = vj
              q(vj) = qm
              qm = vj
  25        continue
          go to 21
  26    if (il(k+1) .le. il(k))  go to 27
          j = jl(irl(k))
          jrl(k) = jrl(j)
          jrl(j) = k
  27    qm = q(np1)
        if (qm .ne. k)  go to 105
        if (luk .eq. 0)  go to 34
        if (lastid .ne. luk)  go to 28
        irul = iru(lasti)
        iju(k) = irul + 1
        if (ju(irul) .ne. k)  iju(k) = iju(k) - 1
        go to 34
  28    if (jumin .gt. juptr)  go to 32
        qm = q(qm)
        do 29 j=jumin,juptr
          if (ju(j) - qm)  29, 30, 32
  29      continue
        go to 32
  30    iju(k) = j
        do 31 i=j,juptr
          if (ju(i) .ne. qm)  go to 32
          qm = q(qm)
          if (qm .gt. n)  go to 34
  31      continue
        juptr = j - 1
  32    jumin = juptr + 1
        iju(k) = jumin
        if (luk .eq. 0)  go to 34
        juptr = juptr + luk
        if (juptr .gt. jumax)  go to 106
          qm = q(np1)
          do 33 j=jumin,juptr
            qm = q(qm)
  33        ju(j) = qm
  34    iru(k) = iju(k)
        iu(k+1) = iu(k) + luk

        i = k
  35      i1 = jru(i)
          if (r(i) .lt. 0)  go to 36
          rend = iju(i) + iu(i+1) - iu(i)
          if (iru(i) .ge. rend)  go to 37
            j = ju(iru(i))
            jru(i) = jru(j)
            jru(j) = i
            go to 37
  36      r(i) = -r(i)
  37      i = i1
          if (i .eq. 0)  go to 38
          iru(i) = iru(i) + 1
          go to 35

  38    i = irac(k)
        if (i .eq. 0)  go to 41
  39      i1 = jra(i)
          ira(i) = ira(i) + 1
          if (ira(i) .ge. ia(r(i)+1))  go to 40
          irai = ira(i)
          jairai = ic(ja(irai))
          if (jairai .gt. i)  go to 40
          jra(i) = irac(jairai)
          irac(jairai) = i
  40      i = i1
          if (i .ne. 0)  go to 39
  41    continue

      ijl(n) = jlptr
      iju(n) = juptr
      flag = 0
      return

 101  flag = n + rk
      return
 102  flag = 2*n + rk
      return
 103  flag = 3*n + k
      return
 105  flag = 5*n + k
      return
 106  flag = 6*n + k
      return
      end
