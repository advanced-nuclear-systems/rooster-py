      subroutine odrv(n, ia,ja,a, p,ip, nsp,isp, path, flag)
      integer  ia(1), ja(1), p(1), ip(1), isp(1), path, flag, v, l, head, tmp, q
      double precision  a(1)
      logical  dflag

      flag = 0
      if (path.lt.1 .or. 5.lt.path)  go to 111

      if ((path-1) * (path-2) * (path-4) .ne. 0)  go to 1
        max = (nsp-n)/2
        v    = 1
        l    = v     +  max
        head = l     +  max
        next = head  +  n
        if (max.lt.n)  go to 110

        call  md(n, ia,ja, max,isp(v),isp(l), isp(head),p,ip, isp(v), flag)
        if (flag.ne.0)  go to 100

   1  if ((path-2) * (path-3) * (path-4) * (path-5) .ne. 0)  go to 2
        tmp = (nsp+1) -      n
        q   = tmp     - (ia(n+1)-1)
        if (q.lt.1)  go to 110

        dflag = path.eq.4 .or. path.eq.5
        call sro(n, ip, ia, ja, a, isp(tmp), isp(q), dflag)

   2  return

 100  return
 110  flag = 10*n + 1
      return
 111  flag = 11*n + 1
      return
      end
