!-----------------------------------------------------------------------
! example problem.
!
! the following is a simple example problem, with the coding
! needed for its solution by lsodes.  the problem is from chemical
! kinetics, and consists of the following 12 rate equations..
!    dy1/dt  = -rk1*y1
!    dy2/dt  = rk1*y1 + rk11*rk14*y4 + rk19*rk14*y5
!                - rk3*y2*y3 - rk15*y2*y12 - rk2*y2
!    dy3/dt  = rk2*y2 - rk5*y3 - rk3*y2*y3 - rk7*y10*y3
!                + rk11*rk14*y4 + rk12*rk14*y6
!    dy4/dt  = rk3*y2*y3 - rk11*rk14*y4 - rk4*y4
!    dy5/dt  = rk15*y2*y12 - rk19*rk14*y5 - rk16*y5
!    dy6/dt  = rk7*y10*y3 - rk12*rk14*y6 - rk8*y6
!    dy7/dt  = rk17*y10*y12 - rk20*rk14*y7 - rk18*y7
!    dy8/dt  = rk9*y10 - rk13*rk14*y8 - rk10*y8
!    dy9/dt  = rk4*y4 + rk16*y5 + rk8*y6 + rk18*y7
!    dy10/dt = rk5*y3 + rk12*rk14*y6 + rk20*rk14*y7
!                + rk13*rk14*y8 - rk7*y10*y3 - rk17*y10*y12
!                - rk6*y10 - rk9*y10
!    dy11/dt = rk10*y8
!    dy12/dt = rk6*y10 + rk19*rk14*y5 + rk20*rk14*y7
!                - rk15*y2*y12 - rk17*y10*y12
!
! with rk1 = rk5 = 0.1,  rk4 = rk8 = rk16 = rk18 = 2.5,
!      rk10 = 5.0,  rk2 = rk6 = 10.0,  rk14 = 30.0,
!      rk3 = rk7 = rk9 = rk11 = rk12 = rk13 = rk19 = rk20 = 50.0,
!      rk15 = rk17 = 100.0.
!
! the t interval is from 0 to 1000, and the initial conditions
! are y1 = 1, y2 = y3 = ... = y12 = 0.  the problem is stiff.
!
! the following coding solves this problem with lsodes, using mf = 121
! and printing results at t = .1, 1., 10., 100., 1000.  it uses
! itol = 1 and mixed relative/absolute tolerance controls.
! during the run and at the end, statistical quantities of interest
! are printed (see optional outputs in the full description below).
!
      double precision atol, rtol, rwork, t, tout, y
      dimension y(12), rwork(500), iwork(30)
      data lrw/500/, liw/30/

      neq = 12
      do i = 1,neq
         y(i) = 0.0d0
      end do
      y(1) = 1.0d0
      t = 0.0d0
      tout = 0.1d0
      itol = 1
      rtol = 1.0d-4
      atol = 1.0d-6
      itask = 1
      istate = 1
      iopt = 0
      mf = 121
      do iout = 1,5
        call lsodes(neq, y, t, tout, itol, rtol, atol, itask, istate, iopt, rwork, lrw, iwork, liw, mf)
        write(6,30)t,iwork(11),rwork(11),(y(i),i=1,neq)
  30    format(//7h at t =,e11.3,4x,
     +    12h no. steps =,i5,4x,12h last step =,e11.3/
     +    13h  y array =  ,4e14.5/13x,4e14.5/13x,4e14.5)
        if(istate .lt. 0)then
           write(6,90)istate
 90        format(///22h error halt.. istate =,i3)
           stop
        end if
        tout = tout*10.0d0
      end do
      lenrw = iwork(17)
      leniw = iwork(18)
      nst = iwork(11)
      nfe = iwork(12)
      nje = iwork(13)
      nlu = iwork(21)
      nnz = iwork(19)
      nnzlu = iwork(25) + iwork(26) + neq
      write (6,70) lenrw,leniw,nst,nfe,nje,nlu,nnz,nnzlu
 70   format(//22h required rwork size =,i4,15h   iwork size =,i4/
     +   12h no. steps =,i4,12h   no. f-s =,i4,12h   no. j-s =,i4,
     +   13h   no. lu-s =,i4/23h no. of nonzeros in j =,i5,
     +   26h   no. of nonzeros in lu =,i5)
      stop
      end
     
      subroutine rhs(neq, t, y, ydot)
      double precision t, y, ydot
      double precision rk1, rk2, rk3, rk4, rk5, rk6, rk7, rk8, rk9, rk10, rk11, rk12, rk13, rk14, rk15, rk16, rk17
      dimension y(12), ydot(12)
      data rk1/0.1d0/, rk2/10.0d0/, rk3/50.0d0/, rk4/2.5d0/, rk5/0.1d0/,
     +   rk6/10.0d0/, rk7/50.0d0/, rk8/2.5d0/, rk9/50.0d0/, rk10/5.0d0/,
     +   rk11/50.0d0/, rk12/50.0d0/, rk13/50.0d0/, rk14/30.0d0/,
     +   rk15/100.0d0/, rk16/2.5d0/, rk17/100.0d0/, rk18/2.5d0/,
     +   rk19/50.0d0/, rk20/50.0d0/
      ydot(1)  = -rk1*y(1)
      ydot(2)  = rk1*y(1) + rk11*rk14*y(4) + rk19*rk14*y(5) - rk3*y(2)*y(3) - rk15*y(2)*y(12) - rk2*y(2)
      ydot(3)  = rk2*y(2) - rk5*y(3) - rk3*y(2)*y(3) - rk7*y(10)*y(3) + rk11*rk14*y(4) + rk12*rk14*y(6)
      ydot(4)  = rk3*y(2)*y(3) - rk11*rk14*y(4) - rk4*y(4)
      ydot(5)  = rk15*y(2)*y(12) - rk19*rk14*y(5) - rk16*y(5)
      ydot(6)  = rk7*y(10)*y(3) - rk12*rk14*y(6) - rk8*y(6)
      ydot(7)  = rk17*y(10)*y(12) - rk20*rk14*y(7) - rk18*y(7)
      ydot(8)  = rk9*y(10) - rk13*rk14*y(8) - rk10*y(8)
      ydot(9)  = rk4*y(4) + rk16*y(5) + rk8*y(6) + rk18*y(7)
      ydot(10) = rk5*y(3) + rk12*rk14*y(6) + rk20*rk14*y(7) + rk13*rk14*y(8) - rk7*y(10)*y(3) - rk17*y(10)*y(12) - rk6*y(10) - rk9*y(10)
      ydot(11) = rk10*y(8)
      ydot(12) = rk6*y(10) + rk19*rk14*y(5) + rk20*rk14*y(7) - rk15*y(2)*y(12) - rk17*y(10)*y(12)
      return
      end
     
      subroutine fjac(neq, t, y, j, ia, ja, pdj)
      double precision t, y, pdj
      double precision rk1, rk2, rk3, rk4, rk5, rk6, rk7, rk8, rk9, rk10, rk11, rk12, rk13, rk14, rk15, rk16, rk17
      dimension y(100), ia(100), ja(100), pdj(100)
      data rk1/0.1d0/, rk2/10.0d0/, rk3/50.0d0/, rk4/2.5d0/, rk5/0.1d0/,
     +   rk6/10.0d0/, rk7/50.0d0/, rk8/2.5d0/, rk9/50.0d0/, rk10/5.0d0/,
     +   rk11/50.0d0/, rk12/50.0d0/, rk13/50.0d0/, rk14/30.0d0/,
     +   rk15/100.0d0/, rk16/2.5d0/, rk17/100.0d0/, rk18/2.5d0/,
     +   rk19/50.0d0/, rk20/50.0d0/
      go to (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), j
 1    pdj(1) = -rk1
      pdj(2) = rk1
      return
 2    pdj(2) = -rk3*y(3) - rk15*y(12) - rk2
      pdj(3) = rk2 - rk3*y(3)
      pdj(4) = rk3*y(3)
      pdj(5) = rk15*y(12)
      pdj(12) = -rk15*y(12)
      return
 3    pdj(2) = -rk3*y(2)
      pdj(3) = -rk5 - rk3*y(2) - rk7*y(10)
      pdj(4) = rk3*y(2)
      pdj(6) = rk7*y(10)
      pdj(10) = rk5 - rk7*y(10)
      return
 4    pdj(2) = rk11*rk14
      pdj(3) = rk11*rk14
      pdj(4) = -rk11*rk14 - rk4
      pdj(9) = rk4
      return
 5    pdj(2) = rk19*rk14
      pdj(5) = -rk19*rk14 - rk16
      pdj(9) = rk16
      pdj(12) = rk19*rk14
      return
 6    pdj(3) = rk12*rk14
      pdj(6) = -rk12*rk14 - rk8
      pdj(9) = rk8
      pdj(10) = rk12*rk14
      return
 7    pdj(7) = -rk20*rk14 - rk18
      pdj(9) = rk18
      pdj(10) = rk20*rk14
      pdj(12) = rk20*rk14
      return
 8    pdj(8) = -rk13*rk14 - rk10
      pdj(10) = rk13*rk14
      pdj(11) = rk10
 9    return
 10   pdj(3) = -rk7*y(3)
      pdj(6) = rk7*y(3)
      pdj(7) = rk17*y(12)
      pdj(8) = rk9
      pdj(10) = -rk7*y(3) - rk17*y(12) - rk6 - rk9
      pdj(12) = rk6 - rk17*y(12)
 11   return
 12   pdj(2) = -rk15*y(2)
      pdj(5) = rk15*y(2)
      pdj(7) = rk17*y(10)
      pdj(10) = -rk17*y(10)
      pdj(12) = -rk15*y(2) - rk17*y(10)
      return
      end
