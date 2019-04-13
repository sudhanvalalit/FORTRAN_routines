  ! cubic hermite interpolation
  ! Written: 02/20/2019
  ! By: Sudhanva Lalit
subroutine splinecbh(n,xin,yin,m,xout,yout)
  implicit none
  double precision xin(n),yin(n),xout(m),yout(m),spl(4)
  integer n, index(4), i1, i2, i3, i4, m, stat, knt,i
  double precision x1,x2,x3,xstart,xend,xstep,x(1000),y

  if (m.eq.1) then
     x(1) = xout(1)
  else
     xstart = xin(1); xend = xin(n); xstep = (xend-xstart)/dble(m)
     print*, "xstart",xstart,"xend",xend,xstep
     do i = 1,m
        x(i) = xstart + dble(i-1)*xstep
     end do
  end if

  do i = 1, m
     call cubhermdh (xin,n,x(i),1,spl,index)
     !
     i1=index(1)
     i2=index(2)
     i3=index(3)
     i4=index(4)
     !
     y = spl(1)*yin(i1)+spl(2)*yin(i2)+spl(3)*yin(i3)+spl(4)*yin(i4)
     !
     xout(i) = x(i); yout(i) = y
  end do

end subroutine splinecbh

subroutine cubhermdh (xold,n,xnew,m,spl,index)

  implicit double precision (a-h,o-z)
  !
  !     This subroutine prepares the interpolation of a function
  !     given at the n grid points xold to the m points xnew.
  !     The interpolating functions are cubic hermitian splines.
  !     The first derivatives at the grid points are taken from a
  !     parabola through the actual grid point and its two neighbours.
  !     For the end points the parabola is taken through the two
  !     right or left neighbours, resp.

  !     In the calling routine you still have to take the sum i=1,4
  !     ynew(j)=sum_i spl(j,i)*yold(index(j,i))

  !     by Dirk Hueber, 08.02.1996
  !     Appendix B in
  !     Few-Body Systems 22, 107 (1997)

  !     parameter (nmax=100,mmax=100)
  dimension xold(n),xnew(m),spl(m,4),index(m,4)
  logical enough

  enough=n.ge.3

  !     if (n.gt.nmax) stop'nmax to small in cubherm'
  !     if (m.gt.mmax) stop'mmax to small in cubherm'

  !     evaluation of the indices

  do j=1,m

     !     If xnew(j) is less than xold(1), we extrapolate using the first
     !     two grid points xold(1) and xold(2).

     index(j,2)=1
  end do
  do i=1,n
     do j=1,m
        if (xnew(j).gt.xold(i)) index(j,2)=i
     end do
  end do

  do j=1,m
     index(j,2)=min(index(j,2),n-1)
     index(j,1)=index(j,2)-1
     index(j,3)=index(j,2)+1
     index(j,4)=index(j,2)+2

     !     Indices 1 and 4 are used only for the ewertcation of the derivatives.
     !     The following settings provide the derivatives at the borders of the
     !     grid xold.

     if (index(j,1).eq.0) index(j,1)=3
     if (index(j,4).eq.n+1) index(j,4)=n-2
  end do
  do j=1,m

     !     We don't extrapolate to the right!

     if (xnew(j).le.xold(n).and.enough) then
        i0=index(j,1)
        i1=index(j,2)
        i2=index(j,3)
        i3=index(j,4)
        x0=xold(i0)
        x1=xold(i1)
        x2=xold(i2)
        x3=xold(i3)

        !      Factors for the derivatives

        d10=x1-x0
        d21=x2-x1
        d32=x3-x2
        d20=x2-x0
        d31=x3-x1
        dfak13=(d21/d10-d10/d21)/d20
        dfak14=-d32/(d21*d31)
        dfak23=d10/(d21*d20)
        dfak24=(d32/d21-d21/d32)/d31
        dfak03=-d21/(d10*d20)
        dfak34=d21/(d32*d31)

        !     the cubic hermitian splines

        xn=xnew(j)
        dn1=xn-x1
        d2n=x2-xn
        phidiv=1./(d21*d21*d21)
        phi1=d2n*d2n*phidiv*(d21+2.*dn1)
        phi2=dn1*dn1*phidiv*(d21+2.*d2n)
        phidiv=phidiv*d21*dn1*d2n
        phi3=phidiv*d2n
        phi4=-phidiv*dn1

        !     combining everything to the final factors

        spl(j,2)=phi1+phi3*dfak13+phi4*dfak14
        spl(j,3)=phi2+phi3*dfak23+phi4*dfak24
        spl(j,1)=phi3*dfak03
        spl(j,4)=phi4*dfak34
     else
        spl(j,2)=0.
        spl(j,3)=0.
        spl(j,1)=0.
        spl(j,4)=0.
     endif
  end do
  
end subroutine cubhermdh
