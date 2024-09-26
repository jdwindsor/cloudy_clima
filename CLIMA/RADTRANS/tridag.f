      SUBROUTINE tridag(a,b,c,r,u,n)
      INTEGER n,NMAX
      REAL a(n),b(n),c(n),r(n),u(n)
      PARAMETER (NMAX=500)
      INTEGER j
      REAL bet,gam(NMAX)
      !do j=1,n
      !  print*, a(j),b(j),c(j),r(j),u(j),j !This seems to break the code
      !enddo
      if (b(1).eq.0.) THEN
        PRINT *, 'tridag: rewrite equations' !EWS - rewrite to remove deprecated pause statement 9/4/2015
        STOP 
      endif 
      bet=b(1)
      u(1)=r(1)/bet
      do 11 j=2,n
        !print*,'bet,gam',bet,c(j-1)/bet
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        if (bet.eq.0.) THEN
          PRINT *,  'tridag failed' !EWS - rewrite to remove deprecated pause statement 9/4/2015
          STOP
        endif
        u(j)=(r(j)-a(j)*u(j-1))/bet
11    continue
      do 12 j=n-1,1,-1
        u(j)=u(j)-gam(j+1)*u(j+1)
12    continue
      return
      END
