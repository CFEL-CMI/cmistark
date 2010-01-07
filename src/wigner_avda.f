      REAL*8 FUNCTION DJM(l, m1, m2, beta)

C     Returns Wigner's "small d" function d^{j}_{m1, m2}(cos(beta))
C     computed by  the  explicit formula [Eq. (3.74)] in
C     Biedenharn and Louck (Encyclopedia of mathematics Vol. 8).
C
C     Angle beta in radians

      implicit real*8(a-h,o-z)
      logical first
      parameter (kmax = 50)
      common /bin/  binom(0:kmax, 0:kmax), jmax
      data first/.true./
      save

      if (first) then
         IF (jmax .EQ. 0) call pascal   ! set jmax and binom
         pi    = acos(-1d0)
         first = .false.
      endif

      if ( 2*l .gt. jmax) then
        write(*,'(/a,i6,/a,/a)')
     .      ' *** Error: array binom too small, size = ', jmax,
     .      '     You can increase kmax in subroutine PASCAL,',
     .      '     but beware of overflow in binomial.'
        stop 16
      endif

      beta_r = beta

      k = min(l + m2, l - m2, l + m1, l - m1)   ! 0 .le. k .le. l
      if (k .eq. l + m2) then
          mu    =  m1 - m2
          nu    = -m1 - m2
          labda =  m1 - m2
      elseif (k .eq. l - m2) then
          mu    = m2 - m1
          nu    = m1 + m2
          labda = 0
      elseif (k .eq. l + m1) then
          mu    =  m2 - m1
          nu    = -m1 - m2
          labda = 0
      elseif (k .eq. l - m1) then
          mu    = m1 - m2
          nu    = m1 + m2
          labda = m1 - m2
      end if

      x  = cos(beta_r)
      call jacobi(k, dble(mu), dble(nu), x, dx)

      prefac = sqrt(binom(2*l-k, k+mu)/binom(k+nu,nu))
      djm    = (-1)**labda * prefac * (sin(beta_r/2d0))**mu
     .        * (cos(beta_r/2d0))**nu * dx

      END

      SUBROUTINE JACOBI ( n, alfa, beta, x, dx )
C
C***********************************************************************
C
C  JACOBI evaluates the Jacobi polynomials at X.
C
C
C  Differential equation:
C
C    (1-X*X) Y'' + (BETA-ALFA-(ALFA+BETA+2) X) Y' + N (N+ALFA+BETA+1) Y = 0
C
C  Recursion:
C
C    P(0,ALFA,BETA) = 1,
C
C    P(1,ALFA,BETA) = (1+0.5*(ALFA+BETA))*X + 0.5*(ALFA-BETA)
C
C    2*N*(N+ALFA+BETA)*(2*N-2+ALFA+BETA) * P(N,ALFA,BETA)  =
C    (2*N+ALFA+BETA-1)*
C    ((ALFA**2-BETA**2)+(2*N+ALFA+BETA)*(2*N+ALFA+BETA-2)*X)
C    * P(N-1,ALFA,BETA)
C    -2*(N-1+ALFA)*(N-1+BETA)*(2*I+ALFA+BETA) * P(N-2,ALFA,BETA)
C
C  Restrictions:
C
C    ALFA > -1
C    BETA > -1
C
C  Special values:
C
C    P(N,ALFA,BETA)(1) = (N+ALFA)!/(N!*ALFA!) for integer ALFA.
C
C  Modified:
C
C    14 April 1999
C
C  Author:
C
C    John Burkardt
C
C  Parameters:
C
C    Input, integer N, the highest order polynomial to compute.  Note
C    that polynomials 0 through N will be computed.
C
C    Input, real ALFA, one of the parameters defining the Jacobi
C    polynomials, ALFA must be greater than -1.
C
C    Input, real BETA, the second parameter defining the Jacobi
C    polynomials, BETA must be greater than -1.
C
C    Input, real X, the point at which the polynomials are to be
C    evaluated.
C
C    Output, real CX(0:N), the values of the first N+1 Jacobi
C    polynomials at the point X.  *** see modification ***

************************************************************************
*     Modification by PES Wormer, 15 July 2004.                        *
*     Routine is from POLPACK, made into double precision              *
*     and Fortran 77.                                                  *
*     Introduction array cx(0:nmax). Returns  element dx = cx(n).      *
*                                                                      *
*     From:                                                            *
*     http://www.psc.edu/~burkardt/src/polpak/polpak.html              *
************************************************************************

      parameter (nmax = 100)
      integer n

      real*8 alfa
      real*8 beta
      real*8 cx(0:nmax)
      real*8 c1
      real*8 c2
      real*8 c3
      real*8 c4
      integer i
      real*8 x
      real*8 dx

      if ( alfa .le. -1.0D+00 ) then
        write ( *, * ) ' '
        write ( *, * ) 'Jacobi - Fatal error!'
        write ( *, * ) '  Illegal input value of ALFA = ', alfa
        write ( *, * ) '  But ALFA must be greater than -1.'
        stop
      end if

      if ( beta .le. -1.0D+00 ) then
        write ( *, * ) ' '
        write ( *, * ) 'Jacobi - Fatal error!'
        write ( *, * ) '  Illegal input value of BETA = ', beta
        write ( *, * ) '  But BETA must be greater than -1.'
        stop
      end if

      if ( n .gt. nmax ) then
        write ( *, * ) ' '
        write ( *, * ) 'Jacobi - Fatal error!'
        write ( *, * ) '  Illegal input value of n = ', n
        write ( *, * ) '  But n must be less  than ', nmax+1
        stop
      end if

      if ( n .lt. 0 ) then
        return
      end if

      cx(0) = 1.0D+00

      if ( n .eq. 0 ) then
        dx = cx(0)
        return
      endif

      cx(1) = ( 1.0D+00 + 0.5D+00 * ( alfa + beta ) ) * x
     .    + 0.5D+00 * ( alfa - beta )

      do i = 2, n

        c1 = 2.D+00 * dble ( i ) * ( dble ( i ) + alfa + beta )
     .    * ( dble ( 2 * i - 2 ) + alfa + beta )

        c2 = ( dble ( 2 * i - 1 ) + alfa + beta )
     .    * ( dble ( 2 * i ) + alfa + beta )
     .    * ( dble ( 2 * i - 2 ) + alfa + beta )

        c3 = ( dble ( 2 * i - 1 ) + alfa + beta )
     .    * ( alfa + beta ) * ( alfa - beta )

        c4 = - dble ( 2 ) * ( dble ( i - 1 ) + alfa )
     .    * ( dble ( i - 1 ) + beta ) * ( dble ( 2 * i ) + alfa + beta)

        cx(i) = ( ( c3 + c2 * x ) * cx(i-1) + c4 * cx(i-2) ) / c1

      enddo

      dx = cx(n)

      END

      SUBROUTINE PASCAL

C     Fill common bin with binomial coefficients (Pascal's triangle)

      implicit double precision(a-h,o-z)
      parameter (kmax = 50)
      common /bin/  binom(0:kmax, 0:kmax), jmax

      if ( kmax .eq. jmax ) return  ! bin has been set

C     Set jmax as sign that binom has been set:
      jmax = kmax

      binom(0,0) = 1.d0
      do  i=1,jmax
         binom(i,0) = 1.d0
         binom(i,i) = 1.d0
         do  j=1,i-1
            binom(i,j) = binom(i-1,j-1) + binom(i-1,j)
         enddo
      enddo

      end
