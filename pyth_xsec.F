*CMZ :  6.00/03 01/12/94  14.48.42  BY  HERMES@KAYAK
*CMZ :  6.00/01 26/08/94  17.26.48  BY  HERMES@KAYAK
*CMZ :  5.01/11 04/03/94  16.26.06  BY  DUE
*CMZ :  5.00/04 30/09/93  16.31.43  BY  MICHAEL DUEREN
*CMZ :  5.00/01 24/01/93  13.23.44  BY  MICHAEL DUEREN
*CMZ :  5.00/00 22/01/93  17.29.08  BY  MICHAEL DUEREN
*-- AUTHOR :    MICHAEL DUEREN   22/01/93
*-- AUTHOR2:    PHILIPP GEIGER   25/01/95

	REAL FUNCTION pyth_xsec(x, Q2,dF1, dF2)

	implicit none
	include "mc_set.inc"

*-----DEKLARATION----------------------

	DOUBLE PRECISION DX, DQ2,dy,dnu,DEBEA    ! kinematic variables
	DOUBLE PRECISION DBRACK   ! intermed. factor
	DOUBLE PRECISION pmass,emass2,alpha2,pi
	DOUBLE PRECISION dxsec  ! sigma
	DOUBLE PRECISION df1,df2
	REAL q2, x

! Load up some convenient variables
        pi=3.1415926d0
        emass2=.261112d-6
        alpha2=.729735d-2**2
        pmass=.938272d0
	DEBEA = dble(mcSet_EneBeam)
	DQ2 = DBLE(Q2)
	DX = DBLE(x)

! Check that the kinematics are reasonable
	DY = dq2/(2D0*pmass*dx)/DEBEA
	IF(DX.LE.0.D0.OR.DX.GE.1.D0.OR.
     1     Dy.LE.0.D0.OR.Dy.GE.1.D0.OR.dq2.le.0) then
         goto 999
	endif
! Determine the UNPOLARIZED xsec

	DBRACK = (4.D0*pi*alpha2)/dq2**2
! ........ that's dSigma/dQ2/dx

	dxsec = DBRACK*(dy**2*(1-2D0*emass2/dq2)*dF1+
     1	          (1-dy-dq2/(4D0*DEBEA**2))*dF2/dx)

	pyth_xsec = real(dxsec)

! Leave by the front door

	RETURN

! Leave by the back door

999	continue
	pyth_xsec = 0.0
	return
	end

