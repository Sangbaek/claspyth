
C-----------------------------------------------------------------------
      SUBROUTINE R1998(DQ2,DX,DR,DELTAR)
C-----------------------------------------------------------------------

C new fit to R  hep-ex/9808028 E143 Collab.
C it is based on the old 3 paramter forms
C 0.005<x<0.86, 0.5<Q2<130 GeV2
C E143 x-section measurement 0.03<x<0.4
C with  overall norm uncertainty 2.5%
C
C U. Stoesslein, October 1998
C

      IMPLICIT NONE

      DOUBLE PRECISION DQ2,DX,DR,DELTAR
      DOUBLE PRECISION Q2,Q2MAX,FAC,RLOG,Q2THR
      DOUBLE PRECISION R_A_NEW,R_A,R_B_NEW,R_B,R_C

      DOUBLE PRECISION A(6),B(6),C(6)

      DATA A/ .0485, 0.5470, 2.0621, -.3804,  0.5090, -.0285/
      DATA B/ .0481, 0.6114, -.3509, -.4611,  0.7172, -.0317/
      DATA C/ .0577, 0.4644, 1.8288,12.3708,-43.1043,41.7415/ 

      DOUBLE PRECISION DR1998
      EXTERNAL DR1998

* use R(0.35) if q2 is below 0.35
      Q2=DQ2
      Q2MAX=0.35
      IF(Q2.LT.Q2MAX) Q2=Q2MAX

      FAC   = 1+12.*(Q2/(1.+Q2))*(.125**2/(DX**2+.125**2))
      RLOG  = FAC/LOG(Q2/.04)
      Q2thr = C(4)*DX+C(5)*DX*DX+C(6)*DX*DX*DX

* new additional terms
      R_A_NEW = (1.+A(4)*DX+A(5)*DX*DX)*DX**(A(6))
      R_A   = A(1)*RLOG + A(2)/SQRT(SQRT(Q2**4+A(3)**4))*R_A_NEW
      R_B_NEW = (1.+B(4)*DX+B(5)*DX*DX)*DX**(B(6))
      R_B   = B(1)*RLOG + (B(2)/Q2 + B(3)/(Q2**2+0.3**2))*R_B_NEW
      R_C   = C(1)*RLOG + C(2)/SQRT((Q2-Q2thr)**2+C(3)**2)
      DR    = (R_A+R_B+R_C)/3.

* straight line fit extrapolation to R(Q2=0)=0
      if (dq2.lt.q2max) DR = DR*DQ2/Q2MAX

* I assume the fit uncertainty only for measured Q2 range
      if (Q2 .GT. 0.5) then
         DELTAR = DR1998(DX,Q2)
      else
         DELTAR=DR
      endif

      RETURN
      END

C--------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION DR1998(DX,DQ2)

* Parameterize uncertainty in R1998 
* associated to the fitting procedure only

      IMPLICIT NONE
      DOUBLE PRECISION DX,DQ2

      DR1998 = 0.0078-0.013*DX+(0.070-0.39*DX+0.7*DX*DX)/(1.7+DQ2)

      RETURN
      END
