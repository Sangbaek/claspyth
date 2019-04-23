*CMZ :  5.00/00 22/01/93  17.29.08  by  Michael Dueren
*-- Author :    Michael Dueren   22/01/93
*-- Author2:    Philipp Geiger   Okt.94
      SUBROUTINE R1990(DQ2,DX,DR)
*:=====================================================================:
*:                                                                     :
*:      Author:    M.Dueren        Last update: 25.01.1995             :
*:                                 Tested: yes                         :
*:                                                                     :
*:      Arguments: DQ2,DX: double prec. input xbj,q2                   :
*:                 DR*: double prec R output                           :
*:                                                                     :
*:      Called by: (evtl.fmf2pl) / rads1g / radw12                  :
*:                                                                     :
*:      Action:    calculate R                                         :
*:                  CREFERENCE:  L.W.WHITLOW, SLAC-REPORT-357          :
*:                  PH.D. THESIS, STANFORD UNIVERSITY,                 :
*:                  MARCH 1990.                                        :
*:                                                                     :
*:=====================================================================:
      IMPLICIT NONE
	include "mc_set.inc"

      DOUBLE PRECISION DQ2, DX
      DOUBLE PRECISION DR

      REAL R
      REAL QQ35, XX
      REAL FAC, RLOG, Q2THR
      REAL R_A, R_B, R_C
C
C Data-Definition of R-Calculation, see
C            L.W.WHITLOW, SLAC-REPORT-357,
C            PH.D. THESIS, STANFORD UNIVERSITY,
C            MARCH 1990.
      REAL AR1990(3), BR1990(3), CR1990(3)
      DATA AR1990  / .06723, .46714, 1.89794 /
      DATA BR1990  / .06347, .57468, -.35342 /
      DATA CR1990  / .05992, .50885, 2.10807 /

      XX=real(DX)
      IF (DQ2.LT.0.35) THEN
        QQ35=0.35
      ELSE
        QQ35=real(DQ2)
      ENDIF
C
C *** If Q2 < 0.35 then variable "R" is calculated at the fixed Q2 of 0.35
C
      FAC   = 1+12.*(QQ35/(1.+QQ35))*(.125**2/(XX**2+.125**2))
      RLOG  = FAC/LOG(QQ35/.04)
      Q2THR = 5.*(1.-XX)**5

      R_A   = AR1990(1)*RLOG +
     &        AR1990(2)/SQRT(SQRT(QQ35**4+AR1990(3)**4))
      R_B   = BR1990(1)*RLOG +
     &        BR1990(2)/QQ35 + BR1990(3)/(QQ35**2+.3**2)
      R_C   = CR1990(1)*RLOG +
     &        CR1990(2)/SQRT((QQ35-Q2THR)**2+CR1990(3)**2)
      R     = (R_A+R_B+R_C)/3.

      IF (DQ2.GE.0.35) THEN
        DR=dble(R)
      ELSE
        DR=dble(R)*DQ2/0.35
      ENDIF

c      print*,'R:',R
      
      END

