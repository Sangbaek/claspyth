      
      SUBROUTINE MKR(DQ2,DX,DR)

*
* select parametrization for R=sigma_L/sigma_T
*
      IMPLICIT NONE
      include "mc_set.inc"
      include "py6strf.inc"

      DOUBLE PRECISION DQ2, DX, DR
      DOUBLE PRECISION DELTAR

      IF ( genSet_R .EQ. '1990' ) THEN
* Whitlow et al.,  Phys.Lett.B 250(1990),193
         CALL R1990(DQ2,DX,DR)
         py6R=DR
      ELSE IF ( genSet_R .EQ. '1998' ) THEN
* E143, hep-ex/9808028
         CALL R1998(DQ2,DX,DR,DELTAR)
         py6R=DR
      ELSE IF ( genSet_R .eq. '0' ) THEN
* pure transverse (sigma_L=0)
         DR = 0.d0
         py6R=0.d0
      ELSE
         write(*,*)( 'MKR: invalid choice for R parametrization' )
      ENDIF

      RETURN
      END
