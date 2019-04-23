*CMZ :  6.00/04 01/12/94  22.35.52  by  hermes@kayak
*CMZ :  6.00/03 01/12/94  14.48.28  by  hermes@kayak
*CMZ :  6.00/00 14/09/94  18.52.09  by  phg
*CMZ :  6.00/00 13/05/94  11.57.28  by  wwc
*CMZ :  5.01/11 13/05/94  11.57.28  by  wwc
*CMZ :  5.00/01 24/01/93  13.23.44  by  Michael Dueren
*CMZ :  5.00/00 22/01/93  17.29.08  by  Michael Dueren
*-- Author :    Michael Dueren   22/01/93
*-- Author2:    Philipp Geiger   Okt.94
* modified and updated for N14 - 14.10.98 abr
* added Kr84 - 15.10.01 abr
* added Ne20 - 15.07.02 abr

      SUBROUTINE MKF2(DQ2,DX,A,Z,DF2,DF1)
*:=====================================================================:
*:                                                                     :
*:      Author:    Philipp Geiger  Last update: 12.02.1995             :
*:                                 Tested: yes                         :
*:                 24.09.98: update of 3He ratio, N14 ratio (abr)      :
*:                                                                     :
*:      Arguments: DQ2,DX: double prec. input xbj,q2                   :
*:                 DF2*,DR*: double prec F  output                     :
*:                                                                     :
*:      Called by: FMDISX                                              :
*:                                                                     :
*:      Action:    calculate f2 structure function, R                  :
*:                 calls different parametrizations of F2              :
*:                   now implemented:                                  :
*:                 TarA  FStruct  Action                               :
*:                 1     'P8'     calls Fctn. DF2H8                    :
*:                 2     'P8'     calls Fctn. DF2D8                    :
*:                 3     'P8'     calls Fctn. DF2H8,DF2D8  *)          :
*:                 1,2   'P15'    calls Subrout. MKF2P15               :
*:                 3     'P15'    calls Subrout. MKF2P15   *)          :
*:                 3     'P8H3'   calls Fctn. DF2D8       **)          :
*:                 3     'P8H4'   calls Fctn. DF2D8      ***)          :
*:                 14    'P8C '   calls Fctn. DF2D8      ***)          :
*:                 14    'P8N '   calls Fctn. DF2D8      ***)          :
*:                 3     '15H3'   calls Subrout. MKF2P15  **)          :
*:                 3     '15H4'   calls Subrout. MKF2P15 ***)          :
*:                 14    'P15C'   calls Subrout. MKF2P15 ***)          :
*:                 14    'P15N'   calls Subrout. MKF2P15 ***)          :
*:                                                                     :
*:                       *): Calc. F2(3He)=1/3(F2p+2*F2d)              :
*:                      **): Calc. F2(3He)=F2d*(F2(3He/F2d))           :
*:                     ***): Calc. F2(3He)=F2d*(F2(4He/F2d))           :
*:=====================================================================:
      IMPLICIT NONE

	include "py6strf.inc"
	include "mc_set.inc"

*-----DEKLARATION----------------------
      DOUBLE PRECISION DX,DQ2, DF2, DF1, df2temp
      DOUBLE PRECISION DF2H, DF2D
      INTEGER A, Z, iflavour
      real rdx,rdq2
      REAL F2ALLM
      EXTERNAL F2ALLM
      DOUBLE PRECISION DNP(5), DF2NF2P
      DOUBLE PRECISION gamma2, dnu, w2, pmass2, pmass, ebeam
      real xpq(-6:6), xdpq(-6:6)
! ... charge of quark flavours
      real qflavour(6)
      save qflavour
      data qflavour/1.,2.,1.,2.,1.,2./

C parameters for ratio F2(n)/F2(p)
C measured at NMC, (NMC Amaudruz et al. CERN-PPE/91-167))
      data DNP   / 0.976D0,    -1.34D0,      1.319D0, 
     &            -2.133D0,     1.533D0/
      DATA PMASS/0.93827/

*-----DEFINITIONS-----------------------

*-------------------------"normal  part"------------------
! If generating from LEPTO/PEPSI parton distribution functions, 
! can deal with all targets on an equal footing

! Fits, fits, and more fits

! ... hydrogen:
      IF(A.EQ.1.and.z.eq.1)THEN
         IF(genSet_FStruct(1:4).EQ.'ALLM') THEN
            rdx=real(dx)
            rdq2=real(dq2)
            DF2=F2ALLM(rdx,rdq2)      
            pmass2=pmass**2
            w2=pmass2+dq2*(1-dx)/dx
            ebeam=mcSet_EneBeam
            dnu=(w2-pmass2+dq2)/(2.*pmass)
            gamma2=dq2/(dnu**2)
            DF1=(1.D0+gamma2)/(2.D0*dx)/(1.D0+py6R)*DF2
         ELSEIF(genSet_FStruct(1:4).EQ.'F2PY') THEN
            call F2PYTH(dx,dq2,df1,df2)
         ELSE
*:           error:
            write(*,*)('invalid parametrisation choice in mkf2')
         ENDIF
*
! ... deuterium:
      ELSEIF(A.EQ.2.and.z.eq.1)THEN
*:        ALLM:
* abr 17.10.01
         IF(genSet_FStruct(1:4).EQ.'ALLM') THEN
            rdx=real(dx)
            rdq2=real(dq2)
            DF1=0.D0
            DF2=F2ALLM(rdx,rdq2)      
            DF2NF2P=DNP(1)+dx*(DNP(2)+dx*(DNP(3)+dx*(DNP(4)+dx*DNP(5))))
            DF2=DF2*0.5*(df2nf2p+1.)
         ELSE
*:           error
            write(*,*)('MKF2: invalid parametrisation choice FStruct')
         ENDIF
*
! ... neutron = 2*(deuterium_per_nucleon) - proton:
      ELSEIF(A.EQ.1.and.z.eq.0)THEN
*:        ALLM:
* abr 17.10.01
         IF(genSet_FStruct(1:4).EQ.'ALLM') THEN
            rdx=real(dx)
            rdq2=real(dq2)
            DF1=0.D0
            DF2=F2ALLM(rdx,rdq2)      
            DF2NF2P=DNP(1)+dx*(DNP(2)+dx*(DNP(3)+dx*(DNP(4)+dx*DNP(5))))
            DF2=DF2*df2nf2p
         ELSE
*:           error
            write(*,*)('MKF2: invalid parametrisation choice FStruct')
         ENDIF
*
      ELSE
*:        error:
         write(*,*)('MKF2: invalid target type')
      ENDIF
      END
