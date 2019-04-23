
      subroutine radgen_init(bUseLUT,bGenLUT)

      implicit none

	include "mc_set.inc"
	include "cmpcom.inc"
	include "tailcom.inc"
	include "radgen.inc"
	include "radgenkeys.inc"

      real ebeam,polbeam,poltarget
      logical bUseLUT,bGenLUT
      INTEGER hermesfile
      CHARACTER*256 tname,fname

! ... force block data modules to be read
      external radata

*     kill_elas_res =0, all rad corrections
*     kill_elas_res =1, ???
*     kill_elas_res =2, no rad corrections
      kill_elas_res=0

*     initialise the rad. correction part !
*     =  2 : fint cernlib extrapolation routine
*     =  1 : 2dim spline
*     =  0 : make lookuptable
*     = -1 : do not lookup table , calculate all events
      if( bUseLUT ) then
        if( bGenLUT ) then
          ixytest = 0
        else
          ixytest = 2
        endif
      else
        ixytest = -1
      endif

*----------------------------
*     ire is target 1,2,3
*     now also 14 (nitrogen) and 84 (krypton)
      if(     mcSet_TarA .eq. 1 .and. mcSet_TarZ .eq. 1 ) then
        ire = 1
      elseif( mcSet_TarA .eq. 2 .and. mcSet_TarZ .eq. 1 ) then
        ire = 2
      elseif( mcSet_TarA .eq. 3 .and. mcSet_TarZ .eq. 2 ) then
        ire = 3
      elseif( mcSet_TarA .eq. 4 .and. mcSet_TarZ .eq. 2 ) then
        ire = 4
      elseif( mcSet_TarA .eq. 14 .and. mcSet_TarZ .eq. 7 ) then
        ire = 14
      elseif( mcSet_TarA .eq. 20 .and. mcSet_TarZ .eq. 10 ) then
        ire = 20
      elseif( mcSet_TarA .eq. 84 .and. mcSet_TarZ .eq. 36 ) then
        ire = 84
      else
        write(*,*)( 'RADGEN_INIT: invalid target selection' )
      endif

*----------------------------
*     plrun : beam polarisation    [-1.,+1.]
*     pnrun : nucleon polarisation     [-1.,0,+1.]

      plrun = 0.
      pnrun = 0.
      if(ire.lt.10) then
            tname='radgen/xytab0unp.dat'
            print *,'radgen_init: open ',ire,tname
      else  
            tname='radgen/xytab00unp.dat'
            print *,'radgen_init: open ',ire,tname
      endif

      if (ire.lt.10) then
          write(tname(13:13),'(i1)')ire
      else
          write(tname(13:14),'(i2)')ire
      endif

*----------------------------
* grid of important regions in theta (7*ntk)
      ntk = 35
*----------------------------
* photonic energy grid
      nrr = 100
*----------------------------
* min. energy in the calo (resolution parameter)
* as the Hermes calorimeter can only see single photons from a minimum
* energy of 0.5 GeV should the parameter not be changed from 0.1 
* to 0.5 GeV (E.C.A)
      demin=0.10

*----------------------------
      ap=2.*amp
      amp2=amp**2
      ap2=2.*amp**2
      if(kill_elas_res.eq.1)amc2=4d0

      if(ire.eq.1)then
        amt=.938272d0
        rtara=1d0
        rtarz=1d0
        fermom=0d0
      elseif(ire.eq.2)then
        amt=1.87561d0
        rtara=2d0
        rtarz=1d0
        fermom=.07d0
      elseif(ire.eq.3)then
        amt=2.80923d0
        rtara=3d0
        rtarz=2d0
        fermom=.164d0
      elseif(ire.eq.4)then
        amt=3.72742d0
        rtara=4d0
        rtarz=2d0
        fermom=.164d0
      elseif(ire.eq.14)then
*abr    amt=14.0067d0     ! wrong value (03.01.02)
        amt=13.0438d0
        rtara=14d0
        rtarz=7d0
        fermom=.221d0
        call fordop
      elseif(ire.eq.20)then
        amt=18.6228d0
        rtara=20d0
        rtarz=10d0
        fermom=.225d0  
        call fordop
      elseif(ire.eq.84)then
*abr    amt=83.80d0       ! wrong value (03.01.02)
        amt=78.1631d0
        rtara=84d0
        rtarz=36d0
        fermom=.260d0
        call fordop
      endif

*-----------------------------
* initialize lookup table
      if(ixytest.ge.0) then
        call xytabl(tname,mcSet_EneBeam,plrun,pnrun,ixytest,ire)
      endif

      end
