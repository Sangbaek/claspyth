C...As an example, consider a main program of the form
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
#include "phydata1234.inc"
      character*500  claspythfiles,fileinp

      INTEGER PYK,PYCHGE,PYCOMP
      REAL pyth_xsec,rccorr,sigobs,sigtrue
      DOUBLE PRECISION pbeam, px, py, pz

C...Input and output strings.
      CHARACTER*12 FRAME,BEAM,TARGET
      CHARACTER*100 PARAM
      LOGICAL UseLut,GenLut,lOK
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
      COMMON/PYINT1/MINT(400),VINT(400)
      Save /PYSUBS/,/PYPARS/,/PYINT1/



	include "mc_set.inc"
	include "py6strf.inc"
	include "mcRadCor.inc"
	include "radgen.inc"
c
#include "pyhtlist.inc"
#include "options.inc"
c
       data v4el0 /0.0,0.0,0.0,0.0/
       data v4pr0 /0.0,0.0,0.0,0.93827/
c      
c ... force block data modules to be read
        external pydata

        px = 0.
        py = 0.
        pz = 0.

       Model=0
       qedrad=0
       UseLUT=.false.
       GenLUT=.false.
       emass=PYMASS(11)
       rccorr =1.
c
       call CL_OPTIONS          ! command line processing
c
       v4el0(3)=cl_beam_energy
       v4el0(4)=cl_beam_energy
c
         if(cl_prom.gt.0) then ! redefine the proton value
            v4pr0(3)=-cl_prom
            v4pr0(4)=sqrt(v4pr0(3)**2+0.88035)
         endif
c
C...Read parameters for PYINIT call (frame, beam particle, target particle,
C...energy).
       FRAME='3MOM'
c
c      defining the 5 mom P(1,1-5) for electron and P(2,1-5) for nucleon
c
        do i_part_2_lund=1,3
         P(1,i_part_2_lund)=v4el0(i_part_2_lund)
         P(2,i_part_2_lund)=v4pr0(i_part_2_lund)
        enddo
c
       BEAM=cl_beam
       TARGET=cl_target
       PBEAM=cl_beam_energy
       PRINT *,'FRAME,BEAM,TARGET,PBEAM',FRAME,BEAM,TARGET,PBEAM
C       
C       
C       
C...Read parameters for PYTHIA-Model = which generation is done
      call getenv('CLASPYTHIA_DECLIST',claspythfiles)
      write(fileinp,'(A,A)') trim(claspythfiles),'/PYTHIA-input.dat'
       print *,'READING PYTHIA input values from',claspythfiles,fileinp
       OPEN(unit=61,file=fileinp,status='unknown')     
       READ(61,*) Model
C...Read parameters of radcorr: do radcorr swith also operated look-up table operation  
       READ(61,*) qedrad
C...Read min/max y of generation range      
       READ(61,*) mcSet_YMin,mcSet_YMax
C...Read min/max Q2 of generation range      
       READ(61,*) mcSet_Q2Min,mcSet_Q2Max
C...Read target type mass and charge
       READ(61,*) mcSet_TarA,mcSet_TarZ
C...Read information for cross section used in radgen
       READ(61,*) genSet_FStruct,genSet_R
C...Loop over reading and setting parameters/switches.
  100  READ(61,'(A)',END=200) PARAM
       CALL PYGIVE(PARAM)
       GOTO 100
       close(61)
C       
C       
C
C...Initialize PYTHIA.      
  200  write(*,*) 'NOW all parameters are read by PYTHIA'
C       call pylist(11)
c       call pylist(12)
       if (Model.eq.0) then
         UseLUT=.false.
         GenLUT=.false.
         qedrad=0
         MSTP(199)=0
       elseif (Model.eq.1) then
         if (qedrad.eq.1) then
           mcSet_EneBeam=real(pbeam)
           UseLUT=.true.
           MSTP(199)=1
           call radgen_init(UseLUT,GenLUT)
           write(*,*) 'I have initialized radgen'
         elseif (qedrad.eq.2) then
           write(*,*) 'radgen lookup table will be generated'
           mcSet_EneBeam=real(pbeam)
           UseLUT=.true.
           GenLUT=.true.
           MSTP(199)=1
           call radgen_init(UseLUT,GenLUT)
           goto 500
         elseif (qedrad.eq.0) then
           mcSet_EneBeam=real(pbeam)
           UseLUT=.false.
           GenLUT=.false.
           MSTP(199)=1
         endif
       endif 
      write(fileinp,'(A,A)') trim(claspythfiles),'/pythdeclist-curr.dat'
       print *,'Reading PYTHIA decay settings',claspythfiles,fileinp
       OPEN (unit=77,file=fileinp,status='UNKNOWN')
       if(cl_dumpdecayfile) then
        print *,'WRITING DECAY file PYUPDA', fileinp
        CALL PYUPDA(1,77)
       else
        print *,'READING DECAY file PYUPDA', fileinp
        CALL PYUPDA(2,77)
       endif 
       call PYLIST(12)
       close(77)
c       stop
c
       call pyinit (FRAME,beam,target, pbeam)
c
c     read the random seed
c
c      dir='/group/gpd/claspith/'
      if(.not.cl_docker.and.cl_rndm) then
         print *,'read random seed from','tmp/'//
     6cl_hispath(1:4)//cl_target(1:2)//'.txti'
         OPEN(unit=42,file='tmp/'//
     6cl_hispath(1:4)//cl_target(1:2)//'.txti',STATUS='OLD',
     6FORM='UNFORMATTED')
        call PYRSET(42,0)
        CLOSE(42)
      endif
c
c      call INIT_PAW            ! init hbook

C...Event generation loop
       DO 300 IEV=1,cl_triggers
         CALL PYEVNT
c         IF(IEV.LE.NPRT) CALL PYLIST(2)
C...Calculate the kinematics of the scattered lepton
        genpprim = sqrt(geneprim**2 - emass**2)
        px = genpprim*sin(genthe)*cos(genphi)
        py = genpprim*sin(genthe)*sin(genphi)
        pz = genpprim*cos(genthe)         

C Calculate all needed to get the info for cross sections and weights
        if(mod(iev,cl_nprint).eq.0) print *, 'event #',iev,nev,(MINT(III),III=1,10)

        IF ((MSTP(199).eq.1).and.(qedrad.ne.0)) then

C P.L. Born cross section at observed kinematics
         call MKF2(dble(genq2),dble(genx),
     +           mcSet_TarA,mcSet_TarZ,py6f2,py6f1)
         IF(genSet_FStruct(1:4).EQ.'ALLM') THEN
          call MKR(dble(genq2),dble(genx),py6r)
          py6f1=(1.+genq2/gennu**2)/(2.*genx*(1.+py6r))*py6f2
         ENDIF
         sigobs=pyth_xsec(genx,genq2,py6f1,py6f2)
     +          *genx*genq2

C P.L. Born cross section at true kinematics
         call MKF2(dble(mcRadCor_Q2True),dble(mcRadCor_XTrue),
     +           mcSet_TarA,mcSet_TarZ,py6f2,py6f1)
          IF(genSet_FStruct(1:4).EQ.'ALLM') THEN
          call MKR(dble(mcRadCor_Q2True),dble(mcRadCor_XTrue),py6r)
          py6f1=(1.+mcRadCor_Q2True/mcRadCor_NuTrue**2)/
     +     (2.*mcRadCor_XTrue*(1.+py6r))*py6f2
          ENDIF
         sigtrue=pyth_xsec(mcRadCor_XTrue, mcRadCor_Q2True,py6f1, py6f2)
     +          *mcRadCor_XTrue*mcRadCor_Q2True

C P.L. Ratio of sigma_born(observed)/sigma_born(true)
          IF(mcRadCor_EBrems.eq.0) THEN
           IF (sig1g.gt.0.) then
            rccorr=(tbor+mcRadCor_TailIne)/sig1g/
     +     (MINT(199)+1)
           ELSE
            rccorr=0.
           ENDIF
          ELSEIF(mcRadCor_EBrems.gt.0) THEN
           IF ((sig1g.gt.0.).and.(sigtrue.gt.0.)) then
             rccorr=(tbor+mcRadCor_TailIne)/sig1g*
     +       sigobs/sigtrue/(MINT(199)+1)
           ELSE
            rccorr=0.
           ENDIF
          ENDIF

        ELSE

         call MKF2(dble(genq2),dble(genx),
     +           mcSet_TarA,mcSet_TarZ,py6f2,py6f1)
         IF(genSet_FStruct(1:4).EQ.'ALLM') THEN
          call MKR(dble(genq2),dble(genx),py6r)
          py6f1=(1.+genq2/gennu**2)/(2.*genx*(1.+py6r))*py6f2
         ENDIF

        ENDIF
C       write(*,*)"F1 =   ",py6f1, "F2 =   ",py6f2, 
C    +            "rccorr = ", rccorr
         call fillevent(lOK)
c
         if(lOK) then
             IEVNTUPLE= IEVNTUPLE+1
c            print *,'---y=',1.0-gelee/cl_beam_energy
            radq2=genq2
            radnu=genx
            gelex=geneprim
            geley=genthe
            gelez=genphi
            gelewgy=mcRadCor_Sigcor
            grade=mcRadCor_EBrems
            gradt=mcRadCor_ThetaBrems
            gradf=mcRadCor_PhiBrems
            radweight=mcRadCor_Redfac
            r_sigrad=mcRadCor_SigRad
            r_sigcor=mcRadCor_Sigcor
            if(mcRadCor_cType.eq.'elas') then 
             iproctpe=2
            else if (mcRadCor_cType.eq.'qela') then
             iproctpe=3
            else if (mcRadCor_cType.eq.'inel') then
             iproctpe=1
            else
             iproctpe=0
            endif
c            
c            if(cl_verb.gt.100) call printntuple
c            call HFNT(22)
            call write2file(i,1)  ! 1-all 2-GSIM-lundc
        endif
  300  CONTINUE

c
c
      close(41)
        if(.not.cl_docker) then
        OPEN(unit=42,file='tmp/'//
     6cl_hispath(1:4)//cl_target(1:2)//'.txt',FORM='UNFORMATTED')
        call PYRGET(42,0)
        CLOSE(42)
	endif
C...Print the Pythia cross section which is needed to get an absolute 
C   normalisation the number is in microbarns
       write(*,*)"Pythia total cross section normalisation",
     +           pari(1)*1000.
      
C...Print cross sections.
       CALL PYSTAT(1)
       CALL PYSTAT(4)

  500  if (qedrad.eq.2) then
         write(*,*) 'lookup table is generated;'
         write(*,*) 'to run now pythia change parameter qedrad to 1'
       endif
   
       END

C***********************************************************************
          subroutine write2file(jj,il)
C***********************************************************************
c write the LUND event:
c
c                          Header is currently dummy---> to be fixed
c
c the header:
c N number of particles in the event
c PARL(1),PARL(2): number of nucleons and protons of target
c cl_pol,elpol: polarization of proton and electron
c LEX,LEYY,LEW2,LEQ2,LEU: x,y,W^2,Q^2, nu
c
c the body:
c i particle number
c plu(i,6) : particle charge
c
c K(I,1) = 3, 13} or \ttt{14}, instead, it  contains special
c K(I,2) : particle code
c K(I,3) : line number of parent particle, where known,otherwise 0. 
c K(I,4) : normally the line number of the first daughter and
c          0 for an undecayed particle or unfragmented parton.
c
c P(I,1) $p_x$, momentum in the $x$ direction, in GeV/$c$.
c P(I,2) $p_y$, momentum in the $y$ direction, in GeV/$c$.
c P(I,3) $p_z$, momentum in the $z$ direction, in GeV/$c$.
c P(I,4) $E$, energy, in GeV.
c P(I,5) $m$, mass, in GeV/$c^2$. 
c       In parton showers, with  space-like virtualities, 
c       i.e.\ where $Q^2 = - m^2 > 0$, one puts \ttt{P(I,5)}$ = -Q$.
c
c V(I,1) $x$ position of production vertex, in cm.
c V(I,2) $y$ position of production vertex, in cm.
c V(I,3) $z$ position of production vertex, in cm.
c V(I,4) time of production, in mm/$c$ ($\approx 3.33 \times 10^{-12}$ s).
c V(I,5) proper lifetime of particle, in mm/$c$
c      ($\approx 3.33 \times 10^{-12}$ s). 
c      If the particle is not expected to decay V(I,5)=0}.
c  A line with K(I,1)=4, i.e.\ a  particle that could have decayed,
c   but did not within the allowed region, has the proper non-zero V(I,5).
cIn the absence of electric or magnetic fields, or other
cdisturbances, the decay vertex \ttt{VP} of an unstable particle
cmay be calculated as VP(j) = V(I,j) + V(I,5)*P(I,j)/P(I,5)},
c
c
c
c
          implicit none
#include "pyhtlist.inc"
#include "options.inc"
#include "phydata1234.inc"
          integer ifile,itotal,iparl23
          DATA ifile /0/
          DATA itotal /0/
          character*200  fname
c
          real plu,random_num,xelex,yeley,zelez
          integer jj,il,i,j,istat,jc,m
          logical lFirst
          DATA lFirst /.TRUE./
          real elpol
          DATA elpol /1.0/
c
          integer  beam_pid,targ_id,process_id
          data beam_pid/11/,targ_id/2212/,process_id/1/
c
          if(itotal.ge.cl_fileev) then   ! each file correspond to 20k events
           close (41)
           lFirst=.TRUE.
           ifile=ifile+1
           itotal=0
          endif
c
          if (lFirst) then   ! open the file
           if(cl_docker) then
             write(fname,'(A)') 'claspyth.dat'
             cl_fileev=cl_triggers+100
           else
             call get_datfilename()

             jc = 0
             do m=1,len(datfilename)
                if ((datfilename(m:m).eq.' ').and.(jc.eq.0)) then
                   jc=m-1
                endif
             enddo
             if (jc.eq.0) jc=len(datfilename)

           if(ifile.le.9) then
             write(fname,'(A,i1,A)') datfilename(1:jc),ifile,'.dat'
             write(6,'(A,i1,A)') datfilename(1:jc),ifile,'.dat'
           else if(ifile.ge.10.and.ifile.le.99) then
            write(fname,'(A,i2,A)') datfilename(1:jc),ifile,'.dat'
            write(6,'(A,i2,A)') datfilename(1:jc),ifile,'.dat'
           else if(ifile.ge.100.and.ifile.le.999) then
            write(fname,'(A,i3,A)') datfilename(1:jc),ifile,'.dat'
            write(6,'(A,i3,A)') datfilename(1:jc),ifile,'.dat'
           else
            write(fname,'(A,i4,A)') datfilename(1:jc),ifile,'.dat'
            write(6,'(A,i4,A)') datfilename(1:jc),ifile,'.dat'
           endif 

           endif ! docker check
             print *, 'open file name ',fname
             OPEN(unit=41,file=fname,status='unknown')
             lFirst=.FALSE.
          endif
c
                itotal=itotal+1
c
               if(il.eq.2) then
                 j=0
                 do i=1,N
                  if(k(i,1).lt.11) j=j+1   ! save space
                 enddo
                 write(41,111) j,1.0,gQ2,gW,1.0*elpol,1.0*elpol
               else
                write(41,101) N,1,1,1.0*elpol,1.0*elpol,
     6     beam_pid,cl_beam_energy,targ_id,process_id,PARI(1)*1000000.0
     6,radweight
               endif
c
 111           FORMAT(2x,I5,10F10.4)
 110           FORMAT(2x,10I10)
 112           FORMAT(2x,10F10.4)
c 101      FORMAT(2x,I10,2F4.0,2I3,5F6.3)
 101      FORMAT(2x,I10,2I6,2F5.1,I5,F10.3,2I5,3E14.7)
 102      FORMAT(2x,I3,F4.0,I5,I7,2I3,5F10.4,2X,5F10.4)
          xelex=(random_num()-0.5)*cl_rast+cl_xpos
          yeley=(random_num()-0.5)*cl_rast+cl_ypos
          zelez=cl_zpos+(random_num()-0.5)*cl_zwidth
              do i=1,N
               if(il.eq.2) then
                 if(k(i,1).lt.11) then   ! save space
                   write(41,110) (k(i,j),j=1,5) 
                   write(41,112) (p(i,j),j=1,5) 
                   write(41,112) cl_xpos,cl_ypos,cl_zpos,0.0,0.0
                 endif 
               else
                 V(i,1)=xelex+V(i,1)*0.1
                 V(i,2)=yeley+V(i,2)*0.1
                 V(i,3)=zelez+V(i,3)*0.1
              write (41,102) i,plu(k(i,2)),(k(i,j),j=1,4),(P(i,j),j=1,5)
     6                                               ,(V(i,j),j=1,3)
               endif
              enddo
          return
          end
c
c
          real function plu(id)
         integer id
         SELECT CASE (id)  
         CASE (11) ! electron
          plu=-1.0
         CASE (2112) ! neutron
          plu=0.0
         CASE (2212) ! proton
          plu=1.0
         CASE (211) ! pi+
          plu=1.0
         CASE (321) ! K+
          plu=1.0
         CASE (-211) ! pi-
          plu=-1.0
         CASE (-321) ! K-
          plu=-1.0
         CASE (111) ! pi0
          plu=0.0
         CASE (22) ! pi0
          plu=0.0
         CASE (113) ! ro0
          plu=0.0
         CASE (213) ! ro+
          plu=1.0
         CASE (223) ! omega
          plu=0.0
         CASE DEFAULT
          plu=0.0
         END SELECT
          return
         end
c 
      subroutine printntuple()
      implicit none
#include "pyhtlist.inc"
c
       print *,'ELEGEN-1 ',gnele,gelebits,ghpar,gproh,geleh,gelee,gelet
       print *,'ELEGEN-2 ',gelef,gelex,geley,gelez,gelewgt,gxb,gq2,gw
c234567890123456789012345678901234567890123456789012345678901234567890123456
*
       print *,'HADGEN-1 ', gnpi0,gnpip,gnpim,gnpro,gngam
       print *,'HADGEN-2 ',gnneu,gnrop,gnro0,gnome
       print *,'HADGEN-3 ',nsproc,nsele,nspi0,nspip,nspim
       print *,'HADGEN-4 ',gpi0e,gpi0t,gpi0f
       print *,'HADGEN-5 ',gpipe,gpipt,gpipf
       print *,'HADGEN-6 ',gpime,gpimt,gpimf
       print *,'HADGEN-7 ',gproe,gprot,gprof
       print *,'HADGEN-8 ',gneue,gneut,gneuf
       print *,'HADGEN-9 ',grope,gropt,gropf
       print *,'HADGEN-10 ',gro0e,gro0t,gro0f
       print *,'HADGEN-11 ',gomee,gomet,gomef
       print *,'HADGEN-12 ',gpipm,gpimm,gpi0m,gro0m,gropm,gomem
       print *,'RADGN-1 ',iproctpe,radq2,radnu,radweight,radweight,r_sigrad,r_sigcor
       print *,'RADGN-2 ',grade,gradt,gradf
       return
       end
c
      subroutine fillevent(lOK)
      implicit none
#include "phydata1234.inc"
#include "pyhtlist.inc"
#include "options.inc"
      integer mint
      double precision VINT
      COMMON/PYINT1/MINT(400),VINT(400)
       real v4pr(4),v4el(4),v4el2(4),v4po(4),v4ne(4)
     6 ,v4pip(4),v4pim(4),v4pi0(4)
     6 ,v4rop(4),v4ro0(4),v4ome(4),tmp4(4),qiu4(4),vmass2,vdotm,tnorm(4)
       real PolarTheta,Azimuthalphi,anu,zmin,mismas
       integer i,part_id
       logical lOK,hOK
       call VZERO(gnele,14)
       call VZERO(gnpos,10)
       call VZERO(gnpi0,44)
       call VZERO(iproctpe,8)
c
       zmin=0.4
       hOK=.FALSE.
       lOK=.FALSE.
c
       do i=1,N
        if(k(i,2).eq.cl_pid) then
           call pylist(1)
           hOK=.TRUE.
        endif
        if(k(i,2).gt.9999.or.k(i,2).lt.-9999) k(i,2)=0 
c
        if(k(i,1).lt.21) then
         part_id=k(i,2)
c
         SELECT CASE (part_id)  
         CASE (11) ! electron
          gnele=gnele+1
          if(gnele.eq.1) then
          v4el(1)=P(i,1)
          v4el(2)=P(i,2)
          v4el(3)=P(i,3)
          v4el(4)=P(i,4)
          gelee=v4el(4)
          gelet=PolarTheta(v4el(1),v4el(2),v4el(3))
          gelef=Azimuthalphi(v4el(1),v4el(2))
          nsele=k(k(i,3),2)
          else
          v4el2(1)=P(i,1)
          v4el2(2)=P(i,2)
          v4el2(3)=P(i,3)
          v4el2(4)=P(i,4)
          gele2=v4el2(4)
          gelt2=PolarTheta(v4el2(1),v4el2(2),v4el2(3))
          gelf2=Azimuthalphi(v4el2(1),v4el2(2))
          nsel2=k(k(i,3),2)
         endif
c
         CASE (-11) ! positron
          gnpos=gnpos+1
          v4po(1)=P(i,1)
          v4po(2)=P(i,2)
          v4po(3)=P(i,3)
          v4po(4)=P(i,4)
          gpose=v4po(4)
          gpost=PolarTheta(v4po(1),v4po(2),v4po(3))
          gposf=Azimuthalphi(v4po(1),v4po(2))
          nspos=k(k(i,3),2)
         CASE (2112) ! neutron
          gnneu=gnneu+1
          v4ne(1)=P(i,1)
          v4ne(2)=P(i,2)
          v4ne(3)=P(i,3)
          v4ne(4)=P(i,4)
          gneue=v4ne(4)
          gneut=PolarTheta(v4ne(1),v4ne(2),v4ne(3))
          gneuf=Azimuthalphi(v4ne(1),v4ne(2))
        CASE (2212) ! proton
          gnpro=gnpro+1
          v4pr(1)=P(i,1)
          v4pr(2)=P(i,2)
          v4pr(3)=P(i,3)
          v4pr(4)=P(i,4)
          gproe=v4pr(4)
          gprot=PolarTheta(v4pr(1),v4pr(2),v4pr(3))
          gprof=Azimuthalphi(v4pr(1),v4pr(2))
         CASE (211) ! pi+
          gnpip=gnpip+1
           if(gnpip.eq.1.or.v4pip(4).lt.P(i,4)) then
           v4pip(1)=P(i,1)
           v4pip(2)=P(i,2)
           v4pip(3)=P(i,3)
           v4pip(4)=P(i,4)
           gpipe=v4pip(4)
           gpipt=PolarTheta(v4pip(1),v4pip(2),v4pip(3))
           gpipf=Azimuthalphi(v4pip(1),v4pip(2))
           nspip=k(k(i,3),2)
           endif
         CASE (-211) ! pi-
          gnpim=gnpim+1
          if(gnpim.eq.1.or.v4pim(4).lt.P(i,4)) then
          v4pim(1)=P(i,1)
          v4pim(2)=P(i,2)
          v4pim(3)=P(i,3)
          v4pim(4)=P(i,4)
          gpime=v4pim(4)
          gpimt=PolarTheta(v4pim(1),v4pim(2),v4pim(3))
          gpimf=Azimuthalphi(v4pim(1),v4pim(2))
          nspim=k(k(i,3),2)
          endif
         CASE (111) ! pi0
          gnpi0=gnpi0+1
          if(gnpi0.eq.1.or.v4pi0(4).lt.P(i,4)) then
          v4pi0(1)=P(i,1)
          v4pi0(2)=P(i,2)
          v4pi0(3)=P(i,3)
          v4pi0(4)=P(i,4)
          gpi0e=v4pi0(4)
          gpi0t=PolarTheta(v4pi0(1),v4pi0(2),v4pi0(3))
          gpi0f=Azimuthalphi(v4pi0(1),v4pi0(2))
          nspi0=k(k(i,3),2)
          endif
         CASE (113) ! ro0
          gnro0=gnro0+1
          v4ro0(1)=P(i,1)
          v4ro0(2)=P(i,2)
          v4ro0(3)=P(i,3)
          v4ro0(4)=P(i,4)
          gro0e=v4ro0(4)
          gro0t=PolarTheta(v4ro0(1),v4ro0(2),v4ro0(3))
          gro0f=Azimuthalphi(v4ro0(1),v4ro0(2))
         CASE (213) ! ro+
          gnrop=gnrop+1
          v4rop(1)=P(i,1)
          v4rop(2)=P(i,2)
          v4rop(3)=P(i,3)
          v4rop(4)=P(i,4)
          grope=v4rop(4)
          gropt=PolarTheta(v4rop(1),v4rop(2),v4rop(3))
          gropf=Azimuthalphi(v4rop(1),v4rop(2))
         CASE (223) ! omega
          gnome=gnome+1
          v4ome(1)=P(i,1)
          v4ome(2)=P(i,2)
          v4ome(3)=P(i,3)
          v4ome(4)=P(i,4)
          gomee=v4ome(4)
          gomet=PolarTheta(v4ome(1),v4ome(2),v4ome(3))
          gomef=Azimuthalphi(v4ome(1),v4ome(2))
         CASE DEFAULT
         END SELECT
        endif       
       enddo

          geleh=MSTI(1)
          ghpar=MINT(2)
          gproh=MSTI(9)
c
         if (gnele.ge.1.and.(gelet*180/3.1415).gt.0.000001         ! take everything
     6        .and.(1.0-gelee/cl_beam_energy).lt.0.999999) then
          if(cl_beam.eq.'gamma') then ! photon beam
            qiu4(4)=cl_beam_energy
            qiu4(3)=cl_beam_energy
            qiu4(1)=0
            qiu4(2)=0
           else
            call vdifm(v4el0,v4el,qiu4,4)
            gq2=vdotm(qiu4,qiu4,4)
c
             if(qiu4(4).gt.0) then
              gxb=gq2/2.0/0.93827/qiu4(4)
c
              if(gxb.gt.0) then
               gw=0.88035+gq2/gxb-gq2
              else
               gw=0
              endif
              anu=qiu4(4)
             else
              gxb=0
              gw=0
              anu=100000.0
             endif
           endif
c          
c        
c
         if(gelet*180/3.1415.gt.0.000001) then      ! write everything
c
c          fill relevant mismases
c
          if(gnpos.gt.0.and.gele2.gt.0) then
           call vsumm(v4el2,v4po,tnorm,4)
           gpmas=vmass2(tnorm)
          endif
c
            if(gnpip.gt.0) gpipm=mismas(qiu4,v4pip)
            if(gnpim.gt.0) gpimm=mismas(qiu4,v4pim)
            if(gnpi0.gt.0) gpi0m=mismas(qiu4,v4pi0)
            if(gnro0.gt.0) gro0m=mismas(qiu4,v4ro0)
            if(gnrop.gt.0) gropm=mismas(qiu4,v4rop)
            if(gnome.gt.0) gomem=mismas(qiu4,v4rop)
c
              if(cl_pid.eq.13) then
               if(gpost*57.3.gt.cl_tmin.and.gelet*57.3.gt.cl_tmin)
c (gpost*57.3.gt.cl_tmin.and.gelt2*57.3.gt.cl_tmin).or.
c     6            (gpost*57.3.gt.cl_tmin.and.gelet*57.3.gt.cl_tmin))  ! e+e- in the forward detector
     6            lOK=.TRUE. 
              else if(cl_pid.eq.11) then 
               if(gpost*57.3.gt.cl_tmin.or.gelet*57.3.gt.cl_tmin.or.
     6            gelt2*57.3.gt.cl_tmin) then ! e+ or e- in the forward detector
                  lOK=.TRUE. 
c              print *,'angles',gelet*57.3,gpost*57.3,gelt2*57.3
                 endif
              else if(cl_pid.eq.-211) then 
               if(gelet*57.3.gt.cl_tmin.or.gelt2*57.3.gt.cl_tmin.or.
     6            gpimt*57.3.gt.cl_tmin) ! pi- or e- in the forward detector
     6            lOK=.TRUE. 
              else if(cl_pid.eq.443) then  ! J/Psi
	        if(hOK)lOK=.TRUE.   
              else if(cl_pid.eq.213) then  ! only rho+
	        if(hOK)lOK=.TRUE.   
              else
               lOK=.TRUE.  ! write everything
              endif
c            endif
cccc            lOK=.TRUE.
              else if(cl_pid.eq.-211) then 
          endif
       endif
c          
      return
      end
c
       real function vmass(a)
       real a(4)
       real vm
        vm= vdotm(a,a,4)
        if (vm.lt.0.0) then
          vmass=sqrt(-vm)
        else
          vmass=-1.0
        endif 
       return
       end
c
       real function vdotm(a,b,n)
       real a(n),b(n),s
       integer i,n
       s=0.0
       do i=1,3
         s=s+a(i)*b(i)
       enddo
       if(n.eq.4) s=s-a(n)*b(n)
       vdotm=s
       return
       end
c   
       subroutine vsumm(a,b,c,n)
       real a(n),b(n),c(n)
       integer i,n
        do i=1,n
          c(i)=a(i)+b(i)
        enddo
       return
       end
c
       subroutine vdifm(a,b,c,n)
       real a(n),b(n),c(n)
       integer i,n
        do i=1,n
          c(i)=a(i)-b(i)
        enddo
       return
       end

c
        real function  PolarTheta(vx,vy,vz)
      implicit none
        real vx,vy,vz,pmod,theta
        pmod=vx*vx+vy*vy+vz*vz
        if(pmod .gt. 0 ) then
         theta=acos(vz/sqrt(pmod))
        else
         theta=-100
        endif
        PolarTheta=theta
        return
        end

        real function Azimuthalphi(vx,vy )
      implicit none
        real pi
        parameter ( PI=3.1415926)
        real vx,vy,pmod,phi,cosf
        pmod=vx*vx+vy*vy
        if(pmod .gt. 0 ) then
         pmod=sqrt(pmod)
         cosf=vx/pmod
        else
         cosf=1.0
        endif
        if(abs(cosf) .le. 1.0) phi=acos(cosf);
        if(vy .lt. 0.0) phi= 2*PI-phi;
        Azimuthalphi=phi
        return
        end
c
       SUBROUTINE CL_OPTIONS
       implicit none
#include "options.inc"
       integer i,iran,irannom,ivalnum
       real VEC1,valnum
       character*20 cnumber
       integer  numopts
       integer*4 now(3)
c
c    command line processing
c
c     defaults
      cl_rndm = .FALSE. 
      cl_xpos = 0.0  
      cl_ypos = 0.0  
      cl_zpos = 0.0  
      cl_zwidth = 5.0
      cl_rast=0.01  
      cl_triggers = 100  
      cl_seed = 0  
      cl_beam_energy=10.6   ! def 
      cl_prom=0.0           ! def 0-clas12
      cl_beam='gamma/e-'     ! electron 
      cl_hispath='pythi'    ! (w9)
      cl_target='p+'     ! (def proton target)
      cl_nprint=1000         ! print every cl_nprint event
      cl_fileev=200         ! change file every cl_fileev
      cl_verb=0
      cl_radgen=.FALSE.
      cl_dis=.FALSE.
      cl_pid=0
      cl_tmin=5.0
      cl_emin=1.0
      cl_dumpdecayfile=.FALSE.
      cl_docker=.FALSE.
c
      numopts  = IARGC()
c      
      if(numopts.gt.0 ) then
        i=0
      else
        goto 2
      endif
 1    i=i+1
        CALL GETARG(i,cnumber)
        if(cnumber.eq.'--trig'.and.i.lt.numopts) then
          i=i+1
          CALL GETARG(i,cnumber)
          cl_triggers=valnum(cnumber)
        elseif(cnumber.eq.'--seed'.and.i.lt.numopts) then
          i=i+1
          CALL GETARG(i,cnumber)
          cl_seed=ivalnum(cnumber)
        elseif(cnumber.eq.'--nmax'.and.i.lt.numopts) then
          i=i+1
          CALL GETARG(i,cnumber)
          cl_fileev=valnum(cnumber)
        elseif(cnumber.eq.'--v'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_verb=valnum(cnumber)
        elseif(cnumber.eq.'--print'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_nprint=valnum(cnumber)
        elseif(cnumber.eq.'--path'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cl_hispath)
        elseif(cnumber.eq.'--targ'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cl_target)
        elseif(cnumber.eq.'--beam'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cl_beam)
        elseif(cnumber.eq.'--tmin'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_tmin=valnum(cnumber)
        elseif(cnumber.eq.'--zwidth'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_zwidth=valnum(cnumber)
        elseif(cnumber.eq.'--rast'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_rast=valnum(cnumber)
        elseif(cnumber.eq.'--radgen') then
           cl_radgen=.TRUE.
        elseif(cnumber.eq.'--dumpdecayfile') then
           cl_dumpdecayfile=.TRUE.
        elseif(cnumber.eq.'--dis') then
           cl_dis=.TRUE.
        elseif(cnumber.eq.'--ebeam'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_beam_energy=valnum(cnumber)           
        elseif(cnumber.eq.'--prom'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_prom=valnum(cnumber)           
        elseif(cnumber.eq.'--pid'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_pid=valnum(cnumber)           
        elseif(cnumber.eq.'--xpos'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_xpos=valnum(cnumber)           
        elseif(cnumber.eq.'--ypos'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_ypos=valnum(cnumber)           
        elseif(cnumber.eq.'--zpos'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_zpos=valnum(cnumber)           
        elseif(cnumber.eq.'--rndm') then
           cl_rndm=.TRUE.
        elseif(cnumber.eq.'--docker') then
           cl_docker=.TRUE.
        else
          call printusage
          stop
        endif
        if(i.lt. numopts) goto 1
 2      continue
        call printvalues
c	
        call RNDMSEED()
c
        return 
        end
c
        subroutine printusage
        print *,' pythia [options]            '
        print *,'      option  value      default    comment  '
        print *,'  --trig nevents         10      number of triggers '
        print *,'  --nmax nevents     200      number of ev in 1 file '
        print *,'  --beam beam_type    gamma/e-     e- beam '
        print *,'     --targ target         p+       n0 for neutron '
        print *,'  --ebeam beam_energy     10.6  e- momentum in GeV '
        print *,'  --prom  proton mom   0-clas12  p+ momentum in GeV '
        print *,'  --xpos  x_position      0.0   beam x pos. in cm '
        print *,'  --ypos  y_position      0.0   beam y pos. in cm '
        print *,'  --zpos  z_position      0.0   target z pos. in cm '
        print *,'  --zwidth  5.0     target length in cm'
        print *,'  --rast    0.01    beam diameter in cm'
        print *,'  --radgen                   include radgen'
        print *,'  --path hisid          w9/p       write dir  '
        print *,'  --print nprint       1000   print nprint event '
        print *,'  --dis                .FALSE.     all dis events in '
        print *,'  --dumpdecayfile  .FALSE.  write the decay file '
        print *,'  --pid      0     write all events '
        print *,'            13      with e- & e+ in forward det.'
        print *,'            11   write e+ or e- in the forward det. '
        print *,'           -211  write pi- and e- in the forward det.'
        print *,'           443  write J/Psi'
        print *,'  --tmin    5      minimum angle in degree for e-,e+ '
        print *,'  --v     verblev       0     verbos. level '
        print *,'     --docker  .FALSE. gen a single file claspyth.dat'
        return
        end

        subroutine printvalues
        implicit none
#include "options.inc"
        print *,'pythia input values:'    
        print *,'triggers:  ',cl_triggers     
        print *,'BEAM type:         ',cl_beam  
        print *,'TARGET:              ', cl_target 
        print *,'BEAM ENERGY:         ',cl_beam_energy  
        print *,'path              ', cl_hispath 
        print *,'nprint            ',cl_nprint
        print *,'verbosity lev.    ',cl_verb
        print *,'RADGEN ON        :   ', cl_radgen 
        return
        end
c
      real function valnum(str)
c     return the real value contained into a string

      implicit integer (k-l)
      character *(*) str
      logical segno
      

      segno=.false.
      valnum=0.00
      lu=length(str)
      
c     check the number sign
      if(str(1:1).eq.'-')then
      segno=.true.
      str=str(2:lu)
      lu=lu-1
      end if

c     check if number is float or integer
      if(index(str,'.').ne.0)then
      iin=index(str,'.')-1
      else
      iin=lu
      end if

      ifr=lu-(iin+1)

c     translate the integer portion
      doi=1,iin  
      k=ichar(str(i:i))-48
      valnum=valnum+float(k)*10.00**float(iin-i)
      end do         

      if(iin.eq.lu)goto 10
      str=str(iin+2:lu)

c     translate the decimal portion
      doi=1,ifr  
      k=ichar(str(i:i))-48
      valnum=valnum+float(k)/10.00**float(i)
      end do

10    if(segno)valnum=-valnum

      return
      end     
c
      integer function ivalnum(str)
c     return the real value contained into a string

      implicit integer (k-l)
      character *(*) str
      logical segno
      double precision valnum
      

      segno=.false.
      valnum=0.00
      lu=length(str)
      
c     check the number sign
      if(str(1:1).eq.'-')then
      segno=.true.
      str=str(2:lu)
      lu=lu-1
      end if

c     check if number is float or integer
      if(index(str,'.').ne.0)then
      iin=index(str,'.')-1
      else
      iin=lu
      end if

      ifr=lu-(iin+1)

c     translate the integer portion
      doi=1,iin  
      k=ichar(str(i:i))-48
      valnum=valnum+float(k)*10.00**float(iin-i)
      end do         

      if(iin.eq.lu)goto 10
      str=str(iin+2:lu)

c     translate the decimal portion
      doi=1,ifr  
      k=ichar(str(i:i))-48
      valnum=valnum+float(k)/10.00**float(i)
      end do

10    if(segno)valnum=-valnum
       ivalnum=int(valnum)
      return
      end     
c
      integer function length(str)
c     return the string length without the blanks characters

      implicit integer (k-l)
      character *(*) str    
      
      lmax=len(str)      
       
c     search the last non blank character
      do i=lmax,1,-1
      if(str(i:i).ne.' ')then
      length=i
      return
      end if
      end do

      length=lmax

      return
      end     
c
       real function mismas(qiu4,v4)
       implicit none
       real qiu4(4),v4(4),pro4(4),tnorm(4),tnorm2(4)
       real Mp,vmass,vmass2,vangle,vdotm,phigstar
       parameter (Mp=0.93827)
c proton
c
        pro4(4)=Mp
        pro4(3)=0
        pro4(2)=0
        pro4(1)=0
c
         call vsumm(qiu4,pro4,tnorm,4)
         call vdifm(tnorm,v4,tnorm2,4)
         mismas=vmass2(tnorm2)
         return
         end
c
       real function vmass2(a)
       implicit none
       real vm,a(4), vdotm
        vm= vdotm(a,a,4)
        vmass2=-vm
       return
       end
c
c
      subroutine get_datfilename()
      implicit none
#include "options.inc"
#include "phydata1234.inc"
      INTEGER MSEL,MSELPD,MSUB,KFIN
      DOUBLE PRECISION CKIN
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
      character*4  cpol
      character*5  cdis
      character*12  cxsec
      character*21  ccuts
      real cl_pol
      data cl_pol/0.0/
c
c
      integer jc,m,ncdis
      write (6,*),'CKIN',(ckin(m),m=69,76)
ccc
      print *,'make the data file name',cl_beam_energy,cl_emin,cl_tmin 
      print *,'************************'
c
      if(cl_beam_energy.lt.10.0) then
      write(ccuts,'(A,F5.3,A,F3.2,A,F3.1)') 'e',cl_beam_energy,
     6'.emin',cl_emin,'tmin',cl_tmin 
      else
      write(ccuts,'(A,F6.3,A,F4.2,A,F3.1)') 'e',cl_beam_energy,
     6'.emn',cl_emin,'tmn',cl_tmin 
      endif
c
      print *,'get the x-section for ',PARI(1),ccuts
c
      if(cl_beam_energy.lt.15.0) then
      write(cxsec,'(A,F7.2,A)') '.xs',PARI(1)*1000000.0,'nb'
      else
      write(cxsec,'(A,F7.1,A)') '.xs',PARI(1)*1000000.0,'nb'
      endif
      print *,'**',ccuts,'++',cxsec
c
      if(cl_pol.eq.1) then
       cpol='.p1.'
      else if (cl_pol.eq.-1) then
       cpol='.m1.'
      else
       cpol='.00.'
      endif
c
      if (cl_pid.eq.0) then 
       cdis='.dis'
       ncdis=4
      else
         if (cl_pid.lt.1000) then
            write (cdis,'(A,I3,A)') '.',cl_pid
            ncdis=4
         else
            write (cdis,'(A,I4,A)') '.',cl_pid
            ncdis=5
         endif
c       write (cdis,'(A,F3.1,I4)') '.',cl_zmin,cl_pid
      endif



      print *,'PARI(1)=',PARI(1)
      write(datfilename,'(A)') 'tmp/claspyt'//
     6cl_target(1:2)//cpol//ccuts//cxsec//cdis(1:ncdis)//'.'
      print *,'filename ',datfilename
      end
C**********************************************************************
C      random number generator
C**********************************************************************
         real function random_num()
         implicit none
         random_num=RAND(0)
         return
         end
c
c
       subroutine RNDMSEED()
       implicit none
#include "options.inc"
       integer*4 now(3),ifilenum
       integer iran,irannom
       real rand,rnd
c
c
          if(cl_seed.eq.0) then
          call itime(now)     ! now(1)=hour, (2)=minute, (3)=second
	  irannom=now(3)+now(1)*3600+now(2)*60
	  else 
	  irannom=cl_seed
	  endif
           CALL SRAND(irannom)
          print *,'RandomSeedActuallyUsed:',irannom,now
c           DO 100 iran=1,irannom
c 100       RND=rand()
c
	   return
	   end

