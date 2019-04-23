C.....Total inclusive Pythia cross section according to the model in  
C.....C.Friberg, T.Sjoestrand, J. High Energy Phys. JHEP 0009, 010 (2000)
C.....(neglecting GVMD cross sections)
C.....P.Liebing, 08/16/2003

      SUBROUTINE F2PYTH(x,q2,f1,f2)

      IMPLICIT NONE
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      INTEGER MSTP,MSTI
      DOUBLE PRECISION PARP,PARI
      SAVE/PYPARS/
      COMMON/PYINT1/MINT(400),VINT(400)
      INTEGER MINT
      DOUBLE PRECISION VINT
      SAVE/PYINT1/
      DOUBLE PRECISION f2dis,f1dis,sdis,svmd1,svmd2,sigvm,x,q2
      DOUBLE PRECISION dipol,rvmd
      DOUBLE PRECISION mrho2,alpha,eps,eta,pmass,pmass2,gevmb
      DOUBLE PRECISION w2,gamma2,nu,convf2,convf1,convf1d,conv,sigh,pi
      DOUBLE PRECISION f1,f2,df2allm,df2allml,df1allm,df1allml,dr,drl
      DOUBLE PRECISION w2l,gamma2l,nul,xl
      REAL f2allm,rx,rq2,rxl
C...Local arrays.
      DOUBLE PRECISION XPQ,XPAR,YPAR
      INTEGER i
      DIMENSION XPQ(-25:25),XPAR(4),YPAR(4)
     
C...X and Y parameters of sigmatot = X * s**epsilon + Y * s**(-eta).
      DATA XPAR/2*13.63D0,10.01D0,0.970D0/
      DATA YPAR/2*31.79D0,-1.51D0,-0.146D0/
      DATA EPS/0.0808D0/,ETA/-0.4525D0/
      DATA MRHO2/0.591822D0/,ALPHA/7.297352533D-3/,PMASS/0.93827D0/,
     &     GEVMB/0.389379292D0/,pi/3.14159265358979324D0/

      EXTERNAL PYPDFU 

c....Kinematic factors and constants
      pmass2=pmass**2
      nul=-1D0
      gamma2l=0D0
      w2l=-1D0
      xl=0D0
      
      if ((x.gt.0D0).and.(x.le.1)) then
       w2=pmass2+q2*(1-x)/x
      else
       f1=0
       f2=0
       return
      endif
      
      if (w2.lt.4) then
       w2l=w2
       w2=4D0
       xl=x
       nul=(w2l-pmass2+q2)/(2D0*pmass)
      endif
      
      nu=(w2-pmass2+q2)/(2D0*pmass)

      if (nu.gt.0D0) then
       gamma2=q2/(nu**2)
       if (nul.gt.0D0) then
        gamma2l=q2/(nul**2)
       endif
      else
       f1=0
       f2=0
       return
      endif
c....For W2<4, we don't trust the PYTHIA F2, so we calculate F2-ALLM(W2,Q2)
c....The real kinematics have and "l" at the end, the kinematics without
c....the "l" are the ones we get by setting W2 to 4
c....Output: f2allm=F2-ALLM(W2=4,Q2),f2allml=F2-ALLM(W2=w2l,Q2)       
      if (w2l.gt.0D0) then
       rq2=real(q2)
       rxl=real(xl)
       call R1990(Q2,Xl,DRl)
       df2allml=dble(f2allm(rxl,rq2))
       df1allml=1D0/(2D0*x*(1+drl))*(1+gamma2l)*df2allml
       x=q2/(2D0*pmass*nu)
       rx=real(x)
       call R1990(Q2,X,DR)
       df2allm=dble(f2allm(rx,rq2))
       df1allm=1D0/(2D0*x*(1+dr))*(1+gamma2)*df2allm
      endif 
c....This factor is needed to convert the Pythia virtual gamma cross 
c....section for VMD to the same level as F2
c....The kinematic factors making the (ep) cross section out of F2 are
c....provided by RADGEN      
      conv=q2*(1D0-x)/(4D0*pi**2*alpha)/gevmb
c....Pythia PDF call, sum PDFs to F2
      call PYPDFU(2212,X,Q2,XPQ)
       f2dis=1D0/9D0*(XPQ(1)+XPQ(-1)+XPQ(3)+XPQ(-3))+
     &       4D0/9D0*(XPQ(2)+XPQ(-2))
c....Suppression factor for DIS
      if (MSTP(19).eq.0) then
       sdis=1.
      else  
       sdis=q2/(q2+mrho2)
       if (MSTP(19).gt.1) then
        sdis=sdis**2
       endif
      endif
C....Sum of Hadronic (Vector Meson) cross sections * Photon couplings
C....const.
      sigh=0.
      do 10 i=1,4
       sigh=sigh+alpha/PARP(160+i)*(XPAR(i)*w2**eps+YPAR(i)*w2**eta) 
   10 continue
C....W2/Q2 suppression of VMD and (1+epsilon R_VMD)
      svmd1=(w2/(w2+q2))**MSTP(20)
      if (MSTP(20).eq.0) then
       dipol=2.575D0
      else
       dipol=2D0
      endif
      if (MSTP(17).eq.3) then
       rvmd=PARP(165)*(q2/mrho2)**0.61
      else
c    ...Attention: This is only good for MSTP(17)=4, i.e., the Pythia
c    ...default       
       rvmd=PARP(165)*(4.*mrho2*q2)/(mrho2+q2)**2
      endif
C  .... Dipole factor for VMD      
      svmd2=(mrho2/(mrho2+q2))**dipol
C.....virtual photon xsec for VMD
      sigvm=svmd1*svmd2*sigh
      convf2=(1D0+rvmd)/(1D0+gamma2)
      convf1=1D0/(2D0*x)
c.....Total "F2"
      f2=sdis*f2dis+conv*convf2*sigvm
      f1dis=(1.D0+gamma2)/(2.D0*x)*f2dis
      f1=sdis*f1dis+conv*convf1*sigvm
      if (w2l.gt.0D0) then
C.....Here we scale F2-ALLM(W2=w2l,Q2) by the factor 
C.....F2-PYTH(W2=4,Q2)/F2-ALLM(W2=4,Q2) (normalize ALLM to PYTHIA model)
       f2=f2/df2allm*df2allml
       f1=f1/df1allm*df1allml
      endif
      RETURN
      END
