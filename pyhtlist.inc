*
      REAL elee,elet,elef,elex,eley,elez,xb,q2,w
      INTEGER nele,eleevtclass,eleh,elerun,eletype,eleevent
      COMMON /ELEREC/ nele,eleevtclass,eleh,elerun,eletype,eleevent,elee
     + ,elet,elef,elex,eley,elez,xb,q2,w
*
      REAL pi0e,pi0t,pi0f,pipe,pipt,pipf,pime,pimt,pimf,proe,prot,prof
     + ,game,gamt,gamf,gam2t,gam2e,gam2f,neue,neut,neuf,pim,pipx,pipy
     + ,pipz,pimx,pimy,pimz,prox,proy,proz
      INTEGER npi0,npip,npim,npro,ngam,nneu
      COMMON /HADREC/ npi0,npip,npim,npro,ngam,nneu,pi0e,pi0t,pi0f,pipe
     + ,pipt,pipf,pime,pimt,pimf,proe,prot,prof,game,gamt,gamf,gam2t
     + ,gam2e,gam2f,neue,neut,neuf,pim,pipx,pipy,pipz,pimx,pimy,pimz
     + ,prox,proy,proz
*
      REAL gelee,gelet,gelef,gelex,geley,gelez,gelewgt,gxb,gq2,gw
      INTEGER gnele,gelebits,ghpar,gproh,geleh
      COMMON /ELEGEN/ gnele,gelebits,ghpar,gproh,geleh,gelee,gelet
     +,gelef,gelex,geley,gelez,gelewgt,gxb,gq2,gw
*
      REAL gpose,gpost,gposf,gpmas,gele2,gelt2,gelf2
      INTEGER gnpos,nspos,nsel2
      COMMON /POSGEN/ gnpos,nspos,nsel2,
     +gpose,gpost,gposf,gpmas,gele2,gelt2,gelf2
c234567890123456789012345678901234567890123456789012345678901234567890123456
*
      REAL gpi0e,gpi0t,gpi0f,gpipe,gpipt,gpipf,gpime,gpimt,gpimf,gproe
     + ,gprot,gprof,ggame,ggamt,ggamf,ggam2t,ggam2e,ggam2f,gneue,gneut
     + ,gneuf,gpipx,gpipy,gpipz,gpimx,gpimy,gpimz,gprox,gproy,gproz
     + ,grope,gropt,gropf,gro0e,gro0t,gro0f
     + ,gpipm,gpimm,gpi0m,gro0m,gropm
     + ,gomee,gomet,gomef,gomem
      INTEGER gnpi0,gnpip,gnpim,gnpro,gngam,gnneu,gnrop,gnro0,gnome
      INTEGER nsproc,nsele,nspi0,nspip,nspim
      COMMON /HADGEN/ gnpi0,gnpip,gnpim,gnpro,gngam
     + ,gnneu,gnrop,gnro0,gnome
     + ,nsproc,nsele,nspi0,nspip,nspim
     + ,gpi0e,gpi0t,gpi0f
     + ,gpipe,gpipt,gpipf
     + ,gpime,gpimt,gpimf
     + ,gproe,gprot,gprof
     + ,gneue,gneut,gneuf
     + ,grope,gropt,gropf
     + ,gro0e,gro0t,gro0f
     + ,gomee,gomet,gomef
     + ,gpipm,gpimm,gpi0m,gro0m,gropm,gomem
*
       real v4el0(4),v4pr0(4)
      COMMON /kinem0/v4el0,v4pr0

*
      INTEGER iproctpe
      REAL radq2,radnu,radweight,r_sigrad,r_sigcor,grade,gradt,gradf
      COMMON /RADGN/iproctpe,radq2,radnu,radweight,r_sigrad,r_sigcor
     + ,grade,gradt,gradf
