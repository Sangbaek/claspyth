       character*120 cl_target,cl_beam,cl_hispath
       character*120 datfilename
       real cl_beam_energy,cl_prom,cl_zpos,cl_zwidth,cl_rast,cl_tmin,cl_emin
       integer cl_triggers,cl_nprint,cl_verb,cl_fileev,cl_pid
       logical clasdisOK,cl_radgen,cl_dis,cl_rndm,cl_ntno,cl_dumpdecayfile
       common /ioptions/ cl_triggers,cl_nprint,cl_verb,cl_fileev,cl_pid
       common /roptions/ cl_beam_energy,cl_prom,cl_zpos,cl_zwidth,cl_rast
     6,cl_tmin,cl_emin
       common /coptions/ cl_target,cl_hispath,cl_beam,datfilename
       common /loptions/ cl_radgen,cl_dis,cl_rndm,cl_ntno,cl_dumpdecayfile
