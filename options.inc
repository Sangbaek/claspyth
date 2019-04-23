       character*120 cl_target,cl_beam,cl_hispath
       character*120 datfilename
	 real cl_beam_energy,cl_zpos
       integer cl_triggers,cl_nprint,cl_verb,cl_fileev
       logical clasdisOK,cl_radgen,cl_dis,cl_rndm
       common /ioptions/ cl_triggers,cl_nprint,cl_verb,cl_fileev
       common /roptions/ cl_beam_energy,cl_zpos
       common /coptions/ cl_target,cl_hispath,cl_beam,datfilename
       common /loptions/ cl_radgen,cl_dis,cl_rndm

