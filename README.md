# claspyth (nocernlib version)
PYTHIA based generator for clas12 (HERMES version with some changes)

To install
`git clone https://github.com/JeffersonLab/claspyth.git`

`cd claspyth` 
`make`

Need to define the path for the pdf sets using env variable CLASPYTHIA_DECLIST to define the path to  directory containing the decay-list file pythdeclist-curr.dat and the input file PYTHIA-input.dat (current directory). In case of running in clas12-mcgen CLASPYTHIA_DECLIST  should be set to $PWD/claspyth

Need to create a tmp directory if running without --docker option

Example
`./claspyth --trig 10000000 --nmax 10000  < input.10.6gev.with-comments`
will write `gemc lund type` data files with 10K events in the directory `tmp` (should be created by hand) with total of 10M events.
A number of control parameters, including JETSET settings are listed in the editable input data file input.10.6gev.with-comments

To get more command line options

`./claspyth --help`
```
        option  value      default    comment
   --trig nevents         10      number of triggers
   --seed seed_number     0      use the local time as seed
   --nmax nevents     200      number of ev in 1 file
   --beam beam_type    gamma/e-     e- beam
   --targ target         p+       n0 for neutron
   --prom  proton mom   0-clas12  p+ momentum in GeV 
   --ebeam beam_energy     10.6  e- momentum in GeV
   --xpos  z_position      0.0   beam x pos. in cm 
   --ypos  y_position      0.0   beam y pos. in cm
   --zpos  z_position      0.0   target z pos. in cm 
   --radgen                   include radgen
   --path hisid          w9/p       write dir
   --print nprint       1000   print nprint event
   --dis                .FALSE.     all dis events in
   --v     verblev       0     verbos. level

 
