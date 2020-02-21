# claspyth (nocernlib version)
PYTHIA based generator for clas12 (HERMES version with some changes)

To install
`git clone https://github.com/JeffersonLab/claspyth-nocernlib.git`

`cd claspyth` 
`make`

Need to define the path for the pdf sets using env variable CLASPYTHIA_DECLIST to define the path to  directory containing the decay-list file pythdeclist-curr.dat

Need to create a tmp directory and define the path for the pdf sets using env variable CLASDIS_PDF to define the path to pdf directory

Example
`./claspyth --trig 10000000 --nmax 10000  < input.10.6gev.with-comments`
will write `gemc lund type` data files with 10K events in the directory `tmp` (should be created by hand) with total of 10M events.
A number of control parameters, including JETSET settings are listed in the editable input data file input.10.6gev.with-comments

To get more command line options

`./claspyth --help`
```
        option  value      default    comment
   --trig nevents         10      number of triggers
   --nmax nevents     200      number of ev in 1 file
   --beam beam_type    gamma/e-     e- beam
      --targ target         p+       n0 for neutron
   --ebeam beam_energy     10.6  e- momentum in GeV
   --zpos  z_position      0.0    e- momentum in GeV
   --radgen                   include radgen
   --path hisid          w9/p       write dir
   --print nprint       1000   print nprint event
   --dis                .FALSE.     all dis events in
   --v     verblev       0     verbos. level

 
