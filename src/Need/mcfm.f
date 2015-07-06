      program mcfm
      implicit none
      character*72 inputfile,workdir
      call determinefilenames(inputfile,workdir)
      call mcfm_main(inputfile,workdir)
      end
