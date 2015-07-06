      subroutine mcfm_main(inputfile,workdir)
************************************************************************
*                                                                      *
*  This is the main program for MCFM                                   *
*                                                                      *
*  The sequence of calls should always be:                             *
*   call mcfm_init          : basic variable initialization, print-out *
*   call mcfm_vegas(warmup) : warm-up the Vegas grid                   *
*   call mcfm_vegas(accum)  : accumulate results                       *
*   call mcfm_exit          : final processing and print-out           *
*                                                                      *
************************************************************************
      implicit none
      include 'constants.f'
      include 'efficiency.f'
      include 'maxwt.f'
      include 'eventbuffer.f'
      include 'gridinfo.f'
      integer itmx1,ncall1,itmx2,ncall2,itmxplots
      double precision integ,integ_err
      logical dryrun
      integer i,pflav,pbarflav
      double precision p(mxpart,4),wt
      character*72 inputfile,workdir
      common/iterat/itmx1,ncall1,itmx2,ncall2
      common/dryrun/dryrun
       
* basic variable initialization, print-out
      call mcfm_init(inputfile,workdir)

* in initial phases, we don't want any unweighting to take place.
* this will be set true in the first call to getevent.
      unweight = .false.
* if we're reading in a grid, there's no need to do any warming-up
      if (readin) dryrun=.true.

* This is the mcfm_vegas(warmup) call
* The Vegas parameters are those read from options.DAT for
* the warm-up stage (itmx1,ncall1) and binning should only take
* place if dryrun is set to true

* There are now 3 modes of operation:
*   dryrun = .false. : warmup, then freeze grid and accumulate
*   dryrun = .true. , readin = .false. : accumulate during warmup
*   dryrun = .true. , readin = .true.  : accumulate with frozen grid

      if ((dryrun .eqv. .false.) .or. 
     .    ((dryrun) .and. (readin .eqv. .false.))) then
* Initialize efficiency variables      
        njetzero=0
        ncutzero=0
        ntotzero=0
        ntotshot=0
        call mcfm_vegas(0,itmx1,ncall1,dryrun,integ,integ_err)
        itmxplots=itmx1
      endif
* This is the mcfm_vegas(accum) call
* This takes place only if dryrun is false
* The Vegas parameters are those read from options.DAT for
* the results stage (itmx2,ncall2) and binning takes place (.true.)
* wtmax may have been set during the dry run, so re-set here :
      wtmax = 0d0
      if ((dryrun .eqv. .false.) .or. 
     .    ((dryrun) .and. (readin .eqv. .true.))) then
* Initialize efficiency variables      
        njetzero=0
        ncutzero=0
        ntotzero=0
        ntotshot=0
        call mcfm_vegas(1,itmx2,ncall2,.true.,integ,integ_err)
        itmxplots=itmx2
      endif
      
* So far we have not used VEGAS to generate any events.
* Make sure future calls to "getevent" are aware of this :
      numstored = 0

      if (evtgen) then
        write(6,*) 'Generate events :'
        do i=1,500
          call mcfm_getevent(p,wt,pflav,pbarflav)
          call fill_stdhep(p,0,0,wt)
c         call write_stdhep(6)
        enddo
      endif

* final processing and print-out
      call mcfm_exit(itmxplots,integ,integ_err)
      
      stop
      end
       
