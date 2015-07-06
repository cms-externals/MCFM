      subroutine userplotter(p,wt,wt2,nd)
c--- Variables passed to this routine:
c
c---      p:  4-momenta of jets in the format p(i,4)
c---          with the particles numbered according to the input file
c---          and components labelled by (px,py,pz,E)??????????
c
c---     wt:  weight of this event
c
c---    wt2:  weight^2 of this event
c
c---     nd:  an integer specifying the dipole number of this contribution
c---          (if applicable), otherwise equal to zero
      implicit none
      include 'constants.f'
      double precision p(mxpart,4),wt,wt2
      integer nd
      return
      end
