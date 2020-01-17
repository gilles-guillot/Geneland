*     Computation of index a 
*     of Rousset J. Evol. Biol. 13 2000 58-62 Eq. 3
*     (after corrections of typos)
*     for genotype at one locus only 
      subroutine rousset1(iloc,nindiv,nloc,nloc2,nall,nallmax,z,a)
      implicit none
*     data
      integer nindiv,nloc,nloc2,nall,nallmax,z
      dimension z(nindiv,nloc2),nall(nloc)
*     output
      double precision a(nindiv,nindiv)
*     working variables 
      integer iloc,iall,iindiv1,iindiv2
      double precision sum,b,w,ssb,ssw

c$$$      write(*,*) 'debut rousset1'
c$$$      write(*,*) 'iloc=',iloc
c$$$      write(*,*) 'nindiv=',nindiv
c$$$      write(*,*) 'nloc=',nloc
c$$$      write(*,*) 'nloc2=',nloc2
c$$$      write(*,*) 'nall=',nall
c$$$      write(*,*) 'nallmax=',nallmax
c$$$      write(*,*) 'z=',z
c$$$      write(*,*) 'a=',a

      sum = 0.d0
      do iindiv1 = 1,nindiv-1
         do iindiv2 = iindiv1+1,nindiv
c            write(*,*) 'iindiv1=',iindiv1
c            write(*,*) 'iindiv2=',iindiv2
            b = ssb(iindiv1,iindiv2,iloc,z,
     &           nindiv,nloc,nloc2,nall,nallmax)
            w = ssw(iindiv1,iindiv2,iloc,z,
     &           nindiv,nloc,nloc2,nall,nallmax)
            a(iindiv1,iindiv2) = 2*b - w 
*            a(iindiv1,iindiv2) = b 
            sum = sum + w
         enddo
      enddo
      do iindiv1 = 1,nindiv-1
         do iindiv2 = iindiv1+1,nindiv
            a(iindiv1,iindiv2) = a(iindiv1,iindiv2)*nindiv*(nindiv-1)/
     &           (4*sum)
         enddo
      enddo
      end subroutine rousset1
 



************************************************************************
*     Computation of index a 
*     of Rousset J. Evol. Biol. 13 2000 58-62 Eq. 3
*     after corrections of typos, (F. Rousset personnal communication)
*     for multilocus genotype
      subroutine areq3(nindiv,nloc,nloc2,nall,nallmax,z,a)
      implicit none
*     data
      integer nindiv,nloc,nloc2,nall,nallmax,z
      dimension z(nindiv,nloc2),nall(nloc)
*     output
      double precision a(nindiv,nindiv)
*     working variables 
      integer iloc,iall,iindiv1,iindiv2
      double precision sum,b,w,ssb,ssw,num,denom

c$$$      write(*,*) 'debut rousset1'
c$$$      write(*,*) 'iloc=',iloc
c$$$      write(*,*) 'nindiv=',nindiv
c$$$      write(*,*) 'nloc=',nloc
c$$$      write(*,*) 'nloc2=',nloc2
c$$$      write(*,*) 'nall=',nall
c$$$      write(*,*) 'nallmax=',nallmax
c$$$      write(*,*) 'z=',z
c$$$      write(*,*) 'a=',a

      sum = 0.d0
      do iindiv1 = 1,nindiv-1
         do iindiv2 = iindiv1+1,nindiv
            a(iindiv1,iindiv2) = 0
         enddo
      enddo

      do iloc = 1,nloc
         do iindiv1 = 1,nindiv-1
            do iindiv2 = iindiv1+1,nindiv
               b = ssb(iindiv1,iindiv2,iloc,z,
     &              nindiv,nloc,nloc2,nall,nallmax)
               w = ssw(iindiv1,iindiv2,iloc,z,
     &              nindiv,nloc,nloc2,nall,nallmax)
               a(iindiv1,iindiv2) = a(iindiv1,iindiv2) + 2*b - w 
               sum = sum + w
            enddo
         enddo
      enddo
      do iindiv1 = 1,nindiv-1
         do iindiv2 = iindiv1+1,nindiv
            a(iindiv1,iindiv2) = a(iindiv1,iindiv2)*
     &           dble(nindiv*(nindiv-1))/(4*sum)
         enddo
      enddo
c$$$      do iloc = 1,nloc
c$$$         do iindiv1 = 1,nindiv-1
c$$$            do iindiv2 = iindiv1+1,nindiv
c$$$               b = ssb(iindiv1,iindiv2,iloc,z,
c$$$     &              nindiv,nloc,nloc2,nall,nallmax)
c$$$               w = ssw(iindiv1,iindiv2,iloc,z,
c$$$     &              nindiv,nloc,nloc2,nall,nallmax)
c$$$               a(iindiv1,iindiv2) = a(iindiv1,iindiv2) + b 
c$$$               sum = sum + w
c$$$            enddo
c$$$         enddo
c$$$      enddo
c$$$      do iindiv1 = 1,nindiv-1
c$$$         do iindiv2 = iindiv1+1,nindiv
c$$$            a(iindiv1,iindiv2) = a(iindiv1,iindiv2)*
c$$$     &           dble(nindiv*(nindiv-1))/(2*sum) -0.5d0
c$$$         enddo
c$$$      enddo
      end subroutine areq3
************************************************************************



 


************************************************************************
*     Computation of index a 
*     of Rousset J. Evol. Biol. 13 2000 58-62 Eq. 3
*     (after corrections of typos)
*     for multilocus genotype
      subroutine areq4(nindiv,nloc,nloc2,nall,nallmax,z,a)
      implicit none
*     data
      integer nindiv,nloc,nloc2,nall,nallmax,z
      dimension z(nindiv,nloc2),nall(nloc)
*     output
      double precision a(nindiv,nindiv)
*     working variables 
      integer iloc,iall,iindiv1,iindiv2
      double precision sum,b,w,ssb,ssw,num,denom

c$$$      write(*,*) 'debut rousset1'
c$$$      write(*,*) 'iloc=',iloc
c$$$      write(*,*) 'nindiv=',nindiv
c$$$      write(*,*) 'nloc=',nloc
c$$$      write(*,*) 'nloc2=',nloc2
c$$$      write(*,*) 'nall=',nall
c$$$      write(*,*) 'nallmax=',nallmax
c$$$      write(*,*) 'z=',z
c$$$      write(*,*) 'a=',a

      sum = 0.d0
      do iindiv1 = 1,nindiv-1
         do iindiv2 = iindiv1+1,nindiv
            a(iindiv1,iindiv2) = 0
         enddo
      enddo

      do iloc = 1,nloc
         do iindiv1 = 1,nindiv-1
            do iindiv2 = iindiv1+1,nindiv
               b = ssb(iindiv1,iindiv2,iloc,z,
     &              nindiv,nloc,nloc2,nall,nallmax)
               w = ssw(iindiv1,iindiv2,iloc,z,
     &              nindiv,nloc,nloc2,nall,nallmax)
               a(iindiv1,iindiv2) = a(iindiv1,iindiv2) + b 
               sum = sum + w
            enddo
         enddo
      enddo
      do iindiv1 = 1,nindiv-1
         do iindiv2 = iindiv1+1,nindiv
            a(iindiv1,iindiv2) = a(iindiv1,iindiv2)*
     &           dble(nindiv*(nindiv-1))/(2*sum) -0.5d0
         enddo
      enddo
      end subroutine areq4
************************************************************************





************************************************************************
      double precision function ssw(iindiv1,iindiv2,iloc,z,
     &     nindiv,nloc,nloc2,nall,nallmax)
      implicit none
      integer iindiv1,iindiv2,iloc,z,nindiv,nloc,nloc2,nall,nallmax
      dimension z(nindiv,nloc2),nall(nloc)
      integer iindiv,igene,listindiv,iall
      double precision x,propindiv
      dimension listindiv(2)

*     propindiv = Xi.:u of Rousset
*     x         = Xij:u

c$$$      write(*,*) ''
c      write(*,*) 'debut ssw'
c$$$      write(*,*)'iindiv1=',iindiv1
c$$$      write(*,*)'iindiv2=',iindiv2
c$$$      write(*,*) 'iloc=',iloc
c$$$      write(*,*) 'nindiv=',nindiv
c$$$      write(*,*) 'nloc=',nloc
c$$$      write(*,*) 'nloc2=',nloc2
c$$$      write(*,*) 'nall=',nall
c$$$      write(*,*) 'nallmax=',nallmax
c$$$      write(*,*) 'z=',z

      listindiv(1) = iindiv1
      listindiv(2) = iindiv2

      ssw = 0.d0
      do iindiv = 1,2
         do iall = 1,nall(iloc)
            propindiv = 0
            if(z(listindiv(iindiv),2*iloc-1) .eq. iall) 
     &           propindiv = propindiv+1
            if(z(listindiv(iindiv),2*iloc)   .eq. iall) 
     &           propindiv = propindiv+1
            propindiv = propindiv/2.d0
c            write(*,*) 'propindiv=',propindiv
            do igene = 1,2
               x = 0d0
               if(z(listindiv(iindiv),2*iloc-1+igene-1) .eq. iall) 
     &              x = 1.d0
               ssw = ssw + (x - propindiv)**2
            enddo
         enddo
      enddo
c      write(*,*) 'ssw=',ssw
      end 
************************************************************************






************************************************************************
      double precision function ssb(iindiv1,iindiv2,iloc,z,
     &     nindiv,nloc,nloc2,nall,nallmax)
      implicit none
      integer iindiv1,iindiv2,iloc,z,nindiv,nloc,nloc2,nall,nallmax
      dimension z(nindiv,nloc2),nall(nloc)
      integer iindiv,igene,listindiv,iall
      double precision x,propindiv,proppair
      dimension listindiv(2)

*     propindiv = Xi.:u of Rousset
*     proppair  = X..:u
*     x         = Xij:u
c$$$      write(*,*) ''
c      write(*,*) 'debut ssb'
c$$$      write(*,*)'iindiv1=',iindiv1
c$$$      write(*,*)'iindiv2=',iindiv2
c$$$      write(*,*) 'iloc=',iloc
c$$$      write(*,*) 'nindiv=',nindiv
c$$$      write(*,*) 'nloc=',nloc
c$$$      write(*,*) 'nloc2=',nloc2
c$$$      write(*,*) 'nall=',nall
c$$$      write(*,*) 'nallmax=',nallmax
c$$$      write(*,*) 'z=',z

      listindiv(1) = iindiv1
      listindiv(2) = iindiv2

      ssb = 0.d0
      do iall = 1,nall(iloc)
c         write(*,*) 'iall=',iall
*     computes proportion of allele iall in current pair
         proppair = 0.d0  
         if(z(iindiv1,2*iloc-1) .eq. iall) 
     &        proppair = proppair+1
         if(z(iindiv1,2*iloc)   .eq. iall) 
     &        proppair = proppair+1
         if(z(iindiv2,2*iloc-1) .eq. iall) 
     &        proppair = proppair+1
         if(z(iindiv2,2*iloc)   .eq. iall) 
     &        proppair = proppair+1
         proppair = proppair/4.d0   
c         write(*,*) 'proppair=',proppair
         do iindiv = 1,2
            propindiv = 0.d0
            if(z(listindiv(iindiv),2*iloc-1) .eq. iall) 
     &           propindiv = propindiv+1
            if(z(listindiv(iindiv),2*iloc)   .eq. iall) 
     &           propindiv = propindiv+1
            propindiv = propindiv/2.d0
c            write(*,*) 'propindiv=',propindiv
            ssb = ssb + (propindiv - proppair)**2
c            write(*,*) 'ssb=',ssb
c            write(*,*) 'propindiv - proppair=',propindiv - proppair
         enddo
      enddo
      ssb = 2*ssb
c      write(*,*) 'ssb=',ssb
      end
************************************************************************

