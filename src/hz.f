      subroutine mcmchz(q,qtmp,zcodom,zdom,zhap,codom,dom,hap,
     &     npop,npopmax,nindiv,f,nlocd,nloch,nalmax,
     &     alphadmix,alphadmixtmp,a,b,c,atmp,btmp,ctmp,amax,bmax,cmax,
     &     dist,nchpath,path,nit,thinning,calluda,calludb,calludc,
     &     calludq,deltab,compar,nitstor,qout,aout,bout,cout)
      implicit none
      integer npop,npopmax,nindiv,nlocd,nloch,nalmax,zcodom,zdom,zhap,
     &     nit,thinning,nchpath,compar,codom,dom,hap,nitstor
      double precision f,q,qtmp,alphadmix,alphadmixtmp,a,b,c,
     &     atmp,btmp,ctmp,amax,bmax,cmax,dist,deltab,
     &     qout,aout,bout,cout
      dimension zcodom(nindiv,2*nlocd),zdom(nindiv,nlocd),
     &     zhap(nindiv,nloch),f(npopmax,nlocd,nalmax),
     &     q(nindiv,npopmax),qtmp(nindiv,npopmax),
     &     alphadmix(nindiv,npopmax), alphadmixtmp(nindiv,npopmax),
     &     dist(nindiv,npopmax),a(npopmax),b(npopmax),c(npopmax),
     &     atmp(npopmax),btmp(npopmax),ctmp(npopmax),
     &     qout(nitstor,nindiv,npopmax),aout(nitstor,npopmax),
     &     bout(nitstor,npopmax),cout(nitstor,npopmax)
      integer iit,ipop,iindiv,calluda,calludb,calludc,calludq,iitstor
      double precision pct,lpriorq,llike
      character*255 path, fileq, filea, fileb, filec, filellike, 
     &     filelpriorq

 2000 format (300(1x,e15.8,1x))


*     init RNG
      call rndstart()
c$$$      do iindiv=1,nindiv
c$$$         write(*,*) 'q(',iindiv,')=',(q(iindiv,ipop),ipop=1,npop)
c$$$         write(*,*) 'qtmp(',iindiv,')=',(qtmp(iindiv,ipop),ipop=1,npop)
c$$$      enddo

      fileq = path(1:nchpath) // "q.txt"
      filea = path(1:nchpath) // "a.txt"
      fileb = path(1:nchpath) // "b.txt"
      filec = path(1:nchpath) // "c.txt"
      filellike = path(1:nchpath) // "llike.txt"
      filelpriorq = path(1:nchpath) // "lpriorq.txt"

c$$$      open(9,file=fileq)
c$$$      open(10,file=filea)
c$$$      open(11,file=fileb)
c$$$      open(12,file=filec)
c$$$      open(17,file=filellike)
c$$$      open(18,file=filelpriorq)

      iitstor = 0
      do iit=1,nit
         if(mod(iit,thinning) .eq. 0) then
            pct = dble(iit)/dble(nit)*100.
            call dblepr('                     ',-1,pct,1)
            iitstor = iitstor + 1 
            do ipop=1,npopmax
               do iindiv=1,nindiv
                  qout(iitstor,iindiv,ipop) = q(iindiv,ipop)
                  aout(iitstor,ipop) = a(ipop)
                  bout(iitstor,ipop) = b(ipop)
                  cout(iitstor,ipop) = c(ipop)
               enddo
            enddo
c$$$            do iindiv=1,nindiv
c$$$               write(9,2000) (sngl(q(iindiv,ipop)),ipop=1,npopmax)
c$$$            enddo
c$$$            write(10,2000) (sngl(a(ipop)),ipop=1,npopmax)
c$$$            write(11,2000) (sngl(b(ipop)),ipop=1,npopmax)
c$$$            write(12,2000) (sngl(c(ipop)),ipop=1,npopmax)
         endif
         if(calludq .eq. 1) then
            call udq(q,qtmp,zcodom,zdom,zhap,codom,dom,hap,
     &     npop,npopmax,nindiv,f,nlocd,nloch,nalmax,alphadmix)
         endif
         if(calluda .eq. 1) then
            call uda(npop,npopmax,nindiv,a,atmp,amax,alphadmix,
     &           alphadmixtmp,q,dist,compar)
         endif
         if(calludb .eq. 1) then
            call udb(npop,npopmax,nindiv,a,b,c,btmp,bmax,alphadmix,
     &     alphadmixtmp,q,dist,compar,deltab)
         endif
         if(calludc .eq. 1) then
            call udc(npop,npopmax,nindiv,a,b,c,ctmp,cmax,alphadmix,
     &           alphadmixtmp,q,dist,compar)
        endif
c     write log- prior and likelihood
c         call ppost(q,zcodom,npop,npopmax,nindiv,f,nlocd,nalmax,alphadmix,
c     &        lpriorq,llike)
      enddo
c$$$      close(9)
c$$$      close(10)
c$$$      close(11)
c$$$      close(12)
c$$$      close(17)
c$$$      close(18)
      call rndend()
      end subroutine mcmchz
      

***********************************************************
*     compute and write prior and likelihood of current state
      subroutine ppost(q,zcodom,npop,npopmax,nindiv,f,nlocd,nalmax,
     &     alphadmix,lpriorq,llike)
      implicit none
      integer npop,npopmax,nindiv,nlocd,nalmax,zcodom
      double precision f,q,alphadmix
      dimension zcodom(nindiv,2*nlocd),f(npop,nlocd,nalmax),
     &     q(nindiv,npopmax),alphadmix(nindiv,npopmax)
      integer ipop,ipop1,ipop2,iindiv,iloc
      double precision llike,lpriorq,gglgamfn,sa,sb
      
c$$$      write(*,*) 'nlocd=',nlocd
c$$$      write(*,*) 'npop=',npop
c$$$      write(*,*) 'nalmax=',nalmax

c$$$      write(*,*) 'q=',q
c$$$      write(*,*) 'f=',f

      lpriorq = 0
      do iindiv=1,nindiv
         sa =0
         do ipop=1,npop
            sa = sa + alphadmix(iindiv,ipop)
            lpriorq = lpriorq - gglgamfn(alphadmix(iindiv,ipop)) + 
     &           (alphadmix(iindiv,ipop)-1)*dlog(q(iindiv,ipop))
         enddo
         lpriorq = lpriorq + gglgamfn(sa)
      enddo

      llike = 0 
      do iindiv = 1,nindiv
c         write(*,*) 'iindiv=',iindiv
         do iloc=1,nlocd
c           write(*,*) 'iloc=',iloc
            sa = 0
            sb = 0
            do ipop=1,npop
c              write(*,*) 'ipop=',ipop
c              write(*,*) 'q(',iindiv,',',ipop,')=',q(iindiv,ipop)
               if(zcodom(iindiv,2*iloc-1) .ne. -999) then 
                  sa = sa + 
     &             q(iindiv,ipop)*f(ipop,iloc,zcodom(iindiv,2*iloc-1))
               endif
               if(zcodom(iindiv,2*iloc) .ne. -999) then 
                  sb = sb + 
     &              q(iindiv,ipop)*f(ipop,iloc,zcodom(iindiv,2*iloc))
               endif
            enddo
            if(zcodom(iindiv,2*iloc-1) .ne. -999) then 
               llike = llike + dlog(sa) 
            endif
            if(zcodom(iindiv,2*iloc) .ne. -999) then 
               llike = llike + dlog(sb)
            endif
c           write(*,*) 'sa=',sa
c           write(*,*) 'sb=',sb
            if(((zcodom(iindiv,2*iloc-1) .ne. zcodom(iindiv,2*iloc)) 
     &        .and.   (zcodom(iindiv,2*iloc-1) .ne. -999)) .and. 
     &          (zcodom(iindiv,2*iloc) .ne. -999)) then
               llike = llike + dlog(2.d0)
            endif
c            write(*,*) 'llike=',llike
         enddo
c           write(*,*) 'sa=',sa
c           write(*,*) 'sb=',sb
c         write(*,*) 'llike=',llike
      enddo
c      write(17,*) llike
c      write(18,*) lpriorq
c      write(*,*) 'fin'
      end subroutine ppost


***********************************************************
      subroutine udq(q,qtmp,zcodom,zdom,zhap,codom,dom,hap,
     &     npop,npopmax,nindiv,f,nlocd,nloch,nalmax,alphadmix)
      implicit none
      integer npop,npopmax,nindiv,nlocd,nloch,nalmax,zcodom,zdom,zhap,
     &     codom,dom,hap
      double precision f,q,qtmp,alphadmix
      dimension zcodom(nindiv,2*nlocd),zdom(nindiv,nlocd),
     &     zhap(nindiv,nloch),f(npop,nlocd,nalmax),
     &     q(nindiv,npopmax),qtmp(nindiv,npopmax),
     &     alphadmix(nindiv,npopmax)
      integer ipop,ipop1,ipop2,iindiv
      double precision ggrunif,ggrbinom,contriblikeq,contprq,delta,
     &     lratio,bern,lratiobis,lpriorq,lpriorqtmp,llike,lliketmp
c      write(*,*) 'debut udq'
      do iindiv=1,nindiv
         do ipop=1,npop
            qtmp(iindiv,ipop) =  q(iindiv,ipop)
         enddo
*     propose new q
*     chose two pop at random
         ipop1 = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
         ipop2 = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
         do while(ipop2 .eq. ipop1)
            ipop2 = 1 + 
     &           idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
         enddo
*     sample small increment
         delta = ggrunif(0.d0,1.d0)/dble(10*npop)
         qtmp(iindiv,ipop1) =  q(iindiv,ipop1) + delta 
         qtmp(iindiv,ipop2) =  q(iindiv,ipop2) - delta 
         lratio = 0 
         if((1-qtmp(iindiv,ipop1) .gt. 1.d-300) .and. 
     &        (qtmp(iindiv,ipop2) .gt. 1.d-300)) then 
            lratio = lratio + 
     &           contriblikeq(q,qtmp,iindiv,
     &           zcodom,zdom,zhap,codom,dom,hap,
     &           npop,npopmax,nindiv,f,nlocd,nloch,nalmax)
            lratio = lratio + 
     &           contprq(q,qtmp,iindiv,npop,npopmax,nindiv,alphadmix)
            bern = 1
            if(lratio .le. 0) then
               bern = ggrbinom(1.d0,dexp(lratio))
            endif
            if(idint(bern) .eq. 1) then 
               do ipop=1,npop
                  q(iindiv,ipop) =  qtmp(iindiv,ipop)
               enddo
            endif
         endif
      enddo
      end subroutine udq

 
      

************************************************************************
*     log of contribution of likelihood to MH ratio in update of 
*     vector of admixture coefficients of individual iindiv
      double precision function contriblikeq(q,qtmp,iindiv,
     &     zcodom,zdom,zhap,codom,dom,hap,
     &     npop,npopmax,nindiv,f,nlocd,nloch,nalmax)
      implicit none
      integer npop,npopmax,nindiv,nlocd,nloch,nalmax,zcodom,zdom,zhap,
     &     iindiv,codom,dom,hap
      double precision f,q,qtmp
      dimension zcodom(nindiv,2*nlocd),zdom(nindiv,nlocd),
     &     zhap(nindiv,nloch),f(npop,nlocd,nalmax),
     &     q(nindiv,npopmax),qtmp(nindiv,npopmax)
      integer iloc,ipop
      double precision c1,c2,c3,c4
      contriblikeq = 0
*     case codominant markers
      if(codom .eq. 1) then
         do iloc=1,nlocd
            c1 = 0
            c2 = 0
            c3 = 0
            c4 = 0
            do ipop=1,npop
               if(zcodom(iindiv,2*iloc-1) .ne. -999) then
                  c1 = c1 +
     &              qtmp(iindiv,ipop)*
     &                 f(ipop,iloc,zcodom(iindiv,2*iloc-1)) 
                  c3 = c3 +
     &                 q(iindiv,ipop)*
     &                 f(ipop,iloc,zcodom(iindiv,2*iloc-1)) 
               endif
               if(zcodom(iindiv,2*iloc) .ne. -999) then
                  c2 = c2 +
     &                 qtmp(iindiv,ipop)*
     &                 f(ipop,iloc,zcodom(iindiv,2*iloc)) 
                  c4 = c4 +
     &                 q(iindiv,ipop)*f(ipop,iloc,zcodom(iindiv,2*iloc)) 
               endif
            enddo
            if(zcodom(iindiv,2*iloc-1) .ne. -999) then
               contriblikeq = contriblikeq + dlog(c1) - dlog(c3)
            endif
            if(zcodom(iindiv,2*iloc) .ne. -999) then
               contriblikeq = contriblikeq + dlog(c2) - dlog(c4) 
            endif
         enddo
      endif
*     case haploid data
      if(hap .eq. 1) then
         do iloc=1,nlocd
            c1 = 0
            c3 = 0
            do ipop=1,npop
               if(zhap(iindiv,iloc) .ne. -999) then
                  c1 = c1 +
     &              qtmp(iindiv,ipop)*
     &                 f(ipop,iloc,zhap(iindiv,iloc)) 
                  c3 = c3 +
     &                 q(iindiv,ipop)*f(ipop,iloc,zhap(iindiv,iloc)) 
               endif
            enddo
            if(zhap(iindiv,iloc) .ne. -999) then
               contriblikeq = contriblikeq + dlog(c1) - dlog(c3)
            endif
         enddo
      endif
*     case dominant markers
      if(dom .eq. 1) then
         do iloc=1,nlocd
            c1 = 0
            c2 = 0
            c3 = 0
            c4 = 0
            if(zhap(iindiv,iloc) .eq. 1) then
               do ipop=1,npop
                  c1 = c1 +
     &                 qtmp(iindiv,ipop)*
     &                 f(ipop,iloc,1) 
                  c3 = c3 +
     &                 q(iindiv,ipop)*f(ipop,iloc,1) 
               enddo
               contriblikeq = contriblikeq + dlog(c1) - dlog(c3)
            endif
            if(zhap(iindiv,iloc) .eq. 2) then
               do ipop=1,npop
                  c1 = c1 +
     &                 qtmp(iindiv,ipop)*f(ipop,iloc,1) 
                  c3 = c3 +
     &                 q(iindiv,ipop)*f(ipop,iloc,1) 
                  c2 = c2 +
     &                 qtmp(iindiv,ipop)* f(ipop,iloc,2) 
                  c4 = c4 +
     &                 q(iindiv,ipop)*f(ipop,iloc,2) 
                  contriblikeq = contriblikeq + 
     &                 dlog(c2**2 +2*c1*c2) - 
     &                 dlog(c4**2 +2*c3*c4)
               enddo
            endif
         enddo
      endif
c      write(*,*) 'contriblikeq = ',contriblikeq
      end function  contriblikeq



************************************************************************
*     log of contribution prior to MH ratio in update of 
*     vector of admixture coefficients of individual iindiv
      double precision function contprq(q,qtmp,iindiv,npop,npopmax,
     &     nindiv,alphadmix)
      implicit none
      integer npop,npopmax,nindiv,iindiv
      double precision q,qtmp,alphadmix
      dimension q(nindiv,npopmax),qtmp(nindiv,npopmax),
     &     alphadmix(nindiv,npopmax)
      integer ipop
      contprq = 0
      do ipop=1,npop
         contprq = contprq + (alphadmix(iindiv,ipop)-1)*
     &        (dlog(qtmp(iindiv,ipop))-dlog(q(iindiv,ipop)))
      enddo
      end function  contprq 


**********************************************************************
*     update a
      subroutine uda(npop,npopmax,nindiv,a,atmp,amax,alphadmix,
     &     alphadmixtmp,q,dist,compar)
      implicit none
      integer npop,npopmax,nindiv,compar
      double precision a,atmp,amax,alphadmix,alphadmixtmp,q,dist
      dimension alphadmix(nindiv,npopmax),alphadmixtmp(nindiv,npopmax),
     &     q(nindiv,npopmax),dist(nindiv,npopmax),a(npopmax),
     &     atmp(npopmax)
      integer ipop,iindiv,k
      double precision ggrnorm,gglgamfn,lratio,s,stmp,ggrbinom,bern
      if(compar .eq. 1) then 
         atmp(1) = a(1) + 0.1*ggrnorm(0.d0,1.d0)
         do ipop = 2,npop
            atmp(ipop) = atmp(1)
         enddo
         lratio = 0
         if((atmp(1) .gt. 1.d-300) .and. (atmp(1) .lt. amax)) then
            do iindiv = 1,nindiv
               s = 0
               stmp = 0
               do ipop = 1,npop
                  s    = s    +  alphadmix(iindiv,ipop)
                  stmp = stmp +  alphadmixtmp(iindiv,ipop)
                  alphadmixtmp(iindiv,ipop)=(atmp(ipop)/a(ipop))*
     &                 alphadmix(iindiv,ipop)
                  lratio =  lratio + 
     &                 (alphadmixtmp(iindiv,ipop)-
     &                 alphadmix(iindiv,ipop))*
     &                 dlog(q(iindiv,ipop)) + 
     &                 gglgamfn(alphadmix(iindiv,ipop)) - 
     &                 gglgamfn(alphadmixtmp(iindiv,ipop))
               enddo
               lratio =  lratio + gglgamfn(stmp) - gglgamfn(s)
            enddo
c            write(*,*) 'lratio=',lratio
            bern = 1
            if(lratio .le. 0) then
               bern = ggrbinom(1.d0,dexp(lratio))
            endif
            if(bern .eq. 1) then 
c               write(*,*) 'bern=',bern
               do ipop = 1,npop
                  a(ipop) = atmp(ipop)
                  do iindiv = 1,nindiv
                     alphadmix(iindiv,ipop) = alphadmixtmp(iindiv,ipop)
                  enddo
               enddo
            endif
         endif
      else
         do k=1,npop
            atmp(k) = a(k) + 0.1*ggrnorm(0.d0,1.d0)
            lratio = 0
            if((atmp(k) .gt. 1.d-300) .and. (atmp(k) .lt. amax)) then
               do iindiv = 1,nindiv
                  s = 0
                  stmp = 0
                  do ipop = 1,npop
                     alphadmixtmp(iindiv,ipop)=(atmp(ipop)/a(ipop))*
     &                    alphadmix(iindiv,ipop)
                     s    = s    +  alphadmix(iindiv,ipop)
                     stmp = stmp +  alphadmixtmp(iindiv,ipop)
                     lratio =  lratio + 
     &                    (alphadmixtmp(iindiv,ipop)-
     &                    alphadmix(iindiv,ipop))*
     &                    dlog(q(iindiv,ipop)) + 
     &                    gglgamfn(alphadmix(iindiv,ipop)) - 
     &                    gglgamfn(alphadmixtmp(iindiv,ipop))
                  enddo
                  lratio =  lratio + gglgamfn(stmp) - gglgamfn(s)
               enddo
c     write(14,*) dexp(lratio)
               bern = 1
               if(lratio .le. 0) then
                  bern = ggrbinom(1.d0,dexp(lratio))
               endif
               if(bern .eq. 1) then 
                  a(k) = atmp(k)
                  do iindiv = 1,nindiv
                     do ipop = 1,npop
                        alphadmix(iindiv,ipop)=alphadmixtmp(iindiv,ipop)
                     enddo
                  enddo
               endif
            endif
         enddo
      endif
      end subroutine uda




**********************************************************************
*     update b
      subroutine udb(npop,npopmax,nindiv,a,b,c,btmp,bmax,alphadmix,
     &     alphadmixtmp,q,dist,compar,deltab)
      implicit none
      integer npop,npopmax,nindiv,compar
      double precision a,b,c,btmp,bmax,alphadmix,alphadmixtmp,q,dist,
     &     deltab
      dimension alphadmix(nindiv,npopmax),alphadmixtmp(nindiv,npopmax),
     &     q(nindiv,npopmax),dist(nindiv,npopmax),a(npopmax),b(npopmax),
     &     c(npopmax),btmp(npopmax)
      integer ipop,iindiv,k
      double precision ggrnorm,gglgamfn,lratio,s,stmp,ggrbinom,bern
      if(compar .eq. 1) then
         btmp(1) = b(1) + deltab*ggrnorm(0.d0,1.d0)
         do ipop = 2,npop
            btmp(ipop) = btmp(1)
         enddo
         s = 0
         lratio = 0
         if((btmp(1) .gt. 1.d-300) .and. (btmp(1) .lt. bmax)) then
            do iindiv = 1,nindiv
               s = 0
               stmp = 0
               do ipop = 1,npop
                  alphadmixtmp(iindiv,ipop) = a(ipop)*
     &                 dexp(-(dist(iindiv,ipop)/btmp(ipop))**c(ipop))
                  s = s +  alphadmix(iindiv,ipop)
                  stmp = stmp +  alphadmixtmp(iindiv,ipop)
                  lratio =  lratio + 
     &                 (alphadmixtmp(iindiv,ipop)-
     &                 alphadmix(iindiv,ipop))*
     &                 dlog(q(iindiv,ipop)) + 
     &                 gglgamfn(alphadmix(iindiv,ipop)) - 
     &                 gglgamfn(alphadmixtmp(iindiv,ipop))
               enddo
               lratio =  lratio + gglgamfn(stmp) - gglgamfn(s) 
            enddo
c     write(15,*) dexp(lratio)
            bern = 1
            if(lratio .le. 0) then
               bern = ggrbinom(1.d0,dexp(lratio))
            endif
            if(bern .eq. 1) then 
               do ipop = 1,npop
                  b(ipop) = btmp(ipop)
                  do iindiv = 1,nindiv
                     alphadmix(iindiv,ipop) = 
     &                    alphadmixtmp(iindiv,ipop)
                  enddo
               enddo
            endif
         endif
      else
         do k = 1,npop
            btmp(k) = b(k) + 0.1*ggrnorm(0.d0,1.d0)
            s = 0
            lratio = 0
            if((btmp(k) .gt. 1.d-300) .and. (btmp(k) .lt. bmax)) then
               do iindiv = 1,nindiv
                  s = 0
                  stmp = 0
                  do ipop = 1,npop
                     alphadmixtmp(iindiv,ipop) = a(ipop)*
     &                    dexp(-(dist(iindiv,ipop)/btmp(ipop))**c(ipop))
                     s = s +  alphadmix(iindiv,ipop)
                     stmp = stmp +  alphadmixtmp(iindiv,ipop)
                     lratio =  lratio + 
     &                    (alphadmixtmp(iindiv,ipop)-
     &                    alphadmix(iindiv,ipop))*
     &                    dlog(q(iindiv,ipop)) + 
     &                    gglgamfn(alphadmix(iindiv,ipop)) - 
     &                    gglgamfn(alphadmixtmp(iindiv,ipop))
                  enddo
                  lratio =  lratio + gglgamfn(stmp) - gglgamfn(s) 
               enddo
c     write(15,*) dexp(lratio)
               bern = 1
               if(lratio .le. 0) then
                  bern = ggrbinom(1.d0,dexp(lratio))
               endif
               if(bern .eq. 1) then 
                  b(k) = btmp(k)
                  do iindiv = 1,nindiv
                     do ipop = 1,npop
                        alphadmix(iindiv,ipop) = 
     &                       alphadmixtmp(iindiv,ipop)
                     enddo
                  enddo
               endif
            endif
         enddo
      endif
      end subroutine udb




**********************************************************************
*     update c
      subroutine udc(npop,npopmax,nindiv,a,b,c,ctmp,cmax,alphadmix,
     &     alphadmixtmp,q,dist,compar)
      implicit none
      integer npop,npopmax,nindiv,compar
      double precision a,b,c,ctmp,cmax,alphadmix,alphadmixtmp,q,dist
      dimension alphadmix(nindiv,npopmax),alphadmixtmp(nindiv,npopmax),
     &     q(nindiv,npopmax),dist(nindiv,npopmax),a(npopmax),b(npopmax),
     &     c(npopmax),ctmp(npopmax)
      integer ipop,iindiv,k
      double precision ggrnorm,gglgamfn,lratio,s,stmp,ggrbinom,bern
      if(compar .eq. 1) then
         ctmp(1) = c(1) + 0.1*ggrnorm(0.d0,1.d0)
         do ipop = 2,npop
            ctmp(ipop) = ctmp(1)
         enddo
         s = 0
         lratio = 0
         if((ctmp(1) .gt. 1.d-300) .and. (ctmp(1) .lt. cmax)) then
            do iindiv = 1,nindiv
               s = 0
               stmp = 0
               do ipop = 1,npop
                  alphadmixtmp(iindiv,ipop) = a(ipop)*
     &                 dexp(-(dist(iindiv,ipop)/b(ipop))**ctmp(ipop))
                  s = s +  alphadmix(iindiv,ipop)
                  stmp = stmp +  alphadmixtmp(iindiv,ipop)
                  lratio =  lratio + 
     &                 (alphadmixtmp(iindiv,ipop)-
     &                 alphadmix(iindiv,ipop))*
     &                 dlog(q(iindiv,ipop)) + 
     &                 gglgamfn(alphadmix(iindiv,ipop)) - 
     &                 gglgamfn(alphadmixtmp(iindiv,ipop))
               enddo
               lratio =  lratio + gglgamfn(stmp) - gglgamfn(s) 
            enddo
            bern = 1
            if(lratio .le. 0) then
               bern = ggrbinom(1.d0,dexp(lratio))
            endif
            if(bern .eq. 1) then 
               do ipop = 1,npop
                  c(ipop) = ctmp(ipop)
                  do iindiv = 1,nindiv
                     alphadmix(iindiv,ipop)=alphadmixtmp(iindiv,ipop)
                  enddo
               enddo
            endif
         endif
         do k = 2,npop
            c(k) = c(1)
         enddo
      else
         do k = 1,npop
            ctmp(k) = c(k) + 0.1*ggrnorm(0.d0,1.d0)
            s = 0
            lratio = 0
            if((ctmp(k) .gt. 1.d-300) .and. (ctmp(k) .lt. cmax)) then
               do iindiv = 1,nindiv
                  s = 0
                  stmp = 0
                  do ipop = 1,npop
                     alphadmixtmp(iindiv,ipop) = a(ipop)*
     &                    dexp(-(dist(iindiv,ipop)/b(ipop))**ctmp(ipop))
                     s = s +  alphadmix(iindiv,ipop)
                     stmp = stmp +  alphadmixtmp(iindiv,ipop)
                     lratio =  lratio + 
     &                    (alphadmixtmp(iindiv,ipop)-
     &                    alphadmix(iindiv,ipop))*
     &                    dlog(q(iindiv,ipop)) + 
     &                    gglgamfn(alphadmix(iindiv,ipop)) - 
     &                    gglgamfn(alphadmixtmp(iindiv,ipop))
                  enddo
                  lratio =  lratio + gglgamfn(stmp) - gglgamfn(s) 
               enddo
               bern = 1
               if(lratio .le. 0) then
                  bern = ggrbinom(1.d0,dexp(lratio))
               endif
               if(bern .eq. 1) then 
                  c(k) = ctmp(k)
                  do iindiv = 1,nindiv
                     do ipop = 1,npop
                        alphadmix(iindiv,ipop)=alphadmixtmp(iindiv,ipop)
                     enddo
                  enddo
               endif
            endif
         enddo
      endif
      end subroutine udc
