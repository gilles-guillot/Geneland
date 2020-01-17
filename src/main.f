      subroutine mcmcgld(s,zz,missloc,z,ql,nql,qtc,nqtc,
     &     path,intpar,dblepar,
     &     nindiv,nlocd,nloch,ncolt,nal,nalmax,nppmax,
     &     npop,npopmin,npopmax,xlim,ylim, indcell,indcelltmp,
     &     distcell,distcelltmp,
     &     t,ttmp,u,utmp,c,ctmp,f,ftmp,fa,drift,drifttmp,
     &     n,ntmp,a,ptmp,
     &     meanqtc,sdqtc,meanqtctmp,sdqtctmp,
     &     nnqtc,sqtc,ssqtc,
     &     ksiqtc,kappaqtc,alphaqtc,betaqtc,gbeta,hbeta,
     &     cellpop,listcell,
     &     yy,fcy,nitsaved,out1,outspace,outfreq,outqtc) 
      implicit none 

*     data
      integer nindiv,nlocd,nlocd2,nloch,nloch2,ncolt,
     &     nal,nalmax,zz,z,ploidy,jcf,nchpath,missloc,
     &     nqtc,ql,nql,nitsaved
      double precision s,qtc

*     hyper parameters
      double precision lambdamax,dt,shape1,shape2,
     &     ksiqtc,kappaqtc,alphaqtc,betaqtc,gbeta,hbeta

*     parameters 
      integer npp,nppmax,npop,npopmin,npopmax,c,ctmp
      double precision lambda,u,utmp,f,t,fa,drift,ftmp,drifttmp,alpha,
     &     meanqtc,sdqtc,meanqtctmp,sdqtctmp,out1,outspace,outfreq,
     &     outqtc

*     modeling/computing options
      integer nit,thinning,fmodel,kfix,spatial,filtna,nudcel,
     &     usegeno1,usegeno2,useqtc,useqtd,useql,intpar
      double precision du,pudcel,dblepar

*     variables de travail
      integer iloc,iindiv,iit,ipp,ipop,ial,indcell,indcelltmp,
     &     n,cellpop,listcell,cellpophost,ntmp,nn,yy,iud,nud,nnqtc,iqtc,
     &     iitstor
      double precision ptmp,xlim,ylim,ggrunif,rpostlamb,
     &     distcell,distcelltmp,a,ttmp,lpriorallvar,llallvar2,
     &     lpriorallvartmp,llallvartmp,fcy,pct,sqtc,ssqtc
      character*255  path,filef,filenpp,filelambda,filenpop,fileu,filec,
     &     filefa,filedrift,filelpp,filell,filet,filesize,
     &     filemq,filesdq,filebetaqtc

*     dimensions
      dimension s(2,nindiv),t(2,nindiv),zz(nindiv,2*nlocd),
     &     z(nindiv,nloch),ql(nindiv,nql),qtc(nindiv,nqtc),
     &     u(2,nppmax),utmp(2,nppmax),
     &     c(nppmax),ctmp(nppmax), f(npopmax,ncolt,nalmax),
     &     xlim(2),ylim(2),nal(ncolt),
     &     indcell(nindiv),indcelltmp(nindiv),
     &     distcell(nindiv),ttmp(2,nindiv),
     &     distcelltmp(nindiv),n(npopmax,ncolt,nalmax),
     &     ntmp(npopmax,ncolt,nalmax),
     &     a(nalmax),ptmp(nalmax),
     &     fa(ncolt,nalmax),drift(npopmax),
     &     ftmp(npopmax,ncolt,nalmax),drifttmp(npopmax),
     &     cellpop(nppmax),listcell(nppmax),cellpophost(nppmax),
     &     yy(nindiv,2*nlocd+2*nloch),fcy(nalmax,2),
     &     missloc(nindiv,nlocd),
     &     meanqtc(npopmax,nqtc),sdqtc(npopmax,nqtc),
     &     meanqtctmp(npopmax,nqtc),sdqtctmp(npopmax,nqtc),
     &     nnqtc(npopmax,nqtc),sqtc(npopmax,nqtc),
     &     ssqtc(npopmax,nqtc),intpar(100),dblepar(100),
     &     ksiqtc(nqtc),kappaqtc(nqtc),alphaqtc(nqtc),betaqtc(nqtc),
     &     gbeta(nqtc),hbeta(nqtc),
     &     out1(nitsaved,5),outspace(nitsaved,5,nppmax),
     &     outfreq(nitsaved,3,npopmax,ncolt,nalmax),
     &     outqtc(nitsaved,2,npopmax,nqtc)


 1000 format (2(1x,e15.8,1x))
 2000 format (300(1x,e15.8,1x))
 3000 format (300(1x,i10,1x))

      call intpr('***************************************',-1,0,0)
      call intpr('***    Starting MCMC simulation     ***',-1,0,0)
      call intpr('***************************************',-1,0,0)

*     init RNG
      call rndstart()
c$$$      write(*,*) "apres rndstart ggrunif=",ggrunif(0.d0,1.d0) 
c$$$      write(*,*) 'nlocd=',nlocd
c$$$      write(*,*) 'nloch=',nloch
c$$$      write(*,*) ''

*     unwrap computing options
      fmodel =     intpar(15+1)
      kfix =       intpar(15+2)
      spatial =    intpar(15+3)
      jcf =        intpar(15+4)
      filtna =     intpar(15+5)
c     parameter ploidy says how to interpret data in matrix z
      ploidy =     intpar(15+6)
      nchpath =    intpar(15+7)
      nit =        intpar(15+8)
      thinning =   intpar(15+9)
      usegeno1 =   intpar(15+10)
      usegeno2 =   intpar(15+11)
      useqtc =     intpar(15+12)
      useqtd =     intpar(15+13)
      useql =      intpar(15+14)

      lambdamax = dblepar(1)
      dt =        dblepar(2)
      shape1 =    dblepar(3)
      shape2 =    dblepar(4)
      pudcel =    dblepar(5)

*     look for smallest rectangle with edges parrallel to axes enclosing the spatial domain
      call limit(nindiv,s,xlim,ylim,dt)

*     Ouverture des fichiers pour l'ecriture des sorties
      filelambda = path(1:nchpath) // "Poisson.process.rate.txt"
      filenpp = path(1:nchpath) // "nuclei.numbers.txt"
      filenpop = path(1:nchpath) // "populations.numbers.txt"
      fileu = path(1:nchpath) // "coord.nuclei.txt"
      filec = path(1:nchpath) // "color.nuclei.txt"
      filef = path(1:nchpath) // "frequencies.txt"
      filefa = path(1:nchpath) // "ancestral.frequencies.txt"
      filedrift = path(1:nchpath) // "drifts.txt"
      filelpp = path(1:nchpath) // "log.posterior.density.txt"
      filell = path(1:nchpath) // "log.likelihood.txt"
      filet = path(1:nchpath) // "hidden.coord.txt"
      filesize = path(1:nchpath) // "size.pop.txt"
      filemq  = path(1:nchpath) // "mean.qtc.txt"
      filesdq  = path(1:nchpath) // "sd.qtc.txt"
      filebetaqtc = path(1:nchpath) // "beta.qtc.txt"


c$$$      if(intpar(1) .eq.1) then
c$$$         open(9,file=filelambda)
c$$$      endif
c$$$      if(intpar(2) .eq.1) then
c$$$         open(10,file=filenpp)
c$$$      endif
c$$$      if(intpar(3) .eq.1) then
c$$$         open(11,file=filenpop)
c$$$      endif   
c$$$      if(intpar(4) .eq.1) then
c$$$         open(12,file=fileu)
c$$$      endif
c$$$      if(intpar(5) .eq.1) then
c$$$         open(13,file=filec)
c$$$      endif
c$$$      if(((usegeno2 .eq. 1) .or. (usegeno1 .eq. 1)) .or. 
c$$$     &     (useql .eq. 1)) then  
c$$$         if(intpar(6) .eq.1) then
c$$$            open(14,file=filef)
c$$$         endif
c$$$         if(intpar(7) .eq.1) then
c$$$            open(15,file=filefa) 
c$$$         endif
c$$$         if(intpar(8) .eq.1) then
c$$$            open(16,file=filedrift)
c$$$         endif
c$$$      endif
c$$$      if(intpar(9) .eq.1) then
c$$$         open(17,file=filelpp)
c$$$      endif
c$$$      if(intpar(10) .eq.1) then
c$$$         open(18,file=filell) 
c$$$      endif
c$$$      if(intpar(11) .eq.1) then
c$$$         open(19,file=filet) 
c$$$      endif
c$$$      if(intpar(12) .eq.1) then
c$$$         open(20,file=filesize)
c$$$      endif
c$$$      if(useqtc .eq. 1) then 
c$$$         if(intpar(13) .eq.1) then
c$$$            open(21,file=filemq)
c$$$         endif
c$$$         if(intpar(14) .eq.1) then
c$$$            open(22,file=filesdq)
c$$$         endif   
c$$$         if(intpar(15) .eq.1) then
c$$$            open(23,file=filebetaqtc)
c$$$         endif
c$$$      endif

      call intpr('Output files have been opened:',-1,0,0)

 
************************
*     Initialization
************************
*     parameter for the Dirichlet model for allele freq. 
*     not used if (fmodel .eq. 1)
      alpha = 1


      du = dsqrt((xlim(2)-xlim(1))*(ylim(2)-ylim(1))/dble(nindiv))
      lambda = lambdamax*ggrunif(0.d0,1.d0) 
      if(spatial .eq. 1) then
         npp = 1 + idint(dint(lambda))
      else
         npp = nindiv
      endif
      call intpr('npp initialised:',-1,0,0)
      if(spatial .eq. 1) then 
         call rprioru(npp,nppmax,xlim,ylim,u)
      else 
         do iindiv=1,nindiv
            u(1,iindiv) = s(1,iindiv)
            u(2,iindiv) = s(2,iindiv) 
         enddo
      endif 
      call intpr('u initialised:',-1,0,0)
      call rpriorc(npp,nppmax,npop,c)
      call intpr('c initialised:',-1,0,0)


c$$$      npp = 2
c$$$      c(1) = 1
c$$$      c(2) = 1
c$$$      u(1,1) = .2
c$$$      u(2,1) = .5
c$$$      u(1,2) = .8
c$$$      u(2,2) = .5
c$$$      npp = 4
c$$$      c(1) = 1
c$$$      c(2) = 2
c$$$      c(3) = 3
c$$$      c(4) = 4
c$$$      u(1,1) = .2
c$$$      u(2,1) = .2
c$$$      u(1,2) = .8
c$$$      u(2,2) = .2
c$$$      u(1,3) = .2
c$$$      u(2,3) = .8
c$$$      u(1,4) = .8
c$$$      u(2,4) = .8
c$$$      npp = 16
c$$$      c(1) = 1
c$$$      c(2) = 1
c$$$      c(3) = 1
c$$$      c(4) = 1
c$$$      c(5) = 1
c$$$      c(6) = 1
c$$$      c(7) = 1
c$$$      c(8) = 1
c$$$      c(9) = 1
c$$$      c(10) = 1
c$$$      c(11) = 1
c$$$      c(12) = 1
c$$$      c(13) = 1
c$$$      c(14) = 1
c$$$      c(15) = 1
c$$$      c(16) = 1
c$$$      u(1,1) = .2
c$$$      u(2,1) = .2
c$$$      u(1,2) = .2
c$$$      u(2,2) = .4
c$$$      u(1,3) = .2
c$$$      u(2,3) = .6
c$$$      u(1,4) = .2
c$$$      u(2,4) = .8
c$$$      u(1,5) = .4
c$$$      u(2,5) = .2
c$$$      u(1,6) = .4
c$$$      u(2,6) = .4
c$$$      u(1,7) = .4
c$$$      u(2,7) = .6
c$$$      u(1,8) = .4
c$$$      u(2,8) = .8
c$$$      u(1,9) = .6
c$$$      u(2,9) = .2
c$$$      u(1,10) = .6
c$$$      u(2,10) = .4
c$$$      u(1,11) = .6
c$$$      u(2,11) = .6
c$$$      u(1,12) = .6
c$$$      u(2,12) = .8
c$$$      u(1,13) = .8
c$$$      u(2,13) = .2
c$$$      u(1,14) = .8
c$$$      u(2,14) = .4
c$$$      u(1,15) = .8
c$$$      u(2,15) = .6
c$$$      u(1,16) = .8
c$$$      u(2,16) = .8

      do iindiv=1,nindiv
         t(1,iindiv) = s(1,iindiv) 
         t(2,iindiv) = s(2,iindiv) 
      enddo
      call intpr('t initialised:',-1,0,0)

      call calccell(nindiv,t,npp,nppmax,u,indcell,distcell)
      call intpr('Voronoi cells initialised:',-1,0,0)

      call rpriordrift(npop,npopmax,drift,fmodel,shape1,shape2)
      call intpr('drift initialised:',-1,0,0)

      call rpriorfaallvar(nlocd,nloch,nql,ncolt,nal,nalmax,fa,
     &     fmodel,ptmp,usegeno2,usegeno1,useql)
      call intpr('fa initialised:',-1,0,0)

      call rpostf4(yy,z,ql,nindiv,nlocd,nloch,nql,ncolt,
     &     nal,nalmax,npop,npopmax,nppmax,c,indcell,n,a,ptmp,
     &     f,fa,drift,usegeno1,usegeno2,useql,ploidy)
      call intpr('f initialised:',-1,0,0)

      if(useqtc .eq.1) then
         call rpostqtc2(qtc,nindiv,nqtc,indcell,c,meanqtc,
     &        sdqtc,nnqtc,sqtc,ssqtc,npop,npopmax,nppmax,
     &        ksiqtc,kappaqtc,alphaqtc,betaqtc)
          call intpr('mean & var quantitative variables initialized:',
     &        -1,0,0)
      endif

      call intpr('All parameters have been initialised:',-1,0,0)
*     end of init
******************************************************     


  
c$$$
c$$$      write(*,*) 'nindiv =', nindiv
c$$$      write(*,*) 'nlocd =',nlocd   
c$$$      write(*,*) 'nloch =',nloch   
c$$$      write(*,*) 'nql =',nql   
c$$$      write(*,*) 'ncolt =',ncolt   
c$$$      write(*,*) 'nal =',nal 
c$$$      write(*,*) 'nalmax =', nalmax
c$$$      write(*,*) 'usegeno2=',usegeno2
c$$$      write(*,*) 'usegeno1=',usegeno1
c$$$      write(*,*) 'useql=',useql
c$$$      write(*,*) 'fmodel=',fmodel
c$$$      write(*,*) 'kfix =',kfix
c$$$      write(*,*) 'spatial =',spatial
c$$$      write(*,*) 'jcf =',jcf
c$$$      write(*,*) 'ploidy=',ploidy
c$$$      write(*,*) 'filtna=',filtna
c$$$      write(*,*) 'dt =', dt
c$$$      write(*,*) 'lambda=',lambda
c$$$      write(*,*) 'lambdamax=',lambdamax  
c$$$      write(*,*) 'npp =', npp 
c$$$      write(*,*) 'nppmax =', nppmax     
c$$$      write(*,*) 'npopmin =',npopmin  
c$$$      write(*,*) 'npopmax =', npopmax  
c$$$      write(*,*) 'nit =', nit  
c$$$      write(*,*) 'thinning =',thinning  
c$$$      write(*,*) 'xlim=', xlim 
c$$$      write(*,*) 'ylim=', ylim 
c$$$      write(*,*) 'zz =', zz
c$$$      write(*,*) 'z =', z
c$$$      write(*,*) 'yy =', yy
c$$$      write(*,*) 'ql =', ql
c$$$      write(*,*) 's =', s
c$$$      write(*,*) 't =', t
c$$$      write(*,*) 'c =', c
c$$$      write(*,*) 'ctmp =',ctmp 
c$$$      write(*,*) 'u =',u  
c$$$      write(*,*) 'utmp =',utmp 
c$$$      write(*,*) 'f =', f  
c$$$      write(*,*) 'ftmp =', ftmp
c$$$      write(*,*) 'indcell =', indcell 
c$$$      write(*,*) 'distcell =', distcell
c$$$      write(*,*) 'indcelltmp =', indcelltmp 
c$$$      write(*,*) 'distcelltmp =',distcelltmp 
c$$$      write(*,*) 'intpar=',intpar
c$$$      write(*,*) 'qtc=',qtc
c$$$      write(*,*) 'meanqtc=',meanqtc
c$$$      write(*,*) 'sdqtc=',sdqtc
c$$$      write(*,*) 'nnqtc=',nnqtc
c$$$      write(*,*) 'intpar=',intpar
c$$$      write(*,*) 'dblepar=',dblepar
c$$$      write(*,*) 'pudcel=',pudcel
c$$$      write(*,*) 'missloc=',missloc
c$$$      write(*,*) 'meanqtc=',sdqtc
c$$$      write(*,*) 'meanqtctmp=',sdqtctmp
c$$$      write(*,*) 'sdqtc=',sdqtc
c$$$      write(*,*) 'sdqtctmp=',sdqtctmp
c$$$      write(*,*) 'ksiqtc=',ksiqtc
c$$$      write(*,*) 'kappaqtc=',kappaqtc
c$$$      write(*,*) 'alphaqtc=',alphaqtc
c$$$      write(*,*) 'betaqtc=',betaqtc
 

*****************************************  
* Main MCMC loop
*****************************************  
      call intpr('Percentage of computations:',-1,0,0)
      iitstor = 0
      do iit=1,nit
         if(mod(iit,thinning) .eq. 0) then
            iitstor = iitstor + 1
*************************************
*     writing results in output files


******************************
*     writing scalar variables 
            pct = dble(iit)/dble(nit)*100.
            call dblepr('                     ',-1,pct,1)
            if(intpar(1) .eq.1) then
c               write(9,*) lambda
               out1(iitstor,1) = lambda
            endif
            if(intpar(2) .eq.1) then
c               write(10,*) npp
               out1(iitstor,2) = dble(npp)
            endif
            if(intpar(3) .eq.1) then
c               write(11,*) npop
               out1(iitstor,3) = dble(npop)
            endif
            if(intpar(10) .eq.1) then
               llallvartmp = llallvar2(yy,z,ql,nindiv,nlocd,nloch,
     &              nql,ncolt,npopmax,nalmax,nppmax,c,f,indcell,qtc,
     &              nqtc,meanqtc,sdqtc,
     &              usegeno2,usegeno1,useql,useqtc,ploidy)
c               write(18,*) llallvartmp
               out1(iitstor,4) = llallvartmp
            endif
            if(intpar(9) .eq.1) then
               lpriorallvartmp = llallvartmp + 
     &              lpriorallvar(lambdamax,lambda,
     &     npop,npp,drift,f,fa,c,nppmax,nindiv,npopmax,nlocd,nloch,nql,
     &     ncolt,nal,nalmax,indcell,fmodel,xlim,ylim,shape1,shape2,
     &     nqtc,meanqtc,sdqtc,ksiqtc,kappaqtc,alphaqtc,betaqtc,
     &     usegeno2,usegeno1,useql,useqtc)
c                write(17,*) lpriorallvartmp + llallvartmp
c     line above seens wrongs corrected as below on 2012/09/06
               out1(iitstor,5) = lpriorallvartmp 
            endif
*     end writing scalar variables 
**********************************             


***************************************
*     writing space variables 
            if(intpar(4) .eq.1) then
c$$$               write(12,1000) 
c$$$     &              (sngl(u(1,ipp)),sngl(u(2,ipp)), ipp=1,nppmax)
               do ipp = 1,nppmax
                  outspace(iitstor,1,ipp) = u(1,ipp)
                  outspace(iitstor,2,ipp) = u(2,ipp)
               enddo
            endif
            if(intpar(5) .eq.1) then
               do ipp=1,nppmax
c                  write(13,*) c(ipp)
                  outspace(iitstor,3,ipp) = dble(c(ipp))
               enddo
            endif
            if(intpar(11) .eq. 1) then
c$$$               write(19,1000) (sngl(t(1,iindiv)),sngl(t(2,iindiv)),
c$$$     &              iindiv=1,nindiv)
               do iindiv = 1,nindiv
                  outspace(iitstor,4,iindiv) = t(1,iindiv)
                  outspace(iitstor,5,iindiv) = t(2,iindiv)
               enddo
            endif
***************************************          

***************************************
*     writing allele frequencies variables 
            if(((usegeno2 .eq. 1) .or. (usegeno1 .eq. 1)) .or. 
     &           (useql .eq. 1)) then  
               if(intpar(6) .eq.1) then
                  do iloc=1,ncolt
                     do ial=1,nalmax
c$$$  write(14,2000) (sngl(f(ipop,iloc,ial)),
c$$$  &                       ipop=1,npopmax)
                        do ipop = 1,npop
                           outfreq(iitstor,1,ipop,iloc,ial) = 
     &                          f(ipop,iloc,ial)
                        enddo
                     enddo
                  enddo  
               endif
               
               if(intpar(7) .eq.1) then
c$$$  write(15,2000) 
c$$$  &                 ((sngl(fa(iloc,ial)),ial=1,nalmax),iloc=1,ncolt)
                  do iloc=1,ncolt
                     do ial=1,nalmax
                        outfreq(iitstor,2,1,iloc,ial) = fa(iloc,ial)
                     enddo
                  enddo
               endif
               if(intpar(8) .eq.1) then
c     write(16,2000) (sngl(drift(ipop)),ipop=1,npopmax)
                  do ipop = 1,npop
                     outfreq(iitstor,3,ipop,1,1) = drift(ipop)
                  enddo
               endif
            endif
*****************************************


*     counting nb of individuals in each pop
c$$$            if(intpar(12) .eq.1) then
c$$$               do ipop = 1,npopmax
c$$$                  n(ipop,1,1) = 0
c$$$               enddo
c$$$               do iindiv = 1,nindiv
c$$$                  n(c(indcell(iindiv)),1,1) =  
c$$$     &                 n(c(indcell(iindiv)),1,1) + 1 
c$$$               enddo
c$$$               write(20,3000) (n(ipop,1,1),ipop=1,npopmax)
c$$$            endif

***************************** 
* writing qtc variables
            if(useqtc .eq. 1) then 
               if(intpar(13) .eq.1) then
                  do ipop = 1,npopmax
c$$$  write(21,2000)  
c$$$  &                    (sngl(meanqtc(ipop,iqtc)),iqtc=1,nqtc)
                     do iqtc = 1,nqtc
                        outqtc(iitstor,1,ipop,iqtc) = meanqtc(ipop,iqtc)
                     enddo
                  enddo
               endif
               if(intpar(14) .eq. 1) then
                  do ipop = 1,npopmax
c$$$  write(22,2000)  
c$$$  &                    (sngl(sdqtc(ipop,iqtc)),iqtc=1,nqtc)
                     do iqtc = 1,nqtc
                        outqtc(iitstor,2,ipop,iqtc) = sdqtc(ipop,iqtc)
                     enddo
                  enddo
               endif
            endif
         endif
*     end writing
*****************************



**************************************
*     defining nb of cells updated (depends on current value of npp)
         nudcel = max0(1,min0(idint(dint(dble(npp)*pudcel)),npp))
c         write(*,*) 'nudcel =',nudcel

**************************************
*     update lambda
c          write(*,*) 'update lambda'
          if(spatial .eq. 1) then
             lambda = rpostlamb(lambdamax,npp)
          endif

          if((fmodel .eq. 1) .and. (((usegeno2 .eq. 1) .or. 
     &         (usegeno1 .eq. 1)) .or. (useql .eq. 1))) then 
**************************************
*     update drift
c              write(*,*) 'update drift'
             call  upddriftallvar(npop,npopmax,nlocd,nloch,nql,ncolt,
     &            nalmax,nal,f,fa,drift,shape1,shape2,
     &            usegeno2,usegeno1,useql)
c              write(*,*) "apres  upddriftallvar=",ggrunif(0.d0,1.d0) 
**************************************
*     update fa
c              write(*,*) 'update fa'
             call updfa2(npop,npopmax,nlocd,nloch,nql,ncolt,nalmax,nal,
     &     f,fa,drift,usegeno2,usegeno1,useql)
c              write(*,*) "apres updfa=",ggrunif(0.d0,1.d0) 
          endif

***************************************
*     update betaqtc
          if(useqtc .eq. 1) then
             call rpostbetaqtc(nqtc,npop,npopmax, 
     &     alphaqtc,betaqtc,gbeta,hbeta,sdqtc)
          endif
     

************************************** 
*     update c and population-specific parameters 
          if(jcf .eq. 1) then 
c             write(*,*) 'update c and f' 
             if(fmodel .eq. 0) then
*     joint update of c and population-specific parameters 
c$$$                write(*,*) ''
c$$$                write(*,*) 'main loop avant udcfallvar3'
c$$$                write(*,*) 'meanqtc=',meanqtc
c$$$                write(*,*) 'sdqtc=',sdqtc
                call udcfallvar3(npop,npopmax,nlocd,nloch,nql,
     &     ncolt,nal,nalmax,indcell,nindiv,npp,nppmax,c,ctmp,nqtc,qtc,
     &     meanqtc,sdqtc,meanqtctmp,sdqtctmp,nnqtc,sqtc,ssqtc,
     &     ksiqtc,kappaqtc,alphaqtc,betaqtc,a,ptmp,f,ftmp,yy,z,ql,n,
     &     ntmp,alpha,nudcel,usegeno2,usegeno1,useql,useqtc,ploidy)
c$$$               write(*,*) 'meanqtc=',meanqtc
c$$$                write(*,*) 'sdqtc=',sdqtc
c$$$                write(*,*) 'main loop apres udcfallvar3'
             else
                call udcfallvarcfm3(npop,npopmax,f,fa,drift,
     &               nlocd,nloch,nql,ncolt,
     &               nal,nalmax,indcell,nindiv,npp,nppmax,c,ctmp,nqtc,
     &               qtc,meanqtc,sdqtc,meanqtctmp,sdqtctmp,nnqtc,
     &               sqtc,ssqtc,ksiqtc,kappaqtc,alphaqtc,betaqtc,
     &               a,ptmp,ftmp,yy,z,ql,n,ntmp,nudcel,
     &               usegeno2,usegeno1,useql,useqtc,ploidy)
             endif
c             write(*,*) "apres c,f=",ggrunif(0.d0,1.d0) 
          else    
*     separate update of c and population-specific parameters     
             if(((usegeno2 .eq. 1)  .or. 
     &           (usegeno1 .eq. 1)) .or. (useql .eq. 1)) then
                call rpostf4(yy,z,ql,nindiv,nlocd,nloch,nql,
     &               ncolt,nal,nalmax,npop,npopmax,nppmax,c,indcell,
     &               n,a,ptmp,f,fa,drift,usegeno1,usegeno2,useql,
     &               ploidy)
             endif
             if(useqtc .eq. 1) then
                call rpostqtc2(qtc,nindiv,nqtc,indcell,c,meanqtc,
     &        sdqtc,nnqtc,sqtc,ssqtc,npop,npopmax,nppmax,
     &        ksiqtc,kappaqtc,alphaqtc,betaqtc)
             endif
             call udcallvar2(npp,nppmax,c,ctmp,yy,z,ql,nindiv,nlocd,
     &            nloch,nql,ncolt,nalmax,npop,npopmax,f,indcell,nudcel,
     &            qtc,nqtc,meanqtc,sdqtc,usegeno2,usegeno1,useql,
     &            useqtc,ploidy)
          endif

          
          if(spatial .eq. 1) then
**************************************
*     update u et mise a jour de indcell et distcell
c             write(*,*) 'update u'
             call uduallvar2(npp,nppmax,c,u,yy,z,ql,nindiv,nlocd,
     &            nloch,nql,ncolt,nalmax,npopmax,f,indcell,distcell,
     &            indcelltmp,distcelltmp,t,xlim,ylim,du,nudcel,
     &            qtc,nqtc,meanqtc,sdqtc,usegeno2,usegeno1,
     &            useql,useqtc,ploidy)
**************************************             
*     update t et mise a jour de indcell et distcell
             if(dt .gt. 1.d-300) then 
                call udtallvar2(npp,nppmax,nindiv,
     &               nlocd,nloch,nql,ncolt,nqtc,nalmax,npopmax,t,ttmp,
     &               dt,s,c,indcell,distcell,indcelltmp,distcelltmp,
     &               u,yy,z,ql,qtc,f,meanqtc,sdqtc,
     &               usegeno2,usegeno1,useql,useqtc,ploidy)
             endif
**************************************
*     birth/death des points du pp
c             write(*,*) 'update npp'
             call bdcellallvar2(nindiv,u,c,utmp,ctmp,npop,npopmax,
     &            nlocd,nloch,nql,ncolt,nalmax,npp,nppmax,yy,z,
     &            ql,f,t,xlim,ylim,indcell,distcell,indcelltmp,
     &            distcelltmp,lambda,qtc,nqtc,meanqtc,sdqtc,
     &            usegeno2,usegeno1,useql,useqtc,ploidy)
          endif

**************************************
*     birth/death de pop
          if(kfix .eq. 0) then 
             if(dble(iit)/dble(nit) .ge. 0.1) then 
c      write(*,*) 'update npop'
                if(fmodel .eq. 0) then 
c$$$                     write(*,*) 'main loop before smdallvar3'
c$$$                     write(*,*) 'meanqtc=',meanqtc
c$$$                     write(*,*) 'sdqtc=',sdqtc
                   call smdallvar3(npop,npopmin,npopmax,f,fa,drift,
     &                  nlocd,nloch,ncolt,nql,
     &                  nal,nalmax,indcell,nindiv,npp,nppmax,c,ctmp,
     &                  a,ptmp,ftmp,drifttmp,yy,z,ql,cellpop,listcell,
     &                  cellpophost,n,ntmp,qtc,nqtc,meanqtc,sdqtc,
     &                  meanqtctmp,sdqtctmp,nnqtc,sqtc,ssqtc,
     &                  ksiqtc,kappaqtc,alphaqtc,betaqtc,
     &                  usegeno2,usegeno1,useql,useqtc,ploidy)
                else
                   call smfallvar3(npop,npopmin,npopmax,f,fa,drift,
     &                  nlocd,nloch,ncolt,nql,
     &                  nal,nalmax,indcell,nindiv,npp,nppmax,c,ctmp,
     &                  a,ptmp,ftmp,drifttmp,yy,z,ql,cellpop,listcell,
     &                  cellpophost,n,ntmp,qtc,nqtc,meanqtc,sdqtc,
     &                  meanqtctmp,sdqtctmp,nnqtc,sqtc,ssqtc,
     &                  ksiqtc,kappaqtc,alphaqtc,betaqtc,shape1,shape2,
     &                  usegeno2,usegeno1,useql,useqtc,ploidy)
                endif
             endif
           endif

**************************************
*     update true unobserved genotypes
*     in case null alleles are suspected
          if(usegeno2 .eq. 1) then 
             if(filtna .eq. 1) then
c                write(*,*) 'udNA3'
                call udNA3(nindiv,nlocd,nloch,nal,nalmax,ncolt,nppmax,
     &               yy,zz,c,indcell,npopmax,f,fcy,npop,missloc) 
             endif    
          endif

**************************************
*     update true unobserved genotypes
*     for dominant markers data 
          if((usegeno1 .eq. 1) .and. (ploidy .eq. 2))then 
c             write(*,*) 'udyDOM3'
             call udyDOM3(nindiv,nlocd,nloch,nal,nalmax,
     &            ncolt,nppmax,yy,z,c,indcell,npopmax,f,fcy,npop)
          endif


***************************

      enddo
*     end of main mcmc loop
***************************



***********************
*     closing all files
c$$$      if(intpar(1) .eq.1) then
c$$$         close(9)
c$$$      endif
c$$$      if(intpar(2) .eq.1) then
c$$$         close(10)
c$$$      endif
c$$$      if(intpar(3) .eq.1) then
c$$$         close(11)
c$$$      endif   
c$$$      if(intpar(4) .eq.1) then
c$$$         close(12)
c$$$      endif
c$$$      if(intpar(5) .eq.1) then
c$$$         close(13)
c$$$      endif
c$$$      if(((usegeno2 .eq. 1) .or. (usegeno1 .eq. 1)) .or. 
c$$$     &     (useql .eq. 1)) then  
c$$$         if(intpar(6) .eq.1) then
c$$$            close(14)
c$$$         endif
c$$$         if(intpar(7) .eq.1) then
c$$$            close(15)
c$$$         endif
c$$$         if(intpar(8) .eq.1) then
c$$$            close(16)
c$$$         endif
c$$$      endif
c$$$      if(intpar(9) .eq.1) then
c$$$         close(17)
c$$$      endif
c$$$      if(intpar(10) .eq.1) then
c$$$         close(18)
c$$$      endif
c$$$      if(intpar(11) .eq.1) then
c$$$         close(19)
c$$$      endif
c$$$      if(intpar(12) .eq.1) then
c$$$         close(20)
c$$$      endif
c$$$      if(useqtc .eq. 1) then 
c$$$         if(intpar(13) .eq.1) then
c$$$            close(21)
c$$$         endif
c$$$         if(intpar(14) .eq.1) then
c$$$            close(22)
c$$$         endif
c$$$         if(intpar(15) .eq.1) then
c$$$            close(23)
c$$$         endif
c$$$      endif

      call intpr('************************************',-1,0,0)
      call intpr('***    End of MCMC simulation    ***',-1,0,0)
      call intpr('************************************',-1,0,0)


c$$$
c$$$      write(*,*) 'nindiv =', nindiv
c$$$      write(*,*) 'nlocd =',nlocd   
c$$$      write(*,*) 'nloch =',nloch   
c$$$      write(*,*) 'nql =',nql   
c$$$      write(*,*) 'ncolt =',ncolt   
c$$$      write(*,*) 'nal =',nal 
c$$$      write(*,*) 'nalmax =', nalmax
c$$$      write(*,*) 'usegeno2=',usegeno2
c$$$      write(*,*) 'usegeno1=',usegeno1
c$$$      write(*,*) 'useql=',useql
c$$$      write(*,*) 'yy =', yy
c$$$      write(*,*) 'z =', z
c$$$      write(*,*) 'ql =', ql
c$$$      write(*,*) 's =', s
c$$$      write(*,*) 't =', t
c$$$      write(*,*) 'dt =', dt
c$$$      write(*,*) 'lambda=',lambda
c$$$      write(*,*) 'lambdamax=',lambdamax  
c$$$      write(*,*) 'npp =', npp 
c$$$      write(*,*) 'nppmax =', nppmax     
c$$$      write(*,*) 'npopmin =',npopmin  
c$$$      write(*,*) 'npopmax =', npopmax  
c$$$      write(*,*) 'c =', c
c$$$      write(*,*) 'ctmp =',ctmp 
c$$$      write(*,*) 'u =',u  
c$$$      write(*,*) 'utmp =',utmp 
c$$$      write(*,*) 'f =', f  
c$$$      write(*,*) 'ftmp =', ftmp
c$$$      write(*,*) 'nit =', nit  
c$$$      write(*,*) 'thinning =',thinning  
c$$$      write(*,*) 'indcell =', indcell 
c$$$      write(*,*) 'distcell =', distcell
c$$$      write(*,*) 'indcelltmp =', indcelltmp 
c$$$      write(*,*) 'distcelltmp =',distcelltmp 
c$$$      write(*,*) 'xlim=', xlim 
c$$$      write(*,*) 'ylim=', ylim 
c$$$      write(*,*) 'intpar=',intpar
c$$$      write(*,*) 'qtc=',qtc
c$$$      write(*,*) 'meanqtc=',meanqtc
c$$$      write(*,*) 'sdqtc=',sdqtc
c$$$      write(*,*) 'nnqtc=',nnqtc
c$$$      write(*,*) 'intpar=',intpar
c$$$      write(*,*) 'dblepar=',dblepar
c$$$      write(*,*) 'pudcel=',pudcel
c$$$      write(*,*) 'missloc=',missloc

c$$$      write(*,*) "avant rndend ggrunif=",ggrunif(0.d0,1.d0) 
      call rndend()
      end subroutine mcmcgld
****  end of main subroutine *********



**************************************************************
*     Update components of parameter beta in distributions 
*     of variances of continuous quantitative variables
*     (as in Richardson and Green JRSS B 1997)
      subroutine rpostbetaqtc(nqtc,npop,npopmax, 
     &     alphaqtc,betaqtc,gbeta,hbeta,sdqtc)
      implicit none
      integer nqtc,npop,npopmax
      double precision alphaqtc,betaqtc,gbeta,hbeta,sdqtc
      dimension alphaqtc(nqtc),betaqtc(nqtc),gbeta(nqtc),hbeta(nqtc),
     &     sdqtc(npopmax,nqtc)
      integer iqtc,ipop
      double precision aa,bb,ggrgam
      do iqtc = 1,nqtc
         aa = gbeta(iqtc) + npop*alphaqtc(iqtc)
         bb = hbeta(iqtc)
         do ipop = 1,npop
            bb = bb + 1/(sdqtc(ipop,iqtc)**2)
         enddo
         betaqtc(iqtc) = ggrgam(aa,1/bb)
      enddo
      end subroutine rpostbetaqtc



**************************************************************
*     sample mean and variance of quantitative variables 
*     from the posterior
*     prior mean : Normal
*     prior 1/sd**2 : Gamma
*     mean and sd a priori independent
*     qtc | mean,sd (likelihood) : Normal
      subroutine rpostqtc(qtc,nindiv,nqtc,indcell,c,meanqtc,sdqtc,
     &     nnqtc,sqtc,ssqtc,npop,npopmax,nppmax,
     &     ksiqtc,kappaqtc,alphaqtc,betaqtc)
      implicit none
      integer nindiv,nqtc,indcell,c,nnqtc,npop,npopmax,nppmax
      double precision qtc,meanqtc,sdqtc,sqtc,ssqtc,
     &     ksiqtc,kappaqtc,alphaqtc,betaqtc
      dimension qtc(nindiv,nqtc),indcell(nindiv),c(nppmax),
     &     meanqtc(npopmax,nqtc),sdqtc(npopmax,nqtc),
     &     nnqtc(npopmax,nqtc),sqtc(npopmax,nqtc),
     &     ssqtc(npopmax,nqtc),
     &     ksiqtc(nqtc),kappaqtc(nqtc),alphaqtc(nqtc),betaqtc(nqtc)
      integer ipop,iindiv,iqtc
      double precision xx,vv,ap,bp,ggrnorm,ggrgam

*     compute empirical sums and sums of squares for quant. variables
      do iqtc = 1,nqtc
         do ipop = 1,npop
            nnqtc(ipop,iqtc) = 0
            sqtc(ipop,iqtc) = 0
            ssqtc(ipop,iqtc) = 0
         enddo
      enddo
      do iqtc = 1,nqtc
         do iindiv = 1,nindiv
            ipop = c(indcell(iindiv))
            if(qtc(iindiv,iqtc) .ne. -999) then
               nnqtc(ipop,iqtc) = nnqtc(ipop,iqtc) + 1
               sqtc(ipop,iqtc) = sqtc(ipop,iqtc) + 
     &              qtc(iindiv,iqtc)
               ssqtc(ipop,iqtc) = ssqtc(ipop,iqtc) + 
     &              qtc(iindiv,iqtc)**2
            endif
         enddo
      enddo

      do iqtc = 1,nqtc
         do ipop = 1,npop
             xx = (sqtc(ipop,iqtc)+
     &           ksiqtc(iqtc)*kappaqtc(iqtc)*sdqtc(ipop,iqtc)**2)/
     &           (nnqtc(ipop,iqtc)+kappaqtc(iqtc)*sdqtc(ipop,iqtc)**2)
             vv = sdqtc(ipop,iqtc)**2 /
     &            (nnqtc(ipop,iqtc)+kappaqtc(iqtc)*sdqtc(ipop,iqtc)**2)
             meanqtc(ipop,iqtc) = ggrnorm(xx,dsqrt(vv))
             ap = alphaqtc(iqtc) + 0.5*nnqtc(ipop,iqtc)
             bp = betaqtc(iqtc) + 0.5*(ssqtc(ipop,iqtc) - 
     &            2*sqtc(ipop,iqtc)*meanqtc(ipop,iqtc) + 
     &            nnqtc(ipop,iqtc)*meanqtc(ipop,iqtc)**2)
             sdqtc(ipop,iqtc) = 1/dsqrt(ggrgam(ap,1/bp))
c$$$             write(6,*) 'in rpostqtc xx=',xx
c$$$             write(6,*) 'in rpostqtc vv=',vv
c$$$             write(6,*) 'in rpostqtc nnqtc=',nnqtc(ipop,iqtc)
c$$$             write(6,*) 'in rpostqtc sqtc=',sqtc(ipop,iqtc)
c$$$             write(6,*) 'in rpostqtc ssqtc=',ssqtc(ipop,iqtc)
c$$$             write(6,*) 'in rpostqtc meanqtc=',meanqtc(ipop,iqtc)
c$$$             write(6,*) 'in rpostqtc ap=',ap
c$$$             write(6,*) 'in rpostqtc bp=',bp
c$$$             write(6,*) 'in rpostqtc sdqtc=',sdqtc(ipop,iqtc)
         enddo
      enddo
      end subroutine rpostqtc
**************************************************************     





**************************************************************
*     sample mean and variance of quantitative variables 
*     from the posterior
*     prior mean : Normal
*     prior 1/sd**2 : Gamma
*     mean and sd a priori independent
*     qtc | mean,sd (likelihood) : Normal
      subroutine rpostqtc2(qtc,nindiv,nqtc,indcell,c,meanqtc,sdqtc,
     &     nnqtc,sqtc,ssqtc,npop,npopmax,nppmax,
     &     ksiqtc,kappaqtc,alphaqtc,betaqtc)
      implicit none
      integer nindiv,nqtc,indcell,c,nnqtc,npop,npopmax,nppmax
      double precision qtc,meanqtc,sdqtc,sqtc,ssqtc,
     &     ksiqtc,kappaqtc,alphaqtc,betaqtc
      dimension qtc(nindiv,nqtc),indcell(nindiv),c(nppmax),
     &     meanqtc(npopmax,nqtc),sdqtc(npopmax,nqtc),
     &     nnqtc(npopmax,nqtc),sqtc(npopmax,nqtc),
     &     ssqtc(npopmax,nqtc),
     &     ksiqtc(nqtc),kappaqtc(nqtc),alphaqtc(nqtc),betaqtc(nqtc)
      integer ipop,iindiv,iqtc
      double precision xx,vv,ap,bp,ggrnorm,ggrgam

c      write(*,*) 'c=',c

*     compute empirical sums and sums of squares for quant. variables
      do iqtc = 1,nqtc
         do ipop = 1,npop
            nnqtc(ipop,iqtc) = 0
            sqtc(ipop,iqtc) = 0
            ssqtc(ipop,iqtc) = 0
         enddo
      enddo
      do iqtc = 1,nqtc
         do iindiv = 1,nindiv
            ipop = c(indcell(iindiv))
            if(qtc(iindiv,iqtc) .ne. -999) then
               nnqtc(ipop,iqtc) = nnqtc(ipop,iqtc) + 1
               sqtc(ipop,iqtc) = sqtc(ipop,iqtc) + 
     &              qtc(iindiv,iqtc)
               ssqtc(ipop,iqtc) = ssqtc(ipop,iqtc) + 
     &              qtc(iindiv,iqtc)**2
            endif
         enddo
      enddo
      do iqtc = 1,nqtc
         do ipop = 1,npop
            if(nnqtc(ipop,iqtc) .ne. 0) then
               ssqtc(ipop,iqtc) = ssqtc(ipop,iqtc) -
     &              (sqtc(ipop,iqtc)**2) / dble(nnqtc(ipop,iqtc))
            else
               ssqtc(ipop,iqtc) = 0
            endif
         enddo
      enddo
      do iqtc = 1,nqtc
         do ipop = 1,npop
             ap = alphaqtc(iqtc) + 0.5*nnqtc(ipop,iqtc)
             if(nnqtc(ipop,iqtc) .ne. 0) then 
                bp = betaqtc(iqtc) + 0.5*ssqtc(ipop,iqtc) + 
     &               0.5*(nnqtc(ipop,iqtc)*kappaqtc(iqtc))/
     &               (nnqtc(ipop,iqtc)+kappaqtc(iqtc)) *
     &               (sqtc(ipop,iqtc)/dble(nnqtc(ipop,iqtc)) - 
     &               ksiqtc(iqtc))**2
             else
                bp = betaqtc(iqtc)
             endif 
             sdqtc(ipop,iqtc) = 1/dsqrt(ggrgam(ap,1/bp))
             xx = (sqtc(ipop,iqtc)+ksiqtc(iqtc)*kappaqtc(iqtc))/
     &        (nnqtc(ipop,iqtc)+kappaqtc(iqtc))
             vv = (sdqtc(ipop,iqtc)**2) / 
     &            (nnqtc(ipop,iqtc)+kappaqtc(iqtc))
             meanqtc(ipop,iqtc) = ggrnorm(xx,dsqrt(vv))
c$$$             write(*,*) 'iqtc,ipop =',iqtc,ipop
c$$$             write(*,*) 'ap =',ap
c$$$             write(*,*) 'bp =',bp
c$$$             write(*,*) 'xx =',xx
c$$$             write(*,*) 'vv =',vv
         enddo
      enddo
c$$$      write(*,*) ''
c$$$      write(*,*) 'fin rpostqtc2'
c$$$      write(*,*) 'meanqtc=',meanqtc
c$$$      write(*,*) 'sdqtc=',sdqtc
      end subroutine rpostqtc2
**************************************************************     



*************************************************************************      
*     update matrix of true genotypes
*     if given genotypes are true corrupted by the presence
*     of null alleles
      subroutine udNA3(nindiv,nlocd,nloch,nal,nalmax,ncolt,nppmax,
     &     yy,zz,c,indcell,npopmax,f,fcy,npop,missloc)
      implicit none
      integer nindiv,nlocd,nloch,nal,nalmax,ncolt,yy,zz,npopmax,npop,
     &     nppmax,c,indcell,missloc
      double precision f,fcy
      dimension nal(ncolt),yy(nindiv,2*nlocd+2*nloch),
     &     zz(nindiv,2*nlocd),
     &     f(npopmax,ncolt,nalmax),fcy(nalmax,2),c(nppmax),
     &     indcell(nindiv),missloc(nindiv,nlocd)
      integer iindiv,iloc,ial1,ipop,ial2,ial,alpha
      double precision u,ggrunif,sp
c      write(*,*) 'debut udNA'

      do iloc = 1,nlocd
         do ipop = 1,npop
*     computes posterior proba of true genotypes 
*     given allele freq and observed genotypes (i.e. true genotypes 
*     blurred by null alleles)
            call postpyNA2(ipop,iloc,f,nal,npopmax,nalmax,ncolt,fcy)
*     sample yy
            do iindiv = 1,nindiv
*     only for indiv in pop ipop 
               if(c(indcell(iindiv)) .eq. ipop) then
*     only for indiv with ambigous genotype (homozzygous)
                  if(zz(iindiv,2*iloc-1) .eq. zz(iindiv,2*iloc)) then
*     case doubly missing data NOT at a missing locus
                     if((zz(iindiv,2*iloc-1) .eq. -999) .and.
     &                    missloc(iindiv,iloc) .eq. 0) then 
                        yy(iindiv,2*iloc-1) = nal(iloc)
                        yy(iindiv,2*iloc)   = nal(iloc)
                     else   
*     case homozzygous (non missing data)
                        u=ggrunif(0.d0,1.d0)
                        alpha = zz(iindiv,2*iloc-1)
                        if(u .le. fcy(alpha,1)) then 
                           yy(iindiv,2*iloc-1) = alpha
                           yy(iindiv,2*iloc)   = alpha
                        else
                           yy(iindiv,2*iloc-1) = alpha
                           yy(iindiv,2*iloc)   = nal(iloc)
                        endif
                     endif
                  endif
               endif
            enddo
         enddo
      enddo
c      write(*,*) 'fin udNA'
      end subroutine udNA3
*************************************************************************


*************************************************************************
*     computes posterior proba of true genotypes 
*     given allele freq and observed genotypes (obs. = true genotypes 
*     corrupted by null alleles)
*     if presence of null allele is assumed, an extra allele 
*     is assumed and information relative to this allele 
*     is stored in the last non empty entry of f,fcy,...
      subroutine postpyNA(ipop,iloc,f,nal,npopmax,nloc,nalmax,fcy)
      implicit none
      integer ipop,iloc,npopmax,nloc,nalmax,nal
      double precision f,fcy
      dimension f(npopmax,nloc,nalmax),nal(nloc),fcy(nalmax,2)
      integer ial
c$$$      write(*,*) 'debut postpyNA'
c$$$      write(*,*) 'nal=',nal
c$$$      write(*,*) 'nalmax=',nalmax 
c$$$      write(*,*) 'fcy=',fcy
c$$$      write(*,*) 'f=',f
*     visit all "genuine" alleles
      do ial = 1,nal(iloc)-1
*     proba to have a genuine homozziguous ial,ial
            fcy(ial,1) = f(ipop,iloc,ial)/
     &           (f(ipop,iloc,ial) + 2*f(ipop,iloc,nal(iloc)))
*     proba to have a false homozziguous
c            fcy(ial,2) = 2*f(ipop,iloc,nal(iloc))/
c     &           (f(ipop,iloc,ial) + 2*f(ipop,iloc,nal(iloc)))
      enddo
c      write(*,*) 'fin postpyNA'
      end subroutine postpyNA
************************************************************************




*************************************************************************
*     computes posterior proba of true genotypes 
*     given allele freq and observed genotypes (obs. = true genotypes 
*     corrupted by null alleles)
*     if presence of null allele is assumed, an extra allele 
*     is assumed and information relative to this allele 
*     is stored in the last non empty entry of f,fcy,...
      subroutine postpyNA2(ipop,iloc,f,nal,npopmax,nalmax,ncolt,fcy)
      implicit none
      integer ipop,iloc,npopmax,nalmax,ncolt,nal
      double precision f,fcy
      dimension f(npopmax,ncolt,nalmax),nal(ncolt),fcy(nalmax,2)
      integer ial
c$$$      write(*,*) 'debut postpyNA'
c$$$      write(*,*) 'nal=',nal
c$$$      write(*,*) 'nalmax=',nalmax 
c$$$      write(*,*) 'fcy=',fcy
c$$$      write(*,*) 'f=',f
*     visit all "genuine" alleles
      do ial = 1,nal(iloc)-1
*     proba to have a genuine homozziguous ial,ial
            fcy(ial,1) = f(ipop,iloc,ial)/
     &           (f(ipop,iloc,ial) + 2*f(ipop,iloc,nal(iloc)))
*     proba to have a false homozziguous
c            fcy(ial,2) = 2*f(ipop,iloc,nal(iloc))/
c     &           (f(ipop,iloc,ial) + 2*f(ipop,iloc,nal(iloc)))
      enddo
c      write(*,*) 'fin postpyNA'
      end subroutine postpyNA2
************************************************************************



*************************************************************************      
*     update part of matrix of true genotypes
*     that corresponds to dominant markers 
*     (stored in columns 2*nlocd+1 to 2*nlocd+2*nloch)
      subroutine udyDOM3(nindiv,nlocd,nloch,nal,nalmax,
     &    ncolt,nppmax,yy,z,c,indcell,npopmax,f,fcy,npop)
      implicit none
      integer nindiv,nlocd,nloch,nal,nalmax,yy,z,
     &     npopmax,npop,nppmax,c,indcell,ncolt
      double precision f,fcy
      dimension nal(ncolt),yy(nindiv,2*nlocd+2*nloch),z(nindiv,nloch),
     &     f(npopmax,ncolt,nalmax),fcy(nalmax,2),c(nppmax),
     &     indcell(nindiv)
      integer iindiv,iloc,ial1,ipop,ial2,ial,alpha
      double precision u,ggrunif,p
c      write(*,*) 'debut udyDOM'
      do iloc = 1,nloch
         do ipop = 1,npop
*     computes posterior proba of genuine homozigous
*     given allele freq and data 
            p = f(ipop,nlocd+iloc,2)/(f(ipop,nlocd+iloc,2) + 
     &           2*f(ipop,nlocd+iloc,1))
*     sample yy
            do iindiv = 1,nindiv
*     only for indiv in pop ipop 
               if(c(indcell(iindiv)) .eq. ipop) then
*     only for indiv with ambiguous obs. (presence of a band)
                  if(z(iindiv,iloc) .eq. 2) then
                     u=ggrunif(0.d0,1.d0)
                     if(u .le. p) then 
                        yy(iindiv,2*nlocd+2*iloc-1) = 2
                        yy(iindiv,2*nlocd+2*iloc)   = 2
                     else
                        yy(iindiv,2*nlocd+2*iloc-1) = 2
                        yy(iindiv,2*nlocd+2*iloc)   = 1
                     endif
                  endif
               endif
            enddo
         enddo
      enddo
c      write(*,*) 'fin udyDOM'
      end subroutine udyDOM3
************************************************************************




************************************************************************
*     Limites du rectangle contenant les coordonnees
      subroutine limit(nindiv,s,xlim,ylim,dt)
      implicit none
      integer nindiv
      double precision s(2,nindiv),xlim(2),ylim(2),dt
      integer iindiv
      xlim(1) = 1.d+300
      xlim(2) = -1.d+300
      ylim(1) = 1.d+300
      ylim(2) = -1.d+300
      do iindiv=1,nindiv
         xlim(1) = dmin1(s(1,iindiv),xlim(1))
         xlim(2) = dmax1(s(1,iindiv),xlim(2))
         ylim(1) = dmin1(s(2,iindiv),ylim(1))
         ylim(2) = dmax1(s(2,iindiv),ylim(2))
      enddo
      xlim(1) = xlim(1) - dt*.5
      xlim(2) = xlim(2) + dt*.5
      ylim(1) = ylim(1) - dt*.5
      ylim(2) = ylim(2) + dt*.5
      end
************************************************************************

************************************************************************
*     points uniformes dans [0,1]x[0,1]
      subroutine rprioru(npp,nppmax,xlim,ylim,u)
      implicit none
      integer npp,nppmax
      double precision u(2,nppmax),ggrunif,xlim(2),ylim(2)
      integer i
c      call intpr('Begin init u ',-1,0,0)
c      write(*,*) 'npp=',npp
c      write(*,*) 'nppmax=',nppmax
c      write(*,*) 'xlim=',xlim
c      write(*,*) 'ylim=',ylim
c      write(*,*) 'u=',u
      do i=1,npp
         u(1,i) = xlim(1)+(xlim(2)-xlim(1))*ggrunif(0.d0,1.d0)
         u(2,i) = ylim(1)+(ylim(2)-ylim(1))*ggrunif(0.d0,1.d0)
      enddo
      if(nppmax .gt. npp) then
         do i=npp+1, nppmax
            u(1,i) = -999.
            u(2,i) = -999.
         enddo
      endif
c      call intpr('End init u ',-1,0,0)
      end
************************************************************************

************************************************************************
*     affectation dans les pops selon une loi uniforme
      subroutine rpriorc(npp,nppmax,npop,c)
      implicit none
      integer npp,nppmax,npop,c(nppmax)
      double precision ggrunif
      integer i
      do i=1,npp
         c(i) = 1+ idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
      enddo
      if(nppmax .gt. npp) then
         do i=npp+1,nppmax
            c(i) = -999
         enddo
      endif
      end
************************************************************************


********************************************************************
*     init of vector of drift parameters
*     prior beta(shape1,shape2)
      subroutine rpriordrift(npop,npopmax,drift,fmodel,shape1,shape2)
      implicit none
      integer npop,npopmax,fmodel
      double precision drift(npopmax),shape1,shape2
      integer ipop
      double precision ggrbet
      if(fmodel .eq. 0) then
         do ipop=1,npop
            drift(ipop) = 0.5d0
         enddo
      else
         do ipop=1,npop
            drift(ipop) = ggrbet(shape1,shape2)
         enddo
      endif
      if(npopmax .gt. npop) then
         do ipop=npop+1,npopmax
            drift(ipop) = -999
         enddo
      endif
      end subroutine rpriordrift
********************************************************************


********************************************************************
*     simulation d'une Dirichlet(1,...,1)
*     (p(1),...,p(k)) uniforme dans l'ensemble {p(1)+...+p(k)=1}
      subroutine dirichlet1(nal,nalmax,p)
      implicit none
      integer nal,nalmax
      double precision p(nalmax)
      integer i
      double precision s,ggrexp
      s = 0.
      do i=1,nal
         p(i) = ggrexp(1.d0)
         s = s + p(i)
      enddo
      do i=1,nal
         p(i) =  p(i)/s
      enddo
      if(nalmax .gt. nal) then
         do i=nal+1,nalmax
            p(i) =  -1
         enddo
      endif
      end
**************************************************************


**************************************************************
*     simulation d'une Dirichlet(a1,...,an)
      subroutine dirichlet(n,nmax,a,p)
      implicit none
      integer n,nmax
      double precision a(nmax),p(nmax)
      integer i
      double precision s,ggrgam
c      write(*,*) 'debut dirichlet'
      s = 0.
      do i=1,n
         p(i) = 0.
         do while(p(i) .lt. 1d-300) 
c            p(i) = ggrgam(1.d0,a(i))
            p(i) = ggrgam(a(i),1.d0)
         enddo
         s = s + p(i)
      enddo
      do i=1,n
         p(i) =  p(i)/s
      enddo
      if(nmax .gt. n) then
         do i=n+1,nmax
            p(i) =  -1
         enddo
      endif
c      write(*,*) 'fin dirichlet'
      end
**************************************************************


**************************************************************
      subroutine rank(n,nmax,x,p)
      implicit none 
      integer n,nmax,p(nmax)
      double precision x(nmax)
      integer i,j
      p(1) = 1
      do i=2,n
         p(i) = 1
         do j=1,i-1
            if(x(i) .lt. x(j)) then 
               p(i) = p(i) + 1
            else 
               p(j) = p(j) + 1
            endif
         enddo
      enddo      
      end
**************************************************************

**************************************************************
*     from numerical recipe p 233
      subroutine indexx(n,nmax,arrin,indx)
      dimension arrin(nmax),indx(nmax),arrtmp(nmax)

      do j=1,n
         arrtmp(j) = - arrin(j)
      enddo
      do j=1,n
         indx(j) = j
      enddo
      if(nmax .gt. n) then 
         do j=n+1, nmax
            indx(j) = j
         enddo
      endif
      if(n .eq. 1) return
      l=n/2+1
      ir = n
 10   continue
      if(l .gt. 1) then
         l = l-1
         indxt = indx(l)
         q = arrtmp(indxt)
      else
         indxt = indx(ir)
         q = arrtmp(indxt)
         indx(ir) = indx(1)
         ir = ir -1
         if(ir .eq. 1) then 
            indx(1) = indxt
            return
         endif
      endif
      i = l
      j = l+l
 20   if(j .le. ir) then
         if(j .lt. ir) then
            if(arrtmp(indx(j)) .lt. arrtmp(indx(j+1))) j = j+1
         endif
         if(q .lt. arrtmp(indx(j))) then
            indx(i) = indx(j)
            i = j
            j = j+j
         else
            j = ir+1
         endif
         go to 20
      endif
      indx(i) = indxt
      go to 10
      end
**************************************************************


**************************************************************
*     tirage des frequences dans toutes les pops
*     a tous les locus
      subroutine rpriorf(npop,npopmax,nloc,nlocmax,nal,nalmax,f,
     &     ptmp)
      implicit none
      integer npop,npopmax,nloc,nlocmax,nal(nlocmax),nalmax
      double precision f(npopmax,nlocmax,nalmax)
      integer k,l,i
      double precision ptmp(nalmax)
c      call intpr('in rpriorf',-1,0,0)
c      call intpr('npop=',-1,npop,1)
c      call intpr('nloc=',-1,nloc,1)
c      call intpr('nalmax=',-1,nalmax,1)

      do  k=1,npop
c           call intpr('k=',-1,k,1)
         do l=1,nloc
c             call intpr('l=',-1,l,1)
            call dirichlet1(nal(l),nalmax,ptmp)
            do i=1,nalmax
c                call intpr('i=',-1,i,1)
               f(k,l,i)  = ptmp(i)
            enddo
         enddo
      enddo
c      call intpr('in rpriorf non dummy entries done',-1,0,0)

      if(npopmax .gt. npop) then
         do k=npop+1, npopmax
            do l=1,nloc
               do i=1,nalmax
                  f(k,l,i)  = -999
               enddo
            enddo
         enddo
      endif
      end subroutine rpriorf
**************************************************************


**************************************************************
*     tirage des frequences dans la pop ancestrale 
*     a tous les locus
      subroutine rpriorfa(nloc,nal,nalmax,fa,fmodel,ptmp)
      implicit none
      integer nloc,nal(nloc),nalmax,fmodel
      double precision fa(nloc,nalmax)
      integer l,i
      double precision ptmp(nalmax)
c      call intpr('in rpriorfa',-1,0,0)
      if(fmodel .eq. 0) then
         do l=1,nloc
            do i=1,nalmax
               fa(l,i)  = 1
            enddo
         enddo
      else
        do l=1,nloc
           call dirichlet1(nal(l),nalmax,ptmp)
           do i=1,nalmax
               fa(l,i)  = ptmp(i)
            enddo
         enddo 
      endif
      end subroutine rpriorfa
**************************************************************



**************************************************************
*     tirage des frequences dans la pop ancestrale 
*     a tous les locus
      subroutine rpriorfaallvar(nlocd,nloch,nql,ncolt,nal,nalmax,fa,
     &     fmodel,ptmp,usegeno2,usegeno1,useql)
      implicit none
      integer nlocd,nloch,nql,ncolt,nal(ncolt),nalmax,fmodel,
     &     usegeno2,usegeno1,useql
      double precision fa(ncolt,nalmax)
      integer iloc,ial
      double precision ptmp(nalmax)
c      call intpr('start rpriorfa',-1,0,0)
c      write(*,*) 'nal=',nal 
c      write(*,*) 'ncolt=',ncolt
c call intpr('nal',-1,nal,0)
      if(fmodel .eq. 0) then
         do iloc=1,ncolt
            do ial=1,nalmax
               fa(iloc,ial)  = 1
            enddo
         enddo
      else
         if(usegeno2 .eq. 1) then
            do iloc=1,nlocd
               call dirichlet1(nal(iloc),nalmax,ptmp)
               do ial=1,nalmax
                  fa(iloc,ial)  = ptmp(ial)
               enddo
            enddo 
         endif
         if(usegeno1 .eq. 1) then
            do iloc=nlocd+1,nlocd+nloch
               call dirichlet1(nal(iloc),nalmax,ptmp)
               do ial=1,nalmax
                  fa(iloc,ial)  = ptmp(ial)
               enddo
            enddo 
         endif
         if(useql .eq. 1) then
            do iloc=nlocd+nloch+1,nlocd+nloch+nql
               call dirichlet1(nal(iloc),nalmax,ptmp)
               do ial=1,nalmax
                  fa(iloc,ial)  = ptmp(ial)
               enddo
            enddo 
         endif
      endif
c      call intpr('end rpriorfa',-1,0,0)
      end subroutine rpriorfaallvar
**************************************************************





**********************************************************************
*     Mise a jour gibbsienne des frequences
*     prior p(f) Dirichlet
*     et une paramétrisation à la Falush (Genetics 2003)
*     p(f|...)  est aussi Dirichlet 
*     ni = nbre de modalites  observees
*     (Cf Falush P. 26)
      subroutine rpostf4(yy,z,ql,nindiv,nlocd,nloch,nql,ncolt,
     &     nal,nalmax,npop,npopmax,nppmax,c,indcell,n,a,ptmp,
     &     f,fa,drift,usegeno1,usegeno2,useql,ploidy)
      implicit none 
      integer nlocd,nlocd2,nloch,nql,npop,npopmax,ncolt,nal(ncolt),
     &     nalmax,nindiv,nppmax,c(nppmax),indcell(nindiv),
     &     usegeno1,usegeno2,useql,ploidy
      double precision f(npopmax,ncolt,nalmax),fa(ncolt,nalmax),
     &     drift(npopmax)
      integer ipop,iloc,iindiv,ial,n(npopmax,ncolt,nalmax),
     &     yy(nindiv,nlocd*2+nloch*2),z(nindiv,nloch),ql(nindiv,nql)
      double precision a(nalmax),ptmp(nalmax)
c$$$      write(*,*) ''
c$$$      write(*,*) ''
c$$$      write(*,*) 'dans rpostf2'
c$$$      write(*,*) 'npop=',npop
c$$$      write(*,*) 'npopmax=',npopmax
c$$$      write(*,*) 'nloc=',nloc
c$$$      write(*,*) 'nlocmax=',nlocmax
c$$$      write(*,*) 'nal=',nal
c$$$      write(*,*) 'nalmax=',nalmax
c$$$      write(*,*) 'f=',f
c$$$      write(*,*) 'fa=',fa
c$$$      write(*,*) 'drift=',drift
c$$$      write(*,*) 'nindiv=',nindiv
c$$$      write(*,*) 'nlocmax2=',nlocmax2
c$$$      write(*,*) 'zz=',zz
c$$$      write(*,*) 'nppmax=',nppmax
c$$$      write(*,*) 'c=',c
c$$$      write(*,*) 'indcell=',indcell
c$$$      write(*,*) 'n=',n
c$$$      write(*,*) 'a=',a
c$$$      write(*,*) 'ptmp=',ptmp
c$$$      
*     comptage des effectifs
      do ipop = 1,npop
         do iloc = 1,ncolt
            do ial =1,nal(iloc)
               n(ipop,iloc,ial)=0
            enddo
         enddo
      enddo
      do iindiv = 1,nindiv
         ipop = c(indcell(iindiv))
         if(usegeno2 .eq. 1) then
c     diploid geno
            do iloc = 1,nlocd
               if(yy(iindiv,2*iloc-1) .ne. -999) then
                  n(ipop,iloc,yy(iindiv,2*iloc-1)) = 
     &                 n(ipop,iloc,yy(iindiv,2*iloc-1)) + 1 
               endif
               if(yy(iindiv,2*iloc) .ne. -999) then 
                  n(ipop,iloc,yy(iindiv,2*iloc)) = 
     &                 n(ipop,iloc,yy(iindiv,2*iloc)) + 1 
               endif
            enddo
         endif
         if(usegeno1 .eq. 1) then
c     data stored in "haploid" genotype matrix
            if(ploidy .eq. 2) then
               do iloc = 1,nloch
                  if(yy(iindiv,2*nlocd+2*iloc-1) .ne. -999) then
                     n(ipop,nlocd+iloc,yy(iindiv,2*nlocd+2*iloc-1)) = 
     &               n(ipop,nlocd+iloc,yy(iindiv,2*nlocd+2*iloc-1)) + 1 
                  endif
                  if(yy(iindiv,2*nlocd+2*iloc) .ne. -999) then 
                     n(ipop,nlocd+iloc,yy(iindiv,2*nlocd+2*iloc)) = 
     &               n(ipop,nlocd+iloc,yy(iindiv,2*nlocd+2*iloc)) + 1 
                  endif
               enddo
            endif
            if(ploidy .eq. 1) then
               do iloc = 1,nloch
                  if(z(iindiv,iloc) .ne. -999) then
                     n(ipop,nlocd+iloc,z(iindiv,iloc)) = 
     &               n(ipop,nlocd+iloc,z(iindiv,iloc)) + 1
                  endif
               enddo
            endif
         endif
         if(useql .eq. 1) then
c     ql variables
            do iloc = 1,nql
               if(ql(iindiv,iloc) .ne. -999) then
                  n(ipop,nlocd+nloch+iloc,ql(iindiv,iloc)) = 
     &            n(ipop,nlocd+nloch+iloc,ql(iindiv,iloc)) + 1
               endif
            enddo
         endif
      enddo
  
*     tirage Dirichlet
      do ipop = 1,npop
         if(usegeno2 .eq. 1) then
c     diploid geno
            do iloc = 1,nlocd
               do ial = 1,nal(iloc)
                  a(ial) = fa(iloc,ial)*(1/drift(ipop)-1) + 
     &                 dble(n(ipop,iloc,ial))
               enddo
               call dirichlet(nal(iloc),nalmax,a,ptmp)
               do ial =1,nal(iloc)
                  f(ipop,iloc,ial) = ptmp(ial)
               enddo
            enddo
         endif
         if(usegeno1 .eq. 1) then
c     haploid geno
            do iloc = nlocd+1,nlocd+nloch
               do ial = 1,nal(iloc)
                  a(ial) = fa(iloc,ial)*(1/drift(ipop)-1) + 
     &                 dble(n(ipop,iloc,ial))
               enddo
               call dirichlet(nal(iloc),nalmax,a,ptmp)
               do ial = 1,nal(iloc)
                  f(ipop,iloc,ial) = ptmp(ial)
               enddo
            enddo
         endif
         if(useql .eq. 1) then
c     ql variables
            do iloc = nlocd+nloch+1,nlocd+nloch+nql
               do ial = 1,nal(iloc)
                  a(ial) = fa(iloc,ial)*(1/drift(ipop)-1) + 
     &                 dble(n(ipop,iloc,ial))
               enddo
               call dirichlet(nal(iloc),nalmax,a,ptmp)
               do ial =1,nal(iloc)
                  f(ipop,iloc,ial) = ptmp(ial)
               enddo
            enddo
         endif
      enddo

c$$$      write(*,*) ''
c$$$      write(*,*) ''
c$$$      write(*,*) 'dans fin rpostf2'
c$$$      write(*,*) 'npop=',npop
c$$$      write(*,*) 'npopmax=',npopmax
c$$$      write(*,*) 'nloc=',nloc
c$$$      write(*,*) 'nlocmax=',nlocmax
c$$$      write(*,*) 'nal=',nal
c$$$      write(*,*) 'nalmax=',nalmax
c$$$      write(*,*) 'f=',f
c$$$      write(*,*) 'fa=',fa
c$$$      write(*,*) 'drift=',drift
c$$$      write(*,*) 'nindiv=',nindiv
c$$$      write(*,*) 'nlocmax2=',nlocmax2
c$$$      write(*,*) 'zz=',zz
c$$$      write(*,*) 'nppmax=',nppmax
c$$$      write(*,*) 'c=',c
c$$$      write(*,*) 'indcell=',indcell
c$$$      write(*,*) 'n=',n
c$$$      write(*,*) 'a=',a
c$$$      write(*,*) 'ptmp=',ptmp
      
      end subroutine rpostf4
**********************************************************************



**********************************************************************
*
*     Mise à jour M-H des freq allelique de la pop ancestrale 
*     fa admet un prior Dirichlet(1,...,1)
      subroutine updfa(npop,npopmax,nlocmax,nalmax,nal,f,fa,drift)
      implicit none 
      integer npop,npopmax,nlocmax,nalmax,nal(nlocmax)
      double precision f(npopmax,nlocmax,nalmax),fa(nlocmax,nalmax),
     &     drift(npopmax)
      integer iloc,ial1,ial2,ipop
      double precision delta,ggrunif,ggrnorm,sigdelta,fa1,fa2,ratio,
     &     lratio,gglgamfn,q,u
      parameter(sigdelta = 0.05) 
      
      do iloc = 1,nlocmax
*     tirage des deux formes alleliques dont les freq seront mises à jour 
         ial1 = 1+ idint(dint(dble(nal(iloc))*ggrunif(0.d0,1.d0)))
         ial2 = ial1 
         do while(ial2 .eq. ial1)
c     write(*,*) 'dans le while'
            ial2 = 1+ idint(dint(dble(nal(iloc))*ggrunif(0.d0,1.d0)))
         enddo
*     tirage de l'increment
         delta = ggrnorm(0.d0,1.d0)*sigdelta
*     perturbation des deux freq
         fa1 = fa(iloc,ial1) + delta
         fa2 = fa(iloc,ial2) - delta
         if(((fa1 .gt. 1d-300) .and. (1-fa1 .gt. 1d-300)) .and.
     &        ((fa2 .gt. 1d-300) .and. (1-fa2 .gt. 1d-300))) then 
*     calcul du log du ratio 
            lratio = 0.
            do ipop = 1,npop
               q = (1.d0-drift(ipop))/drift(ipop)
               lratio = lratio 
     &              + gglgamfn(fa(iloc,ial1)*q)-gglgamfn(fa1*q)
     &              + gglgamfn(fa(iloc,ial2)*q)-gglgamfn(fa2*q)
     &              + delta*q
     &              *log(f(ipop,iloc,ial1)/f(ipop,iloc,ial2))
            enddo
            lratio = dmin1(0.d0,lratio) 
            ratio = dexp(lratio)
            u = ggrunif(0.d0,1.d0)
            if(u .le. ratio) then 
               fa(iloc,ial1) = fa1 
               fa(iloc,ial2) = fa2
            endif
         endif 
      enddo
      end subroutine updfa
************************************************************************   



**********************************************************************
*
*     Mise à jour M-H des freq allelique de la pop ancestrale 
*     fa admet un prior Dirichlet(1,...,1)
      subroutine updfa2(npop,npopmax,nlocd,nloch,nql,ncolt,nalmax,nal,
     &     f,fa,drift,usegeno2,usegeno1,useql)
      implicit none 
      integer npop,npopmax,nlocd,nloch,nql,ncolt,nalmax,nal(ncolt),
     &     usegeno2,usegeno1,useql
      double precision f(npopmax,ncolt,nalmax),fa(ncolt,nalmax),
     &     drift(npopmax)
      integer iloc,ial1,ial2,ipop
      double precision delta,ggrunif,ggrnorm,sigdelta,fa1,fa2,ratio,
     &     lratio,gglgamfn,q,u
      parameter(sigdelta = 0.05) 
      
c      write(*,*) 'debut de updfa2',ggrunif(0.d0,1.d0) 
c      write(*,*) 'nlocd=',nlocd
c      write(*,*) 'nloch=',nloch

c      write(*,*) usegeno2,usegeno1,useql
c      write(*,*) nlocd,nloch,nql,ncolt
      if(usegeno2 .eq. 1) then
c     diploid geno
         do iloc = 1,nlocd
c            write(*,*) 'iloc=',iloc
*     tirage des deux formes alleliques dont les freq seront mises à jour 
            ial1 = 1+ idint(dint(dble(nal(iloc))*ggrunif(0.d0,1.d0)))
            ial2 = ial1 
            do while(ial2 .eq. ial1)
c     write(*,*) 'dans le while'
               ial2 = 1+ idint(dint(dble(nal(iloc))*ggrunif(0.d0,1.d0)))
            enddo
*     tirage de l'increment
            delta = ggrnorm(0.d0,1.d0)*sigdelta
*     perturbation des deux freq
            fa1 = fa(iloc,ial1) + delta
            fa2 = fa(iloc,ial2) - delta
c            write(*,*) 'fa1=',fa1
c            write(*,*) 'fa2=',fa2
            if(((fa1 .gt. 1d-300) .and. (1-fa1 .gt. 1d-300)) .and.
     &           ((fa2 .gt. 1d-300) .and. (1-fa2 .gt. 1d-300))) then 
*     calcul du log du ratio 
               lratio = 0.
               do ipop = 1,npop
                  q = (1.d0-drift(ipop))/drift(ipop)
                  lratio = lratio 
     &                 + gglgamfn(fa(iloc,ial1)*q)-gglgamfn(fa1*q)
     &                 + gglgamfn(fa(iloc,ial2)*q)-gglgamfn(fa2*q)
     &                 + delta*q
     &                 *log(f(ipop,iloc,ial1)/f(ipop,iloc,ial2))
               enddo
c               write(*,*) 'lratio=',lratio
               lratio = dmin1(0.d0,lratio) 
               ratio = dexp(lratio)
               u = ggrunif(0.d0,1.d0)
               if(u .le. ratio) then 
                  fa(iloc,ial1) = fa1 
                  fa(iloc,ial2) = fa2
               endif
            endif 
         enddo
      endif
      if(usegeno1 .eq. 1) then
c     haploid geno
         do iloc = nlocd+1,nlocd+nloch
c            write(*,*) 'iloc=',iloc
*     tirage des deux formes alleliques dont les freq seront mises à jour 
            ial1 = 1+ idint(dint(dble(nal(iloc))*ggrunif(0.d0,1.d0)))
            ial2 = ial1 
            do while(ial2 .eq. ial1)
c     write(*,*) 'dans le while'
               ial2 = 1+ idint(dint(dble(nal(iloc))*ggrunif(0.d0,1.d0)))
            enddo
*     tirage de l'increment
            delta = ggrnorm(0.d0,1.d0)*sigdelta
*     perturbation des deux freq
            fa1 = fa(iloc,ial1) + delta
            fa2 = fa(iloc,ial2) - delta
c            write(*,*) 'fa1=',fa1
c            write(*,*) 'fa2=',fa2
            if(((fa1 .gt. 1d-300) .and. (1-fa1 .gt. 1d-300)) .and.
     &           ((fa2 .gt. 1d-300) .and. (1-fa2 .gt. 1d-300))) then 
*     calcul du log du ratio 
               lratio = 0.
               do ipop = 1,npop
                  q = (1.d0-drift(ipop))/drift(ipop)
                  lratio = lratio 
     &                 + gglgamfn(fa(iloc,ial1)*q)-gglgamfn(fa1*q)
     &                 + gglgamfn(fa(iloc,ial2)*q)-gglgamfn(fa2*q)
     &                 + delta*q
     &                 *log(f(ipop,iloc,ial1)/f(ipop,iloc,ial2))
               enddo
               lratio = dmin1(0.d0,lratio) 
               ratio = dexp(lratio)
               u = ggrunif(0.d0,1.d0)
               if(u .le. ratio) then 
                  fa(iloc,ial1) = fa1 
                  fa(iloc,ial2) = fa2
               endif
            endif 
         enddo
      endif
      if(useql .eq. 1) then
c     haploid geno
         do iloc = nlocd+nloch+1, nlocd+nloch+nql
*     tirage des deux formes alleliques dont les freq seront mises à jour 
            ial1 = 1+ idint(dint(dble(nal(iloc))*ggrunif(0.d0,1.d0)))
            ial2 = ial1 
            do while(ial2 .eq. ial1)
c     write(*,*) 'dans le while'
               ial2 = 1+ idint(dint(dble(nal(iloc))*ggrunif(0.d0,1.d0)))
            enddo
*     tirage de l'increment
            delta = ggrnorm(0.d0,1.d0)*sigdelta
*     perturbation des deux freq
            fa1 = fa(iloc,ial1) + delta
            fa2 = fa(iloc,ial2) - delta
            if(((fa1 .gt. 1d-300) .and. (1-fa1 .gt. 1d-300)) .and.
     &           ((fa2 .gt. 1d-300) .and. (1-fa2 .gt. 1d-300))) then 
*     calcul du log du ratio 
               lratio = 0.
               do ipop = 1,npop
                  q = (1.d0-drift(ipop))/drift(ipop)
                  lratio = lratio 
     &                 + gglgamfn(fa(iloc,ial1)*q)-gglgamfn(fa1*q)
     &                 + gglgamfn(fa(iloc,ial2)*q)-gglgamfn(fa2*q)
     &                 + delta*q
     &                 *log(f(ipop,iloc,ial1)/f(ipop,iloc,ial2))
               enddo
               lratio = dmin1(0.d0,lratio) 
               ratio = dexp(lratio)
               u = ggrunif(0.d0,1.d0)
               if(u .le. ratio) then 
                  fa(iloc,ial1) = fa1 
                  fa(iloc,ial2) = fa2
               endif
            endif 
         enddo
      endif
      end subroutine updfa2
************************************************************************   

************************************************************************
*     Mise à jour M-H du vecteur de dérives génétiques 
*     prior indep. beta sur chaque composante
      subroutine upddrift(npop,npopmax,nlocmax,nalmax,nal,
     &     f,fa,drift,shape1,shape2)
      implicit none 
      integer npop,npopmax,nlocmax,nalmax,nal(nlocmax)
      double precision drift(npopmax),f(npopmax,nlocmax,nalmax),
     &     fa(nlocmax,nalmax)
      integer ipop,iloc,ial
      double precision dtmp,q,qtmp,sigdelta,ratio,lratio,shape1,shape2,
     &     sall,ggrnorm,gglgamfn,u,ggrunif
c      parameter(sigdelta = 0.01)
      sigdelta = 0.5*shape1/(shape1+shape2)

*     boucle sur les pops
      do ipop=1,npop
*     proposition nouvelle valeur
         dtmp = drift(ipop) + ggrnorm(0.d0,1.d0)*sigdelta
         q = (1-drift(ipop))/drift(ipop)
         qtmp = (1-dtmp)/dtmp
         if((dtmp .gt. 1d-300 ) .and. (1-dtmp .gt. 1d-300)) then 

*     calcul du log du ratio
c     prior uniforme
            lratio = 0 
c     prior beta(shape1,shape2)
            lratio = (shape1-1)*dlog(dtmp/drift(ipop)) + 
     &           (shape2-1)*dlog((1-dtmp)/(1-drift(ipop)))
            do iloc=1,nlocmax
               sall = 0.
               do ial = 1,nal(iloc)
                  sall = sall + gglgamfn(fa(iloc,ial)*q)-
     &                 gglgamfn(fa(iloc,ial)*qtmp) +
     &                 fa(iloc,ial)*(qtmp-q)*dlog(f(ipop,iloc,ial))
               enddo
               lratio = lratio + sall + (gglgamfn(qtmp)-gglgamfn(q))
            enddo
 
            lratio = dmin1(0.d0,lratio)
            ratio = dexp(lratio)
            u = ggrunif(0.d0,1.d0)
            if(u .le. ratio) then 
               drift(ipop) = dtmp 
            endif
         endif
      enddo
      end subroutine upddrift
************************************************************************



************************************************************************
*     Mise à jour M-H du vecteur de dérives génétiques 
*     prior indep. beta sur chaque composante
      subroutine upddriftallvar(npop,npopmax,nlocd,nloch,nql,ncolt,
     &    nalmax,nal,f,fa,drift,shape1,shape2,usegeno2,usegeno1,
     &     useql)
      implicit none 
      integer npop,npopmax,nlocd,nloch,nql,ncolt,nalmax,nal(ncolt),
     &     usegeno2,usegeno1,useql
      double precision drift(npopmax),f(npopmax,ncolt,nalmax),
     &     fa(ncolt,nalmax)
      integer ipop,iloc,ial
      double precision dtmp,q,qtmp,sigdelta,ratio,lratio,shape1,shape2,
     &     sall,ggrnorm,gglgamfn,u,ggrunif
c      parameter(sigdelta = 0.01)
      sigdelta = 0.5*shape1/(shape1+shape2)

*     boucle sur les pops
      do ipop=1,npop
*     proposition nouvelle valeur
         dtmp = drift(ipop) + ggrnorm(0.d0,1.d0)*sigdelta
         q = (1-drift(ipop))/drift(ipop)
         qtmp = (1-dtmp)/dtmp
         if((dtmp .gt. 1d-300 ) .and. (1-dtmp .gt. 1d-300)) then 

*     calcul du log du ratio
c     prior uniforme
            lratio = 0 
c     prior beta(shape1,shape2)
            lratio = (shape1-1)*dlog(dtmp/drift(ipop)) + 
     &           (shape2-1)*dlog((1-dtmp)/(1-drift(ipop)))
            if(usegeno2 .eq. 1) then 
               do iloc=1,nlocd
                  sall = 0.
                  do ial = 1,nal(iloc)
                     sall = sall + gglgamfn(fa(iloc,ial)*q)-
     &                    gglgamfn(fa(iloc,ial)*qtmp) +
     &                    fa(iloc,ial)*(qtmp-q)*dlog(f(ipop,iloc,ial))
                  enddo
                  lratio = lratio + sall + (gglgamfn(qtmp)-gglgamfn(q))
               enddo
            endif
            if(usegeno1 .eq. 1) then 
               do iloc=1,nloch
                  sall = 0.
                  do ial = 1,nal(nlocd+iloc)
                     sall = sall + gglgamfn(fa(nlocd+iloc,ial)*q)-
     &                    gglgamfn(fa(nlocd+iloc,ial)*qtmp) +
     &                    fa(nlocd+iloc,ial)*(qtmp-q)*
     &                    dlog(f(ipop,nlocd+iloc,ial))
                  enddo
                  lratio = lratio + sall + (gglgamfn(qtmp)-gglgamfn(q))
               enddo
            endif
            if(useql .eq. 1) then 
               do iloc=1,nql
                  sall = 0.
                  do ial = 1,nal(nlocd+nloch+iloc)
                     sall = sall + gglgamfn(fa(nlocd+nloch+iloc,ial)*q)-
     &                    gglgamfn(fa(nlocd+nloch+iloc,ial)*qtmp) +
     &                    fa(nlocd+nloch+iloc,ial)*(qtmp-q)*
     &                    dlog(f(ipop,nlocd+nloch+iloc,ial))
                  enddo
                  lratio = lratio + sall + (gglgamfn(qtmp)-gglgamfn(q))
               enddo
            endif
            lratio = dmin1(0.d0,lratio)
            ratio = dexp(lratio)
            u = ggrunif(0.d0,1.d0)
            if(u .le. ratio) then 
               drift(ipop) = dtmp 
            endif
         endif
      enddo
      end subroutine upddriftallvar
************************************************************************





************************************************************************
*     recherche la cellule de chaque individu
*     stockage des indices dans indcell
*     stockage des carres des distances dans distcell
      subroutine calccell(nindiv,s,npp,nppmax,u,
     &     indcell,distcell)
      implicit none
      integer nindiv,npp,nppmax,indcell(nindiv)
      double precision distcell(nindiv)
      double precision s(2,nindiv),u(2,nppmax)
      integer iindiv,ipp
      double precision d
c      write(*,*) 's=',s
c      write(*,*) 'u=',u
      do iindiv=1,nindiv
         indcell(iindiv) = -999
         distcell(iindiv) = 1.d+300
         do ipp=1,npp
            d = (s(1,iindiv)-u(1,ipp))**2 + (s(2,iindiv)-u(2,ipp))**2
c           write(*,*) 'd=',d
            if( d .lt. distcell(iindiv) ) then 
               indcell(iindiv) = ipp
               distcell(iindiv) = d
c               write(*,*) 'ipp=',ipp
            endif
         enddo
      enddo
c      write(*,*) 'indcell =',indcell
      end




*     mise a jour de indcell et distcell
*     apres le deplacement d'un point de u (celui d'indice j)
      subroutine vormove(nindiv,s,npp,nppmax,u,
     &     indcell,distcell,indcelltmp,distcelltmp,j)
      implicit none
      integer nindiv,npp,nppmax,indcell(nindiv),
     &     indcelltmp(nindiv),j
      double precision s(2,nindiv),u(2,nppmax),distcell(nindiv),
     &     distcelltmp(nindiv),d
      integer iindiv,ipp

C       write(6,*) 'debut de  vormove'
C       write(6,*) 'j =', j
C       write(6,*) 'indcell',indcell
C       write(6,*)'distcell',distcell
C       write(6,*) 'indcelltmp',indcelltmp
C       write(6,*)'distcelltmp',distcelltmp

      do iindiv=1,nindiv
         if(indcell(iindiv) .eq. j) then 
*     pour les indiv qui etaient dans la cellule j on cherche
*     la nouvelle cellule
            d = 3.d+300
            indcelltmp(iindiv) = -999
            distcelltmp(iindiv) = 3.d+300
            do ipp=1,npp
               d= (s(1,iindiv)-u(1,ipp))**2+(s(2,iindiv)-u(2,ipp))**2
               if( d .lt. distcelltmp(iindiv) ) then 
                  indcelltmp(iindiv) = ipp
                  distcelltmp(iindiv) = d
               endif
            enddo
*     pour les autres indiv on regarde si le nouveau uj s'est intercale
         else
            d = (s(1,iindiv)-u(1,j))**2+(s(2,iindiv)-u(2,j))**2
            if(d .lt. distcell(iindiv)) then
               indcelltmp(iindiv) = j
               distcelltmp(iindiv) = d
            else
               indcelltmp(iindiv) = indcell(iindiv)
               distcelltmp(iindiv) = distcell(iindiv)
            endif
         endif
      enddo

c$$$         call calccell(nindiv,s,npp,nppmax,u,
c$$$     &        indcelltmp2,distcelltmp2)
c$$$         do iindiv=1,nindiv
c$$$            if((indcelltmp2(iindiv) .ne. indcelltmp(iindiv)) .or. 
c$$$     &         (distcelltmp2(iindiv) .ne. distcelltmp(iindiv)))  then
c$$$               write(6,*) 'fin de  vormove'
c$$$               write(6,*) 'j =', j
c$$$               write(6,*) 'iindiv=',iindiv
c$$$               write(6,*) 'indcell',indcell
c$$$               write(6,*)'distcell',distcell
c$$$               write(6,*) 'indcelltmp',indcelltmp
c$$$               write(6,*)'distcelltmp',distcelltmp
c$$$               write(6,*)'indceltmp2',indcelltmp2
c$$$               write(6,*)'distceltmp2',distcelltmp2
c$$$               stop
c$$$            endif
c$$$         enddo
C          write(6,*) 'fin de  vormove'
C          write(6,*) 'j =', j
C          write(6,*) 'indcell',indcell
C          write(6,*)'distcell',distcell
C          write(6,*) 'indcelltmp',indcelltmp
C          write(6,*)'distcelltmp',distcelltmp
      end
************************************************************************

************************************************************************
*     mise a jour de indcell et distcell
*     apres naissance d'un point de u 
      subroutine voradd(s,utmp,
     &     indcell,distcell,indcelltmp,distcelltmp,
     &     nindiv,npp,nppmax)
      implicit none 
      integer nindiv,npp,nppmax,indcell(nindiv),
     &     indcelltmp(nindiv),iindiv
      double precision s(2,nindiv),distcell(nindiv),
     &     distcelltmp(nindiv),d,utmp(2,nppmax)
      
      do iindiv =1,nindiv
*     est-ce que le nouveau point s'est intercale ?
         d = (s(1,iindiv)-utmp(1,npp+1))**2+
     &        (s(2,iindiv)-utmp(2,npp+1))**2
         if(d .lt. distcell(iindiv)) then 
            distcelltmp(iindiv) = d
            indcelltmp(iindiv) = npp+1
         else
            distcelltmp(iindiv) = distcell(iindiv)
            indcelltmp(iindiv) = indcell(iindiv) 
         endif
      enddo
      end 
************************************************************************
 


************************************************************************
*     mise a jour de indcell et distcell
*     apres mort d'un point de u 
      subroutine vorrem(s,utmp,ipprem,
     &     indcell,distcell,indcelltmp,distcelltmp,
     &     nindiv,npp,nppmax)
      implicit none
      integer nindiv,npp,nppmax,indcell(nindiv),
     &     indcelltmp(nindiv),ipprem,iindiv
      double precision s(2,nindiv),utmp(2,nppmax),
     &     distcell(nindiv),distcelltmp(nindiv),d
      integer ipp
      
      do iindiv =1,nindiv
*     est-ce que le site courant dependait de la cellule disparue ?
         if(indcell(iindiv) .eq. ipprem) then
*     si oui on recherche sa nouvelle cellule parmi celles qui restent
*     (les nouvelles)
            distcelltmp(iindiv) = 3.d+300
            do ipp=1,npp-1
               d = (s(1,iindiv)-utmp(1,ipp))**2+
     &              (s(2,iindiv)-utmp(2,ipp))**2
               if( d .lt. distcelltmp(iindiv) ) then 
                  indcelltmp(iindiv) = ipp
                  distcelltmp(iindiv) = d
               endif
            enddo
         else
            if(indcell(iindiv) .lt. ipprem) then
               indcelltmp(iindiv) = indcell(iindiv)
               distcelltmp(iindiv) = distcell(iindiv)
            else
               indcelltmp(iindiv) = indcell(iindiv) - 1
               distcelltmp(iindiv) = distcell(iindiv)
            endif
         endif
      enddo
      end


*     Tirage de lambda selon p(lambda|m) 
*     avec 
*     p(m|lambda) Poisson translatee
*     p(lambda) uniforme dans [0,lambdamax]
*     p(lambda|m) gamma tronquee
      double precision function rpostlamb(lambdamax,m)
      implicit none 
      double precision lambdamax,ggrgam
      integer m
*      write(*,*) 'beg rpostlamb'
      rpostlamb = lambdamax + 1
      do while(rpostlamb .gt. lambdamax)
c         rpostlamb = ggrgam(1.d0,dble(m))
         rpostlamb = ggrgam(dble(m),1.d0)
      enddo
*      write(*,*) 'end rpostlamb'
      end
************************************************************************


************************************************************************     
*     Mise a jour de c sans modif de npp
      subroutine updc(npp,nppmax,c,ctmp,zz,nindiv,nloc,
     &     nlocmax,nlocmax2,nalmax,npop,npopmax,f,indcell,ploidy,
     &     nudcel)
      implicit none 
      integer npp, nppmax,c(nppmax),nindiv,nloc,nlocmax,
     &     nlocmax2,npop,nalmax,npopmax,zz(nindiv,nlocmax2),
     &     indcell(nindiv),ploidy,nudcel
      double precision f(npopmax,nlocmax,nalmax)
      integer ipp,ctmp(nppmax),iud
      double precision ggrunif,r,alpha,ratio,ggrbinom,bern
      
      do ipp=1,npp
         ctmp(ipp) = c(ipp)
      enddo
      if(nppmax .gt. npp) then
         do ipp=npp+1,nppmax
            ctmp(ipp) = -999
         enddo
      endif
c     write(*,*) ''
      do iud=1,nudcel
         ipp = 1 + idint(dint(dble(npp)*ggrunif(0.d0,1.d0)))
         ctmp(ipp) = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
         r = ratio(zz,f,c,ctmp,indcell,indcell,
     &     npopmax,nlocmax,nalmax,nindiv,nloc,nlocmax2,
     &     nppmax,ploidy)
c         write(*,*) 'c=',c
c         write(*,*) 'ctmp=',ctmp
c         write(*,*) 'r=',r
         alpha = dmin1(1.d0,r)
         bern = ggrbinom(1.d0,alpha)
         if(bern .eq. 1) then
            c(ipp) = ctmp(ipp)
         else 
            ctmp(ipp) = c(ipp)
         endif
      enddo
      end subroutine updc
************************************************************************

************************************************************************     
*     Update c when data consist of genotypes and/or quantitative variables
      subroutine updcgq(npp,nppmax,c,ctmp,zz,nindiv,nloc,
     &     nlocmax,nlocmax2,nalmax,npop,npopmax,f,indcell,ploidy,
     &     nudcel,qtc,nqtc,meanqtc,sdqtc,usegeno2,useqtc)
      implicit none 
      integer npp, nppmax,c(nppmax),nindiv,nloc,nlocmax,
     &     nlocmax2,npop,nalmax,npopmax,zz(nindiv,nlocmax2),
     &     indcell(nindiv),ploidy,nudcel,nqtc,usegeno2,useqtc
      double precision f(npopmax,nlocmax,nalmax),qtc(nindiv,nqtc),
     &     meanqtc(npopmax,nqtc),sdqtc(npopmax,nqtc)
      integer ipp,ctmp(nppmax),iud
      double precision ggrunif,r,alpha,lratiogq,ggrbinom,bern
      
      do ipp=1,npp
         ctmp(ipp) = c(ipp)
      enddo
      if(nppmax .gt. npp) then
         do ipp=npp+1,nppmax
            ctmp(ipp) = -999
         enddo
      endif
c     write(*,*) ''
      do iud=1,nudcel
         ipp = 1 + idint(dint(dble(npp)*ggrunif(0.d0,1.d0)))
         ctmp(ipp) = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
         r = dexp(lratiogq(zz,f,c,ctmp,indcell,indcell,
     &     npopmax,nlocmax,nalmax,nindiv,nloc,nlocmax2,
     &     nppmax,ploidy,qtc,nqtc,meanqtc,sdqtc,meanqtc,sdqtc,
     &        usegeno2,useqtc))
c         write(*,*) 'c=',c
c         write(*,*) 'ctmp=',ctmp
c         write(*,*) 'r=',r
         alpha = dmin1(1.d0,r)
         bern = ggrbinom(1.d0,alpha)
         if(bern .eq. 1) then
            c(ipp) = ctmp(ipp)
         else 
            ctmp(ipp) = c(ipp)
         endif
      enddo
      end subroutine updcgq
************************************************************************




************************************************************************     
*     Update c when data consist of genotypes and/or phenotypes
      subroutine udcallvar2(npp,nppmax,c,ctmp,yy,z,ql,nindiv,nlocd,
     &     nloch,nql,ncolt,nalmax,npop,npopmax,f,indcell,nudcel,
     &     qtc,nqtc,meanqtc,sdqtc,usegeno2,usegeno1,useql,useqtc,
     &     ploidy)

      implicit none 
      integer npp, nppmax,c(nppmax),nindiv,nlocd,nloch,nql,ncolt,
     &     npop,nalmax,npopmax,yy(nindiv,2*nlocd+2*nloch),
     &     z(nindiv,nloch),ql(nindiv,nql),indcell(nindiv),nudcel,nqtc,
     &     usegeno2,usegeno1,useql,useqtc,ploidy
      double precision f(npopmax,ncolt,nalmax),qtc(nindiv,nqtc),
     &     meanqtc(npopmax,nqtc),sdqtc(npopmax,nqtc)
      integer ipp,ctmp(nppmax),iud
      double precision ggrunif,r,alpha,lrallvar2,ggrbinom,bern
      
      do ipp=1,npp
         ctmp(ipp) = c(ipp)
      enddo
      if(nppmax .gt. npp) then
         do ipp=npp+1,nppmax
            ctmp(ipp) = -999
         enddo
      endif
c     write(*,*) ''
      do iud=1,nudcel
         ipp = 1 + idint(dint(dble(npp)*ggrunif(0.d0,1.d0)))
         ctmp(ipp) = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
c         r = dexp(lratiogq(yy,f,c,ctmp,indcell,indcell,
c     &     npopmax,nlocmax,nalmax,nindiv,nloc,nlocmax2,
c     &     nppmax,ploidy,qtc,nqtc,meanqtc,sdqtc,meanqtc,sdqtc,
c     &        usegeno2,useqtc))
         r = dexp(lrallvar2(yy,z,ql,f,c,ctmp,indcell,indcell,
     &     npopmax,nlocd,nloch,nql,ncolt,nalmax,nindiv,nppmax,
     &     qtc,nqtc,meanqtc,sdqtc,meanqtc,sdqtc,
     &     usegeno2,usegeno1,useql,useqtc,ploidy))
c         write(*,*) 'c=',c
c         write(*,*) 'ctmp=',ctmp
c         write(*,*) 'r=',r
         alpha = dmin1(1.d0,r)
         bern = ggrbinom(1.d0,alpha)
         if(bern .eq. 1) then
            c(ipp) = ctmp(ipp)
         else 
            ctmp(ipp) = c(ipp)
         endif
      enddo
      end subroutine udcallvar2
************************************************************************





***********************************************************************
*     joint update of c and f under the Dirichlet model
*     single component update of c
*     new f is proposed according to full conditional pi(f*|c*,zz)
      subroutine udcf(npop,npopmax,f,
     &     nloc,nlocmax,nlocmax2,
     &     nal,nalmax,indcell,nindiv,npp,nppmax,c,ctmp,
     &     a,ptmp,ftmp,zz,n,ntmp,ploidy,alpha,nudcel)
      implicit none
      integer npop,npopmax,nloc,nlocmax,nal(nlocmax),
     &     nalmax,nindiv,npp,nppmax,indcell(nindiv),
     &     nlocmax2,c(nppmax),ctmp(nppmax),zz(nindiv,nlocmax2),
     &     n(npopmax,nloc,nalmax),ntmp(npopmax,nloc,nalmax),
     &     ploidy,nudcel
      double precision f(npopmax,nlocmax,nalmax),
     &      ftmp(npopmax,nlocmax,nalmax),
     &     a(nalmax),ptmp(nalmax)
      integer ipop,ipp,ipop1,ipop2,iloc,ial
      double precision ggrunif,lrpf,lratio,ratio,llr6
      integer iipp
      integer n1,n2,ntmp1,ntmp2,iud
      double precision junk,termf9bis,gglgamfn,ggrbinom,bern
      double precision alpha,lrp
c      write(*,*) 'begin udcf'
c      write(*,*) 'npop =',npop
c      write(*,*) 'npp=',npp

*     init. tmp. vector of population membership
      do ipp=1,npp
          ctmp(ipp) = c(ipp)
      enddo
      if(nppmax .gt. npp) then
         do ipp=npp+1,nppmax
            ctmp(ipp) = -999
         enddo
      endif

      do iud=1,nudcel
      ipp = 1 + idint(dint(dble(npp)*ggrunif(0.d0,1.d0)))
*     propose new labeling of a tile
         ctmp(ipp) = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
         ipop1 = c(ipp) 
         ipop2 = ctmp(ipp)
*     counting alleles for states c and ctmp
         call countn(nindiv,nlocmax,nlocmax2,npopmax,
     &        nppmax,nal,nalmax,zz,n,indcell,c,ploidy)
         call countn(nindiv,nlocmax,nlocmax2,npopmax,
     &        nppmax,nal,nalmax,zz,ntmp,indcell,ctmp,ploidy)
c         write(*,*) 'in udcf c=',c(1),c(2)
c         write(*,*) '     ctmp=',ctmp(1),ctmp(2)

         bern = 0 
         if(ipop1 .ne. ipop2) then 
            lratio = 0
            do iloc=1,nloc
               n1 = 0
               n2 = 0
               ntmp1 = 0
               ntmp2 = 0
               do ial=1,nal(iloc)
                  lratio = lratio 
     &                 - gglgamfn(alpha+dble(n(ipop1,iloc,ial)))
     &                 - gglgamfn(alpha+dble(n(ipop2,iloc,ial)))
     &                 + gglgamfn(alpha+dble(ntmp(ipop1,iloc,ial)))
     &                 + gglgamfn(alpha+dble(ntmp(ipop2,iloc,ial)))
                  n1 = n1 + n(ipop1,iloc,ial)
                  n2 = n2 + n(ipop2,iloc,ial)
                  ntmp1 = ntmp1 + ntmp(ipop1,iloc,ial)
                  ntmp2 = ntmp2 + ntmp(ipop2,iloc,ial)
               enddo
               lratio = lratio + 
     &              gglgamfn(alpha*dble(nal(iloc))+dble(n1)) 
     &              + gglgamfn(alpha*dble(nal(iloc))+dble(+n2)) 
     &              - gglgamfn(alpha*dble(nal(iloc))+dble(+ntmp1)) 
     &              - gglgamfn(alpha*dble(nal(iloc))+dble(+ntmp2))
            enddo

c            write(*,*) 'in udcf lratio=',lratio
            lratio = dmin1(0.d0,lratio)
            ratio = dexp(lratio)
            bern = ggrbinom(1.d0,ratio)
         endif
         
         if((bern .eq. 1) .or. (ipop1 .eq. ipop2)) then
*     sample new frequencies
c            write(*,*) 'accept move in udcf'
            call samplef(npop,npopmax,nloc,nlocmax,
     &           nal,nalmax,ipop1,ipop2,ftmp,a,ptmp,ntmp,alpha) 
            c(ipp) = ctmp(ipp)
            do iloc=1,nloc
               do ial=1,nal(iloc)
                  f(ipop1,iloc,ial) = ftmp(ipop1,iloc,ial)
                  f(ipop2,iloc,ial) = ftmp(ipop2,iloc,ial)
               enddo
            enddo
         else 
            ctmp(ipp) = c(ipp)
          endif
      enddo
c      write(*,*) 'end udcf'
      end subroutine udcf
***********************************************************************



***********************************************************************
*     joint update of c, f (under UFM) and mean, variance of quant. variables
*     single component update of c
*     new f is proposed according to full conditional pi(f*|c*,data)
      subroutine udcfparvq(npop,npopmax,f,nloc,nlocmax,nlocmax2,
     &     nal,nalmax,indcell,nindiv,npp,nppmax,c,ctmp,nqtc,qtc,
     &     meanqtc,sdqtc,meanqtctmp,sdqtctmp,nnqtc,sqtc,ssqtc,
     &     ksiqtc,kappaqtc,alphaqtc,betaqtc,
     &     a,ptmp,ftmp,zz,n,ntmp,ploidy,alpha,nudcel,usegeno2,useqtc)
      implicit none
      integer npop,npopmax,nloc,nlocmax,nal(nlocmax),
     &     nalmax,nindiv,npp,nppmax,indcell(nindiv),
     &     nlocmax2,c(nppmax),ctmp(nppmax),zz(nindiv,nlocmax2),
     &     n(npopmax,nloc,nalmax),ntmp(npopmax,nloc,nalmax),
     &     ploidy,nudcel,nnqtc,nqtc,usegeno2,useqtc
      double precision f,ftmp,a,ptmp,qtc,
     &     meanqtc,sdqtc,meanqtctmp,sdqtctmp,sqtc,ssqtc,
     &     ksiqtc,kappaqtc,alphaqtc,betaqtc
      dimension f(npopmax,nlocmax,nalmax),
     &     ftmp(npopmax,nlocmax,nalmax),
     &     a(nalmax),ptmp(nalmax),qtc(nindiv,nqtc),
     &     meanqtc(npopmax,nqtc),sdqtc(npopmax,nqtc),
     &     meanqtctmp(npopmax,nqtc),sdqtctmp(npopmax,nqtc),
     &     ksiqtc(nqtc),kappaqtc(nqtc),alphaqtc(nqtc),betaqtc(nqtc), 
     &     nnqtc(npopmax,nqtc),sqtc(npopmax,nqtc),
     &     ssqtc(npopmax,nqtc)
      integer ipop,ipp,ipop1,ipop2,iloc,ial,iqtc,iipp,n1,n2,ntmp1,
     &     ntmp2,iud
      double precision ggrunif,lrpf,lratio,ratio,llr6,contriblr
      double precision junk,termf9bis,gglgamfn,ggrbinom,bern,alpha,lrp,
     &     lratiogq

c      write(*,*) 'begin udcfparvq'
c      write(*,*) 'npop =',npop
c      write(*,*) 'npp=',npp

*     init. tmp. vector of population membership
      do ipp=1,npp
          ctmp(ipp) = c(ipp)
      enddo
      if(nppmax .gt. npp) then
         do ipp=npp+1,nppmax
            ctmp(ipp) = -999
         enddo
      endif

      do iud=1,nudcel
         ipp = 1 + idint(dint(dble(npp)*ggrunif(0.d0,1.d0)))
*     propose new labeling of a tile
         ctmp(ipp) = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
         ipop1 = c(ipp) 
         ipop2 = ctmp(ipp)
         bern = 0   
         if(ipop1 .ne. ipop2) then
            lratio = 0
            if(usegeno2 .eq. 1) then
*     counting alleles for states c and ctmp
               call countn(nindiv,nlocmax,nlocmax2,npopmax,
     &              nppmax,nal,nalmax,zz,n,indcell,c,ploidy)
               call countn(nindiv,nlocmax,nlocmax2,npopmax,
     &              nppmax,nal,nalmax,zz,ntmp,indcell,ctmp,ploidy)
c     write(*,*) 'in udcf c=',c(1),c(2)
c     write(*,*) '     ctmp=',ctmp(1),ctmp(2)
               do iloc=1,nloc
                  n1 = 0
                  n2 = 0
                  ntmp1 = 0
                  ntmp2 = 0
                  do ial=1,nal(iloc)
                     lratio = lratio 
     &                    - gglgamfn(alpha+dble(n(ipop1,iloc,ial)))
     &                    - gglgamfn(alpha+dble(n(ipop2,iloc,ial)))
     &                    + gglgamfn(alpha+dble(ntmp(ipop1,iloc,ial)))
     &                    + gglgamfn(alpha+dble(ntmp(ipop2,iloc,ial)))
                     n1 = n1 + n(ipop1,iloc,ial)
                     n2 = n2 + n(ipop2,iloc,ial)
                     ntmp1 = ntmp1 + ntmp(ipop1,iloc,ial)
                     ntmp2 = ntmp2 + ntmp(ipop2,iloc,ial)
                  enddo
                  lratio = lratio + 
     &                 gglgamfn(alpha*dble(nal(iloc))+dble(n1)) 
     &                 + gglgamfn(alpha*dble(nal(iloc))+dble(+n2)) 
     &                 - gglgamfn(alpha*dble(nal(iloc))+dble(+ntmp1)) 
     &                 - gglgamfn(alpha*dble(nal(iloc))+dble(+ntmp2))
               enddo
            endif
         
            if(useqtc .eq. 1) then 
c$$$               write(*,*) 'c=',c(1),c(2),c(3),c(4)
c$$$               write(*,*) 'ipp=',ipp
c$$$               write(*,*) 'ipop1=',ipop1
c$$$               write(*,*) 'ipop2=',ipop2
*     propose means and variances         
            call propparqvudc(qtc,nindiv,nqtc,indcell,c,ctmp,
     &           meanqtc,sdqtc,meanqtctmp,sdqtctmp,nnqtc,
     &           sqtc,ssqtc,npop,npopmax,nppmax,
     &           ksiqtc,kappaqtc,alphaqtc,betaqtc,ipop1,ipop2,
     &           contriblr)
c            write(*,*) 'meanqtc=',meanqtc
c            write(*,*) 'sdqtc=',sdqtc
c            write(*,*) 'meanqtctmp=',meanqtctmp
c            write(*,*) 'sdqtctmp=',sdqtctmp
*     contrib likelihood
*     argument usegeno2 is passed as 0 to ignore genotypes 
*     hence avoid having likelihood ratio of genotypes twice 
            lratio = lratio + 
     &           lratiogq(zz,f,c,ctmp,indcell,indcell,
     &           npopmax,nlocmax,nalmax,nindiv,nloc,nlocmax2,
     &           nppmax,ploidy,qtc,nqtc,meanqtc,sdqtc,
     &           meanqtctmp,sdqtctmp,0,useqtc)
c            write(*,*) 'in ucdcfparvq lratio=',lratio-junk      
            lratio = lratio + contriblr 
c     write(*,*) 'in ucdcfparvq lratio=',lratio
         endif
c     write(*,*) 'in udcf lratio=',lratio
         lratio = dmin1(0.d0,lratio)
         ratio = dexp(lratio)
         bern = ggrbinom(1.d0,ratio)
      endif
      
      
      if(bern .eq. 1) then
         c(ipp) = ctmp(ipp)
         if(usegeno2 .eq. 1) then 
*     sample new frequencies
c     write(*,*) 'accept move in udcf'
            call samplef(npop,npopmax,nloc,nlocmax,
     &           nal,nalmax,ipop1,ipop2,ftmp,a,ptmp,ntmp,alpha) 
            do iloc=1,nloc
               do ial=1,nal(iloc)
                  f(ipop1,iloc,ial) = ftmp(ipop1,iloc,ial)
                  f(ipop2,iloc,ial) = ftmp(ipop2,iloc,ial)
               enddo
            enddo
         endif
         if(useqtc .eq. 1) then 
            do iqtc = 1,nqtc
               do ipop = 1,npop
                  meanqtc(ipop,iqtc) =  meanqtctmp(ipop,iqtc)
                  sdqtc(ipop,iqtc) =  sdqtctmp(ipop,iqtc)
               enddo
            enddo
         endif
      else 
         ctmp(ipp) = c(ipp)
      endif
      enddo
c      write(*,*) 'end udcf'
      end subroutine udcfparvq
***********************************************************************




***********************************************************************
*     joint update of c, f (under UFM) and mean, variance of quant. variables
*     single component update of c
*     new f is proposed according to full conditional pi(f*|c*,data)
      subroutine udcfallvar2(npop,npopmax,nlocd,nloch,nql,ncolt, 
     &     nal,nalmax,indcell,nindiv,npp,nppmax,c,ctmp,nqtc,qtc,
     &     meanqtc,sdqtc,meanqtctmp,sdqtctmp,nnqtc,sqtc,ssqtc,
     &     ksiqtc,kappaqtc,alphaqtc,betaqtc,a,ptmp,f,ftmp,yy,z,ql,
     &     n,ntmp,alpha,nudcel,usegeno2,usegeno1,useql,useqtc,
     &     ploidy)
      implicit none
      integer npop,npopmax,nlocd,nloch,nql,ncolt,nal(ncolt),
     &     nalmax,nindiv,npp,nppmax,indcell(nindiv),ql(nindiv,nql),
     &     c(nppmax),ctmp(nppmax),yy(nindiv,2*nlocd+2*nloch),
     &     z(nindiv,nloch),n(npopmax,ncolt,nalmax),
     &     ntmp(npopmax,ncolt,nalmax),nudcel,nnqtc,nqtc,usegeno2,
     &     usegeno1,useql,useqtc,ploidy
      double precision f,ftmp,a,ptmp,qtc,
     &     meanqtc,sdqtc,meanqtctmp,sdqtctmp,sqtc,ssqtc,
     &     ksiqtc,kappaqtc,alphaqtc,betaqtc
      dimension f(npopmax,ncolt,nalmax),
     &     ftmp(npopmax,ncolt,nalmax),
     &     a(nalmax),ptmp(nalmax),qtc(nindiv,nqtc),
     &     meanqtc(npopmax,nqtc),sdqtc(npopmax,nqtc),
     &     meanqtctmp(npopmax,nqtc),sdqtctmp(npopmax,nqtc),
     &     nnqtc(npopmax,nqtc),sqtc(npopmax,nqtc),
     &     ssqtc(npopmax,nqtc),
     &     ksiqtc(nqtc),kappaqtc(nqtc),alphaqtc(nqtc),betaqtc(nqtc)
      integer ipop,ipp,ipop1,ipop2,iloc,ial,iqtc,iipp,n1,n2,ntmp1,
     &     ntmp2,iud
      double precision ggrunif,lrpf,lratio,ratio,llr6,contriblr
      double precision junk,termf9bis,gglgamfn,ggrbinom,bern,alpha,lrp,
     &     lrallvar2

c      write(*,*) 'begin udcfallvar2'
c      write(*,*) 'npop =',npop
c      write(*,*) 'npp=',npp

*     init. tmp. vector of population membership
      do ipp=1,npp
          ctmp(ipp) = c(ipp)
      enddo
      if(nppmax .gt. npp) then
         do ipp=npp+1,nppmax
            ctmp(ipp) = -999
         enddo
      endif

      do iud=1,nudcel
         ipp = 1 + idint(dint(dble(npp)*ggrunif(0.d0,1.d0)))
*     propose new labeling of a tile
         ctmp(ipp) = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
         ipop1 = c(ipp) 
         ipop2 = ctmp(ipp)
         bern = 0   
         if(ipop1 .ne. ipop2) then
            lratio = 0
            if(((usegeno2 .eq. 1) .or. (usegeno1 .eq. 1)) .or. 
     &           useql .eq. 1) then
*     counting alleles for states c and ctmp
c$$$               call countn(nindiv,nlocmax,nlocmax2,npopmax,
c$$$  &              nppmax,nal,nalmax,yy,n,indcell,c,ploidy)
c$$$  call countn(nindiv,nlocmax,nlocmax2,npopmax,
c$$$  &              nppmax,nal,nalmax,yy,ntmp,indcell,ctmp,ploidy)
               call countnallvar2(nindiv,nlocd,nloch,nql,ncolt,
     &              nppmax,nal,nalmax,yy,z,ql,n,indcell,c,npopmax,
     &              usegeno2, usegeno1,useql,ploidy)
               call countnallvar2(nindiv,nlocd,nloch,nql,ncolt,
     &              nppmax,nal,nalmax,yy,z,ql,ntmp,indcell,ctmp,npopmax,
     &              usegeno2, usegeno1,useql,ploidy)
            endif
c      write(*,*) 'in udcf c=',c(1),c(2)
c      write(*,*) '     ctmp=',ctmp(1),ctmp(2)
            if(usegeno2 .eq. 1) then 
               do iloc=1,nlocd
                  n1 = 0
                  n2 = 0
                  ntmp1 = 0
                  ntmp2 = 0
                  do ial=1,nal(iloc)
                     lratio = lratio 
     &                    - gglgamfn(alpha+dble(n(ipop1,iloc,ial)))
     &                    - gglgamfn(alpha+dble(n(ipop2,iloc,ial)))
     &                    + gglgamfn(alpha+dble(ntmp(ipop1,iloc,ial)))
     &                    + gglgamfn(alpha+dble(ntmp(ipop2,iloc,ial)))
                     n1 = n1 + n(ipop1,iloc,ial)
                     n2 = n2 + n(ipop2,iloc,ial)
                     ntmp1 = ntmp1 + ntmp(ipop1,iloc,ial)
                     ntmp2 = ntmp2 + ntmp(ipop2,iloc,ial)
                  enddo
                  lratio = lratio + 
     &                 gglgamfn(alpha*dble(nal(iloc))+dble(n1)) 
     &                 + gglgamfn(alpha*dble(nal(iloc))+dble(+n2)) 
     &                 - gglgamfn(alpha*dble(nal(iloc))+dble(+ntmp1)) 
     &                 - gglgamfn(alpha*dble(nal(iloc))+dble(+ntmp2))
               enddo
            endif
            if(usegeno1 .eq. 1) then 
               do iloc=nlocd+1,nlocd+nloch
                  n1 = 0
                  n2 = 0
                  ntmp1 = 0
                  ntmp2 = 0
                  do ial=1,nal(iloc)
                     lratio = lratio 
     &                    - gglgamfn(alpha+dble(n(ipop1,iloc,ial)))
     &                    - gglgamfn(alpha+dble(n(ipop2,iloc,ial)))
     &                    + gglgamfn(alpha+dble(ntmp(ipop1,iloc,ial)))
     &                    + gglgamfn(alpha+dble(ntmp(ipop2,iloc,ial)))
                     n1 = n1 + n(ipop1,iloc,ial)
                     n2 = n2 + n(ipop2,iloc,ial)
                     ntmp1 = ntmp1 + ntmp(ipop1,iloc,ial)
                     ntmp2 = ntmp2 + ntmp(ipop2,iloc,ial)
                  enddo
                  lratio = lratio + 
     &                 gglgamfn(alpha*dble(nal(iloc))+dble(n1)) 
     &                 + gglgamfn(alpha*dble(nal(iloc))+dble(+n2)) 
     &                 - gglgamfn(alpha*dble(nal(iloc))+dble(+ntmp1)) 
     &                 - gglgamfn(alpha*dble(nal(iloc))+dble(+ntmp2))
               enddo
            endif
            if(useql .eq. 1) then 
               do iloc=nlocd+nloch+1,nlocd+nloch+nql
                  n1 = 0
                  n2 = 0
                  ntmp1 = 0
                  ntmp2 = 0
                  do ial=1,nal(iloc)
                     lratio = lratio 
     &                    - gglgamfn(alpha+dble(n(ipop1,iloc,ial)))
     &                    - gglgamfn(alpha+dble(n(ipop2,iloc,ial)))
     &                    + gglgamfn(alpha+dble(ntmp(ipop1,iloc,ial)))
     &                    + gglgamfn(alpha+dble(ntmp(ipop2,iloc,ial)))
                     n1 = n1 + n(ipop1,iloc,ial)
                     n2 = n2 + n(ipop2,iloc,ial)
                     ntmp1 = ntmp1 + ntmp(ipop1,iloc,ial)
                     ntmp2 = ntmp2 + ntmp(ipop2,iloc,ial)
                  enddo
                  lratio = lratio + 
     &                 gglgamfn(alpha*dble(nal(iloc))+dble(n1)) 
     &                 + gglgamfn(alpha*dble(nal(iloc))+dble(+n2)) 
     &                 - gglgamfn(alpha*dble(nal(iloc))+dble(+ntmp1)) 
     &                 - gglgamfn(alpha*dble(nal(iloc))+dble(+ntmp2))
               enddo
            endif
            
            if(useqtc .eq. 1) then 
c$$$  write(*,*) 'c=',c(1),c(2),c(3),c(4)
c$$$  write(*,*) 'ipp=',ipp
c$$$  write(*,*) 'ipop1=',ipop1
c$$$  write(*,*) 'ipop2=',ipop2
*     propose means and variances         
               call propparqvudc(qtc,nindiv,nqtc,indcell,c,ctmp,
     &              meanqtc,sdqtc,meanqtctmp,sdqtctmp,nnqtc,
     &              sqtc,ssqtc,npop,npopmax,nppmax,
     &              ksiqtc,kappaqtc,alphaqtc,betaqtc,ipop1,ipop2,
     &              contriblr)
c     write(*,*) 'meanqtc=',meanqtc
c     write(*,*) 'sdqtc=',sdqtc
c     write(*,*) 'meanqtctmp=',meanqtctmp
c     write(*,*) 'sdqtctmp=',sdqtctmp
*     contrib likelihood
*     argument usegeno2,usegeno1,useql passed as 0 to disregard these variables 
*     hence avoid having their contrib to likelihood ratio  twice 
c$$$  lratio = lratio + 
c$$$  &           lratiogq(yy,f,c,ctmp,indcell,indcell,
c$$$  &           npopmax,nlocmax,nalmax,nindiv,nloc,nlocmax2,
c$$$  &           nppmax,ploidy,qtc,nqtc,meanqtc,sdqtc,
c$$$  &           meanqtctmp,sdqtctmp,0,useqtc)
c               write(*,*) 'c=',c(1),c(2)
c               write(*,*) 'ctmp=',ctmp(1),ctmp(2)
c               write(*,*) 'in ucdcfparvq lratio=',lratio
               lratio = lratio + lrallvar2(yy,z,ql,
     &              f,c,ctmp,indcell,indcell,
     &              npopmax,nlocd,nloch,nql,ncolt,nalmax,nindiv,
     &              nppmax,qtc,nqtc,meanqtc,sdqtc,meanqtctmp,sdqtctmp,
     &              0,0,0,useqtc,ploidy)
c                write(*,*) 'in ucdcfparvq lratio=',lratio
c     write(*,*) 'in ucdcfparvq lratio=',lratio-junk      
               lratio = lratio + contriblr 
c               write(*,*) 'in ucdcfparvq lratio=',lratio
            endif
c      write(*,*) 'in udcf lratio=',lratio
            lratio = dmin1(0.d0,lratio)
            ratio = dexp(lratio)
            bern = ggrbinom(1.d0,ratio)
c            write(*,*) 'in udcf ratio=',ratio
         endif
         
         
         if(bern .eq. 1) then
            c(ipp) = ctmp(ipp)
            if(((usegeno2 .eq. 1) .or. (usegeno1 .eq. 1)) .or. 
     &           useql .eq. 1) then
*     sample new frequencies
c     write(*,*) 'accept move in udcf'
c$$$               call samplef(npop,npopmax,nloc,nlocmax,
c$$$     &              nal,nalmax,ipop1,ipop2,ftmp,a,ptmp,ntmp,alpha) 
               call samplefallvar(npopmax,nlocd,nloch,nql,ncolt,
     &              nal,nalmax,ipop1,ipop2,ftmp,a,ptmp,ntmp,alpha,
     &              usegeno2,usegeno1,useql)
               do iloc=1,ncolt
                  do ial=1,nal(iloc)
                     f(ipop1,iloc,ial) = ftmp(ipop1,iloc,ial)
                     f(ipop2,iloc,ial) = ftmp(ipop2,iloc,ial)
                  enddo
               enddo
            endif
            if(useqtc .eq. 1) then 
               do iqtc = 1,nqtc
                  do ipop = 1,npop
                     meanqtc(ipop,iqtc) =  meanqtctmp(ipop,iqtc)
                     sdqtc(ipop,iqtc) =  sdqtctmp(ipop,iqtc)
                  enddo
               enddo
            endif
         else 
            ctmp(ipp) = c(ipp)
         endif
      enddo
c      write(*,*) 'end udcf'
      end subroutine udcfallvar2
***********************************************************************






***********************************************************************
*     joint update of c, f (under UFM) and mean, variance of quant. variables
*     single component update of c
*     new f  proposed according to full conditional 
*     new mean and variances proposed according to full conditional 
*     prior mean-variance: Normal Scaled Inverse Gamma
      subroutine udcfallvar3(npop,npopmax,nlocd,nloch,nql,ncolt, 
     &     nal,nalmax,indcell,nindiv,npp,nppmax,c,ctmp,nqtc,qtc,
     &     meanqtc,sdqtc,meanqtctmp,sdqtctmp,nnqtc,sqtc,ssqtc,
     &     ksiqtc,kappaqtc,alphaqtc,betaqtc,a,ptmp,f,ftmp,yy,z,ql,
     &     n,ntmp,alpha,nudcel,usegeno2,usegeno1,useql,useqtc,
     &     ploidy)
      implicit none
      integer npop,npopmax,nlocd,nloch,nql,ncolt,nal(ncolt),
     &     nalmax,nindiv,npp,nppmax,indcell(nindiv),ql(nindiv,nql),
     &     c(nppmax),ctmp(nppmax),yy(nindiv,2*nlocd+2*nloch),
     &     z(nindiv,nloch),n(npopmax,ncolt,nalmax),
     &     ntmp(npopmax,ncolt,nalmax),nudcel,nnqtc,nqtc,usegeno2,
     &     usegeno1,useql,useqtc,ploidy
      double precision f,ftmp,a,ptmp,qtc,
     &     meanqtc,sdqtc,meanqtctmp,sdqtctmp,sqtc,ssqtc,
     &     ksiqtc,kappaqtc,alphaqtc,betaqtc
      dimension f(npopmax,ncolt,nalmax),
     &     ftmp(npopmax,ncolt,nalmax),
     &     a(nalmax),ptmp(nalmax),qtc(nindiv,nqtc),
     &     meanqtc(npopmax,nqtc),sdqtc(npopmax,nqtc),
     &     meanqtctmp(npopmax,nqtc),sdqtctmp(npopmax,nqtc),
     &     nnqtc(npopmax,nqtc),sqtc(npopmax,nqtc),
     &     ssqtc(npopmax,nqtc),
     &     ksiqtc(nqtc),kappaqtc(nqtc),alphaqtc(nqtc),betaqtc(nqtc)
      integer ipop,ipp,ipop1,ipop2,iloc,ial,iqtc,iipp,n1,n2,ntmp1,
     &     ntmp2,iud,i,j
      double precision ggrunif,lrpf,lratio,ratio,llr6,contriblr
      double precision junk,termf9bis,gglgamfn,ggrbinom,bern,alpha,lrp,
     &     lrallvar2

c$$$      write(*,*) 'begin udcfallvar3'
c$$$      write(*,*) 'npop =',npop
c$$$      write(*,*) 'npp=',npp

*     init. tmp. vector of population membership
      do ipp=1,npp
          ctmp(ipp) = c(ipp)
      enddo
      if(nppmax .gt. npp) then
         do ipp=npp+1,nppmax
            ctmp(ipp) = -999
         enddo
      endif

      if(useqtc .eq. 1) then
         do ipop=1,npop
            do iqtc=1,nqtc
               meanqtctmp(ipop,iqtc) = meanqtc(ipop,iqtc)
               sdqtctmp(ipop,iqtc) = sdqtc(ipop,iqtc)
            enddo
         enddo
      endif

      do iud=1,nudcel
         ipp = 1 + idint(dint(dble(npp)*ggrunif(0.d0,1.d0)))
*     propose new labeling of a tile
         ctmp(ipp) = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
         ipop1 = c(ipp) 
         ipop2 = ctmp(ipp)
         bern = 0   
         if(ipop1 .ne. ipop2) then
            lratio = 0
            if(((usegeno2 .eq. 1) .or. (usegeno1 .eq. 1)) .or. 
     &           useql .eq. 1) then
*     counting alleles for states c and ctmp
c$$$               call countn(nindiv,nlocmax,nlocmax2,npopmax,
c$$$  &              nppmax,nal,nalmax,yy,n,indcell,c,ploidy)
c$$$  call countn(nindiv,nlocmax,nlocmax2,npopmax,
c$$$  &              nppmax,nal,nalmax,yy,ntmp,indcell,ctmp,ploidy)
               call countnallvar2(nindiv,nlocd,nloch,nql,ncolt,
     &              nppmax,nal,nalmax,yy,z,ql,n,indcell,c,npopmax,
     &              usegeno2, usegeno1,useql,ploidy)
               call countnallvar2(nindiv,nlocd,nloch,nql,ncolt,
     &              nppmax,nal,nalmax,yy,z,ql,ntmp,indcell,ctmp,npopmax,
     &              usegeno2, usegeno1,useql,ploidy)
            endif
 
            if(usegeno2 .eq. 1) then 
               do iloc=1,nlocd
                  n1 = 0
                  n2 = 0
                  ntmp1 = 0
                  ntmp2 = 0
                  do ial=1,nal(iloc)
                     lratio = lratio 
     &                    - gglgamfn(alpha+dble(n(ipop1,iloc,ial)))
     &                    - gglgamfn(alpha+dble(n(ipop2,iloc,ial)))
     &                    + gglgamfn(alpha+dble(ntmp(ipop1,iloc,ial)))
     &                    + gglgamfn(alpha+dble(ntmp(ipop2,iloc,ial)))
                     n1 = n1 + n(ipop1,iloc,ial)
                     n2 = n2 + n(ipop2,iloc,ial)
                     ntmp1 = ntmp1 + ntmp(ipop1,iloc,ial)
                     ntmp2 = ntmp2 + ntmp(ipop2,iloc,ial)
                  enddo
                  lratio = lratio + 
     &                 gglgamfn(alpha*dble(nal(iloc))+dble(n1)) 
     &                 + gglgamfn(alpha*dble(nal(iloc))+dble(+n2)) 
     &                 - gglgamfn(alpha*dble(nal(iloc))+dble(+ntmp1)) 
     &                 - gglgamfn(alpha*dble(nal(iloc))+dble(+ntmp2))
               enddo
            endif
            if(usegeno1 .eq. 1) then 
               do iloc=nlocd+1,nlocd+nloch
                  n1 = 0
                  n2 = 0
                  ntmp1 = 0
                  ntmp2 = 0
                  do ial=1,nal(iloc)
                     lratio = lratio 
     &                    - gglgamfn(alpha+dble(n(ipop1,iloc,ial)))
     &                    - gglgamfn(alpha+dble(n(ipop2,iloc,ial)))
     &                    + gglgamfn(alpha+dble(ntmp(ipop1,iloc,ial)))
     &                    + gglgamfn(alpha+dble(ntmp(ipop2,iloc,ial)))
                     n1 = n1 + n(ipop1,iloc,ial)
                     n2 = n2 + n(ipop2,iloc,ial)
                     ntmp1 = ntmp1 + ntmp(ipop1,iloc,ial)
                     ntmp2 = ntmp2 + ntmp(ipop2,iloc,ial)
                  enddo
                  lratio = lratio + 
     &                 gglgamfn(alpha*dble(nal(iloc))+dble(n1)) 
     &                 + gglgamfn(alpha*dble(nal(iloc))+dble(+n2)) 
     &                 - gglgamfn(alpha*dble(nal(iloc))+dble(+ntmp1)) 
     &                 - gglgamfn(alpha*dble(nal(iloc))+dble(+ntmp2))
               enddo
            endif
            if(useql .eq. 1) then 
               do iloc=nlocd+nloch+1,nlocd+nloch+nql
                  n1 = 0
                  n2 = 0
                  ntmp1 = 0
                  ntmp2 = 0
                  do ial=1,nal(iloc)
                     lratio = lratio 
     &                    - gglgamfn(alpha+dble(n(ipop1,iloc,ial)))
     &                    - gglgamfn(alpha+dble(n(ipop2,iloc,ial)))
     &                    + gglgamfn(alpha+dble(ntmp(ipop1,iloc,ial)))
     &                    + gglgamfn(alpha+dble(ntmp(ipop2,iloc,ial)))
                     n1 = n1 + n(ipop1,iloc,ial)
                     n2 = n2 + n(ipop2,iloc,ial)
                     ntmp1 = ntmp1 + ntmp(ipop1,iloc,ial)
                     ntmp2 = ntmp2 + ntmp(ipop2,iloc,ial)
                  enddo
                  lratio = lratio + 
     &                 gglgamfn(alpha*dble(nal(iloc))+dble(n1)) 
     &                 + gglgamfn(alpha*dble(nal(iloc))+dble(+n2)) 
     &                 - gglgamfn(alpha*dble(nal(iloc))+dble(+ntmp1)) 
     &                 - gglgamfn(alpha*dble(nal(iloc))+dble(+ntmp2))
               enddo
            endif
            
            if(useqtc .eq. 1) then 
c$$$  write(*,*) 'c=',c(1),c(2),c(3),c(4)
c$$$  write(*,*) 'ipp=',ipp
c$$$  write(*,*) 'ipop1=',ipop1
c$$$  write(*,*) 'ipop2=',ipop2
*     propose means and variances         
c$$$      write(*,*) 'udcfallvar3 meanqtc=',sdqtc
c$$$      write(*,*) 'meanqtctmp=',sdqtctmp
c$$$      write(*,*) 'sdqtc=',sdqtc
c$$$      write(*,*) 'sdqtctmp=',sdqtctmp
c$$$      write(*,*) 'ksiqtc=',ksiqtc
c$$$      write(*,*) 'kappaqtc=',kappaqtc
c$$$      write(*,*) 'alphaqtc=',alphaqtc
c$$$      write(*,*) 'betaqtc=',betaqtc
c               write(*,*) 'lratio=',lratio
c$$$               write(*,*) 'in udcfallvar3 c(ipp)=',c(ipp)
c$$$               write(*,*) 'in udcfallvar3 ctmp(ipp)=',ctmp(ipp)
c$$$               write(*,*) 'in udcfallvar3  ipop1,ipop2=',ipop1,ipop2
               call propparqvudc3(qtc,nindiv,nqtc,indcell,c,ctmp,
     &              meanqtc,sdqtc,meanqtctmp,sdqtctmp,nnqtc,
     &              sqtc,ssqtc,npop,npopmax,nppmax,
     &              ksiqtc,kappaqtc,alphaqtc,betaqtc,ipop1,ipop2,
     &              contriblr)
               lratio = lratio + contriblr 
c               write(*,*) 'udcfallvar3 aft propparqvucd3 lratio=',lratio
*     contrib likelihood
*     argument usegeno2,usegeno1,useql passed as 0 to disregard these variables 
*     hence avoid having their contrib to likelihood ratio  twice 
               lratio = lratio + lrallvar2(yy,z,ql,
     &              f,c,ctmp,indcell,indcell,
     &              npopmax,nlocd,nloch,nql,ncolt,nalmax,nindiv,
     &              nppmax,qtc,nqtc,meanqtc,sdqtc,meanqtctmp,sdqtctmp,
     &              0,0,0,useqtc,ploidy)
c               write(*,*) 'lratio=',lratio
c$$$               write(6,*) 'udcfallvar3  apres propparqvudc3='
c$$$               write(*,*) 'c=',c(1),c(2),c(3),c(4)
c$$$               write(*,*) 'ctmp=',ctmp(1),ctmp(2),ctmp(3),ctmp(4)   
c$$$               do i=1,npop
c$$$                  write(*,*) 'meanqtc=',(meanqtc(i,j),j=1,nqtc)
c$$$               enddo
c$$$               do i=1,npop
c$$$                  write(*,*) 'sdqtc=',(sdqtc(i,j),j=1,nqtc)
c$$$               enddo
c$$$               do i=1,npop
c$$$                  write(*,*) 'meanqtctmp=',(meanqtctmp(i,j),j=1,nqtc)
c$$$               enddo
c$$$               do i=1,npop
c$$$                  write(*,*) 'sdqtctmp=',(meanqtctmp(i,j),j=1,nqtc)
c$$$               enddo
c               write(*,*) 'udcfallvar3 after lrallvar2 lratio=',lratio
            endif
c$$$      write(*,*) 'in udcfallvar3 lratio=',lratio
            lratio = dmin1(0.d0,lratio)
            ratio = dexp(lratio)
            bern = ggrbinom(1.d0,ratio)
c            write(*,*) 'in udcf ratio=',ratio
         endif
         
         
         if(bern .eq. 1) then
            c(ipp) = ctmp(ipp)
            if(((usegeno2 .eq. 1) .or. (usegeno1 .eq. 1)) .or. 
     &           useql .eq. 1) then
*     sample new frequencies
c     write(*,*) 'accept move in udcf'
c$$$               call samplef(npop,npopmax,nloc,nlocmax,
c$$$     &              nal,nalmax,ipop1,ipop2,ftmp,a,ptmp,ntmp,alpha) 
               call samplefallvar(npopmax,nlocd,nloch,nql,ncolt,
     &              nal,nalmax,ipop1,ipop2,ftmp,a,ptmp,ntmp,alpha,
     &              usegeno2,usegeno1,useql)
               do iloc=1,ncolt
                  do ial=1,nal(iloc)
                     f(ipop1,iloc,ial) = ftmp(ipop1,iloc,ial)
                     f(ipop2,iloc,ial) = ftmp(ipop2,iloc,ial)
                  enddo
               enddo
            endif
            if(useqtc .eq. 1) then 
               do iqtc = 1,nqtc
                     meanqtc(ipop1,iqtc) =  meanqtctmp(ipop1,iqtc)
                     meanqtc(ipop2,iqtc) =  meanqtctmp(ipop2,iqtc)
                     sdqtc(ipop1,iqtc) =  sdqtctmp(ipop1,iqtc)
                     sdqtc(ipop2,iqtc) =  sdqtctmp(ipop2,iqtc)
               enddo
            endif
         else 
            ctmp(ipp) = c(ipp)
         endif
      enddo
c      write(*,*) 'end udcfallvar3'
      end subroutine udcfallvar3
***********************************************************************










***********************************************************************
*     joint update of c, f (under UFM) and mean, variance of quant. variables
*     multiple component update of c
*     new f  proposed according to full conditional 
*     new mean and variances proposed according to full conditional 
*     prior mean-variance: Normal Scaled Inverse Gamma
      subroutine udcfallvar4(npop,npopmax,nlocd,nloch,nql,ncolt, 
     &     nal,nalmax,indcell,nindiv,npp,nppmax,c,ctmp,nqtc,qtc,
     &     meanqtc,sdqtc,meanqtctmp,sdqtctmp,nnqtc,sqtc,ssqtc,
     &     ksiqtc,kappaqtc,alphaqtc,betaqtc,a,ptmp,f,ftmp,yy,z,ql,
     &     n,ntmp,alpha,nudcel,usegeno2,usegeno1,useql,useqtc,
     &     ploidy)
      implicit none
      integer npop,npopmax,nlocd,nloch,nql,ncolt,nal(ncolt),
     &     nalmax,nindiv,npp,nppmax,indcell(nindiv),ql(nindiv,nql),
     &     c(nppmax),ctmp(nppmax),yy(nindiv,2*nlocd+2*nloch),
     &     z(nindiv,nloch),n(npopmax,ncolt,nalmax),
     &     ntmp(npopmax,ncolt,nalmax),nudcel,nnqtc,nqtc,usegeno2,
     &     usegeno1,useql,useqtc,ploidy
      double precision f,ftmp,a,ptmp,qtc,
     &     meanqtc,sdqtc,meanqtctmp,sdqtctmp,sqtc,ssqtc,
     &     ksiqtc,kappaqtc,alphaqtc,betaqtc
      dimension f(npopmax,ncolt,nalmax),
     &     ftmp(npopmax,ncolt,nalmax),
     &     a(nalmax),ptmp(nalmax),qtc(nindiv,nqtc),
     &     meanqtc(npopmax,nqtc),sdqtc(npopmax,nqtc),
     &     meanqtctmp(npopmax,nqtc),sdqtctmp(npopmax,nqtc),
     &     nnqtc(npopmax,nqtc),sqtc(npopmax,nqtc),
     &     ssqtc(npopmax,nqtc),
     &     ksiqtc(nqtc),kappaqtc(nqtc),alphaqtc(nqtc),betaqtc(nqtc)
      integer ipop,ipp,iloc,ial,iqtc,iipp,n1,n2,ntmp1,
     &     ntmp2,iud,i,j,nnudcel
      double precision ggrunif,lrpf,lratio,ratio,llr6,contriblr
      double precision junk,termf9bis,gglgamfn,ggrbinom,bern,alpha,lrp,
     &     lrallvar2
*     init. tmp. vector of population membership
      do ipp=1,npp
          ctmp(ipp) = c(ipp)
      enddo
      if(nppmax .gt. npp) then
         do ipp=npp+1,nppmax
            ctmp(ipp) = -999
         enddo
      endif
      if(useqtc .eq. 1) then
         do ipop=1,npop
            do iqtc=1,nqtc
               meanqtctmp(ipop,iqtc) = meanqtc(ipop,iqtc)
               sdqtctmp(ipop,iqtc) = sdqtc(ipop,iqtc)
            enddo
         enddo
      endif

*     propose new labeling of some tiles
      nnudcel = 2 + idint(dint(dble(nudcel-1)*ggrunif(0.d0,1.d0)))
      do iud=1,nnudcel
         ipp = 1 + idint(dint(dble(npp)*ggrunif(0.d0,1.d0)))
         ctmp(ipp) = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
      enddo
      bern = 0   
      lratio = 0
      if(((usegeno2 .eq. 1) .or. (usegeno1 .eq. 1)) .or. 
     &     useql .eq. 1) then
*     counting alleles for states c and ctmp
         call countnallvar2(nindiv,nlocd,nloch,nql,ncolt,
     &        nppmax,nal,nalmax,yy,z,ql,n,indcell,c,npopmax,
     &        usegeno2, usegeno1,useql,ploidy)
         call countnallvar2(nindiv,nlocd,nloch,nql,ncolt,
     &        nppmax,nal,nalmax,yy,z,ql,ntmp,indcell,ctmp,npopmax,
     &        usegeno2, usegeno1,useql,ploidy)
      endif
 
      do ipop=1,npop
         if(usegeno2 .eq. 1) then 
            do iloc=1,nlocd
               n1 = 0
               ntmp1 = 0
               do ial=1,nal(iloc)
                  lratio = lratio 
     &                 - gglgamfn(alpha+dble(n(ipop,iloc,ial)))
     &                 + gglgamfn(alpha+dble(ntmp(ipop,iloc,ial)))
                  n1 = n1 + n(ipop,iloc,ial)
                  ntmp1 = ntmp1 + ntmp(ipop,iloc,ial)
               enddo
               lratio = lratio + 
     &              gglgamfn(alpha*dble(nal(iloc))+dble(n1)) 
     &              - gglgamfn(alpha*dble(nal(iloc))+dble(+ntmp1)) 
            enddo
         endif
         if(usegeno1 .eq. 1) then 
            do iloc=nlocd+1,nlocd+nloch
               n1 = 0
               ntmp1 = 0
               do ial=1,nal(iloc)
                  lratio = lratio 
     &                 - gglgamfn(alpha+dble(n(ipop,iloc,ial)))
     &                 + gglgamfn(alpha+dble(ntmp(ipop,iloc,ial)))
                  n1 = n1 + n(ipop,iloc,ial)
                  ntmp1 = ntmp1 + ntmp(ipop,iloc,ial)
               enddo
               lratio = lratio + 
     &              gglgamfn(alpha*dble(nal(iloc))+dble(n1)) 
     &              - gglgamfn(alpha*dble(nal(iloc))+dble(+ntmp1)) 
            enddo
         endif
         if(useql .eq. 1) then 
            do iloc=nlocd+nloch+1,nlocd+nloch+nql
               n1 = 0
               ntmp1 = 0
               do ial=1,nal(iloc)
                  lratio = lratio 
     &                 - gglgamfn(alpha+dble(n(ipop,iloc,ial)))
     &                 + gglgamfn(alpha+dble(ntmp(ipop,iloc,ial)))
                  n1 = n1 + n(ipop,iloc,ial)
                  ntmp1 = ntmp1 + ntmp(ipop,iloc,ial)
               enddo
               lratio = lratio + 
     &              gglgamfn(alpha*dble(nal(iloc))+dble(n1)) 
     &              - gglgamfn(alpha*dble(nal(iloc))+dble(+ntmp1)) 
            enddo
         endif
      enddo
            
      if(useqtc .eq. 1) then 
         call propparqvudc4(qtc,nindiv,nqtc,indcell,c,ctmp,
     &        meanqtc,sdqtc,meanqtctmp,sdqtctmp,nnqtc,
     &        sqtc,ssqtc,npop,npopmax,nppmax,
     &        ksiqtc,kappaqtc,alphaqtc,betaqtc,contriblr)
         lratio = lratio + contriblr 
         lratio = lratio + lrallvar2(yy,z,ql,
     &        f,c,ctmp,indcell,indcell,
     &        npopmax,nlocd,nloch,nql,ncolt,nalmax,nindiv,
     &        nppmax,qtc,nqtc,meanqtc,sdqtc,meanqtctmp,sdqtctmp,
     &              0,0,0,useqtc,ploidy)
      endif
       if(lratio .gt. 0) then 
c          write(*,*) 'lratio=',lratio
       endif
      lratio = dmin1(0.d0,lratio)
     
      ratio = dexp(lratio)
      bern = ggrbinom(1.d0,ratio)
         
      if(bern .eq. 1) then
         do ipp=1,npp
            c(ipp) = ctmp(ipp)
         enddo
         if(((usegeno2 .eq. 1) .or. (usegeno1 .eq. 1)) .or. 
     &           useql .eq. 1) then
*     sample new frequencies
            call samplefallvar4(npop,npopmax,nlocd,nloch,nql,ncolt,
     &           nal,nalmax,ftmp,a,ptmp,ntmp,alpha,
     &           usegeno2,usegeno1,useql)
            do ipop=1,npop
               do iloc=1,ncolt
                  do ial=1,nal(iloc)
                     f(ipop,iloc,ial) = ftmp(ipop,iloc,ial)
                  enddo
               enddo
            enddo
         endif
         if(useqtc .eq. 1) then 
            do ipop=1,npop
               do iqtc = 1,nqtc
                  meanqtc(ipop,iqtc) =  meanqtctmp(ipop,iqtc)
                  sdqtc(ipop,iqtc) =  sdqtctmp(ipop,iqtc)
               enddo
            enddo
         endif
      else 
         do ipp=1,npp
            ctmp(ipp) = c(ipp)
         enddo
      endif
         
c      write(*,*) 'end udcfallvar3'
      end subroutine udcfallvar4
***********************************************************************








***********************************************************************
*     joint update of c and f under the CFM
*     single component update of c
*     new f is proposed according to full conditional pi(f*|c*,zz)
      subroutine udcf2(npop,npopmax,f,fa,drift,
     &     nloc,nlocmax,nlocmax2,
     &     nal,nalmax,indcell,nindiv,npp,nppmax,c,ctmp,
     &     a,ptmp,ftmp,zz,n,ntmp,ploidy,nudcel)
      implicit none
      integer npop,npopmax,nloc,nlocmax,nal(nlocmax),
     &     nalmax,nindiv,npp,nppmax,indcell(nindiv),
     &     nlocmax2,c(nppmax),ctmp(nppmax),zz(nindiv,nlocmax2),
     &     n(npopmax,nloc,nalmax),ntmp(npopmax,nloc,nalmax),
     &     ploidy,nudcel
      double precision f(npopmax,nlocmax,nalmax),drift(npopmax),
     &      ftmp(npopmax,nlocmax,nalmax),
     &     a(nalmax),ptmp(nalmax),fa(nlocmax,nalmax)
      integer ipop,ipp,ipop1,ipop2,iloc,ial,iud
      double precision alpha,ggrunif,lratio,lTf
      double precision bern,ggrbinom

*     init. tmp. vector of population membership
      do ipp=1,npp
          ctmp(ipp) = c(ipp)
      enddo
      if(nppmax .gt. npp) then
         do ipp=npp+1,nppmax
            ctmp(ipp) = -999
         enddo
      endif

      do iud=1,nudcel
      ipp = 1 + idint(dint(dble(npp)*ggrunif(0.d0,1.d0)))
c         write(*,*) 'ipp=',ipp
*     propose new labeling of a tile
         ctmp(ipp) = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
         ipop1 = c(ipp)
         ipop2 = ctmp(ipp)
*     counting alleles for both states of c
         call countn(nindiv,nlocmax,nlocmax2,npopmax,
     &        nppmax,nal,nalmax,zz,n,indcell,c,ploidy)
         call countn(nindiv,nlocmax,nlocmax2,npopmax,
     &        nppmax,nal,nalmax,zz,ntmp,indcell,ctmp,ploidy)
         
         bern = 0
         if(ipop1 .ne. ipop2) then 
            lratio = lTf(ipop1,ntmp,fa,drift,npopmax,nloc,nal,nalmax)
     &           + lTf(ipop2,ntmp,fa,drift,npopmax,nloc,nal,nalmax)
     &           - lTf(ipop1,n,fa,drift,npopmax,nloc,nal,nalmax)
     &           - lTf(ipop2,n,fa,drift,npopmax,nloc,nal,nalmax)
            
            lratio = dmin1(0.d0,lratio)
            alpha = dexp(lratio)
            bern = ggrbinom(1.d0,alpha)
         endif
c         write(*,*) 'bern=',bern

         if((bern .eq. 1) .or. (ipop1 .eq. ipop2)) then
            call samplef2(npop,npopmax,nloc,nlocmax,
     &           nal,nalmax,ipop1,ipop2,f,ftmp,
     &           fa,drift,a,ptmp,ntmp)
            c(ipp) = ctmp(ipp)
            do iloc=1,nloc
               do ial=1,nal(iloc)
                  f(ipop1,iloc,ial) = ftmp(ipop1,iloc,ial)
                  f(ipop2,iloc,ial) = ftmp(ipop2,iloc,ial)
               enddo
            enddo
         else 
            ctmp(ipp) = c(ipp)
            do iloc=1,nloc
               do ial=1,nal(iloc)
                  ftmp(ipop1,iloc,ial) = f(ipop1,iloc,ial)
                  ftmp(ipop2,iloc,ial) = f(ipop2,iloc,ial)
               enddo
            enddo
         endif        
      enddo
      end subroutine udcf2
************************************************************************



c$$$
c$$$***********************************************************************
c$$$*     joint update of c, f (under the CFM) 
c$$$*     and mean, variance of quant. variables
c$$$*     single component update of c
c$$$*     new f is proposed according to full conditional pi(f*|c*,zz)
c$$$      subroutine udcf2parvq(npop,npopmax,f,fa,drift,
c$$$     &     nloc,nlocmax,nlocmax2,
c$$$     &     nal,nalmax,indcell,nindiv,npp,nppmax,c,ctmp,nqtc,qtc,
c$$$     &     meanqtc,sdqtc,meanqtctmp,sdqtctmp,nnqtc,sqtc,ssqtc,
c$$$     &     ksiqtc,kappaqtc,alphaqtc,betaqtc,
c$$$     &     a,ptmp,ftmp,zz,n,ntmp,ploidy,nudcel,usegeno2,useqtc)
c$$$      implicit none
c$$$      integer npop,npopmax,nloc,nlocmax,nal(nlocmax),
c$$$     &     nalmax,nindiv,npp,nppmax,indcell(nindiv),
c$$$     &     nlocmax2,c(nppmax),ctmp(nppmax),zz(nindiv,nlocmax2),
c$$$     &     n(npopmax,nloc,nalmax),ntmp(npopmax,nloc,nalmax),
c$$$     &     ploidy,nudcel,nnqtc,nqtc,usegeno2,useqtc
c$$$      double precision f,drift,ftmp,a,ptmp,fa,qtc,
c$$$     &     meanqtc,sdqtc,meanqtctmp,sdqtctmp,sqtc,ssqtc,
c$$$     &     ksiqtc,kappaqtc,alphaqtc,betaqtc
c$$$      dimension f(npopmax,nlocmax,nalmax),drift(npopmax),
c$$$     &      ftmp(npopmax,nlocmax,nalmax),
c$$$     &     a(nalmax),ptmp(nalmax),fa(nlocmax,nalmax),
c$$$     &     qtc(nindiv,nqtc),meanqtc(npopmax,nqtc),
c$$$     &     sdqtc(npopmax,nqtc),ksiqtc(nqtc),
c$$$     &     meanqtctmp(npopmax,nqtc),sdqtctmp(npopmax,nqtc),
c$$$     &     nnqtc(npopmax,nqtc),sqtc(npopmax,nqtc),
c$$$     &     ssqtc(npopmax,nqtc)
c$$$      integer ipop,ipp,ipop1,ipop2,iloc,ial,iud,iqtc
c$$$      double precision alpha,ggrunif,lratio,lTf,contriblr,bern,ggrbinom,
c$$$     &     lratiogq
c$$$
c$$$*     init. tmp. vector of population membership
c$$$      do ipp=1,npp
c$$$          ctmp(ipp) = c(ipp)
c$$$      enddo
c$$$      if(nppmax .gt. npp) then
c$$$         do ipp=npp+1,nppmax
c$$$            ctmp(ipp) = -999
c$$$         enddo
c$$$      endif
c$$$
c$$$      do iud=1,nudcel
c$$$      ipp = 1 + idint(dint(dble(npp)*ggrunif(0.d0,1.d0)))
c$$$c         write(*,*) 'ipp=',ipp
c$$$*     propose new labeling of a tile
c$$$         ctmp(ipp) = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
c$$$         ipop1 = c(ipp)
c$$$         ipop2 = ctmp(ipp)
c$$$         bern = 0
c$$$         if(ipop1 .ne. ipop2) then 
c$$$            lratio = 0
c$$$            if(usegeno2 .eq. 1) then
c$$$*     counting alleles for both states of c
c$$$            call countn(nindiv,nlocmax,nlocmax2,npopmax,
c$$$     &           nppmax,nal,nalmax,zz,n,indcell,c,ploidy)
c$$$            call countn(nindiv,nlocmax,nlocmax2,npopmax,
c$$$     &           nppmax,nal,nalmax,zz,ntmp,indcell,ctmp,ploidy)
c$$$               lratio = lratio + 
c$$$     &              lTf(ipop1,ntmp,fa,drift,npopmax,nloc,nal,nalmax)
c$$$     &              + lTf(ipop2,ntmp,fa,drift,npopmax,nloc,nal,nalmax)
c$$$     &              - lTf(ipop1,n,fa,drift,npopmax,nloc,nal,nalmax)
c$$$     &              - lTf(ipop2,n,fa,drift,npopmax,nloc,nal,nalmax)
c$$$            endif
c$$$            if(useqtc .eq. 1) then 
c$$$*     propose means and variances 
c$$$               call propparqvudc(qtc,nindiv,nqtc,indcell,c,ctmp,
c$$$     &           meanqtc,sdqtc,meanqtctmp,sdqtctmp,nnqtc,
c$$$     &              sqtc,ssqtc,npop,npopmax,nppmax,
c$$$     &              ksiqtc,kappaqtc,alphaqtc,betaqtc,ipop1,ipop2,
c$$$     &              contriblr)
c$$$*     contrib likelihood
c$$$*     argument usegeno2 is passed as 0 to ignore genotypes 
c$$$*     hence avoid having likelihood ratio of genotypes twice 
c$$$               lratio = lratio + 
c$$$     &              lratiogq(zz,f,c,ctmp,indcell,indcell,
c$$$     &              npopmax,nlocmax,nalmax,nindiv,nloc,nlocmax2,
c$$$     &              nppmax,ploidy,qtc,nqtc,meanqtc,sdqtc,
c$$$     &              meanqtctmp,sdqtctmp,0,useqtc)
c$$$               
c$$$               lratio = lratio + contriblr 
c$$$            endif
c$$$            lratio = dmin1(0.d0,lratio)
c$$$            alpha = dexp(lratio)
c$$$            bern = ggrbinom(1.d0,alpha)
c$$$         endif
c$$$            
c$$$c     write(*,*) 'bern=',bern
c$$$            
c$$$         if(bern .eq. 1) then
c$$$            c(ipp) = ctmp(ipp)
c$$$            if(usegeno2 .eq. 1) then 
c$$$               call samplef2(npop,npopmax,nloc,nlocmax,
c$$$     &              nal,nalmax,ipop1,ipop2,f,ftmp,
c$$$     &              fa,drift,a,ptmp,ntmp)
c$$$               do iloc=1,nloc
c$$$                  do ial=1,nal(iloc)
c$$$                     f(ipop1,iloc,ial) = ftmp(ipop1,iloc,ial)
c$$$                     f(ipop2,iloc,ial) = ftmp(ipop2,iloc,ial)
c$$$                  enddo
c$$$               enddo
c$$$            endif
c$$$            if(useqtc .eq. 1) then 
c$$$               do iqtc = 1,nqtc
c$$$                  do ipop = 1,npop
c$$$                     meanqtc(ipop,iqtc) =  meanqtctmp(ipop,iqtc)
c$$$                     sdqtc(ipop,iqtc) =  sdqtctmp(ipop,iqtc)
c$$$                  enddo
c$$$               enddo
c$$$            endif
c$$$         else 
c$$$            ctmp(ipp) = c(ipp)
c$$$c$$$  do iloc=1,nloc
c$$$c$$$  do ial=1,nal(iloc)
c$$$c$$$  ftmp(ipop1,iloc,ial) = f(ipop1,iloc,ial)
c$$$c$$$  ftmp(ipop2,iloc,ial) = f(ipop2,iloc,ial)
c$$$c$$$  enddo
c$$$c$$$  enddo
c$$$         endif        
c$$$      enddo
c$$$      end subroutine udcf2parvq
c$$$***********************************************************************





***********************************************************************
*     joint update of c, f (under the CFM) 
*     and mean, variance of quant. variables
*     single component update of c
*     new f is proposed according to full conditional pi(f*|c*,yy)
      subroutine udcfallvarcfm3(npop,npopmax,f,fa,drift,
     &     nlocd,nloch,nql,ncolt,
     &     nal,nalmax,indcell,nindiv,npp,nppmax,c,ctmp,nqtc,qtc,
     &     meanqtc,sdqtc,meanqtctmp,sdqtctmp,nnqtc,sqtc,ssqtc,
     &     ksiqtc,kappaqtc,alphaqtc,betaqtc,
     &     a,ptmp,ftmp,yy,z,ql,n,ntmp,nudcel,
     &     usegeno2,usegeno1,useql,useqtc,ploidy)
      implicit none
      integer npop,npopmax,nlocd,nloch,nql,ncolt,nal(ncolt),
     &     nalmax,nindiv,npp,nppmax,indcell(nindiv),
     &     c(nppmax),ctmp(nppmax),yy(nindiv,2*nlocd+2*nloch),
     &     z(nindiv,nloch),ql(nindiv,nql),n(npopmax,ncolt,nalmax),
     &     ntmp(npopmax,ncolt,nalmax),nudcel,nnqtc,nqtc,
     &     usegeno2,usegeno1,useql,useqtc,ploidy
      double precision f,drift,ftmp,a,ptmp,fa,qtc,
     &     meanqtc,sdqtc,meanqtctmp,sdqtctmp,sqtc,ssqtc,
     &     ksiqtc,kappaqtc,alphaqtc,betaqtc
      dimension f(npopmax,ncolt,nalmax),drift(npopmax),
     &      ftmp(npopmax,ncolt,nalmax),
     &     a(nalmax),ptmp(nalmax),fa(ncolt,nalmax),
     &     qtc(nindiv,nqtc),meanqtc(npopmax,nqtc),
     &     sdqtc(npopmax,nqtc),
     &     meanqtctmp(npopmax,nqtc),sdqtctmp(npopmax,nqtc),
     &     nnqtc(npopmax,nqtc),sqtc(npopmax,nqtc),
     &     ssqtc(npopmax,nqtc),
     &     ksiqtc(nqtc),kappaqtc(nqtc),alphaqtc(nqtc),betaqtc(nqtc) 
      integer ipop,ipp,ipop1,ipop2,iloc,ial,iud,iqtc
      double precision alpha,ggrunif,lratio,lTfallvar,contriblr,bern,
     &     ggrbinom,lrallvar2

*     init. tmp. vector of population membership
      do ipp=1,npp
          ctmp(ipp) = c(ipp)
      enddo
      if(nppmax .gt. npp) then
         do ipp=npp+1,nppmax
            ctmp(ipp) = -999
         enddo
      endif

      do iud=1,nudcel
      ipp = 1 + idint(dint(dble(npp)*ggrunif(0.d0,1.d0)))
c         write(*,*) 'ipp=',ipp
*     propose new labeling of a tile
         ctmp(ipp) = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
         ipop1 = c(ipp)
         ipop2 = ctmp(ipp)
         bern = 0
         if(ipop1 .ne. ipop2) then 
            lratio = 0
            if(((usegeno2 .eq. 1) .or. (usegeno1 .eq. 1)) .or. 
     &           useql .eq. 1) then
*     counting alleles for both states of c
c$$$            call countn(nindiv,nlocmax,nlocmax2,npopmax,
c$$$     &           nppmax,nal,nalmax,yy,n,indcell,c,ploidy)
c$$$            call countn(nindiv,nlocmax,nlocmax2,npopmax,
c$$$     &           nppmax,nal,nalmax,yy,ntmp,indcell,ctmp,ploidy)
               call countnallvar2(nindiv,nlocd,nloch,nql,ncolt,
     &              nppmax,nal,nalmax,yy,z,ql,n,indcell,c,npopmax,
     &              usegeno2,usegeno1,useql,ploidy)
               call countnallvar2(nindiv,nlocd,nloch,nql,ncolt,
     &              nppmax,nal,nalmax,yy,z,ql,ntmp,indcell,ctmp,npopmax,
     &              usegeno2,usegeno1,useql,ploidy)
c$$$               lratio = lratio + 
c$$$     &              lTf(ipop1,ntmp,fa,drift,npopmax,nloc,nal,nalmax)
c$$$     &              + lTf(ipop2,ntmp,fa,drift,npopmax,nloc,nal,nalmax)
c$$$     &              - lTf(ipop1,n,fa,drift,npopmax,nloc,nal,nalmax)
c$$$     &              - lTf(ipop2,n,fa,drift,npopmax,nloc,nal,nalmax)
               lratio = lratio + 
     &              lTfallvar(ipop1,ntmp,fa,drift,npopmax,nlocd,nloch,
     &              nql,ncolt,nal,nalmax,usegeno2,usegeno1,useql) +
     &              lTfallvar(ipop2,ntmp,fa,drift,npopmax,nlocd,nloch,
     &              nql,ncolt,nal,nalmax,usegeno2,usegeno1,useql) -
     &              lTfallvar(ipop1,n,fa,drift,npopmax,nlocd,nloch,
     &              nql,ncolt,nal,nalmax,usegeno2,usegeno1,useql) - 
     &              lTfallvar(ipop2,n,fa,drift,npopmax,nlocd,nloch,
     &              nql,ncolt,nal,nalmax,usegeno2,usegeno1,useql)
            endif
            if(useqtc .eq. 1) then 
*     propose means and variances 
               call propparqvudc3(qtc,nindiv,nqtc,indcell,c,ctmp,
     &           meanqtc,sdqtc,meanqtctmp,sdqtctmp,nnqtc,
     &              sqtc,ssqtc,npop,npopmax,nppmax,
     &              ksiqtc,kappaqtc,alphaqtc,betaqtc,ipop1,ipop2,
     &              contriblr)
*     contrib likelihood
*     argument usegeno2,usegeno1,useql passed as 0 to disregard these variables 
*     hence avoid having their contrib to likelihood ratio  twice 
c$$$               lratio = lratio + 
c$$$     &              lratiogq(yy,f,c,ctmp,indcell,indcell,
c$$$     &              npopmax,nlocmax,nalmax,nindiv,nloc,nlocmax2,
c$$$     &              nppmax,ploidy,qtc,nqtc,meanqtc,sdqtc,
c$$$     &              meanqtctmp,sdqtctmp,0,useqtc)
               lratio = lratio + lrallvar2(yy,z,ql,
     &              f,c,ctmp,indcell,indcell,
     &              npopmax,nlocd,nloch,nql,ncolt,nalmax,nindiv,
     &              nppmax,qtc,nqtc,meanqtc,sdqtc,meanqtctmp,sdqtctmp,
     &              0,0,0,useqtc,ploidy)               
               lratio = lratio + contriblr 
            endif
            lratio = dmin1(0.d0,lratio)
            alpha = dexp(lratio)
            bern = ggrbinom(1.d0,alpha)
         endif
            
c     write(*,*) 'bern=',bern
            
         if(bern .eq. 1) then
            c(ipp) = ctmp(ipp)
            if((usegeno2 .eq. 1) .or. (usegeno1 .eq. 1)) then 
c               write(*,*) 'update c,f: bern=1 usegeno2=',usegeno2
c$$$               call samplef2(npop,npopmax,nloc,nlocmax,
c$$$     &              nal,nalmax,ipop1,ipop2,f,ftmp,
c$$$  &              fa,drift,a,ptmp,ntmp)
               call samplef2allvar(npop,npopmax,nlocd,nloch,nql,ncolt,
     &              nal,nalmax,ipop1,ipop2,f,ftmp,fa,drift,a,ptmp,ntmp,
     &              usegeno2,usegeno1,useql)
               do iloc=1,ncolt
                  do ial=1,nal(iloc)
                     f(ipop1,iloc,ial) = ftmp(ipop1,iloc,ial)
                     f(ipop2,iloc,ial) = ftmp(ipop2,iloc,ial)
                  enddo
               enddo
            endif
            if(useqtc .eq. 1) then 
               do iqtc = 1,nqtc
                  do ipop = 1,npop
                     meanqtc(ipop,iqtc) =  meanqtctmp(ipop,iqtc)
                     sdqtc(ipop,iqtc) =  sdqtctmp(ipop,iqtc)
                  enddo
               enddo
            endif
         else 
            ctmp(ipp) = c(ipp)
c$$$  do iloc=1,nloc
c$$$  do ial=1,nal(iloc)
c$$$  ftmp(ipop1,iloc,ial) = f(ipop1,iloc,ial)
c$$$  ftmp(ipop2,iloc,ial) = f(ipop2,iloc,ial)
c$$$  enddo
c$$$  enddo
         endif        
      enddo
      end subroutine udcfallvarcfm3
************************************************************************




c$$$***********************************************************************
c$$$*     joint update of c, f (under the CFM) 
c$$$*     and mean, variance of quant. variables
c$$$*     single component update of c
c$$$*     new f is proposed according to full conditional pi(f*|c*,yy)
c$$$      subroutine udcfallvarcfm2(npop,npopmax,f,fa,drift,
c$$$     &     nlocd,nloch,nql,ncolt,
c$$$     &     nal,nalmax,indcell,nindiv,npp,nppmax,c,ctmp,nqtc,qtc,
c$$$     &     meanqtc,sdqtc,meanqtctmp,sdqtctmp,nnqtc,sqtc,ssqtc,
c$$$     &     ksiqtc,kappaqtc,alphaqtc,betaqtc,
c$$$     &     a,ptmp,ftmp,yy,z,ql,n,ntmp,nudcel,
c$$$     &     usegeno2,usegeno1,useql,useqtc,ploidy)
c$$$      implicit none
c$$$      integer npop,npopmax,nlocd,nloch,nql,ncolt,nal(ncolt),
c$$$     &     nalmax,nindiv,npp,nppmax,indcell(nindiv),
c$$$     &     c(nppmax),ctmp(nppmax),yy(nindiv,2*nlocd+2*nloch),
c$$$     &     z(nindiv,nloch),ql(nindiv,nql),n(npopmax,ncolt,nalmax),
c$$$     &     ntmp(npopmax,ncolt,nalmax),nudcel,nnqtc,nqtc,
c$$$     &     usegeno2,usegeno1,useql,useqtc,ploidy
c$$$      double precision f,drift,ftmp,a,ptmp,fa,qtc,
c$$$     &     meanqtc,sdqtc,meanqtctmp,sdqtctmp,sqtc,ssqtc,
c$$$     &     ksiqtc,kappaqtc,alphaqtc,betaqtc
c$$$      dimension f(npopmax,ncolt,nalmax),drift(npopmax),
c$$$     &      ftmp(npopmax,ncolt,nalmax),
c$$$     &     a(nalmax),ptmp(nalmax),fa(ncolt,nalmax),
c$$$     &     qtc(nindiv,nqtc),meanqtc(npopmax,nqtc),
c$$$     &     sdqtc(npopmax,nqtc),
c$$$     &     meanqtctmp(npopmax,nqtc),sdqtctmp(npopmax,nqtc),
c$$$     &     nnqtc(npopmax,nqtc),sqtc(npopmax,nqtc),
c$$$     &     ssqtc(npopmax,nqtc),
c$$$     &     ksiqtc(nqtc),kappaqtc(nqtc),alphaqtc(nqtc),betaqtc(nqtc) 
c$$$      integer ipop,ipp,ipop1,ipop2,iloc,ial,iud,iqtc
c$$$      double precision alpha,ggrunif,lratio,lTfallvar,contriblr,bern,
c$$$     &     ggrbinom,lrallvar2
c$$$
c$$$*     init. tmp. vector of population membership
c$$$      do ipp=1,npp
c$$$          ctmp(ipp) = c(ipp)
c$$$      enddo
c$$$      if(nppmax .gt. npp) then
c$$$         do ipp=npp+1,nppmax
c$$$            ctmp(ipp) = -999
c$$$         enddo
c$$$      endif
c$$$
c$$$      do iud=1,nudcel
c$$$      ipp = 1 + idint(dint(dble(npp)*ggrunif(0.d0,1.d0)))
c$$$c         write(*,*) 'ipp=',ipp
c$$$*     propose new labeling of a tile
c$$$         ctmp(ipp) = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
c$$$         ipop1 = c(ipp)
c$$$         ipop2 = ctmp(ipp)
c$$$         bern = 0
c$$$         if(ipop1 .ne. ipop2) then 
c$$$            lratio = 0
c$$$            if(((usegeno2 .eq. 1) .or. (usegeno1 .eq. 1)) .or. 
c$$$     &           useql .eq. 1) then
c$$$*     counting alleles for both states of c
c$$$c$$$            call countn(nindiv,nlocmax,nlocmax2,npopmax,
c$$$c$$$     &           nppmax,nal,nalmax,yy,n,indcell,c,ploidy)
c$$$c$$$            call countn(nindiv,nlocmax,nlocmax2,npopmax,
c$$$c$$$     &           nppmax,nal,nalmax,yy,ntmp,indcell,ctmp,ploidy)
c$$$               call countnallvar2(nindiv,nlocd,nloch,nql,ncolt,
c$$$     &              nppmax,nal,nalmax,yy,z,ql,n,indcell,c,npopmax,
c$$$     &              usegeno2,usegeno1,useql,ploidy)
c$$$               call countnallvar2(nindiv,nlocd,nloch,nql,ncolt,
c$$$     &              nppmax,nal,nalmax,yy,z,ql,ntmp,indcell,ctmp,npopmax,
c$$$     &              usegeno2,usegeno1,useql,ploidy)
c$$$c$$$               lratio = lratio + 
c$$$c$$$     &              lTf(ipop1,ntmp,fa,drift,npopmax,nloc,nal,nalmax)
c$$$c$$$     &              + lTf(ipop2,ntmp,fa,drift,npopmax,nloc,nal,nalmax)
c$$$c$$$     &              - lTf(ipop1,n,fa,drift,npopmax,nloc,nal,nalmax)
c$$$c$$$     &              - lTf(ipop2,n,fa,drift,npopmax,nloc,nal,nalmax)
c$$$               lratio = lratio + 
c$$$     &              lTfallvar(ipop1,ntmp,fa,drift,npopmax,nlocd,nloch,
c$$$     &              nql,ncolt,nal,nalmax,usegeno2,usegeno1,useql) +
c$$$     &              lTfallvar(ipop2,ntmp,fa,drift,npopmax,nlocd,nloch,
c$$$     &              nql,ncolt,nal,nalmax,usegeno2,usegeno1,useql) -
c$$$     &              lTfallvar(ipop1,n,fa,drift,npopmax,nlocd,nloch,
c$$$     &              nql,ncolt,nal,nalmax,usegeno2,usegeno1,useql) - 
c$$$     &              lTfallvar(ipop2,n,fa,drift,npopmax,nlocd,nloch,
c$$$     &              nql,ncolt,nal,nalmax,usegeno2,usegeno1,useql)
c$$$            endif
c$$$            if(useqtc .eq. 1) then 
c$$$*     propose means and variances 
c$$$               call propparqvudc(qtc,nindiv,nqtc,indcell,c,ctmp,
c$$$     &           meanqtc,sdqtc,meanqtctmp,sdqtctmp,nnqtc,
c$$$     &              sqtc,ssqtc,npop,npopmax,nppmax,
c$$$     &              ksiqtc,kappaqtc,alphaqtc,betaqtc,ipop1,ipop2,
c$$$     &              contriblr)
c$$$*     contrib likelihood
c$$$*     argument usegeno2,usegeno1,useql passed as 0 to disregard these variables 
c$$$*     hence avoid having their contrib to likelihood ratio  twice 
c$$$c$$$               lratio = lratio + 
c$$$c$$$     &              lratiogq(yy,f,c,ctmp,indcell,indcell,
c$$$c$$$     &              npopmax,nlocmax,nalmax,nindiv,nloc,nlocmax2,
c$$$c$$$     &              nppmax,ploidy,qtc,nqtc,meanqtc,sdqtc,
c$$$c$$$     &              meanqtctmp,sdqtctmp,0,useqtc)
c$$$               lratio = lratio + lrallvar2(yy,z,ql,
c$$$     &              f,c,ctmp,indcell,indcell,
c$$$     &              npopmax,nlocd,nloch,nql,ncolt,nalmax,nindiv,
c$$$     &              nppmax,qtc,nqtc,meanqtc,sdqtc,meanqtctmp,sdqtctmp,
c$$$     &              0,0,0,useqtc,ploidy)               
c$$$               lratio = lratio + contriblr 
c$$$            endif
c$$$            lratio = dmin1(0.d0,lratio)
c$$$            alpha = dexp(lratio)
c$$$            bern = ggrbinom(1.d0,alpha)
c$$$         endif
c$$$            
c$$$c     write(*,*) 'bern=',bern
c$$$            
c$$$         if(bern .eq. 1) then
c$$$            c(ipp) = ctmp(ipp)
c$$$            if((usegeno2 .eq. 1) .or. (usegeno1 .eq. 1)) then 
c$$$c               write(*,*) 'update c,f: bern=1 usegeno2=',usegeno2
c$$$c$$$               call samplef2(npop,npopmax,nloc,nlocmax,
c$$$c$$$     &              nal,nalmax,ipop1,ipop2,f,ftmp,
c$$$c$$$  &              fa,drift,a,ptmp,ntmp)
c$$$               call samplef2allvar(npop,npopmax,nlocd,nloch,nql,ncolt,
c$$$     &              nal,nalmax,ipop1,ipop2,f,ftmp,fa,drift,a,ptmp,ntmp,
c$$$     &              usegeno2,usegeno1,useql)
c$$$               do iloc=1,ncolt
c$$$                  do ial=1,nal(iloc)
c$$$                     f(ipop1,iloc,ial) = ftmp(ipop1,iloc,ial)
c$$$                     f(ipop2,iloc,ial) = ftmp(ipop2,iloc,ial)
c$$$                  enddo
c$$$               enddo
c$$$            endif
c$$$            if(useqtc .eq. 1) then 
c$$$               do iqtc = 1,nqtc
c$$$                  do ipop = 1,npop
c$$$                     meanqtc(ipop,iqtc) =  meanqtctmp(ipop,iqtc)
c$$$                     sdqtc(ipop,iqtc) =  sdqtctmp(ipop,iqtc)
c$$$                  enddo
c$$$               enddo
c$$$            endif
c$$$         else 
c$$$            ctmp(ipp) = c(ipp)
c$$$c$$$  do iloc=1,nloc
c$$$c$$$  do ial=1,nal(iloc)
c$$$c$$$  ftmp(ipop1,iloc,ial) = f(ipop1,iloc,ial)
c$$$c$$$  ftmp(ipop2,iloc,ial) = f(ipop2,iloc,ial)
c$$$c$$$  enddo
c$$$c$$$  enddo
c$$$         endif        
c$$$      enddo
c$$$      end subroutine udcfallvarcfm2
c$$$************************************************************************
c$$$
c$$$

************************************************************************
*     Modification de u
*     composante par composante 
*     avec proposal uniforme sur un carre de cote du 
*     centre sur le point courant (random walk)
      subroutine updu(npp,nppmax,c,u,zz,nindiv,nloc,nlocmax,
     &     nlocmax2,nalmax,npopmax,f,indcell,distcell,
     &     indcelltmp,distcelltmp,
     &     s,xlim,ylim,du,ploidy,nudcel)
      implicit none 
      integer npp, nppmax,c(nppmax),nindiv,nloc,nlocmax,
     &     nlocmax2,nalmax,npopmax,zz(nindiv,nlocmax2),
     &     indcell(nindiv),ploidy,nudcel
      double precision u(2,nppmax),f(npopmax,nlocmax,nalmax),
     &     distcell(nindiv),s(2,nindiv),xlim(2),ylim(2),du
      integer ipp,iindiv,indcelltmp(nindiv),iud
      double precision utmp(2,nppmax),ggrunif,r,alpha,
     &     distcelltmp(nindiv),surf,surftmp,dx,dy,ratio
      double precision bern,ggrbinom
*     initialisation du tableau tmporaire
      do ipp=1,npp
         utmp(1,ipp) = u(1,ipp)
         utmp(2,ipp) = u(2,ipp)
      enddo
      if(nppmax .gt. npp) then
         do ipp=npp+1,nppmax
            utmp(1,ipp) = -999.
            utmp(2,ipp) = -999.
         enddo
      endif

c      write(*,*) 'npp=', npp
c      write(*,*) 'u=', u

      do iud=1,nudcel
      ipp = 1 + idint(dint(dble(npp)*ggrunif(0.d0,1.d0)))
*     proposition d un deplacement d un point de u
         utmp(1,ipp) = max(u(1,ipp)-du/2.,xlim(1)) + ggrunif(0.d0,1.d0)*
     &        (min(u(1,ipp)+du/2.,xlim(2))-max(u(1,ipp)-du/2.,xlim(1)))
         utmp(2,ipp) = max(u(2,ipp)-du/2.,ylim(1)) + ggrunif(0.d0,1.d0)*
     &        (min(u(2,ipp)+du/2.,ylim(2))-max(u(2,ipp)-du/2.,ylim(1)))

*     calcul de l aire du domaine ou il pouvait aller 
         dx = min(du/2.,u(1,ipp)-xlim(1),xlim(2)-u(1,ipp))
         dy = min(du/2.,u(2,ipp)-ylim(1),ylim(2)-u(2,ipp))
         surf = (dx+du/2.)*(dy+du/2.)
         dx = min(du/2.,utmp(1,ipp)-xlim(1),xlim(2)-utmp(1,ipp))
         dy = min(du/2.,utmp(2,ipp)-ylim(1),ylim(2)-utmp(2,ipp))
         surftmp = (dx+du/2.)*(dy+du/2.)

*     modif de indcell et distcell
         call vormove(nindiv,s,npp,nppmax,utmp,
     &        indcell,distcell,indcelltmp,distcelltmp,ipp)
c         write(*,*) 'apres vormove'
         r = ratio(zz,f,c,c,indcell,indcelltmp,
     &     npopmax,nlocmax,nalmax,nindiv,nloc,nlocmax2,
     &     nppmax,ploidy)
c         write(*,*) 'r=',r
         r = r*surf/surftmp

         alpha = dmin1(1.d0,r)
c         write(*,*) 'alpha=',alpha
         bern = ggrbinom(1.d0,alpha)
         if(bern .eq. 1) then
            u(1,ipp) = utmp(1,ipp)
            u(2,ipp) = utmp(2,ipp)
            do iindiv=1,nindiv
               indcell(iindiv) = indcelltmp(iindiv)
               distcell(iindiv) = distcelltmp(iindiv)
            enddo
         else 
            utmp(1,ipp) = u(1,ipp)
            utmp(2,ipp) = u(2,ipp)
         endif
      enddo
      end subroutine updu
************************************************************************

************************************************************************
*     Update u component-wise
*     data consist of genotypes and/or quantitative variables
      subroutine updugq(npp,nppmax,c,u,zz,nindiv,nloc,nlocmax,
     &     nlocmax2,nalmax,npopmax,f,indcell,distcell,
     &     indcelltmp,distcelltmp,
     &     s,xlim,ylim,du,ploidy,nudcel,
     &     qtc,nqtc,meanqtc,sdqtc,usegeno2,useqtc)
      implicit none 
      integer npp, nppmax,c(nppmax),nindiv,nloc,nlocmax,
     &     nlocmax2,nalmax,npopmax,zz(nindiv,nlocmax2),
     &     indcell(nindiv),ploidy,nudcel,nqtc,usegeno2,useqtc
      double precision u(2,nppmax),f(npopmax,nlocmax,nalmax),
     &     distcell(nindiv),s(2,nindiv),xlim(2),ylim(2),du,
     &     qtc(nindiv,nqtc),
     &     meanqtc(npopmax,nqtc),sdqtc(npopmax,nqtc)
      integer ipp,iindiv,indcelltmp(nindiv),iud
      double precision utmp(2,nppmax),ggrunif,r,alpha,
     &     distcelltmp(nindiv),surf,surftmp,dx,dy,lratiogq
      double precision bern,ggrbinom
*     initialisation du tableau tmporaire
      do ipp=1,npp
         utmp(1,ipp) = u(1,ipp)
         utmp(2,ipp) = u(2,ipp)
      enddo
      if(nppmax .gt. npp) then
         do ipp=npp+1,nppmax
            utmp(1,ipp) = -999.
            utmp(2,ipp) = -999.
         enddo
      endif
      do iud=1,nudcel
      ipp = 1 + idint(dint(dble(npp)*ggrunif(0.d0,1.d0)))
*     proposition d un deplacement d un point de u
         utmp(1,ipp) = max(u(1,ipp)-du/2.,xlim(1)) + ggrunif(0.d0,1.d0)*
     &        (min(u(1,ipp)+du/2.,xlim(2))-max(u(1,ipp)-du/2.,xlim(1)))
         utmp(2,ipp) = max(u(2,ipp)-du/2.,ylim(1)) + ggrunif(0.d0,1.d0)*
     &        (min(u(2,ipp)+du/2.,ylim(2))-max(u(2,ipp)-du/2.,ylim(1)))

*     calcul de l aire du domaine ou il pouvait aller 
         dx = min(du/2.,u(1,ipp)-xlim(1),xlim(2)-u(1,ipp))
         dy = min(du/2.,u(2,ipp)-ylim(1),ylim(2)-u(2,ipp))
         surf = (dx+du/2.)*(dy+du/2.)
         dx = min(du/2.,utmp(1,ipp)-xlim(1),xlim(2)-utmp(1,ipp))
         dy = min(du/2.,utmp(2,ipp)-ylim(1),ylim(2)-utmp(2,ipp))
         surftmp = (dx+du/2.)*(dy+du/2.)

*     modif de indcell et distcell
         call vormove(nindiv,s,npp,nppmax,utmp,
     &        indcell,distcell,indcelltmp,distcelltmp,ipp)
c         write(*,*) 'apres vormove'
         r = dexp(lratiogq(zz,f,c,c,indcell,indcelltmp,
     &     npopmax,nlocmax,nalmax,nindiv,nloc,nlocmax2,
     &     nppmax,ploidy,qtc,nqtc,meanqtc,sdqtc,
     &        meanqtc,sdqtc,usegeno2,useqtc))
c         write(*,*) 'r=',r
         r = r*surf/surftmp

         alpha = dmin1(1.d0,r)
c         write(*,*) 'alpha=',alpha
         bern = ggrbinom(1.d0,alpha)
         if(bern .eq. 1) then
            u(1,ipp) = utmp(1,ipp)
            u(2,ipp) = utmp(2,ipp)
            do iindiv=1,nindiv
               indcell(iindiv) = indcelltmp(iindiv)
               distcell(iindiv) = distcelltmp(iindiv)
            enddo
         else 
            utmp(1,ipp) = u(1,ipp)
            utmp(2,ipp) = u(2,ipp)
         endif
      enddo
      end subroutine updugq
************************************************************************

************************************************************************
*     Update u component-wise
*     data consist of diploid and/or haploid genotypes 
*     and/or quantitative variables
      subroutine uduallvar2(npp,nppmax,c,u,yy,z,ql,nindiv,nlocd,
     &     nloch,nql,ncolt,nalmax,npopmax,f,indcell,distcell,
     &     indcelltmp,distcelltmp,s,xlim,ylim,du,nudcel,
     &     qtc,nqtc,meanqtc,sdqtc,usegeno2,usegeno1,useql,useqtc,
     &     ploidy)
      implicit none 
      integer npp, nppmax,c(nppmax),nindiv,nlocd,nloch,nql,ncolt,
     &     nalmax,npopmax,yy(nindiv,2*nlocd+2*nloch),z(nindiv,nloch),
     &     ql(nindiv,nql),indcell(nindiv),nudcel,nqtc,
     &     usegeno2,usegeno1,useql,useqtc,ploidy
      double precision u(2,nppmax),f(npopmax,ncolt,nalmax),
     &     distcell(nindiv),s(2,nindiv),xlim(2),ylim(2),du,
     &     qtc(nindiv,nqtc),
     &     meanqtc(npopmax,nqtc),sdqtc(npopmax,nqtc)
      integer ipp,iindiv,indcelltmp(nindiv),iud
      double precision utmp(2,nppmax),ggrunif,r,alpha,
     &     distcelltmp(nindiv),surf,surftmp,dx,dy,lrallvar2
      double precision bern,ggrbinom
*     initialisation du tableau tmporaire
      do ipp=1,npp
         utmp(1,ipp) = u(1,ipp)
         utmp(2,ipp) = u(2,ipp)
      enddo
      if(nppmax .gt. npp) then
         do ipp=npp+1,nppmax
            utmp(1,ipp) = -999.
            utmp(2,ipp) = -999.
         enddo
      endif
      do iud=1,nudcel
      ipp = 1 + idint(dint(dble(npp)*ggrunif(0.d0,1.d0)))
*     proposition d un deplacement d un point de u
         utmp(1,ipp) = max(u(1,ipp)-du/2.,xlim(1)) + ggrunif(0.d0,1.d0)*
     &        (min(u(1,ipp)+du/2.,xlim(2))-max(u(1,ipp)-du/2.,xlim(1)))
         utmp(2,ipp) = max(u(2,ipp)-du/2.,ylim(1)) + ggrunif(0.d0,1.d0)*
     &        (min(u(2,ipp)+du/2.,ylim(2))-max(u(2,ipp)-du/2.,ylim(1)))

*     calcul de l aire du domaine ou il pouvait aller 
         dx = min(du/2.,u(1,ipp)-xlim(1),xlim(2)-u(1,ipp))
         dy = min(du/2.,u(2,ipp)-ylim(1),ylim(2)-u(2,ipp))
         surf = (dx+du/2.)*(dy+du/2.)
         dx = min(du/2.,utmp(1,ipp)-xlim(1),xlim(2)-utmp(1,ipp))
         dy = min(du/2.,utmp(2,ipp)-ylim(1),ylim(2)-utmp(2,ipp))
         surftmp = (dx+du/2.)*(dy+du/2.)

*     modif de indcell et distcell
         call vormove(nindiv,s,npp,nppmax,utmp,
     &        indcell,distcell,indcelltmp,distcelltmp,ipp)
c         write(*,*) 'apres vormove'
c$$$         r = dexp(lratiogq(yy,f,c,c,indcell,indcelltmp,
c$$$     &     npopmax,nlocmax,nalmax,nindiv,nloc,nlocmax2,
c$$$     &     nppmax,ploidy,qtc,nqtc,meanqtc,sdqtc,
c$$$     &        meanqtc,sdqtc,usegeno2,useqtc))
         r = dexp(lrallvar2(yy,z,ql,f,c,c,indcell,indcelltmp,
     &     npopmax,nlocd,nloch,nql,ncolt,nalmax,nindiv,nppmax,
     &     qtc,nqtc,meanqtc,sdqtc,meanqtc,sdqtc,
     &     usegeno2,usegeno1,useql,useqtc,ploidy))
c         write(*,*) 'r=',r
         r = r*surf/surftmp

         alpha = dmin1(1.d0,r)
c         write(*,*) 'alpha=',alpha
         bern = ggrbinom(1.d0,alpha)
         if(bern .eq. 1) then
            u(1,ipp) = utmp(1,ipp)
            u(2,ipp) = utmp(2,ipp)
            do iindiv=1,nindiv
               indcell(iindiv) = indcelltmp(iindiv)
               distcell(iindiv) = distcelltmp(iindiv)
            enddo
         else 
            utmp(1,ipp) = u(1,ipp)
            utmp(2,ipp) = u(2,ipp)
         endif
      enddo
      end subroutine uduallvar2
************************************************************************




************************************************************************
*     mise a jour de t    
* 
      subroutine updt(npp,nppmax,nindiv,
     &     nloc,nlocmax,nlocmax2,nalmax,npopmax,
     &     t,ttmp,dt,s,c,indcell,distcell,indcelltmp,distcelltmp,
     &     u,zz,f,ploidy)
      implicit none 
      integer npp,nppmax,nindiv,nloc,nlocmax,nlocmax2,nalmax,
     &     npopmax,c(nppmax),indcell(nindiv),zz(nindiv,nlocmax2),
     &     ploidy
      double precision t(2,nindiv),s(2,nindiv),distcell(nindiv),
     &     u(2,nppmax),f(npopmax,nlocmax,nalmax),dt
      integer iindiv,ipp,indcelltmp(nindiv)
      double precision ggrunif,d,ttmp(2,nindiv),r,alpha,
     &     distcelltmp(nindiv),ratio
      double precision bern,ggrbinom,accept

*     initialisation
      do iindiv = 1,nindiv
         ttmp(1,iindiv) = t(1,iindiv)
         ttmp(2,iindiv) = t(2,iindiv)
         indcelltmp(iindiv) = indcell(iindiv)
         distcelltmp(iindiv) = distcell(iindiv)
      enddo

c      do iindiv = 1,nindiv
      iindiv= 1 + idint(dint(dble(nindiv)*ggrunif(0.d0,1.d0)))
*     proposition d'une modif de t
         ttmp(1,iindiv) = s(1,iindiv) + dt*(ggrunif(0.d0,1.d0)-.5)
         ttmp(2,iindiv) = s(2,iindiv) + dt*(ggrunif(0.d0,1.d0)-.5)

*     modif de indcell et distcell
         distcelltmp(iindiv) = 3.d+300
         do ipp = 1,npp
            d = (ttmp(1,iindiv)-u(1,ipp))**2+
     &           (ttmp(2,iindiv)-u(2,ipp))**2
            if(d .lt. distcelltmp(iindiv)) then 
               indcelltmp(iindiv)  = ipp
               distcelltmp(iindiv) = d
            endif
         enddo

*     proba d'acceptation
         if(indcelltmp(iindiv) .ne. indcell(iindiv)) then 
            r = ratio(zz,f,c,c,indcell,indcelltmp,
     &           npopmax,nlocmax,nalmax,nindiv,nloc,
     &           nlocmax2,nppmax,ploidy)
         else 
            r = 1.
         endif
         alpha = dmin1(1.d0,r)
         accept = ggrbinom(1.d0,alpha)
*     mise a jour en cas d'acceptation
         if(accept .eq. 1) then 
            indcell(iindiv) = indcelltmp(iindiv)
            distcell(iindiv) = distcelltmp(iindiv)
            t(1,iindiv) = ttmp(1,iindiv) 
            t(2,iindiv) = ttmp(2,iindiv)
         endif
c      enddo
      end subroutine updt
************************************************************************





************************************************************************
*     update true un-observed spatial coordinates of individuals  
      subroutine udtallvar2(npp,nppmax,nindiv,
     &     nlocd,nloch,nql,ncolt,nqtc,nalmax,npopmax,
     &     t,ttmp,dt,s,c,indcell,distcell,indcelltmp,distcelltmp,
     &     u,yy,z,ql,qtc,f,meanqtc,sdqtc,
     &     usegeno2,usegeno1,useql,useqtc,ploidy)
      implicit none 
      integer npp,nppmax,nindiv,nlocd,nlocd2,nloch,nql,ncolt,nqtc,
     &     nalmax,npopmax,c(nppmax),indcell(nindiv),
     &     yy(nindiv,2*nlocd+2*nloch),z(nindiv,nloch),ql(nindiv,nql)
      double precision t(2,nindiv),s(2,nindiv),distcell(nindiv),
     &     u(2,nppmax),f(npopmax,ncolt,nalmax),dt,qtc(nindiv,nqtc),
     &     meanqtc(npopmax,nqtc),sdqtc(npopmax,nqtc)
      integer iindiv,ipp,indcelltmp(nindiv),
     &     usegeno2,usegeno1,useql,useqtc,ploidy
      double precision ggrunif,d,ttmp(2,nindiv),r,alpha,
     &     distcelltmp(nindiv),lrallvar2
      double precision bern,ggrbinom,accept

*     initialisation
      do iindiv = 1,nindiv
         ttmp(1,iindiv) = t(1,iindiv)
         ttmp(2,iindiv) = t(2,iindiv)
         indcelltmp(iindiv) = indcell(iindiv)
         distcelltmp(iindiv) = distcell(iindiv)
      enddo

c      do iindiv = 1,nindiv
      iindiv= 1 + idint(dint(dble(nindiv)*ggrunif(0.d0,1.d0)))
*     proposition d'une modif de t
         ttmp(1,iindiv) = s(1,iindiv) + dt*(ggrunif(0.d0,1.d0)-.5)
         ttmp(2,iindiv) = s(2,iindiv) + dt*(ggrunif(0.d0,1.d0)-.5)

*     modif de indcell et distcell
         distcelltmp(iindiv) = 3.d+300
         do ipp = 1,npp
            d = (ttmp(1,iindiv)-u(1,ipp))**2+
     &           (ttmp(2,iindiv)-u(2,ipp))**2
            if(d .lt. distcelltmp(iindiv)) then 
               indcelltmp(iindiv)  = ipp
               distcelltmp(iindiv) = d
            endif
         enddo

*     proba d'acceptation
         if(indcelltmp(iindiv) .ne. indcell(iindiv)) then 
c$$$            r = ratio(yy,f,c,c,indcell,indcelltmp,
c$$$     &           npopmax,nlocmax,nalmax,nindiv,nloc,
c$$$     &           nlocmax2,nppmax,ploidy)
            r = dexp(lrallvar2(yy,z,ql,f,c,c,indcell,indcelltmp,
     &     npopmax,nlocd,nloch,nql,ncolt,nalmax,nindiv,nppmax,
     &     qtc,nqtc,meanqtc,sdqtc,meanqtc,sdqtc,
     &     usegeno2,usegeno1,useql,useqtc,ploidy))
         else 
            r = 1.
         endif
         alpha = dmin1(1.d0,r)
         accept = ggrbinom(1.d0,alpha)
*     mise a jour en cas d'acceptation
         if(accept .eq. 1) then 
            indcell(iindiv) = indcelltmp(iindiv)
            distcell(iindiv) = distcelltmp(iindiv)
            t(1,iindiv) = ttmp(1,iindiv) 
            t(2,iindiv) = ttmp(2,iindiv)
         endif
c      enddo
      end subroutine udtallvar2
************************************************************************




************************************************************************
*     naissance ou mort d'une cellule
*     avec prior Poisson(lambda) tronquée :   0 < m < nppmax
      subroutine bdpp(nindiv,u,c,utmp,ctmp,npop,npopmax,
     &     nloc,nlocmax,nlocmax2,nalmax,npp,nppmax,zz,f,s,xlim,ylim,
     &     indcell,distcell,indcelltmp,distcelltmp,lambda,ploidy)
      implicit none 
      integer nindiv,nloc,nlocmax,nlocmax2,
     &     npop,npopmax,
     &     nalmax,npp,nppmax,zz(nindiv,nlocmax2),c(nppmax),
     &     indcell(nindiv),ploidy
      double precision u(2,nindiv),f(npopmax,nlocmax,nalmax),xlim(2),
     &     ylim(2),s(2,nindiv),distcell(nindiv),lambda

      integer ctmp(nppmax),indcelltmp(nindiv),ipp,npptmp,
     &     iindiv,ipprem
      double precision utmp(2,nppmax),distcelltmp(nindiv),ggrunif,
     &     ratio,r,alpha,ggrbinom,b,accept
      
*     naissance ou mort ?
      b = ggrbinom(1.d0,0.5d0)

      if(b .eq. 1) then
         if(npp .ne. nppmax) then 
*     naissance
            do ipp = 1,npp
               utmp (1,ipp) = u(1,ipp)
               utmp (2,ipp) = u(2,ipp)
               ctmp(ipp) = c(ipp)
            enddo
            npptmp = npp + 1
            ctmp(npptmp) = 1+ idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
            utmp(1,npptmp) = xlim(1)+(xlim(2)-xlim(1))*
     &           ggrunif(0.d0,1.d0)
            utmp(2,npptmp) = ylim(1)+(ylim(2)-ylim(1))*
     &           ggrunif(0.d0,1.d0)
            if(nppmax .gt. npptmp) then
               do ipp=npptmp+1,nppmax
                  ctmp(ipp) = -999
                  utmp(1,ipp) = -999.
                  utmp(2,ipp) = -999.
               enddo
            endif
            
            call voradd(s,utmp,indcell,distcell,indcelltmp,
     &           distcelltmp,nindiv,npp,nppmax)
            r = ratio(zz,f,c,ctmp,indcell,indcelltmp,npopmax,nlocmax,
     &           nalmax,nindiv,nloc,nlocmax2,nppmax,ploidy)
            r = r*lambda/dble(npp+1)
            alpha = dmin1(1.d0,r)
            accept = ggrbinom(1.d0,alpha)
            if(accept .eq. 1) then 
               npp = npptmp
               do iindiv=1,nindiv
                  indcell(iindiv) = indcelltmp(iindiv)
                  distcell(iindiv) = distcelltmp(iindiv)
               enddo
               do ipp = 1,nppmax
                  u (1,ipp) = utmp(1,ipp)
                  u (2,ipp) = utmp(2,ipp)
                  c(ipp) = ctmp(ipp)
               enddo
            endif
         endif
      else
*     mort
         if(npp .ne. 1) then 
            ipprem = 1+ aint(dble(npp)*ggrunif(0.d0,1.d0))
            if(ipprem .ne. 1) then 
               do ipp = 1,ipprem-1
                  utmp (1,ipp) = u(1,ipp)
                  utmp (2,ipp) = u(2,ipp)
                  ctmp(ipp) = c(ipp)
               enddo
            endif
            if(ipprem .ne. npp) then 
               do ipp = ipprem,npp-1
                  utmp (1,ipp) = u(1,ipp+1)
                  utmp (2,ipp) = u(2,ipp+1)
                  ctmp(ipp) = c(ipp+1)
               enddo
            endif
            do ipp=npp,nppmax
               utmp (1,ipp) = -999.
               utmp (2,ipp) = -999.
               ctmp(ipp) = -999
            enddo

            call vorrem(s,utmp,ipprem,indcell,distcell,
     &           indcelltmp,distcelltmp,nindiv,npp,nppmax)

            r = ratio(zz,f,c,ctmp,indcell,indcelltmp,npopmax,nlocmax,
     &           nalmax,nindiv,nloc,nlocmax2,nppmax,ploidy)
            r = r*dble(npp)/lambda
            alpha = dmin1(1.d0,r)
            accept = ggrbinom(1.d0,alpha)
            if(accept .eq. 1) then 
               npp = npp-1
               do iindiv=1,nindiv
                  indcell(iindiv) = indcelltmp(iindiv)
                  distcell(iindiv) = distcelltmp(iindiv)
               enddo
               do ipp = 1,nppmax
                  u (1,ipp) = utmp(1,ipp)
                  u (2,ipp) = utmp(2,ipp)
                  c(ipp) = ctmp(ipp)
               enddo
            endif
         endif
      endif
      end subroutine bdpp
************************************************************************



************************************************************************
*     birth and death of tiles
*     double sided truncated Poisson(lambda)  prior:   0 < m < nppmax
      subroutine bdppgq(nindiv,u,c,utmp,ctmp,npop,npopmax,
     &     nloc,nlocmax,nlocmax2,nalmax,npp,nppmax,zz,f,s,xlim,ylim,   
     &     indcell,distcell,indcelltmp,distcelltmp,lambda,ploidy,
     &     qtc,nqtc,meanqtc,sdqtc,usegeno2,useqtc)
      implicit none 
      integer nindiv,nloc,nlocmax,nlocmax2,
     &     npop,npopmax,
     &     nalmax,npp,nppmax,zz(nindiv,nlocmax2),c(nppmax),
     &     indcell(nindiv),ploidy,nqtc,usegeno2,useqtc
      double precision u(2,nindiv),f(npopmax,nlocmax,nalmax),xlim(2),
     &     ylim(2),s(2,nindiv),distcell(nindiv),lambda,
     &     qtc(nindiv,nqtc),
     &     meanqtc(npopmax,nqtc),sdqtc(npopmax,nqtc)

      integer ctmp(nppmax),indcelltmp(nindiv),ipp,npptmp,
     &     iindiv,ipprem
      double precision utmp(2,nppmax),distcelltmp(nindiv),ggrunif,
     &     lratiogq,r,alpha,ggrbinom,accept,b
      
*     naissance ou mort ?
      b = ggrbinom(1.d0,0.5d0)

      if(b .eq. 1) then
         if(npp .ne. nppmax) then 
*     naissance
            do ipp = 1,npp
               utmp (1,ipp) = u(1,ipp)
               utmp (2,ipp) = u(2,ipp)
               ctmp(ipp) = c(ipp)
            enddo
            npptmp = npp + 1
            ctmp(npptmp) = 1+ idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
            utmp(1,npptmp) = xlim(1)+(xlim(2)-xlim(1))*
     &           ggrunif(0.d0,1.d0)
            utmp(2,npptmp) = ylim(1)+(ylim(2)-ylim(1))*
     &           ggrunif(0.d0,1.d0)
            if(nppmax .gt. npptmp) then
               do ipp=npptmp+1,nppmax
                  ctmp(ipp) = -999
                  utmp(1,ipp) = -999.
                  utmp(2,ipp) = -999.
               enddo
            endif
            
            call voradd(s,utmp,indcell,distcell,indcelltmp,
     &           distcelltmp,nindiv,npp,nppmax)
            r = dexp(lratiogq(zz,f,c,ctmp,indcell,indcelltmp,npopmax,
     &           nlocmax,nalmax,nindiv,nloc,nlocmax2,nppmax,ploidy,
     &           qtc,nqtc,meanqtc,sdqtc,meanqtc,sdqtc,
     &           usegeno2,useqtc))
            r = r*lambda/dble(npp+1)
            alpha = dmin1(1.d0,r)
            accept = ggrbinom(1.d0,alpha)
            if(accept .eq. 1) then 
               npp = npptmp
               do iindiv=1,nindiv
                  indcell(iindiv) = indcelltmp(iindiv)
                  distcell(iindiv) = distcelltmp(iindiv)
               enddo
               do ipp = 1,nppmax
                  u (1,ipp) = utmp(1,ipp)
                  u (2,ipp) = utmp(2,ipp)
                  c(ipp) = ctmp(ipp)
               enddo
            endif
         endif
      else
*     mort
         if(npp .ne. 1) then 
            ipprem = 1+ aint(dble(npp)*ggrunif(0.d0,1.d0))
            if(ipprem .ne. 1) then 
               do ipp = 1,ipprem-1
                  utmp (1,ipp) = u(1,ipp)
                  utmp (2,ipp) = u(2,ipp)
                  ctmp(ipp) = c(ipp)
               enddo
            endif
            if(ipprem .ne. npp) then 
               do ipp = ipprem,npp-1
                  utmp (1,ipp) = u(1,ipp+1)
                  utmp (2,ipp) = u(2,ipp+1)
                  ctmp(ipp) = c(ipp+1)
               enddo
            endif
            do ipp=npp,nppmax
               utmp (1,ipp) = -999.
               utmp (2,ipp) = -999.
               ctmp(ipp) = -999
            enddo

            call vorrem(s,utmp,ipprem,indcell,distcell,
     &           indcelltmp,distcelltmp,nindiv,npp,nppmax)

            r = dexp(lratiogq(zz,f,c,ctmp,indcell,indcelltmp,npopmax,
     &           nlocmax,nalmax,nindiv,nloc,nlocmax2,nppmax,ploidy,
     &           qtc,nqtc,meanqtc,sdqtc,meanqtc,sdqtc,
     &           usegeno2,useqtc))
            r = r*dble(npp)/lambda
            alpha = dmin1(1.d0,r)
            accept = ggrbinom(1.d0,alpha)
            if(accept .eq. 1) then 
               npp = npp-1
               do iindiv=1,nindiv
                  indcell(iindiv) = indcelltmp(iindiv)
                  distcell(iindiv) = distcelltmp(iindiv)
               enddo
               do ipp = 1,nppmax
                  u (1,ipp) = utmp(1,ipp)
                  u (2,ipp) = utmp(2,ipp)
                  c(ipp) = ctmp(ipp)
               enddo
            endif
         endif
      endif
      end subroutine bdppgq
************************************************************************




************************************************************************
*     birth and death of cells
*     double sided truncated Poisson(lambda)  prior:   0 < m < nppmax
      subroutine bdcellallvar2(nindiv,u,c,utmp,ctmp,npop,npopmax,
     &     nlocd,nloch,nql,ncolt,nalmax,npp,nppmax,yy,z,ql,f,s,  
     &     xlim,ylim,indcell,distcell,indcelltmp,distcelltmp,lambda,
     &     qtc,nqtc,meanqtc,sdqtc,usegeno2,usegeno1,useql,useqtc,
     &     ploidy)
      implicit none 
      integer nindiv,nlocd,nloch,nql,ncolt,npop,npopmax,
     &     nalmax,npp,nppmax,yy(nindiv,2*nlocd+2*nloch),z(nindiv,nloch),
     &     ql(nindiv,nql),c(nppmax),
     &     indcell(nindiv),nqtc,usegeno2,usegeno1,useql,useqtc,
     &     ploidy
      double precision u(2,nindiv),f(npopmax,ncolt,nalmax),xlim(2),
     &     ylim(2),s(2,nindiv),distcell(nindiv),lambda,qtc(nindiv,nqtc),
     &     meanqtc(npopmax,nqtc),sdqtc(npopmax,nqtc)

      integer ctmp(nppmax),indcelltmp(nindiv),ipp,npptmp,
     &     iindiv,ipprem
      double precision utmp(2,nppmax),distcelltmp(nindiv),ggrunif,
     &     lrallvar2,r,alpha,ggrbinom,accept,b
      
*     naissance ou mort ?
      b = ggrbinom(1.d0,0.5d0)

      if(b .eq. 1) then
         if(npp .ne. nppmax) then 
*     naissance
            do ipp = 1,npp
               utmp (1,ipp) = u(1,ipp)
               utmp (2,ipp) = u(2,ipp)
               ctmp(ipp) = c(ipp)
            enddo
            npptmp = npp + 1
            ctmp(npptmp) = 1+ idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
            utmp(1,npptmp) = xlim(1)+(xlim(2)-xlim(1))*
     &           ggrunif(0.d0,1.d0)
            utmp(2,npptmp) = ylim(1)+(ylim(2)-ylim(1))*
     &           ggrunif(0.d0,1.d0)
            if(nppmax .gt. npptmp) then
               do ipp=npptmp+1,nppmax
                  ctmp(ipp) = -999
                  utmp(1,ipp) = -999.
                  utmp(2,ipp) = -999.
               enddo
            endif
            
            call voradd(s,utmp,indcell,distcell,indcelltmp,
     &           distcelltmp,nindiv,npp,nppmax)
c$$$            r = dexp(lratiogq(yy,f,c,ctmp,indcell,indcelltmp,npopmax,
c$$$     &           nlocmax,nalmax,nindiv,nloc,nlocmax2,nppmax,ploidy,
c$$$     &           qtc,nqtc,meanqtc,sdqtc,meanqtc,sdqtc,
c$$$     &           usegeno2,useqtc))
            r = dexp(lrallvar2(yy,z,ql,f,c,ctmp,indcell,indcelltmp,
     &           npopmax,nlocd,nloch,nql,ncolt,nalmax,
     &           nindiv,nppmax,qtc,nqtc,meanqtc,sdqtc,meanqtc,sdqtc,
     &           usegeno2,usegeno1,useql,useqtc,ploidy))
            r = r*lambda/dble(npp+1)
            alpha = dmin1(1.d0,r)
            accept = ggrbinom(1.d0,alpha)
            if(accept .eq. 1) then 
               npp = npptmp
               do iindiv=1,nindiv
                  indcell(iindiv) = indcelltmp(iindiv)
                  distcell(iindiv) = distcelltmp(iindiv)
               enddo
               do ipp = 1,nppmax
                  u (1,ipp) = utmp(1,ipp)
                  u (2,ipp) = utmp(2,ipp)
                  c(ipp) = ctmp(ipp)
               enddo
            endif
         endif
      else
*     mort
         if(npp .ne. 1) then 
            ipprem = 1+ aint(dble(npp)*ggrunif(0.d0,1.d0))
            if(ipprem .ne. 1) then 
               do ipp = 1,ipprem-1
                  utmp (1,ipp) = u(1,ipp)
                  utmp (2,ipp) = u(2,ipp)
                  ctmp(ipp) = c(ipp)
               enddo
            endif
            if(ipprem .ne. npp) then 
               do ipp = ipprem,npp-1
                  utmp (1,ipp) = u(1,ipp+1)
                  utmp (2,ipp) = u(2,ipp+1)
                  ctmp(ipp) = c(ipp+1)
               enddo
            endif
            do ipp=npp,nppmax
               utmp (1,ipp) = -999.
               utmp (2,ipp) = -999.
               ctmp(ipp) = -999
            enddo

            call vorrem(s,utmp,ipprem,indcell,distcell,
     &           indcelltmp,distcelltmp,nindiv,npp,nppmax)

c$$$            r = dexp(lratiogq(yy,f,c,ctmp,indcell,indcelltmp,npopmax,
c$$$     &           nlocmax,nalmax,nindiv,nloc,nlocmax2,nppmax,ploidy,
c$$$     &           qtc,nqtc,meanqtc,sdqtc,meanqtc,sdqtc,
c$$$     &           usegeno2,useqtc))
            r = dexp(lrallvar2(yy,z,ql,f,c,ctmp,indcell,indcelltmp,
     &           npopmax,nlocd,nloch,nql,ncolt,nalmax,
     &           nindiv,nppmax,qtc,nqtc,meanqtc,sdqtc,meanqtc,sdqtc,
     &           usegeno2,usegeno1,useql,useqtc,ploidy))
            r = r*dble(npp)/lambda
            alpha = dmin1(1.d0,r)
            accept = ggrbinom(1.d0,alpha)
            if(accept .eq. 1) then 
               npp = npp-1
               do iindiv=1,nindiv
                  indcell(iindiv) = indcelltmp(iindiv)
                  distcell(iindiv) = distcelltmp(iindiv)
               enddo
               do ipp = 1,nppmax
                  u (1,ipp) = utmp(1,ipp)
                  u (2,ipp) = utmp(2,ipp)
                  c(ipp) = ctmp(ipp)
               enddo
            endif
         endif
      endif
      end subroutine bdcellallvar2
************************************************************************




************************************************************************
*     Likelihood ratio p(zz|theta*)/p(zz|theta)
*     when data consist of genotypes only
      double precision function ratio(zz,f,c,ctmp,indcell,indcelltmp,
     &     npopmax,nlocmax,nalmax,nindiv,nloc,nlocmax2,
     &     nppmax,ploidy)
      implicit none
      integer npopmax,nlocmax,nalmax,nindiv,nloc,nlocmax2,
     &     nppmax,zz(nindiv,nlocmax2),c(nppmax),ctmp(nppmax),
     &     indcell(nindiv),indcelltmp(nindiv),ploidy
      double precision f(npopmax,nlocmax,nalmax)
      integer iindiv,iloc,ial1,ial2,ipop,ipoptmp

      ratio = 1.
      do iindiv=1,nindiv
         ipop = c(indcell(iindiv))
         ipoptmp = ctmp(indcelltmp(iindiv))
         do iloc=1,nloc
            ial1 = zz(iindiv,2*iloc-1)
            ial2 = zz(iindiv,2*iloc)
            if(ial1 .ne. -999) then 
               ratio = ratio*
     &              (f(ipoptmp,iloc,ial1)/f(ipop,iloc,ial1))
            endif
            if(ial2 .ne. -999) then 
               ratio = ratio*
     &              (f(ipoptmp,iloc,ial2)/f(ipop,iloc,ial2))
            endif
         enddo
      enddo
      if(ploidy .eq. 1) then 
         ratio = dsqrt(ratio)
      endif
      end function ratio


************************************************************************
*     Log of likelihood ratio p(zz|theta*)/p(zz|theta)
*     when data consist of genotypes and/or quantitative variables
      double precision function lratiogq(zz,f,c,ctmp,indcell,indcelltmp,
     &     npopmax,nlocmax,nalmax,nindiv,nloc,nlocmax2,nppmax,ploidy,
     &     qtc,nqtc,meanqtc,sdqtc,meanqtctmp,sdqtctmp,
     &     usegeno2,useqtc)
      implicit none
      integer npopmax,nlocmax,nalmax,nindiv,nloc,nlocmax2,
     &     nppmax,zz(nindiv,nlocmax2),c(nppmax),ctmp(nppmax),
     &     indcell(nindiv),indcelltmp(nindiv),ploidy,nqtc,
     &     usegeno2,useqtc
      double precision f(npopmax,nlocmax,nalmax),qtc(nindiv,nqtc),
     &     meanqtc(npopmax,nqtc),sdqtc(npopmax,nqtc),
     &     meanqtctmp(npopmax,nqtc),sdqtctmp(npopmax,nqtc)
      integer iindiv,iloc,ial1,ial2,ipop,ipoptmp,iqtc
c      write(*,*) 'inside lratiogq meanqtc=',meanqtc
c      write(*,*) 'inside lratiogq meanqtctmp=',meanqtctmp
c      write(*,*) 'inside lratiogq sdqtc=',sdqtc
c      write(*,*) 'inside lratiogq sdqtctmp=',sdqtctmp
c      write(*,*) 'inside lratiogq c=',c
c      write(*,*) 'inside lratiogq ctmp=',ctmp

      lratiogq = 0.
      do iindiv=1,nindiv
         ipop = c(indcell(iindiv))
         ipoptmp = ctmp(indcelltmp(iindiv))
         if(usegeno2 .eq. 1) then 
*     contrib genotypes
            do iloc=1,nloc
               ial1 = zz(iindiv,2*iloc-1)
               ial2 = zz(iindiv,2*iloc)
               if(ial1 .ne. -999) then 
                  lratiogq = lratiogq + 
     &                 dlog(f(ipoptmp,iloc,ial1)) - 
     &                 dlog(f(ipop,iloc,ial1))
               endif
               if(ploidy .eq. 2) then 
                  if(ial2 .ne. -999) then 
                     lratiogq = lratiogq + 
     &                    dlog(f(ipoptmp,iloc,ial2)) - 
     &                    dlog(f(ipop,iloc,ial2))
                  endif
               endif
            enddo
         endif
         if(useqtc .eq. 1) then
*     contrib quant. var.
            do iqtc = 1,nqtc
               lratiogq = lratiogq + 
     &              dlog(sdqtc(ipop,iqtc)) -
     &              dlog(sdqtctmp(ipoptmp,iqtc)) -
     &              0.5*(
     &              ((qtc(iindiv,iqtc)-meanqtctmp(ipoptmp,iqtc))/
     &              sdqtctmp(ipoptmp,iqtc))**2 -
     &              ((qtc(iindiv,iqtc)-meanqtc(ipop,iqtc))/
     &              sdqtc(ipop,iqtc))**2)
c               write(*,*) 'inside lratiogq=',lratiogq
            enddo
         endif
      enddo
      end function lratiogq
************************************************************************


************************************************************************
*     Log of likelihood ratio p(data|theta*)/p(data|theta)
*     when data consist of genotypes 
*     and/or quantitative variables
      double precision function lrallvar2(yy,z,ql,
     &     f,c,ctmp,indcell,indcelltmp,
     &     npopmax,nlocd,nloch,nql,ncolt,nalmax,nindiv,nppmax,
     &     qtc,nqtc,meanqtc,sdqtc,meanqtctmp,sdqtctmp,
     &     usegeno2,usegeno1,useql,useqtc,ploidy)
      implicit none
      integer npopmax,nlocd,nlocd2,nloch,nql,ncolt,nalmax,nindiv,
     &     nppmax,yy(nindiv,2*nlocd+2*nloch),z(nindiv,nloch),
     &     ql(nindiv,nql),c(nppmax),ctmp(nppmax),
     &     indcell(nindiv),indcelltmp(nindiv),nqtc,
     &     usegeno2,usegeno1,useql,useqtc,ploidy
      double precision f(npopmax,ncolt,nalmax),qtc(nindiv,nqtc),
     &     meanqtc(npopmax,nqtc),sdqtc(npopmax,nqtc),
     &     meanqtctmp(npopmax,nqtc),sdqtctmp(npopmax,nqtc)
      integer iindiv,iloc,ial1,ial2,ipop,ipoptmp,iqtc
c      write(*,*) 'inside lratiogq meanqtc=',meanqtc
c      write(*,*) 'inside lratiogq meanqtctmp=',meanqtctmp
c      write(*,*) 'inside lratiogq sdqtc=',sdqtc
c      write(*,*) 'inside lratiogq sdqtctmp=',sdqtctmp
c      write(*,*) 'inside lratiogq c=',c
c      write(*,*) 'inside lratiogq ctmp=',ctmp

      lrallvar2 = 0.
      do iindiv=1,nindiv
         ipop = c(indcell(iindiv))
         ipoptmp = ctmp(indcelltmp(iindiv))
         if(usegeno2 .eq. 1) then 
*     contrib diploid genotypes 
            do iloc=1,nlocd
               ial1 = yy(iindiv,2*iloc-1)
               ial2 = yy(iindiv,2*iloc)
               if(ial1 .ne. -999) then 
                  lrallvar2 = lrallvar2 + 
     &                 dlog(f(ipoptmp,iloc,ial1)) - 
     &                 dlog(f(ipop,iloc,ial1))
               endif
               if(ial2 .ne. -999) then 
                  lrallvar2 = lrallvar2 + 
     &                 dlog(f(ipoptmp,iloc,ial2)) - 
     &                 dlog(f(ipop,iloc,ial2))
               endif
            enddo
         endif
         if(usegeno1 .eq. 1) then 
*     contrib haploid genotypes 
            if(ploidy .eq. 2) then 
               do iloc=1,nloch
                  ial1 = yy(iindiv,2*nlocd+2*iloc-1)
                  ial2 = yy(iindiv,2*nlocd+2*iloc)
                  if(ial1 .ne. -999) then 
                     lrallvar2 = lrallvar2 + 
     &                    dlog(f(ipoptmp,nlocd+iloc,ial1)) - 
     &                    dlog(f(ipop,nlocd+iloc,ial1))
                  endif
                  if(ial2 .ne. -999) then 
                     lrallvar2 = lrallvar2 + 
     &                    dlog(f(ipoptmp,nlocd+iloc,ial2)) - 
     &                    dlog(f(ipop,nlocd+iloc,ial2))
                  endif
               enddo
            endif
            if(ploidy .eq. 1) then 
               do iloc=1,nloch
                  ial1 = z(iindiv,iloc)
                  if(ial1 .ne. -999) then 
                     lrallvar2 = lrallvar2 + 
     &                    dlog(f(ipoptmp,nlocd+iloc,ial1)) - 
     &                    dlog(f(ipop,nlocd+iloc,ial1))
                  endif
               enddo
            endif
         endif
         if(useql .eq. 1) then 
*     contrib qualit. variables
            do iloc=1,nql
               ial1 = ql(iindiv,iloc)
               if(ial1 .ne. -999) then 
                  lrallvar2 = lrallvar2 + 
     &                 dlog(f(ipoptmp,nlocd+nloch+iloc,ial1)) - 
     &                 dlog(f(ipop,nlocd+nloch+iloc,ial1))
               endif
            enddo
         endif
         if(useqtc .eq. 1) then
*     contrib quant. var.
            do iqtc = 1,nqtc
               if(qtc(iindiv,iqtc) .ne. -999) then 
                  lrallvar2 = lrallvar2 + 
     &                 dlog(sdqtc(ipop,iqtc)) -
     &                 dlog(sdqtctmp(ipoptmp,iqtc)) -
     &                 0.5*(
     &                 ((qtc(iindiv,iqtc)-meanqtctmp(ipoptmp,iqtc))/
     &                 sdqtctmp(ipoptmp,iqtc))**2 -
     &                 ((qtc(iindiv,iqtc)-meanqtc(ipop,iqtc))/
     &                 sdqtc(ipop,iqtc))**2)
c     write(*,*) 'inside lrallvar2=',lrallvar2
               endif
            enddo
         endif
      enddo
      end function lrallvar2
************************************************************************




************************************************************************
*     Indice des cellules dans une pop
      subroutine who(c,ipop,npp,nppmax,cellpop,
     &     ncellpop)
      implicit none
      integer npp,nppmax,c(nppmax),ipop,cellpop(nppmax),ncellpop
      integer ipp,ii
c      write(*,*) 'who'
      ii = 1
      ncellpop = 0
      do ipp=1,npp
         if(c(ipp) .eq. ipop) then
           cellpop(ii) = ipp
           ncellpop = ncellpop + 1
           ii = ii + 1
        endif
      enddo
      if(nppmax .gt. ncellpop) then
         do ipp= ncellpop+1, nppmax
            cellpop(ipp) = -999
         enddo
      endif
      end subroutine who
************************************************************************

************************************************************************
*     Tirage de nu cellules parmi ncellpop cellules
*      
      subroutine sample(cellpop,nppmax,nu,ncellpop,listcell)
      implicit none
      integer nppmax,cellpop(nppmax),nu,ncellpop,listcell(nppmax)
      integer isamp,ii,jj
      double precision ggrunif
c      write(*,*) 'sample'

*     init
      ii = 1 + idint(dint(dble(ncellpop)*ggrunif(0.d0,1.d0)))
      listcell(1) = cellpop(ii)
      if(nu .gt. 1) then
         do isamp = 2,nu
c             write(*,*) 'cellpop=',cellpop
*     translation
            if(ii .eq. 1) then 
               do jj=2,ncellpop
                  cellpop(jj-1) = cellpop(jj)
               enddo
            else
               if(ii .ne. ncellpop) then
                  do jj=ii+1,ncellpop
                     cellpop(jj-1) = cellpop(jj)
                  enddo
               endif
            endif
            cellpop(ncellpop-isamp+1) = -999
*     tirage parmi les ncellpop-isamp cellules restantes
            ii = 1 + 
     &           idint(dint(dble(ncellpop-isamp)*ggrunif(0.d0,1.d0)))
            listcell(isamp) = cellpop(ii)
         enddo
      endif
c      write(*,*) 'nu=',nu
c      write(*,*) 'ncellpop=',ncellpop
c      write(*,*) 'cellpop=',cellpop
c      write(*,*) 'listcell=',listcell
      end subroutine sample
************************************************************************


***********************************************************************
*
*     Tirage de nu cellules parmi ncellpop cellules
*      version corrigée de sample apres un bug 
*     trouvé en septembre 2005 à Göteborg
      subroutine sample2(cellpop,nppmax,nu,ncellpop,listcell)
      implicit none
      integer nppmax,cellpop(nppmax),nu,ncellpop,listcell(nppmax)
      integer isamp,ii,jj
      double precision ggrunif
c      call rndstart()
c      write(*,*) 'sample2'
c      write(*,*) 'nu=',nu
c      write(*,*) 'ncellpop=',ncellpop
c      write(*,*) 'cellpop=',cellpop
*     init
      ii = 1 + idint(dint(dble(ncellpop)*ggrunif(0.d0,1.d0)))
      listcell(1) = cellpop(ii)
c      write(*,*) 'listcell(1)=',listcell(1)
      if(nu .gt. 1) then
         do isamp = 2,nu
c             write(*,*) 'cellpop=',cellpop
*     translation
            if(ii .eq. 1) then 
               do jj=2,ncellpop
                  cellpop(jj-1) = cellpop(jj)
               enddo
            else
               if(ii .ne. ncellpop) then
                  do jj=ii+1,ncellpop
                     cellpop(jj-1) = cellpop(jj)
                  enddo
               endif
            endif
            cellpop(ncellpop-(isamp-1)+1) = -999
c             write(*,*) 'cellpop=',cellpop
*     tirage parmi les ncellpop-isamp cellules restantes
            ii = 1 + 
     &           idint(dint(dble(ncellpop-isamp)*ggrunif(0.d0,1.d0)))
c            write(*,*) 'ii=',ii
            listcell(isamp) = cellpop(ii)
c            write(*,*) 'listcell(isamp)=',listcell(isamp)
         enddo
      endif
c      write(*,*) 'nu=',nu
c      write(*,*) 'ncellpop=',ncellpop
c      write(*,*) 'cellpop=',cellpop
c      write(*,*) 'listcell=',listcell
c      call rndend()
      end subroutine sample2
************************************************************************

*******************************************************************
*     split d'une pop en deux
*     reallocation de nu cellules dont les indices
*     sont dans listcell
*     dans la pop ipop
      subroutine split(ipop,c,ctmp,nppmax,nu,listcell)
      implicit none
      integer ipop,nppmax,c(nppmax),ctmp(nppmax),nu,
     &     listcell(nppmax)
      integer ipp,ii
c      write(*,*) 'debut de split'
c      write(*,*) 'nu=',nu
c      write(*,*) 'ipop=',ipop
      do ipp=1,nppmax
         ctmp(ipp) = c(ipp)
      enddo
      if(nu .gt. 0) then
         do ii=1,nu
            ctmp(listcell(ii)) = ipop
         enddo
      endif
c      write(*,*)'c=',c
c      write(*,*)'ctmp=',ctmp
c      write(*,*) 'fin de split'
      end subroutine split
***********************************************************************


*******************************************************************
*     merge de deux  pops en une : 
*     reallocation des nu cellules de la pop ipoprem 
*     dont les indices sont dans listcell
*     dans la pop ipophost
      subroutine merging(ipoprem,ipophost,
     &     c,ctmp,nppmax,nu,listcell)
      implicit none
      integer ipoprem,ipophost,nppmax,
     &     c(nppmax),ctmp(nppmax),nu,listcell(nppmax)
      integer ipp,ii
c      write(*,*) 'debut de merge'
c      write(*,*) 'nu=',nu
c      write(*,*) 'ipoprem=',ipoprem
      do ipp=1,nppmax
         ctmp(ipp) = c(ipp)
      enddo
      if(ipoprem .gt. ipophost) then 
         if(nu .gt. 0) then
            do ii=1,nu
               ctmp(listcell(ii)) = ipophost
            enddo
         endif
      else
         if(nu .gt. 0) then
            do ii=1,nu
               ctmp(listcell(ii)) = ipophost - 1
            enddo
         endif
      endif
      do ipp=1,nppmax
         if(c(ipp) .gt. ipoprem) ctmp(ipp) = c(ipp)-1
      enddo
c      write(*,*)'c=',c
c      write(*,*)'ctmp=',ctmp
c      write(*,*) 'fin de merge'
      end subroutine merging
*******************************************************************
 

********************************************************************
*     Mise a jour de c et f en cas d acceptation d'un split/merge
      subroutine accept5(nppmax,npopmax,nlocmax,nalmax,
     &     nal,c,ctmp,f,ftmp,drift,drifttmp)
      implicit none
      integer nppmax,npopmax,nlocmax,nalmax,
     &     nal(nlocmax),c(nppmax),ctmp(nppmax)
      double precision f(npopmax,nlocmax,nalmax),
     &     ftmp(npopmax,nlocmax,nalmax),
     &     drift(npopmax),drifttmp(npopmax)
      integer ipop,iloc,ial,ipp
c      write(*,*) 'debut de accept5'
c      write(*,*) 'f=',f
c      write(*,*) 'ftmp=',ftmp
      do ipp=1,nppmax
         c(ipp) = ctmp(ipp)
      enddo
      do ipop = 1,npopmax
         do iloc= 1,nlocmax
            do ial=1,nal(iloc)
               f(ipop,iloc,ial) = ftmp(ipop,iloc,ial)
            enddo
         enddo
         drift(ipop) = drifttmp(ipop)
      enddo
c      write(*,*) 'f=',f
c      write(*,*) 'fin de accept5'
      end subroutine accept5
********************************************************************



********************************************************************
*     coefficients du binome C_n^p
*
      double precision function bico(n,p)
      implicit none
      integer n,p
      double precision gglgamfn
      bico = dexp(gglgamfn(dble(n+1))-gglgamfn(dble(p+1))-
     &     gglgamfn(dble(n-p+1)))
c      write(*,*) 'in bico '
c$$$      write(*,*) 'n=', n
c$$$      write(*,*) 'p=', p
c$$$      write(*,*) 'gglgamfn(dble(n+1))=',gglgamfn(dble(n+1))
c$$$      write(*,*) 'gglgamfn(dble(p+1))=',gglgamfn(dble(p+1))
c$$$      write(*,*) 'gglgamfn(dble(n-p+1)))=',gglgamfn(dble(n-p+1))
c$$$      write(*,*) 'dexp()=',dexp(gglgamfn(dble(n+1))-gglgamfn(dble(p+1))-
c$$$     &     gglgamfn(dble(n-p+1)))
c$$$      write(*,*) 'bico =', nint(dexp(gglgamfn(dble(n+1))-
c$$$     &     gglgamfn(dble(p+1))-
c$$$     &     gglgamfn(dble(n-p+1))))
c      write(*,*) 'bico =',bico
      
      end function bico
***********************************************************************




***********************************************************************
*     ln du coefficient du binome C_n^p
*
      double precision function lbico(n,p)
      implicit none
      integer n,p
      double precision gglgamfn
      lbico = gglgamfn(dble(n+1))-gglgamfn(dble(p+1))-
     &     gglgamfn(dble(n-p+1))
      end function lbico
***********************************************************************





************************************************************************
*     log du ratio des vraisemblances dans bdpop6
*
      double precision function llr6(zz,f,ftmp,c,ctmp,indcell,
     &     indcelltmp,npopmax,nlocmax,nalmax,nindiv,nloc,nlocmax2,
     &     nppmax,ploidy)
      implicit none
      integer npopmax,nlocmax,nalmax,nindiv,nloc,nlocmax2,
     &     nppmax,zz(nindiv,nlocmax2),c(nppmax),ctmp(nppmax),
     &     indcell(nindiv),indcelltmp(nindiv),ploidy
      double precision f(npopmax,nlocmax,nalmax),
     &     ftmp(npopmax,nlocmax,nalmax)
      integer iindiv,iloc,ial1,ial2,ipop,ipoptmp

      llr6 = 0

*     log du rapport des vraisemblances
      do iindiv=1,nindiv
         ipop = c(indcell(iindiv))
         ipoptmp = ctmp(indcelltmp(iindiv))
         do iloc=1,nloc
            ial1 = zz(iindiv,2*iloc-1)
            ial2 = zz(iindiv,2*iloc)
            if(ial1 .ne. -999) then 
               llr6 = llr6 + 
     &              dlog(ftmp(ipoptmp,iloc,ial1)) - 
     &              dlog(f(ipop,iloc,ial1))
            endif
            if(ial2 .ne. -999) then 
               llr6 = llr6 + 
     &              dlog(ftmp(ipoptmp,iloc,ial2)) - 
     &              dlog(f(ipop,iloc,ial2))
            endif
            
c$$$            if(ipoptmp .eq. 3) then 
c$$$               write(*,*) 'ipoptmp=',ipoptmp
c$$$               write(*,*) 'iinidiv=',iindiv
c$$$               write(*,*) 'c=', c(indcell(iindiv))
c$$$               write(*,*) 'ctmp=', ctmp(indcelltmp(iindiv))
c$$$               write(*,*) 'zz=',zz(iindiv,1),zz(iindiv,2)
c$$$               write(*,*) 'llr6 =',llr6
c$$$               write(*,*) 'ftmp(ipoptmp,iloc,ial1)=',
c$$$     &              ftmp(ipoptmp,iloc,ial1) 
c$$$               write(*,*) 'f(ipop,iloc,ial1)=',
c$$$     &              f(ipop,iloc,ial1)
c$$$               write(*,*) 'ftmp(ipoptmp,iloc,ial2)=',
c$$$     &              ftmp(ipoptmp,iloc,ial2)
c$$$               write(*,*) 'f(ipop,iloc,ial2)=',
c$$$     &              f(ipop,iloc,ial2)
c$$$            endif

         enddo
      enddo
      if(ploidy .eq. 1) llr6 = 0.5d0*llr6 
      end function llr6
************************************************************************




************************************************************************
*     comptage des alleles dans chaque pop pour c
      subroutine countn(nindiv,nlocmax,nlocmax2,npopmax,
     &     nppmax,nal,nalmax,zz,n,indcell,c,ploidy)
      implicit none
      integer  nindiv,nlocmax,nlocmax2,npopmax,nppmax,nalmax,ploidy,
     &     zz(nindiv,nlocmax2),nal(nlocmax),
     &     n(npopmax,nlocmax,nalmax),c(nppmax),indcell(nindiv)
      integer ipop,iloc,ial,iindiv
*     init du tableau
      do ipop = 1,npopmax
         do iloc = 1,nlocmax
            do ial =1,nal(iloc)
               n(ipop,iloc,ial)=0
            enddo
         enddo
      enddo
*     comptage
      do iindiv = 1,nindiv
         do iloc = 1,nlocmax
            if(zz(iindiv,2*iloc-1) .ne. -999) then
               n(c(indcell(iindiv)),iloc,zz(iindiv,2*iloc-1)) = 
     &              n(c(indcell(iindiv)),iloc,zz(iindiv,2*iloc-1))+ 1 
            endif
            if(ploidy .eq. 2) then
               if(zz(iindiv,2*iloc) .ne. -999) then 
                  n(c(indcell(iindiv)),iloc,zz(iindiv,2*iloc)) = 
     &                 n(c(indcell(iindiv)),iloc,zz(iindiv,2*iloc)) + 1 
               endif
            endif
         enddo
      enddo
      end subroutine countn
************************************************************************




************************************************************************
*     comptage des alleles dans chaque pop pour c
      subroutine countnallvar2(nindiv,nlocd,nloch,nql,ncolt,
     &     nppmax,nal,nalmax,yy,z,ql,n,indcell,c,npopmax,
     &     usegeno2, usegeno1,useql,ploidy)
      implicit none
      integer  nindiv,nlocd,nloch,nql,ncolt,nppmax,nalmax,
     &     npopmax,yy(nindiv,2*nlocd+2*nloch),z(nindiv,nloch),
     &     ql(nindiv,nql),nal(ncolt),n(npopmax,ncolt,nalmax),c(nppmax),
     &     indcell(nindiv),usegeno2, usegeno1,useql,ploidy
      integer ipop,iloc,ial,iindiv
c$$$      write(*,*) 'debut countnallvar'
c$$$      write(*,*) 'npopmax=',npopmax
c$$$      write(*,*) 'nindiv =', nindiv,"\n"
c$$$      write(*,*) 'nlocd =',nlocd   ,"\n"
c$$$      write(*,*) 'nloch =',nloch   ,"\n"
c$$$      write(*,*) 'nql =',nql   ,"\n"
c$$$      write(*,*) 'ncolt =',ncolt   ,"\n"
c$$$      write(*,*) 'nal =',nal ,"\n"
c$$$      write(*,*) 'nalmax =', nalmax,"\n"
c$$$      write(*,*) 'usegeno2=',usegeno2,"\n"
c$$$      write(*,*) 'usegeno1=',usegeno1,"\n"
c$$$      write(*,*) 'useql=',useql,"\n"
c$$$      write(*,*) 'yy =', yy,"\n"
c$$$      write(*,*) 'z =', z,"\n"
c$$$      write(*,*) 'ql =', ql,"\n"
c$$$      write(*,*) 'npopmax =', npopmax ,"\n" 
c$$$      write(*,*) 'indcell =', indcell ,"\n"
*     init du tableau
      do ipop = 1,npopmax
         do iloc = 1,ncolt
            do ial =1,nal(iloc)
               n(ipop,iloc,ial) = 0
            enddo
         enddo
      enddo
*     comptage
      do iindiv = 1,nindiv
c         write(*,*) 'iindiv=',iindiv
         if(usegeno2 .eq. 1) then
*     diploid geno
            do iloc = 1,nlocd
c               write(*,*) 'iloc=',iloc
               if(yy(iindiv,2*iloc-1) .ne. -999) then
                  n(c(indcell(iindiv)),iloc,yy(iindiv,2*iloc-1)) = 
     &            n(c(indcell(iindiv)),iloc,yy(iindiv,2*iloc-1))+ 1 
               endif
               if(yy(iindiv,2*iloc) .ne. -999) then 
                  n(c(indcell(iindiv)),iloc,yy(iindiv,2*iloc)) = 
     &            n(c(indcell(iindiv)),iloc,yy(iindiv,2*iloc)) + 1 
               endif
            enddo
         endif
         if(usegeno1 .eq. 1) then
*     haploid geno
            if(ploidy .eq. 2) then 
               do iloc = 1,nloch
c     write(*,*) 'iloc=',iloc
                  if(yy(iindiv,2*nlocd+2*iloc-1) .ne. -999) then
                     n(c(indcell(iindiv)),nlocd+iloc,
     &                    yy(iindiv,2*nlocd+2*iloc-1)) = 
     &               n(c(indcell(iindiv)),nlocd+iloc,
     &                    yy(iindiv,2*nlocd+2*iloc-1))+ 1 
                  endif
                  if(yy(iindiv,2*nlocd+2*iloc) .ne. -999) then 
                     n(c(indcell(iindiv)),nlocd+iloc,
     &                    yy(iindiv,2*nlocd+2*iloc)) = 
     &               n(c(indcell(iindiv)),nlocd+iloc,
     &                    yy(iindiv,2*nlocd+2*iloc)) + 1 
                  endif
               enddo
            endif
            if(ploidy .eq. 1) then 
               do iloc = 1,nloch
                  if(z(iindiv,iloc) .ne. -999) then
                     n(c(indcell(iindiv)),nlocd+iloc,z(iindiv,iloc)) = 
     &               n(c(indcell(iindiv)),nlocd+iloc,z(iindiv,iloc)) + 1 
                  endif
               enddo
            endif
         endif
         if(useql .eq. 1) then
*     ql variable
            do iloc = 1,nql
               if(ql(iindiv,iloc) .ne. -999) then
c                  write(*,*) 'ql(iindiv,iloc)=',ql(iindiv,iloc)
                  n(c(indcell(iindiv)),nlocd+nloch+iloc,
     &                 ql(iindiv,iloc)) = 
     &            n(c(indcell(iindiv)),nlocd+nloch+iloc,
     &                 ql(iindiv,iloc))+ 1 
               endif
            enddo
         endif
      enddo
c$$$      write(*,*) 'ds countnallvar c=',c(1),c(2),c(3),c(4)
c$$$      write(*,*) 'ds countnallvar n=', (n(ipop,1,1),ipop=1,npopmax)
c$$$      write(*,*) 'n=',n
      end subroutine countnallvar2
************************************************************************


************************************************************************
*     log du ratio (prob cond. complete)/prior
*     pour les frequences
*     dans un split de la pop ipop
      double precision function lrf(ipop,npopmax,nlocmax,nal,nalmax,
     &     f,fa,drift,n)
      implicit none
      integer ipop,npopmax,nlocmax,nal(nlocmax),nalmax,
     &     n(npopmax,nlocmax,nalmax)
      double precision f(npopmax,nlocmax,nalmax),
     &     fa(nlocmax,nalmax),
     &     drift(npopmax)
      integer iloc,ial,nn
      double precision ss,gglgamfn,q

      lrf = 0.
      q = (1-drift(ipop))/drift(ipop)
      do iloc = 1,nlocmax
         ss = 0
         nn = 0
         do ial = 1,nal(iloc)
            ss = ss + gglgamfn(fa(iloc,ial) * q +
     &           dble(n(ipop,iloc,ial))) +
     &           (1 - fa(iloc,ial) * q - dble(n(ipop,iloc,ial)))*
     &           dlog(f(ipop,iloc,ial))
            nn = nn + n(ipop,iloc,ial)
         enddo
c         write(*,*) 'nn=',nn
         lrf = lrf + gglgamfn(dble(nal(iloc))) -
     &        gglgamfn(q + nn) + ss

      enddo
      end function lrf
***********************************************************************


***********************************************************************
*
*     Naissance et mort de pop avec réallocations 
*     (split/merge)
*     proposition de drift* selon prior
*     proposition de f* selon conditionelle complète 
*     dans les deux sens
*
      subroutine bdpop9(npop,npopmin,npopmax,f,fa,drift,
     &     nloc,nlocmax,nlocmax2,
     &     nal,nalmax,indcell,nindiv,npp,nppmax,c,ctmp,
     &     a,ptmp,ftmp,drifttmp,zz,cellpop,listcell,
     &     cellpophost,n,ntmp,ploidy)
      implicit none
      integer npop,npopmin,npopmax,nloc,nlocmax,nal(nlocmax),
     &     nalmax,nindiv,npp,nppmax,indcell(nindiv),
     &     nlocmax2,c(nppmax),ctmp(nppmax),zz(nindiv,nlocmax2),
     &     n(npopmax,nloc,nalmax),ntmp(npopmax,nloc,nalmax),
     &     ploidy
      double precision f(npopmax,nlocmax,nalmax),drift(npopmax),
     &      ftmp(npopmax,nlocmax,nalmax),drifttmp(npopmax),
     &     a(nalmax),ptmp(nalmax),fa(nlocmax,nalmax)
      integer ipoprem,ipp,isplit,
     &     cellpop(nppmax),ncellpop,nu,listcell(nppmax),
     &     ipophost,ncellpophost,cellpophost(nppmax),ii
      double precision alpha,ggrunif,lbico,lratio,llr6,termf9
      double precision b,bern,ggrbinom      
       do ipp=1,nppmax
          cellpop(ipp) = -999
          listcell(ipp) = -999
       enddo
*     naissance ou mort ?
       b = ggrbinom(1.d0,0.5d0)
       
       if(b .eq. 1) then
          if(npop .lt. npopmax) then 
c             write(*,*) 'naissance'
*     split

*     choix de la pop qui split
             isplit = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
             
*     recherche des cellules affectees a cette pop
             call who(c,isplit,npp,nppmax,cellpop,ncellpop)
             
             if(ncellpop .gt. 0) then
*     tirage du nombre de cellules reallouees
                nu = idint(dint(dble(ncellpop)*ggrunif(0.d0,1.d0)))
                if(nu .gt. 0) then
                   
*     tirage des cellules reallouees
                   call sample2(cellpop,nppmax,nu,ncellpop,
     &                  listcell)
                   
*     proposition de reallocation dans la pop npop+1
                   call split(npop+1,c,ctmp,nppmax,nu,
     &                  listcell)
                else 
                   do ipp = 1,nppmax
                      ctmp(ipp) = c(ipp)
                   enddo
                endif
             else
                nu = 0
                do ipp = 1,nppmax
                   ctmp(ipp) = c(ipp)
                enddo
             endif

*     comptage des alleles sur chaque locus pour c puis ctmp
             call countn(nindiv,nlocmax,nlocmax2,npopmax,
     &            nppmax,nal,nalmax,zz,n,indcell,c,ploidy)
             call countn(nindiv,nlocmax,nlocmax2,npopmax,
     &            nppmax,nal,nalmax,zz,ntmp,indcell,ctmp,ploidy)
             
*     proposition nouvelle freq et derive 
c     call addfreq5(isplit,npop,npopmax,nloc,nlocmax,
c     &     nal,nalmax,f,ftmp,fa,drift,drifttmp,a,ptmp)
             call addfreq7(npop,npopmax,nloc,nlocmax,
     &            nal,nalmax,isplit,
     &            f,ftmp,fa,drift,drifttmp,a,ptmp,ntmp)
             
*     calcul du log du ratio
*     terme des vraisemblances
             lratio =  llr6(zz,f,ftmp,c,ctmp,indcell,
     &            indcell,npopmax,nlocmax,nalmax,
     &            nindiv,nloc,nlocmax2,nppmax,ploidy)
*     terme des freq.
c$$$             lratio = lratio + termfsplit(isplit,npop,npopmax,
c$$$     &            nlocmax,nal,nalmax,
c$$$     &            f,ftmp,n,ntmp,fa,drift,drifttmp) 
             lratio = lratio 
     &           + termf9(npopmax,nloc,nal,nalmax,n,f,fa,drift,isplit) 
     &            -termf9(npopmax,nloc,nal,nalmax,ntmp,ftmp,fa,
     &            drifttmp,isplit) 
     &            -termf9(npopmax,nloc,nal,nalmax,ntmp,ftmp,fa,
     &            drifttmp,npop+1) 

             
*     terme des proposal sur c
             lratio = lratio + dlog(2*dble(ncellpop+1)) + 
     &            lbico(ncellpop,nu) - dlog(dble(npop+1)) 
             
*     terme des priors sur c
             lratio = lratio + 
     &            dble(npp)*(dlog(dble(npop)) - 
     &            dlog(dble(npop+1)))
             
             lratio = dmin1(0.d0,lratio)
             alpha = dexp(lratio)
             bern = ggrbinom(1.d0,alpha)

c$$$                   write(*,*) 'npop=',npop
c$$$                   write(*,*) 'npp=',npp
c$$$                   write(*,*) 'isplit=',isplit
c$$$                   write(*,*) 'c=',c
c$$$                   write(*,*) 'ncellpop=',ncellpop  
c$$$                   write(*,*) 'nu=',nu
c$$$                   write(*,*) 'listcell=',listcell
c$$$                   write(*,*) 'ctmp=',ctmp 
c$$$                   do ipop=1,npopmax
c$$$                      do iloc=1,nlocmax
c$$$                         do ial=1,nal(iloc)
c$$$                            write(*,*) 
c$$$     &                           'n(',ipop,',',iloc,',',ial,')=',
c$$$     &                           n(ipop,iloc,ial)
c$$$                         enddo
c$$$                      enddo
c$$$                   enddo
*c$$$                   do ipop=1,npopmax
c$$$                      do iloc=1,nlocmax
c$$$                         do ial=1,nal(iloc)
c$$$                            write(*,*) 
c$$$     &                           'ntmp(',ipop,',',iloc,',',ial,')=',
c$$$     &                           ntmp(ipop,iloc,ial)
c$$$                         enddo
c$$$                      enddo
c$$$                   enddo
c$$$                   do ipop=1,npopmax
c$$$                      do iloc=1,nlocmax
c$$$                         do ial=1,nal(iloc)
c$$$                            write(*,*) 
c$$$     &                           'f(',ipop,',',iloc,',',ial,')=',
c$$$     &                       f(ipop,iloc,ial)
c$$$                         enddo
c$$$                      enddo
c$$$                   enddo
c$$$                   do ipop=1,npopmax
c$$$                      do iloc=1,nlocmax
c$$$                         do ial=1,nal(iloc)
c$$$                            write(*,*) 
c$$$     &                           'ftmp(',ipop,',',iloc,',',ial,')=',
c$$$     &                           ftmp(ipop,iloc,ial)
c$$$                         enddo
c$$$                      enddo
c$$$                   enddo
c$$$                   write(*,*) 'terme vrais. =',
c$$$     &                 llr6(zz,f,ftmp,c,ctmp,indcell,
c$$$     &                  indcell,npopmax,nlocmax,nalmax,
c$$$     &                  nindiv,nloc,nlocmax2,nppmax) 
c$$$                   write(*,*) 'terme freq =',
c$$$     &                  termfsplit(isplit,npop,npopmax,
c$$$     &                  nlocmax,nal,nalmax,
c$$$     &                  f,ftmp,n,ntmp,fa,drift,drifttmp)
c$$$                   write(*,*) ' term prop c=',
c$$$     &                  dlog(2*dble(ncellpop+1)) + 
c$$$     &                  lbico(ncellpop,nu) - dlog(dble(npop+1)) 
c$$$                   write(*,*) ' term prior c=',
c$$$     &                   dble(npp)*(dlog(dble(npop)) - 
c$$$     &                  dlog(dble(npop+1)))
c                    write(*,*) 'alpha=',alpha 

             
             if(bern .eq. 1) then
                call accept5(nppmax,npopmax,nlocmax,
     &               nalmax,nal,c,ctmp,f,ftmp,drift,drifttmp)
                npop = npop + 1
             endif
          endif

*     merge
      else
         if(npop .gt. npopmin) then 
c             write(*,*) 'mort'
*     tirage de la pop qui meurt
            ipoprem = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
            
*     tirage de la pop hote
            ipophost = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
            do while(ipophost .eq. ipoprem)
               ipophost = 1 + 
     &              idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
            enddo

*     on range dans la pop d'indice le plus petit
            if(ipophost .gt. ipoprem) then
               ii = ipophost
               ipophost = ipoprem
               ipoprem = ii
            endif
            
*     recherche des cellules qui vont etre reallouees
            call who(c,ipoprem,npp,nppmax,cellpop,ncellpop)
            
*     recherche des cellules de la pop hote
            call who(c,ipophost,npp,nppmax,cellpophost,
     &           ncellpophost)
 
*     next line corrected by Gilles on 05/01/08
*            if(ncellpop .gt. 0) then
*     proposition de reallocation dans la pop ipophost
               call merging(ipoprem,ipophost,c,ctmp,nppmax,
     &              ncellpop,cellpop)

*     comptage des alleles sur chaque locus pour c puis ctmp
               call countn(nindiv,nlocmax,nlocmax2,npopmax,
     &              nppmax,nal,nalmax,zz,n,indcell,c,ploidy)
               call countn(nindiv,nlocmax,nlocmax2,npopmax,
     &              nppmax,nal,nalmax,zz,ntmp,indcell,ctmp,ploidy)

*     propostion du nouveau tableau de freq et de derives
c               call remfreq5(ipoprem,ipophost,npop,npopmax,
c     &              nloc,nlocmax,nal,nalmax,f,ftmp,drift,drifttmp,
c     &              a,fa)
               call remfreq7(ipoprem,ipophost,
     &     npop,npopmax,nloc,nlocmax,nal,
     &     nalmax,f,ftmp,drift,drifttmp,fa,a,ptmp,ntmp)
               
*     calcul du log du ratio  
*     terme des vraisemblances
               lratio =  llr6(zz,f,ftmp,c,ctmp,indcell,
     &              indcell,npopmax,nlocmax,nalmax,
     &              nindiv,nloc,nlocmax2,nppmax,ploidy)
               
*     terme des freq.
c$$$               lratio = lratio + 
c$$$     &              termfmerge(ipophost,ipoprem,
c$$$     &              npopmax,nlocmax,
c$$$     &              nal,nalmax,
c$$$     &              f,ftmp,n,ntmp,fa,drift,drifttmp)
            lratio = lratio  
     &     + termf9(npopmax,nloc,nal,nalmax,n,f,fa,drift,ipophost) 
     &      + termf9(npopmax,nloc,nal,nalmax,n,f,fa,drift,ipoprem) 
     & -termf9(npopmax,nloc,nal,nalmax,ntmp,ftmp,fa,
     &              drifttmp,ipophost)


*     terme des proposal sur c
               lratio = lratio + dlog(dble(npop)) - 
     &              dlog(2*dble(ncellpop+ncellpophost+1)) -
     &              lbico(ncellpop+ncellpophost,ncellpop) 

*     terme des priors sur c
               lratio = lratio + 
     &              dble(npp)*(dlog(dble(npop)) - 
     &              dlog(dble(npop-1)))
               lratio = dmin1(0.d0,lratio)
               alpha  = dexp(lratio)
               bern = ggrbinom(1.d0,alpha)      
      
         
c$$$               write(*,*) 'npop=',npop
c$$$               write(*,*) 'npp=',npp
c$$$               write(*,*) 'ipoprem=',ipoprem
c$$$               write(*,*) 'cellpop=',cellpop
c$$$               write(*,*) 'ipophost=',ipophost
c$$$               write(*,*) 'c=',c
c$$$               write(*,*) 'ctmp=',ctmp 
c$$$               do ipop=1,npopmax
c$$$                  do iloc=1,nlocmax
c$$$                     do ial=1,nal(iloc)
c$$$                        write(*,*) 
c$$$     &                           'n(',ipop,',',iloc,',',ial,')=',
c$$$     &                       n(ipop,iloc,ial)
c$$$                     enddo
c$$$                  enddo
c$$$               enddo
c$$$               do ipop=1,npopmax
c$$$                  do iloc=1,nlocmax
c$$$                     do ial=1,nal(iloc)
c$$$                        write(*,*) 
c$$$     &                       'ntmp(',ipop,',',iloc,',',ial,')=',
c$$$     &                           ntmp(ipop,iloc,ial)
c$$$                     enddo
c$$$                  enddo
c$$$               enddo
c$$$               do ipop=1,npopmax
c$$$                  do iloc=1,nlocmax
c$$$                     do ial=1,nal(iloc)
c$$$                        write(*,*) 
c$$$     &                           'f(',ipop,',',iloc,',',ial,')=',
c$$$     &                       f(ipop,iloc,ial)
c$$$                     enddo
c$$$                  enddo
c$$$               enddo
c$$$               do ipop=1,npopmax
c$$$                  do iloc=1,nlocmax
c$$$                     do ial=1,nal(iloc)
c$$$                        write(*,*) 
c$$$     &                       'ftmp(',ipop,',',iloc,',',ial,')=',
c$$$     &                           ftmp(ipop,iloc,ial)
c$$$                     enddo
c$$$                  enddo
c$$$               enddo
c$$$               write(*,*) 'terme vrais. =',
c$$$     &              llr6(zz,f,ftmp,c,ctmp,indcell,
c$$$     &              indcell,npopmax,nlocmax,nalmax,
c$$$     &              nindiv,nloc,nlocmax2,nppmax)
c$$$               write(*,*) 'terme freq =',
c$$$     &              termfmerge(ipophost,ipoprem,
c$$$     &              npopmax,nlocmax,
c$$$     &              nal,nalmax,
c$$$     &              f,ftmp,n,ntmp,fa,drift,drifttmp)
c$$$               write(*,*) ' term prop c=',
c$$$     &              dlog(dble(npop)) - 
c$$$     &              dlog(2*dble(ncellpop+ncellpophost+1)) -
c$$$     &              lbico(ncellpop+ncellpophost,ncellpop) 
c$$$               write(*,*) ' term prior c=',
c$$$     &              dble(npp)*(dlog(dble(npop)) - 
c$$$     &              dlog(dble(npop-1)))
c               write(*,*) 'alpha=',alpha 


               if(bern .eq. 1) then
                  call accept5(nppmax,npopmax,nlocmax,
     &                 nalmax,nal,c,ctmp,f,ftmp,drift,drifttmp)
                  npop = npop - 1
               endif
*            endif
         endif
      endif
      end subroutine bdpop9
***********************************************************************



***********************************************************************
*     split/merge populations in the spatial D-model
*     changes from bdpop7bis:
*     - process populations whatever the number of tiles or individuals
*       they have
      subroutine bdpop9bis(npop,npopmin,npopmax,f,fa,drift,
     &     nloc,nlocmax,nlocmax2,
     &     nal,nalmax,indcell,nindiv,npp,nppmax,c,ctmp,
     &     a,ptmp,ftmp,drifttmp,zz,cellpop,listcell,
     &     cellpophost,n,ntmp,ploidy)
      implicit none
      integer npop,npopmin,npopmax,nloc,nlocmax,nal(nlocmax),
     &     nalmax,nindiv,npp,nppmax,indcell(nindiv),
     &     nlocmax2,c(nppmax),ctmp(nppmax),zz(nindiv,nlocmax2),
     &     n(npopmax,nloc,nalmax),ntmp(npopmax,nloc,nalmax),
     &     ploidy
      double precision f(npopmax,nlocmax,nalmax),drift(npopmax),
     &      ftmp(npopmax,nlocmax,nalmax),drifttmp(npopmax),
     &     a(nalmax),ptmp(nalmax),fa(nlocmax,nalmax)
      integer ipoprem,ipp,isplit,
     &     cellpop(nppmax),ncellpop,nu,listcell(nppmax),
     &     ipophost,ncellpophost,cellpophost(nppmax),ii
      double precision alpha,ggrunif,lbico,lratio,llr6,termf9bis,
     &     gglgamfn
      integer iloc
      double precision bern,ggrbinom,b

c      write(*,*) ''

      do ipp=1,nppmax
         cellpop(ipp) = -999
         listcell(ipp) = -999
      enddo
*     naissance ou mort ?
      b = ggrbinom(1.d0,0.5d0)
      
      if(b .eq. 1) then
         if(npop .lt. npopmax) then 
c            write(*,*) 'naissance'
*     split
            
*     choix de la pop qui split
            isplit = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
c            write(*,*) 'isplit=',isplit

*     recherche des cellules affectees a cette pop
            call who(c,isplit,npp,nppmax,cellpop,ncellpop)

            if(ncellpop .gt. 0) then
*     tirage du nombre de cellules reallouees 
               nu = idint(dint(dble(ncellpop+1)*ggrunif(0.d0,1.d0)))
               if(nu .gt. 0) then                 

*     tirage des cellules reallouees
                  call sample2(cellpop,nppmax,nu,ncellpop,
     &                 listcell)
 
*     proposition de reallocation dans la pop npop+1
                  call split(npop+1,c,ctmp,nppmax,nu,
     &                 listcell)
               else 
                  do ipp = 1,nppmax
                     ctmp(ipp) = c(ipp)
                  enddo
               endif
            else
               nu = 0
               do ipp = 1,nppmax
                  ctmp(ipp) = c(ipp)
               enddo
            endif
            
*     comptage des alleles sur chaque locus pour c puis ctmp
            call countn(nindiv,nlocmax,nlocmax2,npopmax,
     &           nppmax,nal,nalmax,zz,n,indcell,c,ploidy)
            call countn(nindiv,nlocmax,nlocmax2,npopmax,
     &           nppmax,nal,nalmax,zz,ntmp,indcell,ctmp,ploidy)

            
*     proposition nouvelle freq et derive 
            call addfreq7bis(npop,npopmax,nloc,
     &           nlocmax,nal,nalmax,isplit,
     &           f,ftmp,fa,drift,drifttmp,a,ptmp,ntmp)

*     calcul du log du ratio
*     terme des vraisemblances
c            write(*,*) 'calcul du log du ratio'
            lratio =  llr6(zz,f,ftmp,c,ctmp,indcell,
     &           indcell,npopmax,nlocmax,nalmax,
     &           nindiv,nloc,nlocmax2,nppmax,ploidy)

*     term proposal freq.
            lratio = lratio 
     &           + termf9bis(npopmax,nloc,nal,nalmax,n,f,isplit) 
     &         - termf9bis(npopmax,nloc,nal,nalmax,ntmp,ftmp,isplit) 
     &         - termf9bis(npopmax,nloc,nal,nalmax,ntmp,ftmp,npop+1)  

* term prior freq
            do iloc = 1,nloc
               lratio = lratio + gglgamfn(dble(nal(iloc)))
            enddo

*     terme des proposal sur c
            lratio = lratio + dlog(2*dble(ncellpop+1)) + 
     &           lbico(ncellpop,nu) - dlog(dble(npop+1)) 
*     terme des priors sur c
            lratio = lratio + 
     &           dble(npp)*(dlog(dble(npop)) - 
     &           dlog(dble(npop+1)))

            lratio = dmin1(0.d0,lratio)
            alpha = dexp(lratio)
            bern = ggrbinom(1.d0,alpha)

            if(bern .eq. 1) then
               call accept5(nppmax,npopmax,nlocmax,
     &              nalmax,nal,c,ctmp,f,ftmp,drift,drifttmp)
               npop = npop + 1
            endif
         endif

      else
         if(npop .gt. npopmin) then 
c      write(*,*) 'mort'
*     tirage de la pop qui meurt
            ipoprem = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
            
*     tirage de la pop hote
            ipophost = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
            do while(ipophost .eq. ipoprem)
               ipophost = 1 
     &              + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
            enddo
            
*     on range dans la pop d'indice le plus petit
            if(ipophost .gt. ipoprem) then
               ii = ipophost
               ipophost = ipoprem
               ipoprem = ii
            endif
            
*     recherche des cellules qui vont etre reallouees
            call who(c,ipoprem,npp,nppmax,cellpop,ncellpop)
            
*     recherche des cellules de la pop hote
            call who(c,ipophost,npp,nppmax,cellpophost,
     &           ncellpophost)

*     proposition de reallocation dans la pop ipophost
            call merging(ipoprem,ipophost,c,ctmp,nppmax,
     &           ncellpop,cellpop)
            
*     comptage des alleles sur chaque locus pour c puis ctmp
            call countn(nindiv,nlocmax,nlocmax2,npopmax,
     &           nppmax,nal,nalmax,zz,n,indcell,c,ploidy)
            call countn(nindiv,nlocmax,nlocmax2,npopmax,
     &           nppmax,nal,nalmax,zz,ntmp,indcell,ctmp,ploidy)

*     propostion du nouveau tableau de freq et de derives
            call remfreq7bis(ipoprem,ipophost,
     &           npop,npopmax,nloc,nlocmax,nal,
     &           nalmax,f,ftmp,drift,drifttmp,fa,a,ptmp,
     &           ntmp)
            
*     calcul du log du ratio  
*     terme des vraisemblances
            lratio =  llr6(zz,f,ftmp,c,ctmp,indcell,
     &           indcell,npopmax,nlocmax,nalmax,
     &           nindiv,nloc,nlocmax2,nppmax,ploidy)

            lratio = lratio  
     &           + termf9bis(npopmax,nloc,nal,nalmax,n,f,ipophost) 
     &           + termf9bis(npopmax,nloc,nal,nalmax,n,f,ipoprem) 
     &       -termf9bis(npopmax,nloc,nal,nalmax,ntmp,ftmp,ipophost)

* term prior freq
            do iloc = 1,nloc
               lratio = lratio - gglgamfn(dble(nal(iloc)))
            enddo
     
*     terme des proposal sur c
            lratio = lratio + dlog(dble(npop)) - 
     &           dlog(2*dble(ncellpop+ncellpophost+1)) -
     &           lbico(ncellpop+ncellpophost,ncellpop) 

*     terme des priors sur c
            lratio = lratio + 
     &           dble(npp)*(dlog(dble(npop)) - 
     &           dlog(dble(npop-1)))

            lratio = dmin1(0.d0,lratio)
            alpha  = dexp(lratio)
            bern = ggrbinom(1.d0,alpha) 

            if(bern .eq. 1) then
               call accept5(nppmax,npopmax,nlocmax,
     &              nalmax,nal,c,ctmp,f,ftmp,drift,drifttmp)
               npop = npop - 1
            endif
         endif
      endif 
      
      end subroutine bdpop9bis
***********************************************************************
 


************************************************************************
*     ajoute une pop 
*     dans  le tableau des dérives selon le prior
*     et dans le tableau des frequences 
*     selon la conditionelle complète pour les deux nouveaux groupes
*     (sans modifier les tableaux en entrée)
      subroutine addfreq7(npop,npopmax,nloc,nlocmax,
     &     nal,nalmax,isplit,f,ftmp,fa,
     &     drift,drifttmp,a,ptmp,ntmp)
      implicit none
      integer npop,npopmax,nloc,nlocmax,nal(nlocmax),
     &     nalmax,ntmp(npopmax,nloc,nalmax),
     &     isplit
      double precision f(npopmax,nlocmax,nalmax),drift(npopmax),
     &     ftmp(npopmax,nlocmax,nalmax),drifttmp(npopmax),
     &     fa(nlocmax,nalmax),a(nalmax)
      integer iloc,ial,ipop
      double precision ptmp(nalmax),ggrunif
      double precision bern,ggrbinom

c      write(*,*) 'debut addfreq7'
c      write(*,*) 'fa=',fa
c      write(*,*) 'drift=',drift
c      write(*,*) 'ntmp=',ntmp

*     remplissage de f et drift pour le pops pre-existantes
      do ipop = 1,npop
         drifttmp(ipop) = drift(ipop)
         do iloc=1,nloc
            do ial=1,nal(iloc)
               ftmp(ipop,iloc,ial) = f(ipop,iloc,ial)
            enddo 
         enddo
      enddo

*     nouvelles derives
      drifttmp(isplit) = ggrunif(0.d0,1.d0)
      drifttmp(npop+1) = ggrunif(0.d0,1.d0)

*     nouvelles frequences 
*     pour celle qui reste
      do iloc=1,nloc
         do ial = 1,nal(iloc)
            a(ial) = fa(iloc,ial)*(1-drifttmp(isplit))/
     &           drifttmp(isplit)+dble(ntmp(isplit,iloc,ial))
         enddo
         call dirichlet(nal(iloc),nalmax,a,ptmp)
         do ial=1,nal(iloc)
            ftmp(isplit,iloc,ial)  = ptmp(ial)
         enddo
      enddo
*     pour la nouvelle
      do iloc=1,nloc
         do ial = 1,nal(iloc)
            a(ial) = fa(iloc,ial)*(1-drifttmp(npop+1))/
     &           drifttmp(npop+1)+dble(ntmp(npop+1,iloc,ial))
         enddo
         call dirichlet(nal(iloc),nalmax,a,ptmp)
         do ial=1,nal(iloc)
            ftmp(npop+1,iloc,ial)  = ptmp(ial)
         enddo
      enddo

c$$$      write(*,*) 'f(',isplit,1,1,')=',f(isplit,1,1)
c$$$      write(*,*) 'f(',isplit,1,2,')=',f(isplit,1,2)
c$$$      write(*,*) 'ftmp(',isplit,1,1,')=',ftmp(isplit,1,1)
c$$$      write(*,*) 'ftmp(',isplit,1,2,')=',ftmp(isplit,1,2)
c$$$      write(*,*) 'f(',npop+1,1,1,')=',f(npop+1,1,1)
c$$$      write(*,*) 'f(',npop+1,1,2,')=',f(npop+1,1,2)
c$$$      write(*,*) 'ftmp(',npop+1,1,1,')=',ftmp(npop+1,1,1)
c$$$      write(*,*) 'ftmp(',npop+1,1,2,')=',ftmp(npop+1,1,2)
c$$$
c$$$      write(*,*) 'f(',isplit,2,1,')=',f(isplit,2,1)
c$$$      write(*,*) 'f(',isplit,2,2,')=',f(isplit,2,2)
c$$$      write(*,*) 'ftmp(',isplit,2,1,')=',ftmp(isplit,2,1)
c$$$      write(*,*) 'ftmp(',isplit,2,2,')=',ftmp(isplit,2,2)
c$$$      write(*,*) 'f(',npop+1,2,1,')=',f(npop+1,2,1)
c$$$      write(*,*) 'f(',npop+1,2,2,')=',f(npop+1,2,2)
c$$$      write(*,*) 'ftmp(',npop+1,2,1,')=',ftmp(npop+1,2,1)
c$$$      write(*,*) 'ftmp(',npop+1,2,2,')=',ftmp(npop+1,2,2)
c$$$
      end subroutine addfreq7
***********************************************************************



***********************************************************************
*     ajoute une pop 
*     dans le tableau des frequences 
*     selon la conditionelle complète pour les deux nouveaux groupes
*     (sans modifier les tableaux en entrée)
      subroutine addfreq8(npop,npopmax,nloc,nlocmax,
     &     nal,nalmax,isplit,f,ftmp,fa,
     &     drift,drifttmp,a,ptmp,ntmp)
      implicit none
      integer npop,npopmax,nloc,nlocmax,nal(nlocmax),
     &     nalmax,ntmp(npopmax,nloc,nalmax),
     &     isplit
      double precision f(npopmax,nlocmax,nalmax),drift(npopmax),
     &     ftmp(npopmax,nlocmax,nalmax),drifttmp(npopmax),
     &     fa(nlocmax,nalmax),a(nalmax)
      integer iloc,ial,ipop
      double precision ptmp(nalmax),ggrunif
*     remplissage de f et drift pour le pops pre-existantes
      do ipop = 1,npop
         do iloc=1,nloc
            do ial=1,nal(iloc)
               ftmp(ipop,iloc,ial) = f(ipop,iloc,ial)
            enddo 
         enddo
      enddo
*     nouvelles frequences 
*     pour celle qui reste
      do iloc=1,nloc
         do ial = 1,nal(iloc)
            a(ial) = fa(iloc,ial)*(1-drifttmp(isplit))/
     &           drifttmp(isplit)+dble(ntmp(isplit,iloc,ial))
         enddo
         call dirichlet(nal(iloc),nalmax,a,ptmp)
         do ial=1,nal(iloc)
            ftmp(isplit,iloc,ial)  = ptmp(ial)
         enddo
      enddo
*     pour la nouvelle
      do iloc=1,nloc
         do ial = 1,nal(iloc)
            a(ial) = fa(iloc,ial)*(1-drifttmp(npop+1))/
     &           drifttmp(npop+1)+dble(ntmp(npop+1,iloc,ial))
         enddo
         call dirichlet(nal(iloc),nalmax,a,ptmp)
         do ial=1,nal(iloc)
            ftmp(npop+1,iloc,ial)  = ptmp(ial)
         enddo
      enddo
      end subroutine addfreq8
***********************************************************************



***********************************************************************
*     ajoute une pop 
*     dans le tableau des frequences 
*     selon la conditionelle complète pour les deux nouveaux groupes
*     (sans modifier les tableaux en entrée)
      subroutine addfall(npop,npopmax,nlocd,nlocd2,nloch,nql,
     &     ncolt,nal,nalmax,isplit,f,ftmp,fa,drift,drifttmp,a,ptmp,ntmp,
     &     usegeno2,usegeno1,useql)
      implicit none
      integer npop,npopmax,nlocd,nlocd2,nloch,nql,ncolt,nal(ncolt),
     &     nalmax,ntmp(npopmax,ncolt,nalmax),
     &     isplit,usegeno2,usegeno1,useql
      double precision f(npopmax,ncolt,nalmax),drift(npopmax),
     &     ftmp(npopmax,ncolt,nalmax),drifttmp(npopmax),
     &     fa(ncolt,nalmax),a(nalmax)
      integer iloc,ial,ipop
      double precision ptmp(nalmax),ggrunif
*     remplissage de f et drift pour le pops pre-existantes
      do ipop = 1,npop
         do iloc=1,ncolt
            do ial=1,nal(iloc)
               ftmp(ipop,iloc,ial) = f(ipop,iloc,ial)
            enddo 
         enddo
      enddo
*     nouvelles frequences 
      if(usegeno2 .eq. 1) then 
*     pour celle qui reste
         do iloc=1,nlocd
            do ial = 1,nal(iloc)
               a(ial) = fa(iloc,ial)*(1-drifttmp(isplit))/
     &              drifttmp(isplit)+dble(ntmp(isplit,iloc,ial))
            enddo
            call dirichlet(nal(iloc),nalmax,a,ptmp)
            do ial=1,nal(iloc)
               ftmp(isplit,iloc,ial)  = ptmp(ial)
            enddo
         enddo
*     pour la nouvelle
         do iloc=1,nlocd
            do ial = 1,nal(iloc)
               a(ial) = fa(iloc,ial)*(1-drifttmp(npop+1))/
     &              drifttmp(npop+1)+dble(ntmp(npop+1,iloc,ial))
            enddo
            call dirichlet(nal(iloc),nalmax,a,ptmp)
            do ial=1,nal(iloc)
               ftmp(npop+1,iloc,ial)  = ptmp(ial)
            enddo
         enddo
      endif
      if(usegeno1 .eq. 1) then 
*     pour celle qui reste
         do iloc=1,nloch
            do ial = 1,nal(nlocd+iloc)
               a(ial) = fa(nlocd+iloc,ial)*(1-drifttmp(isplit))/
     &              drifttmp(isplit)+dble(ntmp(isplit,nlocd+iloc,ial))
            enddo
            call dirichlet(nal(nlocd+iloc),nalmax,a,ptmp)
            do ial=1,nal(nlocd+iloc)
               ftmp(isplit,nlocd+iloc,ial)  = ptmp(ial)
            enddo
         enddo
*     pour la nouvelle
         do iloc=1,nloch
            do ial = 1,nal(nlocd+iloc)
               a(ial) = fa(nlocd+iloc,ial)*(1-drifttmp(npop+1))/
     &              drifttmp(npop+1)+dble(ntmp(npop+1,nlocd+iloc,ial))
            enddo
            call dirichlet(nal(nlocd+iloc),nalmax,a,ptmp)
            do ial=1,nal(nlocd+iloc)
               ftmp(npop+1,nlocd+iloc,ial)  = ptmp(ial)
            enddo
         enddo
      endif
      if(useql .eq. 1) then 
*     pour celle qui reste
         do iloc=1,nql
            do ial = 1,nal(nlocd+nloch+iloc)
               a(ial) = fa(nlocd+nloch+iloc,ial)*(1-drifttmp(isplit))/
     &              drifttmp(isplit)+
     &              dble(ntmp(isplit,nlocd+nloch+iloc,ial))
            enddo
            call dirichlet(nal(nlocd+nloch+iloc),nalmax,a,ptmp)
            do ial=1,nal(nlocd+nloch+iloc)
               ftmp(isplit,nlocd+nloch+iloc,ial)  = ptmp(ial)
            enddo
         enddo
*     pour la nouvelle
         do iloc=1,nql
            do ial = 1,nal(nlocd+nloch+iloc)
               a(ial) = fa(nlocd+nloch+iloc,ial)*(1-drifttmp(npop+1))/
     &              drifttmp(npop+1)+
     &              dble(ntmp(npop+1,nlocd+nloch+iloc,ial))
            enddo
            call dirichlet(nal(nlocd+nloch+iloc),nalmax,a,ptmp)
            do ial=1,nal(nlocd+nloch+iloc)
               ftmp(npop+1,nlocd+nloch+iloc,ial)  = ptmp(ial)
            enddo
         enddo
      endif
      end subroutine addfall
***********************************************************************





***********************************************************************
*     ajoute une pop 
*     dans  le tableau des frequences  selon le prior
*     selon la conditionelle complète pour les deux nouveaux groupes
*     et une valeur 0.5d0 dans le tableau des dérives 
*     (sans modifier les tableaux en entrée)
*     pour court-circuiter le F-model
      subroutine addfreq7bis(npop,npopmax,nloc,nlocmax,
     &     nal,nalmax,isplit,f,ftmp,
     &     fa,drift,drifttmp,a,ptmp,ntmp)
      implicit none
      integer npop,npopmax,nloc,nlocmax,nal(nlocmax),
     &     nalmax,ntmp(npopmax,nloc,nalmax),
     &     isplit
      double precision f(npopmax,nlocmax,nalmax),drift(npopmax),
     &     ftmp(npopmax,nlocmax,nalmax),drifttmp(npopmax),
     &     fa(nlocmax,nalmax),a(nalmax)
      integer iloc,ial,ipop
      double precision ptmp(nalmax)
*     remplissage de f et drift pour le pops pre-existantes
      do ipop = 1,npop
         drifttmp(ipop) = drift(ipop)
         do iloc=1,nloc
            do ial=1,nal(iloc)
               ftmp(ipop,iloc,ial) = f(ipop,iloc,ial)
            enddo 
         enddo
      enddo
*     nouvelles derives
      drifttmp(isplit) = 0.5d0
      drifttmp(npop+1) = 0.5d0
*     nouvelles frequences 
*     pour celle qui reste
      do iloc=1,nloc
         do ial = 1,nal(iloc)
            a(ial) = fa(iloc,ial)*(1-drifttmp(isplit))/
     &           drifttmp(isplit)+dble(ntmp(isplit,iloc,ial))
         enddo
         call dirichlet(nal(iloc),nalmax,a,ptmp)
         do ial=1,nal(iloc)
            ftmp(isplit,iloc,ial)  = ptmp(ial)
         enddo
      enddo
*     pour la nouvelle
      do iloc=1,nloc
         do ial = 1,nal(iloc)
            a(ial) = fa(iloc,ial)*(1-drifttmp(npop+1))/
     &           drifttmp(npop+1)+dble(ntmp(npop+1,iloc,ial))
         enddo
         call dirichlet(nal(iloc),nalmax,a,ptmp)
         do ial=1,nal(iloc)
            ftmp(npop+1,iloc,ial)  = ptmp(ial)
         enddo
      enddo
      end subroutine addfreq7bis
***********************************************************************




***********************************************************************
*     ajoute une pop 
*     dans  le tableau des frequences selon le prior
*     selon la conditionelle complète pour les deux nouveaux groupes
*     et une valeur 0.5d0 dans le tableau des dérives 
*     (sans modifier les tableaux en entrée)
*     pour court-circuiter le F-model
      subroutine addfallbis(npop,npopmax,nlocd,nlocd2,nloch,nql,ncolt,
     &     nal,nalmax,isplit,f,ftmp,fa,drift,drifttmp,a,ptmp,ntmp,
     &     usegeno2,usegeno1,useql)
      implicit none
      integer npop,npopmax,nlocd,nlocd2,nloch,nql,ncolt,nal(ncolt),
     &     nalmax,ntmp(npopmax,ncolt,nalmax),isplit,
     &     usegeno2,usegeno1,useql
      double precision f(npopmax,ncolt,nalmax),drift(npopmax),
     &     ftmp(npopmax,ncolt,nalmax),drifttmp(npopmax),
     &     fa(ncolt,nalmax),a(nalmax)
      integer iloc,ial,ipop
      double precision ptmp(nalmax)
*     remplissage de f et drift pour le pops pre-existantes
      do ipop = 1,npop
         drifttmp(ipop) = drift(ipop)
         do iloc=1,ncolt
            do ial=1,nal(iloc)
               ftmp(ipop,iloc,ial) = f(ipop,iloc,ial)
            enddo 
         enddo
      enddo
*     nouvelles derives
      drifttmp(isplit) = 0.5d0
      drifttmp(npop+1) = 0.5d0
*     nouvelles frequences 
      if(usegeno2 .eq. 1) then
*     pour celle qui reste
         do iloc=1,nlocd
            do ial = 1,nal(iloc)
               a(ial) = fa(iloc,ial)*(1-drifttmp(isplit))/
     &              drifttmp(isplit)+dble(ntmp(isplit,iloc,ial))
            enddo
            call dirichlet(nal(iloc),nalmax,a,ptmp)
            do ial=1,nal(iloc)
               ftmp(isplit,iloc,ial)  = ptmp(ial)
            enddo
         enddo
*     pour la nouvelle
         do iloc=1,nlocd
            do ial = 1,nal(iloc)
               a(ial) = fa(iloc,ial)*(1-drifttmp(npop+1))/
     &              drifttmp(npop+1)+dble(ntmp(npop+1,iloc,ial))
            enddo
            call dirichlet(nal(iloc),nalmax,a,ptmp)
            do ial=1,nal(iloc)
               ftmp(npop+1,iloc,ial)  = ptmp(ial)
            enddo
         enddo
      endif
      if(usegeno1 .eq. 1) then
*     pour celle qui reste
         do iloc=1,nloch
            do ial = 1,nal(nlocd+iloc)
               a(ial) = fa(nlocd+iloc,ial)*(1-drifttmp(isplit))/
     &              drifttmp(isplit)+dble(ntmp(isplit,nlocd+iloc,ial))
            enddo
            call dirichlet(nal(nlocd+iloc),nalmax,a,ptmp)
            do ial=1,nal(nlocd+iloc)
               ftmp(isplit,nlocd+iloc,ial)  = ptmp(ial)
            enddo
         enddo
*     pour la nouvelle
         do iloc=1,nloch
            do ial = 1,nal(nlocd+iloc)
               a(ial) = fa(nlocd+iloc,ial)*(1-drifttmp(npop+1))/
     &              drifttmp(npop+1)+dble(ntmp(npop+1,nlocd+iloc,ial))
            enddo
            call dirichlet(nal(nlocd+iloc),nalmax,a,ptmp)
            do ial=1,nal(nlocd+iloc)
               ftmp(npop+1,nlocd+iloc,ial)  = ptmp(ial)
            enddo
         enddo
      endif
      if(useql .eq. 1) then
*     pour celle qui reste
         do iloc=1,nql
            do ial = 1,nal(nlocd+nloch+iloc)
               a(ial) = fa(nlocd+nloch+iloc,ial)*(1-drifttmp(isplit))/
     &              drifttmp(isplit)+
     &              dble(ntmp(isplit,nlocd+nloch+iloc,ial))
            enddo
            call dirichlet(nal(nlocd+nloch+iloc),nalmax,a,ptmp)
            do ial=1,nal(nlocd+nloch+iloc)
               ftmp(isplit,nlocd+nloch+iloc,ial)  = ptmp(ial)
            enddo
         enddo
*     pour la nouvelle
         do iloc=1,nql
            do ial = 1,nal(nlocd+nloch+iloc)
               a(ial) = fa(nlocd+nloch+iloc,ial)*(1-drifttmp(npop+1))/
     &              drifttmp(npop+1)+
     &              dble(ntmp(npop+1,nlocd+nloch+iloc,ial))
            enddo
            call dirichlet(nal(nlocd+nloch+iloc),nalmax,a,ptmp)
            do ial=1,nal(nlocd+nloch+iloc)
               ftmp(npop+1,nlocd+nloch+iloc,ial)  = ptmp(ial)
            enddo
         enddo
      endif
      end subroutine addfallbis
***********************************************************************



***********************************************************************
*     enleve une pop des tableau des frequences et des derives
*     tirage d'une freq selon posterior apres un merge de deux pops
*     tirage d'une derive selon prior 
*     sans modifier des tableaux en entrée
      subroutine remfreq7(ipoprem,ipophost,
     &     npop,npopmax,nloc,nlocmax,nal,
     &     nalmax,f,ftmp,drift,drifttmp,fa,a,ptmp,ntmp)
      implicit none
      integer ipoprem,ipophost,
     &     npop,npopmax,nloc,nlocmax,nlocmax2,nal(nlocmax),
     &     nalmax,nppmax,nindiv,
     &     ntmp(npopmax,nlocmax,nalmax)
      double precision f(npopmax,nlocmax,nalmax),drift(npopmax),
     &     ftmp(npopmax,nlocmax,nalmax),drifttmp(npopmax),
     &     fa(nlocmax,nalmax),a(nalmax) 
      integer ipop,iloc,ial
      double precision ptmp(nalmax),ggrunif

*     supprime la pop ipoprem
      if(ipoprem .eq. 1) then
         do ipop =2,npop
            do iloc=1,nloc
               do ial=1,nal(iloc)
                  ftmp(ipop-1,iloc,ial) = f(ipop,iloc,ial)
               enddo 
            enddo
            drifttmp(ipop-1) = drift(ipop)
         enddo
      else
         do ipop =1,ipoprem-1
            do iloc=1,nloc
               do ial=1,nal(iloc)
                  ftmp(ipop,iloc,ial) = f(ipop,iloc,ial)
               enddo 
            enddo
            drifttmp(ipop) = drift(ipop)
         enddo
         if(ipoprem .ne. npop) then
            do ipop =ipoprem+1,npop
               do iloc=1,nloc
                  do ial=1,nal(iloc)
                     ftmp(ipop-1,iloc,ial) = f(ipop,iloc,ial)
                  enddo 
               enddo
               drifttmp(ipop-1) = drift(ipop)
            enddo
         endif
      endif
      do ipop=npop,npopmax
         do iloc=1,nloc
            do ial=1,nal(iloc)
               ftmp(ipop,iloc,ial) = -999
            enddo 
         enddo
         drifttmp(ipop) = -999
      enddo

*     tirage de la derive et de la freq
      
*     derive pour la nouvelle pop
      drifttmp(ipophost) = ggrunif(0.d0,1.d0)
      
*     frequences pour la nouvelle pop
      do iloc=1,nloc
         do ial = 1,nal(iloc)
            a(ial) = fa(iloc,ial)*(1-drifttmp(ipophost))/
     &           drifttmp(ipophost)+
     &           dble(ntmp(ipophost,iloc,ial))
         enddo
         call dirichlet(nal(iloc),nalmax,a,ptmp)
         do ial=1,nal(iloc)
            ftmp(ipophost,iloc,ial)  = ptmp(ial)
         enddo
      enddo

c$$$      write(*,*) 'f(',ipophost,1,1,')=',f(ipophost,1,1)
c$$$      write(*,*) 'f(',ipophost,1,2,')=',f(ipophost,1,2)
c$$$      write(*,*) 'f(',ipophost,2,1,')=',f(ipophost,2,1)
c$$$      write(*,*) 'f(',ipophost,2,2,')=',f(ipophost,2,2)
c$$$
c$$$      write(*,*) 'ftmp(',ipophost,1,1,')=',ftmp(ipophost,1,1)
c$$$      write(*,*) 'ftmp(',ipophost,1,2,')=',ftmp(ipophost,1,2)
c$$$      write(*,*) 'ftmp(',ipophost,2,1,')=',ftmp(ipophost,2,1)
c$$$      write(*,*) 'ftmp(',ipophost,2,2,')=',ftmp(ipophost,2,2)

      end subroutine remfreq7
******************************************************************


******************************************************************
*     enleve une pop des tableau des frequences
*     tirage d'une freq selon posterior apres un merge de deux pops
*     sans modifier des tableaux en entrée
      subroutine remfreq8(ipoprem,ipophost,
     &     npop,npopmax,nloc,nlocmax,nal,
     &     nalmax,f,ftmp,drift,drifttmp,fa,a,ptmp,ntmp)
      implicit none
      integer ipoprem,ipophost,
     &     npop,npopmax,nloc,nlocmax,nlocmax2,nal(nlocmax),
     &     nalmax,nppmax,nindiv,
     &     ntmp(npopmax,nlocmax,nalmax)
      double precision f(npopmax,nlocmax,nalmax),drift(npopmax),
     &     ftmp(npopmax,nlocmax,nalmax),drifttmp(npopmax),
     &     fa(nlocmax,nalmax),a(nalmax) 
      integer ipop,iloc,ial
      double precision ptmp(nalmax),ggrunif

*     supprime la pop ipoprem
      if(ipoprem .eq. 1) then
         do ipop =2,npop
            do iloc=1,nloc
               do ial=1,nal(iloc)
                  ftmp(ipop-1,iloc,ial) = f(ipop,iloc,ial)
               enddo 
            enddo
         enddo
      else
         do ipop =1,ipoprem-1
            do iloc=1,nloc
               do ial=1,nal(iloc)
                  ftmp(ipop,iloc,ial) = f(ipop,iloc,ial)
               enddo 
            enddo
         enddo
         if(ipoprem .ne. npop) then
            do ipop =ipoprem+1,npop
               do iloc=1,nloc
                  do ial=1,nal(iloc)
                     ftmp(ipop-1,iloc,ial) = f(ipop,iloc,ial)
                  enddo 
               enddo
            enddo
         endif
      endif
      do ipop=npop,npopmax
         do iloc=1,nloc
            do ial=1,nal(iloc)
               ftmp(ipop,iloc,ial) = -999
            enddo 
         enddo
      enddo

*     frequences pour la nouvelle pop
      do iloc=1,nloc
         do ial = 1,nal(iloc)
            a(ial) = fa(iloc,ial)*(1-drifttmp(ipophost))/
     &           drifttmp(ipophost)+
     &           dble(ntmp(ipophost,iloc,ial))
         enddo
         call dirichlet(nal(iloc),nalmax,a,ptmp)
         do ial=1,nal(iloc)
            ftmp(ipophost,iloc,ial)  = ptmp(ial)
         enddo
      enddo
      end subroutine remfreq8
******************************************************************



******************************************************************
*     enleve une pop des tableau des frequences
*     tirage d'une freq selon posterior apres un merge de deux pops
*     sans modifier des tableaux en entrée
      subroutine remfall(ipoprem,ipophost,
     &     npop,npopmax,nlocd,nloch,nql,ncolt,nal,
     &     nalmax,f,ftmp,drift,drifttmp,fa,a,ptmp,ntmp,
     &     usegeno2,usegeno1,useql)
      implicit none
      integer ipoprem,ipophost,
     &     npop,npopmax,nlocd,nloch,nql,ncolt,nal(ncolt),
     &     nalmax,nppmax,nindiv,ntmp(npopmax,ncolt,nalmax),
     &      usegeno2,usegeno1,useql
      double precision f(npopmax,ncolt,nalmax),drift(npopmax),
     &     ftmp(npopmax,ncolt,nalmax),drifttmp(npopmax),
     &     fa(ncolt,nalmax),a(nalmax) 
      integer ipop,iloc,ial
      double precision ptmp(nalmax),ggrunif
*     supprime la pop ipoprem
      if(ipoprem .eq. 1) then
         do ipop =2,npop
            do iloc=1,ncolt
               do ial=1,nal(iloc)
                  ftmp(ipop-1,iloc,ial) = f(ipop,iloc,ial)
               enddo 
            enddo
         enddo
      else
         do ipop =1,ipoprem-1
            do iloc=1,ncolt
               do ial=1,nal(iloc)
                  ftmp(ipop,iloc,ial) = f(ipop,iloc,ial)
               enddo 
            enddo
         enddo
         if(ipoprem .ne. npop) then
            do ipop =ipoprem+1,npop
               do iloc=1,ncolt
                  do ial=1,nal(iloc)
                     ftmp(ipop-1,iloc,ial) = f(ipop,iloc,ial)
                  enddo 
               enddo
            enddo
         endif
      endif
      do ipop=npop,npopmax
         do iloc=1,ncolt
            do ial=1,nal(iloc)
               ftmp(ipop,iloc,ial) = -999
            enddo 
         enddo
      enddo
*     frequences pour la nouvelle pop
      if(usegeno2 .eq. 1) then 
         do iloc=1,nlocd
            do ial = 1,nal(iloc)
               a(ial) = fa(iloc,ial)*(1-drifttmp(ipophost))/
     &              drifttmp(ipophost)+
     &              dble(ntmp(ipophost,iloc,ial))
            enddo
            call dirichlet(nal(iloc),nalmax,a,ptmp)
            do ial=1,nal(iloc)
               ftmp(ipophost,iloc,ial)  = ptmp(ial)
            enddo
         enddo
      endif
      if(usegeno1 .eq. 1) then 
         do iloc=1,nloch
            do ial = 1,nal(nlocd+iloc)
               a(ial) = fa(nlocd+iloc,ial)*(1-drifttmp(ipophost))/
     &              drifttmp(ipophost)+
     &              dble(ntmp(ipophost,nlocd+iloc,ial))
            enddo
            call dirichlet(nal(nlocd+iloc),nalmax,a,ptmp)
            do ial=1,nal(nlocd+iloc)
               ftmp(ipophost,nlocd+iloc,ial)  = ptmp(ial)
            enddo
         enddo
      endif
      if(useql .eq. 1) then 
         do iloc=1,nql
            do ial = 1,nal(nlocd+nloch+iloc)
               a(ial) = fa(nlocd+nloch+iloc,ial)*(1-drifttmp(ipophost))/
     &              drifttmp(ipophost)+
     &              dble(ntmp(ipophost,nlocd+nloch+iloc,ial))
            enddo
            call dirichlet(nal(nlocd+nloch+iloc),nalmax,a,ptmp)
            do ial=1,nal(nlocd+nloch+iloc)
               ftmp(ipophost,nlocd+nloch+iloc,ial)  = ptmp(ial)
            enddo
         enddo
      endif
      end subroutine remfall
******************************************************************




******************************************************************
*     enleve une pop des tableau des derives
*     sans modifier des tableaux en entrée
      subroutine remdrift(ipoprem,ipophost,npop,npopmax,drift,drifttmp,
     &     shape1,shape2)
      implicit none
      integer ipoprem,ipophost,
     &     npop,npopmax
      double precision drift(npopmax),drifttmp(npopmax),shape1,shape2
      integer ipop
      double precision ggrbet
*     supprime la pop ipoprem
      if(ipoprem .eq. 1) then
         do ipop =2,npop
            drifttmp(ipop-1) = drift(ipop)
         enddo
      else
         do ipop =1,ipoprem-1
            drifttmp(ipop) = drift(ipop)
         enddo
         if(ipoprem .ne. npop) then
            do ipop =ipoprem+1,npop
               drifttmp(ipop-1) = drift(ipop)
            enddo
         endif
      endif
      do ipop=npop,npopmax
         drifttmp(ipop) = -999
      enddo
*     derive pour la nouvelle pop
      drifttmp(ipophost) =  ggrbet(shape1,shape2)
      end subroutine remdrift


***********************************************
*     enleve une pop des tableau des derives
*     sans modifier des tableaux en entrée
*     d* = (d1+d2)/2
      subroutine remdrift2(ipoprem,ipophost,npop,npopmax,drift,drifttmp)
      implicit none
      integer ipoprem,ipophost,
     &     npop,npopmax
      double precision drift(npopmax),drifttmp(npopmax)
      integer ipop
      double precision ggrbet
*     supprime la pop ipoprem
      if(ipoprem .eq. 1) then
         do ipop =2,npop
            drifttmp(ipop-1) = drift(ipop)
         enddo
      else
         do ipop =1,ipoprem-1
            drifttmp(ipop) = drift(ipop)
         enddo
         if(ipoprem .ne. npop) then
            do ipop =ipoprem+1,npop
               drifttmp(ipop-1) = drift(ipop)
            enddo
         endif
      endif
      do ipop=npop,npopmax
         drifttmp(ipop) = -999
      enddo
*     derive pour la nouvelle pop
      drifttmp(ipophost) = (drift(ipophost) + drift(ipoprem))/2
      end subroutine remdrift2



*****************************************************************
*     enleve une pop des tableaux des frequences et des derives
*     tirage d'une freq selon posterior apres un merge de deux pops
*     la nouvelle derive est mise à 0.5d0
*     (sans modifier les tableaux en entrée)
*     c'est pour court-circuiter le F-model 
      subroutine remfreq7bis(ipoprem,ipophost,
     &     npop,npopmax,nloc,nlocmax,nal,
     &     nalmax,f,ftmp,drift,drifttmp,fa,a,ptmp,ntmp)
      implicit none
      integer ipoprem,ipophost,
     &     npop,npopmax,nloc,nlocmax,nal(nlocmax),
     &     nalmax,nppmax,nindiv,
     &     ntmp(npopmax,nlocmax,nalmax)
      double precision f(npopmax,nlocmax,nalmax),drift(npopmax),
     &     ftmp(npopmax,nlocmax,nalmax),drifttmp(npopmax),
     &     fa(nlocmax,nalmax),a(nalmax) 
      integer ipop,iloc,ial
      double precision ptmp(nalmax)

*     supprime la pop ipoprem
      if(ipoprem .eq. 1) then
         do ipop =2,npop
            do iloc=1,nloc
               do ial=1,nal(iloc)
                  ftmp(ipop-1,iloc,ial) = f(ipop,iloc,ial)
               enddo 
            enddo
            drifttmp(ipop-1) = drift(ipop)
         enddo
      else
         do ipop =1,ipoprem-1
            do iloc=1,nloc
               do ial=1,nal(iloc)
                  ftmp(ipop,iloc,ial) = f(ipop,iloc,ial)
               enddo 
            enddo
            drifttmp(ipop) = drift(ipop)
         enddo
         if(ipoprem .ne. npop) then
            do ipop =ipoprem+1,npop
               do iloc=1,nloc
                  do ial=1,nal(iloc)
                     ftmp(ipop-1,iloc,ial) = f(ipop,iloc,ial)
                  enddo 
               enddo
               drifttmp(ipop-1) = drift(ipop)
            enddo
         endif
      endif
      do ipop=npop,npopmax
         do iloc=1,nloc
            do ial=1,nal(iloc)
               ftmp(ipop,iloc,ial) = -999
            enddo 
         enddo
         drifttmp(ipop) = -999
      enddo

*     tirage de la derive et de la freq
      
*     derive pour la nouvelle pop
      drifttmp(ipophost) = 0.5d0
      
*     frequences pour la nouvelle pop
      do iloc=1,nloc
         do ial = 1,nal(iloc)
            a(ial) = fa(iloc,ial)*(1-drifttmp(ipophost))/
     &           drifttmp(ipophost)+
     &           dble(ntmp(ipophost,iloc,ial))
         enddo
         call dirichlet(nal(iloc),nalmax,a,ptmp)
         do ial=1,nal(iloc)
            ftmp(ipophost,iloc,ial)  = ptmp(ial)
         enddo
      enddo

c$$$      write(*,*) 'f(',ipophost,1,1,')=',f(ipophost,1,1)
c$$$      write(*,*) 'f(',ipophost,1,2,')=',f(ipophost,1,2)
c$$$      write(*,*) 'f(',ipophost,2,1,')=',f(ipophost,2,1)
c$$$      write(*,*) 'f(',ipophost,2,2,')=',f(ipophost,2,2)
c$$$
c$$$      write(*,*) 'ftmp(',ipophost,1,1,')=',ftmp(ipophost,1,1)
c$$$      write(*,*) 'ftmp(',ipophost,1,2,')=',ftmp(ipophost,1,2)
c$$$      write(*,*) 'ftmp(',ipophost,2,1,')=',ftmp(ipophost,2,1)
c$$$      write(*,*) 'ftmp(',ipophost,2,2,')=',ftmp(ipophost,2,2)

      end subroutine remfreq7bis
************************************************************************




*****************************************************************
*     enleve une pop des tableaux des frequences et des derives
*     tirage d'une freq selon posterior apres un merge de deux pops
*     la nouvelle derive est mise à 0.5d0
*     (sans modifier les tableaux en entrée)
*     c'est pour court-circuiter le F-model 
      subroutine remfallbis(ipoprem,ipophost,
     &     npop,npopmax,nlocd,nlocd2,nloch,nql,ncolt,nal,
     &     nalmax,f,ftmp,drift,drifttmp,fa,a,ptmp,ntmp,
     &     usegeno2,usegeno1,useql)
      implicit none
      integer ipoprem,ipophost,
     &     npop,npopmax,nlocd,nlocd2,nloch,nql,ncolt,nal(ncolt),
     &     nalmax,nppmax,nindiv, ntmp(npopmax,ncolt,nalmax),
     &    usegeno2,usegeno1,useql
      double precision f(npopmax,ncolt,nalmax),drift(npopmax),
     &     ftmp(npopmax,ncolt,nalmax),drifttmp(npopmax),
     &     fa(ncolt,nalmax),a(nalmax) 
      integer ipop,iloc,ial
      double precision ptmp(nalmax)

*     supprime la pop ipoprem
      if(ipoprem .eq. 1) then
         do ipop =2,npop
            do iloc=1,ncolt
               do ial=1,nal(iloc)
                  ftmp(ipop-1,iloc,ial) = f(ipop,iloc,ial)
               enddo 
            enddo
            drifttmp(ipop-1) = drift(ipop)
         enddo
      else
         do ipop =1,ipoprem-1
            do iloc=1,ncolt
               do ial=1,nal(iloc)
                  ftmp(ipop,iloc,ial) = f(ipop,iloc,ial)
               enddo 
            enddo
            drifttmp(ipop) = drift(ipop)
         enddo
         if(ipoprem .ne. npop) then
            do ipop =ipoprem+1,npop
               do iloc=1,ncolt
                  do ial=1,nal(iloc)
                     ftmp(ipop-1,iloc,ial) = f(ipop,iloc,ial)
                  enddo 
               enddo
               drifttmp(ipop-1) = drift(ipop)
            enddo
         endif
      endif
      do ipop=npop,npopmax
         do iloc=1,ncolt
            do ial=1,nal(iloc)
               ftmp(ipop,iloc,ial) = -999
            enddo 
         enddo
         drifttmp(ipop) = -999
      enddo

*     tirage de la derive et de la freq
      
*     derive pour la nouvelle pop
      drifttmp(ipophost) = 0.5d0
      
*     frequences pour la nouvelle pop
      if(usegeno2 .eq. 1) then 
         do iloc=1,nlocd
            do ial = 1,nal(iloc)
               a(ial) = fa(iloc,ial)*(1-drifttmp(ipophost))/
     &              drifttmp(ipophost)+
     &              dble(ntmp(ipophost,iloc,ial))
            enddo
            call dirichlet(nal(iloc),nalmax,a,ptmp)
            do ial=1,nal(iloc)
               ftmp(ipophost,iloc,ial)  = ptmp(ial)
            enddo
         enddo
      endif
      if(usegeno1 .eq. 1) then 
         do iloc=1,nloch
            do ial = 1,nal(nlocd+iloc)
               a(ial) = fa(nlocd+iloc,ial)*(1-drifttmp(ipophost))/
     &              drifttmp(ipophost)+
     &              dble(ntmp(ipophost,nlocd+iloc,ial))
            enddo
            call dirichlet(nal(nlocd+iloc),nalmax,a,ptmp)
            do ial=1,nal(nlocd+iloc)
               ftmp(ipophost,nlocd+iloc,ial)  = ptmp(ial)
            enddo
         enddo
      endif
      if(useql .eq. 1) then 
         do iloc=1,nql
            do ial = 1,nal(nlocd+nloch+iloc)
               a(ial) = fa(nlocd+nloch+iloc,ial)*(1-drifttmp(ipophost))/
     &              drifttmp(ipophost)+
     &              dble(ntmp(ipophost,nlocd+nloch+iloc,ial))
            enddo
            call dirichlet(nal(nlocd+nloch+iloc),nalmax,a,ptmp)
            do ial=1,nal(nlocd+nloch+iloc)
               ftmp(ipophost,nlocd+nloch+iloc,ial)  = ptmp(ial)
            enddo
         enddo
      endif
c$$$      write(*,*) 'f(',ipophost,1,1,')=',f(ipophost,1,1)
c$$$      write(*,*) 'f(',ipophost,1,2,')=',f(ipophost,1,2)
c$$$      write(*,*) 'f(',ipophost,2,1,')=',f(ipophost,2,1)
c$$$      write(*,*) 'f(',ipophost,2,2,')=',f(ipophost,2,2)
c$$$
c$$$      write(*,*) 'ftmp(',ipophost,1,1,')=',ftmp(ipophost,1,1)
c$$$      write(*,*) 'ftmp(',ipophost,1,2,')=',ftmp(ipophost,1,2)
c$$$      write(*,*) 'ftmp(',ipophost,2,1,')=',ftmp(ipophost,2,1)
c$$$      write(*,*) 'ftmp(',ipophost,2,2,')=',ftmp(ipophost,2,2)

      end subroutine remfallbis
************************************************************************




************************************************************************
*     terme des freq dans le log ratio pour un split 
      double precision function termfsplit(isplit,npop,npopmax,
     &     nlocmax,nal,nalmax,
     &     f,ftmp,n,ntmp,fa,drift,drifttmp)
      implicit none
      integer npopmax,nlocmax,nalmax,
     &     n(npopmax,nlocmax,nalmax),
     &     ntmp(npopmax,nlocmax,nalmax),nal(nlocmax),
     &     isplit,npop
      double precision f(npopmax,nlocmax,nalmax),drift(npopmax),
     &     ftmp(npopmax,nlocmax,nalmax),drifttmp(npopmax),
     &     fa(nlocmax,nalmax)
      integer iloc, ial,nn
      double precision q,gglgamfn,ss,tt

      
      termfsplit  = 0.
*    ln( pi[ f_isplit |...] / pi[f_isplit| fa, drift] ) 
      tt = 0
      q = (1-drift(isplit))/drift(isplit)
      do iloc = 1,nlocmax
         ss = 0
         nn = 0
         do ial = 1,nal(iloc)
            ss = ss + gglgamfn(fa(iloc,ial) * q) -
     &           gglgamfn(fa(iloc,ial) * q +
     &           dble(n(isplit,iloc,ial))) +
     &           dble(n(isplit,iloc,ial))*
     &           dlog(f(isplit,iloc,ial))
            nn = nn + n(isplit,iloc,ial)
         enddo
         tt = gglgamfn(q + dble(nn)) - gglgamfn(q) + ss
      enddo
      termfsplit = termfsplit + tt

*    ln( pi[ f_isplit^* |fa, drift] / pi[f_isplit^*| ...] ) 
      tt = 0
      q = (1-drifttmp(isplit))/drifttmp(isplit)
      do iloc = 1,nlocmax
         ss = 0
         nn = 0
         do ial = 1,nal(iloc)
            ss = ss + gglgamfn(fa(iloc,ial) * q) -
     &           gglgamfn(fa(iloc,ial) * q +
     &           dble(ntmp(isplit,iloc,ial))) +
     &           dble(ntmp(isplit,iloc,ial))*
     &           dlog(ftmp(isplit,iloc,ial))
            nn = nn + ntmp(isplit,iloc,ial)
         enddo
         tt = gglgamfn(q + dble(nn)) - gglgamfn(q) + ss
      enddo
      termfsplit = termfsplit - tt

*    ln( pi[ f_{npop+1}^* |fa, drift] / pi[f_{npop+1}^*| ...] ) 
      tt = 0
      q = (1-drifttmp(npop+1))/drifttmp(npop+1)
      do iloc = 1,nlocmax
         ss = 0
         nn = 0
         do ial = 1,nal(iloc)
            ss = ss + gglgamfn(fa(iloc,ial) * q) -
     &           gglgamfn(fa(iloc,ial) * q +
     &           dble(ntmp(npop+1,iloc,ial))) +
     &           dble(ntmp(npop+1,iloc,ial))*
     &           dlog(ftmp(npop+1,iloc,ial))
            nn = nn + ntmp(npop+1,iloc,ial)
         enddo
         tt = gglgamfn(q + dble(nn)) - gglgamfn(q) + ss
      enddo
      termfsplit = termfsplit - tt
      end function termfsplit



*******************************************************
*
*     terme des freq dans le log ratio pour un split 
*     quand f_Alj =1 pour tout j et drift(ipop) =0.5d0
      double precision function termfsplitbis(isplit,npop,npopmax,
     &     nlocmax,nal,nalmax,
     &     f,ftmp,n,ntmp,fa,drift,drifttmp)
      implicit none
      integer npopmax,nlocmax,nalmax,
     &     n(npopmax,nlocmax,nalmax),
     &     ntmp(npopmax,nlocmax,nalmax),nal(nlocmax),
     &     isplit,npop
      double precision f(npopmax,nlocmax,nalmax),drift(npopmax),
     &     ftmp(npopmax,nlocmax,nalmax),drifttmp(npopmax),
     &     fa(nlocmax,nalmax)
      integer iloc, ial,nn
      double precision q,gglgamfn,ss,tt

      
      termfsplitbis  = 0.
*    ln( pi[ f_isplit |...] / pi[f_isplit| fa, drift] ) 
      tt = 0
      q = (1-drift(isplit))/drift(isplit)
      do iloc = 1,nlocmax
         ss = 0
         nn = 0
         do ial = 1,nal(iloc)
            ss = ss + gglgamfn(fa(iloc,ial) * q) -
     &           gglgamfn(fa(iloc,ial) * q +
     &           dble(n(isplit,iloc,ial))) +
     &           dble(n(isplit,iloc,ial))*
     &           dlog(f(isplit,iloc,ial))
            nn = nn + n(isplit,iloc,ial)
         enddo
         tt = gglgamfn(dble(nal(iloc)) + dble(nn)) - 
     &        gglgamfn(dble(nal(iloc))) + ss
      enddo
      termfsplitbis = termfsplitbis + tt

*    ln( pi[ f_isplit^* |fa, drift] / pi[f_isplit^*| ...] ) 
      tt = 0
      q = (1-drifttmp(isplit))/drifttmp(isplit)
      do iloc = 1,nlocmax
         ss = 0
         nn = 0
         do ial = 1,nal(iloc)
            ss = ss + gglgamfn(fa(iloc,ial) * q) -
     &           gglgamfn(fa(iloc,ial) * q +
     &           dble(ntmp(isplit,iloc,ial))) +
     &           dble(ntmp(isplit,iloc,ial))*
     &           dlog(ftmp(isplit,iloc,ial))
            nn = nn + ntmp(isplit,iloc,ial)
         enddo
         tt = gglgamfn(dble(nal(iloc)) + dble(nn)) - 
     &        gglgamfn(dble(nal(iloc))) + ss
      enddo
      termfsplitbis = termfsplitbis - tt

*    ln( pi[ f_{npop+1}^* |fa, drift] / pi[f_{npop+1}^*| ...] ) 
      tt = 0
      q = (1-drifttmp(npop+1))/drifttmp(npop+1)
      do iloc = 1,nlocmax
         ss = 0
         nn = 0
         do ial = 1,nal(iloc)
            ss = ss + gglgamfn(fa(iloc,ial) * q) -
     &           gglgamfn(fa(iloc,ial) * q +
     &           dble(ntmp(npop+1,iloc,ial))) +
     &           dble(ntmp(npop+1,iloc,ial))*
     &           dlog(ftmp(npop+1,iloc,ial))
            nn = nn + ntmp(npop+1,iloc,ial)
         enddo
         tt = gglgamfn(dble(nal(iloc)) + dble(nn)) - 
     &        gglgamfn(dble(nal(iloc))) + ss
      enddo
      termfsplitbis = termfsplitbis - tt

      end function termfsplitbis


*********************************************************************
*     terme des freq dans le log ratio pour un split 
      double precision function termfmerge(ihost,irem,npopmax,
     &     nlocmax,nal,nalmax,
     &     f,ftmp,n,ntmp,fa,drift,drifttmp)
      implicit none
      integer npopmax,nlocmax,nalmax,
     &     n(npopmax,nlocmax,nalmax),
     &     ntmp(npopmax,nlocmax,nalmax),nal(nlocmax),
     &     ihost,irem
      double precision f(npopmax,nlocmax,nalmax),drift(npopmax),
     &     ftmp(npopmax,nlocmax,nalmax),drifttmp(npopmax),
     &     fa(nlocmax,nalmax)
      integer iloc, ial,nn
      double precision q,gglgamfn,ss,tt

      
      termfmerge  = 0.
*     ln( pi[f_host| ...]/pi[f_host|fa,drift])
      tt = 0
      q = (1-drift(ihost))/drift(ihost)
      do iloc = 1,nlocmax
         ss = 0
         nn = 0
         do ial = 1,nal(iloc)
            ss = ss + gglgamfn(fa(iloc,ial) * q) -
     &           gglgamfn(fa(iloc,ial) * q +
     &           dble(n(ihost,iloc,ial))) +
     &           dble(n(ihost,iloc,ial))*
     &           dlog(f(ihost,iloc,ial))
            nn = nn + n(ihost,iloc,ial)
         enddo
         tt = gglgamfn(q + dble(nn)) - gglgamfn(q) + ss
      enddo
      termfmerge = termfmerge + tt

*     ln( pi[f_rem| ...]/pi[f_rem|fa,drift])
      tt = 0
      q = (1-drift(irem))/drift(irem)
      do iloc = 1,nlocmax
         ss = 0
         nn = 0
         do ial = 1,nal(iloc)
            ss = ss + gglgamfn(fa(iloc,ial) * q) -
     &           gglgamfn(fa(iloc,ial) * q +
     &           dble(n(irem,iloc,ial))) +
     &           dble(n(irem,iloc,ial))*
     &           dlog(f(irem,iloc,ial))
            nn = nn + n(irem,iloc,ial)
         enddo
         tt = gglgamfn(q + dble(nn)) - gglgamfn(q) + ss
      enddo
      termfmerge = termfmerge + tt

*     ln( pi[f_host^*| ...]/pi[f_host^*|fa,drift])
      tt = 0
      q = (1-drifttmp(ihost))/drifttmp(ihost)
      do iloc = 1,nlocmax
         ss = 0
         nn = 0
         do ial = 1,nal(iloc)
            ss = ss + gglgamfn(fa(iloc,ial) * q) -
     &           gglgamfn(fa(iloc,ial) * q +
     &           dble(ntmp(ihost,iloc,ial))) +
     &           dble(ntmp(ihost,iloc,ial))*
     &           dlog(ftmp(ihost,iloc,ial))
            nn = nn + ntmp(ihost,iloc,ial)
         enddo
         tt = gglgamfn(q + dble(nn)) - gglgamfn(q) + ss
      enddo
      termfmerge = termfmerge - tt
      end function termfmerge


*
*     terme des freq dans le log ratio pour un split 
*     quand f_Alj =1 pour tout j et drift(ipop) =0.5d0
      double precision function termfmergebis(ihost,irem,npopmax,
     &     nlocmax,nal,nalmax,
     &     f,ftmp,n,ntmp,fa,drift,drifttmp)
      implicit none
      integer npopmax,nlocmax,nalmax,
     &     n(npopmax,nlocmax,nalmax),
     &     ntmp(npopmax,nlocmax,nalmax),nal(nlocmax),
     &     ihost,irem
      double precision f(npopmax,nlocmax,nalmax),drift(npopmax),
     &     ftmp(npopmax,nlocmax,nalmax),drifttmp(npopmax),
     &     fa(nlocmax,nalmax)
      integer iloc, ial,nn
      double precision q,gglgamfn,ss,tt

      
      termfmergebis  = 0.
*     ln( pi[f_host| ...]/pi[f_host|fa,drift])
      tt = 0
      q = (1-drift(ihost))/drift(ihost)
      do iloc = 1,nlocmax
         ss = 0
         nn = 0
         do ial = 1,nal(iloc)
            ss = ss + gglgamfn(fa(iloc,ial) * q) -
     &           gglgamfn(fa(iloc,ial) * q +
     &           dble(n(ihost,iloc,ial))) +
     &           dble(n(ihost,iloc,ial))*
     &           dlog(f(ihost,iloc,ial))
            nn = nn + n(ihost,iloc,ial)
         enddo
         tt = gglgamfn(dble(nal(iloc)) + dble(nn)) - 
     &        gglgamfn(dble(nal(iloc))) + ss
      enddo
      termfmergebis = termfmergebis + tt

*     ln( pi[f_rem| ...]/pi[f_rem|fa,drift])
      tt = 0
      q = (1-drift(irem))/drift(irem)
      do iloc = 1,nlocmax
         ss = 0
         nn = 0
         do ial = 1,nal(iloc)
            ss = ss + gglgamfn(fa(iloc,ial) * q) -
     &           gglgamfn(fa(iloc,ial) * q +
     &           dble(n(irem,iloc,ial))) +
     &           dble(n(irem,iloc,ial))*
     &           dlog(f(irem,iloc,ial))
            nn = nn + n(irem,iloc,ial)
         enddo
         tt = gglgamfn(dble(nal(iloc)) + dble(nn)) - 
     &        gglgamfn(dble(nal(iloc))) + ss
      enddo
      termfmergebis = termfmergebis + tt

*     ln( pi[f_host^*| ...]/pi[f_host^*|fa,drift])
      tt = 0
      q = (1-drifttmp(ihost))/drifttmp(ihost)
      do iloc = 1,nlocmax
         ss = 0
         nn = 0
         do ial = 1,nal(iloc)
            ss = ss + gglgamfn(fa(iloc,ial) * q) -
     &           gglgamfn(fa(iloc,ial) * q +
     &           dble(ntmp(ihost,iloc,ial))) +
     &           dble(ntmp(ihost,iloc,ial))*
     &           dlog(ftmp(ihost,iloc,ial))
            nn = nn + ntmp(ihost,iloc,ial)
         enddo
         tt = gglgamfn(dble(nal(iloc)) + dble(nn)) - 
     &        gglgamfn(dble(nal(iloc))) + ss
      enddo
      termfmergebis = termfmergebis - tt
      end function termfmergebis



*******************************************************
*     term from proposal of frequencies in a split in bdpop9bis
      double precision function termf9(npopmax,nloc,nal,nalmax,n,
     &     f,fa,drift,ipop)
      implicit none 
      integer npopmax,nloc,nal,nalmax,n,ipop,ipop2
      double precision f,fa,drift
      dimension nal(nloc),n(npopmax,nloc,nalmax),
     &     f(npopmax,nloc,nalmax),fa(nloc,nalmax),drift(npopmax)
      integer iloc,ial,nn
      double precision gglgamfn,q
*     next line corrected by gilles on 5/01/08
*      q = drift(ipop) / (1- drift(ipop))
      q = (1-drift(ipop)) / drift(ipop)
      termf9 = 0.d0
      do iloc=1,nloc
         nn = 0
         do ial=1,nal(iloc)
            termf9 = termf9 + 
     &           dble(n(ipop,iloc,ial))*dlog(f(ipop,iloc,ial)) - 
     &           gglgamfn(q*fa(iloc,ial)+dble(n(ipop,iloc,ial))) + 
     &           gglgamfn(q*fa(iloc,ial))
            nn = nn + n(ipop,iloc,ial)
         enddo
*     next line corrected by gilles on 5/01/08
*         termf9 = termf9 + 
*     &        gglgamfn(q+dble(nn)) -  gglgamfn(dble(nn))
         termf9 = termf9 + 
     &        gglgamfn(q+dble(nn)) - gglgamfn(q)
      enddo
      end function termf9

 
*******************************************************
*     term from proposal of frequencies in a split in bdpop9bis
      double precision function termf9bis(npopmax,nloc,nal,nalmax,n,
     &     f,ipop)
      implicit none 
      integer npopmax,nloc,nal,nalmax,n,ipop,ipop2
      double precision f
      dimension nal(nloc),n(npopmax,nloc,nalmax),
     &     f(npopmax,nloc,nalmax)
      integer iloc,ial,nn
      double precision gglgamfn
      termf9bis = 0.d0
      do iloc=1,nloc
         nn = 0
         do ial=1,nal(iloc)
            termf9bis = termf9bis + 
     &           dble(n(ipop,iloc,ial))*dlog(f(ipop,iloc,ial)) - 
     &           gglgamfn(1+dble(n(ipop,iloc,ial))) 
            nn = nn + n(ipop,iloc,ial)
         enddo
         termf9bis = termf9bis + 
     &        gglgamfn(dble(nal(iloc)+nn))
      enddo
      end function termf9bis

 

********************************************************
*     log vraisemblance
      double precision function ll(zz,nindiv,nlocmax,nlocmax2,npopmax,
     &     nalmax,nppmax,c,f,indcell,ploidy)
      implicit none
      integer nindiv,nlocmax,nlocmax2,npopmax,
     &     zz(nindiv,nlocmax2),nppmax,c(nppmax),nalmax,
     &     indcell(nindiv),ploidy
      double precision f(npopmax,nlocmax,nalmax)
      integer iindiv,iloc,ial1,ial2,ipop
      ll = 0
      do iindiv = 1,nindiv
         ipop = c(indcell(iindiv))
         do iloc = 1,nlocmax
            ial1 = zz(iindiv,2*iloc-1)
            ial2 = zz(iindiv,2*iloc)
            if(ial1 .ne. -999) then 
               ll = ll + dlog(f(ipop,iloc,ial1))
            endif
            if((ial2 .ne. -999) .and. (ploidy .eq.2)) then 
               ll = ll + dlog(f(ipop,iloc,ial2))
            endif
            if(((ial1 .ne. ial2)  .and. (ploidy .eq. 2)) .and. 
     &          ((ial1 .ne. -999) .and. (ial2 .ne. -999))) then 
               ll = ll + dlog(2.d0)
            endif
         enddo
      enddo
      end function ll
************************************************************************



************************************************************************
*     log likelihood for the case where 
*     data consist of genotypes and/or quantitative variables
      double precision function llgq(zz,nindiv,nlocmax,nlocmax2,
     &     npopmax,nalmax,nppmax,c,f,indcell,ploidy,
     &     qtc,nqtc,meanqtc,sdqtc,usegeno2,useqtc)
      implicit none
      integer nindiv,nlocmax,nlocmax2,npopmax,
     &     zz(nindiv,nlocmax2),nppmax,c(nppmax),nalmax,
     &     indcell(nindiv),ploidy,nqtc,usegeno2,useqtc
      double precision f(npopmax,nlocmax,nalmax),qtc(nindiv,nqtc),
     &     meanqtc(npopmax,nqtc),sdqtc(npopmax,nqtc)
      integer iindiv,iloc,ial1,ial2,ipop,iqtc
      double precision ggdnorm
      llgq = 0
      do iindiv = 1,nindiv
         ipop = c(indcell(iindiv))
*     contrib genotypes
         if(usegeno2 .eq. 1) then 
            do iloc = 1,nlocmax
               ial1 = zz(iindiv,2*iloc-1)
               ial2 = zz(iindiv,2*iloc)
               if(ial1 .ne. -999) then 
                  llgq = llgq + dlog(f(ipop,iloc,ial1))
               endif
               if((ial2 .ne. -999) .and. (ploidy .eq.2)) then 
                  llgq = llgq + dlog(f(ipop,iloc,ial2))
               endif
               if(((ial1 .ne. ial2)  .and. (ploidy .eq. 2)) .and. 
     &              ((ial1 .ne. -999) .and. (ial2 .ne. -999))) then 
                  llgq = llgq + dlog(2.d0)
               endif
            enddo
         endif
*     contrib quantit. variables
         if(useqtc .eq. 1) then 
            do iqtc = 1,nqtc
               llgq = llgq + ggdnorm(qtc(iindiv,iqtc),
     &              meanqtc(ipop,iqtc),sdqtc(ipop,iqtc),1)
            enddo
         endif
      enddo
      end function llgq
************************************************************************
 



************************************************************************
*     log likelihood for the case where 
*     data consist of diploid genotypes and/or haploid genotypes 
*     and/or  quantitative variables
      double precision function llallvar2(yy,z,ql,nindiv,nlocd,
     &     nloch,nql,ncolt,npopmax,nalmax,nppmax,c,f,indcell,
     &     qtc,nqtc,meanqtc,sdqtc,usegeno2,usegeno1,useql,useqtc,
     &     ploidy)
      implicit none
      integer nindiv,nlocd,nlocd2,nloch,nql,ncolt,npopmax,
     &     yy(nindiv,2*nlocd+2*nloch),z(nindiv,nloch),ql(nindiv,nql),
     &     nppmax,c(nppmax),nalmax,indcell(nindiv),ploidy,nqtc,
     &     usegeno2,usegeno1,useql,useqtc
      double precision f(npopmax,ncolt,nalmax),qtc(nindiv,nqtc),
     &     meanqtc(npopmax,nqtc),sdqtc(npopmax,nqtc)
      integer iindiv,iloc,ial1,ial2,ipop,iqtc
      double precision ggdnorm
      llallvar2 = 0
      do iindiv = 1,nindiv
         ipop = c(indcell(iindiv))
*     contrib diploid genotypes
         if(usegeno2 .eq. 1) then 
            do iloc = 1,nlocd
               ial1 = yy(iindiv,2*iloc-1)
               ial2 = yy(iindiv,2*iloc)
               if(ial1 .ne. -999) then 
                  llallvar2 = llallvar2 + dlog(f(ipop,iloc,ial1))
               endif
               if(ial2 .ne. -999) then 
                  llallvar2 = llallvar2 + dlog(f(ipop,iloc,ial2))
               endif
               if((ial1 .ne. ial2) .and. 
     &              ((ial1 .ne. -999) .and. (ial2 .ne. -999))) then 
                  llallvar2 = llallvar2 + dlog(2.d0)
               endif
            enddo
         endif
*     contrib haploid genotypes
         if(usegeno1 .eq. 1) then 
            if(ploidy .eq. 2) then 
               do iloc = 1,nloch
                  ial1 = yy(iindiv,2*nlocd+2*iloc-1)
                  ial2 = yy(iindiv,2*nlocd+2*iloc)
                  if(ial1 .ne. -999) then 
                     llallvar2 = llallvar2 
     &                    + dlog(f(ipop,nlocd+iloc,ial1))
                  endif
                  if(ial2 .ne. -999) then 
                     llallvar2 = llallvar2 
     &                    + dlog(f(ipop,nlocd+iloc,ial2))
                  endif
                  if((ial1 .ne. ial2) .and. 
     &                 ((ial1 .ne. -999) .and. (ial2 .ne. -999))) then 
                     llallvar2 = llallvar2 + dlog(2.d0)
                  endif
               enddo
            endif
            if(ploidy .eq. 1) then
               do iloc = 1,nloch
                  ial1 = z(iindiv,iloc)
                  if(ial1 .ne. -999) then 
                     llallvar2 = llallvar2 
     &                    + dlog(f(ipop,nlocd+iloc,ial1))
                  endif
               enddo
            endif
         endif
*     contrib qualit. variables
         if(useql .eq. 1) then 
            do iloc = 1,nql
               ial1 = ql(iindiv,iloc)
               if(ial1 .ne. -999) then 
                  llallvar2 = llallvar2 + 
     &                 dlog(f(ipop,nlocd+nloch+iloc,ial1))
               endif
            enddo
         endif
*     contrib quantit. variables
         if(useqtc .eq. 1) then 
            do iqtc = 1,nqtc
               llallvar2 = llallvar2 + ggdnorm(qtc(iindiv,iqtc),
     &              meanqtc(ipop,iqtc),sdqtc(ipop,iqtc),1)
            enddo
         endif
      enddo
      end function llallvar2
************************************************************************
 



************************************************************************
*     log de la proba a posteriori du vecteur de parametres
      double precision function lpp(lambdamax,lambda,zz,npop,npp,drift,
     &     f,fa,c,nppmax,
     &     nindiv,nlocmax2,npopmax,nlocmax,nalmax,indcell,nal,
     &     fmodel,xlim,ylim,shape1,shape2,ploidy)
      implicit none
      integer nindiv,nlocmax2,npop,npopmax,nlocmax,nalmax,
     &     npp,nppmax,zz(nindiv,nlocmax2),indcell(nindiv),
     &     c(nppmax),nal(nlocmax),fmodel,ploidy
      double precision drift(npopmax),f(npopmax,nlocmax,nalmax),
     &     fa(nlocmax,nalmax),lambdamax,lambda,xlim(2),ylim(2)
      integer ipp,ipop,iloc,ial
      double precision gglgamfn,shape1,shape2,ll,lg
      lpp = 0
*     contrib npop
      lpp = lpp - dlog(dble(npopmax))
*     contrib lambda
      lpp = - dlog(lambdamax)
*     contrib npp
      lpp = lpp - lambda + dble(npp)*dlog(lambda) - 
     &     gglgamfn(dble(npp+1))
*     contrib c
      lpp = lpp - dble(npp)*dlog(dble(npop))
*     contrib u
      lpp = lpp - dble(npp)*dlog((xlim(2)-xlim(1))*(ylim(2)-ylim(1)))
*     contrib freq
      if(fmodel .eq. 0) then
         lg = 0
         do iloc = 1,nlocmax
            lg = lg + gglgamfn(dble(nal(iloc)))
         enddo
         lpp= lpp + dble(npop) * lg
      endif
      if(fmodel .eq. 1) then
*     contrib ancestral freq
         lg = 0
         do iloc = 1,nlocmax
            lg = lg + gglgamfn(dble(nal(iloc)))
         enddo
         lpp= lpp + lg
*     contrib drifts
         do ipop = 1,npop
            lpp = lpp + 
     &           gglgamfn(shape1+shape2) - gglgamfn(shape1) 
     &           - gglgamfn(shape2) + (shape1-1)*dlog(drift(ipop)) + 
     &           (shape2-1)*dlog(1-drift(ipop))
         enddo
*     contrib freq
         do ipop = 1,npop
            do iloc = 1,nlocmax
               lpp = lpp + gglgamfn((1-drift(ipop))/drift(ipop))
               do ial=1,nal(iloc)
                  lpp = lpp -gglgamfn(fa(iloc,ial)*
     &                 (1-drift(ipop))/drift(ipop)) + 
     &                 (fa(iloc,ial)* 
     &                 (1-drift(ipop))/drift(ipop)-1)*
     &              dlog(f(ipop,iloc,ial))
               enddo
            enddo
         enddo
      endif
      
*     contrib likelihood
      lpp = lpp + ll(zz,nindiv,nlocmax,nlocmax2,npopmax,
     &     nalmax,nppmax,c,f,indcell,ploidy)
      end function lpp
************************************************************************


c$$$************************************************************************
c$$$*     log de la proba a posteriori du vecteur de parametres
c$$$      double precision function lppgq(lambdamax,lambda,zz,npop,npp,
c$$$     &     drift,f,fa,c,nppmax,
c$$$     &     nindiv,nlocmax2,npopmax,nlocmax,nalmax,indcell,nal,
c$$$     &     fmodel,xlim,ylim,shape1,shape2,ploidy,
c$$$     &     qtc,nqtc,meanqtc,sdqtc,usegeno2,useqtc,
c$$$     &     ksiqtc,kappaqtc,alphaqtc,betaqtc)
c$$$      implicit none
c$$$      integer nindiv,nlocmax2,npop,npopmax,nlocmax,nalmax,
c$$$     &     npp,nppmax,zz(nindiv,nlocmax2),indcell(nindiv),
c$$$     &     c(nppmax),nal(nlocmax),fmodel,ploidy,nqtc,usegeno2,useqtc
c$$$      double precision drift(npopmax),f(npopmax,nlocmax,nalmax),
c$$$     &     fa(nlocmax,nalmax),lambdamax,lambda,xlim(2),ylim(2),
c$$$     &     qtc(nindiv,nqtc),
c$$$     &     meanqtc(npopmax,nqtc),sdqtc(npopmax,nqtc),
c$$$     &     ksiqtc,kappaqtc,alphaqtc,betaqtc
c$$$      integer ipp,ipop,iloc,ial,iqtc
c$$$      double precision gglgamfn,shape1,shape2,llgq,lg,ggdnorm,ggdgamma
c$$$*     contrib npop
c$$$      lppgq = lppgq - dlog(dble(npopmax))
c$$$*     contrib lambda
c$$$      lppgq = - dlog(lambdamax)
c$$$*     contrib npp
c$$$      lppgq = lppgq - lambda + dble(npp)*dlog(lambda) - 
c$$$     &     gglgamfn(dble(npp+1))
c$$$*     contrib c
c$$$      lppgq = lppgq - dble(npp)*dlog(dble(npop))
c$$$*     contrib u
c$$$      lppgq = lppgq - npp*dlog((xlim(2)-xlim(1))*(ylim(2)-ylim(1)))
c$$$      if(usegeno2 .eq. 1) then 
c$$$*     contrib freq
c$$$         if(fmodel .eq. 0) then
c$$$            lg = 0
c$$$            do iloc = 1,nlocmax
c$$$               lg = lg + gglgamfn(dble(nal(iloc)))
c$$$            enddo
c$$$            lppgq= lppgq + dble(npop) * lg
c$$$         endif
c$$$         if(fmodel .eq. 1) then
c$$$*     contrib ancestral freq
c$$$            lg = 0
c$$$            do iloc = 1,nlocmax
c$$$               lg = lg + gglgamfn(dble(nal(iloc)))
c$$$            enddo
c$$$            lppgq= lppgq + lg
c$$$*     contrib drifts
c$$$            do ipop = 1,npop
c$$$               lppgq = lppgq + 
c$$$     &              gglgamfn(shape1+shape2) - gglgamfn(shape1) 
c$$$     &              - gglgamfn(shape2) + (shape1-1)*dlog(drift(ipop)) + 
c$$$     &              (shape2-1)*dlog(1-drift(ipop))
c$$$            enddo
c$$$*     contrib freq
c$$$            do ipop = 1,npop
c$$$               do iloc = 1,nlocmax
c$$$                  lppgq = lppgq + gglgamfn((1-drift(ipop))/drift(ipop))
c$$$                  do ial=1,nal(iloc)
c$$$                     lppgq = lppgq -gglgamfn(fa(iloc,ial)*
c$$$     &                    (1-drift(ipop))/drift(ipop)) + 
c$$$     &                    (fa(iloc,ial)* 
c$$$     &                    (1-drift(ipop))/drift(ipop)-1)*
c$$$     &                    dlog(f(ipop,iloc,ial))
c$$$                  enddo
c$$$               enddo
c$$$            enddo
c$$$         endif
c$$$      endif
c$$$      if(useqtc .eq. 1) then
c$$$*     contrib mean of quantit. variables
c$$$         do iqtc = 1,nqtc
c$$$            do ipop = 1,npop
c$$$               lppgq = lppgq + 
c$$$     &            ggdnorm(meanqtc(ipop,iqtc),ksiqtc,
c$$$     &              dsqrt(1/kappaqtc),1)
c$$$            enddo
c$$$         enddo
c$$$*     contrib precision of quantit. variables
c$$$         do iqtc = 1,nqtc
c$$$            do ipop = 1,npop
c$$$               lppgq = lppgq + 
c$$$     &              ggdgamma(1/sdqtc(ipop,iqtc)**2,
c$$$     &              alphaqtc,1/betaqtc,1)
c$$$            enddo
c$$$         enddo
c$$$      endif
c$$$      
c$$$*     contrib likelihood
c$$$      lppgq = lppgq + llgq(zz,nindiv,nlocmax,nlocmax2,npopmax,
c$$$     &     nalmax,nppmax,c,f,indcell,ploidy,
c$$$     &     qtc,nqtc,meanqtc,sdqtc,usegeno2,useqtc)
c$$$      end function lppgq
c$$$************************************************************************
c$$$



************************************************************************
*     log de la proba a priori du vecteur de parametres
      double precision function lpriorallvar(lambdamax,lambda,
     &     npop,npp,drift,f,fa,c,nppmax,nindiv,npopmax,nlocd,nloch,nql,
     &     ncolt,nal,nalmax,indcell,fmodel,xlim,ylim,shape1,shape2,
     &     nqtc,meanqtc,sdqtc,ksiqtc,kappaqtc,alphaqtc,betaqtc,
     &     usegeno2,usegeno1,useql,useqtc)
      implicit none
      integer nindiv,npop,npopmax,nlocd,nloch,nql,ncolt,nalmax,
     &     npp,nppmax,indcell(nindiv),
     &     c(nppmax),nal(ncolt),fmodel,ploidy,nqtc,usegeno2,
     &     usegeno1,useql,useqtc
      double precision drift(npopmax),f(npopmax,ncolt,nalmax),
     &     fa(ncolt,nalmax),lambdamax,lambda,xlim(2),ylim(2),
     &     meanqtc(npopmax,nqtc),sdqtc(npopmax,nqtc),
     &     ksiqtc,kappaqtc,alphaqtc,betaqtc
      integer ipp,ipop,iloc,ial,iqtc
      double precision gglgamfn,shape1,shape2,llgq,lg,ggdnorm,ggdgamma
      dimension  ksiqtc(nqtc),kappaqtc(nqtc),alphaqtc(nqtc),
     &     betaqtc(nqtc)

      lpriorallvar = 0
*     contrib npop
      lpriorallvar = lpriorallvar - dlog(dble(npopmax))
*     contrib lambda
      lpriorallvar = lpriorallvar - dlog(lambdamax)
*     contrib npp
      lpriorallvar = lpriorallvar - lambda + dble(npp)*dlog(lambda) - 
     &     gglgamfn(dble(npp+1))
*     contrib c
      lpriorallvar = lpriorallvar - dble(npp)*dlog(dble(npop))
*     contrib u
      lpriorallvar = lpriorallvar - 
     &     dble(npp)*dlog((xlim(2)-xlim(1))*(ylim(2)-ylim(1)))
*     contrib drifts
      if((((usegeno2 .eq. 1) .or. (usegeno1 .eq. 1)) 
     &     .or. (useql .eq. 1)) .and. (fmodel .eq. 1)) then 
         do ipop = 1,npop
            lpriorallvar = lpriorallvar + 
     &           gglgamfn(shape1+shape2) - gglgamfn(shape1) - 
     &           gglgamfn(shape2) + (shape1-1)*dlog(drift(ipop)) + 
     &           (shape2-1)*dlog(1-drift(ipop))
         enddo
      endif
*     contrib ancestral freq
      if(fmodel .eq. 1) then
         if(usegeno2 .eq. 1) then
            lg = 0
            do iloc = 1,nlocd
               lg = lg + gglgamfn(dble(nal(iloc)))
            enddo
            lpriorallvar= lpriorallvar + lg
         endif
         if(usegeno1 .eq. 1) then
            lg = 0
            do iloc = 1,nloch
               lg = lg + gglgamfn(dble(nal(nlocd+iloc)))
            enddo
            lpriorallvar= lpriorallvar + lg
         endif
         if(useql .eq. 1) then
            lg = 0
            do iloc = 1,nql
               lg = lg + gglgamfn(dble(nal(nlocd+nloch+iloc)))
            enddo
            lpriorallvar= lpriorallvar + lg
         endif
      endif
*     contrib freq
      if(fmodel .eq. 0) then
         if(usegeno2 .eq. 1) then 
            lg = 0
            do iloc = 1,nlocd
               lg = lg + gglgamfn(dble(nal(iloc)))
            enddo
            lpriorallvar= lpriorallvar + dble(npop) * lg
         endif
         if(usegeno1 .eq. 1) then 
            lg = 0
            do iloc = 1,nloch
               lg = lg + gglgamfn(dble(nal(nlocd+iloc)))
            enddo
            lpriorallvar= lpriorallvar + dble(npop) * lg
         endif
        if(useql .eq. 1) then 
            lg = 0
            do iloc = 1,nql
               lg = lg + gglgamfn(dble(nal(nlocd+nloch+iloc)))
            enddo
            lpriorallvar= lpriorallvar + dble(npop) * lg
         endif
      endif
      if(fmodel .eq. 1) then
         if(usegeno2 .eq. 1) then 
            do ipop = 1,npop
               do iloc = 1,nlocd
                  lpriorallvar = lpriorallvar + 
     &                 gglgamfn((1-drift(ipop))/drift(ipop))
                  do ial=1,nal(iloc)
                     lpriorallvar = lpriorallvar -gglgamfn(fa(iloc,ial)*
     &                    (1-drift(ipop))/drift(ipop)) + (fa(iloc,ial)*
     &                    (1-drift(ipop))/drift(ipop)-1)*
     &                    dlog(f(ipop,iloc,ial))
                  enddo
               enddo
            enddo
         endif
         if(usegeno1 .eq. 1) then 
            do ipop = 1,npop
               do iloc = 1,nloch
                  lpriorallvar = lpriorallvar + 
     &                 gglgamfn((1-drift(ipop))/drift(ipop))
                  do ial=1,nal(nlocd+iloc)
                     lpriorallvar = lpriorallvar - 
     &                    gglgamfn(fa(nlocd+iloc,ial)*
     &                    (1-drift(ipop))/drift(ipop)) + 
     &                    (fa(nlocd+iloc,ial)*
     &                    (1-drift(ipop))/drift(ipop)-1)*
     &                    dlog(f(ipop,nlocd+iloc,ial))
                  enddo
               enddo
            enddo
         endif
         if(useql .eq. 1) then 
            do ipop = 1,npop
               do iloc = 1,nql
                  lpriorallvar = lpriorallvar + 
     &                 gglgamfn((1-drift(ipop))/drift(ipop))
                  do ial=1,nal(nlocd+nloch+iloc)
                     lpriorallvar = lpriorallvar - 
     &                    gglgamfn(fa(nlocd+nloch+iloc,ial)*
     &                    (1-drift(ipop))/drift(ipop)) + 
     &                    (fa(nlocd+nloch+iloc,ial)*
     &                    (1-drift(ipop))/drift(ipop)-1)*
     &                    dlog(f(ipop,nlocd+nloch+iloc,ial))
                  enddo
               enddo
            enddo
         endif
      endif
      if(useqtc .eq. 1) then
*     contrib precision of quantit. variables
         do iqtc = 1,nqtc
            do ipop = 1,npop
               lpriorallvar = lpriorallvar + 
     &              ggdgamma(1/sdqtc(ipop,iqtc)**2,
     &              alphaqtc(iqtc),1/betaqtc(iqtc),1)
            enddo
         enddo
*     contrib mean of quantit. variables
         do iqtc = 1,nqtc
            do ipop = 1,npop
               lpriorallvar = lpriorallvar + 
     &            ggdnorm(meanqtc(ipop,iqtc),ksiqtc(iqtc),
     &              sdqtc(ipop,iqtc)/dsqrt(kappaqtc(iqtc)),1)
            enddo
         enddo
      endif
      end function lpriorallvar
************************************************************************





***********************************************************************
*     sample freq from full contitionnal pi(f|u,c,zz)
*     drifts and ancestral are not changed
      subroutine samplef(npop,npopmax,nloc,nlocmax,
     &     nal,nalmax,ipop1,ipop2,ftmp,a,ptmp,ntmp,alpha)
      implicit none
      integer npop,npopmax,nloc,nlocmax,nal,
     &     nalmax,ntmp,
     &     ipop1,ipop2
      double precision ftmp,a
      dimension nal(nlocmax),
     &     ntmp(npopmax,nloc,nalmax),
     &     ftmp(npopmax,nlocmax,nalmax),
     &     a(nalmax),ptmp(nalmax)
      integer iloc,ial,ipop
      double precision ptmp
      double precision alpha

c      write(*,*) 'debut sample f'

*     new freq for pop ipop1
      do iloc=1,nloc
         do ial = 1,nal(iloc)
            a(ial) = alpha + dble(ntmp(ipop1,iloc,ial))
         enddo
         call dirichlet(nal(iloc),nalmax,a,ptmp)
         do ial=1,nal(iloc)
            ftmp(ipop1,iloc,ial)  = ptmp(ial)
         enddo
      enddo

*     for pop ipop2
      if(ipop2 .ne. ipop1) then
         do iloc=1,nloc
            do ial = 1,nal(iloc)
               a(ial) = 1. + dble(ntmp(ipop2,iloc,ial))
            enddo
            call dirichlet(nal(iloc),nalmax,a,ptmp)
            do ial=1,nal(iloc)
               ftmp(ipop2,iloc,ial)  = ptmp(ial)
            enddo
         enddo
      endif
c      write(*,*) 'end  sample f'
      end subroutine samplef 
***********************************************************************






***********************************************************************
*     sample freq from full contitionnal pi(f|u,c,zz)
*     drifts and ancestral are not changed
      subroutine samplefallvar(npopmax,nlocd,nloch,nql,ncolt,
     &     nal,nalmax,ipop1,ipop2,ftmp,a,ptmp,ntmp,alpha,
     &     usegeno2,usegeno1,useql)
      implicit none
      integer npopmax,nlocd,nloch,nql,ncolt,nal,nalmax,ntmp,ipop1,ipop2,
     &     usegeno2,usegeno1,useql
      double precision ftmp,a
      integer iloc,ial,ipop
      double precision ptmp,alpha
      dimension nal(ncolt),ntmp(npopmax,ncolt,nalmax),
     &     ftmp(npopmax,ncolt,nalmax),a(nalmax),ptmp(nalmax)
c      write(*,*) 'debut sample f'
      if(usegeno2 .eq. 1) then 
*     new freq for pop ipop1
         do iloc=1,nlocd
            do ial = 1,nal(iloc)
               a(ial) = alpha + dble(ntmp(ipop1,iloc,ial))
            enddo
            call dirichlet(nal(iloc),nalmax,a,ptmp)
            do ial=1,nal(iloc)
               ftmp(ipop1,iloc,ial)  = ptmp(ial)
            enddo
         enddo
*     for pop ipop2
         if(ipop2 .ne. ipop1) then
            do iloc=1,nlocd
               do ial = 1,nal(iloc)
                  a(ial) = 1. + dble(ntmp(ipop2,iloc,ial))
               enddo
               call dirichlet(nal(iloc),nalmax,a,ptmp)
               do ial=1,nal(iloc)
                  ftmp(ipop2,iloc,ial)  = ptmp(ial)
               enddo
            enddo
         endif
      endif
      if(usegeno1 .eq. 1) then 
*     new freq for pop ipop1
         do iloc=1,nloch
            do ial = 1,nal(nlocd+iloc)
               a(ial) = alpha + dble(ntmp(ipop1,nlocd+iloc,ial))
            enddo
            call dirichlet(nal(nlocd+iloc),nalmax,a,ptmp)
            do ial=1,nal(nlocd+iloc)
               ftmp(ipop1,nlocd+iloc,ial)  = ptmp(ial)
            enddo
         enddo
*     for pop ipop2
         if(ipop2 .ne. ipop1) then
            do iloc=1,nloch
               do ial = 1,nal(nlocd+iloc)
                  a(ial) = alpha + dble(ntmp(ipop2,nlocd+iloc,ial))
               enddo
               call dirichlet(nal(nlocd+iloc),nalmax,a,ptmp)
               do ial=1,nal(nlocd+iloc)
                  ftmp(ipop2,nlocd+iloc,ial)  = ptmp(ial)
               enddo
            enddo
         endif
      endif
      if(useql .eq. 1) then 
*     new freq for pop ipop1
         do iloc=1,nql
            do ial = 1,nal(nlocd+nloch+iloc)
               a(ial) = alpha + dble(ntmp(ipop1,nlocd+nloch+iloc,ial))
            enddo
            call dirichlet(nal(nlocd+nloch+iloc),nalmax,a,ptmp)
            do ial=1,nal(nlocd+nloch+iloc)
               ftmp(ipop1,nlocd+nloch+iloc,ial)  = ptmp(ial)
            enddo
         enddo
*     for pop ipop2
         if(ipop2 .ne. ipop1) then
            do iloc=1,nql
               do ial = 1,nal(nlocd+nloch+iloc)
                  a(ial) = alpha+ dble(ntmp(ipop2,nlocd+nloch+iloc,ial))
               enddo
               call dirichlet(nal(nlocd+nloch+iloc),nalmax,a,ptmp)
               do ial=1,nal(nlocd+nloch+iloc)
                  ftmp(ipop2,nlocd+nloch+iloc,ial)  = ptmp(ial)
               enddo
            enddo
         endif
      endif
c      write(*,*) 'end  sample f'
      end subroutine samplefallvar 
***********************************************************************







***********************************************************************
*     sample freq from full contitionnal pi(f|u,c,zz)
*     drifts and ancestral are not changed
      subroutine samplefallvar4(npop,npopmax,nlocd,nloch,nql,ncolt,
     &     nal,nalmax,ftmp,a,ptmp,ntmp,alpha,
     &     usegeno2,usegeno1,useql)
      implicit none
      integer npop,npopmax,nlocd,nloch,nql,ncolt,nal,nalmax,ntmp,
     &     usegeno2,usegeno1,useql
      double precision ftmp,a
      integer ipop,iloc,ial
      double precision ptmp,alpha
      dimension nal(ncolt),ntmp(npopmax,ncolt,nalmax),
     &     ftmp(npopmax,ncolt,nalmax),a(nalmax),ptmp(nalmax)
c      write(*,*) 'debut sample f'
      do ipop=1,npop
         if(usegeno2 .eq. 1) then 
*     new freq for pop ipop
            do iloc=1,nlocd
               do ial = 1,nal(iloc)
                  a(ial) = alpha + dble(ntmp(ipop,iloc,ial))
               enddo
               call dirichlet(nal(iloc),nalmax,a,ptmp)
               do ial=1,nal(iloc)
                  ftmp(ipop,iloc,ial)  = ptmp(ial)
               enddo
            enddo
         endif
         if(usegeno1 .eq. 1) then 
*     new freq for pop ipop
            do iloc=1,nloch
               do ial = 1,nal(nlocd+iloc)
                  a(ial) = alpha + dble(ntmp(ipop,nlocd+iloc,ial))
               enddo
               call dirichlet(nal(nlocd+iloc),nalmax,a,ptmp)
               do ial=1,nal(nlocd+iloc)
                  ftmp(ipop,nlocd+iloc,ial)  = ptmp(ial)
               enddo
            enddo
         endif
         if(useql .eq. 1) then 
*     new freq for pop ipop
            do iloc=1,nql
               do ial = 1,nal(nlocd+nloch+iloc)
                  a(ial) = alpha + dble(ntmp(ipop,nlocd+nloch+iloc,ial))
               enddo
               call dirichlet(nal(nlocd+nloch+iloc),nalmax,a,ptmp)
               do ial=1,nal(nlocd+nloch+iloc)
                  ftmp(ipop,nlocd+nloch+iloc,ial)  = ptmp(ial)
               enddo
            enddo
         endif
      enddo
      end subroutine samplefallvar4
***********************************************************************



***********************************************************************
*     sample freq from full contitionnal pi(f|u,c,zz)
*     drifts and ancestral are not changed
*     if CFM for allele frequencies
      subroutine samplef2(npop,npopmax,nloc,nlocmax,
     &     nal,nalmax,ipop1,ipop2,f,ftmp,
     &     fa,drift,a,ptmp,ntmp)
      implicit none
      integer npop,npopmax,nloc,nlocmax,nal,
     &     nalmax,ntmp,
     &     ipop1,ipop2
      double precision f,drift,ftmp,fa,a
      dimension nal(nlocmax),
     &     ntmp(npopmax,nloc,nalmax),
     &     f(npopmax,nlocmax,nalmax),drift(npopmax),
     &     ftmp(npopmax,nlocmax,nalmax),
     &     fa(nlocmax,nalmax),a(nalmax),ptmp(nalmax)
      integer iloc,ial,ipop
      double precision ptmp

*     init ftmp 
c$$$      do ipop = 1,npop
c$$$         do iloc=1,nloc
c$$$            do ial=1,nal(iloc)
c$$$               ftmp(ipop,iloc,ial) = f(ipop,iloc,ial)
c$$$            enddo 
c$$$         enddo
c$$$      enddo

      do iloc=1,nloc
         do ial=1,nal(iloc)
            ftmp(ipop1,iloc,ial) = f(ipop1,iloc,ial)
         enddo 
      enddo
      do iloc=1,nloc
         do ial=1,nal(iloc)
            ftmp(ipop2,iloc,ial) = f(ipop2,iloc,ial)
         enddo 
      enddo


*     new freq
*     for pop ipop1
      do iloc=1,nloc
         do ial = 1,nal(iloc)
            a(ial) = fa(iloc,ial)*(1-drift(ipop1))/
     &           drift(ipop1)+dble(ntmp(ipop1,iloc,ial))
         enddo
         call dirichlet(nal(iloc),nalmax,a,ptmp)
         do ial=1,nal(iloc)
            ftmp(ipop1,iloc,ial)  = ptmp(ial)
         enddo
      enddo

*     for pop ipop2
      if(ipop2 .ne. ipop1) then
         do iloc=1,nloc
            do ial = 1,nal(iloc)
               a(ial) = fa(iloc,ial)*(1-drift(ipop2))/
     &              drift(ipop2)+dble(ntmp(ipop2,iloc,ial))
            enddo
            call dirichlet(nal(iloc),nalmax,a,ptmp)
            do ial=1,nal(iloc)
               ftmp(ipop2,iloc,ial)  = ptmp(ial)
            enddo
         enddo
      endif
      end subroutine samplef2 
***********************************************************************




***********************************************************************
*     sample freq from full contitionnal pi(f|u,c,zz)
*     drifts and ancestral are not changed
*     if CFM for allele frequencies
      subroutine samplef2allvar(npop,npopmax,nlocd,nloch,nql,ncolt,
     &     nal,nalmax,ipop1,ipop2,f,ftmp,fa,drift,a,ptmp,ntmp,
     &     usegeno2,usegeno1,useql)
      implicit none
      integer npop,npopmax,nlocd,nloch,nql,ncolt,nal,nalmax,ntmp,
     &     ipop1,ipop2,usegeno2,usegeno1,useql
      double precision f,drift,ftmp,fa,a
      dimension nal(ncolt),
     &     ntmp(npopmax,ncolt,nalmax),
     &     f(npopmax,ncolt,nalmax),drift(npopmax),
     &     ftmp(npopmax,ncolt,nalmax),
     &     fa(ncolt,nalmax),a(nalmax),ptmp(nalmax)
      integer iloc,ial,ipop
      double precision ptmp

      if(usegeno2 .eq. 1) then 
*     for pop ipop1
         do iloc=1,nlocd
            do ial = 1,nal(iloc)
               a(ial) = fa(iloc,ial)*(1-drift(ipop1))/
     &              drift(ipop1)+dble(ntmp(ipop1,iloc,ial))
            enddo
            call dirichlet(nal(iloc),nalmax,a,ptmp)
            do ial=1,nal(iloc)
               ftmp(ipop1,iloc,ial)  = ptmp(ial)
            enddo
         enddo
*     for pop ipop2
         if(ipop2 .ne. ipop1) then
            do iloc=1,nlocd
               do ial = 1,nal(iloc)
                  a(ial) = fa(iloc,ial)*(1-drift(ipop2))/
     &                 drift(ipop2)+dble(ntmp(ipop2,iloc,ial))
               enddo
               call dirichlet(nal(iloc),nalmax,a,ptmp)
               do ial=1,nal(iloc)
                  ftmp(ipop2,iloc,ial)  = ptmp(ial)
               enddo
            enddo
         endif
      endif
      if(usegeno1 .eq. 1) then 
*     for pop ipop1
         do iloc=1,nloch
            do ial = 1,nal(nlocd+iloc)
               a(ial) = fa(nlocd+iloc,ial)*(1-drift(ipop1))/
     &              drift(ipop1)+dble(ntmp(ipop1,nlocd+iloc,ial))
            enddo
            call dirichlet(nal(nlocd+iloc),nalmax,a,ptmp)
            do ial=1,nal(nlocd+iloc)
               ftmp(ipop1,nlocd+iloc,ial)  = ptmp(ial)
            enddo
         enddo
*     for pop ipop2
         if(ipop2 .ne. ipop1) then
            do iloc=1,nloch
               do ial = 1,nal(nlocd+iloc)
                  a(ial) = fa(nlocd+iloc,ial)*(1-drift(ipop2))/
     &                 drift(ipop2)+dble(ntmp(ipop2,nlocd+iloc,ial))
               enddo
               call dirichlet(nal(nlocd+iloc),nalmax,a,ptmp)
               do ial=1,nal(nlocd+iloc)
                  ftmp(ipop2,nlocd+iloc,ial)  = ptmp(ial)
               enddo
            enddo
         endif
      endif
      if(useql .eq. 1) then 
*     for pop ipop1
         do iloc=1,nql
            do ial = 1,nal(nlocd+nloch+iloc)
               a(ial) = fa(nlocd+nloch+iloc,ial)*(1-drift(ipop1))/
     &              drift(ipop1)+dble(ntmp(ipop1,nlocd+nloch+iloc,ial))
            enddo
            call dirichlet(nal(nlocd+nloch+iloc),nalmax,a,ptmp)
            do ial=1,nal(nlocd+nloch+iloc)
               ftmp(ipop1,nlocd+nloch+iloc,ial)  = ptmp(ial)
            enddo
         enddo
*     for pop ipop2
         if(ipop2 .ne. ipop1) then
            do iloc=1,nql
               do ial = 1,nal(nlocd+nloch+iloc)
                  a(ial) = fa(nlocd+nloch+iloc,ial)*(1-drift(ipop2))/
     &                 drift(ipop2)+
     &                 dble(ntmp(ipop2,nlocd+nloch+iloc,ial))
               enddo
               call dirichlet(nal(nlocd+nloch+iloc),nalmax,a,ptmp)
               do ial=1,nal(nlocd+nloch+iloc)
                  ftmp(ipop2,nlocd+nloch+iloc,ial)  = ptmp(ial)
               enddo
            enddo
         endif
      endif
      end subroutine samplef2allvar 
***********************************************************************


***********************************************************************
*     log of ratio of proposals in a joint update of c and f
*     (in subroutine udcf)
      double precision function lrpf(npopmax,nloc,nal,nalmax,n,
     &     ntmp,f,ftmp,ipop1,ipop2)
      implicit none 
      integer npopmax,nloc,nal,nalmax,n,ntmp,ipop1,ipop2
      double precision f,ftmp
      dimension nal(nloc),n(npopmax,nloc,nalmax),
     &     ntmp(npopmax,nloc,nalmax),f(npopmax,nloc,nalmax),
     &     ftmp(npopmax,nloc,nalmax)
      integer iloc,ial,n1,n2,ntmp1,ntmp2
      double precision gglgamfn
      double precision alpha
      alpha = 1

      lrpf = 0.d0
      do iloc=1,nloc
c     write(*,*) 'iloc=',iloc
         n1 = 0
         n2 = 0
         ntmp1 = 0
         ntmp2 = 0
         do ial=1,nal(iloc)
c     write(*,*) 'ial=',ial
            lrpf = lrpf + 
     &     (n(ipop1,iloc,ial)+alpha-1)*dlog(f(ipop1,iloc,ial)) 
     &   + (n(ipop2,iloc,ial)+alpha-1)*dlog(f(ipop2,iloc,ial))
     &   - (ntmp(ipop1,iloc,ial)+alpha-1)*dlog(ftmp(ipop1,iloc,ial)) 
     &   - (ntmp(ipop2,iloc,ial)+alpha-1)*dlog(ftmp(ipop2,iloc,ial)) 
     &           - gglgamfn(alpha+dble(n(ipop1,iloc,ial)))
     &           - gglgamfn(alpha+dble(n(ipop2,iloc,ial)))
     &           + gglgamfn(alpha+dble(ntmp(ipop1,iloc,ial)))
     &           + gglgamfn(alpha+dble(ntmp(ipop2,iloc,ial)))
            n1 = n1 + n(ipop1,iloc,ial)
            n2 = n2 + n(ipop2,iloc,ial)
            ntmp1 = ntmp1 + ntmp(ipop1,iloc,ial)
            ntmp2 = ntmp2 + ntmp(ipop2,iloc,ial)
c     write(*,*) 'lrpf=',lrpf
         enddo
         lrpf = lrpf + 
     &          gglgamfn(alpha*dble(nal(iloc))+dble(n1)) 
     &        + gglgamfn(alpha*dble(nal(iloc))+dble(n2)) 
     &        - gglgamfn(alpha*dble(nal(iloc))+dble(ntmp1)) 
     &        - gglgamfn(alpha*dble(nal(iloc))+dble(ntmp2))
c     write(*,*) 'lrpf=',lrpf
      enddo
c$$$      if(ipop1 .ne. ipop2) then
c$$$         do iloc=1,nloc
c$$$c         write(*,*) 'iloc=',iloc
c$$$            n1 = 0
c$$$            n2 = 0
c$$$            ntmp1 = 0
c$$$            ntmp2 = 0
c$$$            do ial=1,nal(iloc)
c$$$c     write(*,*) 'ial=',ial
c$$$               lrpf = lrpf + 
c$$$     &              n(ipop1,iloc,ial)*dlog(f(ipop1,iloc,ial)) 
c$$$     &            + n(ipop2,iloc,ial)*dlog(f(ipop2,iloc,ial))
c$$$     &             - ntmp(ipop1,iloc,ial)*dlog(ftmp(ipop1,iloc,ial)) 
c$$$     &             - ntmp(ipop2,iloc,ial)*dlog(ftmp(ipop2,iloc,ial)) 
c$$$     &              - gglgamfn(1+dble(n(ipop1,iloc,ial)))
c$$$     &              - gglgamfn(1+dble(n(ipop2,iloc,ial)))
c$$$     &              + gglgamfn(1+dble(ntmp(ipop1,iloc,ial)))
c$$$     &              + gglgamfn(1+dble(ntmp(ipop2,iloc,ial)))
c$$$               n1 = n1 + n(ipop1,iloc,ial)
c$$$               n2 = n2 + n(ipop2,iloc,ial)
c$$$               ntmp1 = ntmp1 + ntmp(ipop1,iloc,ial)
c$$$               ntmp2 = ntmp2 + ntmp(ipop2,iloc,ial)
c$$$c         write(*,*) 'lrpf=',lrpf
c$$$            enddo
c$$$            lrpf = lrpf + gglgamfn(dble(nal(iloc)+n1)) 
c$$$     &           + gglgamfn(dble(nal(iloc)+n2)) 
c$$$     &           - gglgamfn(dble(nal(iloc)+ntmp1)) 
c$$$     &           - gglgamfn(dble(nal(iloc)+ntmp2))
c$$$c         write(*,*) 'lrpf=',lrpf
c$$$         enddo
c$$$      else
c$$$         do iloc=1,nloc
c$$$c         write(*,*) 'iloc=',iloc
c$$$            n1 = 0
c$$$            ntmp1 = 0
c$$$            do ial=1,nal(iloc)
c$$$c     write(*,*) 'ial=',ial
c$$$               lrpf = lrpf + 
c$$$     &              n(ipop1,iloc,ial)*dlog(f(ipop1,iloc,ial)) 
c$$$     &             - ntmp(ipop1,iloc,ial)*dlog(ftmp(ipop1,iloc,ial)) 
c$$$     &             - gglgamfn(1+dble(n(ipop1,iloc,ial)))
c$$$     &             + gglgamfn(1+dble(ntmp(ipop1,iloc,ial)))
c$$$               n1 = n1 + n(ipop1,iloc,ial)
c$$$               ntmp1 = ntmp1 + ntmp(ipop1,iloc,ial)
c$$$c         write(*,*) 'lrpf=',lrpf
c$$$            enddo
c$$$            lrpf = lrpf + gglgamfn(dble(nal(iloc)+n1)) 
c$$$     &           - gglgamfn(dble(nal(iloc)+ntmp1)) 
c$$$c         write(*,*) 'lrpf=',lrpf
c$$$         enddo
c$$$      endif                     
      end 




***********************************************************************
*
*     log of ratio of proposals x priors in a joint update of c and f
*     (in subroutine udcf2)
*     
      double precision function lrppf(npopmax,nloc,nal,nalmax,n,
     &     ntmp,fa,drift,f,ftmp,ipop1,ipop2)
      implicit none 
      integer npopmax,nloc,nal,nalmax,n,ntmp,ipop1,ipop2
      double precision fa,drift,f,ftmp
      dimension nal(nloc),n(npopmax,nloc,nalmax),
     &     ntmp(npopmax,nloc,nalmax),f(npopmax,nloc,nalmax),
     &     fa(nloc,nalmax),ftmp(npopmax,nloc,nalmax),
     &     drift(npopmax)
      integer iloc,ial,n1,n2,ntmp1,ntmp2
      double precision gglgamfn,q1,q2
      lrppf = 0.d0
      q1 = drift(ipop1)/(1-drift(ipop1))
      q2 = drift(ipop2)/(1-drift(ipop2))
      
      do iloc=1,nloc
         n1 = 0
         n2 = 0
         ntmp1 = 0
         ntmp2 = 0
         do ial=1,nal(iloc)
            lrppf = lrppf + 
     &         n(ipop1,iloc,ial)*dlog(f(ipop1,iloc,ial)) 
     &       + n(ipop2,iloc,ial)*dlog(f(ipop2,iloc,ial) )
     &       - ntmp(ipop1,iloc,ial)*dlog(ftmp(ipop1,iloc,ial)) 
     &       - ntmp(ipop2,iloc,ial)*dlog(ftmp(ipop2,iloc,ial)) 
     &       + gglgamfn(fa(iloc,ial)*q1+
     &           dble(ntmp(ipop1,iloc,ial)))
     &       + gglgamfn(fa(iloc,ial)*q2+
     &           dble(ntmp(ipop2,iloc,ial)))
     &       - gglgamfn(fa(iloc,ial)*q1+
     &        dble(n(ipop1,iloc,ial)))
     &       - gglgamfn(fa(iloc,ial)*q2+
     &        dble(n(ipop2,iloc,ial)))
            n1 = n1 + n(ipop1,iloc,ial)
            n2 = n2 + n(ipop2,iloc,ial)
            ntmp1 = ntmp1 + ntmp(ipop1,iloc,ial)
            ntmp2 = ntmp2 + ntmp(ipop2,iloc,ial)
         enddo
         lrppf = lrppf + gglgamfn(q1+dble(n1)) 
     &           + gglgamfn(q2+dble(n2)) 
     &           - gglgamfn(q1+dble(ntmp1)) 
     &           - gglgamfn(q2+dble(ntmp2))
      enddo
      end function lrppf



************************************************************************
*     split and merge of populations with acceptance according to MH ratio
*     MH ratio does not depend on proposed frequencies
*     D-model
      subroutine smd(npop,npopmin,npopmax,f,fa,drift,
     &     nloc,nloc2,
     &     nal,nalmax,indcell,nindiv,npp,nppmax,c,ctmp,
     &     a,ptmp,ftmp,drifttmp,zz,cellpop,listcell,
     &     cellpophost,n,ntmp,ploidy)
      implicit none 
      integer npop,npopmin,npopmax,nloc,nal(nloc),
     &     nalmax,nindiv,npp,nppmax,indcell(nindiv),
     &     nloc2,c(nppmax),ctmp(nppmax),zz(nindiv,nloc2),
     &     n(npopmax,nloc,nalmax),ntmp(npopmax,nloc,nalmax),
     &     ploidy
      double precision f(npopmax,nloc,nalmax),drift(npopmax),
     &      ftmp(npopmax,nloc,nalmax),drifttmp(npopmax),
     &     a(nalmax),ptmp(nalmax),fa(nloc,nalmax)
      integer ipoprem,ipp,isplit,
     &     cellpop(nppmax),ncellpop,nu,listcell(nppmax),
     &     ipophost,ncellpophost,cellpophost(nppmax),ii
      double precision alpha,ggrunif,lbico,lratio,lTfd
      integer iloc
      double precision llr6,termf9bis,gglgamfn
      double precision bern,ggrbinom,b
      
c      write(*,*) 'debut smd' 
c      write(*,*) 'npop =',npop

      do ipp=1,nppmax
         cellpop(ipp) = -999
         listcell(ipp) = -999
      enddo
*     naissance ou mort ?
      b = ggrbinom(1.d0,0.5d0)
      if(b .eq. 1) then
         if(npop .lt. npopmax) then 

c            write(*,*) 'naissance'
            
*     split
*     choix de la pop qui split
            isplit = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
*     recherche des cellules affectees a cette pop
            call who(c,isplit,npp,nppmax,cellpop,ncellpop)
            if(ncellpop .gt. 0) then
*     tirage du nombre de cellules reallouees
C     ligne suivante corrigee le 17/01/08
C     nu = idint(dint(dble(ncellpop)*ggrunif(0.d0,1.d0)))
               nu = idint(dint(dble(ncellpop+1)*ggrunif(0.d0,1.d0)))
c               write(*,*) 'nu=',nu
               if(nu .gt. 0) then
*     tirage des cellules reallouees
                  call sample2(cellpop,nppmax,nu,ncellpop,listcell)
c                   write(*,*) 'listcell=',listcell
*     proposition de reallocation dans la pop npop+1
                  call split(npop+1,c,ctmp,nppmax,nu,listcell)
               else 
                  do ipp = 1,nppmax
                     ctmp(ipp) = c(ipp)
                  enddo
               endif
            else
               nu = 0
               do ipp = 1,nppmax
                  ctmp(ipp) = c(ipp)
               enddo
            endif
*     comptage des alleles sur chaque locus pour c puis ctmp
            call countn(nindiv,nloc,nloc2,npopmax,
     &           nppmax,nal,nalmax,zz,n,indcell,c,ploidy)
            call countn(nindiv,nloc,nloc2,npopmax,
     &           nppmax,nal,nalmax,zz,ntmp,indcell,ctmp,ploidy)
            
*     calcul du log du ratio
*     terme des proposal sur c
            lratio = dlog(2*dble(ncellpop+1)) + 
     &           lbico(ncellpop,nu) - dlog(dble(npop+1)) 
            
*     terme des priors sur c
            lratio = lratio + dble(npp)*(dlog(dble(npop)) - 
     &           dlog(dble(npop+1)))
            
c            write(*,*) 'in smd   c= ',c(1),c(2)
c            write(*,*) '         ctmp= ',ctmp(1),ctmp(2)
c            write(*,*) '       Rc=',lratio 
            
*     terme des frequences
            lratio = lratio +
     &           lTfd(isplit,ntmp,npopmax,nloc,nal,nalmax)  +
     &           lTfd(npop+1,ntmp,npopmax,nloc,nal,nalmax) -
     &           lTfd(isplit,n,npopmax,nloc,nal,nalmax) 

c            write(*,*) 'in smd Tf=',
c     &           lTfd(isplit,ntmp,npopmax,nloc,nal,nalmax)  
c     &           + lTfd(npop+1,ntmp,npopmax,nloc,nal,nalmax) 
c     &           - lTfd(isplit,n,npopmax,nloc,nal,nalmax) 
c            write(*,*) 'in smd lratio=',lratio 

c            write(*,*) 'lratio=',lratio
            lratio = dmin1(0.d0,lratio)
            alpha = dexp(lratio)
            bern = ggrbinom(1.d0,alpha)
             
            if(bern .eq. 1) then
*     proposition nouvelle freq et derive 
               call addfreq7bis(npop,npopmax,nloc,nloc,
     &              nal,nalmax,isplit,
     &              f,ftmp,fa,drift,drifttmp,a,ptmp,ntmp)
               call accept5(nppmax,npopmax,nloc,
     &              nalmax,nal,c,ctmp,f,ftmp,drift,drifttmp)
               npop = npop + 1
            endif
         endif
*     merge
      else
         if(npop .gt. npopmin) then 
c            write(*,*) 'mort'
*     tirage de la pop qui meurt
            ipoprem = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
            
*     tirage de la pop hote
            ipophost = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
            do while(ipophost .eq. ipoprem)
               ipophost = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
            enddo
            
*     on range dans la pop d'indice le plus petit
            if(ipophost .gt. ipoprem) then
               ii = ipophost
               ipophost = ipoprem
               ipoprem = ii
            endif
            
*     recherche des cellules qui vont etre reallouees
            call who(c,ipoprem,npp,nppmax,cellpop,ncellpop)
            
*     recherche des cellules de la pop hote
            call who(c,ipophost,npp,nppmax,cellpophost,
     &           ncellpophost)
            
*     proposition de reallocation dans la pop ipophost
            call merging(ipoprem,ipophost,c,ctmp,nppmax,
     &           ncellpop,cellpop)
            
*     comptage des alleles sur chaque locus pour c puis ctmp
            call countn(nindiv,nloc,nloc2,npopmax,
     &           nppmax,nal,nalmax,zz,n,indcell,c,ploidy)
            call countn(nindiv,nloc,nloc2,npopmax,
     &           nppmax,nal,nalmax,zz,ntmp,indcell,ctmp,ploidy)
            
*     calcul du log du ratio  
*     terme des proposal sur c
            lratio = dlog(dble(npop)) - 
     &           dlog(2*dble(ncellpop+ncellpophost+1)) -
     &           lbico(ncellpop+ncellpophost,ncellpop) 
*     terme des priors sur c
            lratio = lratio + 
     &           dble(npp)*(dlog(dble(npop)) - 
     &           dlog(dble(npop-1)))
c            write(*,*) 'terme en c lratio =',lratio
*     term des freq
            lratio = lratio + 
     &           lTfd(ipophost,ntmp,npopmax,nloc,nal,nalmax) -  
     &           lTfd(ipophost,n,npopmax,nloc,nal,nalmax) - 
     &           lTfd(ipoprem,n,npopmax,nloc,nal,nalmax) 
c            write(*,*) 'terme en f=',
c     &           lTfd(ipophost,ntmp,npopmax,nloc,nal,nalmax)  
c     &           - lTfd(ipophost,n,npopmax,nloc,nal,nalmax) 
c     &           - lTfd(ipoprem,n,npopmax,nloc,nal,nalmax) 
c            write(*,*) 'lratio =',lratio
             lratio = dmin1(0.d0,lratio)
             alpha  = dexp(lratio)
             bern = ggrbinom(1.d0,alpha)      
                          
             if(bern .eq. 1) then
*     propostion du nouveau tableau de freq et de derives
                call remfreq7bis(ipoprem,ipophost,
     &               npop,npopmax,nloc,nloc,nal,
     &               nalmax,f,ftmp,drift,drifttmp,fa,a,ptmp,ntmp)
                call accept5(nppmax,npopmax,nloc,
     &               nalmax,nal,c,ctmp,f,ftmp,drift,drifttmp)
                npop = npop - 1
             endif
          endif
       endif

c       write(*,*) 'fin smd' 
       end subroutine smd
************************************************************************     




c$$$
c$$$
c$$$************************************************************************
c$$$*     split and merge of populations 
c$$$*     D-model for allele frequencies
c$$$*     MH ratio does not depend on proposed frequencies
c$$$*     extended for quantitative variables 
c$$$      subroutine smdgq(npop,npopmin,npopmax,f,fa,drift,
c$$$     &     nloc,nloc2,
c$$$     &     nal,nalmax,indcell,nindiv,npp,nppmax,c,ctmp,
c$$$     &     a,ptmp,ftmp,drifttmp,zz,cellpop,listcell,
c$$$     &     cellpophost,n,ntmp,ploidy,qtc,nqtc,meanqtc,sdqtc,
c$$$     &     meanqtctmp,sdqtctmp,nnqtc,sqtc,ssqtc,
c$$$     &     ksiqtc,kappaqtc,alphaqtc,betaqtc,usegeno2,useqtc)
c$$$      implicit none 
c$$$      integer npop,npopmin,npopmax,nloc,nal(nloc),
c$$$     &     nalmax,nindiv,npp,nppmax,indcell(nindiv),
c$$$     &     nloc2,c(nppmax),ctmp(nppmax),zz(nindiv,nloc2),
c$$$     &     n(npopmax,nloc,nalmax),ntmp(npopmax,nloc,nalmax),
c$$$     &     ploidy,nqtc,nnqtc,usegeno2,useqtc
c$$$      double precision f(npopmax,nloc,nalmax),drift(npopmax),
c$$$     &      ftmp(npopmax,nloc,nalmax),drifttmp(npopmax),
c$$$     &     a(nalmax),ptmp(nalmax),fa(nloc,nalmax),qtc(nindiv,nqtc),
c$$$     &     meanqtc(npopmax,nqtc),sdqtc(npopmax,nqtc),
c$$$     &     meanqtctmp(npopmax,nqtc),sdqtctmp(npopmax,nqtc),
c$$$     &     sqtc(npopmax,nqtc),ssqtc(npopmax,nqtc),
c$$$     &     ksiqtc,kappaqtc,alphaqtc,betaqtc
c$$$      integer ipoprem,ipp,isplit,
c$$$     &     cellpop(nppmax),ncellpop,nu,listcell(nppmax),
c$$$     &     ipophost,ncellpophost,cellpophost(nppmax),ii
c$$$      double precision alpha,ggrunif,lbico,lratio,lTfd
c$$$      integer iloc,iqtc,ipop,ipoptmp,iindiv
c$$$      double precision llr6,termf9bis,gglgamfn,bern,ggrbinom,ggdnorm,
c$$$     &     ggdgamma,lratiogq,contriblr,junk,b
c$$$      
c$$$c      write(*,*) 'debut smd' 
c$$$c      write(*,*) 'npop =',npop
c$$$
c$$$      do ipp=1,nppmax
c$$$         cellpop(ipp) = -999
c$$$         listcell(ipp) = -999
c$$$      enddo
c$$$*     naissance ou mort ?
c$$$      b = ggrbinom(1.d0,0.5d0)
c$$$      if(b .eq. 1) then
c$$$         if(npop .lt. npopmax) then 
c$$$c     write(*,*) 'naissance'
c$$$c     write(*,*) 'npop=',npop           
c$$$*     split
c$$$*     choix de la pop qui split
c$$$            isplit = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
c$$$*     recherche des cellules affectees a cette pop
c$$$            call who(c,isplit,npp,nppmax,cellpop,ncellpop)
c$$$            if(ncellpop .gt. 0) then
c$$$*     tirage du nombre de cellules reallouees
c$$$C     ligne suivante corrigee le 17/01/08
c$$$C     nu = idint(dint(dble(ncellpop)*ggrunif(0.d0,1.d0)))
c$$$               nu = idint(dint(dble(ncellpop+1)*ggrunif(0.d0,1.d0)))
c$$$c     write(*,*) 'nu=',nu
c$$$               if(nu .gt. 0) then
c$$$*     tirage des cellules reallouees
c$$$                  call sample2(cellpop,nppmax,nu,ncellpop,listcell)
c$$$c                   write(*,*) 'listcell=',listcell
c$$$*     proposition de reallocation dans la pop npop+1
c$$$                  call split(npop+1,c,ctmp,nppmax,nu,listcell)
c$$$               else 
c$$$                  do ipp = 1,nppmax
c$$$                     ctmp(ipp) = c(ipp)
c$$$                  enddo
c$$$               endif
c$$$            else
c$$$               nu = 0
c$$$               do ipp = 1,nppmax
c$$$                  ctmp(ipp) = c(ipp)
c$$$               enddo
c$$$            endif
c$$$
c$$$            if(usegeno2 .eq. 1) then 
c$$$*     comptage des alleles sur chaque locus pour c puis ctmp
c$$$               call countn(nindiv,nloc,nloc2,npopmax,
c$$$     &              nppmax,nal,nalmax,zz,n,indcell,c,ploidy)
c$$$               call countn(nindiv,nloc,nloc2,npopmax,
c$$$     &              nppmax,nal,nalmax,zz,ntmp,indcell,ctmp,ploidy)
c$$$            endif
c$$$
c$$$*     calcul du log du ratio
c$$$*     terme des proposal sur c
c$$$            lratio = dlog(2*dble(ncellpop+1)) + 
c$$$     &           lbico(ncellpop,nu) - dlog(dble(npop+1)) 
c$$$            
c$$$*     terme des priors sur c
c$$$            lratio = lratio + dble(npp)*(dlog(dble(npop)) - 
c$$$     &           dlog(dble(npop+1)))
c$$$            
c$$$c            write(*,*) 'in smdgq   c= ',c(1),c(2),c(3),c(4)
c$$$c            write(*,*) '        ctmp= ',ctmp(1),ctmp(2),ctmp(3),ctmp(4)
c$$$c            write(*,*) '       Rc=',lratio 
c$$$            
c$$$            if(usegeno2 .eq. 1) then 
c$$$*     contribition of frequencies, this term includes likelihood ratio
c$$$               lratio = lratio +
c$$$     &              lTfd(isplit,ntmp,npopmax,nloc,nal,nalmax)  +
c$$$     &              lTfd(npop+1,ntmp,npopmax,nloc,nal,nalmax) -
c$$$     &              lTfd(isplit,n,npopmax,nloc,nal,nalmax) 
c$$$c     write(*,*) 'in smd Tf=',
c$$$c     &           lTfd(isplit,ntmp,npopmax,nloc,nal,nalmax)  
c$$$c     &           + lTfd(npop+1,ntmp,npopmax,nloc,nal,nalmax) 
c$$$c     &           - lTfd(isplit,n,npopmax,nloc,nal,nalmax) 
c$$$c     write(*,*) 'in smd lratio=',lratio 
c$$$            endif
c$$$            
c$$$            
c$$$            if(useqtc .eq. 1) then 
c$$$*     propose mean and variance of quantitative variables
c$$$               call propparqvsplit(qtc,nindiv,nqtc,indcell,c,ctmp,
c$$$     &              meanqtc,sdqtc,meanqtctmp,sdqtctmp,nnqtc,sqtc,
c$$$     &              ssqtc,npop,npopmax,nppmax,
c$$$     &              ksiqtc,kappaqtc,alphaqtc,betaqtc,isplit,contriblr)
c$$$*     contrib prior and proposal 
c$$$               lratio = lratio + contriblr
c$$$*     contrib likelihood quant. variable
c$$$*     argument usegeno2 is passed as 0 to ignore genotypes 
c$$$*     hence avoid having likelihood ratio of genotypes twice 
c$$$               junk = lratio
c$$$c      write(*,*) 'in smdgq lratiogq c=',c
c$$$c      write(*,*) 'in smdgq ctmp=',ctmp
c$$$               lratio = lratio + 
c$$$     &              lratiogq(zz,f,c,ctmp,indcell,indcell,
c$$$     &              npopmax,nloc,nalmax,nindiv,nloc,nloc2,
c$$$     &              nppmax,ploidy,qtc,nqtc,meanqtc,sdqtc,
c$$$     &              meanqtctmp,sdqtctmp,0,useqtc)
c$$$c               write(*,*) 'likelihood ratio=',lratio-junk
c$$$            endif
c$$$
c$$$c            write(*,*) 'lratio for spit=',lratio
c$$$            lratio = dmin1(0.d0,lratio)
c$$$            alpha = dexp(lratio)
c$$$            bern = ggrbinom(1.d0,alpha)
c$$$             
c$$$            if(bern .eq. 1) then
c$$$               if(usegeno2 .eq. 1) then 
c$$$*     proposition nouvelle freq et derive 
c$$$                  call addfreq7bis(npop,npopmax,nloc,nloc,
c$$$     &                 nal,nalmax,isplit,
c$$$     &                 f,ftmp,fa,drift,drifttmp,a,ptmp,ntmp)
c$$$               endif
c$$$               if(useqtc .eq. 1) then 
c$$$                  do iqtc = 1,nqtc
c$$$                     do ipop = 1,npopmax
c$$$                        meanqtc(ipop,iqtc) =  meanqtctmp(ipop,iqtc)
c$$$                        sdqtc(ipop,iqtc) =  sdqtctmp(ipop,iqtc)
c$$$                     enddo
c$$$                  enddo
c$$$               endif
c$$$               call accept5(nppmax,npopmax,nloc,
c$$$     &              nalmax,nal,c,ctmp,f,ftmp,drift,drifttmp)
c$$$               npop = npop + 1
c$$$            endif
c$$$         endif
c$$$*     merge
c$$$      else
c$$$         if(npop .gt. npopmin) then 
c$$$c      write(*,*) 'mort'
c$$$c      write(*,*) 'npop=',npop
c$$$*     tirage de la pop qui meurt
c$$$            ipoprem = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
c$$$            
c$$$*     tirage de la pop hote
c$$$            ipophost = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
c$$$            do while(ipophost .eq. ipoprem)
c$$$               ipophost = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
c$$$            enddo
c$$$            
c$$$*     on range dans la pop d'indice le plus petit
c$$$            if(ipophost .gt. ipoprem) then
c$$$               ii = ipophost
c$$$               ipophost = ipoprem
c$$$               ipoprem = ii
c$$$            endif
c$$$            
c$$$*     recherche des cellules qui vont etre reallouees
c$$$            call who(c,ipoprem,npp,nppmax,cellpop,ncellpop)
c$$$            
c$$$*     recherche des cellules de la pop hote
c$$$            call who(c,ipophost,npp,nppmax,cellpophost,
c$$$     &           ncellpophost)
c$$$            
c$$$*     proposition de reallocation dans la pop ipophost
c$$$            call merging(ipoprem,ipophost,c,ctmp,nppmax,
c$$$     &           ncellpop,cellpop)
c$$$            
c$$$            if(usegeno2 .eq. 1) then 
c$$$*     comptage des alleles sur chaque locus pour c puis ctmp
c$$$               call countn(nindiv,nloc,nloc2,npopmax,
c$$$     &              nppmax,nal,nalmax,zz,n,indcell,c,ploidy)
c$$$               call countn(nindiv,nloc,nloc2,npopmax,
c$$$     &              nppmax,nal,nalmax,zz,ntmp,indcell,ctmp,ploidy)
c$$$            endif
c$$$            
c$$$*     calcul du log du ratio  
c$$$*     terme des proposal sur c
c$$$            lratio = dlog(dble(npop)) - 
c$$$     &           dlog(2*dble(ncellpop+ncellpophost+1)) -
c$$$     &           lbico(ncellpop+ncellpophost,ncellpop) 
c$$$*     terme des priors sur c
c$$$            lratio = lratio + 
c$$$     &           dble(npp)*(dlog(dble(npop)) - 
c$$$     &           dlog(dble(npop-1)))
c$$$c     write(*,*) 'terme en c lratio =',lratio
c$$$            if(usegeno2 .eq. 1) then 
c$$$*     term des freq
c$$$               lratio = lratio + 
c$$$     &              lTfd(ipophost,ntmp,npopmax,nloc,nal,nalmax) -  
c$$$     &              lTfd(ipophost,n,npopmax,nloc,nal,nalmax) - 
c$$$     &              lTfd(ipoprem,n,npopmax,nloc,nal,nalmax) 
c$$$c     write(*,*) 'terme en f=',
c$$$c     &           lTfd(ipophost,ntmp,npopmax,nloc,nal,nalmax)  
c$$$c     &           - lTfd(ipophost,n,npopmax,nloc,nal,nalmax) 
c$$$c     &           - lTfd(ipoprem,n,npopmax,nloc,nal,nalmax) 
c$$$c     write(*,*) 'lratio =',lratio
c$$$            endif
c$$$            
c$$$            if(useqtc .eq. 1) then 
c$$$*     propose mean and variance of quantitative variables
c$$$               call propparqvmerge(qtc,nindiv,nqtc,indcell,c,ctmp,
c$$$     &              meanqtc,sdqtc,meanqtctmp,sdqtctmp,nnqtc,sqtc,
c$$$     &              ssqtc,npop,npopmax,nppmax,
c$$$     &              ksiqtc,kappaqtc,alphaqtc,betaqtc,ipophost,ipoprem,
c$$$     &              contriblr)
c$$$*     contrib prior and proposal 
c$$$               lratio = lratio + contriblr   
c$$$*     contrib likelihood quant. variable
c$$$*     argument usegeno2 is passed as 0 to ignore genotypes 
c$$$*     hence avoid having likelihood ratio of genotypes twice 
c$$$               lratio = lratio + 
c$$$     &              lratiogq(zz,f,c,ctmp,indcell,indcell,
c$$$     &              npopmax,nloc,nalmax,nindiv,nloc,nloc2,
c$$$     &              nppmax,ploidy,qtc,nqtc,meanqtc,sdqtc,
c$$$     &              meanqtctmp,sdqtctmp,0,useqtc)    
c$$$        
c$$$            endif
c$$$
c$$$           lratio = dmin1(0.d0,lratio)
c$$$           alpha  = dexp(lratio)
c$$$           bern = ggrbinom(1.d0,alpha)      
c$$$                          
c$$$           if(bern .eq. 1) then
c$$$              if(usegeno2 .eq. 1) then 
c$$$*     propostion du nouveau tableau de freq et de derives
c$$$                 call remfreq7bis(ipoprem,ipophost,
c$$$     &                npop,npopmax,nloc,nloc,nal,
c$$$     &                nalmax,f,ftmp,drift,drifttmp,fa,a,ptmp,ntmp)
c$$$              endif
c$$$              if(useqtc .eq. 1) then 
c$$$                 do iqtc = 1,nqtc
c$$$                    do ipop = 1,npopmax
c$$$                       meanqtc(ipop,iqtc) =  meanqtctmp(ipop,iqtc)
c$$$                       sdqtc(ipop,iqtc) =  sdqtctmp(ipop,iqtc)
c$$$                    enddo
c$$$                 enddo
c$$$              endif
c$$$              call accept5(nppmax,npopmax,nloc,
c$$$     &             nalmax,nal,c,ctmp,f,ftmp,drift,drifttmp)
c$$$              npop = npop - 1
c$$$           endif
c$$$        endif
c$$$      endif
c$$$
c$$$c       write(*,*) 'fin smd' 
c$$$      end subroutine smdgq
c$$$************************************************************************





************************************************************************
*     split and merge of populations 
*     uncorrelated model for frequencies
*     MH ratio does not depend on proposed frequencies
*     extended for quantitative variables 
      subroutine smdallvar3(npop,npopmin,npopmax,f,fa,drift,
     &     nlocd,nloch,ncolt,nql,
     &     nal,nalmax,indcell,nindiv,npp,nppmax,c,ctmp,
     &     a,ptmp,ftmp,drifttmp,yy,z,ql,cellpop,listcell,
     &     cellpophost,n,ntmp,qtc,nqtc,meanqtc,sdqtc,
     &     meanqtctmp,sdqtctmp,nnqtc,sqtc,ssqtc,
     &     ksiqtc,kappaqtc,alphaqtc,betaqtc,
     &     usegeno2,usegeno1,useql,useqtc,ploidy)
      implicit none 
      integer npop,npopmin,npopmax,ncolt,nal(ncolt),
     &     nalmax,nindiv,npp,nppmax,indcell(nindiv),nlocd,nlocd2,
     &     nloch,nql,c(nppmax),ctmp(nppmax),yy(nindiv,2*nlocd+2*nloch),
     &     z(nindiv,nloch),ql(nindiv,nql),
     &     n(npopmax,ncolt,nalmax),ntmp(npopmax,ncolt,nalmax),
     &     nqtc,nnqtc(npopmax,nqtc),usegeno2,usegeno1,useql,useqtc,
     &     ploidy
      double precision f(npopmax,ncolt,nalmax),drift(npopmax),
     &      ftmp(npopmax,ncolt,nalmax),drifttmp(npopmax),
     &     a(nalmax),ptmp(nalmax),fa(ncolt,nalmax),qtc(nindiv,nqtc),
     &     meanqtc(npopmax,nqtc),sdqtc(npopmax,nqtc),
     &     meanqtctmp(npopmax,nqtc),sdqtctmp(npopmax,nqtc),
     &     sqtc(npopmax,nqtc),ssqtc(npopmax,nqtc),
     &     ksiqtc,kappaqtc,alphaqtc,betaqtc
      integer ipoprem,ipp,isplit,
     &     cellpop(nppmax),ncellpop,nu,listcell(nppmax),
     &     ipophost,ncellpophost,cellpophost(nppmax),ii,i,j
      double precision alpha,ggrunif,lbico,lratio,lTfdallvar
      integer iloc,iqtc,ipop,ipoptmp,iindiv
      double precision llr6,termf9bis,gglgamfn,bern,ggrbinom,ggdnorm,
     &     ggdgamma,lrallvar2,contriblr,b,junk
      dimension ksiqtc(nqtc),kappaqtc(nqtc),alphaqtc(nqtc),betaqtc(nqtc)
      
c      write(*,*) 'debut smdallvar3' 
c      write(*,*) 'npop =',npop


      do ipp=1,nppmax
         cellpop(ipp) = -999
         listcell(ipp) = -999
      enddo
*     naissance ou mort ?
      b = ggrbinom(1.d0,0.5d0)
c      write(*,*) 'b=',b
      if(b .eq. 1) then
         if(npop .lt. npopmax) then 
c$$$            write(*,*) 'naissance'
c$$$            write(*,*) 'npop=',npop       
c$$$             write(*,*) 'npp=',npp       
*     split
*     choix de la pop qui split
            isplit = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
*     recherche des cellules affectees a cette pop
            call who(c,isplit,npp,nppmax,cellpop,ncellpop)
            if(ncellpop .gt. 0) then
*     tirage du nombre de cellules reallouees
               nu = idint(dint(dble(ncellpop+1)*ggrunif(0.d0,1.d0)))
c               write(*,*) 'nu=',nu
               if(nu .gt. 0) then
*     tirage des cellules reallouees
                  call sample2(cellpop,nppmax,nu,ncellpop,listcell)
c                   write(*,*) 'listcell=',listcell
*     proposition de reallocation dans la pop npop+1
                  call split(npop+1,c,ctmp,nppmax,nu,listcell)
               else 
                  do ipp = 1,nppmax
                     ctmp(ipp) = c(ipp)
                  enddo
               endif
            else
               nu = 0
               do ipp = 1,nppmax
                  ctmp(ipp) = c(ipp)
               enddo
            endif
c$$$            write(*,*) 'npp=',npp
c$$$            write(*,*) 'isplit=',isplit
c$$$            write(*,*) 'c=',c(1),c(2),c(3),c(4)
c$$$            write(*,*) 'ctmp=',ctmp(1),ctmp(2),ctmp(3),ctmp(4)

            if(((usegeno2 .eq. 1) .or. (usegeno1 .eq. 1)) .or. 
     &           useql .eq. 1) then
*     comptage des alleles sur chaque locus pour c puis ctmp
               call countnallvar2(nindiv,nlocd,nloch,nql,ncolt,
     &              nppmax,nal,nalmax,yy,z,ql,n,indcell,c,npopmax,
     &              usegeno2,usegeno1,useql,ploidy)
               call countnallvar2(nindiv,nlocd,nloch,nql,ncolt,
     &              nppmax,nal,nalmax,yy,z,ql,ntmp,indcell,ctmp,npopmax,
     &              usegeno2,usegeno1,useql,ploidy)
            endif

*     calcul du log du ratio
*     terme des proposal sur c
            lratio = dlog(2*dble(ncellpop+1)) + 
     &           lbico(ncellpop,nu) - dlog(dble(npop+1)) 
c            write(6,*) 'split apres proposal c lratio=',lratio 
            
            
*     terme des priors sur c
            lratio = lratio + dble(npp)*(dlog(dble(npop)) - 
     &           dlog(dble(npop+1)))
c            write(6,*) 'split apres prior c lratio=',lratio 
c            write(*,*) 'in smdgq   c= ',c(1),c(2),c(3),c(4)
c            write(*,*) '        ctmp= ',ctmp(1),ctmp(2),ctmp(3),ctmp(4)
c            write(*,*) '       Rc=',lratio 
            if(((usegeno2 .eq. 1) .or. (usegeno1 .eq. 1)) .or. 
     &           useql .eq. 1) then
c               write(*,*) 'avant appel lTfdallvar'
*     contribition of frequencies, this term includes likelihood ratio
c$$$               lratio = lratio +
c$$$     &              lTfd(isplit,ntmp,npopmax,nloc,nal,nalmax)  +
c$$$     &              lTfd(npop+1,ntmp,npopmax,nloc,nal,nalmax) -
c$$$     &              lTfd(isplit,n,npopmax,nloc,nal,nalmax) 
               lratio = lratio + 
     &              lTfdallvar(isplit,ntmp,npopmax,nlocd,nloch,
     &              nql,ncolt,nal,nalmax,usegeno2,usegeno1,useql) + 
     &              lTfdallvar(npop+1,ntmp,npopmax,nlocd,nloch,
     &              nql,ncolt,nal,nalmax,usegeno2,usegeno1,useql) - 
     &              lTfdallvar(isplit,n,npopmax,nlocd,nloch,
     &              nql,ncolt,nal,nalmax,usegeno2,usegeno1,useql)
c               write(*,*) 'apres appel lTfdallvar'
c               write(6,*) 'split apres prior contrib freq =',lratio 
            endif
            
            
            if(useqtc .eq. 1) then 
*     propose mean and variance of quantitative variables
               call propqvsplit3(qtc,nindiv,nqtc,indcell,c,ctmp,
     &              meanqtc,sdqtc,meanqtctmp,sdqtctmp,nnqtc,sqtc,
     &              ssqtc,npop,npopmax,nppmax,
     &              ksiqtc,kappaqtc,alphaqtc,betaqtc,isplit,contriblr)
*     contrib prior and proposal 
               lratio = lratio + contriblr
c$$$               write(6,*) 'smdallvar3  apres propqvsplit3=',lratio 
c$$$               do i=1,npop
c$$$                  write(*,*) 'meanqtc=',(meanqtc(i,j),j=1,nqtc)
c$$$               enddo
c$$$               do i=1,npop
c$$$                  write(*,*) 'sdqtc=',(sdqtc(i,j),j=1,nqtc)
c$$$               enddo
c$$$               do i=1,npop
c$$$                  write(*,*) 'meanqtctmp=',(meanqtctmp(i,j),j=1,nqtc)
c$$$               enddo
c$$$               do i=1,npop
c$$$                  write(*,*) 'sdqtctmp=',(meanqtctmp(i,j),j=1,nqtc)
c$$$               enddo
*     contrib likelihood quant. variable
*     argument usegeno2 is passed as 0 to ignore genotypes 
*     hence avoid having likelihood ratio of genotypes twice 
c               junk = lratio
c      write(*,*) 'in smdgq lratiogq c=',c
c      write(*,*) 'in smdgq ctmp=',ctmp
c$$$               lratio = lratio + 
c$$$     &              lratiogq(yy,f,c,ctmp,indcell,indcell,
c$$$     &              npopmax,nloc,nalmax,nindiv,nloc,nloc2,
c$$$     &              nppmax,ploidy,qtc,nqtc,meanqtc,sdqtc,
c$$$     &              meanqtctmp,sdqtctmp,0,useqtc)
               lratio = lratio + lrallvar2(yy,z,ql,
     &              f,c,ctmp,indcell,indcell,
     &              npopmax,nlocd,nloch,nql,ncolt,nalmax,nindiv,
     &              nppmax,qtc,nqtc,meanqtc,sdqtc,meanqtctmp,sdqtctmp,
     &              0,0,0,useqtc,ploidy)
c               write(6,*) 'in smdallvar3 split apres llhood=',lratio 
            endif

            lratio = dmin1(0.d0,lratio)
            alpha = dexp(lratio)
            bern = ggrbinom(1.d0,alpha)
             
            if(bern .eq. 1) then
               if(((usegeno2 .eq. 1) .or. (usegeno1 .eq. 1)) .or. 
     &              useql .eq. 1) then
*     proposition nouvelle freq et derive 
                  call addfallbis(npop,npopmax,nlocd,nlocd2,nloch,nql,
     &                 ncolt,nal,nalmax,isplit,f,ftmp,fa,drift,drifttmp,
     &                 a,ptmp,ntmp,usegeno2,usegeno1,useql)
               endif
               if(useqtc .eq. 1) then 
                  do iqtc = 1,nqtc
                     do ipop = 1,npopmax
                        meanqtc(ipop,iqtc) =  meanqtctmp(ipop,iqtc)
                        sdqtc(ipop,iqtc) =  sdqtctmp(ipop,iqtc)
                     enddo
                  enddo
               endif
               call accept5(nppmax,npopmax,ncolt,
     &              nalmax,nal,c,ctmp,f,ftmp,drift,drifttmp)
               npop = npop + 1
            endif
         endif
*     merge
      else
         if(npop .gt. npopmin) then 
c      write(*,*) 'mort'
c      write(*,*) 'npop=',npop
*     tirage de la pop qui meurt
            ipoprem = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
            
*     tirage de la pop hote
            ipophost = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
            do while(ipophost .eq. ipoprem)
               ipophost = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
            enddo
            
*     on range dans la pop d'indice le plus petit
            if(ipophost .gt. ipoprem) then
               ii = ipophost
               ipophost = ipoprem
               ipoprem = ii
            endif
            
*     recherche des cellules qui vont etre reallouees
            call who(c,ipoprem,npp,nppmax,cellpop,ncellpop)
            
*     recherche des cellules de la pop hote
            call who(c,ipophost,npp,nppmax,cellpophost,
     &           ncellpophost)
            
*     proposition de reallocation dans la pop ipophost
            call merging(ipoprem,ipophost,c,ctmp,nppmax,
     &           ncellpop,cellpop)
c$$$            write(*,*) 'npp=',npp
c$$$            write(*,*) 'ipoprem=',ipoprem
c$$$            write(*,*) 'ipophost=',ipophost
c$$$            write(*,*) 'c=',c(1),c(2),c(3),c(4)
c$$$            write(*,*) 'ctmp=',ctmp(1),ctmp(2),ctmp(3),ctmp(4)   
            
            if(((usegeno2 .eq. 1) .or. (usegeno1 .eq. 1)) .or. 
     &           useql .eq. 1) then
*     comptage des alleles sur chaque locus pour c puis ctmp
c     write(*,*) 'avant appel countnallvar'
               call countnallvar2(nindiv,nlocd,nloch,nql,ncolt,
     &              nppmax,nal,nalmax,yy,z,ql,n,indcell,c,npopmax,
     &              usegeno2, usegeno1,useql,ploidy)
               call countnallvar2(nindiv,nlocd,nloch,nql,ncolt,
     &              nppmax,nal,nalmax,yy,z,ql,ntmp,indcell,ctmp,npopmax,
     &              usegeno2, usegeno1,useql,ploidy)
c     write(*,*) 'apres appel countnallvar'
            endif
        

*     calcul du log du ratio  
*     terme des proposal sur c
            lratio = dlog(dble(npop)) - 
     &           dlog(2*dble(ncellpop+ncellpophost+1)) -
     &           lbico(ncellpop+ncellpophost,ncellpop) 
c            write(6,*) 'merge apres proposal c lratio=',lratio 
*     terme des priors sur c
            lratio = lratio + 
     &           dble(npp)*(dlog(dble(npop)) - 
     &           dlog(dble(npop-1)))
c            write(6,*) 'merge apres prior c lratio=',lratio 
c     write(*,*) 'terme en c lratio =',lratio
            if(((usegeno2 .eq. 1) .or. (usegeno1 .eq. 1)) .or. 
     &           useql .eq. 1) then
*     term des freq
c$$$  lratio = lratio + 
c$$$     &              lTfd(ipophost,ntmp,npopmax,nloc,nal,nalmax) -  
c$$$  &              lTfd(ipophost,n,npopmax,nloc,nal,nalmax) - 
c$$$  &              lTfd(ipoprem,n,npopmax,nloc,nal,nalmax) 
c               write(*,*) 'avant appel lTfdallvar'
               lratio = lratio + 
     &              lTfdallvar(ipophost,ntmp,npopmax,nlocd,nloch,
     &              nql,ncolt,nal,nalmax,usegeno2,usegeno1,useql) -
     &              lTfdallvar(ipophost,n,npopmax,nlocd,nloch,
     &              nql,ncolt,nal,nalmax,usegeno2,usegeno1,useql) - 
     &              lTfdallvar(ipoprem,n,npopmax,nlocd,nloch,
     &              nql,ncolt,nal,nalmax,usegeno2,usegeno1,useql)
c               write(*,*) 'apres appel lTfdallvar'
c               write(6,*) 'merge apres prior contrib freq =',lratio 
            endif
            
            if(useqtc .eq. 1) then 
*     propose mean and variance of quantitative variables
               call propqvmrg3(qtc,nindiv,nqtc,indcell,c,ctmp,
     &              meanqtc,sdqtc,meanqtctmp,sdqtctmp,nnqtc,sqtc,
     &              ssqtc,npop,npopmax,nppmax,
     &              ksiqtc,kappaqtc,alphaqtc,betaqtc,ipophost,ipoprem,
     &              contriblr)
*     contrib prior and proposal 
               lratio = lratio + contriblr   
c$$$               write(6,*) 'smdallvar3  apres propqvmrg3=',lratio 
c$$$               do i=1,npop
c$$$                  write(*,*) 'meanqtc=',(meanqtc(i,j),j=1,nqtc)
c$$$               enddo
c$$$               do i=1,npop
c$$$                  write(*,*) 'sdqtc=',(sdqtc(i,j),j=1,nqtc)
c$$$               enddo
c$$$               do i=1,npop
c$$$                  write(*,*) 'meanqtctmp=',(meanqtctmp(i,j),j=1,nqtc)
c$$$               enddo
c$$$               do i=1,npop
c$$$                  write(*,*) 'sdqtctmp=',(meanqtctmp(i,j),j=1,nqtc)
c$$$               enddo
*     contrib likelihood quant. variable
*     argument usegeno2 is passed as 0 to ignore genotypes 
*     hence avoid having likelihood ratio of genotypes twice 
c$$$  lratio = lratio + 
c$$$  &              lratiogq(yy,f,c,ctmp,indcell,indcell,
c$$$     &              npopmax,nloc,nalmax,nindiv,nloc,nloc2,
c$$$  &              nppmax,ploidy,qtc,nqtc,meanqtc,sdqtc,
c$$$  &              meanqtctmp,sdqtctmp,0,useqtc)    
               lratio = lratio + lrallvar2(yy,z,ql,
     &              f,c,ctmp,indcell,indcell,
     &              npopmax,nlocd,nloch,nql,ncolt,nalmax,nindiv,
     &              nppmax,qtc,nqtc,meanqtc,sdqtc,meanqtctmp,sdqtctmp,
     &              0,0,0,useqtc,ploidy)   
c               write(6,*) 'in smdallvar3 merge apres llhood=',lratio 
            endif

            lratio = dmin1(0.d0,lratio)
            alpha  = dexp(lratio)
            bern = ggrbinom(1.d0,alpha)      
            
            if(bern .eq. 1) then
               if(((usegeno2 .eq. 1) .or. (usegeno1 .eq. 1)) .or. 
     &              useql .eq. 1) then
*     propostion du nouveau tableau de freq et de derives
c$$$  call remfreq7bis(ipoprem,ipophost,
c$$$  &                  npop,npopmax,nloc,nal,
c$$$  &                  nalmax,f,ftmp,drift,drifttmp,fa,a,ptmp,ntmp)
                  call  remfallbis(ipoprem,ipophost,
     &                 npop,npopmax,nlocd,nlocd2,nloch,nql,ncolt,nal,
     &                 nalmax,f,ftmp,drift,drifttmp,fa,a,ptmp,ntmp,
     &                 usegeno2,usegeno1,useql)
               endif
               if(useqtc .eq. 1) then 
                  do iqtc = 1,nqtc
                     do ipop = 1,npopmax
                        meanqtc(ipop,iqtc) =  meanqtctmp(ipop,iqtc)
                        sdqtc(ipop,iqtc) =  sdqtctmp(ipop,iqtc)
                     enddo
                  enddo
               endif
               call accept5(nppmax,npopmax,ncolt,
     &              nalmax,nal,c,ctmp,f,ftmp,drift,drifttmp)
               npop = npop - 1
            endif
         endif
      endif

c$$$      write(*,*) 'z =', z
c$$$      write(*,*) 'yy =', yy
c$$$      write(*,*) 'ql =', ql
c$$$      write(*,*) 'c =', c
c$$$      write(*,*) 'ctmp =',ctmp 
c$$$      write(*,*) 'f =', f  
c$$$      write(*,*) 'ftmp =', ftmp
c$$$      write(*,*) 'qtc=',qtc
c$$$      write(*,*) 'meanqtc=',meanqtc
c$$$      write(*,*) 'sdqtc=',sdqtc
c$$$      write(*,*) 'nnqtc=',nnqtc
c       write(*,*) 'fin smdallvar3' 
      end subroutine smdallvar3
************************************************************************


c$$$
c$$$************************************************************************
c$$$*     split and merge of populations 
c$$$*     uncorrelated model for frequencies
c$$$*     MH ratio does not depend on proposed frequencies
c$$$*     extended for quantitative variables 
c$$$      subroutine smdallvar2(npop,npopmin,npopmax,f,fa,drift,
c$$$     &     nlocd,nloch,ncolt,nql,
c$$$     &     nal,nalmax,indcell,nindiv,npp,nppmax,c,ctmp,
c$$$     &     a,ptmp,ftmp,drifttmp,yy,z,ql,cellpop,listcell,
c$$$     &     cellpophost,n,ntmp,qtc,nqtc,meanqtc,sdqtc,
c$$$     &     meanqtctmp,sdqtctmp,nnqtc,sqtc,ssqtc,
c$$$     &     ksiqtc,kappaqtc,alphaqtc,betaqtc,
c$$$     &     usegeno2,usegeno1,useql,useqtc,ploidy)
c$$$      implicit none 
c$$$      integer npop,npopmin,npopmax,ncolt,nal(ncolt),
c$$$     &     nalmax,nindiv,npp,nppmax,indcell(nindiv),nlocd,nlocd2,
c$$$     &     nloch,nql,c(nppmax),ctmp(nppmax),yy(nindiv,2*nlocd+2*nloch),
c$$$     &     z(nindiv,nloch),ql(nindiv,nql),
c$$$     &     n(npopmax,ncolt,nalmax),ntmp(npopmax,ncolt,nalmax),
c$$$     &     nqtc,nnqtc,usegeno2,usegeno1,useql,useqtc,ploidy
c$$$      double precision f(npopmax,ncolt,nalmax),drift(npopmax),
c$$$     &      ftmp(npopmax,ncolt,nalmax),drifttmp(npopmax),
c$$$     &     a(nalmax),ptmp(nalmax),fa(ncolt,nalmax),qtc(nindiv,nqtc),
c$$$     &     meanqtc(npopmax,nqtc),sdqtc(npopmax,nqtc),
c$$$     &     meanqtctmp(npopmax,nqtc),sdqtctmp(npopmax,nqtc),
c$$$     &     sqtc(npopmax,nqtc),ssqtc(npopmax,nqtc),
c$$$     &     ksiqtc,kappaqtc,alphaqtc,betaqtc
c$$$      integer ipoprem,ipp,isplit,
c$$$     &     cellpop(nppmax),ncellpop,nu,listcell(nppmax),
c$$$     &     ipophost,ncellpophost,cellpophost(nppmax),ii
c$$$      double precision alpha,ggrunif,lbico,lratio,lTfdallvar
c$$$      integer b,iloc,iqtc,ipop,ipoptmp,iindiv
c$$$      double precision llr6,termf9bis,gglgamfn,bern,ggrbinom,ggdnorm,
c$$$     &     ggdgamma,lrallvar2,contriblr,junk
c$$$      dimension ksiqtc(nqtc),kappaqtc(nqtc),alphaqtc(nqtc),betaqtc(nqtc)
c$$$      
c$$$c      write(*,*) 'debut smdallvar2' 
c$$$c      write(*,*) 'npop =',npop
c$$$
c$$$
c$$$      do ipp=1,nppmax
c$$$         cellpop(ipp) = -999
c$$$         listcell(ipp) = -999
c$$$      enddo
c$$$*     naissance ou mort ?
c$$$      b = ggrbinom(1.d0,0.5d0)
c$$$      if(b .eq. 1) then
c$$$         if(npop .lt. npopmax) then 
c$$$c            write(*,*) 'naissance'
c$$$c            write(*,*) 'npop=',npop       
c$$$c             write(*,*) 'npp=',npp       
c$$$*     split
c$$$*     choix de la pop qui split
c$$$            isplit = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
c$$$*     recherche des cellules affectees a cette pop
c$$$            call who(c,isplit,npp,nppmax,cellpop,ncellpop)
c$$$            if(ncellpop .gt. 0) then
c$$$*     tirage du nombre de cellules reallouees
c$$$               nu = idint(dint(dble(ncellpop+1)*ggrunif(0.d0,1.d0)))
c$$$c               write(*,*) 'nu=',nu
c$$$               if(nu .gt. 0) then
c$$$*     tirage des cellules reallouees
c$$$                  call sample2(cellpop,nppmax,nu,ncellpop,listcell)
c$$$c                   write(*,*) 'listcell=',listcell
c$$$*     proposition de reallocation dans la pop npop+1
c$$$                  call split(npop+1,c,ctmp,nppmax,nu,listcell)
c$$$               else 
c$$$                  do ipp = 1,nppmax
c$$$                     ctmp(ipp) = c(ipp)
c$$$                  enddo
c$$$               endif
c$$$            else
c$$$               nu = 0
c$$$               do ipp = 1,nppmax
c$$$                  ctmp(ipp) = c(ipp)
c$$$               enddo
c$$$            endif
c$$$c            write(*,*) 'c=',c(1),c(2),c(3),c(4)
c$$$c            write(*,*) 'ctmp=',ctmp(1),ctmp(2),ctmp(3),ctmp(4)
c$$$
c$$$            if(((usegeno2 .eq. 1) .or. (usegeno1 .eq. 1)) .or. 
c$$$     &           useql .eq. 1) then
c$$$*     comptage des alleles sur chaque locus pour c puis ctmp
c$$$               call countnallvar2(nindiv,nlocd,nloch,nql,ncolt,
c$$$     &              nppmax,nal,nalmax,yy,z,ql,n,indcell,c,npopmax,
c$$$     &              usegeno2,usegeno1,useql,ploidy)
c$$$               call countnallvar2(nindiv,nlocd,nloch,nql,ncolt,
c$$$     &              nppmax,nal,nalmax,yy,z,ql,ntmp,indcell,ctmp,npopmax,
c$$$     &              usegeno2,usegeno1,useql,ploidy)
c$$$            endif
c$$$
c$$$*     calcul du log du ratio
c$$$*     terme des proposal sur c
c$$$            lratio = dlog(2*dble(ncellpop+1)) + 
c$$$     &           lbico(ncellpop,nu) - dlog(dble(npop+1)) 
c$$$c            write(6,*) 'split apres proposal c lratio=',lratio 
c$$$            
c$$$            
c$$$*     terme des priors sur c
c$$$            lratio = lratio + dble(npp)*(dlog(dble(npop)) - 
c$$$     &           dlog(dble(npop+1)))
c$$$c            write(6,*) 'split apres prior c lratio=',lratio 
c$$$c            write(*,*) 'in smdgq   c= ',c(1),c(2),c(3),c(4)
c$$$c            write(*,*) '        ctmp= ',ctmp(1),ctmp(2),ctmp(3),ctmp(4)
c$$$c            write(*,*) '       Rc=',lratio 
c$$$            if(((usegeno2 .eq. 1) .or. (usegeno1 .eq. 1)) .or. 
c$$$     &           useql .eq. 1) then
c$$$c               write(*,*) 'avant appel lTfdallvar'
c$$$*     contribition of frequencies, this term includes likelihood ratio
c$$$c$$$               lratio = lratio +
c$$$c$$$     &              lTfd(isplit,ntmp,npopmax,nloc,nal,nalmax)  +
c$$$c$$$     &              lTfd(npop+1,ntmp,npopmax,nloc,nal,nalmax) -
c$$$c$$$     &              lTfd(isplit,n,npopmax,nloc,nal,nalmax) 
c$$$               lratio = lratio + 
c$$$     &              lTfdallvar(isplit,ntmp,npopmax,nlocd,nloch,
c$$$     &              nql,ncolt,nal,nalmax,usegeno2,usegeno1,useql) + 
c$$$     &              lTfdallvar(npop+1,ntmp,npopmax,nlocd,nloch,
c$$$     &              nql,ncolt,nal,nalmax,usegeno2,usegeno1,useql) - 
c$$$     &              lTfdallvar(isplit,n,npopmax,nlocd,nloch,
c$$$     &              nql,ncolt,nal,nalmax,usegeno2,usegeno1,useql)
c$$$c               write(*,*) 'apres appel lTfdallvar'
c$$$c               write(6,*) 'split apres prior contrib freq =',lratio 
c$$$            endif
c$$$            
c$$$            
c$$$            if(useqtc .eq. 1) then 
c$$$*     propose mean and variance of quantitative variables
c$$$               call propparqvsplit(qtc,nindiv,nqtc,indcell,c,ctmp,
c$$$     &              meanqtc,sdqtc,meanqtctmp,sdqtctmp,nnqtc,sqtc,
c$$$     &              ssqtc,npop,npopmax,nppmax,
c$$$     &              ksiqtc,kappaqtc,alphaqtc,betaqtc,isplit,contriblr)
c$$$*     contrib prior and proposal 
c$$$               lratio = lratio + contriblr
c$$$*     contrib likelihood quant. variable
c$$$*     argument usegeno2 is passed as 0 to ignore genotypes 
c$$$*     hence avoid having likelihood ratio of genotypes twice 
c$$$c               junk = lratio
c$$$c      write(*,*) 'in smdgq lratiogq c=',c
c$$$c      write(*,*) 'in smdgq ctmp=',ctmp
c$$$c$$$               lratio = lratio + 
c$$$c$$$     &              lratiogq(yy,f,c,ctmp,indcell,indcell,
c$$$c$$$     &              npopmax,nloc,nalmax,nindiv,nloc,nloc2,
c$$$c$$$     &              nppmax,ploidy,qtc,nqtc,meanqtc,sdqtc,
c$$$c$$$     &              meanqtctmp,sdqtctmp,0,useqtc)
c$$$               lratio = lratio + lrallvar2(yy,z,ql,
c$$$     &              f,c,ctmp,indcell,indcell,
c$$$     &              npopmax,nlocd,nloch,nql,ncolt,nalmax,nindiv,
c$$$     &              nppmax,qtc,nqtc,meanqtc,sdqtc,meanqtctmp,sdqtctmp,
c$$$     &              0,0,0,useqtc,ploidy)
c$$$c               write(6,*) 'split apres llhood=',lratio 
c$$$            endif
c$$$
c$$$c            write(*,*) 'lratio for split=',lratio
c$$$            lratio = dmin1(0.d0,lratio)
c$$$            alpha = dexp(lratio)
c$$$            bern = ggrbinom(1.d0,alpha)
c$$$             
c$$$            if(bern .eq. 1) then
c$$$               if(((usegeno2 .eq. 1) .or. (usegeno1 .eq. 1)) .or. 
c$$$     &              useql .eq. 1) then
c$$$*     proposition nouvelle freq et derive 
c$$$                  call addfallbis(npop,npopmax,nlocd,nlocd2,nloch,nql,
c$$$     &                 ncolt,nal,nalmax,isplit,f,ftmp,fa,drift,drifttmp,
c$$$     &                 a,ptmp,ntmp,usegeno2,usegeno1,useql)
c$$$               endif
c$$$               if(useqtc .eq. 1) then 
c$$$                  do iqtc = 1,nqtc
c$$$                     do ipop = 1,npopmax
c$$$                        meanqtc(ipop,iqtc) =  meanqtctmp(ipop,iqtc)
c$$$                        sdqtc(ipop,iqtc) =  sdqtctmp(ipop,iqtc)
c$$$                     enddo
c$$$                  enddo
c$$$               endif
c$$$               call accept5(nppmax,npopmax,ncolt,
c$$$     &              nalmax,nal,c,ctmp,f,ftmp,drift,drifttmp)
c$$$               npop = npop + 1
c$$$            endif
c$$$         endif
c$$$*     merge
c$$$      else
c$$$         if(npop .gt. npopmin) then 
c$$$c      write(*,*) 'mort'
c$$$c      write(*,*) 'npop=',npop
c$$$*     tirage de la pop qui meurt
c$$$            ipoprem = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
c$$$            
c$$$*     tirage de la pop hote
c$$$            ipophost = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
c$$$            do while(ipophost .eq. ipoprem)
c$$$               ipophost = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
c$$$            enddo
c$$$            
c$$$*     on range dans la pop d'indice le plus petit
c$$$            if(ipophost .gt. ipoprem) then
c$$$               ii = ipophost
c$$$               ipophost = ipoprem
c$$$               ipoprem = ii
c$$$            endif
c$$$            
c$$$*     recherche des cellules qui vont etre reallouees
c$$$            call who(c,ipoprem,npp,nppmax,cellpop,ncellpop)
c$$$            
c$$$*     recherche des cellules de la pop hote
c$$$            call who(c,ipophost,npp,nppmax,cellpophost,
c$$$     &           ncellpophost)
c$$$            
c$$$*     proposition de reallocation dans la pop ipophost
c$$$            call merging(ipoprem,ipophost,c,ctmp,nppmax,
c$$$     &           ncellpop,cellpop)
c$$$c            write(*,*) 'c=',c(1),c(2),c(3),c(4)
c$$$c            write(*,*) 'ctmp=',ctmp(1),ctmp(2),ctmp(3),ctmp(4)   
c$$$            
c$$$            if(((usegeno2 .eq. 1) .or. (usegeno1 .eq. 1)) .or. 
c$$$     &           useql .eq. 1) then
c$$$*     comptage des alleles sur chaque locus pour c puis ctmp
c$$$c     write(*,*) 'avant appel countnallvar'
c$$$               call countnallvar2(nindiv,nlocd,nloch,nql,ncolt,
c$$$     &              nppmax,nal,nalmax,yy,z,ql,n,indcell,c,npopmax,
c$$$     &              usegeno2, usegeno1,useql,ploidy)
c$$$               call countnallvar2(nindiv,nlocd,nloch,nql,ncolt,
c$$$     &              nppmax,nal,nalmax,yy,z,ql,ntmp,indcell,ctmp,npopmax,
c$$$     &              usegeno2, usegeno1,useql,ploidy)
c$$$c     write(*,*) 'apres appel countnallvar'
c$$$            endif
c$$$        
c$$$
c$$$*     calcul du log du ratio  
c$$$*     terme des proposal sur c
c$$$            lratio = dlog(dble(npop)) - 
c$$$     &           dlog(2*dble(ncellpop+ncellpophost+1)) -
c$$$     &           lbico(ncellpop+ncellpophost,ncellpop) 
c$$$c            write(6,*) 'merge apres proposal c lratio=',lratio 
c$$$*     terme des priors sur c
c$$$            lratio = lratio + 
c$$$     &           dble(npp)*(dlog(dble(npop)) - 
c$$$     &           dlog(dble(npop-1)))
c$$$c            write(6,*) 'merge apres prior c lratio=',lratio 
c$$$c     write(*,*) 'terme en c lratio =',lratio
c$$$            if(((usegeno2 .eq. 1) .or. (usegeno1 .eq. 1)) .or. 
c$$$     &           useql .eq. 1) then
c$$$*     term des freq
c$$$c$$$  lratio = lratio + 
c$$$c$$$     &              lTfd(ipophost,ntmp,npopmax,nloc,nal,nalmax) -  
c$$$c$$$  &              lTfd(ipophost,n,npopmax,nloc,nal,nalmax) - 
c$$$c$$$  &              lTfd(ipoprem,n,npopmax,nloc,nal,nalmax) 
c$$$c               write(*,*) 'avant appel lTfdallvar'
c$$$               lratio = lratio + 
c$$$     &              lTfdallvar(ipophost,ntmp,npopmax,nlocd,nloch,
c$$$     &              nql,ncolt,nal,nalmax,usegeno2,usegeno1,useql) -
c$$$     &              lTfdallvar(ipophost,n,npopmax,nlocd,nloch,
c$$$     &              nql,ncolt,nal,nalmax,usegeno2,usegeno1,useql) - 
c$$$     &              lTfdallvar(ipoprem,n,npopmax,nlocd,nloch,
c$$$     &              nql,ncolt,nal,nalmax,usegeno2,usegeno1,useql)
c$$$c               write(*,*) 'apres appel lTfdallvar'
c$$$c               write(6,*) 'merge apres prior contrib freq =',lratio 
c$$$            endif
c$$$            
c$$$            if(useqtc .eq. 1) then 
c$$$*     propose mean and variance of quantitative variables
c$$$               call propparqvmerge(qtc,nindiv,nqtc,indcell,c,ctmp,
c$$$     &              meanqtc,sdqtc,meanqtctmp,sdqtctmp,nnqtc,sqtc,
c$$$     &              ssqtc,npop,npopmax,nppmax,
c$$$     &              ksiqtc,kappaqtc,alphaqtc,betaqtc,ipophost,ipoprem,
c$$$     &              contriblr)
c$$$*     contrib prior and proposal 
c$$$               lratio = lratio + contriblr   
c$$$*     contrib likelihood quant. variable
c$$$*     argument usegeno2 is passed as 0 to ignore genotypes 
c$$$*     hence avoid having likelihood ratio of genotypes twice 
c$$$c$$$  lratio = lratio + 
c$$$c$$$  &              lratiogq(yy,f,c,ctmp,indcell,indcell,
c$$$c$$$     &              npopmax,nloc,nalmax,nindiv,nloc,nloc2,
c$$$c$$$  &              nppmax,ploidy,qtc,nqtc,meanqtc,sdqtc,
c$$$c$$$  &              meanqtctmp,sdqtctmp,0,useqtc)    
c$$$               lratio = lratio + lrallvar2(yy,z,ql,
c$$$     &              f,c,ctmp,indcell,indcell,
c$$$     &              npopmax,nlocd,nloch,nql,ncolt,nalmax,nindiv,
c$$$     &              nppmax,qtc,nqtc,meanqtc,sdqtc,meanqtctmp,sdqtctmp,
c$$$     &              0,0,0,useqtc,ploidy)   
c$$$c               write(6,*) 'merge apres llhood=',lratio 
c$$$            endif
c$$$c            write(*,*) 'lratio for merge=',lratio
c$$$
c$$$            lratio = dmin1(0.d0,lratio)
c$$$            alpha  = dexp(lratio)
c$$$            bern = ggrbinom(1.d0,alpha)      
c$$$            
c$$$            if(bern .eq. 1) then
c$$$               if(((usegeno2 .eq. 1) .or. (usegeno1 .eq. 1)) .or. 
c$$$     &              useql .eq. 1) then
c$$$*     propostion du nouveau tableau de freq et de derives
c$$$c$$$  call remfreq7bis(ipoprem,ipophost,
c$$$c$$$  &                  npop,npopmax,nloc,nal,
c$$$c$$$  &                  nalmax,f,ftmp,drift,drifttmp,fa,a,ptmp,ntmp)
c$$$                  call  remfallbis(ipoprem,ipophost,
c$$$     &                 npop,npopmax,nlocd,nlocd2,nloch,nql,ncolt,nal,
c$$$     &                 nalmax,f,ftmp,drift,drifttmp,fa,a,ptmp,ntmp,
c$$$     &                 usegeno2,usegeno1,useql)
c$$$               endif
c$$$               if(useqtc .eq. 1) then 
c$$$                  do iqtc = 1,nqtc
c$$$                     do ipop = 1,npopmax
c$$$                        meanqtc(ipop,iqtc) =  meanqtctmp(ipop,iqtc)
c$$$                        sdqtc(ipop,iqtc) =  sdqtctmp(ipop,iqtc)
c$$$                     enddo
c$$$                  enddo
c$$$               endif
c$$$               call accept5(nppmax,npopmax,ncolt,
c$$$     &              nalmax,nal,c,ctmp,f,ftmp,drift,drifttmp)
c$$$               npop = npop - 1
c$$$            endif
c$$$         endif
c$$$      endif
c$$$
c$$$c$$$      write(*,*) 'z =', z
c$$$c$$$      write(*,*) 'yy =', yy
c$$$c$$$      write(*,*) 'ql =', ql
c$$$c$$$      write(*,*) 'c =', c
c$$$c$$$      write(*,*) 'ctmp =',ctmp 
c$$$c$$$      write(*,*) 'f =', f  
c$$$c$$$      write(*,*) 'ftmp =', ftmp
c$$$c$$$      write(*,*) 'qtc=',qtc
c$$$c$$$      write(*,*) 'meanqtc=',meanqtc
c$$$c$$$      write(*,*) 'sdqtc=',sdqtc
c$$$c$$$      write(*,*) 'nnqtc=',nnqtc
c$$$c       write(*,*) 'fin smdallvar2' 
c$$$      end subroutine smdallvar2
c$$$************************************************************************



            
            
      

************************************************************************
*     propose parameters (means and variances) of quantitative variables 
*     in the update of a cell's membership  
*     and returns contrib of prior and proposal to log ratio
      subroutine propparqvudc(qtc,nindiv,nqtc,indcell,c,ctmp,
     &     meanqtc,sdqtc,meanqtctmp,sdqtctmp,nnqtc,sqtc,ssqtc,
     &     npop,npopmax,nppmax,
     &     ksiqtc,kappaqtc,alphaqtc,betaqtc,ipop1,ipop2,contriblr)
      implicit none
      integer nindiv,nqtc,indcell,c,ctmp,nnqtc,npop,npopmax,nppmax,
     &     ipop1,ipop2
      double precision qtc,meanqtc,sdqtc,
     &     meanqtctmp,sdqtctmp,sqtc,ssqtc,
     &     ksiqtc,kappaqtc,alphaqtc,betaqtc,contriblr
      dimension qtc(nindiv,nqtc),indcell(nindiv),c(nppmax),
     &     ctmp(nppmax),meanqtc(npopmax,nqtc),sdqtc(npopmax,nqtc),
     &     meanqtctmp(npopmax,nqtc),sdqtctmp(npopmax,nqtc),
     &     nnqtc(npopmax,nqtc),sqtc(npopmax,nqtc),
     &     ssqtc(npopmax,nqtc),
     &     ksiqtc(nqtc),kappaqtc(nqtc),alphaqtc(nqtc),betaqtc(nqtc)
      integer ipop,iindiv,iqtc
      double precision xx,vv,ap,bp,ggrnorm,ggrgam,ggdnorm,ggdgamma,junk
       contriblr = 0
*     compute empirical  sums  and sums of squares of quant. variables
*     for current state 
      do iqtc = 1,nqtc
         do ipop = 1,npopmax
            nnqtc(ipop,iqtc) = 0
            sqtc(ipop,iqtc) = 0
            ssqtc(ipop,iqtc) = 0
         enddo
      enddo
      do iqtc = 1,nqtc
         do iindiv = 1,nindiv
            ipop = c(indcell(iindiv))
            if(qtc(iindiv,iqtc) .ne. -999) then
               nnqtc(ipop,iqtc) = nnqtc(ipop,iqtc) + 1
               sqtc(ipop,iqtc) = sqtc(ipop,iqtc) + 
     &              qtc(iindiv,iqtc)
               ssqtc(ipop,iqtc) = ssqtc(ipop,iqtc) + 
     &              qtc(iindiv,iqtc)**2
            endif
         enddo
      enddo

*     contrib proposal  reverse move to log ratio
      do iqtc = 1,nqtc
*     mean 
         xx = (sqtc(ipop1,iqtc)+ksiqtc(iqtc)*kappaqtc(iqtc))/
     &        (nnqtc(ipop1,iqtc)+kappaqtc(iqtc))
         vv = 1 / (nnqtc(ipop1,iqtc)+kappaqtc(iqtc))
         contriblr = contriblr +
     &        ggdnorm(meanqtc(ipop1,iqtc),xx,dsqrt(vv),1)
c         write(*,*) 'ggdnorm(meanqtc(ipop1,iqtc),xx,dsqrt(vv),1)', 
c     &        ggdnorm(meanqtc(ipop1,iqtc),xx,dsqrt(vv),1)   
         xx = (sqtc(ipop2,iqtc)+ksiqtc(iqtc)*kappaqtc(iqtc))/
     &        (nnqtc(ipop2,iqtc)+kappaqtc(iqtc))
         vv = 1 / (nnqtc(ipop2,iqtc)+kappaqtc(iqtc))
         contriblr = contriblr +
     &        ggdnorm(meanqtc(ipop2,iqtc),xx,dsqrt(vv),1)
c        write(*,*) 'ggdnorm(meanqtc(ipop2,iqtc),xx,dsqrt(vv),1)', 
c    &        ggdnorm(meanqtc(ipop2,iqtc),xx,dsqrt(vv),1)          
*     variance
         ap = alphaqtc(iqtc) + 0.5*nnqtc(ipop1,iqtc)
         bp = betaqtc(iqtc) + 0.5*(ssqtc(ipop1,iqtc) - 
     &        2*sqtc(ipop1,iqtc)*meanqtc(ipop1,iqtc) + 
     &        nnqtc(ipop1,iqtc)*meanqtc(ipop1,iqtc)**2)
         contriblr = contriblr + 
     &        ggdgamma(1/sdqtc(ipop1,iqtc)**2,ap,1/bp,1)
c         write(*,*) 'ggdgamma(1/sdqtc(ipop1,iqtc)**2,ap,1/bp,1)',  
c     &        ggdgamma(1/sdqtc(ipop1,iqtc)**2,ap,1/bp,1)  
         ap = alphaqtc(iqtc) + 0.5*nnqtc(ipop2,iqtc)
         bp = betaqtc(iqtc) + 0.5*(ssqtc(ipop2,iqtc) - 
     &        2*sqtc(ipop2,iqtc)*meanqtc(ipop2,iqtc) + 
     &        nnqtc(ipop2,iqtc)*meanqtc(ipop2,iqtc)**2)
         contriblr = contriblr + 
     &        ggdgamma(1/sdqtc(ipop2,iqtc)**2,ap,1/bp,1)
c        write(*,*) 'ggdgamma(1/sdqtc(ipop2,iqtc)**2,ap,1/bp,1)',  
c     &        ggdgamma(1/sdqtc(ipop2,iqtc)**2,ap,1/bp,1)
      enddo
c      write(*,*) 'contrib proposal  reverse move=',contriblr
      junk = contriblr

*     compute empirical  sums  and sums of squares of quant. variables
*     for proposed state
      do iqtc = 1,nqtc
         do ipop = 1,npopmax
            nnqtc(ipop,iqtc) = 0
            sqtc(ipop,iqtc) = 0
            ssqtc(ipop,iqtc) = 0
            meanqtctmp(ipop,iqtc) = meanqtc(ipop,iqtc)
            sdqtctmp(ipop,iqtc) = sdqtc(ipop,iqtc)
         enddo
      enddo
      do iqtc = 1,nqtc
         do iindiv = 1,nindiv
            ipop = ctmp(indcell(iindiv))
            if(qtc(iindiv,iqtc) .ne. -999) then
               nnqtc(ipop,iqtc) = nnqtc(ipop,iqtc) + 1
               sqtc(ipop,iqtc) = sqtc(ipop,iqtc) + 
     &              qtc(iindiv,iqtc)
               ssqtc(ipop,iqtc) = ssqtc(ipop,iqtc) + 
     &              qtc(iindiv,iqtc)**2
            endif
         enddo
      enddo

*     sample mean from the full conditional
*     assuming current value for variance is 1
      do iqtc = 1,nqtc
*     propose mean pop ipop1 direct move 
*     and compute contrib proposal to log ratio
         xx = (sqtc(ipop1,iqtc)+ksiqtc(iqtc)*kappaqtc(iqtc))/
     &        (nnqtc(ipop1,iqtc)+kappaqtc(iqtc))
         vv = 1 / (nnqtc(ipop1,iqtc)+kappaqtc(iqtc))
         meanqtctmp(ipop1,iqtc) = ggrnorm(xx,dsqrt(vv))
         contriblr = contriblr - 
     &        ggdnorm(meanqtctmp(ipop1,iqtc),xx,dsqrt(vv),1)
c         write(*,*) ' ggdnorm(meanqtctmp(ipop1,iqtc),xx,dsqrt(vv),1)',
c     &        ggdnorm(meanqtctmp(ipop1,iqtc),xx,dsqrt(vv),1)
*     propose mean pop ipop2 direct move 
*     and compute contrib proposal to log ratio
         xx = (sqtc(ipop2,iqtc)+ksiqtc(iqtc)*kappaqtc(iqtc))/
     &        (nnqtc(ipop2,iqtc)+kappaqtc(iqtc))
         vv = 1 / (nnqtc(ipop2,iqtc)+kappaqtc(iqtc))
         meanqtctmp(ipop2,iqtc) = ggrnorm(xx,dsqrt(vv))
         contriblr = contriblr - 
     &        ggdnorm(meanqtctmp(ipop2,iqtc),xx,dsqrt(vv),1)
c         write(*,*) ' ggdnorm(meanqtctmp(ipop2,iqtc),xx,dsqrt(vv),1) ',
c     &     ggdnorm(meanqtctmp(ipop2,iqtc),xx,dsqrt(vv),1)
*     propose variance pop ipop1 direct move 
*     and compute contrib proposal to log ratio
*     sample variance from full conditional
         ap = alphaqtc(iqtc) + 0.5*nnqtc(ipop1,iqtc)
         bp = betaqtc(iqtc) + 0.5*(ssqtc(ipop1,iqtc) - 
     &        2*sqtc(ipop1,iqtc)*meanqtctmp(ipop1,iqtc) + 
     &        nnqtc(ipop1,iqtc)*meanqtctmp(ipop1,iqtc)**2)
         sdqtctmp(ipop1,iqtc) = 1/dsqrt(ggrgam(ap,1/bp))
         contriblr = contriblr - 
     &        ggdgamma(1/sdqtctmp(ipop1,iqtc)**2,ap,1/bp,1)
c         write(*,*) 'ggdgamma(1/sdqtctmp(ipop1,iqtc)**2,ap,1/bp,1) ',
c     &     ggdgamma(1/sdqtctmp(ipop1,iqtc)**2,ap,1/bp,1)
*     propose variance pop ipop2 direct move 
*     and compute contrib proposal to log ratio
*     sample variance from full conditional
         ap = alphaqtc(iqtc) + 0.5*nnqtc(ipop2,iqtc)
         bp = betaqtc(iqtc) + 0.5*(ssqtc(ipop2,iqtc) - 
     &        2*sqtc(ipop2,iqtc)*meanqtctmp(ipop2,iqtc) + 
     &        nnqtc(ipop2,iqtc)*meanqtctmp(ipop2,iqtc)**2)
         sdqtctmp(ipop2,iqtc) = 1/dsqrt(ggrgam(ap,1/bp))
         contriblr = contriblr - 
     &        ggdgamma(1/sdqtctmp(ipop2,iqtc)**2,ap,1/bp,1)
c         write(*,*) 'ggdgamma(1/sdqtctmp(ipop2,iqtc)**2,ap,1/bp,1) ',
c     &     ggdgamma(1/sdqtctmp(ipop2,iqtc)**2,ap,1/bp,1)
      enddo 
c      write(*,*) 'contrib proposal  direct move=',contriblr-junk
      junk = contriblr

      do iqtc = 1,nqtc
*     contrib prior of means to log ratio
         contriblr = contriblr + 
     &        ggdnorm(meanqtctmp(ipop1,iqtc),
     &        ksiqtc(iqtc),1/dsqrt(kappaqtc(iqtc)),1) + 
     &        ggdnorm(meanqtctmp(ipop2,iqtc),
     &        ksiqtc(iqtc),1/dsqrt(kappaqtc(iqtc)),1)  - 
     &        ggdnorm(meanqtc(ipop1,iqtc),
     &        ksiqtc(iqtc),1/dsqrt(kappaqtc(iqtc)),1)  - 
     &        ggdnorm(meanqtc(ipop2,iqtc),
     &        ksiqtc(iqtc),1/dsqrt(kappaqtc(iqtc)),1) 
      enddo
c      write(*,*) 'contrib prior means=',contriblr-junk
      junk = contriblr
      do iqtc = 1,nqtc
*     contrib prior of variances to log ratio
         contriblr = contriblr +
     &        ggdgamma(1/sdqtctmp(ipop1,iqtc)**2,
     &        alphaqtc(iqtc),1/betaqtc(iqtc),1) + 
     &        ggdgamma(1/sdqtctmp(ipop2,iqtc)**2,
     &        alphaqtc(iqtc),1/betaqtc(iqtc),1) -
     &        ggdgamma(1/sdqtc(ipop1,iqtc)**2,
     &        alphaqtc(iqtc),1/betaqtc(iqtc),1) -
     &        ggdgamma(1/sdqtc(ipop2,iqtc)**2,
     &        alphaqtc(iqtc),1/betaqtc(iqtc),1)
c$$$         write(*,*) 'contrib prior variances:'
c$$$         write(*,*) 'ggdgamma(1/sdqtctmp(ipop1,iqtc)**2,
c$$$     & alphaqtc(iqtc),1/betaqtc(iqtc),1=',
c$$$     &     ggdgamma(1/sdqtctmp(ipop1,iqtc)**2,
c$$$     &        alphaqtc(iqtc),1/betaqtc(iqtc),1)
c$$$         write(*,*) 'ggdgamma(1/sdqtctmp(ipop2,iqtc)**2,
c$$$     & alphaqtc(iqtc),1/betaqtc(iqtc),1)=',
c$$$     &        ggdgamma(1/sdqtctmp(ipop2,iqtc)**2,
c$$$     &        alphaqtc(iqtc),1/betaqtc(iqtc),1)
c$$$      write(*,*) 'ggdgamma(1/sdqtc(ipop1,iqtc)**2,
c$$$     & alphaqtc(iqtc),1/betaqtc(iqtc),1)=',
c$$$     &        ggdgamma(1/sdqtc(ipop1,iqtc)**2,
c$$$     &        alphaqtc(iqtc),1/betaqtc(iqtc),1)
c$$$      write(*,*) 'ggdgamma(1/sdqtc(ipop2,iqtc)**2,
c$$$     & alphaqtc(iqtc),1/betaqtc(iqtc),1)=',
c$$$     &     ggdgamma(1/sdqtc(ipop2,iqtc)**2,
c$$$     &     alphaqtc(iqtc),1/betaqtc(iqtc),1)
c$$$      write(*,*) 'global contrib prior variances=',contriblr-junk
      enddo

c$$$      write(*,*) 'in propparqvudc, meanqtc='
c$$$      write(*,*) (meanqtc(1,iqtc), iqtc=1,3)
c$$$      write(*,*) (meanqtc(2,iqtc), iqtc=1,3)
c$$$      write(*,*) 'in propparqvudc, sdqtc='
c$$$      write(*,*) (sdqtc(1,iqtc), iqtc=1,3)
c$$$      write(*,*) (sdqtc(2,iqtc), iqtc=1,3)
c$$$      write(*,*) 'in propparqvudc, meanqtctmp='
c$$$      write(*,*) (meanqtctmp(1,iqtc), iqtc=1,3)
c$$$      write(*,*) (meanqtctmp(2,iqtc), iqtc=1,3)
c$$$      write(*,*) 'in propparqvudc, sdqtctmp='
c$$$      write(*,*) (sdqtctmp(1,iqtc), iqtc=1,3)
c$$$      write(*,*) (sdqtctmp(2,iqtc), iqtc=1,3)
c$$$      write(*,*) 'contriblr=',contriblr

      end subroutine propparqvudc



            
      

************************************************************************
*     propose parameters (means and variances) of quantitative variables 
*     in the update of a cell's membership  
*     and returns contrib of prior and proposal to log ratio
*     new mean and variances proposed according to full conditional
*     prior mean-variance: Normal Scaled Inverse Gamma
      subroutine propparqvudc3(qtc,nindiv,nqtc,indcell,c,ctmp,
     &     meanqtc,sdqtc,meanqtctmp,sdqtctmp,nnqtc,sqtc,ssqtc,
     &     npop,npopmax,nppmax,
     &     ksiqtc,kappaqtc,alphaqtc,betaqtc,ipop1,ipop2,contriblr)
      implicit none
      integer nindiv,nqtc,indcell,c,ctmp,nnqtc,npop,npopmax,nppmax,
     &     ipop1,ipop2
      double precision qtc,meanqtc,sdqtc,
     &     meanqtctmp,sdqtctmp,sqtc,ssqtc,
     &     ksiqtc,kappaqtc,alphaqtc,betaqtc,contriblr
      dimension qtc(nindiv,nqtc),indcell(nindiv),c(nppmax),
     &     ctmp(nppmax),meanqtc(npopmax,nqtc),sdqtc(npopmax,nqtc),
     &     meanqtctmp(npopmax,nqtc),sdqtctmp(npopmax,nqtc),
     &     nnqtc(npopmax,nqtc),sqtc(npopmax,nqtc),
     &     ssqtc(npopmax,nqtc),
     &     ksiqtc(nqtc),kappaqtc(nqtc),alphaqtc(nqtc),betaqtc(nqtc)
      integer ipop,iindiv,iqtc
      double precision xx,vv,ap,bp,ggrnorm,ggrgam,ggdnorm,ggdgamma,res,
     &     junk
c$$$      write(*,*) 'debut propparqvudc3' 
c$$$      write(*,*) 'meanqtc=',meanqtc
c$$$      write(*,*) 'meanqtctmp=',meanqtctmp
c$$$      write(*,*) 'sdqtc=',sdqtc
c$$$      write(*,*) 'sdqtctmp=',sdqtctmp
c$$$      write(*,*) 'ipop1 ipop2=',ipop1, ipop2
c$$$      write(*,*) 'ksiqtc=',ksiqtc
c$$$      write(*,*) 'kappaqtc=',kappaqtc
c$$$      write(*,*) 'alphaqtc=',alphaqtc
c$$$      write(*,*) 'betaqtc=',betaqtc
      call contproprev(nqtc,ipop1,ipop2,npopmax,nindiv,nppmax,
     &     nnqtc,meanqtc,sdqtc,c,
     &     indcell,sqtc,ssqtc,alphaqtc,betaqtc,qtc,
     &     ksiqtc,kappaqtc,res)
      contriblr = res 
c      write(*,*) 'after controprev contriblr=',contriblr

      call contpropdir(nqtc,ipop1,ipop2,nindiv,nnqtc,npopmax,
     &     nppmax,ctmp,indcell,qtc,meanqtctmp,sdqtctmp,sqtc,ssqtc,
     &     alphaqtc,betaqtc,ksiqtc,kappaqtc,res)
      contriblr = contriblr + res 
c      write(*,*) 'after contrpropdir contriblr=',contriblr

      call contprior(ipop1,ipop2,npopmax,nqtc,meanqtc,sdqtc,
     &     meanqtctmp,sdqtctmp,ksiqtc,kappaqtc,alphaqtc,betaqtc,res)
      contriblr = contriblr + res 
c      write(*,*) 'after contrprior contriblr=',contriblr

c$$$      write(*,*) 'in propparqvudc, meanqtc='
c$$$      write(*,*) (meanqtc(1,iqtc), iqtc=1,3)
c$$$      write(*,*) (meanqtc(2,iqtc), iqtc=1,3)
c$$$      write(*,*) 'in propparqvudc, sdqtc='
c$$$      write(*,*) (sdqtc(1,iqtc), iqtc=1,3)
c$$$      write(*,*) (sdqtc(2,iqtc), iqtc=1,3)
c$$$      write(*,*) 'in propparqvudc, meanqtctmp='
c$$$      write(*,*) (meanqtctmp(1,iqtc), iqtc=1,3)
c$$$      write(*,*) (meanqtctmp(2,iqtc), iqtc=1,3)
c$$$      write(*,*) 'in propparqvudc, sdqtctmp='
c$$$      write(*,*) (sdqtctmp(1,iqtc), iqtc=1,3)
c$$$      write(*,*) (sdqtctmp(2,iqtc), iqtc=1,3)
c$$$      write(*,*) 'contriblr=',contriblr
c$$$      write(*,*) 'fin propparqvudc3' 
c$$$      write(*,*) 'meanqtc=',meanqtc
c$$$      write(*,*) 'meanqtctmp=',meanqtctmp
c$$$      write(*,*) 'sdqtc=',sdqtc
c$$$      write(*,*) 'sdqtctmp=',sdqtctmp
      end subroutine propparqvudc3


************************************************************************
*     propose parameters (means and variances) of quantitative variables 
*     in the update of a cell's membership  
*     and returns contrib of prior and proposal to log ratio
*     new mean and variances proposed according to full conditional
*     prior mean-variance: Normal Scaled Inverse Gamma
      subroutine propparqvudc4(qtc,nindiv,nqtc,indcell,c,ctmp,
     &     meanqtc,sdqtc,meanqtctmp,sdqtctmp,nnqtc,sqtc,ssqtc,
     &     npop,npopmax,nppmax,ksiqtc,kappaqtc,alphaqtc,betaqtc,
     &     contriblr)
      implicit none
      integer nindiv,nqtc,indcell,c,ctmp,nnqtc,npop,npopmax,nppmax,
     &     ipop1,ipop2
      double precision qtc,meanqtc,sdqtc,
     &     meanqtctmp,sdqtctmp,sqtc,ssqtc,
     &     ksiqtc,kappaqtc,alphaqtc,betaqtc,contriblr
      dimension qtc(nindiv,nqtc),indcell(nindiv),c(nppmax),
     &     ctmp(nppmax),meanqtc(npopmax,nqtc),sdqtc(npopmax,nqtc),
     &     meanqtctmp(npopmax,nqtc),sdqtctmp(npopmax,nqtc),
     &     nnqtc(npopmax,nqtc),sqtc(npopmax,nqtc),
     &     ssqtc(npopmax,nqtc),
     &     ksiqtc(nqtc),kappaqtc(nqtc),alphaqtc(nqtc),betaqtc(nqtc)
      integer ipop,iindiv,iqtc
      double precision xx,vv,ap,bp,ggrnorm,ggrgam,ggdnorm,ggdgamma,res
      call contproprev4(nqtc,npop,npopmax,nindiv,nppmax,
     &     nnqtc,meanqtc,sdqtc,c,
     &     indcell,sqtc,ssqtc,alphaqtc,betaqtc,qtc,
     &     ksiqtc,kappaqtc,res)
      contriblr = res 
      call contpropdir4(nqtc,npop,nindiv,nnqtc,npopmax,
     &     nppmax,ctmp,indcell,qtc,meanqtctmp,sdqtctmp,sqtc,ssqtc,
     &     alphaqtc,betaqtc,ksiqtc,kappaqtc,res)
      contriblr = contriblr + res 
      call contprior4(npop,npopmax,nqtc,meanqtc,sdqtc,
     &     meanqtctmp,sdqtctmp,ksiqtc,kappaqtc,alphaqtc,betaqtc,res)
      contriblr = contriblr + res 
      end subroutine propparqvudc4





******************************************************
*     contrib proposal of reverse move in MH ratio 
*     in joint update of c,mu,sigma
*     proposal: full conditionnal p(mu,sigma|c,x,...)
      subroutine contproprev(nqtc,ipop1,ipop2,npopmax,nindiv,nppmax,
     &     nnqtc,meanqtc,sdqtc,c,
     &     indcell,sqtc,ssqtc,alphaqtc,betaqtc,qtc,
     &     ksiqtc,kappaqtc,res)
      implicit none
      integer nqtc,ipop1,ipop2,nnqtc,npopmax,nindiv,nppmax,c,indcell
      double precision meanqtc,sdqtc,sqtc,ssqtc,res,qtc,
     &     alphaqtc,betaqtc,ksiqtc,kappaqtc
      integer iqtc,ipop,iindiv
      double precision xx,vv,ap,bp,ggdnorm,ggdgamma
      dimension meanqtc(npopmax,nqtc),sdqtc(npopmax,nqtc),
     &     nnqtc(npopmax,nqtc),sqtc(npopmax,nqtc),
     &     ssqtc(npopmax,nqtc),ksiqtc(nqtc),kappaqtc(nqtc),
     &     alphaqtc(nqtc),betaqtc(nqtc),c(nppmax),indcell(nindiv),
     &     qtc(nindiv,nqtc)
      res = 0
*     nnqtc sample size in the various clusters
*     sqtc: sum 
*     ssqtc:  Sum_i (x_i - xbar)^2
      do iqtc = 1,nqtc
         do ipop = 1,npopmax
            nnqtc(ipop,iqtc) = 0
            sqtc(ipop,iqtc) = 0
            ssqtc(ipop,iqtc) = 0
         enddo
      enddo
      do iqtc = 1,nqtc
         do iindiv = 1,nindiv
            ipop = c(indcell(iindiv))
            if(qtc(iindiv,iqtc) .ne. -999) then
               nnqtc(ipop,iqtc) = nnqtc(ipop,iqtc) + 1
               sqtc(ipop,iqtc) = sqtc(ipop,iqtc) + 
     &              qtc(iindiv,iqtc)
               ssqtc(ipop,iqtc) = ssqtc(ipop,iqtc) + 
     &              qtc(iindiv,iqtc)**2
            endif
         enddo
      enddo
      do iqtc = 1,nqtc
         do ipop = 1,npopmax
            if(nnqtc(ipop,iqtc) .ne. 0) then
               ssqtc(ipop,iqtc) = ssqtc(ipop,iqtc) -
     &              (sqtc(ipop,iqtc)**2) / dble(nnqtc(ipop,iqtc))
            else
               ssqtc(ipop,iqtc) = 0
            endif
         enddo
      enddo

*     contrib proposal  reverse move to log ratio
      do iqtc = 1,nqtc
*     mean 
         xx = (sqtc(ipop1,iqtc)+ksiqtc(iqtc)*kappaqtc(iqtc))/
     &        (nnqtc(ipop1,iqtc)+kappaqtc(iqtc))
         vv = sdqtc(ipop1,iqtc)**2 / (nnqtc(ipop1,iqtc)+kappaqtc(iqtc))
         res = res +
     &        ggdnorm(meanqtc(ipop1,iqtc),xx,dsqrt(vv),1)
c$$$         write(*,*) 'ggdnorm(meanqtc(ipop1,iqtc),xx,dsqrt(vv),1)', 
c$$$     &        ggdnorm(meanqtc(ipop1,iqtc),xx,dsqrt(vv),1)   
         xx = (sqtc(ipop2,iqtc)+ksiqtc(iqtc)*kappaqtc(iqtc))/
     &        (nnqtc(ipop2,iqtc)+kappaqtc(iqtc))
         vv = sdqtc(ipop2,iqtc)**2 / (nnqtc(ipop2,iqtc)+kappaqtc(iqtc))
         res = res +
     &        ggdnorm(meanqtc(ipop2,iqtc),xx,dsqrt(vv),1)
c$$$         write(*,*) 'ggdnorm(meanqtc(ipop2,iqtc),xx,dsqrt(vv),1)', 
c$$$     &        ggdnorm(meanqtc(ipop2,iqtc),xx,dsqrt(vv),1)          
*     variance
         ap = alphaqtc(iqtc) + 0.5*nnqtc(ipop1,iqtc)
         if(nnqtc(ipop1,iqtc) .ne. 0) then
            bp = betaqtc(iqtc) + 0.5*ssqtc(ipop1,iqtc) + 
     &           0.5*(nnqtc(ipop1,iqtc)*kappaqtc(iqtc))/
     &           (nnqtc(ipop1,iqtc)+kappaqtc(iqtc)) *
     &           (sqtc(ipop1,iqtc)/dble(nnqtc(ipop1,iqtc)) - 
     &           ksiqtc(iqtc))**2
         else 
            bp = betaqtc(iqtc)
         endif
         res = res + ggdgamma(1/sdqtc(ipop1,iqtc)**2,ap,1/bp,1)
c$$$         write(*,*) 'ggdgamma(1/sdqtc(ipop1,iqtc)**2,ap,1/bp,1)',  
c$$$     &        ggdgamma(1/sdqtc(ipop1,iqtc)**2,ap,1/bp,1)  
         ap = alphaqtc(iqtc) + 0.5*nnqtc(ipop2,iqtc)
         if(nnqtc(ipop2,iqtc) .ne. 0) then
            bp = betaqtc(iqtc) + 0.5*ssqtc(ipop2,iqtc) + 
     &           0.5*(nnqtc(ipop2,iqtc)*kappaqtc(iqtc))/
     &           (nnqtc(ipop2,iqtc)+kappaqtc(iqtc)) *
     &           (sqtc(ipop2,iqtc)/dble(nnqtc(ipop2,iqtc)) - 
     &           ksiqtc(iqtc))**2
         else
            bp = betaqtc(iqtc)
         endif
            res = res + ggdgamma(1/sdqtc(ipop2,iqtc)**2,ap,1/bp,1)
c$$$         write(*,*) 'ggdgamma(1/sdqtc(ipop2,iqtc)**2,ap,1/bp,1)',  
c$$$     &        ggdgamma(1/sdqtc(ipop2,iqtc)**2,ap,1/bp,1)
      enddo
c$$$      write(*,*) 'contrib proposal  reverse move=',res
      end subroutine contproprev


******************************************************
*     contrib proposal of reverse move in MH ratio 
*     in joint update of c,mu,sigma
*     proposal: full conditionnal p(mu,sigma|c,x,...)
      subroutine contproprev4(nqtc,npop,npopmax,nindiv,nppmax,
     &     nnqtc,meanqtc,sdqtc,c,
     &     indcell,sqtc,ssqtc,alphaqtc,betaqtc,qtc,
     &     ksiqtc,kappaqtc,res)
      implicit none
      integer nqtc,nnqtc,npop,npopmax,nindiv,nppmax,c,indcell
      double precision meanqtc,sdqtc,sqtc,ssqtc,res,qtc,
     &     alphaqtc,betaqtc,ksiqtc,kappaqtc
      integer iqtc,ipop,iindiv
      double precision xx,vv,ap,bp,ggdnorm,ggdgamma
      dimension meanqtc(npopmax,nqtc),sdqtc(npopmax,nqtc),
     &     nnqtc(npopmax,nqtc),sqtc(npopmax,nqtc),
     &     ssqtc(npopmax,nqtc),ksiqtc(nqtc),kappaqtc(nqtc),
     &     alphaqtc(nqtc),betaqtc(nqtc),c(nppmax),indcell(nindiv),
     &     qtc(nindiv,nqtc)
      res = 0
*     nnqtc sample size in the various clusters
*     sqtc: sum 
*     ssqtc:  Sum_i (x_i - xbar)^2
      do iqtc = 1,nqtc
         do ipop = 1,npopmax
            nnqtc(ipop,iqtc) = 0
            sqtc(ipop,iqtc) = 0
            ssqtc(ipop,iqtc) = 0
         enddo
      enddo
      do iqtc = 1,nqtc
         do iindiv = 1,nindiv
            ipop = c(indcell(iindiv))
            if(qtc(iindiv,iqtc) .ne. -999) then
               nnqtc(ipop,iqtc) = nnqtc(ipop,iqtc) + 1
               sqtc(ipop,iqtc) = sqtc(ipop,iqtc) + 
     &              qtc(iindiv,iqtc)
               ssqtc(ipop,iqtc) = ssqtc(ipop,iqtc) + 
     &              qtc(iindiv,iqtc)**2
            endif
         enddo
      enddo
      do iqtc = 1,nqtc
         do ipop = 1,npopmax
            if(nnqtc(ipop,iqtc) .ne. 0) then
               ssqtc(ipop,iqtc) = ssqtc(ipop,iqtc) -
     &              (sqtc(ipop,iqtc)**2) / dble(nnqtc(ipop,iqtc))
            else
               ssqtc(ipop,iqtc) = 0
            endif
         enddo
      enddo

*     contrib proposal  reverse move to log ratio
      do ipop=1,npop
         do iqtc = 1,nqtc
*     mean 
            xx = (sqtc(ipop,iqtc)+ksiqtc(iqtc)*kappaqtc(iqtc))/
     &           (nnqtc(ipop,iqtc)+kappaqtc(iqtc))
            vv = sdqtc(ipop,iqtc)**2 / 
     &           (nnqtc(ipop,iqtc)+kappaqtc(iqtc))
            res = res +
     &           ggdnorm(meanqtc(ipop,iqtc),xx,dsqrt(vv),1)
*     variance
            ap = alphaqtc(iqtc) + 0.5*nnqtc(ipop,iqtc)
            if(nnqtc(ipop,iqtc) .ne. 0) then
               bp = betaqtc(iqtc) + 0.5*ssqtc(ipop,iqtc) + 
     &              0.5*(nnqtc(ipop,iqtc)*kappaqtc(iqtc))/
     &              (nnqtc(ipop,iqtc)+kappaqtc(iqtc)) *
     &              (sqtc(ipop,iqtc)/dble(nnqtc(ipop,iqtc)) - 
     &              ksiqtc(iqtc))**2
            else 
               bp = betaqtc(iqtc)
            endif
            res = res + ggdgamma(1/sdqtc(ipop,iqtc)**2,ap,1/bp,1)
         enddo
      enddo
      end subroutine contproprev4




******************************************************
*     propose mu,sigma in joint update of c,mu,sigma
*     compute contrib proposal of direct move in MH ratio 
*     proposal: full conditionnal p(mu,sigma|c,x,...)
      subroutine contpropdir(nqtc,ipop1,ipop2,nindiv,nnqtc,npopmax,
     &     nppmax,ctmp,indcell,qtc,meanqtctmp,sdqtctmp,sqtc,ssqtc,
     &     alphaqtc,betaqtc,ksiqtc,kappaqtc,res)
      implicit none
      integer nqtc,ipop1,ipop2,nnqtc,npopmax,nindiv,ctmp,indcell,nppmax
      double precision meanqtctmp,sdqtctmp,sqtc,ssqtc,res,
     &     alphaqtc,betaqtc,ksiqtc,kappaqtc,qtc
      integer iqtc,ipop,iindiv
      double precision xx,vv,ap,bp,ggdnorm,ggrnorm,ggdgamma,ggrgam,junk
      dimension meanqtctmp(npopmax,nqtc),sdqtctmp(npopmax,nqtc),
     &     nnqtc(npopmax,nqtc),sqtc(npopmax,nqtc),
     &     ssqtc(npopmax,nqtc),ksiqtc(nqtc),kappaqtc(nqtc),
     &     alphaqtc(nqtc),betaqtc(nqtc),ctmp(nppmax),indcell(nindiv),
     &     qtc(nindiv,nqtc)
c$$$      write(*,*) 'debut contpropdir meanqtc='
c$$$      write(*,*) 'meanqtctmp=',meanqtctmp
c$$$      write(*,*) 'sdqtctmp=',sdqtctmp
c$$$      write(*,*) 'ksiqtc=',ksiqtc
c$$$      write(*,*) 'kappaqtc=',kappaqtc
c$$$      write(*,*) 'alphaqtc=',alphaqtc
c$$$      write(*,*) 'betaqtc=',betaqtc
      res = 0
*     compute empirical  sums  and sums of squares of quant. variables
*     for proposed state
      do iqtc = 1,nqtc
         do ipop = 1,npopmax
            nnqtc(ipop,iqtc) = 0
            sqtc(ipop,iqtc) = 0
            ssqtc(ipop,iqtc) = 0
         enddo
      enddo
      do iqtc = 1,nqtc
         do iindiv = 1,nindiv
            ipop = ctmp(indcell(iindiv))
            if(qtc(iindiv,iqtc) .ne. -999) then
               nnqtc(ipop,iqtc) = nnqtc(ipop,iqtc) + 1
               sqtc(ipop,iqtc) = sqtc(ipop,iqtc) + 
     &              qtc(iindiv,iqtc)
               ssqtc(ipop,iqtc) = ssqtc(ipop,iqtc) + 
     &              qtc(iindiv,iqtc)**2
            endif
         enddo
      enddo
      do iqtc = 1,nqtc
         do ipop = 1,npopmax
            if(nnqtc(ipop,iqtc) .ne. 0) then
               ssqtc(ipop,iqtc) = ssqtc(ipop,iqtc) -
     &              (sqtc(ipop,iqtc)**2) / dble(nnqtc(ipop,iqtc))
            else
               ssqtc(ipop,iqtc) = 0
            endif
         enddo
      enddo
c$$$      write(*,*) 'nnqtc=',nnqtc
c$$$      write(*,*) 'sqtc=',sqtc
c$$$      write(*,*) 'ssqtc=',ssqtc


*     sample mean from the full conditional
      do iqtc = 1,nqtc
*     propose mean pop ipop1 direct move 
*     and compute contrib proposal to log ratio

c$$$         write(*,*) ' ggdnorm(meanqtctmp(ipop2,iqtc),xx,dsqrt(vv),1) ',
c$$$     &     ggdnorm(meanqtctmp(ipop2,iqtc),xx,dsqrt(vv),1)
*     propose variance pop ipop1 direct move 
*     and compute contrib proposal to log ratio
*     sample variance from full conditional
         ap = alphaqtc(iqtc) + 0.5*nnqtc(ipop1,iqtc)
         if(nnqtc(ipop1,iqtc) .ne. 0) then 
            bp = betaqtc(iqtc) + 0.5*ssqtc(ipop1,iqtc) + 
     &           0.5*(nnqtc(ipop1,iqtc)*kappaqtc(iqtc))/
     &           (nnqtc(ipop1,iqtc)+kappaqtc(iqtc)) *
     &           (sqtc(ipop1,iqtc)/dble(nnqtc(ipop1,iqtc)) - 
     &           ksiqtc(iqtc))**2
         else
            bp = betaqtc(iqtc)
         endif
c         write(*,*) 'ap=',ap
c         write(*,*) 'bp=',bp
         sdqtctmp(ipop1,iqtc) = 1/dsqrt(ggrgam(ap,1/bp))
c         write(*,*) 'sdqtctmp(ipop1,iqtc)=',sdqtctmp(ipop1,iqtc)
         res = res - 
     &        ggdgamma(1/sdqtctmp(ipop1,iqtc)**2,ap,1/bp,1)
 
c$$$         write(*,*) 'ggdgamma(1/sdqtctmp(ipop1,iqtc)**2,ap,1/bp,1) ',
c$$$     &     ggdgamma(1/sdqtctmp(ipop1,iqtc)**2,ap,1/bp,1)
*     propose variance pop ipop2 direct move 
*     and compute contrib proposal to log ratio
*     sample variance from full conditional
         ap = alphaqtc(iqtc) + 0.5*nnqtc(ipop2,iqtc)
         if(nnqtc(ipop2,iqtc) .ne. 0) then 
            bp = betaqtc(iqtc) + 0.5*ssqtc(ipop2,iqtc) + 
     &           0.5*(nnqtc(ipop2,iqtc)*kappaqtc(iqtc))/
     &           (nnqtc(ipop2,iqtc)+kappaqtc(iqtc)) *
     &           (sqtc(ipop2,iqtc)/dble(nnqtc(ipop2,iqtc)) - 
     &           ksiqtc(iqtc))**2
         else
            bp = betaqtc(iqtc)
         endif
c         write(*,*) 'ap=',ap
c         write(*,*) 'bp=',bp
         sdqtctmp(ipop2,iqtc) = 1/dsqrt(ggrgam(ap,1/bp))
c         write(*,*) 'sdqtctmp(ipop2,iqtc)=',sdqtctmp(ipop2,iqtc)
         res = res - 
     &        ggdgamma(1/sdqtctmp(ipop2,iqtc)**2,ap,1/bp,1)
         
         xx = (sqtc(ipop1,iqtc)+ksiqtc(iqtc)*kappaqtc(iqtc))/
     &        (nnqtc(ipop1,iqtc)+kappaqtc(iqtc))
         vv = (sdqtctmp(ipop1,iqtc)**2)/
     &        (nnqtc(ipop1,iqtc)+kappaqtc(iqtc))
c$$$         write(*,*) 'xx=',xx
c$$$         write(*,*) 'vv=',vv
         meanqtctmp(ipop1,iqtc) = ggrnorm(xx,dsqrt(vv))
         res = res - 
     &        ggdnorm(meanqtctmp(ipop1,iqtc),xx,dsqrt(vv),1)

*     propose mean pop ipop2 direct move 
*     and compute contrib proposal to log ratio
         xx = (sqtc(ipop2,iqtc)+ksiqtc(iqtc)*kappaqtc(iqtc))/
     &        (nnqtc(ipop2,iqtc)+kappaqtc(iqtc))
         vv = (sdqtctmp(ipop2,iqtc)**2)/
     &        (nnqtc(ipop2,iqtc)+kappaqtc(iqtc))
c$$$         write(*,*) 'xx=',xx
c$$$         write(*,*) 'vv=',vv
         meanqtctmp(ipop2,iqtc) = ggrnorm(xx,dsqrt(vv))
         res = res - 
     &        ggdnorm(meanqtctmp(ipop2,iqtc),xx,dsqrt(vv),1)
      enddo 
c$$$      write(*,*) 'contrib proposal  direct move=',res-junk

c$$$      write(*,*) 'fin contpropdir meanqtc='
c$$$      write(*,*) 'meanqtctmp=',meanqtctmp
c$$$      write(*,*) 'sdqtctmp=',sdqtctmp
c$$$      write(*,*) 'ksiqtc=',ksiqtc
c$$$      write(*,*) 'kappaqtc=',kappaqtc
c$$$      write(*,*) 'alphaqtc=',alphaqtc
c$$$      write(*,*) 'betaqtc=',betaqtc
      end subroutine contpropdir



******************************************************
*     propose mu,sigma in joint update of c,mu,sigma
*     compute contrib proposal of direct move in MH ratio 
*     proposal: full conditionnal p(mu,sigma|c,x,...)
      subroutine contpropdir4(nqtc,npop,nindiv,nnqtc,npopmax,
     &     nppmax,ctmp,indcell,qtc,meanqtctmp,sdqtctmp,sqtc,ssqtc,
     &     alphaqtc,betaqtc,ksiqtc,kappaqtc,res)
      implicit none
      integer nqtc,nnqtc,npop,npopmax,nindiv,ctmp,indcell,nppmax
      double precision meanqtctmp,sdqtctmp,sqtc,ssqtc,res,
     &     alphaqtc,betaqtc,ksiqtc,kappaqtc,qtc
      integer iqtc,ipop,iindiv
      double precision xx,vv,ap,bp,ggdnorm,ggrnorm,ggdgamma,ggrgam,junk
      dimension meanqtctmp(npopmax,nqtc),sdqtctmp(npopmax,nqtc),
     &     nnqtc(npopmax,nqtc),sqtc(npopmax,nqtc),
     &     ssqtc(npopmax,nqtc),ksiqtc(nqtc),kappaqtc(nqtc),
     &     alphaqtc(nqtc),betaqtc(nqtc),ctmp(nppmax),indcell(nindiv),
     &     qtc(nindiv,nqtc)
      res = 0
*     compute empirical  sums  and sums of squares of quant. variables
*     for proposed state
      do iqtc = 1,nqtc
         do ipop = 1,npopmax
            nnqtc(ipop,iqtc) = 0
            sqtc(ipop,iqtc) = 0
            ssqtc(ipop,iqtc) = 0
         enddo
      enddo
      do iqtc = 1,nqtc
         do iindiv = 1,nindiv
            ipop = ctmp(indcell(iindiv))
            if(qtc(iindiv,iqtc) .ne. -999) then
               nnqtc(ipop,iqtc) = nnqtc(ipop,iqtc) + 1
               sqtc(ipop,iqtc) = sqtc(ipop,iqtc) + 
     &              qtc(iindiv,iqtc)
               ssqtc(ipop,iqtc) = ssqtc(ipop,iqtc) + 
     &              qtc(iindiv,iqtc)**2
            endif
         enddo
      enddo
      do iqtc = 1,nqtc
         do ipop = 1,npopmax
            if(nnqtc(ipop,iqtc) .ne. 0) then
               ssqtc(ipop,iqtc) = ssqtc(ipop,iqtc) -
     &              (sqtc(ipop,iqtc)**2) / dble(nnqtc(ipop,iqtc))
            else
               ssqtc(ipop,iqtc) = 0
            endif
         enddo
      enddo

*     sample means from the full conditional
      do ipop=1,npop
         do iqtc = 1,nqtc
*     propose mean pop ipop direct move 
*     and compute contrib proposal to log ratio
            
*     propose variance pop ipop direct move 
*     and compute contrib proposal to log ratio
*     sample variance from full conditional
            ap = alphaqtc(iqtc) + 0.5*nnqtc(ipop,iqtc)
            if(nnqtc(ipop,iqtc) .ne. 0) then 
               bp = betaqtc(iqtc) + 0.5*ssqtc(ipop,iqtc) + 
     &              0.5*(nnqtc(ipop,iqtc)*kappaqtc(iqtc))/
     &              (nnqtc(ipop,iqtc)+kappaqtc(iqtc)) *
     &              (sqtc(ipop,iqtc)/dble(nnqtc(ipop,iqtc)) - 
     &              ksiqtc(iqtc))**2
            else
               bp = betaqtc(iqtc)
            endif
            sdqtctmp(ipop,iqtc) = 1/dsqrt(ggrgam(ap,1/bp))
c     write(*,*) 'sdqtctmp(ipop,iqtc)=',sdqtctmp(ipop,iqtc)
            res = res - 
     &           ggdgamma(1/sdqtctmp(ipop,iqtc)**2,ap,1/bp,1)
            
            xx = (sqtc(ipop,iqtc)+ksiqtc(iqtc)*kappaqtc(iqtc))/
     &           (nnqtc(ipop,iqtc)+kappaqtc(iqtc))
            vv = (sdqtctmp(ipop,iqtc)**2)/
     &           (nnqtc(ipop,iqtc)+kappaqtc(iqtc))
            meanqtctmp(ipop,iqtc) = ggrnorm(xx,dsqrt(vv))
            res = res - 
     &           ggdnorm(meanqtctmp(ipop,iqtc),xx,dsqrt(vv),1)
         enddo 
      enddo
      end subroutine contpropdir4



*********************************************************
*     contrib to log MH ratio of prior of mu and sigma 
*     in joint update of c, mu, sigma
*     prior: NSIG
      subroutine contprior(ipop1,ipop2,npopmax,nqtc,
     &     meanqtc,sdqtc,meanqtctmp,sdqtctmp,
     &     ksiqtc,kappaqtc,alphaqtc,betaqtc,res)
      implicit none
      integer ipop1,ipop2,npopmax,nqtc
      double precision meanqtc,sdqtc,meanqtctmp,sdqtctmp,
     &     ksiqtc,kappaqtc,alphaqtc,betaqtc,res
      integer iqtc
      double precision ggdnorm,ggdgamma
      dimension meanqtc(npopmax,nqtc),sdqtc(npopmax,nqtc),
     &     meanqtctmp(npopmax,nqtc),sdqtctmp(npopmax,nqtc),
     &     ksiqtc(nqtc),kappaqtc(nqtc),alphaqtc(nqtc),betaqtc(nqtc)
     
c$$$      write(*,*) 'contprior: meanqtc=',sdqtc
c$$$      write(*,*) 'contprior: meanqtctmp=',sdqtctmp
c$$$      write(*,*) 'contprior: sdqtc=',sdqtc
c$$$      write(*,*) 'contprior: sdqtctmp=',sdqtctmp
c$$$      write(*,*) 'ksiqtc=',ksiqtc
c$$$      write(*,*) 'kappaqtc=',kappaqtc
c$$$      write(*,*) 'alphaqtc=',alphaqtc
c$$$      write(*,*) 'betaqtc=',betaqtc

      res = 0 
      do iqtc = 1,nqtc
*     contrib prior of means to log ratio
         res = res + 
     &        ggdnorm(meanqtctmp(ipop1,iqtc),
     &        ksiqtc(iqtc),
     &        sdqtctmp(ipop1,iqtc)/dsqrt(kappaqtc(iqtc)),1) + 
     &        ggdnorm(meanqtctmp(ipop2,iqtc),
     &        ksiqtc(iqtc),
     &        sdqtctmp(ipop2,iqtc)/dsqrt(kappaqtc(iqtc)),1)  - 
     &        ggdnorm(meanqtc(ipop1,iqtc),
     &        ksiqtc(iqtc),
     &        sdqtc(ipop1,iqtc)/dsqrt(kappaqtc(iqtc)),1)  - 
     &        ggdnorm(meanqtc(ipop2,iqtc),
     &        ksiqtc(iqtc),
     &        sdqtc(ipop2,iqtc)/dsqrt(kappaqtc(iqtc)),1) 
      enddo
c$$$      write(*,*) 'contprior: res=',res

c$$$      junk = res
      do iqtc = 1,nqtc
*     contrib prior of variances to log ratio
         res = res +
     &        ggdgamma(1/sdqtctmp(ipop1,iqtc)**2,
     &        alphaqtc(iqtc),1/betaqtc(iqtc),1) + 
     &        ggdgamma(1/sdqtctmp(ipop2,iqtc)**2,
     &        alphaqtc(iqtc),1/betaqtc(iqtc),1) -
     &        ggdgamma(1/sdqtc(ipop1,iqtc)**2,
     &        alphaqtc(iqtc),1/betaqtc(iqtc),1) -
     &        ggdgamma(1/sdqtc(ipop2,iqtc)**2,
     &        alphaqtc(iqtc),1/betaqtc(iqtc),1)

c$$$         write(*,*) 'contprior: res=',res
c$$$ 

c$$$         write(*,*) 'contrib prior variances:'
c$$$         write(*,*) 'ggdgamma(1/sdqtctmp(ipop1,iqtc)**2,
c$$$     & alphaqtc(iqtc),1/betaqtc(iqtc),1=',
c$$$     &     ggdgamma(1/sdqtctmp(ipop1,iqtc)**2,
c$$$     &        alphaqtc(iqtc),1/betaqtc(iqtc),1)
c$$$         write(*,*) 'ggdgamma(1/sdqtctmp(ipop2,iqtc)**2,
c$$$     & alphaqtc(iqtc),1/betaqtc(iqtc),1)=',
c$$$     &        ggdgamma(1/sdqtctmp(ipop2,iqtc)**2,
c$$$     &        alphaqtc(iqtc),1/betaqtc(iqtc),1)
c$$$      write(*,*) 'ggdgamma(1/sdqtc(ipop1,iqtc)**2,
c$$$     & alphaqtc(iqtc),1/betaqtc(iqtc),1)=',
c$$$     &        ggdgamma(1/sdqtc(ipop1,iqtc)**2,
c$$$     &        alphaqtc(iqtc),1/betaqtc(iqtc),1)
c$$$      write(*,*) 'ggdgamma(1/sdqtc(ipop2,iqtc)**2,
c$$$     & alphaqtc(iqtc),1/betaqtc(iqtc),1)=',
c$$$     &     ggdgamma(1/sdqtc(ipop2,iqtc)**2,
c$$$     &     alphaqtc(iqtc),1/betaqtc(iqtc),1)
c$$$      write(*,*) 'global contrib prior variances=',res-junk
      enddo      
      end subroutine contprior





*********************************************************
*     contrib to log MH ratio of prior of mu and sigma 
*     in joint update of c, mu, sigma
*     prior: NSIG
      subroutine contprior4(npop,npopmax,nqtc,
     &     meanqtc,sdqtc,meanqtctmp,sdqtctmp,
     &     ksiqtc,kappaqtc,alphaqtc,betaqtc,res)
      implicit none
      integer npop,npopmax,nqtc
      double precision meanqtc,sdqtc,meanqtctmp,sdqtctmp,
     &     ksiqtc,kappaqtc,alphaqtc,betaqtc,res
      integer iqtc,ipop
      double precision ggdnorm,ggdgamma
      dimension meanqtc(npopmax,nqtc),sdqtc(npopmax,nqtc),
     &     meanqtctmp(npopmax,nqtc),sdqtctmp(npopmax,nqtc),
     &     ksiqtc(nqtc),kappaqtc(nqtc),alphaqtc(nqtc),betaqtc(nqtc)
      res = 0 
      do ipop=1,npop
         do iqtc = 1,nqtc
*     contrib prior of means to log ratio
            res = res + 
     &           ggdnorm(meanqtctmp(ipop,iqtc),
     &           ksiqtc(iqtc),
     &           sdqtctmp(ipop,iqtc)/dsqrt(kappaqtc(iqtc)),1) -
     &           ggdnorm(meanqtc(ipop,iqtc),
     &           ksiqtc(iqtc),
     &           sdqtc(ipop,iqtc)/dsqrt(kappaqtc(iqtc)),1)  
         enddo
         do iqtc = 1,nqtc
*     contrib prior of variances to log ratio
            res = res +
     &           ggdgamma(1/sdqtctmp(ipop,iqtc)**2,
     &           alphaqtc(iqtc),1/betaqtc(iqtc),1) -
     &           ggdgamma(1/sdqtc(ipop,iqtc)**2,
     &           alphaqtc(iqtc),1/betaqtc(iqtc),1) 
         enddo      
      enddo
      end subroutine contprior4


************************************************************************
*     propose mean and variance of quantitative variables in a split 
*     and returns contrib of prior and proposal to log ratio
      subroutine propqvsplit3(qtc,nindiv,nqtc,indcell,c,ctmp,
     &     meanqtc,sdqtc,meanqtctmp,sdqtctmp,nnqtc,
     &     sqtc,ssqtc,npop,npopmax,nppmax,
     &     ksiqtc,kappaqtc,alphaqtc,betaqtc,isplit,contriblr)
      implicit none
      integer nindiv,nqtc,indcell,c,ctmp,nnqtc,npop,npopmax,nppmax,
     &     isplit
      double precision qtc,meanqtc,sdqtc,
     &     meanqtctmp,sdqtctmp,sqtc,ssqtc,
     &     ksiqtc,kappaqtc,alphaqtc,betaqtc,contriblr
      dimension qtc(nindiv,nqtc),indcell(nindiv),c(nppmax),
     &     ctmp(nppmax),meanqtc(npopmax,nqtc),sdqtc(npopmax,nqtc),
     &     meanqtctmp(npopmax,nqtc),sdqtctmp(npopmax,nqtc),
     &     nnqtc(npopmax,nqtc),sqtc(npopmax,nqtc),
     &     ssqtc(npopmax,nqtc),
     &     ksiqtc(nqtc),kappaqtc(nqtc),alphaqtc(nqtc),betaqtc(nqtc)
      integer ipop,iindiv,iqtc
      double precision xx,vv,ap,bp,ggrnorm,ggrgam,ggdnorm,ggdgamma,junk

      contriblr = 0
*     compute empirical sums and sums of squares of quant. variables
*     for current state 
      do iqtc = 1,nqtc
         do ipop = 1,npopmax
            nnqtc(ipop,iqtc) = 0
            sqtc(ipop,iqtc) = 0
            ssqtc(ipop,iqtc) = 0
            meanqtctmp(ipop,iqtc) = meanqtc(ipop,iqtc)
            sdqtctmp(ipop,iqtc) = sdqtc(ipop,iqtc)
         enddo
      enddo
      do iqtc = 1,nqtc
         do iindiv = 1,nindiv
            ipop = c(indcell(iindiv))
            if(qtc(iindiv,iqtc) .ne. -999) then
               nnqtc(ipop,iqtc) = nnqtc(ipop,iqtc) + 1
               sqtc(ipop,iqtc) = sqtc(ipop,iqtc) + 
     &              qtc(iindiv,iqtc)
               ssqtc(ipop,iqtc) = ssqtc(ipop,iqtc) + 
     &              qtc(iindiv,iqtc)**2
            endif
         enddo
      enddo
      do iqtc = 1,nqtc
         do ipop = 1,npopmax
            if(nnqtc(ipop,iqtc) .ne. 0) then
               ssqtc(ipop,iqtc) = ssqtc(ipop,iqtc) -
     &              (sqtc(ipop,iqtc)**2) / dble(nnqtc(ipop,iqtc))
            else
               ssqtc(ipop,iqtc) = 0
            endif
         enddo
      enddo


*     contrib proposal pop isplit reverse move to log ratio
      do iqtc = 1,nqtc
*     mean 
         xx = (sqtc(isplit,iqtc)+ksiqtc(iqtc)*kappaqtc(iqtc))/
     &        (nnqtc(isplit,iqtc)+kappaqtc(iqtc))
         vv = sdqtc(isplit,iqtc)**2 / 
     &        (nnqtc(isplit,iqtc)+kappaqtc(iqtc))
         contriblr = contriblr +
     &        ggdnorm(meanqtc(isplit,iqtc),xx,dsqrt(vv),1)
*     variance
         ap = alphaqtc(iqtc) + 0.5*nnqtc(isplit,iqtc)
         if(nnqtc(isplit,iqtc) .ne. 0) then
            bp = betaqtc(iqtc) + 0.5*ssqtc(isplit,iqtc) + 
     &           0.5*(nnqtc(isplit,iqtc)*kappaqtc(iqtc))/
     &           (nnqtc(isplit,iqtc)+kappaqtc(iqtc)) *
     &           (sqtc(isplit,iqtc)/dble(nnqtc(isplit,iqtc)) - 
     &           ksiqtc(iqtc))**2
         else 
            bp = betaqtc(iqtc)
         endif
         contriblr = contriblr + 
     &        ggdgamma(1/sdqtc(isplit,iqtc)**2,ap,1/bp,1)
      enddo

*     compute empirical  sums  and sums of squares of quant. variables
*     for proposed state
      do iqtc = 1,nqtc
         do ipop = 1,npopmax
            nnqtc(ipop,iqtc) = 0
            sqtc(ipop,iqtc) = 0
            ssqtc(ipop,iqtc) = 0
            meanqtctmp(ipop,iqtc) = meanqtc(ipop,iqtc)
            sdqtctmp(ipop,iqtc) = sdqtc(ipop,iqtc)
         enddo
      enddo
      do iqtc = 1,nqtc
         do iindiv = 1,nindiv
            ipop = ctmp(indcell(iindiv))
            if(qtc(iindiv,iqtc) .ne. -999) then
               nnqtc(ipop,iqtc) = nnqtc(ipop,iqtc) + 1
               sqtc(ipop,iqtc) = sqtc(ipop,iqtc) + 
     &              qtc(iindiv,iqtc)
               ssqtc(ipop,iqtc) = ssqtc(ipop,iqtc) + 
     &              qtc(iindiv,iqtc)**2
            endif
         enddo
      enddo
      do iqtc = 1,nqtc
         do ipop = 1,npopmax
            if(nnqtc(ipop,iqtc) .ne. 0) then
               ssqtc(ipop,iqtc) = ssqtc(ipop,iqtc) -
     &              (sqtc(ipop,iqtc)**2) / dble(nnqtc(ipop,iqtc))
            else
               ssqtc(ipop,iqtc) = 0
            endif
         enddo
      enddo


*     sample mean from the full conditional, prior is NSIG
*     and compute contrib proposal to log ratio
*     propose mean,var  pop isplit direct move 
      do iqtc = 1,nqtc
*     propose variance pop isplit direct move 
*     and compute contrib proposal to log ratio
*     sample variance from full conditional
         ap = alphaqtc(iqtc) + 0.5*nnqtc(isplit,iqtc)
         if(nnqtc(isplit,iqtc) .ne. 0) then 
            bp = betaqtc(iqtc) + 0.5*ssqtc(isplit,iqtc) + 
     &           0.5*(nnqtc(isplit,iqtc)*kappaqtc(iqtc))/
     &           (nnqtc(isplit,iqtc)+kappaqtc(iqtc)) *
     &           (sqtc(isplit,iqtc)/dble(nnqtc(isplit,iqtc)) - 
     &           ksiqtc(iqtc))**2
         else
            bp = betaqtc(iqtc)
         endif
         sdqtctmp(isplit,iqtc) = 1/dsqrt(ggrgam(ap,1/bp))
         contriblr = contriblr - 
     &        ggdgamma(1/sdqtctmp(isplit,iqtc)**2,ap,1/bp,1)
 
*     propose mean,variance pop npop+1 direct move 
*     and compute contrib proposal to log ratio
*     sample variance from full conditional
         ap = alphaqtc(iqtc) + 0.5*nnqtc(npop+1,iqtc)
         if(nnqtc(npop+1,iqtc) .ne. 0) then 
            bp = betaqtc(iqtc) + 0.5*ssqtc(npop+1,iqtc) + 
     &           0.5*(nnqtc(npop+1,iqtc)*kappaqtc(iqtc))/
     &           (nnqtc(npop+1,iqtc)+kappaqtc(iqtc)) *
     &           (sqtc(npop+1,iqtc)/dble(nnqtc(npop+1,iqtc)) - 
     &           ksiqtc(iqtc))**2
         else
            bp = betaqtc(iqtc)
         endif
         sdqtctmp(npop+1,iqtc) = 1/dsqrt(ggrgam(ap,1/bp))
         contriblr = contriblr - 
     &        ggdgamma(1/sdqtctmp(npop+1,iqtc)**2,ap,1/bp,1)
         
         xx = (sqtc(isplit,iqtc)+ksiqtc(iqtc)*kappaqtc(iqtc))/
     &        (nnqtc(isplit,iqtc)+kappaqtc(iqtc))
         vv = (sdqtctmp(isplit,iqtc)**2)/
     &        (nnqtc(isplit,iqtc)+kappaqtc(iqtc))
         meanqtctmp(isplit,iqtc) = ggrnorm(xx,dsqrt(vv))
         contriblr = contriblr - 
     &        ggdnorm(meanqtctmp(isplit,iqtc),xx,dsqrt(vv),1)

*     propose mean pop npop+1 direct move 
*     and compute contrib proposal to log ratio
         xx = (sqtc(npop+1,iqtc)+ksiqtc(iqtc)*kappaqtc(iqtc))/
     &        (nnqtc(npop+1,iqtc)+kappaqtc(iqtc))
         vv = (sdqtctmp(npop+1,iqtc)**2)/
     &        (nnqtc(npop+1,iqtc)+kappaqtc(iqtc))
         meanqtctmp(npop+1,iqtc) = ggrnorm(xx,dsqrt(vv))
         contriblr = contriblr - 
     &        ggdnorm(meanqtctmp(npop+1,iqtc),xx,dsqrt(vv),1)
      enddo 

*     contrib prior
      do iqtc = 1,nqtc
*     contrib prior of means to log ratio
         contriblr = contriblr + 
     &        ggdnorm(meanqtctmp(isplit,iqtc),ksiqtc(iqtc),
     &        sdqtctmp(isplit,iqtc)/dsqrt(kappaqtc(iqtc)),1) + 
     &        ggdnorm(meanqtctmp(npop+1,iqtc),ksiqtc(iqtc),
     &        sdqtctmp(npop+1,iqtc)/dsqrt(kappaqtc(iqtc)),1)  - 
     &        ggdnorm(meanqtc(isplit,iqtc),ksiqtc(iqtc),
     &        sdqtc(isplit,iqtc)/dsqrt(kappaqtc(iqtc)),1) 
*     contrib prior of variances to log ratio
         contriblr = contriblr +
     &        ggdgamma(1/sdqtctmp(isplit,iqtc)**2,
     &        alphaqtc(iqtc),1/betaqtc(iqtc),1) + 
     &        ggdgamma(1/sdqtctmp(npop+1,iqtc)**2,
     &        alphaqtc(iqtc),1/betaqtc(iqtc),1) -
     &        ggdgamma(1/sdqtc(isplit,iqtc)**2,
     &        alphaqtc(iqtc),1/betaqtc(iqtc),1)
      enddo
c$$$      write(*,*) ''
c$$$      write(*,*) 'fin propqvsplit3'
c$$$      write(*,*) 'meanqtctmp=',meanqtctmp
c$$$      write(*,*) 'sdqtctmp=',sdqtctmp
c$$$      write(*,*) 'contriblr=',contriblr
      end subroutine propqvsplit3


************************************************************************
*     propose mean and variance of quantitative variables in a split 
*     and returns contrib of prior and proposal to log ratio
      subroutine propparqvsplit(qtc,nindiv,nqtc,indcell,c,ctmp,
     &     meanqtc,sdqtc,meanqtctmp,sdqtctmp,nnqtc,
     &     sqtc,ssqtc,npop,npopmax,nppmax,
     &     ksiqtc,kappaqtc,alphaqtc,betaqtc,isplit,contriblr)
      implicit none
      integer nindiv,nqtc,indcell,c,ctmp,nnqtc,npop,npopmax,nppmax,
     &     isplit
      double precision qtc,meanqtc,sdqtc,
     &     meanqtctmp,sdqtctmp,sqtc,ssqtc,
     &     ksiqtc,kappaqtc,alphaqtc,betaqtc,contriblr
      dimension qtc(nindiv,nqtc),indcell(nindiv),c(nppmax),
     &     ctmp(nppmax),meanqtc(npopmax,nqtc),sdqtc(npopmax,nqtc),
     &     meanqtctmp(npopmax,nqtc),sdqtctmp(npopmax,nqtc),
     &     nnqtc(npopmax,nqtc),sqtc(npopmax,nqtc),
     &     ssqtc(npopmax,nqtc),
     &     ksiqtc(nqtc),kappaqtc(nqtc),alphaqtc(nqtc),betaqtc(nqtc)
      integer ipop,iindiv,iqtc
      double precision xx,vv,ap,bp,ggrnorm,ggrgam,ggdnorm,ggdgamma,junk

      contriblr = 0
*     compute empirical sums and sums of squares of quant. variables
*     for current state 
      do iqtc = 1,nqtc
         do ipop = 1,npopmax
            nnqtc(ipop,iqtc) = 0
            sqtc(ipop,iqtc) = 0
            ssqtc(ipop,iqtc) = 0
            meanqtctmp(ipop,iqtc) = meanqtc(ipop,iqtc)
            sdqtctmp(ipop,iqtc) = sdqtc(ipop,iqtc)
         enddo
      enddo
      do iqtc = 1,nqtc
         do iindiv = 1,nindiv
            ipop = c(indcell(iindiv))
            if(qtc(iindiv,iqtc) .ne. -999) then
               nnqtc(ipop,iqtc) = nnqtc(ipop,iqtc) + 1
               sqtc(ipop,iqtc) = sqtc(ipop,iqtc) + 
     &              qtc(iindiv,iqtc)
               ssqtc(ipop,iqtc) = ssqtc(ipop,iqtc) + 
     &              qtc(iindiv,iqtc)**2
            endif
         enddo
      enddo


*     contrib proposal pop isplit reverse move to log ratio
      do iqtc = 1,nqtc
*     mean 
         junk = contriblr
         xx = (sqtc(isplit,iqtc)+ksiqtc(iqtc)*kappaqtc(iqtc))/
     &        (nnqtc(isplit,iqtc)+kappaqtc(iqtc))
         vv = 1 / (nnqtc(isplit,iqtc)+kappaqtc(iqtc))
         contriblr = contriblr +
     &        ggdnorm(meanqtc(isplit,iqtc),xx,dsqrt(vv),1)
*     variance
         junk = contriblr
         ap = alphaqtc(iqtc) + 0.5*nnqtc(isplit,iqtc)
         bp = betaqtc(iqtc) + 0.5*(ssqtc(isplit,iqtc) - 
     &        2*sqtc(isplit,iqtc)*meanqtc(isplit,iqtc) + 
     &        nnqtc(isplit,iqtc)*meanqtc(isplit,iqtc)**2)
          contriblr = contriblr + 
     &        ggdgamma(1/sdqtc(isplit,iqtc)**2,ap,1/bp,1)
      enddo

*     compute empirical  sums  and sums of squares of quant. variables
*     for proposed state
      do iqtc = 1,nqtc
         do ipop = 1,npopmax
            nnqtc(ipop,iqtc) = 0
            sqtc(ipop,iqtc) = 0
            ssqtc(ipop,iqtc) = 0
            meanqtctmp(ipop,iqtc) = meanqtc(ipop,iqtc)
            sdqtctmp(ipop,iqtc) = sdqtc(ipop,iqtc)
         enddo
      enddo
      do iqtc = 1,nqtc
         do iindiv = 1,nindiv
            ipop = ctmp(indcell(iindiv))
            if(qtc(iindiv,iqtc) .ne. -999) then
               nnqtc(ipop,iqtc) = nnqtc(ipop,iqtc) + 1
               sqtc(ipop,iqtc) = sqtc(ipop,iqtc) + 
     &              qtc(iindiv,iqtc)
               ssqtc(ipop,iqtc) = ssqtc(ipop,iqtc) + 
     &              qtc(iindiv,iqtc)**2
            endif
         enddo
      enddo

*     sample mean from the full conditional
*     assuming current value for variance is 1
      do iqtc = 1,nqtc
*     propose mean pop isplit direct move 
*     and compute contrib proposal to log ratio
         junk = contriblr
         xx = (sqtc(isplit,iqtc)+ksiqtc(iqtc)*kappaqtc(iqtc))/
     &        (nnqtc(isplit,iqtc)+kappaqtc(iqtc))
         vv = 1 / (nnqtc(isplit,iqtc)+kappaqtc(iqtc))
         meanqtctmp(isplit,iqtc) = ggrnorm(xx,dsqrt(vv))
         contriblr = contriblr - 
     &        ggdnorm(meanqtctmp(isplit,iqtc),xx,dsqrt(vv),1)
*     propose mean pop npop+1 direct move 
*     and compute contrib proposal to log ratio
         xx = (sqtc(npop+1,iqtc)+ksiqtc(iqtc)*kappaqtc(iqtc))/
     &        (nnqtc(npop+1,iqtc)+kappaqtc(iqtc))
         vv = 1 / (nnqtc(npop+1,iqtc)+kappaqtc(iqtc))
         meanqtctmp(npop+1,iqtc) = ggrnorm(xx,dsqrt(vv))
         contriblr = contriblr - 
     &        ggdnorm(meanqtctmp(npop+1,iqtc),xx,dsqrt(vv),1)
*     propose variance pop isplit direct move 
*     and compute contrib proposal to log ratio
*     sample variance from full conditional
         ap = alphaqtc(iqtc) + 0.5*nnqtc(isplit,iqtc)
         bp = betaqtc(iqtc) + 0.5*(ssqtc(isplit,iqtc) - 
     &        2*sqtc(isplit,iqtc)*meanqtctmp(isplit,iqtc) + 
     &        nnqtc(isplit,iqtc)*meanqtctmp(isplit,iqtc)**2)
         sdqtctmp(isplit,iqtc) = 1/dsqrt(ggrgam(ap,1/bp))
         contriblr = contriblr - 
     &        ggdgamma(1/sdqtctmp(isplit,iqtc)**2,ap,1/bp,1)
*     propose variance pop npop+1 direct move 
*     and compute contrib proposal to log ratio
*     sample variance from full conditional
         ap = alphaqtc(iqtc) + 0.5*nnqtc(npop+1,iqtc)
         bp = betaqtc(iqtc) + 0.5*(ssqtc(npop+1,iqtc) - 
     &        2*sqtc(npop+1,iqtc)*meanqtctmp(npop+1,iqtc) + 
     &        nnqtc(npop+1,iqtc)*meanqtctmp(npop+1,iqtc)**2)
         sdqtctmp(npop+1,iqtc) = 1/dsqrt(ggrgam(ap,1/bp))
         contriblr = contriblr - 
     &        ggdgamma(1/sdqtctmp(npop+1,iqtc)**2,ap,1/bp,1)
      enddo 

      do iqtc = 1,nqtc
*     contrib prior of means to log ratio
         contriblr = contriblr + 
     &        ggdnorm(meanqtctmp(isplit,iqtc),
     &        ksiqtc(iqtc),1/dsqrt(kappaqtc(iqtc)),1) + 
     &        ggdnorm(meanqtctmp(npop+1,iqtc),
     &        ksiqtc(iqtc),1/dsqrt(kappaqtc(iqtc)),1)  - 
     &        ggdnorm(meanqtc(isplit,iqtc),
     &        ksiqtc(iqtc),1/dsqrt(kappaqtc(iqtc)),1) 
*     contrib prior of variances to log ratio
         contriblr = contriblr +
     &        ggdgamma(1/sdqtctmp(isplit,iqtc)**2,
     &        alphaqtc(iqtc),1/betaqtc(iqtc),1) + 
     &        ggdgamma(1/sdqtctmp(npop+1,iqtc)**2,
     &        alphaqtc(iqtc),1/betaqtc(iqtc),1) -
     &        ggdgamma(1/sdqtc(isplit,iqtc)**2,
     &        alphaqtc(iqtc),1/betaqtc(iqtc),1)
      enddo
      end subroutine propparqvsplit

c$$$
c$$$
c$$$************************************************************************
c$$$*     propose mean and variance of quantitative variables in a merge
c$$$*     and returns contrib of prior and proposal to log ratio
c$$$      subroutine propparqvmerge(qtc,nindiv,nqtc,indcell,c,ctmp,
c$$$     &     meanqtc,sdqtc,meanqtctmp,sdqtctmp,nnqtc,
c$$$     &     sqtc,ssqtc,npop,npopmax,nppmax,
c$$$     &     ksiqtc,kappaqtc,alphaqtc,betaqtc,ipophost,ipoprem,contriblr)
c$$$      implicit none
c$$$      integer nindiv,nqtc,indcell,c,ctmp,nnqtc,npop,npopmax,nppmax,
c$$$     &     ipophost,ipoprem
c$$$      double precision qtc,meanqtc,sdqtc,
c$$$     &     meanqtctmp,sdqtctmp,sqtc,ssqtc,
c$$$     &     ksiqtc,kappaqtc,alphaqtc,betaqtc,contriblr
c$$$      dimension qtc(nindiv,nqtc),indcell(nindiv),c(nppmax),
c$$$     &     ctmp(nppmax),meanqtc(npopmax,nqtc),sdqtc(npopmax,nqtc),
c$$$     &     meanqtctmp(npopmax,nqtc),sdqtctmp(npopmax,nqtc),
c$$$     &     nnqtc(npopmax,nqtc),sqtc(npopmax,nqtc),
c$$$     &     ssqtc(npopmax,nqtc)
c$$$      integer ipop,iindiv,iqtc
c$$$      double precision xx,vv,ap,bp,ggrnorm,ggrgam,ggdnorm,ggdgamma
c$$$
c$$$      contriblr = 0
c$$$*     compute empirical  sums  and sums of squares of quant. variables
c$$$*     for current state 
c$$$      do iqtc = 1,nqtc
c$$$         do ipop = 1,npopmax
c$$$            nnqtc(ipop,iqtc) = 0
c$$$            sqtc(ipop,iqtc) = 0
c$$$            ssqtc(ipop,iqtc) = 0
c$$$            meanqtctmp(ipop,iqtc) = meanqtc(ipop,iqtc)
c$$$            sdqtctmp(ipop,iqtc) = sdqtc(ipop,iqtc)
c$$$         enddo
c$$$      enddo
c$$$      do iqtc = 1,nqtc
c$$$         do iindiv = 1,nindiv
c$$$            ipop = c(indcell(iindiv))
c$$$            if(qtc(iindiv,iqtc) .ne. -999) then
c$$$               nnqtc(ipop,iqtc) = nnqtc(ipop,iqtc) + 1
c$$$               sqtc(ipop,iqtc) = sqtc(ipop,iqtc) + 
c$$$     &              qtc(iindiv,iqtc)
c$$$               ssqtc(ipop,iqtc) = ssqtc(ipop,iqtc) + 
c$$$     &              qtc(iindiv,iqtc)**2
c$$$            endif
c$$$         enddo
c$$$      enddo
c$$$
c$$$      do iqtc = 1,nqtc
c$$$*     contrib proposal pop ipophost to log ratio in reverse move 
c$$$*     mean 
c$$$         xx = (sqtc(ipophost,iqtc)+ksiqtc*kappaqtc)/
c$$$     &        (nnqtc(ipophost,iqtc)+kappaqtc)
c$$$         vv = 1 / (nnqtc(ipophost,iqtc)+kappaqtc)
c$$$         contriblr = contriblr +
c$$$     &        ggdnorm(meanqtc(ipophost,iqtc),xx,dsqrt(vv),1)
c$$$*     variance
c$$$         ap = alphaqtc + 0.5*nnqtc(ipophost,iqtc)
c$$$         bp = betaqtc + 0.5*(ssqtc(ipophost,iqtc) - 
c$$$     &        2*sqtc(ipophost,iqtc)*meanqtc(ipophost,iqtc) + 
c$$$     &        nnqtc(ipophost,iqtc)*meanqtc(ipophost,iqtc)**2)
c$$$         contriblr = contriblr + 
c$$$     &        ggdgamma(1/sdqtc(ipophost,iqtc)**2,ap,1/bp,1)
c$$$*     contrib proposal pop ipoprem to log ratio in reverse move 
c$$$*     mean 
c$$$         xx = (sqtc(ipoprem,iqtc)+ksiqtc*kappaqtc)/
c$$$     &        (nnqtc(ipoprem,iqtc)+kappaqtc)
c$$$         vv = 1 / (nnqtc(ipoprem,iqtc)+kappaqtc)
c$$$         contriblr = contriblr +
c$$$     &        ggdnorm(meanqtc(ipoprem,iqtc),xx,dsqrt(vv),1)
c$$$*     variance
c$$$         ap = alphaqtc + 0.5*nnqtc(ipoprem,iqtc)
c$$$         bp = betaqtc + 0.5*(ssqtc(ipoprem,iqtc) - 
c$$$     &        2*sqtc(ipoprem,iqtc)*meanqtc(ipoprem,iqtc) + 
c$$$     &        nnqtc(ipoprem,iqtc)*meanqtc(ipoprem,iqtc)**2)
c$$$         contriblr = contriblr + 
c$$$     &        ggdgamma(1/sdqtc(ipoprem,iqtc)**2,ap,1/bp,1)
c$$$      enddo
c$$$
c$$$
c$$$*     compute empirical  sums  and sums of squares of quant. variables
c$$$*     for proposed state
c$$$      do iqtc = 1,nqtc
c$$$         do ipop = 1,npopmax
c$$$            nnqtc(ipop,iqtc) = 0
c$$$            sqtc(ipop,iqtc) = 0
c$$$            ssqtc(ipop,iqtc) = 0
c$$$            meanqtctmp(ipop,iqtc) = meanqtc(ipop,iqtc)
c$$$            sdqtctmp(ipop,iqtc) = sdqtc(ipop,iqtc)
c$$$         enddo
c$$$      enddo
c$$$      do iqtc = 1,nqtc
c$$$         do iindiv = 1,nindiv
c$$$            ipop = ctmp(indcell(iindiv))
c$$$            if(qtc(iindiv,iqtc) .ne. -999) then
c$$$               nnqtc(ipop,iqtc) = nnqtc(ipop,iqtc) + 1
c$$$               sqtc(ipop,iqtc) = sqtc(ipop,iqtc) + 
c$$$     &              qtc(iindiv,iqtc)
c$$$               ssqtc(ipop,iqtc) = ssqtc(ipop,iqtc) + 
c$$$     &              qtc(iindiv,iqtc)**2
c$$$            endif
c$$$         enddo
c$$$      enddo
c$$$
c$$$*     sample mean from the full conditional
c$$$*     assuming current value for variance is 1
c$$$      do iqtc = 1,nqtc
c$$$*     propose mean pop ipophost direct move 
c$$$*     and compute contrib proposal to log ratio
c$$$         xx = (sqtc(ipophost,iqtc)+ksiqtc*kappaqtc)/
c$$$     &        (nnqtc(ipophost,iqtc)+kappaqtc)
c$$$         vv = 1 / (nnqtc(ipophost,iqtc)+kappaqtc)
c$$$         meanqtctmp(ipophost,iqtc) = ggrnorm(xx,dsqrt(vv))
c$$$         contriblr = contriblr - 
c$$$     &        ggdnorm(meanqtctmp(ipophost,iqtc),xx,dsqrt(vv),1)
c$$$*     propose variance pop ipophost direct move 
c$$$*     and compute contrib proposal to log ratio
c$$$*     sample variance from full conditional
c$$$         ap = alphaqtc + 0.5*nnqtc(ipophost,iqtc)
c$$$         bp = betaqtc + 0.5*(ssqtc(ipophost,iqtc) - 
c$$$     &        2*sqtc(ipophost,iqtc)*meanqtctmp(ipophost,iqtc) + 
c$$$     &        nnqtc(ipophost,iqtc)*meanqtctmp(ipophost,iqtc)**2)
c$$$         sdqtctmp(ipophost,iqtc) = 1/dsqrt(ggrgam(ap,1/bp))
c$$$         contriblr = contriblr - 
c$$$     &        ggdgamma(1/sdqtctmp(ipophost,iqtc)**2,ap,1/bp,1)
c$$$      enddo 
c$$$*     arrange values for other populations
c$$$      do iqtc = 1,nqtc
c$$$         if(ipoprem .ne. npop) then
c$$$            do ipop =ipoprem+1,npop
c$$$               meanqtctmp(ipoprem-1,iqtc) = meanqtc(ipoprem,iqtc)
c$$$               sdqtctmp(ipoprem-1,iqtc) = sdqtc(ipoprem,iqtc)
c$$$            enddo
c$$$         endif
c$$$         do ipop = npop,npopmax
c$$$            meanqtctmp(ipop,iqtc) = -999
c$$$            sdqtctmp(ipop,iqtc) = -999
c$$$         enddo
c$$$      enddo 
c$$$
c$$$      do iqtc = 1,nqtc
c$$$*     contrib prior of means to log ratio
c$$$         contriblr = contriblr + 
c$$$     &        ggdnorm(meanqtctmp(ipophost,iqtc),
c$$$     &        ksiqtc,1/dsqrt(kappaqtc),1) -
c$$$     &        ggdnorm(meanqtc(ipophost,iqtc),
c$$$     &        ksiqtc,1/dsqrt(kappaqtc),1) -
c$$$     &        ggdnorm(meanqtc(ipoprem,iqtc),
c$$$     &        ksiqtc,1/dsqrt(kappaqtc),1) 
c$$$*     contrib prior of variances to log ratio
c$$$         contriblr = contriblr +
c$$$     &        ggdgamma(1/sdqtctmp(ipophost,iqtc)**2,
c$$$     &        alphaqtc,1/betaqtc,1) - 
c$$$     &        ggdgamma(1/sdqtc(ipophost,iqtc)**2,
c$$$     &        alphaqtc,1/betaqtc,1) - 
c$$$     &        ggdgamma(1/sdqtc(ipoprem,iqtc)**2,
c$$$     &        alphaqtc,1/betaqtc,1)
c$$$      enddo
c$$$      end subroutine propparqvmerge



************************************************************************
*     propose mean and variance of quantitative variables in a merge
*     and returns contrib of prior and proposal to log ratio
      subroutine propqvmrg3(qtc,nindiv,nqtc,indcell,c,ctmp,
     &     meanqtc,sdqtc,meanqtctmp,sdqtctmp,nnqtc,
     &     sqtc,ssqtc,npop,npopmax,nppmax,
     &     ksiqtc,kappaqtc,alphaqtc,betaqtc,ipophost,ipoprem,contriblr)
      implicit none
      integer nindiv,nqtc,indcell,c,ctmp,nnqtc,npop,npopmax,nppmax,
     &     ipophost,ipoprem
      double precision qtc,meanqtc,sdqtc,
     &     meanqtctmp,sdqtctmp,sqtc,ssqtc,
     &     ksiqtc,kappaqtc,alphaqtc,betaqtc,contriblr
      dimension qtc(nindiv,nqtc),indcell(nindiv),c(nppmax),
     &     ctmp(nppmax),meanqtc(npopmax,nqtc),sdqtc(npopmax,nqtc),
     &     meanqtctmp(npopmax,nqtc),sdqtctmp(npopmax,nqtc),
     &     nnqtc(npopmax,nqtc),sqtc(npopmax,nqtc),
     &     ssqtc(npopmax,nqtc),ksiqtc(nqtc),kappaqtc(nqtc),
     &     alphaqtc(nqtc),betaqtc(nqtc)
      integer ipop,iindiv,iqtc
      double precision xx,vv,ap,bp,ggrnorm,ggrgam,ggdnorm,ggdgamma
c$$$      write(*,*) ''
c$$$      write(*,*) 'debut propqvmrg3'
c$$$      write(*,*) 'ipophost,ipoprem=',ipophost,ipoprem
c$$$      write(*,*) 'meanqtc=',meanqtc
c$$$      write(*,*) 'sdqtc=',sdqtc
      contriblr = 0
*     compute empirical  sums  and sums of squares of quant. variables
*     for current state 
      do iqtc = 1,nqtc
         do ipop = 1,npopmax
            nnqtc(ipop,iqtc) = 0
            sqtc(ipop,iqtc) = 0
            ssqtc(ipop,iqtc) = 0
            meanqtctmp(ipop,iqtc) = meanqtc(ipop,iqtc)
            sdqtctmp(ipop,iqtc) = sdqtc(ipop,iqtc)
         enddo
      enddo
      do iqtc = 1,nqtc
         do iindiv = 1,nindiv
            ipop = c(indcell(iindiv))
            if(qtc(iindiv,iqtc) .ne. -999) then
               nnqtc(ipop,iqtc) = nnqtc(ipop,iqtc) + 1
               sqtc(ipop,iqtc) = sqtc(ipop,iqtc) + 
     &              qtc(iindiv,iqtc)
               ssqtc(ipop,iqtc) = ssqtc(ipop,iqtc) + 
     &              qtc(iindiv,iqtc)**2
            endif
         enddo
      enddo
      do iqtc = 1,nqtc
         do ipop = 1,npopmax
            if(nnqtc(ipop,iqtc) .ne. 0) then
               ssqtc(ipop,iqtc) = ssqtc(ipop,iqtc) -
     &              (sqtc(ipop,iqtc)**2) / dble(nnqtc(ipop,iqtc))
            else
               ssqtc(ipop,iqtc) = 0
            endif
         enddo
      enddo

      do iqtc = 1,nqtc
*     contrib proposal pop ipophost to log ratio in reverse move 
*     mean 
         xx = (sqtc(ipophost,iqtc)+ksiqtc(iqtc)*kappaqtc(iqtc))/
     &        (nnqtc(ipophost,iqtc)+kappaqtc(iqtc))
         vv = sdqtc(ipophost,iqtc)**2 / 
     &        (nnqtc(ipophost,iqtc)+kappaqtc(iqtc))
         contriblr = contriblr +
     &        ggdnorm(meanqtc(ipophost,iqtc),xx,dsqrt(vv),1)
*     variance
         ap = alphaqtc(iqtc) + 0.5*nnqtc(ipophost,iqtc)
         if(nnqtc(ipophost,iqtc) .ne. 0) then 
            bp = betaqtc(iqtc) + 0.5*ssqtc(ipophost,iqtc) + 
     &           0.5*(nnqtc(ipophost,iqtc)*kappaqtc(iqtc))/
     &           (nnqtc(ipophost,iqtc)+kappaqtc(iqtc)) *
     &           (sqtc(ipophost,iqtc)/dble(nnqtc(ipophost,iqtc)) - 
     &           ksiqtc(iqtc))**2
         else
            bp = betaqtc(iqtc)
         endif
         contriblr = contriblr + 
     &        ggdgamma(1/sdqtc(ipophost,iqtc)**2,ap,1/bp,1)

*     contrib proposal pop ipoprem to log ratio in reverse move 
*     mean 
         xx = (sqtc(ipoprem,iqtc)+ksiqtc(iqtc)*kappaqtc(iqtc))/
     &        (nnqtc(ipoprem,iqtc)+kappaqtc(iqtc))
         vv = sdqtc(ipoprem,iqtc)**2 / 
     &        (nnqtc(ipoprem,iqtc)+kappaqtc(iqtc))
         contriblr = contriblr +
     &        ggdnorm(meanqtc(ipoprem,iqtc),xx,dsqrt(vv),1)
*     variance
         ap = alphaqtc(iqtc) + 0.5*nnqtc(ipoprem,iqtc)
         if(nnqtc(ipoprem,iqtc) .ne. 0) then 
            bp = betaqtc(iqtc) + 0.5*ssqtc(ipoprem,iqtc) + 
     &           0.5*(nnqtc(ipoprem,iqtc)*kappaqtc(iqtc))/
     &           (nnqtc(ipoprem,iqtc)+kappaqtc(iqtc)) *
     &           (sqtc(ipoprem,iqtc)/dble(nnqtc(ipoprem,iqtc)) - 
     &           ksiqtc(iqtc))**2
         else
            bp = betaqtc(iqtc)
         endif
         contriblr = contriblr + 
     &        ggdgamma(1/sdqtc(ipoprem,iqtc)**2,ap,1/bp,1)
      enddo

*     compute empirical  sums  and sums of squares of quant. variables
*     for proposed state
      do iqtc = 1,nqtc
         do ipop = 1,npopmax
            nnqtc(ipop,iqtc) = 0
            sqtc(ipop,iqtc) = 0
            ssqtc(ipop,iqtc) = 0
            meanqtctmp(ipop,iqtc) = meanqtc(ipop,iqtc)
            sdqtctmp(ipop,iqtc) = sdqtc(ipop,iqtc)
         enddo
      enddo
      do iqtc = 1,nqtc
         do iindiv = 1,nindiv
            ipop = ctmp(indcell(iindiv))
            if(qtc(iindiv,iqtc) .ne. -999) then
               nnqtc(ipop,iqtc) = nnqtc(ipop,iqtc) + 1
               sqtc(ipop,iqtc) = sqtc(ipop,iqtc) + 
     &              qtc(iindiv,iqtc)
               ssqtc(ipop,iqtc) = ssqtc(ipop,iqtc) + 
     &              qtc(iindiv,iqtc)**2
            endif
         enddo
      enddo
      do iqtc = 1,nqtc
         do ipop = 1,npopmax
            if(nnqtc(ipop,iqtc) .ne. 0) then
               ssqtc(ipop,iqtc) = ssqtc(ipop,iqtc) -
     &              (sqtc(ipop,iqtc)**2) / dble(nnqtc(ipop,iqtc))
            else
               ssqtc(ipop,iqtc) = 0
            endif
         enddo
      enddo

*     sample mean from the full conditional
      do iqtc = 1,nqtc
*     propose variance pop ipophost direct move 
*     and compute contrib proposal to log ratio
*     sample variance from full conditional
         ap = alphaqtc(iqtc) + 0.5*nnqtc(ipophost,iqtc)
         if(nnqtc(ipophost,iqtc) .ne. 0) then 
            bp = betaqtc(iqtc) + 0.5*ssqtc(ipophost,iqtc) + 
     &           0.5*(nnqtc(ipophost,iqtc)*kappaqtc(iqtc))/
     &           (nnqtc(ipophost,iqtc)+kappaqtc(iqtc)) *
     &           (sqtc(ipophost,iqtc)/dble(nnqtc(ipophost,iqtc)) - 
     &           ksiqtc(iqtc))**2
         else
            bp = betaqtc(iqtc)
         endif
         sdqtctmp(ipophost,iqtc) = 1/dsqrt(ggrgam(ap,1/bp))
         contriblr = contriblr - 
     &        ggdgamma(1/sdqtctmp(ipophost,iqtc)**2,ap,1/bp,1)
*     propose mean pop ipophost direct move 
*     and compute contrib proposal to log ratio
         xx = (sqtc(ipophost,iqtc)+ksiqtc(iqtc)*kappaqtc(iqtc))/
     &        (nnqtc(ipophost,iqtc)+kappaqtc(iqtc))
         vv = sdqtctmp(ipophost,iqtc)**2 
     &        / (nnqtc(ipophost,iqtc)+kappaqtc(iqtc))
         meanqtctmp(ipophost,iqtc) = ggrnorm(xx,dsqrt(vv))
         contriblr = contriblr - 
     &        ggdnorm(meanqtctmp(ipophost,iqtc),xx,dsqrt(vv),1)

      enddo 
*     arrange values for other populations
      do iqtc = 1,nqtc
         if(ipoprem .ne. npop) then
            do ipop =ipoprem+1,npop
               meanqtctmp(ipoprem-1,iqtc) = meanqtc(ipoprem,iqtc)
               sdqtctmp(ipoprem-1,iqtc) = sdqtc(ipoprem,iqtc)
            enddo
         endif
         do ipop = npop,npopmax
            meanqtctmp(ipop,iqtc) = -999
            sdqtctmp(ipop,iqtc) = -999
         enddo
      enddo 


      do iqtc = 1,nqtc
*     contrib prior of means to log ratio
         contriblr = contriblr + 
     &        ggdnorm(meanqtctmp(ipophost,iqtc),ksiqtc(iqtc),
     &        sdqtctmp(ipophost,iqtc)/dsqrt(kappaqtc(iqtc)),1) -
     &        ggdnorm(meanqtc(ipophost,iqtc),ksiqtc(iqtc),
     &        sdqtc(ipophost,iqtc)/dsqrt(kappaqtc(iqtc)),1) -
     &        ggdnorm(meanqtc(ipoprem,iqtc),ksiqtc(iqtc),
     &        sdqtc(ipoprem,iqtc)/dsqrt(kappaqtc(iqtc)),1) 
*     contrib prior of variances to log ratio
         contriblr = contriblr +
     &        ggdgamma(1/sdqtctmp(ipophost,iqtc)**2,
     &        alphaqtc,1/betaqtc(iqtc),1) - 
     &        ggdgamma(1/sdqtc(ipophost,iqtc)**2,
     &        alphaqtc,1/betaqtc(iqtc),1) - 
     &        ggdgamma(1/sdqtc(ipoprem,iqtc)**2,
     &        alphaqtc,1/betaqtc(iqtc),1)
      enddo
c$$$      write(*,*) 'fin propqvmrg3'
c$$$      write(*,*) 'meanqtctmp=',meanqtctmp
c$$$      write(*,*) 'sdqtctmp=',sdqtctmp
c$$$      write(*,*) 'contriblr=',contriblr
      end subroutine propqvmrg3


c$$$
c$$$************************************************************************
c$$$*     propose mean and variance of quantitative variables in a merge
c$$$      subroutine propparqvmrg(qtc,nindiv,nqtc,indcell,ctmp,
c$$$     &     meanqtctmp,sdqtctmp,nnqtc,sqtc,ssqtc,npop,npopmax,
c$$$     &     nppmax,ksiqtc,kappaqtc,alphaqtc,betaqtc,ihost)
c$$$      implicit none
c$$$      integer nindiv,nqtc,indcell,ctmp,nnqtc,npop,npopmax,nppmax,
c$$$     &     ihost
c$$$      double precision qtc,meanqtctmp,sdqtctmp,sqtc,ssqtc,
c$$$     &     ksiqtc,kappaqtc,alphaqtc,betaqtc
c$$$      dimension qtc(nindiv,nqtc),indcell(nindiv),ctmp(nppmax),
c$$$     &     meanqtctmp(npopmax,nqtc),sdqtctmp(npopmax,nqtc),
c$$$     &     nnqtc(npopmax,nqtc),sqtc(npopmax,nqtc),
c$$$     &     ssqtc(npopmax,nqtc)
c$$$      integer ipop,iindiv,iqtc
c$$$      double precision xx,vv,ap,bp,ggrnorm,ggrgam
c$$$*     compute empirical  sums  and sums of squares for quant. variables
c$$$      do iqtc = 1,nqtc
c$$$            nnqtc(ihost,iqtc) = 0
c$$$            sqtc(ihost,iqtc) = 0
c$$$            ssqtc(ihost,iqtc) = 0
c$$$      enddo
c$$$      do iqtc = 1,nqtc
c$$$         do iindiv = 1,nindiv
c$$$            ipop = ctmp(indcell(iindiv))
c$$$            if(qtc(iindiv,iqtc) .ne. -999) then
c$$$               nnqtc(ipop,iqtc) = nnqtc(ipop,iqtc) + 1
c$$$               sqtc(ipop,iqtc) = sqtc(ipop,iqtc) + 
c$$$     &              qtc(iindiv,iqtc)
c$$$               ssqtc(ipop,iqtc) = ssqtc(ipop,iqtc) + 
c$$$     &              qtc(iindiv,iqtc)**2
c$$$            endif
c$$$         enddo
c$$$      enddo
c$$$*     sample mean from the full conditional
c$$$*     assuming current value for variance is 1
c$$$      do iqtc = 1,nqtc
c$$$         xx = (sqtc(ihost,iqtc)+ksiqtc*kappaqtc)/
c$$$     &        (nnqtc(ihost,iqtc)+kappaqtc)
c$$$         vv = 1 / (nnqtc(ihost,iqtc)+kappaqtc)
c$$$         meanqtctmp(ihost,iqtc) = ggrnorm(xx,dsqrt(vv))
c$$$*     sample variance from full conditional
c$$$         ap = alphaqtc + 0.5*nnqtc(ihost,iqtc)
c$$$         bp = betaqtc + 0.5*(ssqtc(ihost,iqtc) - 
c$$$     &        2*sqtc(ihost,iqtc)*meanqtctmp(ihost,iqtc) + 
c$$$     &        nnqtc(ihost,iqtc)*meanqtctmp(ihost,iqtc)**2)
c$$$         sdqtctmp(ihost,iqtc) = 1/dsqrt(ggrgam(ap,1/bp))
c$$$      enddo 
c$$$      end subroutine propparqvmrg

************************************************************************
*     split and merge of populations with acceptance according to MH ratio
*     CFM
*     d* drawn form prior
      subroutine sm(npop,npopmin,npopmax,f,fa,drift,
     &     nloc,nloc2,
     &     nal,nalmax,indcell,nindiv,npp,nppmax,c,ctmp,
     &     a,ptmp,ftmp,drifttmp,zz,cellpop,listcell,
     &     cellpophost,n,ntmp,ploidy,shape1,shape2)
      implicit none 
      integer npop,npopmin,npopmax,nloc,nal(nloc),
     &     nalmax,nindiv,npp,nppmax,indcell(nindiv),
     &     nloc2,c(nppmax),ctmp(nppmax),zz(nindiv,nloc2),
     &     n(npopmax,nloc,nalmax),ntmp(npopmax,nloc,nalmax),
     &     ploidy
      double precision f(npopmax,nloc,nalmax),drift(npopmax),
     &      ftmp(npopmax,nloc,nalmax),drifttmp(npopmax),
     &     a(nalmax),ptmp(nalmax),fa(nloc,nalmax),shape1,shape2
      integer ipoprem,ipp,isplit,ipop,
     &     cellpop(nppmax),ncellpop,nu,listcell(nppmax),
     &     ipophost,ncellpophost,cellpophost(nppmax),ii
      double precision alpha,ggrunif,lbico,lratio,lTf
      integer iloc
      double precision gglgamfn,ggrbet,rr
      double precision bern,b,ggrbinom

c         write(*,*) 'debut sm' 
      do ipp=1,nppmax
         cellpop(ipp) = -999
         listcell(ipp) = -999
      enddo
*     naissance ou mort ?
      b = ggrbinom(1.d0,0.5d0)
      
      if(b .eq. 1) then
         if(npop .lt. npopmax) then 
c            write(*,*) 'naissance'
*     split
            
*     choix de la pop qui split
            isplit = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
            
*     recherche des cellules affectees a cette pop
            call who(c,isplit,npp,nppmax,cellpop,ncellpop)
            
            if(ncellpop .gt. 0) then
*     tirage du nombre de cellules reallouees
C     ligne suivante corrigee le 17/01/08
C     nu = idint(dint(dble(ncellpop)*ggrunif(0.d0,1.d0)))
               nu = idint(dint(dble(ncellpop+1)*ggrunif(0.d0,1.d0)))
               if(nu .gt. 0) then
*     tirage des cellules reallouees
                  call sample2(cellpop,nppmax,nu,ncellpop,
     &                 listcell)
*     proposition de reallocation dans la pop npop+1
                  call split(npop+1,c,ctmp,nppmax,nu,listcell)
               else 
                  do ipp = 1,nppmax
                     ctmp(ipp) = c(ipp)
                  enddo
               endif
            else
               nu = 0
               do ipp = 1,nppmax
                  ctmp(ipp) = c(ipp)
               enddo
            endif
*     comptage des alleles sur chaque locus pour c puis ctmp
            call countn(nindiv,nloc,nloc2,npopmax,
     &           nppmax,nal,nalmax,zz,n,indcell,c,ploidy)
            call countn(nindiv,nloc,nloc2,npopmax,
     &           nppmax,nal,nalmax,zz,ntmp,indcell,ctmp,ploidy)
*     nouvelles derives      
            do ipop = 1,npop
               drifttmp(ipop) = drift(ipop)
            enddo
            drifttmp(isplit) = ggrbet(shape1,shape2)
            drifttmp(npop+1) = ggrbet(shape1,shape2)
*     terme des proposal sur c
            lratio = dlog(2*dble(ncellpop+1)) + 
     &           lbico(ncellpop,nu) - dlog(dble(npop+1)) 
*     terme des priors sur c
            lratio = lratio + dble(npp)*(dlog(dble(npop)) - 
     &           dlog(dble(npop+1)))
*     terme des frequences
            lratio = lratio 
     &         + lTf(isplit,ntmp,fa,drifttmp,npopmax,nloc,nal,nalmax)  
     &         + lTf(npop+1,ntmp,fa,drifttmp,npopmax,nloc,nal,nalmax) 
     &         - lTf(isplit,n,fa,drift,npopmax,nloc,nal,nalmax) 
c            write(*,*) 'lratio=',lratio
            lratio = dmin1(0.d0,lratio)
            alpha = dexp(lratio)
            bern = ggrbinom(1.d0,alpha)
            if(bern .eq. 1) then
*     proposition nouvelle freq 
               call addfreq8(npop,npopmax,nloc,nloc,
     &              nal,nalmax,isplit,
     &              f,ftmp,fa,drift,drifttmp,a,ptmp,ntmp)
               call accept5(nppmax,npopmax,nloc,
     &              nalmax,nal,c,ctmp,f,ftmp,drift,drifttmp)
               npop = npop + 1
            endif
         endif
         
*     merge
      else
         if(npop .gt. npopmin) then 
c            write(*,*) 'mort'
*     tirage de la pop qui meurt
            ipoprem = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
*     tirage de la pop hote
            ipophost = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
            do while(ipophost .eq. ipoprem)
               ipophost = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
            enddo
*     on range dans la pop d'indice le plus petit
            if(ipophost .gt. ipoprem) then
               ii = ipophost
               ipophost = ipoprem
               ipoprem = ii
            endif
*     recherche des cellules qui vont etre reallouees
            call who(c,ipoprem,npp,nppmax,cellpop,ncellpop)
*     recherche des cellules de la pop hote
            call who(c,ipophost,npp,nppmax,cellpophost,
     &           ncellpophost)
*     proposition de reallocation dans la pop ipophost
            call merging(ipoprem,ipophost,c,ctmp,nppmax,
     &           ncellpop,cellpop)
*     comptage des alleles sur chaque locus pour c puis ctmp
            call countn(nindiv,nloc,nloc2,npopmax,
     &           nppmax,nal,nalmax,zz,n,indcell,c,ploidy)
            call countn(nindiv,nloc,nloc2,npopmax,
     &           nppmax,nal,nalmax,zz,ntmp,indcell,ctmp,ploidy)
*     enleve une pop dans le tableau tmporaire des derives            
            call remdrift(ipoprem,ipophost,npop,npopmax,drift,
     &           drifttmp,shape1,shape2)
*     calcul du log du ratio  
*     terme des proposal sur c
            lratio = dlog(dble(npop)) - 
     &           dlog(2*dble(ncellpop+ncellpophost+1)) -
     &           lbico(ncellpop+ncellpophost,ncellpop) 
*     terme des priors sur c
            lratio = lratio + 
     &           dble(npp)*(dlog(dble(npop)) - 
     &           dlog(dble(npop-1)))
c     write(*,*) 'term en c lratio =',lratio
*     term des freq
            lratio = lratio + 
     &         lTf(ipophost,ntmp,fa,drifttmp,npopmax,nloc,nal,nalmax) - 
     &         lTf(ipophost,n,fa,drift,npopmax,nloc,nal,nalmax) -
     &        lTf(ipoprem,n,fa,drift,npopmax,nloc,nal,nalmax) 
c$$$            write(*,*) 'term en f=',
c$$$     &          lTf(ipophost,ntmp,fa,drifttmp,npopmax,nloc,nal,nalmax)  
c$$$     &        - lTf(ipophost,n,fa,drift,npopmax,nloc,nal,nalmax) 
c$$$     &        - lTf(ipoprem,n,fa,drift,npopmax,nloc,nal,nalmax) 
c$$$            write(*,*) 'drifttmp(ipophost)=',drifttmp(ipophost)
c$$$            write(*,*) 'term en f*_k1=',
c$$$     &          lTf(ipophost,ntmp,fa,drifttmp,npopmax,nloc,nal,nalmax) 
c$$$            write(*,*) 'drift(ipophost)=',drift(ipophost)
c$$$             write(*,*) 'term en f_k1=',
c$$$     &        - lTf(ipophost,n,fa,drift,npopmax,nloc,nal,nalmax) 
c$$$             write(*,*) 'drift(ipoprem)=',drift(ipoprem)
c$$$            write(*,*) 'term en f_k2=',
c$$$     &        - lTf(ipoprem,n,fa,drift,npopmax,nloc,nal,nalmax) 
c$$$            write(*,*) 'lratio =',lratio
c            write(*,*) 'lratio=',lratio
            lratio = dmin1(0.d0,lratio)
            alpha  = dexp(lratio)
            bern = ggrbinom(1.d0,alpha)      
                          
            if(bern .eq. 1) then
*     propostion du nouveau tableau de freq et de derives
               call remfreq8(ipoprem,ipophost,
     &              npop,npopmax,nloc,nloc,nal,
     &              nalmax,f,ftmp,drift,drifttmp,fa,a,ptmp,ntmp)
               call accept5(nppmax,npopmax,nloc,
     &              nalmax,nal,c,ctmp,f,ftmp,drift,drifttmp)
               npop = npop - 1
            endif
         endif
      endif
      end subroutine sm
************************************************************************      



************************************************************************
*     split and merge of populations with acceptance according to MH ratio
*     CFM
*     (d1*,d2*) = d1-u,d2+u
      subroutine sm2(npop,npopmin,npopmax,f,fa,drift,
     &     nloc,nloc2,
     &     nal,nalmax,indcell,nindiv,npp,nppmax,c,ctmp,
     &     a,ptmp,ftmp,drifttmp,zz,cellpop,listcell,
     &     cellpophost,n,ntmp,ploidy,shape1,shape2)
      implicit none 
      integer npop,npopmin,npopmax,nloc,nal(nloc),
     &     nalmax,nindiv,npp,nppmax,indcell(nindiv),
     &     nloc2,c(nppmax),ctmp(nppmax),zz(nindiv,nloc2),
     &     n(npopmax,nloc,nalmax),ntmp(npopmax,nloc,nalmax),
     &     ploidy
      double precision f(npopmax,nloc,nalmax),drift(npopmax),
     &      ftmp(npopmax,nloc,nalmax),drifttmp(npopmax),
     &     a(nalmax),ptmp(nalmax),fa(nloc,nalmax),shape1,shape2
      integer ipoprem,ipp,isplit,ipop,
     &     cellpop(nppmax),ncellpop,nu,listcell(nppmax),
     &     ipophost,ncellpophost,cellpophost(nppmax),ii
      double precision alpha,ggrunif,lbico,lratio,lTf
      integer iloc
      double precision gglgamfn,rr,deltad
      double precision bern,b,ggrbinom,ggrnorm

      deltad = shape1/(shape1+shape2)

c         write(*,*) 'debut sm' 
      do ipp=1,nppmax
         cellpop(ipp) = -999
         listcell(ipp) = -999
      enddo
*     naissance ou mort ?
      b = ggrbinom(1.d0,0.5d0)
      if(b .eq. 1) then
         if(npop .lt. npopmax) then 
c            write(*,*) 'naissance'
*     split
*     choix de la pop qui split
            isplit = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
*     recherche des cellules affectees a cette pop
            call who(c,isplit,npp,nppmax,cellpop,ncellpop)
            if(ncellpop .gt. 0) then
*     tirage du nombre de cellules reallouees
               nu = idint(dint(dble(ncellpop+1)*ggrunif(0.d0,1.d0)))
               if(nu .gt. 0) then
*     tirage des cellules reallouees
                  call sample2(cellpop,nppmax,nu,ncellpop,listcell)
*     proposition de reallocation dans la pop npop+1
                  call split(npop+1,c,ctmp,nppmax,nu,listcell)
               else 
                  do ipp = 1,nppmax
                     ctmp(ipp) = c(ipp)
                  enddo
               endif
            else
               nu = 0
               do ipp = 1,nppmax
                  ctmp(ipp) = c(ipp)
               enddo
            endif
*     comptage des alleles sur chaque locus pour c puis ctmp
            call countn(nindiv,nloc,nloc2,npopmax,
     &           nppmax,nal,nalmax,zz,n,indcell,c,ploidy)
            call countn(nindiv,nloc,nloc2,npopmax,
     &           nppmax,nal,nalmax,zz,ntmp,indcell,ctmp,ploidy)
*     nouvelles derives      
            do ipop = 1,npop
               drifttmp(ipop) = drift(ipop)
            enddo
c     next line modified 29/02/08
c            rr = ggrunif(0.d0,1.d0)*deltad
            rr = ggrnorm(0.d0,1.d0)*deltad
            drifttmp(isplit) = drift(isplit) - rr
            drifttmp(npop+1) = drift(isplit) + rr
            if(((drifttmp(isplit)   .gt. 1d-300) .and. 
     &          (1-drifttmp(isplit) .gt. 1d-300)).and.
     &         ((drifttmp(npop+1)   .gt. 1d-300) .and. 
     &          (1-drifttmp(npop+1) .gt. 1d-300))) then

*     terme des proposal sur c
               lratio = dlog(2*dble(ncellpop+1)) + 
     &              lbico(ncellpop,nu) - dlog(dble(npop+1)) 
*     terme des priors sur c
               lratio = lratio + dble(npp)*(dlog(dble(npop)) - 
     &              dlog(dble(npop+1)))
*     terme des frequences
               lratio = lratio 
     &         + lTf(isplit,ntmp,fa,drifttmp,npopmax,nloc,nal,nalmax)  
     &         + lTf(npop+1,ntmp,fa,drifttmp,npopmax,nloc,nal,nalmax) 
     &         - lTf(isplit,n,fa,drift,npopmax,nloc,nal,nalmax) 
*     term proposal drift
               lratio = lratio + dlog(2*deltad) 
     &              + .5*rr**2 + dlog(dsqrt(2*3.141593d0))
*     term prior drift
               lratio = lratio 
     &              + (shape1-1)*dlog(drifttmp(isplit))
     &              + (shape2-1)*dlog(1-drifttmp(isplit))
     &              + (shape1-1)*dlog(drifttmp(npop+1))
     &              + (shape2-1)*dlog(1-drifttmp(npop+1))
     &              - (shape1-1)*dlog(drift(isplit))
     &              - (shape2-1)*dlog(1-drift(isplit))
     &              + gglgamfn(shape1+shape2)
     &              - gglgamfn(shape1)-gglgamfn(shape2)
c            write(*,*) 'lratio=',lratio
            
               lratio = dmin1(0.d0,lratio)
               alpha = dexp(lratio)
               bern = ggrbinom(1.d0,alpha)
             
               if(bern .eq. 1) then
*     proposition nouvelle freq 
                  call addfreq8(npop,npopmax,nloc,nloc,
     &                 nal,nalmax,isplit,
     &                 f,ftmp,fa,drift,drifttmp,a,ptmp,ntmp)
                  call accept5(nppmax,npopmax,nloc,
     &                 nalmax,nal,c,ctmp,f,ftmp,drift,drifttmp)
                  npop = npop + 1
               endif
            endif
         endif
         
*     merge
      else
         if(npop .gt. npopmin) then 
c      write(*,*) 'mort'
*     tirage de la pop qui meurt
            ipoprem = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
*     tirage de la pop hote
            ipophost = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
            do while(ipophost .eq. ipoprem)
               ipophost = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
            enddo
            
c            if(abs(drift(ipoprem)-drift(ipophost)) .lt. 2*deltad) then
*     on range dans la pop d'indice le plus petit
               if(ipophost .gt. ipoprem) then
                  ii = ipophost
                  ipophost = ipoprem
                  ipoprem = ii
               endif
*     recherche des cellules qui vont etre reallouees
               call who(c,ipoprem,npp,nppmax,cellpop,ncellpop)
*     recherche des cellules de la pop hote
               call who(c,ipophost,npp,nppmax,cellpophost,
     &              ncellpophost)
*     proposition de reallocation dans la pop ipophost
               call merging(ipoprem,ipophost,c,ctmp,nppmax,
     &              ncellpop,cellpop)
*     comptage des alleles sur chaque locus pour c puis ctmp
               call countn(nindiv,nloc,nloc2,npopmax,
     &              nppmax,nal,nalmax,zz,n,indcell,c,ploidy)
               call countn(nindiv,nloc,nloc2,npopmax,
     &              nppmax,nal,nalmax,zz,ntmp,indcell,ctmp,ploidy)
*     enleve une pop dans le tableau tmporaire des derives            
               call remdrift2(ipoprem,ipophost,npop,npopmax,drift,
     &           drifttmp)
               rr = (drift(ipophost) - drift(ipoprem))/(2*deltad)
c$$$               write(*,*) 'drifttmp=',drifttmp
*     calcul du log du ratio  
*     terme des proposal sur c
               lratio = dlog(dble(npop)) - 
     &              dlog(2*dble(ncellpop+ncellpophost+1)) -
     &              lbico(ncellpop+ncellpophost,ncellpop) 
*     terme des priors sur c
               lratio = lratio + 
     &              dble(npp)*(dlog(dble(npop)) - 
     &              dlog(dble(npop-1)))
c$$$               write(*,*) 'term en c lratio =',lratio
*     term des freq
               lratio = lratio + 
     &        lTf(ipophost,ntmp,fa,drifttmp,npopmax,nloc,nal,nalmax) -  
     &        lTf(ipophost,n,fa,drift,npopmax,nloc,nal,nalmax) - 
     &        lTf(ipoprem,n,fa,drift,npopmax,nloc,nal,nalmax) 
*     terme proposal d
               lratio = lratio - dlog(2*deltad)
     &              - .5*rr**2 - dlog(dsqrt(2*3.141593d0))
*     term prior drift
               lratio = lratio 
     &              + (shape1-1)*dlog(drifttmp(ipophost)) 
     &              + (shape2-1)*dlog(1-drifttmp(ipophost))
     &              - (shape1-1)*dlog(drift(ipophost))
     &              - (shape2-1)*dlog(1-drift(ipophost))
     &              - (shape1-1)*dlog(drift(ipoprem))
     &              - (shape2-1)*dlog(1-drift(ipoprem))
     &              - gglgamfn(shape1+shape2)
     &              + gglgamfn(shape1)+gglgamfn(shape2)               

               lratio = dmin1(0.d0,lratio)
               alpha  = dexp(lratio)
               bern = ggrbinom(1.d0,alpha)      
               
               if(bern .eq. 1) then
*     propostion du nouveau tableau de freq et de derives
                  call remfreq8(ipoprem,ipophost,
     &                 npop,npopmax,nloc,nloc,nal,
     &                 nalmax,f,ftmp,drift,drifttmp,fa,a,ptmp,ntmp)
                  call accept5(nppmax,npopmax,nloc,
     &                 nalmax,nal,c,ctmp,f,ftmp,drift,drifttmp)
                  npop = npop - 1
               endif
            endif
         endif
       end subroutine sm2
***********************************************************************



************************************************************************
*     split and merge of populations with acceptance according to MH ratio
*     CFM
*     (d1*,d2*) = d1-u,d2+u
*     MH ratio does not depend on proposed frequencies
*     extended for quantitative variables 
      subroutine smfallvar3(npop,npopmin,npopmax,f,fa,drift,
     &     nlocd,nloch,ncolt,nql,
     &     nal,nalmax,indcell,nindiv,npp,nppmax,c,ctmp,
     &     a,ptmp,ftmp,drifttmp,yy,z,ql,cellpop,listcell,
     &     cellpophost,n,ntmp,qtc,nqtc,meanqtc,sdqtc,
     &     meanqtctmp,sdqtctmp,nnqtc,sqtc,ssqtc,
     &     ksiqtc,kappaqtc,alphaqtc,betaqtc,shape1,shape2,
     &     usegeno2,usegeno1,useql,useqtc,ploidy)
      implicit none 
      integer npop,npopmin,npopmax,ncolt,nal(ncolt),
     &     nalmax,nindiv,npp,nppmax,indcell(nindiv),nlocd,nlocd2,
     &     nloch,nql,c(nppmax),ctmp(nppmax),yy(nindiv,2*nlocd+2*nloch),
     &     z(nindiv,nloch),ql(nindiv,nql),
     &     n(npopmax,ncolt,nalmax),ntmp(npopmax,ncolt,nalmax),
     &     nqtc,nnqtc(npopmax,nqtc),usegeno2,usegeno1,useql,useqtc,
     &     ploidy
      double precision f(npopmax,ncolt,nalmax),drift(npopmax),
     &      ftmp(npopmax,ncolt,nalmax),drifttmp(npopmax),
     &     a(nalmax),ptmp(nalmax),fa(ncolt,nalmax),qtc(nindiv,nqtc),
     &     meanqtc(npopmax,nqtc),sdqtc(npopmax,nqtc),
     &     meanqtctmp(npopmax,nqtc),sdqtctmp(npopmax,nqtc),
     &     sqtc(npopmax,nqtc),ssqtc(npopmax,nqtc),
     &     ksiqtc,kappaqtc,alphaqtc,betaqtc,shape1,shape2
      integer ipoprem,ipp,isplit,
     &     cellpop(nppmax),ncellpop,nu,listcell(nppmax),
     &     ipophost,ncellpophost,cellpophost(nppmax),ii
      double precision alpha,ggrunif,lbico,lratio,lTfallvar
      integer iloc,iqtc,ipop,ipoptmp,iindiv
      double precision llr6,termf9bis,gglgamfn,bern,b,ggrbinom,ggdnorm,
     &     ggdgamma,ggrnorm,lrallvar2,contriblr,deltad,rr,junk
      dimension ksiqtc(nqtc),kappaqtc(nqtc),alphaqtc(nqtc),betaqtc(nqtc)
 
      deltad = shape1/(shape1+shape2)

c         write(*,*) 'debut sm' 
      do ipp=1,nppmax
         cellpop(ipp) = -999
         listcell(ipp) = -999
      enddo
*     naissance ou mort ?
      b = ggrbinom(1.d0,0.5d0)
      if(b .eq. 1) then
         if(npop .lt. npopmax) then 
c            write(*,*) 'naissance'
c            write(*,*) 'npop=',npop    
*     split
*     choix de la pop qui split
            isplit = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
*     recherche des cellules affectees a cette pop
            call who(c,isplit,npp,nppmax,cellpop,ncellpop)
            if(ncellpop .gt. 0) then
*     tirage du nombre de cellules reallouees
               nu = idint(dint(dble(ncellpop+1)*ggrunif(0.d0,1.d0)))
               if(nu .gt. 0) then
*     tirage des cellules reallouees
                  call sample2(cellpop,nppmax,nu,ncellpop,listcell)
*     proposition de reallocation dans la pop npop+1
                  call split(npop+1,c,ctmp,nppmax,nu,listcell)
               else 
                  do ipp = 1,nppmax
                     ctmp(ipp) = c(ipp)
                  enddo
               endif
            else
               nu = 0
               do ipp = 1,nppmax
                  ctmp(ipp) = c(ipp)
               enddo
            endif
c            write(*,*) 'c=',c(1),c(2),c(3),c(4)
c            write(*,*) 'ctmp=',ctmp(1),ctmp(2),ctmp(3),ctmp(4)

            if(((usegeno2 .eq. 1) .or. (usegeno1 .eq. 1)) .or. 
     &           useql .eq. 1) then
*     comptage des alleles sur chaque locus pour c puis ctmp
c               write(*,*) 'dans smfallvar usegeno2=',usegeno2
               call countnallvar2(nindiv,nlocd,nloch,nql,ncolt,
     &              nppmax,nal,nalmax,yy,z,ql,n,indcell,c,npopmax,
     &              usegeno2,usegeno1,useql,ploidy)
               call countnallvar2(nindiv,nlocd,nloch,nql,ncolt,
     &              nppmax,nal,nalmax,yy,z,ql,ntmp,indcell,ctmp,npopmax,
     &              usegeno2,usegeno1,useql,ploidy)
c               write(*,*) 'n(',isplit,',1,1)=',n(isplit,1,1)
c               write(*,*) 'ntmp(',isplit,',1,1)=',ntmp(isplit,1,1)
c               write(*,*) 'ntmp(',npop+1,',1,1)=',ntmp(npop+1,1,1)
*     nouvelles derives      
               do ipop = 1,npop
                  drifttmp(ipop) = drift(ipop)
               enddo
               rr = ggrnorm(0.d0,1.d0)*deltad
               drifttmp(isplit) = drift(isplit) - rr
               drifttmp(npop+1) = drift(isplit) + rr
            endif
c            write(*,*) 'drift = ',drift
c            write(*,*) 'drifttmp = ',drifttmp
            if(((drifttmp(isplit)   .gt. 1d-300) .and. 
     &           (1-drifttmp(isplit) .gt. 1d-300)).and.
     &           ((drifttmp(npop+1)   .gt. 1d-300) .and. 
     &           (1-drifttmp(npop+1) .gt. 1d-300))) then
               
*     terme des proposal sur c
               lratio = dlog(2*dble(ncellpop+1)) + 
     &              lbico(ncellpop,nu) - dlog(dble(npop+1)) 
c               write(6,*) 'split apres proposal c lratio=',lratio 
*     terme des priors sur c
               lratio = lratio + dble(npp)*(dlog(dble(npop)) - 
     &              dlog(dble(npop+1)))
c               write(6,*) 'split apres prior c lratio=',lratio 
c               write(*,*) 'term en c lratio =',lratio               
               if(((usegeno2 .eq. 1) .or. (usegeno1 .eq. 1)) .or. 
     &              useql .eq. 1) then
*     contribition of frequencies, this term includes likelihood ratio
                  lratio = lratio + 
     &                 lTfallvar(isplit,ntmp,fa,drifttmp,npopmax,nlocd,
     &                 nloch,nql,ncolt,nal,nalmax,usegeno2,usegeno1,
     &                 useql) + 
     &                 lTfallvar(npop+1,ntmp,fa,drifttmp,npopmax,nlocd,
     &                 nloch,nql,ncolt,nal,nalmax,usegeno2,usegeno1,
     &                 useql) - 
     &                 lTfallvar(isplit,n,fa,drift,npopmax,nlocd,
     &                 nloch,nql,ncolt,nal,nalmax,usegeno2,usegeno1,
     &                 useql)
c$$$         write(*,*) lTfallvar(isplit,ntmp,fa,drifttmp,npopmax,nlocd,
c$$$     &                 nloch,nql,ncolt,nal,nalmax,usegeno2,usegeno1,
c$$$     &                 useql)
c$$$         write(*,*) lTfallvar(npop+1,ntmp,fa,drifttmp,npopmax,nlocd,
c$$$     &                 nloch,nql,ncolt,nal,nalmax,usegeno2,usegeno1,
c$$$     &                 useql)
c$$$         write(*,*) -lTfallvar(isplit,n,fa,drift,npopmax,nlocd,
c$$$     &                 nloch,nql,ncolt,nal,nalmax,usegeno2,usegeno1,
c$$$     &                 useql)
c$$$         write(*,*) 'apres Tf lratio =',lratio     
c                  write(6,*) 'split apres prior contrib freq =',lratio 
*     term proposal drift
                  lratio = lratio + dlog(2*deltad) 
     &                 + .5*rr**2 + dlog(dsqrt(2*3.141593d0))
*     term prior drift
                  lratio = lratio 
     &                 + (shape1-1)*dlog(drifttmp(isplit))
     &                 + (shape2-1)*dlog(1-drifttmp(isplit))
     &                 + (shape1-1)*dlog(drifttmp(npop+1))
     &                 + (shape2-1)*dlog(1-drifttmp(npop+1))
     &                 - (shape1-1)*dlog(drift(isplit))
     &                 - (shape2-1)*dlog(1-drift(isplit))
     &                 + gglgamfn(shape1+shape2)
     &                 - gglgamfn(shape1)-gglgamfn(shape2)
c                  write(*,*) 'apres T drift lratio=',lratio
c                  write(6,*) 'split apres prior contrib drift =',lratio 
               endif

               if(useqtc .eq. 1) then 
*     propose mean and variance of quantitative variables
                  call propqvsplit3(qtc,nindiv,nqtc,indcell,c,ctmp,
     &                 meanqtc,sdqtc,meanqtctmp,sdqtctmp,nnqtc,sqtc,
     &                 ssqtc,npop,npopmax,nppmax,ksiqtc,kappaqtc,
     &                 alphaqtc,betaqtc,isplit,contriblr)
*     contrib prior and proposal 
                  lratio = lratio + contriblr
*     contrib likelihood quant. variable
*     argument usegeno2 is passed as 0 to ignore genotypes 
*     hence avoid having likelihood ratio of genotypes twice 
                  lratio = lratio + lrallvar2(yy,z,ql,
     &                 f,c,ctmp,indcell,indcell,npopmax,nlocd,
     &                 nloch,nql,ncolt,nalmax,nindiv,nppmax,qtc,nqtc,
     &                 meanqtc,sdqtc,meanqtctmp,sdqtctmp,0,0,0,useqtc,
     &                 ploidy)
               endif
               lratio = dmin1(0.d0,lratio)
               alpha = dexp(lratio)
               bern = ggrbinom(1.d0,alpha)
             
               if(bern .eq. 1) then
                  if(((usegeno2 .eq. 1) .or. (usegeno1 .eq. 1)) .or. 
     &                 useql .eq. 1) then
*     proposition nouvelle freq 
                     call addfall(npop,npopmax,nlocd,nlocd2,nloch,nql,
     &                    ncolt,nal,nalmax,isplit,f,ftmp,fa,drift,
     &                    drifttmp,a,ptmp,ntmp, 
     &                    usegeno2,usegeno1,useql)
                  endif
                  if(useqtc .eq. 1) then 
                     do iqtc = 1,nqtc
                        do ipop = 1,npopmax
                           meanqtc(ipop,iqtc) =  meanqtctmp(ipop,iqtc)
                           sdqtc(ipop,iqtc) =  sdqtctmp(ipop,iqtc)
                        enddo
                     enddo
                  endif
                  call accept5(nppmax,npopmax,ncolt,
     &                 nalmax,nal,c,ctmp,f,ftmp,drift,drifttmp)
                  npop = npop + 1
               endif
            endif
         endif
*     merge
      else
         if(npop .gt. npopmin) then 
c      write(*,*) 'mort'
c      write(*,*) 'npop=',npop
*     tirage de la pop qui meurt
            ipoprem = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
*     tirage de la pop hote
            ipophost = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
            do while(ipophost .eq. ipoprem)
               ipophost = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
            enddo
*     on range dans la pop d'indice le plus petit
            if(ipophost .gt. ipoprem) then
               ii = ipophost
               ipophost = ipoprem
               ipoprem = ii
            endif
*     recherche des cellules qui vont etre reallouees
            call who(c,ipoprem,npp,nppmax,cellpop,ncellpop)
*     recherche des cellules de la pop hote
            call who(c,ipophost,npp,nppmax,cellpophost,
     &           ncellpophost)
*     proposition de reallocation dans la pop ipophost
            call merging(ipoprem,ipophost,c,ctmp,nppmax,
     &           ncellpop,cellpop)
c            write(*,*) 'c=',c(1),c(2),c(3),c(4)
c            write(*,*) 'ctmp=',ctmp(1),ctmp(2),ctmp(3),ctmp(4) 

*     comptage des alleles sur chaque locus pour c puis ctmp
            if(((usegeno2 .eq. 1) .or. (usegeno1 .eq. 1)) .or. 
     &           useql .eq. 1) then
c               write(*,*) 'dans smfallvar usegeno2=',usegeno2
               call countnallvar2(nindiv,nlocd,nloch,nql,ncolt,
     &              nppmax,nal,nalmax,yy,z,ql,n,indcell,c,npopmax,
     &              usegeno2, usegeno1,useql,ploidy)
               call countnallvar2(nindiv,nlocd,nloch,nql,ncolt,
     &              nppmax,nal,nalmax,yy,z,ql,ntmp,indcell,ctmp,npopmax,
     &              usegeno2, usegeno1,useql,ploidy)
c               write(*,*) 'n(',ipoprem,',1,1)=',n(ipoprem,1,1)
c               write(*,*) 'n(',ipophost,',1,1)=',n(ipophost,1,1)
c               write(*,*) 'ntmp(',ipophost,',1,1)=',ntmp(ipophost,1,1)

*     enleve une pop dans le tableau tmporaire des derives            
               call remdrift2(ipoprem,ipophost,npop,npopmax,drift,
     &              drifttmp)
               rr = (drift(ipophost) - drift(ipoprem))/(2*deltad)
c               write(*,*) 'drift=',drift
c               write(*,*) 'drifttmp=',drifttmp
            endif
*     calcul du log du ratio  
*     terme des proposal sur c
            lratio = dlog(dble(npop)) - 
     &           dlog(2*dble(ncellpop+ncellpophost+1)) -
     &           lbico(ncellpop+ncellpophost,ncellpop) 
c            write(6,*) 'merge apres proposal c lratio=',lratio 
*     terme des priors sur c
            lratio = lratio + 
     &           dble(npp)*(dlog(dble(npop)) - 
     &           dlog(dble(npop-1)))
c            write(6,*) 'merge apres prior c lratio=',lratio 
c           write(*,*) 'term en c lratio =',lratio
            if(((usegeno2 .eq. 1) .or. (usegeno1 .eq. 1)) .or. 
     &           useql .eq. 1) then
*     term des freq
c$$$  *     term des freq
               lratio = lratio + 
     &              lTfallvar(ipophost,ntmp,fa,drifttmp,npopmax,nlocd,
     &              nloch,nql,ncolt,nal,nalmax,usegeno2,usegeno1,
     &              useql) -
     &              lTfallvar(ipophost,n,fa,drift,npopmax,nlocd,
     &              nloch,nql,ncolt,nal,nalmax,usegeno2,usegeno1,
     &              useql) - 
     &              lTfallvar(ipoprem,n,fa,drift,npopmax,nlocd,
     &              nloch,nql,ncolt,nal,nalmax,usegeno2,usegeno1,
     &              useql)
c               write(6,*) 'merge apres prior contrib freq =',lratio 
c$$$           write(*,*) lTfallvar(ipophost,ntmp,fa,drifttmp,npopmax,nlocd,
c$$$     &              nloch,nql,ncolt,nal,nalmax,usegeno2,usegeno1,
c$$$     &              useql)
c$$$           write(*,*) -lTfallvar(ipophost,n,fa,drift,npopmax,nlocd,
c$$$     &              nloch,nql,ncolt,nal,nalmax,usegeno2,usegeno1,
c$$$     &              useql)
c$$$           write(*,*) -lTfallvar(ipoprem,n,fa,drift,npopmax,nlocd,
c$$$     &              nloch,nql,ncolt,nal,nalmax,usegeno2,usegeno1,
c$$$     &              useql)

c          write(*,*) 'apres Tf lratio =',lratio 
*     terme proposal d
               lratio = lratio - dlog(2*deltad)
     &              - .5*rr**2 - dlog(dsqrt(2*3.141593d0))
*     term prior drift
               lratio = lratio 
     &              + (shape1-1)*dlog(drifttmp(ipophost)) 
     &              + (shape2-1)*dlog(1-drifttmp(ipophost))
     &              - (shape1-1)*dlog(drift(ipophost))
     &              - (shape2-1)*dlog(1-drift(ipophost))
     &              - (shape1-1)*dlog(drift(ipoprem))
     &              - (shape2-1)*dlog(1-drift(ipoprem))
     &              - gglgamfn(shape1+shape2)
     &              + gglgamfn(shape1)+gglgamfn(shape2)  
c               write(*,*) 'apres T drift lratio=',lratio
c               write(6,*) 'merge apres prior contrib drift =',lratio 
            endif

            if(useqtc .eq. 1) then 
*     propose mean and variance of quantitative variables
               call propqvmrg3(qtc,nindiv,nqtc,indcell,c,ctmp,
     &              meanqtc,sdqtc,meanqtctmp,sdqtctmp,nnqtc,sqtc,
     &              ssqtc,npop,npopmax,nppmax,
     &              ksiqtc,kappaqtc,alphaqtc,betaqtc,ipophost,
     &              ipoprem,contriblr)
*     contrib prior and proposal 
               lratio = lratio + contriblr   
*     contrib likelihood quant. variable
*     argument usegeno2 is passed as 0 to ignore genotypes 
*     hence avoid having likelihood ratio of genotypes twice 
               lratio = lratio + lrallvar2(yy,z,ql,
     &              f,c,ctmp,indcell,indcell,npopmax,nlocd,
     &              nloch,nql,ncolt,nalmax,nindiv,nppmax,qtc,nqtc,
     &              meanqtc,sdqtc,meanqtctmp,sdqtctmp,0,0,0,useqtc,
     &              ploidy)  
            endif
c            write(*,*) 'lratio for merge=',lratio
            lratio = dmin1(0.d0,lratio)
            alpha  = dexp(lratio)
            bern = ggrbinom(1.d0,alpha)      
            
            if(bern .eq. 1) then
               if(((usegeno2 .eq. 1) .or. (usegeno1 .eq. 1)) .or. 
     &              useql .eq. 1) then
*     propostion du nouveau tableau de freq et de derives
                  call  remfall(ipoprem,ipophost,
     &                    npop,npopmax,nlocd,nloch,nql,ncolt,nal,
     &                 nalmax,f,ftmp,drift,drifttmp,fa,a,ptmp,ntmp,
     &                 usegeno2,usegeno1,useql)
               endif
               if(useqtc .eq. 1) then 
                  do iqtc = 1,nqtc
                     do ipop = 1,npopmax
                        meanqtc(ipop,iqtc) =  meanqtctmp(ipop,iqtc)
                        sdqtc(ipop,iqtc) =  sdqtctmp(ipop,iqtc)
                     enddo
                  enddo
               endif
               call accept5(nppmax,npopmax,ncolt,
     &              nalmax,nal,c,ctmp,f,ftmp,drift,drifttmp)
               npop = npop - 1
            endif
         endif
      endif
      end subroutine smfallvar3
***********************************************************************     



c$$$
c$$$************************************************************************
c$$$*     split and merge of populations with acceptance according to MH ratio
c$$$*     CFM
c$$$*     (d1*,d2*) = d1-u,d2+u
c$$$*     MH ratio does not depend on proposed frequencies
c$$$*     extended for quantitative variables 
c$$$      subroutine smfallvar2(npop,npopmin,npopmax,f,fa,drift,
c$$$     &     nlocd,nloch,ncolt,nql,
c$$$     &     nal,nalmax,indcell,nindiv,npp,nppmax,c,ctmp,
c$$$     &     a,ptmp,ftmp,drifttmp,yy,z,ql,cellpop,listcell,
c$$$     &     cellpophost,n,ntmp,qtc,nqtc,meanqtc,sdqtc,
c$$$     &     meanqtctmp,sdqtctmp,nnqtc,sqtc,ssqtc,
c$$$     &     ksiqtc,kappaqtc,alphaqtc,betaqtc,shape1,shape2,
c$$$     &     usegeno2,usegeno1,useql,useqtc,ploidy)
c$$$      implicit none 
c$$$      integer npop,npopmin,npopmax,ncolt,nal(ncolt),
c$$$     &     nalmax,nindiv,npp,nppmax,indcell(nindiv),nlocd,nlocd2,
c$$$     &     nloch,nql,c(nppmax),ctmp(nppmax),yy(nindiv,2*nlocd+2*nloch),
c$$$     &     z(nindiv,nloch),ql(nindiv,nql),
c$$$     &     n(npopmax,ncolt,nalmax),ntmp(npopmax,ncolt,nalmax),
c$$$     &     nqtc,nnqtc,usegeno2,usegeno1,useql,useqtc,ploidy
c$$$      double precision f(npopmax,ncolt,nalmax),drift(npopmax),
c$$$     &      ftmp(npopmax,ncolt,nalmax),drifttmp(npopmax),
c$$$     &     a(nalmax),ptmp(nalmax),fa(ncolt,nalmax),qtc(nindiv,nqtc),
c$$$     &     meanqtc(npopmax,nqtc),sdqtc(npopmax,nqtc),
c$$$     &     meanqtctmp(npopmax,nqtc),sdqtctmp(npopmax,nqtc),
c$$$     &     sqtc(npopmax,nqtc),ssqtc(npopmax,nqtc),
c$$$     &     ksiqtc,kappaqtc,alphaqtc,betaqtc,shape1,shape2
c$$$      integer ipoprem,ipp,isplit,
c$$$     &     cellpop(nppmax),ncellpop,nu,listcell(nppmax),
c$$$     &     ipophost,ncellpophost,cellpophost(nppmax),ii
c$$$      double precision alpha,ggrunif,lbico,lratio,lTfallvar
c$$$      integer iloc,iqtc,ipop,ipoptmp,iindiv
c$$$      double precision llr6,termf9bis,gglgamfn,bern,b,ggrbinom,ggdnorm,
c$$$     &     ggdgamma,ggrnorm,lrallvar2,contriblr,deltad,rr,junk
c$$$      dimension ksiqtc(nqtc),kappaqtc(nqtc),alphaqtc(nqtc),betaqtc(nqtc)
c$$$ 
c$$$      deltad = shape1/(shape1+shape2)
c$$$
c$$$c         write(*,*) 'debut sm' 
c$$$      do ipp=1,nppmax
c$$$         cellpop(ipp) = -999
c$$$         listcell(ipp) = -999
c$$$      enddo
c$$$*     naissance ou mort ?
c$$$      b = ggrbinom(1.d0,0.5d0)
c$$$      if(b .eq. 1) then
c$$$         if(npop .lt. npopmax) then 
c$$$c            write(*,*) 'naissance'
c$$$c            write(*,*) 'npop=',npop    
c$$$*     split
c$$$*     choix de la pop qui split
c$$$            isplit = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
c$$$*     recherche des cellules affectees a cette pop
c$$$            call who(c,isplit,npp,nppmax,cellpop,ncellpop)
c$$$            if(ncellpop .gt. 0) then
c$$$*     tirage du nombre de cellules reallouees
c$$$               nu = idint(dint(dble(ncellpop+1)*ggrunif(0.d0,1.d0)))
c$$$               if(nu .gt. 0) then
c$$$*     tirage des cellules reallouees
c$$$                  call sample2(cellpop,nppmax,nu,ncellpop,listcell)
c$$$*     proposition de reallocation dans la pop npop+1
c$$$                  call split(npop+1,c,ctmp,nppmax,nu,listcell)
c$$$               else 
c$$$                  do ipp = 1,nppmax
c$$$                     ctmp(ipp) = c(ipp)
c$$$                  enddo
c$$$               endif
c$$$            else
c$$$               nu = 0
c$$$               do ipp = 1,nppmax
c$$$                  ctmp(ipp) = c(ipp)
c$$$               enddo
c$$$            endif
c$$$c            write(*,*) 'c=',c(1),c(2),c(3),c(4)
c$$$c            write(*,*) 'ctmp=',ctmp(1),ctmp(2),ctmp(3),ctmp(4)
c$$$
c$$$            if(((usegeno2 .eq. 1) .or. (usegeno1 .eq. 1)) .or. 
c$$$     &           useql .eq. 1) then
c$$$*     comptage des alleles sur chaque locus pour c puis ctmp
c$$$c               write(*,*) 'dans smfallvar usegeno2=',usegeno2
c$$$               call countnallvar2(nindiv,nlocd,nloch,nql,ncolt,
c$$$     &              nppmax,nal,nalmax,yy,z,ql,n,indcell,c,npopmax,
c$$$     &              usegeno2,usegeno1,useql,ploidy)
c$$$               call countnallvar2(nindiv,nlocd,nloch,nql,ncolt,
c$$$     &              nppmax,nal,nalmax,yy,z,ql,ntmp,indcell,ctmp,npopmax,
c$$$     &              usegeno2,usegeno1,useql,ploidy)
c$$$c               write(*,*) 'n(',isplit,',1,1)=',n(isplit,1,1)
c$$$c               write(*,*) 'ntmp(',isplit,',1,1)=',ntmp(isplit,1,1)
c$$$c               write(*,*) 'ntmp(',npop+1,',1,1)=',ntmp(npop+1,1,1)
c$$$*     nouvelles derives      
c$$$               do ipop = 1,npop
c$$$                  drifttmp(ipop) = drift(ipop)
c$$$               enddo
c$$$               rr = ggrnorm(0.d0,1.d0)*deltad
c$$$               drifttmp(isplit) = drift(isplit) - rr
c$$$               drifttmp(npop+1) = drift(isplit) + rr
c$$$            endif
c$$$c            write(*,*) 'drift = ',drift
c$$$c            write(*,*) 'drifttmp = ',drifttmp
c$$$            if(((drifttmp(isplit)   .gt. 1d-300) .and. 
c$$$     &           (1-drifttmp(isplit) .gt. 1d-300)).and.
c$$$     &           ((drifttmp(npop+1)   .gt. 1d-300) .and. 
c$$$     &           (1-drifttmp(npop+1) .gt. 1d-300))) then
c$$$               
c$$$*     terme des proposal sur c
c$$$               lratio = dlog(2*dble(ncellpop+1)) + 
c$$$     &              lbico(ncellpop,nu) - dlog(dble(npop+1)) 
c$$$c               write(6,*) 'split apres proposal c lratio=',lratio 
c$$$*     terme des priors sur c
c$$$               lratio = lratio + dble(npp)*(dlog(dble(npop)) - 
c$$$     &              dlog(dble(npop+1)))
c$$$c               write(6,*) 'split apres prior c lratio=',lratio 
c$$$c               write(*,*) 'term en c lratio =',lratio               
c$$$               if(((usegeno2 .eq. 1) .or. (usegeno1 .eq. 1)) .or. 
c$$$     &              useql .eq. 1) then
c$$$*     contribition of frequencies, this term includes likelihood ratio
c$$$                  lratio = lratio + 
c$$$     &                 lTfallvar(isplit,ntmp,fa,drifttmp,npopmax,nlocd,
c$$$     &                 nloch,nql,ncolt,nal,nalmax,usegeno2,usegeno1,
c$$$     &                 useql) + 
c$$$     &                 lTfallvar(npop+1,ntmp,fa,drifttmp,npopmax,nlocd,
c$$$     &                 nloch,nql,ncolt,nal,nalmax,usegeno2,usegeno1,
c$$$     &                 useql) - 
c$$$     &                 lTfallvar(isplit,n,fa,drift,npopmax,nlocd,
c$$$     &                 nloch,nql,ncolt,nal,nalmax,usegeno2,usegeno1,
c$$$     &                 useql)
c$$$c$$$         write(*,*) lTfallvar(isplit,ntmp,fa,drifttmp,npopmax,nlocd,
c$$$c$$$     &                 nloch,nql,ncolt,nal,nalmax,usegeno2,usegeno1,
c$$$c$$$     &                 useql)
c$$$c$$$         write(*,*) lTfallvar(npop+1,ntmp,fa,drifttmp,npopmax,nlocd,
c$$$c$$$     &                 nloch,nql,ncolt,nal,nalmax,usegeno2,usegeno1,
c$$$c$$$     &                 useql)
c$$$c$$$         write(*,*) -lTfallvar(isplit,n,fa,drift,npopmax,nlocd,
c$$$c$$$     &                 nloch,nql,ncolt,nal,nalmax,usegeno2,usegeno1,
c$$$c$$$     &                 useql)
c$$$c$$$         write(*,*) 'apres Tf lratio =',lratio     
c$$$c                  write(6,*) 'split apres prior contrib freq =',lratio 
c$$$*     term proposal drift
c$$$                  lratio = lratio + dlog(2*deltad) 
c$$$     &                 + .5*rr**2 + dlog(dsqrt(2*3.141593d0))
c$$$*     term prior drift
c$$$                  lratio = lratio 
c$$$     &                 + (shape1-1)*dlog(drifttmp(isplit))
c$$$     &                 + (shape2-1)*dlog(1-drifttmp(isplit))
c$$$     &                 + (shape1-1)*dlog(drifttmp(npop+1))
c$$$     &                 + (shape2-1)*dlog(1-drifttmp(npop+1))
c$$$     &                 - (shape1-1)*dlog(drift(isplit))
c$$$     &                 - (shape2-1)*dlog(1-drift(isplit))
c$$$     &                 + gglgamfn(shape1+shape2)
c$$$     &                 - gglgamfn(shape1)-gglgamfn(shape2)
c$$$c                  write(*,*) 'apres T drift lratio=',lratio
c$$$c                  write(6,*) 'split apres prior contrib drift =',lratio 
c$$$               endif
c$$$
c$$$               if(useqtc .eq. 1) then 
c$$$*     propose mean and variance of quantitative variables
c$$$                  call propparqvsplit(qtc,nindiv,nqtc,indcell,c,ctmp,
c$$$     &                 meanqtc,sdqtc,meanqtctmp,sdqtctmp,nnqtc,sqtc,
c$$$     &                 ssqtc,npop,npopmax,nppmax,ksiqtc,kappaqtc,
c$$$     &                 alphaqtc,betaqtc,isplit,contriblr)
c$$$*     contrib prior and proposal 
c$$$                  lratio = lratio + contriblr
c$$$*     contrib likelihood quant. variable
c$$$*     argument usegeno2 is passed as 0 to ignore genotypes 
c$$$*     hence avoid having likelihood ratio of genotypes twice 
c$$$                  lratio = lratio + lrallvar2(yy,z,ql,
c$$$     &                 f,c,ctmp,indcell,indcell,npopmax,nlocd,
c$$$     &                 nloch,nql,ncolt,nalmax,nindiv,nppmax,qtc,nqtc,
c$$$     &                 meanqtc,sdqtc,meanqtctmp,sdqtctmp,0,0,0,useqtc,
c$$$     &                 ploidy)
c$$$               endif
c$$$               lratio = dmin1(0.d0,lratio)
c$$$               alpha = dexp(lratio)
c$$$               bern = ggrbinom(1.d0,alpha)
c$$$             
c$$$               if(bern .eq. 1) then
c$$$                  if(((usegeno2 .eq. 1) .or. (usegeno1 .eq. 1)) .or. 
c$$$     &                 useql .eq. 1) then
c$$$*     proposition nouvelle freq 
c$$$                     call addfall(npop,npopmax,nlocd,nlocd2,nloch,nql,
c$$$     &                    ncolt,nal,nalmax,isplit,f,ftmp,fa,drift,
c$$$     &                    drifttmp,a,ptmp,ntmp, 
c$$$     &                    usegeno2,usegeno1,useql)
c$$$                  endif
c$$$                  if(useqtc .eq. 1) then 
c$$$                     do iqtc = 1,nqtc
c$$$                        do ipop = 1,npopmax
c$$$                           meanqtc(ipop,iqtc) =  meanqtctmp(ipop,iqtc)
c$$$                           sdqtc(ipop,iqtc) =  sdqtctmp(ipop,iqtc)
c$$$                        enddo
c$$$                     enddo
c$$$                  endif
c$$$                  call accept5(nppmax,npopmax,ncolt,
c$$$     &                 nalmax,nal,c,ctmp,f,ftmp,drift,drifttmp)
c$$$                  npop = npop + 1
c$$$               endif
c$$$            endif
c$$$         endif
c$$$*     merge
c$$$      else
c$$$         if(npop .gt. npopmin) then 
c$$$c      write(*,*) 'mort'
c$$$c      write(*,*) 'npop=',npop
c$$$*     tirage de la pop qui meurt
c$$$            ipoprem = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
c$$$*     tirage de la pop hote
c$$$            ipophost = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
c$$$            do while(ipophost .eq. ipoprem)
c$$$               ipophost = 1 + idint(dint(dble(npop)*ggrunif(0.d0,1.d0)))
c$$$            enddo
c$$$*     on range dans la pop d'indice le plus petit
c$$$            if(ipophost .gt. ipoprem) then
c$$$               ii = ipophost
c$$$               ipophost = ipoprem
c$$$               ipoprem = ii
c$$$            endif
c$$$*     recherche des cellules qui vont etre reallouees
c$$$            call who(c,ipoprem,npp,nppmax,cellpop,ncellpop)
c$$$*     recherche des cellules de la pop hote
c$$$            call who(c,ipophost,npp,nppmax,cellpophost,
c$$$     &           ncellpophost)
c$$$*     proposition de reallocation dans la pop ipophost
c$$$            call merging(ipoprem,ipophost,c,ctmp,nppmax,
c$$$     &           ncellpop,cellpop)
c$$$c            write(*,*) 'c=',c(1),c(2),c(3),c(4)
c$$$c            write(*,*) 'ctmp=',ctmp(1),ctmp(2),ctmp(3),ctmp(4) 
c$$$
c$$$*     comptage des alleles sur chaque locus pour c puis ctmp
c$$$            if(((usegeno2 .eq. 1) .or. (usegeno1 .eq. 1)) .or. 
c$$$     &           useql .eq. 1) then
c$$$c               write(*,*) 'dans smfallvar usegeno2=',usegeno2
c$$$               call countnallvar2(nindiv,nlocd,nloch,nql,ncolt,
c$$$     &              nppmax,nal,nalmax,yy,z,ql,n,indcell,c,npopmax,
c$$$     &              usegeno2, usegeno1,useql,ploidy)
c$$$               call countnallvar2(nindiv,nlocd,nloch,nql,ncolt,
c$$$     &              nppmax,nal,nalmax,yy,z,ql,ntmp,indcell,ctmp,npopmax,
c$$$     &              usegeno2, usegeno1,useql,ploidy)
c$$$c               write(*,*) 'n(',ipoprem,',1,1)=',n(ipoprem,1,1)
c$$$c               write(*,*) 'n(',ipophost,',1,1)=',n(ipophost,1,1)
c$$$c               write(*,*) 'ntmp(',ipophost,',1,1)=',ntmp(ipophost,1,1)
c$$$
c$$$*     enleve une pop dans le tableau tmporaire des derives            
c$$$               call remdrift2(ipoprem,ipophost,npop,npopmax,drift,
c$$$     &              drifttmp)
c$$$               rr = (drift(ipophost) - drift(ipoprem))/(2*deltad)
c$$$c               write(*,*) 'drift=',drift
c$$$c               write(*,*) 'drifttmp=',drifttmp
c$$$            endif
c$$$*     calcul du log du ratio  
c$$$*     terme des proposal sur c
c$$$            lratio = dlog(dble(npop)) - 
c$$$     &           dlog(2*dble(ncellpop+ncellpophost+1)) -
c$$$     &           lbico(ncellpop+ncellpophost,ncellpop) 
c$$$c            write(6,*) 'merge apres proposal c lratio=',lratio 
c$$$*     terme des priors sur c
c$$$            lratio = lratio + 
c$$$     &           dble(npp)*(dlog(dble(npop)) - 
c$$$     &           dlog(dble(npop-1)))
c$$$c            write(6,*) 'merge apres prior c lratio=',lratio 
c$$$c           write(*,*) 'term en c lratio =',lratio
c$$$            if(((usegeno2 .eq. 1) .or. (usegeno1 .eq. 1)) .or. 
c$$$     &           useql .eq. 1) then
c$$$*     term des freq
c$$$c$$$  *     term des freq
c$$$               lratio = lratio + 
c$$$     &              lTfallvar(ipophost,ntmp,fa,drifttmp,npopmax,nlocd,
c$$$     &              nloch,nql,ncolt,nal,nalmax,usegeno2,usegeno1,
c$$$     &              useql) -
c$$$     &              lTfallvar(ipophost,n,fa,drift,npopmax,nlocd,
c$$$     &              nloch,nql,ncolt,nal,nalmax,usegeno2,usegeno1,
c$$$     &              useql) - 
c$$$     &              lTfallvar(ipoprem,n,fa,drift,npopmax,nlocd,
c$$$     &              nloch,nql,ncolt,nal,nalmax,usegeno2,usegeno1,
c$$$     &              useql)
c$$$c               write(6,*) 'merge apres prior contrib freq =',lratio 
c$$$c$$$           write(*,*) lTfallvar(ipophost,ntmp,fa,drifttmp,npopmax,nlocd,
c$$$c$$$     &              nloch,nql,ncolt,nal,nalmax,usegeno2,usegeno1,
c$$$c$$$     &              useql)
c$$$c$$$           write(*,*) -lTfallvar(ipophost,n,fa,drift,npopmax,nlocd,
c$$$c$$$     &              nloch,nql,ncolt,nal,nalmax,usegeno2,usegeno1,
c$$$c$$$     &              useql)
c$$$c$$$           write(*,*) -lTfallvar(ipoprem,n,fa,drift,npopmax,nlocd,
c$$$c$$$     &              nloch,nql,ncolt,nal,nalmax,usegeno2,usegeno1,
c$$$c$$$     &              useql)
c$$$
c$$$c          write(*,*) 'apres Tf lratio =',lratio 
c$$$*     terme proposal d
c$$$               lratio = lratio - dlog(2*deltad)
c$$$     &              - .5*rr**2 - dlog(dsqrt(2*3.141593d0))
c$$$*     term prior drift
c$$$               lratio = lratio 
c$$$     &              + (shape1-1)*dlog(drifttmp(ipophost)) 
c$$$     &              + (shape2-1)*dlog(1-drifttmp(ipophost))
c$$$     &              - (shape1-1)*dlog(drift(ipophost))
c$$$     &              - (shape2-1)*dlog(1-drift(ipophost))
c$$$     &              - (shape1-1)*dlog(drift(ipoprem))
c$$$     &              - (shape2-1)*dlog(1-drift(ipoprem))
c$$$     &              - gglgamfn(shape1+shape2)
c$$$     &              + gglgamfn(shape1)+gglgamfn(shape2)  
c$$$c               write(*,*) 'apres T drift lratio=',lratio
c$$$c               write(6,*) 'merge apres prior contrib drift =',lratio 
c$$$            endif
c$$$
c$$$            if(useqtc .eq. 1) then 
c$$$*     propose mean and variance of quantitative variables
c$$$               call propparqvmerge(qtc,nindiv,nqtc,indcell,c,ctmp,
c$$$     &              meanqtc,sdqtc,meanqtctmp,sdqtctmp,nnqtc,sqtc,
c$$$     &              ssqtc,npop,npopmax,nppmax,
c$$$     &              ksiqtc,kappaqtc,alphaqtc,betaqtc,ipophost,
c$$$     &              ipoprem,contriblr)
c$$$*     contrib prior and proposal 
c$$$               lratio = lratio + contriblr   
c$$$*     contrib likelihood quant. variable
c$$$*     argument usegeno2 is passed as 0 to ignore genotypes 
c$$$*     hence avoid having likelihood ratio of genotypes twice 
c$$$               lratio = lratio + lrallvar2(yy,z,ql,
c$$$     &              f,c,ctmp,indcell,indcell,npopmax,nlocd,
c$$$     &              nloch,nql,ncolt,nalmax,nindiv,nppmax,qtc,nqtc,
c$$$     &              meanqtc,sdqtc,meanqtctmp,sdqtctmp,0,0,0,useqtc,
c$$$     &              ploidy)  
c$$$            endif
c$$$c            write(*,*) 'lratio for merge=',lratio
c$$$            lratio = dmin1(0.d0,lratio)
c$$$            alpha  = dexp(lratio)
c$$$            bern = ggrbinom(1.d0,alpha)      
c$$$            
c$$$            if(bern .eq. 1) then
c$$$               if(((usegeno2 .eq. 1) .or. (usegeno1 .eq. 1)) .or. 
c$$$     &              useql .eq. 1) then
c$$$*     propostion du nouveau tableau de freq et de derives
c$$$                  call  remfall(ipoprem,ipophost,
c$$$     &                    npop,npopmax,nlocd,nloch,nql,ncolt,nal,
c$$$     &                 nalmax,f,ftmp,drift,drifttmp,fa,a,ptmp,ntmp,
c$$$     &                 usegeno2,usegeno1,useql)
c$$$               endif
c$$$               if(useqtc .eq. 1) then 
c$$$                  do iqtc = 1,nqtc
c$$$                     do ipop = 1,npopmax
c$$$                        meanqtc(ipop,iqtc) =  meanqtctmp(ipop,iqtc)
c$$$                        sdqtc(ipop,iqtc) =  sdqtctmp(ipop,iqtc)
c$$$                     enddo
c$$$                  enddo
c$$$               endif
c$$$               call accept5(nppmax,npopmax,ncolt,
c$$$     &              nalmax,nal,c,ctmp,f,ftmp,drift,drifttmp)
c$$$               npop = npop - 1
c$$$            endif
c$$$         endif
c$$$      endif
c$$$      end subroutine smfallvar2
c$$$***********************************************************************     
c$$$



***********************************************************************
*     log of term coming from frequencies in a split-merge
*     under CFM
      double precision function lTf(ipop,n,fa,drift,npopmax,nloc,nal,
     &     nalmax)
      implicit none
      integer ipop,n,npopmin,npopmax,nloc,nal,nalmax
      double precision fa,drift
      dimension n(npopmax,nloc,nalmax),fa(nloc,nalmax),nal(nloc), 
     &     drift(npopmax)
      integer iloc,ial,sumn
      double precision gglgamfn
      lTf = 0
      do iloc = 1,nloc
         sumn = 0
         do ial = 1,nal(iloc)
            lTf = lTf + gglgamfn(dble(n(ipop,iloc,ial)) + 
     &           fa(iloc,ial)*(1-drift(ipop))/drift(ipop)) - 
     &           gglgamfn(fa(iloc,ial)*(1-drift(ipop))/drift(ipop))
            sumn = sumn + n(ipop,iloc,ial)
         enddo
         lTf = lTf + gglgamfn((1-drift(ipop))/drift(ipop)) - 
     &        gglgamfn(sumn+(1-drift(ipop))/drift(ipop))
      enddo
      end function lTf
***********************************************************************



***********************************************************************
*     log of term coming from frequencies in a split-merge
*     under CFM
      double precision function lTfallvar(ipop,n,fa,drift,npopmax,nlocd,
     &     nloch,nql,ncolt,nal,nalmax,usegeno2,usegeno1,useql)
      implicit none
      integer ipop,n,npopmin,npopmax,nlocd,nloch,nql,ncolt,nal,nalmax,
     &     usegeno2,usegeno1,useql
      double precision fa,drift
      dimension n(npopmax,ncolt,nalmax),fa(ncolt,nalmax),nal(ncolt), 
     &     drift(npopmax)
      integer iloc,ial,sumn
      double precision gglgamfn
      lTfallvar = 0
      if(usegeno2 .eq. 1) then 
         do iloc = 1,nlocd
            sumn = 0
            do ial = 1,nal(iloc)
               lTfallvar = lTfallvar + gglgamfn(dble(n(ipop,iloc,ial)) + 
     &              fa(iloc,ial)*(1-drift(ipop))/drift(ipop)) - 
     &              gglgamfn(fa(iloc,ial)*(1-drift(ipop))/drift(ipop))
               sumn = sumn + n(ipop,iloc,ial)
            enddo
            lTfallvar = lTfallvar + 
     &           gglgamfn((1-drift(ipop))/drift(ipop)) - 
     &           gglgamfn(sumn+(1-drift(ipop))/drift(ipop))
c            write(*,*) 'lTfallvar=',lTfallvar
         enddo
      endif
      if(usegeno1 .eq. 1) then 
         do iloc = 1,nloch
            sumn = 0
            do ial = 1,nal(nlocd+iloc)
               lTfallvar = lTfallvar + 
     &              gglgamfn(dble(n(ipop,nlocd+iloc,ial)) + 
     &              fa(nlocd+iloc,ial)*(1-drift(ipop))/drift(ipop)) - 
     &              gglgamfn(fa(nlocd+iloc,ial)*
     &              (1-drift(ipop))/drift(ipop))
               sumn = sumn + n(ipop,nlocd+iloc,ial)
            enddo
            lTfallvar = lTfallvar + 
     &           gglgamfn((1-drift(ipop))/drift(ipop)) - 
     &           gglgamfn(sumn+(1-drift(ipop))/drift(ipop))
         enddo
      endif
      if(useql .eq. 1) then 
         do iloc = 1,nql
            sumn = 0
            do ial = 1,nal(nlocd+nloch+iloc)
               lTfallvar = lTfallvar + 
     &              gglgamfn(dble(n(ipop,nlocd+nloch+iloc,ial)) + 
     &              fa(nlocd+nloch+iloc,ial)*
     &              (1-drift(ipop))/drift(ipop)) - 
     &              gglgamfn(fa(nlocd+nloch+iloc,ial)*
     &              (1-drift(ipop))/drift(ipop))
               sumn = sumn + n(ipop,nlocd+nloch+iloc,ial)
            enddo
            lTfallvar = lTfallvar + 
     &           gglgamfn((1-drift(ipop))/drift(ipop)) - 
     &           gglgamfn(sumn+(1-drift(ipop))/drift(ipop))
         enddo
      endif
c      write(*,*) 'usegeno2=',usegeno2
      end function lTfallvar
***********************************************************************




***********************************************************************
*     log of term coming from frequencies in a split-merge
*     under the Dirichlet model
      double precision function lTfd(ipop,n,npopmax,nloc,nal,nalmax)
      implicit none
      integer ipop,n,npopmin,npopmax,nloc,nal,nalmax
      dimension n(npopmax,nloc,nalmax),nal(nloc)
      integer iloc,ial,sumn
      double precision gglgamfn
      lTfd = 0
      do iloc = 1,nloc
         sumn = 0
         do ial = 1,nal(iloc)
            lTfd = lTfd + gglgamfn(dble(n(ipop,iloc,ial)) + 1)
            sumn = sumn + n(ipop,iloc,ial)
         enddo
         lTfd = lTfd + gglgamfn(dble(nal(iloc))) - 
     &        gglgamfn(sumn+dble(nal(iloc)))
      enddo
      end function lTfd
***********************************************************************




***********************************************************************
*     log of term coming from frequencies in a split-merge
*     under the Dirichlet model
      double precision function lTfdallvar(ipop,n,npopmax,nlocd,nloch,
     &     nql,ncolt,nal,nalmax,usegeno2,usegeno1,useql)
      implicit none
      integer ipop,n,npopmin,npopmax,nlocd,nloch,nql,ncolt,nal,nalmax,
     &     usegeno2,usegeno1,useql
      dimension n(npopmax,ncolt,nalmax),nal(ncolt)
      integer iloc,ial,sumn
      double precision gglgamfn
      lTfdallvar = 0
      if(usegeno2 .eq. 1) then
         do iloc = 1,nlocd
            sumn = 0
            do ial = 1,nal(iloc)
               lTfdallvar = lTfdallvar + 
     &              gglgamfn(dble(n(ipop,iloc,ial)) + 1)
               sumn = sumn + n(ipop,iloc,ial)
            enddo
            lTfdallvar = lTfdallvar + gglgamfn(dble(nal(iloc)))- 
     &           gglgamfn(sumn+dble(nal(iloc)))
         enddo
      endif
      if(usegeno1 .eq. 1) then
         do iloc = 1,nloch
            sumn = 0
            do ial = 1,nal(nlocd+iloc)
               lTfdallvar = lTfdallvar + 
     &              gglgamfn(dble(n(ipop,nlocd+iloc,ial)) + 1)
               sumn = sumn + n(ipop,nlocd+iloc,ial)
            enddo
            lTfdallvar = lTfdallvar + gglgamfn(dble(nal(nlocd+iloc)))- 
     &           gglgamfn(sumn+dble(nal(nlocd+iloc)))
         enddo
      endif
      if(useql .eq. 1) then
         do iloc = 1,nql
            sumn = 0
            do ial = 1,nal(nlocd+nloch+iloc)
               lTfdallvar = lTfdallvar + 
     &              gglgamfn(dble(n(ipop,nlocd+nloch+iloc,ial)) + 1)
               sumn = sumn + n(ipop,nlocd+nloch+iloc,ial)
            enddo
            lTfdallvar = lTfdallvar + 
     &           gglgamfn(dble(nal(nlocd+nloch+iloc)))- 
     &           gglgamfn(sumn+dble(nal(nlocd+nloch+iloc)))
         enddo
      endif
      end function lTfdallvar
***********************************************************************




***********************************************************************
      subroutine  postprocesschain2(nxdommax,nydommax,burnin,ninrub,
     &     npopmax,nppmax,nindiv,nlocd,nloch,nql,ncolt,
     &     nal,nalmax,xlim,ylim,dt,nit,
     &     thinning,filenpop,filenpp,fileu,filec,filef,fileperm,filedom,
     &     filemeanqtc,filemeanf,s,u,c,f,pivot,fpiv,fmean,dom,coorddom,
     &     indcel,distcel,order,ordertmp,npopest,usegeno2,usegeno1,
     &     useql,useqtc,nqtc,meanqtc,meanqtcpiv,nitsaved,outorderf)
      implicit none
      character*255 fileu,filec,filenpp,filenpop,filedom,filef,fileperm,
     &      filemeanqtc,filemeanf
      integer nit,thinning,npp,npop,iit,nindiv,nxdommax,
     &     nydommax,npopmax,ipp,nppmax,c,ixdom,iydom,idom,indcel,
     &     ipop,nlocd,nloch,nql,ncolt,nal,nalmax,ijunk,order,ordertmp,
     &     ipopperm,burnin,ninrub,npopest,nnit,iloc,ial,pivot,
     &     usegeno2,usegeno1,useql,useqtc,nqtc,nitsaved,outorderf
      double precision s,u,xlim,ylim,coorddom,dom,domperm,distcel,dt,
     &     f,fpiv,fmean,meanqtc,meanqtcpiv
      integer iitsub,iqtc,iitmodK
*     dimensionnement 
      dimension s(2,nindiv),u(2,nppmax),c(nppmax),
     &     dom(nxdommax*nydommax,npopmax),xlim(2),ylim(2),
     &     domperm(nxdommax*nydommax,npopmax),
     &     coorddom(2,nxdommax*nydommax),indcel(nxdommax*nydommax),
     &     distcel(nxdommax*nydommax),order(npopmax),ordertmp(npopmax),
     &     f(npopmax,ncolt,nalmax),fpiv(npopmax,ncolt,nalmax),
     &     fmean(npopmax,ncolt,nalmax),
     &     nal(ncolt),meanqtc(npopmax,nqtc),meanqtcpiv(npopmax,nqtc),
     &     outorderf(nitsaved,npopmax)

      open(9,file=filenpop)
      open(10,file=filenpp)
      open(11,file=fileu)
      open(12,file=filec)
      if(((usegeno2 .eq. 1) .or. (usegeno1 .eq. 1)) .or. 
     &     (useql .eq. 1)) then  
         open(13,file=filef)
         open(17,file=filemeanf)
      endif
      open(14,file=fileperm)
      open(15,file=filedom)
      if(useqtc .eq. 1) then
         open(16,file=filemeanqtc)
      endif

c      write(6,*) 'debut postproc order=',order

c      write(6,*) 'npopest=', npopest

*     coordonnées de la grille 
      call limit(nindiv,s,xlim,ylim,dt)
      idom = 1
      do ixdom =1,nxdommax
c         write(6,*) 'ixdom=',ixdom
         do iydom=1,nydommax
c            write(6,*) 'iydom=',iydom
            coorddom(1,idom) = xlim(1) + 
     &           dble(ixdom-1)*(xlim(2) - xlim(1))/dble(nxdommax-1)
            coorddom(2,idom) = ylim(1) +
     &           dble(iydom-1)*(ylim(2) - ylim(1))/dble(nydommax-1)
            do ipop=1,npopmax
               dom(idom,ipop) = 0.
               domperm(idom,ipop) = 0.
            enddo
            idom = idom + 1
         enddo
      enddo

****************
*     read value of pivot state   
      nnit = 0
      do iit=1,int(dble(nit)/dble(thinning))
         read(9,*) npop
         if(((usegeno2 .eq. 1) .or. (usegeno1 .eq. 1)) .or. 
     &        (useql .eq. 1)) then
            do iloc=1,ncolt
               do ial=1,nalmax
                  read(13,*) (f(ipop,iloc,ial),ipop=1,npopmax)
               enddo
            enddo
         endif
         if(useqtc .eq. 1) then
            do ipop=1,npopmax
                read(16,*) (meanqtc(ipop,iqtc),iqtc=1,nqtc)
            enddo
         endif
         if((npop .eq. npopest) .and. (iit .gt. burnin)) then 
            nnit = nnit + 1
         endif
         if(nnit .eq. pivot) then
            if(((usegeno2 .eq. 1) .or. (usegeno1 .eq. 1)) .or. 
     &           (useql .eq. 1)) then
               do ipop = 1,npopmax
                  do iloc=1,ncolt
                     do ial=1,nalmax
                        fpiv(ipop,iloc,ial) = f(ipop,iloc,ial)
                     enddo
                  enddo
               enddo
            endif
            if(useqtc .eq. 1) then
               do ipop = 1,npopmax
                  do iqtc = 1,nqtc
                     meanqtcpiv(ipop,iqtc) =  meanqtc(ipop,iqtc)
                  enddo
               enddo
            endif
         endif
      enddo
      rewind 9
      if(((usegeno2 .eq. 1) .or. (usegeno1 .eq. 1)) .or. 
     &     (useql .eq. 1)) then
         rewind 13
      endif
      if(useqtc .eq. 1) then
         rewind 16
      endif


**************
*     relabel wrt to pivot or take pivot as estimator
c      write(6,*) 'relabel'
      nnit = 0
      iitmodK = 0
      do iit=1,int(dble(nit)/dble(thinning))
         read(9,*) npop
         read(10,*) npp
         do ipp=1,nppmax
            read(11,*) u(1,ipp),u(2,ipp)
            read(12,*) c(ipp)
         enddo

*     read current state
         if(((usegeno2 .eq. 1) .or. (usegeno1 .eq. 1)) .or. 
     &        (useql .eq. 1)) then
            do iloc=1,ncolt
               do ial=1,nalmax
                  read(13,*) (f(ipop,iloc,ial),ipop=1,npopmax)
               enddo
            enddo  
         endif
         if(useqtc .eq. 1) then
            do ipop = 1,npopmax
               read(16,*) (meanqtc(ipop,iqtc),iqtc=1,nqtc)
            enddo
         endif

         if((npop .eq. npopest) .and. (iit .gt. burnin)) then 
            nnit = nnit + 1 
            if(npopest .lt. 10) then 
               iitmodK = iitmodK + 1
               call relaballvar(npopmax,nlocd,nloch,nql,ncolt,nalmax,
     &              nal,npop,f,fpiv,meanqtc,meanqtcpiv,nqtc,usegeno2,
     &              usegeno1,useql,useqtc,order,ordertmp)
c               write(14,*) (order(ipop),ipop = 1,npopmax)
               do ipop = 1,npopmax
                  outorderf(iitmodK,ipop) = order(ipop)
               enddo
               if(((usegeno2 .eq. 1) .or. (usegeno1 .eq. 1)) .or. 
     &              (useql .eq. 1)) then
                  do ipop=1,npop
                     do iloc=1,ncolt
                        do ial=1,nalmax
                           fmean(ipop,iloc,ial) = fmean(ipop,iloc,ial) + 
     &                      f(order(ipop),iloc,ial)
                        enddo
                     enddo
                  enddo
               endif
            endif
            if((npopest .lt. 10) .or. (nnit .eq. pivot)) then 
               call calccell(nxdommax*nydommax,coorddom,
     &              npp,nppmax,u,indcel,distcel)
               do idom=1,nxdommax*nydommax
                  ipop = order(c(indcel(idom)))
                  dom(idom,ipop) = dom(idom,ipop) + 1.
               enddo
            endif
         endif
      enddo

*     divide by the right number of iterations
      if(npopest .lt. 10) then 
         do idom=1,nxdommax*nydommax
            do ipop=1,npopmax
               dom(idom,ipop) = dom(idom,ipop)/dble(nnit)
            enddo
         enddo
         if(((usegeno2 .eq. 1) .or. (usegeno1 .eq. 1)) .or. 
     &        (useql .eq. 1)) then
            do ipop=1,npopest
               do iloc=1,ncolt
                  do ial=1,nalmax
                     fmean(ipop,iloc,ial) = 
     &                    fmean(ipop,iloc,ial)/dble(nnit)
                  enddo
               enddo
            enddo
         endif
      endif

c$$$      if(((usegeno2 .eq. 1) .or. (usegeno1 .eq. 1)) .or. 
c$$$     &     (useql .eq. 1)) then
c$$$ 1000    format (300(1x,e15.8,1x))
c$$$         do iloc=1,ncolt
c$$$            do ial=1,nalmax
c$$$               write(17,1000) (sngl(fmean(ipop,iloc,ial)),
c$$$     &              ipop=1,npopmax)
c$$$            enddo
c$$$         enddo  
c$$$      endif

c$$$ 2000 format (1000(e15.5,1x))
c$$$      do idom=1,nxdommax*nydommax
c$$$         write(15,2000) coorddom(1,idom),  coorddom(2,idom), 
c$$$     &        (dom(idom,ipop), ipop=1,npopmax)
c$$$      enddo
c      write(*,*) coorddom
      close(9)
      close(10)
      close(11)
      close(12)
      if(((usegeno2 .eq. 1) .or. (usegeno1 .eq. 1)) .or. 
     &     (useql .eq. 1)) then  
         close(13)
         close(17)
      endif
      close(14)
      close(15)
      if(useqtc .eq. 1) then 
         close(16)
      endif                  
      end subroutine postprocesschain2

c$$$
c$$$***********************************************************************
c$$$      subroutine  postprocessmultchain(nxdommax,nydommax,burnin,ninrub,
c$$$     &     npopmax,nppmax,nindiv,nloc,nal,nalmax,xlim,ylim,dt,nit,
c$$$     &     thinning,nrun,pathall,nchpathall,
c$$$     &     s,u,c,f,pivot,chirunpiv,nchirunpiv,fpiv,dom,coorddom,indcel,
c$$$     &     distcel,order,ordertmp,npopest)
c$$$      implicit none
c$$$      character*255 fileu,filec,filenpp,filenpop,filedom,filef,fileperm,
c$$$     &     pathall,path,chirunpiv,chirun,fileftmp     
c$$$      integer nit,thinning,npp,npop,iit,nindiv,nxdommax,
c$$$     &     nydommax,npopmax,ipp,nppmax,c,ixdom,iydom,idom,indcel,
c$$$     &     ipop,nloc,nal,nalmax,ijunk,order,ordertmp,ipopperm,burnin,
c$$$     &     ninrub,npopest,nnit,iloc,ial,pivot,nrun,irun,nchpathall,
c$$$     &     irunpiv,nchirunpiv,resirun,nchpath
c$$$      double precision s,u,xlim,ylim,coorddom,dom,domperm,distcel,dt,
c$$$     &     f,fpiv
c$$$      integer iitsub
c$$$*     dimensionnement 
c$$$      dimension s(2,nindiv),u(2,nppmax),c(nppmax),
c$$$     &     dom(nxdommax*nydommax,npopmax),xlim(2),ylim(2),
c$$$     &     domperm(nxdommax*nydommax,npopmax),
c$$$     &     coorddom(2,nxdommax*nydommax),indcel(nxdommax*nydommax),
c$$$     &     distcel(nxdommax*nydommax),order(npopmax),ordertmp(npopmax),
c$$$     &     f(npopmax,nloc,nalmax),fpiv(npopmax,nloc,nalmax),nal(nloc)
c$$$
c$$$c      write(6,*) 'debut postproc order=',order
c$$$      fileperm = pathall(1:nchpathall) // "/perm.txt"
c$$$      filedom = pathall(1:nchpathall) // "/proba.pop.membership.txt"
c$$$      open(14,file=fileperm)
c$$$      open(15,file=filedom)
c$$$******************************
c$$$*     coordonnées de la grille 
c$$$      call limit(nindiv,s,xlim,ylim,dt)
c$$$      idom = 1
c$$$      do ixdom =1,nxdommax
c$$$c         write(6,*) 'ixdom=',ixdom
c$$$         do iydom=1,nydommax
c$$$c            write(6,*) 'iydom=',iydom
c$$$            coorddom(1,idom) = xlim(1) + 
c$$$     &           float(ixdom-1)*(xlim(2) - xlim(1))/float(nxdommax-1)
c$$$            coorddom(2,idom) = ylim(1) +
c$$$     &           float(iydom-1)*(ylim(2) - ylim(1))/float(nydommax-1)
c$$$            do ipop=1,npopmax
c$$$               dom(idom,ipop) = 0.
c$$$               domperm(idom,ipop) = 0.
c$$$            enddo
c$$$            idom = idom + 1
c$$$         enddo
c$$$      enddo
c$$$
c$$$****************************************
c$$$*     read frequencies for pivot state  
c$$$*     index of the run containing pivot state
c$$$c      write(*,*) 'en fortran pathall=',pathall
c$$$c      write(*,*) 'en fortran pivot=',pivot
c$$$c      write(*,*) 'en fortran nchpathall=',nchpathall
c$$$c      write(*,*) 'en fortran chirunpiv=',chirunpiv
c$$$c      write(*,*) 'en fortran nchirunpiv=',nchirunpiv
c$$$c      nchpathall = len_trim(pathall)
c$$$c      nchirunpiv =  len_trim(chirunpiv)
c$$$      irunpiv = 1 + 
c$$$     &     int(aint(float(pivot-1)/(float(nit)/float(thinning))))
c$$$      fileftmp = pathall(1:nchpathall) // chirunpiv(1:nchirunpiv)
c$$$c      write(*,*) 'en fortran fileftmp=',fileftmp
c$$$      filef = fileftmp(1:(nchpathall+nchirunpiv)) // "/frequencies.txt"
c$$$c      write(*,*) 'en fortran filef=',filef
c$$$      open(13,file=filef)
c$$$*     index of saved iteration containing pivot  in this file
c$$$      iit = pivot - (irunpiv-1)*float(nit)/float(thinning)
c$$$      do iitsub = 1,iit
c$$$         do iloc=1,nloc
c$$$c     write(6,*) 'iloc=',iloc
c$$$            do ial=1,nalmax
c$$$c     write(6,*) 'ial=',ial
c$$$               read(13,*) (fpiv(ipop,iloc,ial),ipop=1,npopmax)
c$$$            enddo
c$$$         enddo
c$$$      enddo
c$$$      close(13)
c$$$
c$$$
c$$$      nnit = 0
c$$$      do iit=1,int(float(nrun*nit)/float(thinning))
c$$$         irun = 1 + aint(float(iit-1)/(float(nit)/float(thinning)))
c$$$         if(mod(iit-1,int(float(nit)/float(thinning))) .eq. 0) then 
c$$$c            write(*,*) 'iit=',iit
c$$$c            write(*,*) 'irun=',irun
c$$$*     define path to file and  open files (Ref.: book Maryse Ain, p.340) 
c$$$c           write(*,*) 'opening files'
c$$$c            write(*,*) 'pathall=',pathall
c$$$            resirun = irun
c$$$            if(irun .gt. 999) then 
c$$$               path = pathall(1:nchpathall) // 
c$$$     &              char(int(aint(float(irun)/1000)) + ichar('0'))
c$$$               nchpath = nchpathall + 1
c$$$               resirun = resirun - 1000*int(aint(float(irun)/1000))
c$$$               path = path(1:nchpath) // 
c$$$     &              char(int(aint(float(resirun)/100)) + ichar('0'))
c$$$               nchpath = nchpath + 1
c$$$               resirun = resirun - 100*int(aint(float(resirun)/100))
c$$$               path = path(1:nchpath) // 
c$$$     &              char(int(aint(float(resirun)/10)) + ichar('0'))
c$$$               path = path(1:nchpath) //
c$$$     &              char(resirun + ichar('0'))
c$$$               nchpath = nchpath + 1
c$$$c               write(*,*) 'path=',path
c$$$            endif
c$$$            if((irun .gt. 99) .and. (irun .le. 999))then 
c$$$               path = pathall(1:nchpathall) // 
c$$$     &              char(int(aint(float(irun)/100)) + ichar('0'))
c$$$               nchpath = nchpathall + 1
c$$$               resirun = resirun - 100*int(aint(float(irun)/100))
c$$$               path = path(1:nchpath) // 
c$$$     &              char(int(aint(float(resirun)/10)) + ichar('0'))
c$$$               resirun = resirun - 10*int(aint(float(resirun)/10))
c$$$               path = path(1:nchpath) //
c$$$     &              char(resirun + ichar('0'))
c$$$               nchpath = nchpath + 1
c$$$c               write(*,*) 'path=',path
c$$$            endif
c$$$            if((irun .gt. 9) .and. (irun .le. 99)) then 
c$$$               path = pathall(1:nchpathall) // 
c$$$     &              char(int(aint(float(irun)/10)) + ichar('0'))
c$$$               nchpath = nchpathall + 1
c$$$               resirun = resirun - 10*int(aint(float(irun)/10))
c$$$               path = path(1:nchpath) // 
c$$$     &              char(resirun + ichar('0'))
c$$$               nchpath = nchpath + 1
c$$$c               write(*,*) 'path=',path
c$$$            endif
c$$$            if(irun .le. 9) then 
c$$$               path = pathall(1:nchpathall) // char(irun + ichar('0'))
c$$$               nchpath = nchpathall + 1
c$$$            endif
c$$$
c$$$            filenpp = path(1:nchpath) // "/nuclei.numbers.txt"
c$$$            filenpop = path(1:nchpath)// "/populations.numbers.txt"
c$$$            fileu = path(1:nchpath)   // "/coord.nuclei.txt"
c$$$            filec = path(1:nchpath)   // "/color.nuclei.txt"
c$$$            filef = path(1:nchpath)   // "/frequencies.txt"
c$$$            open(9,file=filenpop)
c$$$            open(10,file=filenpp)
c$$$            open(11,file=fileu)
c$$$            open(12,file=filec)
c$$$            open(13,file=filef)
c$$$         endif
c$$$
c$$$
c$$$*     read and relabel  wrt to pivot or take pivot as estimator
c$$$         read(9,*) npop
c$$$         read(10,*) npp
c$$$         do ipp=1,nppmax
c$$$            read(11,*) u(1,ipp),u(2,ipp)
c$$$            read(12,*) c(ipp)
c$$$         enddo
c$$$         do iloc=1,nloc
c$$$            do ial=1,nalmax
c$$$               read(13,*) (f(ipop,iloc,ial),ipop=1,npopmax)
c$$$            enddo
c$$$         enddo  
c$$$         if((npop .eq. npopest) .and. (iit .gt. burnin)) then 
c$$$            nnit = nnit + 1 
c$$$            if(npopest .lt. 10) then 
c$$$               call relabel(npopmax,nloc,nalmax,nal,npopest,f,fpiv,
c$$$     &              order,ordertmp)
c$$$               write(14,*) (order(ipop),ipop = 1,npopmax)
c$$$            endif
c$$$            if((npopest .lt. 10) .or. (iit .eq. pivot)) then 
c$$$               call calccell(nxdommax*nydommax,coorddom,
c$$$     &              npp,nppmax,u,indcel,distcel)
c$$$               do idom=1,nxdommax*nydommax
c$$$                  ipop = order(c(indcel(idom)))
c$$$                  dom(idom,ipop) = dom(idom,ipop) + 1.
c$$$               enddo
c$$$            endif
c$$$         endif
c$$$      
c$$$         if(mod(iit,int(float(nit)/float(thinning))) .eq. 0) then 
c$$$            close(9)
c$$$            close(10)
c$$$            close(11)
c$$$            close(12)
c$$$            close(13)
c$$$         endif
c$$$      enddo
c$$$      if(npopest .lt. 10) then 
c$$$         do idom=1,nxdommax*nydommax
c$$$            do ipop=1,npopmax
c$$$               dom(idom,ipop) = dom(idom,ipop)/float(nnit)
c$$$            enddo
c$$$         enddo
c$$$      endif
c$$$
c$$$ 2000 format (1000(e15.5,1x))
c$$$      do idom=1,nxdommax*nydommax
c$$$         write(15,2000) coorddom(1,idom),  coorddom(2,idom), 
c$$$     &        (dom(idom,ipop), ipop=1,npopmax)
c$$$      enddo  
c$$$      close(14)
c$$$      close(15)     
c$$$      
c$$$      end subroutine postprocessmultchain
c$$$
c$$$


*********************************************************************
*     posterior probability of population membership for individuals
*
      subroutine  pppmindiv2(nindiv,s,npopmax,nppmax,
     &     indcell,distcell,u,c,pmp,filenpop,filenpp,fileu,filec,
     &     fileperm,nit,thinning,burnin,order,npopest,pivot)
      implicit none
 
      integer npopmax,nppmax,nindiv,indcell,npop,npp,c,
     &     nit,burnin,npopest,order(npopmax),thinning,pivot
      double precision pmp,distcell,u,s

      integer iit,ipp,iindiv,ipop,nnit,iitsub
      character*255 filenpop,fileu,filec,filenpp,fileperm
      

      dimension indcell(nindiv),distcell(nindiv),
     &     pmp(nindiv,npopmax),u(2,nppmax),c(nppmax),s(2,nindiv)
c$$$      write(6,*) '      **********************************************'
c$$$      write(6,*) '      *  Computing posterior probabilities          '
c$$$      write(6,*) '      *  of population membership for individuals   '
c$$$      write(6,*) '      **********************************************'
c$$$      write(*,*) 'nindiv=',nindiv
c$$$      write(*,*) 's=',s
c$$$      write(*,*) 'npopmax=',npopmax
c$$$      write(*,*) 'nppmax=',nppmax
c$$$c$$$      write(*,*) 'indcell=',indcell
c$$$c$$$      write(*,*) 'distcell=',distcell
c$$$c$$$      write(*,*) 'u=',u
c$$$c$$$      write(*,*) 'c=',c
c$$$c$$$      write(*,*) 'pmp=',pmp
c$$$      write(*,*) 'filenpp=',filenpp
c$$$      write(*,*) 'fileu=',fileu
c$$$      write(*,*) 'filec=',filec
c$$$      write(*,*) 'nit=',nit
c$$$      write(*,*) 'burnin=',burnin,'\n'
c$$$      write(*,*) 'iit=',iit
c$$$      write(*,*) 'ipp=',ipp
c$$$      write(*,*) 'npp=',npp, '\n'
      open(9,file=filenpop)
      open(10,file=filenpp)
      open(11,file=fileu)
      open(12,file=filec)
      open(13,file=fileperm)
      nnit = 0 
      do iit=1,int(float(nit)/float(thinning))
         read(9,*) npop
         read(10,*) npp
         do ipp=1,nppmax
            read(11,*) u(1,ipp),u(2,ipp)
            read(12,*) c(ipp)
         enddo
         if((npop .eq. npopest) .and. (iit .gt. burnin)) then 
            nnit = nnit + 1 
            if(npopest .lt. 10) then
               read(13,*) (order(ipop),ipop=1,npopmax)
            endif
            if((npopest .lt. 10) .or. (nnit .eq. pivot)) then 
               call calccell(nindiv,s,npp,nppmax,u,indcell,distcell)
               do iindiv=1,nindiv
                  ipop = order(c(indcell(iindiv)))
                  pmp(iindiv,ipop) =  pmp(iindiv,ipop) + 1.
               enddo
            endif
         endif
      enddo
      if(npopest .lt. 10) then
         do iindiv=1,nindiv 
            do ipop=1,npopmax
               pmp(iindiv,ipop) = pmp(iindiv,ipop)/float(nnit)
            enddo
         enddo
      endif
      close(9)
      close(10)
      close(11)
      close(12)
      close(13)

c$$$      write(*,*) 'nindiv=',nindiv
c$$$      write(*,*) 's=',s
c$$$      write(*,*) 'npopmax=',npopmax
c$$$      write(*,*) 'nppmax=',nppmax
c$$$c$$$      write(*,*) 'indcell=',indcell
c$$$c$$$      write(*,*) 'distcell=',distcell
c$$$c$$$      write(*,*) 'u=',u
c$$$c$$$      write(*,*) 'c=',c
c$$$c$$$      write(*,*) 'pmp=',pmp
c$$$      write(*,*) 'filenpp=',filenpp
c$$$      write(*,*) 'fileu=',fileu
c$$$      write(*,*) 'filec=',filec
c$$$      write(*,*) 'nit=',nit
c$$$      write(*,*) 'burnin=',burnin,'\n'
c$$$      write(*,*) 'iit=',iit
c$$$      write(*,*) 'ipp=',ipp
c$$$      write(*,*) 'npp=',npp, '\n'

      end subroutine  pppmindiv2



*********************************************************************
*     posterior probability of population membership for individuals
*
      subroutine  pppmindivmultchain(nindiv,s,npopmax,nppmax,
     &     indcell,distcell,u,c,pmp,pathall,nchpathall,nrun,
     &     nit,thinning,burnin,order,npopest,pivot)
      implicit none
 
      integer npopmax,nppmax,nindiv,indcell,npop,npp,c,nrun,
     &     nit,burnin,npopest,order(npopmax),thinning,pivot,
     &     nchpathall
      double precision pmp,distcell,u,s
      character*255 pathall

      integer iit,ipp,iindiv,ipop,nnit,iitsub,irun,nchpath,resirun
      character*255 path,filenpop,fileu,filec,filenpp,fileperm
      

      dimension indcell(nindiv),distcell(nindiv),
     &     pmp(nindiv,npopmax),u(2,nppmax),c(nppmax),s(2,nindiv)
c$$$      write(6,*) '      **********************************************'
c$$$      write(6,*) '      *  Computing posterior probabilities          '
c$$$      write(6,*) '      *  of population membership for individuals   '
c$$$      write(6,*) '      **********************************************'
c$$$      write(*,*) 'nindiv=',nindiv
c$$$      write(*,*) 's=',s
c$$$      write(*,*) 'npopmax=',npopmax
c$$$      write(*,*) 'nppmax=',nppmax
c$$$c$$$      write(*,*) 'indcell=',indcell
c$$$c$$$      write(*,*) 'distcell=',distcell
c$$$c$$$      write(*,*) 'u=',u
c$$$c$$$      write(*,*) 'c=',c
c$$$c$$$      write(*,*) 'pmp=',pmp
c$$$      write(*,*) 'filenpp=',filenpp
c$$$      write(*,*) 'fileu=',fileu
c$$$      write(*,*) 'filec=',filec
c$$$      write(*,*) 'nit=',nit
c$$$      write(*,*) 'burnin=',burnin,'\n'
c$$$      write(*,*) 'iit=',iit
c$$$      write(*,*) 'ipp=',ipp
c$$$      write(*,*) 'npp=',npp, '\n'

      fileperm = pathall(1:nchpathall) // "/perm.txt"
      open(13,file=fileperm)

      nnit = 0
      do iit=1,int(float(nrun*nit)/float(thinning))
         irun = 1 + aint(float(iit-1)/(float(nit)/float(thinning)))
         if(mod(iit-1,int(float(nit)/float(thinning))) .eq. 0) then 
c            write(*,*) 'iit=',iit
c            write(*,*) 'irun=',irun
*     open files
c           write(*,*) 'opening files'
c            write(*,*) 'pathall=',pathall
*     cf book M Ain, p. 340 
            resirun = irun
            if(irun .gt. 999)then 
               path = pathall(1:nchpathall) // 
     &              char(int(aint(float(irun)/1000)) + ichar('0'))
               nchpath = nchpathall + 1
               resirun = resirun - 1000*int(aint(float(irun)/1000))
               path = path(1:nchpath) // 
     &              char(int(aint(float(resirun)/100)) + ichar('0'))
               nchpath = nchpath + 1
               resirun = resirun - 100*int(aint(float(resirun)/100))
               path = path(1:nchpath) // 
     &              char(int(aint(float(resirun)/10)) + ichar('0'))
               path = path(1:nchpath) //
     &              char(resirun + ichar('0'))
               nchpath = nchpath + 1
c               write(*,*) 'path=',path
            endif
            if((irun .gt. 99) .and. (irun .le. 999))then 
               path = pathall(1:nchpathall) // 
     &              char(int(aint(float(irun)/100)) + ichar('0'))
               nchpath = nchpathall + 1
               resirun = resirun - 100*int(aint(float(irun)/100))
               path = path(1:nchpath) // 
     &              char(int(aint(float(resirun)/10)) + ichar('0'))
               resirun = resirun - 10*int(aint(float(resirun)/10))
               path = path(1:nchpath) //
     &              char(resirun + ichar('0'))
               nchpath = nchpath + 1
c               write(*,*) 'path=',path
            endif
            if((irun .gt. 9) .and. (irun .le. 99)) then 
               path = pathall(1:nchpathall) // 
     &              char(int(aint(float(irun)/10)) + ichar('0'))
               nchpath = nchpathall + 1
               resirun = resirun - 10*int(aint(float(irun)/10))
               path = path(1:nchpath) // 
     &              char(resirun + ichar('0'))
               nchpath = nchpath + 1
c               write(*,*) 'path=',path
            endif
            if(irun .le. 9) then 
               path = pathall(1:nchpathall) // char(irun + ichar('0'))
               nchpath = nchpathall + 1
            endif
            filenpp = path(1:nchpath) // "/nuclei.numbers.txt"
            filenpop = path(1:nchpath) // "/populations.numbers.txt"
            fileu = path(1:nchpath) // "/coord.nuclei.txt"
            filec = path(1:nchpath) // "/color.nuclei.txt"           
            open(9,file=filenpop)
            open(10,file=filenpp)
            open(11,file=fileu)
            open(12,file=filec)
         endif
*     do the real job now
         read(9,*) npop
         read(10,*) npp
         do ipp=1,nppmax
            read(11,*) u(1,ipp),u(2,ipp)
            read(12,*) c(ipp)
         enddo
         if((npop .eq. npopest) .and. (iit .gt. burnin)) then 
            nnit = nnit + 1 
            if(npopest .lt. 10) then
               read(13,*) (order(ipop),ipop=1,npopmax)
            endif
            if((npopest .lt. 10) .or. (iit .eq. pivot)) then 
               call calccell(nindiv,s,npp,nppmax,u,indcell,distcell)
               do iindiv=1,nindiv
                  ipop = order(c(indcell(iindiv)))
                  pmp(iindiv,ipop) =  pmp(iindiv,ipop) + 1.
               enddo
            endif
         endif
*     close files
         if(mod(iit,int(float(nit)/float(thinning))) .eq. 0) then 
            close(9)
            close(10)
            close(11)
            close(12)
         endif
      enddo
      close(13)
c$$$CCCCCCCCCCCCCCCCCCCCCCCCCCCC
c$$$      open(9,file=filenpop)
c$$$      open(10,file=filenpp)
c$$$      open(11,file=fileu)
c$$$      open(12,file=filec)
c$$$      open(13,file=fileperm)
c$$$      nnit = 0 
c$$$      do iit=1,int(float(nit)/float(thinning))
c$$$         read(9,*) npop
c$$$         read(10,*) npp
c$$$         do ipp=1,nppmax
c$$$            read(11,*) u(1,ipp),u(2,ipp)
c$$$            read(12,*) c(ipp)
c$$$         enddo
c$$$         if((npop .eq. npopest) .and. (iit .gt. burnin)) then 
c$$$            nnit = nnit + 1 
c$$$            if(npopest .lt. 10) then
c$$$               read(13,*) (order(ipop),ipop=1,npopmax)
c$$$            endif
c$$$            if((npopest .lt. 10) .or. (iit .eq. pivot)) then 
c$$$               call calccell(nindiv,s,npp,nppmax,u,indcell,distcell)
c$$$               do iindiv=1,nindiv
c$$$                  ipop = order(c(indcell(iindiv)))
c$$$                  pmp(iindiv,ipop) =  pmp(iindiv,ipop) + 1.
c$$$               enddo
c$$$            endif
c$$$         endif
c$$$      enddo
c$$$      if(npopest .lt. 10) then
c$$$         do iindiv=1,nindiv 
c$$$            do ipop=1,npopmax
c$$$               pmp(iindiv,ipop) = pmp(iindiv,ipop)/float(nnit)
c$$$            enddo
c$$$         enddo
c$$$      endif
c$$$      close(9)
c$$$      close(10)
c$$$      close(11)
c$$$      close(12)
c$$$      close(13)
c$$$CCCCCCCCCCCCCCCCCCCCCCCCCCCc

c$$$      write(*,*) 'nindiv=',nindiv
c$$$      write(*,*) 's=',s
c$$$      write(*,*) 'npopmax=',npopmax
c$$$      write(*,*) 'nppmax=',nppmax
c$$$c$$$      write(*,*) 'indcell=',indcell
c$$$c$$$      write(*,*) 'distcell=',distcell
c$$$c$$$      write(*,*) 'u=',u
c$$$c$$$      write(*,*) 'c=',c
c$$$c$$$      write(*,*) 'pmp=',pmp
c$$$      write(*,*) 'filenpp=',filenpp
c$$$      write(*,*) 'fileu=',fileu
c$$$      write(*,*) 'filec=',filec
c$$$      write(*,*) 'nit=',nit
c$$$      write(*,*) 'burnin=',burnin,'\n'
c$$$      write(*,*) 'iit=',iit
c$$$      write(*,*) 'ipp=',ipp
c$$$      write(*,*) 'npp=',npp, '\n'
      end subroutine  pppmindivmultchain


************************************************
      Subroutine relabel(npopmax,nloc,nalmax,nal,npop,f,fpiv,order,
     &     ordertmp)
*     find partition that minimizzes scalar product between f and fpiv
      implicit none
      Integer npopmax,nloc,nalmax,nal,npop,order,ordertmp
      double precision f,fpiv
      Dimension f(npopmax,nloc,nalmax),fpiv(npopmax,nloc,nalmax),
     &     order(npopmax),ordertmp(npopmax),nal(nloc)
      double precision sp,sptmp,spf
      integer ipop
      integer I,I1,J,G,H

      sp = 0 
      do ipop=1,npop
         ordertmp(ipop) = ipop
      enddo
      If (npop.Gt.1) Go To 10
C      Call Sum(npop,ordertmp)
c      Write(*,*) 'ordertmp=',ordertmp
      sptmp = spf(npopmax,nloc,nalmax,nal,npop,f,fpiv,ordertmp)
c      Write(*,*) 'ordertmp=',ordertmp
c      Write(*,*) 'sptmp=',sptmp
      if(sptmp .gt. sp) then
         sp = sptmp
         do ipop=1,npop
            order(ipop) = ordertmp(ipop)
         enddo
      endif
90    Return
10    Continue
      I=npop-2
C      Call Sum(npop,ordertmp)
c      Write(*,*) 'ordertmp=',ordertmp
      sptmp = spf(npopmax,nloc,nalmax,nal,npop,f,fpiv,ordertmp)
c      Write(*,*) 'ordertmp=',ordertmp
c      Write(*,*) 'sptmp=',sptmp
      if(sptmp .gt. sp) then
         sp = sptmp
         do ipop=1,npop
            order(ipop) = ordertmp(ipop)
         enddo
      endif
      G=ordertmp(npop-1)
      H=ordertmp(npop)
      If (G .eq. H) Go To 20
      ordertmp(npop)=G
      ordertmp(npop-1)=H
C      Call Sum(npop,order)
c      Write(*,*) 'order=',order
c      Write(*,*) 'ordertmp=',ordertmp
      sptmp = spf(npopmax,nloc,nalmax,nal,npop,f,fpiv,ordertmp)
c      Write(*,*) 'ordertmp=',ordertmp
c      Write(*,*) 'sptmp=',sptmp
      if(sptmp .gt. sp) then
         sp = sptmp
         do ipop=1,npop
            order(ipop) = ordertmp(ipop)
         enddo
      endif
      ordertmp(npop-1)=G
      ordertmp(npop)=H
20    Continue
      If (I.Eq.0) Go To 90
      H=ordertmp(I)
      I1=I+1
      Do 30 J=I1,npop
      If (ordertmp(J) .Le. H) Go To 30
      ordertmp(I)=ordertmp(J)
      ordertmp(J)=H
      Go To 10
30    Continue
31    Continue
      Do 40 J=I1,npop
      ordertmp(J-1)=ordertmp(J)
40    Continue
      ordertmp(npop)=H
      I=(I-1)
      Go To 20
      End Subroutine relabel
************************************************************************





************************************************************************
*     scalar product of two arrays of frequencies
      double precision function spf(npopmax,nloc,nalmax,nal,npop,
     &     f,fpiv,order)
      implicit none
      Integer npopmax,nloc,nalmax,nal,npop,order
      double precision f,fpiv
      Dimension f(npopmax,nloc,nalmax),fpiv(npopmax,nloc,nalmax),
     &     order(npopmax),nal(nloc)
      integer ipop,iloc,ial
      spf = 0
      do ipop=1,npop
         do iloc = 1,nloc
            do ial=1,nal(iloc)
               spf = spf + f(order(ipop),iloc,ial)*fpiv(ipop,iloc,ial)
            enddo
         enddo
      enddo
      end function spf
************************************************************************
      



************************************************************************
      Subroutine relabelgq(npopmax,nloc,nalmax,nal,npop,f,fpiv,
     &     meanqtc,meanqtcpiv,nqtc,usegeno2,useqtc,order,ordertmp)
*     find partition that minimizzes scalar product 
*     between pivot and current parameter
      implicit none
      Integer npopmax,nloc,nalmax,nal,npop,order,ordertmp,
     &     nqtc,usegeno2,useqtc
      double precision f,fpiv,meanqtc,meanqtcpiv
      Dimension f(npopmax,nloc,nalmax),fpiv(npopmax,nloc,nalmax),
     &     order(npopmax),ordertmp(npopmax),nal(nloc),
     &     meanqtc(npopmax,nqtc),meanqtcpiv(npopmax,nqtc)
      double precision sp,sptmp,spfgq
      integer ipop
      integer I,I1,J,G,H

      sp = 0 
      do ipop=1,npop
         ordertmp(ipop) = ipop
      enddo
      If (npop .Gt. 1) Go To 10
C      Call Sum(npop,ordertmp)
c      Write(*,*) 'ordertmp=',ordertmp
      sptmp = spfgq(npopmax,nloc,nalmax,nal,npop,f,fpiv,
     &     meanqtc,meanqtcpiv,nqtc,usegeno2,useqtc,ordertmp)
c      Write(*,*) 'ordertmp=',ordertmp
c      Write(*,*) 'sptmp=',sptmp
      if(sptmp .gt. sp) then
         sp = sptmp
         do ipop=1,npop
            order(ipop) = ordertmp(ipop)
         enddo
      endif
90    Return
10    Continue
      I=npop-2
C      Call Sum(npop,ordertmp)
c      Write(*,*) 'ordertmp=',ordertmp
      sptmp = spfgq(npopmax,nloc,nalmax,nal,npop,f,fpiv,
     &     meanqtc,meanqtcpiv,nqtc,usegeno2,useqtc,ordertmp)
c      Write(*,*) 'ordertmp =',ordertmp
c      Write(*,*) 'sptmp=',sptmp
      if(sptmp .gt. sp) then
         sp = sptmp
         do ipop=1,npop
            order(ipop) = ordertmp(ipop)
         enddo
      endif
      G=ordertmp(npop-1)
      H=ordertmp(npop)
      If (G .eq. H) Go To 20
      ordertmp(npop)=G
      ordertmp(npop-1)=H
C      Call Sum(npop,order)
c      Write(*,*) 'order=',order
c      Write(*,*) 'ordertmp=',ordertmp
      sptmp = spfgq(npopmax,nloc,nalmax,nal,npop,f,fpiv,
     &     meanqtc,meanqtcpiv,nqtc,usegeno2,useqtc,ordertmp)
c      Write(*,*) 'ordertmp=',ordertmp
c      Write(*,*) 'sptmp=',sptmp
      if(sptmp .gt. sp) then
         sp = sptmp
         do ipop=1,npop
            order(ipop) = ordertmp(ipop)
         enddo
      endif
      ordertmp(npop-1)=G
      ordertmp(npop)=H
20    Continue
      If (I.Eq.0) Go To 90
      H=ordertmp(I)
      I1=I+1
      Do 30 J=I1,npop
      If (ordertmp(J) .Le. H) Go To 30
      ordertmp(I)=ordertmp(J)
      ordertmp(J)=H
      Go To 10
30    Continue
31    Continue
      Do 40 J=I1,npop
      ordertmp(J-1)=ordertmp(J)
40    Continue
      ordertmp(npop)=H
      I=(I-1)
      Go To 20
      End Subroutine relabelgq
************************************************************************




************************************************************************
      Subroutine relaballvar(npopmax,nlocd,nloch,nql,ncolt,nalmax,nal,
     &     npop,f,fpiv,meanqtc,meanqtcpiv,nqtc,usegeno2,usegeno1,
     &     useql,useqtc,order,ordertmp)
*     find partition that minimizzes scalar product 
*     between pivot and current parameter
      implicit none
      Integer npopmax,nlocd,nloch,nql,ncolt,nalmax,nal,npop,order,
     &     ordertmp,nqtc,usegeno2,usegeno1,useql,useqtc
      double precision f,fpiv,meanqtc,meanqtcpiv
      Dimension f(npopmax,ncolt,nalmax),fpiv(npopmax,ncolt,nalmax),
     &     order(npopmax),ordertmp(npopmax),nal(ncolt),
     &     meanqtc(npopmax,nqtc),meanqtcpiv(npopmax,nqtc)
      double precision sp,sptmp,spfallvar
      integer ipop
      integer I,I1,J,G,H

      sp = 0 
      do ipop=1,npop
         ordertmp(ipop) = ipop
      enddo
      If (npop .Gt. 1) Go To 10
C      Call Sum(npop,ordertmp)
c      Write(*,*) 'ordertmp=',ordertmp
      sptmp = spfallvar(npopmax,nlocd,nloch,nql,ncolt,
     &     nalmax,nal,npop,f,fpiv,meanqtc,meanqtcpiv,nqtc,
     &     usegeno2,usegeno1,useql,useqtc,ordertmp)
c      Write(*,*) 'ordertmp=',ordertmp
c      Write(*,*) 'sptmp=',sptmp
      if(sptmp .gt. sp) then
         sp = sptmp
         do ipop=1,npop
            order(ipop) = ordertmp(ipop)
         enddo
      endif
90    Return
10    Continue
      I=npop-2
C      Call Sum(npop,ordertmp)
c      Write(*,*) 'ordertmp=',ordertmp
      sptmp = spfallvar(npopmax,nlocd,nloch,nql,ncolt,
     &     nalmax,nal,npop,f,fpiv,meanqtc,meanqtcpiv,nqtc,
     &     usegeno2,usegeno1,useql,useqtc,ordertmp)
c      Write(*,*) 'ordertmp =',ordertmp
c      Write(*,*) 'sptmp=',sptmp
      if(sptmp .gt. sp) then
         sp = sptmp
         do ipop=1,npop
            order(ipop) = ordertmp(ipop)
         enddo
      endif
      G=ordertmp(npop-1)
      H=ordertmp(npop)
      If (G .eq. H) Go To 20
      ordertmp(npop)=G
      ordertmp(npop-1)=H
C      Call Sum(npop,order)
c      Write(*,*) 'order=',order
c      Write(*,*) 'ordertmp=',ordertmp
      sptmp = spfallvar(npopmax,nlocd,nloch,nql,ncolt,
     &     nalmax,nal,npop,f,fpiv,meanqtc,meanqtcpiv,nqtc,
     &     usegeno2,usegeno1,useql,useqtc,ordertmp)
c      Write(*,*) 'ordertmp=',ordertmp
c      Write(*,*) 'sptmp=',sptmp
      if(sptmp .gt. sp) then
         sp = sptmp
         do ipop=1,npop
            order(ipop) = ordertmp(ipop)
         enddo
      endif
      ordertmp(npop-1)=G
      ordertmp(npop)=H
20    Continue
      If (I.Eq.0) Go To 90
      H=ordertmp(I)
      I1=I+1
      Do 30 J=I1,npop
      If (ordertmp(J) .Le. H) Go To 30
      ordertmp(I)=ordertmp(J)
      ordertmp(J)=H
      Go To 10
30    Continue
31    Continue
      Do 40 J=I1,npop
      ordertmp(J-1)=ordertmp(J)
40    Continue
      ordertmp(npop)=H
      I=(I-1)
      Go To 20
      End Subroutine relaballvar
************************************************************************




************************************************************************
*     scalar product of two vectors of parameters
*     parameters = freq and/or mean of quantit. variables
      double precision function spfgq(npopmax,nloc,nalmax,nal,npop,
     &     f,fpiv,meanqtc,meanqtcpiv,nqtc,usegeno2,useqtc,
     &     order)
      implicit none
      Integer npopmax,nloc,nalmax,nal,npop,order,nqtc,usegeno2,
     &     useqtc
      double precision f,fpiv,meanqtc,meanqtcpiv
      Dimension f(npopmax,nloc,nalmax),fpiv(npopmax,nloc,nalmax),
     &     order(npopmax),nal(nloc),meanqtc(npopmax,nqtc),
     &     meanqtcpiv(npopmax,nqtc)
      integer ipop,iloc,ial,iqtc
      spfgq = 0
     
      if(usegeno2 .eq. 1) then
         do ipop=1,npop
            do iloc = 1,nloc
               do ial=1,nal(iloc)
                  spfgq = spfgq + f(order(ipop),iloc,ial)*
     &                 fpiv(ipop,iloc,ial)
               enddo
            enddo
         enddo
      endif
      if(useqtc .eq. 1) then 
         do ipop=1,npop
            do iqtc = 1,nqtc
               spfgq = spfgq + meanqtc(order(ipop),iqtc)*
     &              meanqtcpiv(ipop,iqtc)
            enddo
         enddo
      endif
     
      end function spfgq
************************************************************************





************************************************************************
*     scalar product of two vectors of parameters
*     parameters = freq and/or mean of quantit. variables
      double precision function spfallvar(npopmax,nlocd,nloch,nql,ncolt,
     &     nalmax,nal,npop,f,fpiv,meanqtc,meanqtcpiv,nqtc,
     &     usegeno2,usegeno1,useql,useqtc,order)
      implicit none
      Integer npopmax,nlocd,nloch,nql,ncolt,nalmax,nal,npop,order,nqtc,
     &     usegeno2,usegeno1,useql,useqtc
      double precision f,fpiv,meanqtc,meanqtcpiv
      Dimension f(npopmax,ncolt,nalmax),fpiv(npopmax,ncolt,nalmax),
     &     order(npopmax),nal(ncolt),meanqtc(npopmax,nqtc),
     &     meanqtcpiv(npopmax,nqtc)
      integer ipop,iloc,ial,iqtc
      spfallvar = 0
      if(usegeno2 .eq. 1) then
         do ipop=1,npop
            do iloc = 1,nlocd
               do ial=1,nal(iloc)
                  spfallvar = spfallvar + f(order(ipop),iloc,ial)*
     &                 fpiv(ipop,iloc,ial)
               enddo
            enddo
         enddo
      endif
      if(usegeno1 .eq. 1) then
         do ipop=1,npop
            do iloc = 1,nloch
               do ial=1,nal(nlocd+iloc)
                  spfallvar = spfallvar + f(order(ipop),nlocd+iloc,ial)*
     &                 fpiv(ipop,nlocd+iloc,ial)
               enddo
            enddo
         enddo
      endif
      if(useql .eq. 1) then
         do ipop=1,npop
            do iloc = 1,nql
               do ial=1,nal(nlocd+nloch+iloc)
                  spfallvar = spfallvar + 
     &                 f(order(ipop),nlocd+nloch+iloc,ial)*
     &                 fpiv(ipop,nlocd+nloch+iloc,ial)
               enddo
            enddo
         enddo
      endif
      if(useqtc .eq. 1) then 
         do ipop=1,npop
            do iqtc = 1,nqtc
               spfallvar = spfallvar + meanqtc(order(ipop),iqtc)*
     &              meanqtcpiv(ipop,iqtc)
            enddo
         enddo
      endif
      end function spfallvar
************************************************************************



************************************************************************
*     computes number of similar entries in a vector of ctrue
*     and permuted entries of vector cest
      integer function simil(nindiv,kest,ctrue,cest,perm)
      implicit none
      integer ctrue,cest,perm,nindiv,kest,c
      dimension ctrue(nindiv),cest(nindiv),perm(kest)
      integer i
      simil=0
      do i=1,nindiv
         if (perm(cest(i)) .eq. ctrue(i)) then
            simil = simil + 1
         endif
      enddo
      end function simil
************************************************************************


************************************************************************
*     find permutation that maximizze the number of similar entries
*     in ctrue and cest
      SUBROUTINE maxsim(N,nindiv,ctrue,cest,permopt,sim)
	implicit none
      INTEGER E,permopt,N
      integer ctrue,cest,nindiv
      dimension cest(nindiv),ctrue(nindiv)
      DIMENSION E(N),permopt(N)
      integer sim,simtmp,simil
      INTEGER G,H,I,J,I1,k
      sim = 0
      simtmp = 0
      do k = 1,N
         E(k) = k 
         permopt(k) = k
      enddo
      IF (N.GT.1) GO TO 10
c      CALL SUM(N,E)
      simtmp = simil(nindiv,N,ctrue,cest,E)
      if(simtmp .gt. sim) then 
         sim = simtmp
         do k = 1,N
            permopt(k) = E(k)
         enddo
      endif
90    RETURN
10    CONTINUE
      I=N-2
c      CALL SUM(N,E)
      simtmp = simil(nindiv,N,ctrue,cest,E)
      if(simtmp .gt. sim) then 
         sim = simtmp
         do k = 1,N
            permopt(k) = E(k)
         enddo
      endif

      G=E(N-1)
      H=E(N)
      IF (G.EQ.H) GO TO 20
      E(N)=G
      E(N-1)=H
c      CALL SUM(N,E)
       simtmp = simil(nindiv,N,ctrue,cest,E)
      if(simtmp .gt. sim) then 
         sim = simtmp
         do k = 1,N
            permopt(k) = E(k)
         enddo
      endif
      E(N-1)=G
      E(N)=H
20    CONTINUE
      IF (I.EQ.0) GO TO 90
      H=E(I)
      I1=I+1
      DO 30 J=I1,N
      IF (E(J) .LE. H) GO TO 30
      E(I)=E(J)
      E(J)=H
      GO TO 10
30    CONTINUE
31    CONTINUE
      DO 40 J=I1,N
      E(J-1)=E(J)
40    CONTINUE
      E(N)=H
      I=(I-1)
      GO TO 20
      END subroutine maxsim
******************************************************************



******************************************************************
*     post processing output of MCMC
*     to plot a sequence of tessellations
      subroutine  tessdyn(npopmax,nppmax,
     &     nxgrid,nygrid,indgrid,npp,c,ccur,nit,burnin,ninrub,
     &     dom,distgrid,coorddom,u,ucur,nindiv,s,xlim,ylim,dt)
      implicit none
 
      integer npopmax,nppmax,nxgrid,nygrid,indgrid,npp,c,
     &     nit,burnin,nindiv,ninrub
      double precision dom,distgrid,coorddom,u,s

      integer iit,ipp,idom,ixdom,iydom,nppcur,ccur,ipop
      double precision xlim(2),ylim(2),ucur,dt,pct

      dimension indgrid(nxgrid*nygrid),distgrid(nxgrid*nygrid),
     &     dom(nxgrid*nygrid,npopmax),coorddom(2,nxgrid*nygrid),
     &  npp(nit),u(2,nppmax,nit),ucur(2,nppmax),
     &     c(nppmax,nit),ccur(nppmax),s(2,nindiv)
      


*     look for smallest rectangle enclosing the spatial domain
      call limit(nindiv,s,xlim,ylim,dt)

*     coordinates of grid
      idom = 1
      do ixdom =1,nxgrid
c         write(6,*) 'ixdom=',ixdom
         do iydom=1,nygrid
c            write(6,*) 'iydom=',iydom
            coorddom(1,idom) = xlim(1) + 
     &           dble(ixdom-1)*(xlim(2) - xlim(1))/dble(nxgrid-1)
            coorddom(2,idom) = ylim(1) +
     &           dble(iydom-1)*(ylim(2) - ylim(1))/dble(nygrid-1)
            do ipop=1,npopmax
               dom(idom,ipop) = 0.d0
            enddo
            idom = idom + 1
         enddo
      enddo

*     sequentially processes states of the chain
      do iit=1,nit
10000    format(f7.3,' %')
         pct = dble(iit)/dble(nit)*100.
         call dblepr('                     ',-1,pct,1)
c         write(6,10000)dble(iit)/dble(nit)*100.
         
         if((iit .gt. burnin) .and. (iit .le. nit-ninrub)) then
            nppcur = npp(iit)
c            write(*,*) 'iit=',iit
c            write(*,*) 'burnin = ', burnin
c            write(*,*) 'ninrub = ', ninrub
            do ipp=1,nppcur
               ucur(1,ipp) = u(1,ipp,iit)
               ucur(2,ipp) = u(2,ipp,iit)
               ccur(ipp) = c(ipp,iit)
            enddo
            
             
c            write(*,*) 'nppcur=',nppcur
c            write(*,*) 'ccur=',ccur
c            write(*,*) 'ucur=',ucur
c            write(*,*) 'avant calccell'
            call calccell(nxgrid*nygrid,coorddom,
     &           nppcur,nppmax,ucur,indgrid,distgrid)
c         write(*,*) 'apres calccell'
         
            do idom=1,nxgrid*nygrid
               ipop = ccur(indgrid(idom))
               
               dom(idom,ipop) =  dom(idom,ipop) + 1.
            enddo
         endif 
      enddo
       

      do idom=1,nxgrid*nygrid 
         do ipop=1,npopmax
            dom(idom,ipop) = dom(idom,ipop)/dble(nit-burnin-ninrub)
         enddo 
      enddo
      end subroutine tessdyn
***********************************************************************
 
  





