*     Calcul du Fst d'apres programme en Turbo Pascal d'Arnaud Estoup
*     le vecteur c contient la variable de classe
*     les effectifs des pop doivent etre donnees en entree
 
      subroutine ggfst(nindiv,nloc,nloc2,nall,npop,effcl,z,c,
     &     tabindiv,kk,Fistot,Fsttot,Fittot)
      implicit none
      integer nindiv,nloc,nloc2,nall,npop,effcl,z,c,tabindiv,kk

      double precision Fsttot,Fittot,Fistot
      
      dimension nall(nloc),z(nindiv,nloc2),c(nindiv),effcl(npop),
     &     tabindiv(nindiv,npop),kk(npop)

      integer iloc,ipop,iall,iindiv
      integer g1,g2,k

      double precision s1,s2,s3,s1l,s2l,s3l,ni,sni,sni2,sniA,sniAA,
     &     s2A,nA,AA,nc,MSG,MSI,MSP,s2G,s2I,s2P,Fst,Fit,Fis
  

*     Recherche des indices des indiv de chaque pop 
      do ipop=1,npop
         kk(ipop) = 1
      enddo
      do iindiv=1,nindiv
         tabindiv(kk(c(iindiv)),c(iindiv)) = iindiv
         kk(c(iindiv)) = kk(c(iindiv)) + 1
      enddo
c$$$      write(*,*) (tabindiv(iindiv,1),iindiv=1,nindiv)
c$$$      write(*,*) (tabindiv(iindiv,2),iindiv=1,nindiv)
c$$$      write(*,*) 'z=',z
c$$$      write(*,*) 'c=',c         

      s1 = 0.
      s2 = 0.
      s3 = 0.
      do iloc=1,nloc
c         write(*,*) 'iloc=', iloc
         s1l = 0.
         s2l = 0.
         s3l = 0.
         do iall=1,nall(iloc)
c            write(*,*) 'iall=', iall
            sni = 0.
            sni2 = 0.
            sniA = 0.
            sniAA = 0.
            s2A = 0.
c            k = 1
            do ipop=1,npop
               ni = 0.
               nA = 0.
               AA = 0.
               do iindiv=1,idint(dble(effcl(ipop)))
                  k = tabindiv(iindiv,ipop)
                  g1 = z(k,2*(iloc-1)+1)
                  g2 = z(k,2*(iloc-1)+2)
                  if((g1 .eq. iall) .and. (g2 .eq. iall)) then 
                     AA = AA + 1.
                  endif
                  if(g1 .eq. iall) then 
                     nA = nA + 1. 
                  endif
                  if(g2 .eq. iall) then 
                     nA = nA + 1. 
                  endif
C     test on missing values added by Gilles on 7/4/08
C     corrected on 29/07/08
                  if((z(k,2*(iloc-1)+1) .ne. -999) .and. 
     &                 (z(k,2*(iloc-1)+2) .ne. -999)) then 
                     ni = ni + 1
                  endif
c                    write(*,*) 'AA,nA=',AA,nA
               enddo 
c               k = k+ni
               sniA = sniA + nA
               sniAA = sniAA + AA
               sni = sni + ni 
               sni2 = sni2 + ni*ni
               s2A = s2A + nA*nA/(2*ni)
c            write(*,*) 'sniA,sniAA,sni,sni2,s2=',sniA,sniAA,sni,sni2,s2A
            enddo
            if(sni .gt. 0) then 
               nc = (sni-sni2/sni)
            else
               nc = 0
            endif
            if(sni*nc .gt. 0) then
               MSG = (0.5*sniA-sniAA)/sni
               MSI = (0.5*sniA+sniAA-s2A)/(sni-dble(npop))
               MSP = (s2A-0.5*(sniA**2)/sni)/(dble(npop)-1)
               s2G = MSG
               s2I = 0.5*(MSI-MSG)
               s2P = (MSP-MSI)/(2*nc)
               s1l = s1l + s2P
               s2l = s2l + s2P + s2I
               s3l = s3l + s2P + s2I + s2G
            endif
c            write(*,*) 'variables=',nc,MSG,MSI,s2G,s2I,s2P,s1l,s2l,s3l
         enddo
         Fst = s1l/s3l
         Fit = s2l/s3l
         Fis = (Fit-Fst)/(1-Fst)
         s1 = s1 + s1l
         s2 = s2 + s2l
         s3 = s3 + s3l
c         write(*,*) 'variables=',Fst,Fit,Fis,s1,s2,s3
      enddo
      Fsttot = s1/s3
      Fittot = S2/s3
      Fistot = (Fittot-Fsttot)/(1-Fsttot)
      
      end subroutine ggfst
      




 
c$$$
c$$$
c$$$
c$$$*     Calcul du Fst d'apres programme en Turbo Pascal d'Arnaud Estoup
c$$$*     le vecteur c contient la variable de classe
c$$$*     les effectifs des pop doivent etre donnees en entree
c$$$ 
c$$$      subroutine ggfst(nindiv,nloc,nloc2,nall,npop,effcl,z,c,
c$$$     &     tabindiv,kk,Fistot,Fsttot,Fittot)
c$$$      implicit none
c$$$      integer nindiv,nloc,nloc2,nall,npop,effcl,z,c,tabindiv,kk
c$$$
c$$$      double precision Fsttot,Fittot,Fistot
c$$$      
c$$$      dimension nall(nloc),z(nindiv,nloc2),c(nindiv),effcl(npop),
c$$$     &     tabindiv(nindiv,npop),kk(npop)
c$$$
c$$$      integer iloc,ipop,iall,iindiv
c$$$      integer g1,g2,k
c$$$
c$$$      double precision s1,s2,s3,s1l,s2l,s3l,ni,sni,sni2,sniA,sniAA,
c$$$     &     s2A,nA,AA,nc,MSG,MSI,MSP,s2G,s2I,s2P,Fst,Fit,Fis
c$$$  
c$$$
c$$$*     Recherche des indices des indiv de chaque pop 
c$$$      do ipop=1,npop
c$$$         kk(ipop) = 1
c$$$      enddo
c$$$      do iindiv=1,nindiv
c$$$         tabindiv(kk(c(iindiv)),c(iindiv)) = iindiv
c$$$         kk(c(iindiv)) = kk(c(iindiv)) + 1
c$$$      enddo
c$$$c$$$      write(*,*) (tabindiv(iindiv,1),iindiv=1,nindiv)
c$$$c$$$      write(*,*) (tabindiv(iindiv,2),iindiv=1,nindiv)
c$$$c$$$      write(*,*) 'z=',z
c$$$c$$$      write(*,*) 'c=',c         
c$$$
c$$$      s1 = 0.
c$$$      s2 = 0.
c$$$      s3 = 0.
c$$$      do iloc=1,nloc
c$$$c         write(*,*) 'iloc=', iloc
c$$$         s1l = 0.
c$$$         s2l = 0.
c$$$         s3l = 0.
c$$$         do iall=1,nall(iloc)
c$$$c            write(*,*) 'iall=', iall
c$$$            sni = 0.
c$$$            sni2 = 0.
c$$$            sniA = 0.
c$$$            sniAA = 0.
c$$$            s2A = 0.
c$$$c            k = 1
c$$$            do ipop=1,npop
c$$$               ni = dble(effcl(ipop))
c$$$               nA = 0.
c$$$               AA = 0.
c$$$               do iindiv=1,idint(ni) 
c$$$                  k = tabindiv(iindiv,ipop)
c$$$                  g1 = z(k,2*(iloc-1)+1)
c$$$                  g2 = z(k,2*(iloc-1)+2)
c$$$                  if((g1 .eq. iall) .and. (g2 .eq. iall)) then 
c$$$                     AA = AA + 1.
c$$$                  endif
c$$$                  if(g1 .eq. iall) then 
c$$$                     nA = nA + 1. 
c$$$                  endif
c$$$                  if(g2 .eq. iall) then 
c$$$                     nA = nA + 1. 
c$$$                  endif
c$$$C     test on missing values added by Gilles on 7/4/08
c$$$                  if((z(k,2*(iloc-1)+1) .eq. -999) .and. 
c$$$     &                 (z(k,2*(iloc-1)+2) .eq. -999)) then 
c$$$                     ni = ni + 1
c$$$                  endif
c$$$c                    write(*,*) 'AA,nA=',AA,nA
c$$$               enddo 
c$$$c               k = k+ni
c$$$               sniA = sniA + nA
c$$$               sniAA = sniAA + AA
c$$$               sni = sni + ni 
c$$$               sni2 = sni2 + ni*ni
c$$$               s2A = s2A + nA*nA/(2*ni)
c$$$c            write(*,*) 'sniA,sniAA,sni,sni2,s2=',sniA,sniAA,sni,sni2,s2A
c$$$            enddo
c$$$c$$$            nc = (sni-sni2/sni)/(npop-1)
c$$$c$$$            MSG = (0.5*sniA-sniAA)/sni
c$$$c$$$            MSI = (0.5*sniA+sniAA-s2A)/(sni-npop)
c$$$c$$$            MSP = (s2A-0.5*(sniA**2)/sni)/(npop-1) 
c$$$            nc = (sni-sni2/sni)/(dble(npop)-1)
c$$$            MSG = (0.5*sniA-sniAA)/sni
c$$$            MSI = (0.5*sniA+sniAA-s2A)/(sni-dble(npop))
c$$$            MSP = (s2A-0.5*(sniA**2)/sni)/(dble(npop)-1)
c$$$            s2G = MSG
c$$$            s2I = 0.5*(MSI-MSG)
c$$$            s2P = (MSP-MSI)/(2*nc)
c$$$            s1l = s1l + s2P
c$$$            s2l = s2l + s2P + s2I
c$$$            s3l = s3l + s2P + s2I + s2G
c$$$c            write(*,*) 'variables=',nc,MSG,MSI,s2G,s2I,s2P,s1l,s2l,s3l
c$$$         enddo
c$$$         Fst = s1l/s3l
c$$$         Fit = s2l/s3l
c$$$         Fis = (Fit-Fst)/(1-Fst)
c$$$         s1 = s1 + s1l
c$$$         s2 = s2 + s2l
c$$$         s3 = s3 + s3l
c$$$c         write(*,*) 'variables=',Fst,Fit,Fis,s1,s2,s3
c$$$      enddo
c$$$      Fsttot = s1/s3
c$$$      Fittot = S2/s3
c$$$      Fistot = (Fittot-Fsttot)/(1-Fsttot)
c$$$      
c$$$      end subroutine ggfst
c$$$      
c$$$
c$$$
c$$$
c$$$
c$$$ 
c$$$
c$$$
c$$$
c$$$
c$$$
c$$$
c$$$
c$$$
c$$$
c$$$
c$$$ 
c$$$
c$$$
c$$$
c$$$
c$$$
c$$$
c$$$
c$$$
c$$$
c$$$
c$$$ 
c$$$
c$$$
c$$$
c$$$
c$$$
c$$$
c$$$
c$$$
c$$$
c$$$
c$$$
c$$$
c$$$
c$$$
c$$$
c$$$ 
c$$$
c$$$



