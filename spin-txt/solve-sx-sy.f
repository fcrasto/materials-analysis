c********************************************
c Last update 05/2018 Felipe Crasto de Lima
c********************************************
      Program spintxt 

      implicit none

      real max_length, length, points
      parameter (max_length = 800, length = 1000, points = 70)
      character*3 nada
      integer l, ii
      integer i, j, k, sud
      integer svec, clin, nskip, LORBIT
      parameter (svec = 1500)
      integer iion(svec), iband(svec,svec)
      real a(svec,svec)   ! numero de auto-valores up
      real orb(10,svec,svec)
      real norm(svec,svec)
cccccccccccccccccccccccccccccccc
      real orbx(10,svec,svec)
      real orby(10,svec,svec)
      real orbz(10,svec,svec)
cccccccccccccccccccccccccccccccc
      real e_fermi 
      real ene    (svec,svec), sPnlm  (svec,svec)
      real pyPnlm (svec,svec), pzPnlm (svec,svec), pxPnlm (svec,svec)
      real dxyPnlm(svec,svec), dyzPnlm(svec,svec), dz2Pnlm(svec,svec)
      real dxzPnlm(svec,svec), dx2Pnlm(svec,svec), totPnlm(svec,svec)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real sPnlmx  (svec,svec)
      real pyPnlmx (svec,svec), pzPnlmx (svec,svec), pxPnlmx (svec,svec)
      real dxyPnlmx(svec,svec), dyzPnlmx(svec,svec), dz2Pnlmx(svec,svec)
      real dxzPnlmx(svec,svec), dx2Pnlmx(svec,svec), totPnlmx(svec,svec)
      real sPnlmy  (svec,svec)
      real pyPnlmy (svec,svec), pzPnlmy (svec,svec), pxPnlmy (svec,svec)
      real dxyPnlmy(svec,svec), dyzPnlmy(svec,svec), dz2Pnlmy(svec,svec)
      real dxzPnlmy(svec,svec), dx2Pnlmy(svec,svec), totPnlmy(svec,svec)
      real sPnlmz  (svec,svec)
      real pyPnlmz (svec,svec), pzPnlmz (svec,svec), pxPnlmz (svec,svec)
      real dxyPnlmz(svec,svec), dyzPnlmz(svec,svec), dz2Pnlmz(svec,svec)
      real dxzPnlmz(svec,svec), dx2Pnlmz(svec,svec), totPnlmz(svec,svec)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real XNada(svec),YNada(svec), ZNada(svec)
      real AA(3,3), XXNada, YYNada, ZZNada
      real kp(3), kpt(svec)
      real maxtot, maxcs, maxcpx, maxcpy, maxcpz
      real maxcdxy, maxcdyz, maxcdz2, maxcdxz, maxcdx2
      real maxnorm
      real tartot, targs, targpx, targpy, targpz, perc
      real targdxy, targdyz, targdz2, targdxz, targdx2
      real targnorm
      character c1*15, c2*20, c3*19, c5*19
      character c7*4,  c8*9,  c9*7, kkk*1
      integer   nkpt, ibin, sk
      integer ni, nk, nb, ij, naqui, ispin
      integer n0, nf, n02, nf2, n03, nf3, n04, nf4
      integer n05, nf5
      logical orb_s
      logical orb_px, orb_py, orb_pz
      logical orb_dxy, orb_dyz, orb_dz2, orb_dxz, orb_dx2

      open(90, file = 'inp-spin')

      read(90,*) ispin
      read(90,*) clin  
      read(90,*) LORBIT
      read(90,*) e_fermi
      read(90,*) n0
      read(90,*) nf
      read(90,*) n02
      read(90,*) nf2
      read(90,*) n03
      read(90,*) nf3
      read(90,*) n04
      read(90,*) nf4
      read(90,*) n05
      read(90,*) nf5
      read(90,*) perc
      read(90,*) orb_s
      read(90,*) orb_px
      read(90,*) orb_py
      read(90,*) orb_pz
      read(90,*) orb_dxy
      read(90,*) orb_dyz
      read(90,*) orb_dz2
      read(90,*) orb_dxz
      read(90,*) orb_dx2
      do i = 1, 3
        read(90,*) (AA(i,j), j = 1, 3)
      enddo
   
      do i = 1, 1000
        iion(i) = 1000
      enddo
       
      open (10, file = 'PROCAR', ACCESS = 'SEQUENTIAL')
      open (13, file = 'grid-sx-sy.dat')
cccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do i = n0, nf
        iion(i) = i
      enddo

      if (n02 .ne. 0 .and. nf2 .ne. 0) then
        do i = n02, nf2
          iion(i) = i
        enddo
      endif
      if (n03 .ne. 0 .and. nf3 .ne. 0) then
        do i = n03, nf3
          iion(i) = i
        enddo
      endif
      if (n04 .ne. 0 .and. nf4 .ne. 0) then
        do i = n04, nf4
          iion(i) = i
        enddo
      endif
      if (n05 .ne. 0 .and. nf5 .ne. 0) then
        do i = n05, nf5
          iion(i) = i
        enddo
      endif
      
      read(10,*) nada

      DO SUD = 1, ISPIN ! Spin up (1) or down (2)

      read(10,100) c1, nk, c2, nb, c3, ni
      write(6,100) c1, nk, c2, nb, c3, ni

      if (clin .eq. 1) then
        nskip = 5*(ni + 1)
       else
        nskip = 2*(ni + 1)  
      endif

      orb(:,:,:) = 0.
        
      do ij = 1, nk
        read(10,*)
        read(10,*) c9, nkpt, kkk, kp(1), kp(2), kp(3)
c       write(6,*) c9, nkpt, kkk, kp(1), kp(2), kp(3)            
        if (mod(ij, 5) .eq. 0)  write(*,*) c9, nkpt
        read(10,*)
        XNada(ij) = kp(1)*AA(1,1) + kp(2)*AA(2,1) + kp(3)*AA(3,1)
        YNada(ij) = kp(1)*AA(1,2) + kp(2)*AA(2,2) + kp(3)*AA(3,2)
        ZNada(ij) = kp(1)*AA(1,3) + kp(2)*AA(2,3) + kp(3)*AA(3,3)
        if (ij .eq. 1) then
            kpt(ij) = 0.0
           else
            XXNada = XNada(ij-1) - XNada(ij)
            YYNada = YNada(ij-1) - YNada(ij)
            ZZNada = ZNada(ij-1) - ZNada(ij)
            kpt(ij) = sqrt(XXNada**2 + YYNada**2 + ZZNada**2)
            kpt(ij) = kpt(ij-1) + kpt(ij)
        endif
c       write(6,*) kpt(ij)
c Loop das bandas.
        do j = 1, nb
          read(10,104) c7, iband(j,ij), c8, ene(j,ij), c5
          !write(6,104) c7, iband(j,ij), c8, ene(j,ij), c5
          do l = 1, 2
            read(10,*)
          enddo ! l
c          if(clin .eq. 1)then
c          do sk = 1, 3*(ni+1)
c          read(10,*)
c          enddo
c          endif
          do k = 1, ni
c Este eh o projetado Pnlm|Ylm>
            read(10,*) naqui, sPnlm(j,ij),
     :         pyPnlm (j,ij), pzPnlm (j,ij), pxPnlm (j,ij),
     :         dxyPnlm(j,ij), dyzPnlm(j,ij), dz2Pnlm(j,ij),
     :         dxzPnlm(j,ij), dx2Pnlm(j,ij), totPnlm(j,ij)

            norm(j,ij) = norm(j,ij) + totPnlm(j,ij)


            if (k .eq. iion(k)) then
              if (ij .eq. 1 .and. j .eq. 1) write(6,*)'ion = ', k
              a(j,ij) = ene(j,ij) - e_fermi
              orb(1,j,ij) = orb(1,j,ij) + totPnlm(j,ij) 
              if (orb_s) then
                orb(2,j,ij) = orb(2,j,ij) + sPnlm(j,ij) 
              endif
              if (orb_px) then
                orb(3,j,ij) = orb(3,j,ij) + pxPnlm(j,ij) 
              endif
              if (orb_py) then
                orb(4,j,ij) = orb(4,j,ij) + pyPnlm(j,ij) 
              endif
              if (orb_pz) then
                orb(5,j,ij) = orb(5,j,ij) + pzPnlm(j,ij) 
              endif
              if (orb_dxy) then
                orb(6,j,ij) = orb(6, j,ij) + dxyPnlm(j,ij) 
              endif
              if (orb_dxz) then
                orb(7,j,ij) = orb(7, j,ij) + dyzPnlm(j,ij) 
              endif
              if (orb_dz2) then
                orb(8,j,ij) = orb(8, j,ij) + dz2Pnlm(j,ij) 
              endif
              if (orb_dxz) then
                orb(9,j,ij) = orb(9, j,ij) + dxzPnlm(j,ij) 
              endif
              if (orb_dx2) then
                orb(10,j,ij) = orb(10,j,ij) + dx2Pnlm(j,ij) 
              endif 
            endif
          enddo ! k
          if (LORBIT .eq. 12) then
            do k = 1, nskip 
              read(10,*)
            enddo
           else
            read(10,*)
          endif 
ccccccccccccc Felipe cccc
c          elseif(clin .eq. 1)then
c          do sk = 1, 3*(ni+1)
c          read(10,*)
c          enddo
c          endif
cccccccccccccccccccccccc spin texture cccccccccccccccccccc
          do k = 1, ni
c Este eh o projetado Pnlmx|Ylm>
            read(10,*) naqui, sPnlmx(j,ij),
     :         pyPnlmx (j,ij), pzPnlmx (j,ij), pxPnlmx (j,ij),
     :         dxyPnlmx(j,ij), dyzPnlmx(j,ij), dz2Pnlmx(j,ij),
     :         dxzPnlmx(j,ij), dx2Pnlmx(j,ij), totPnlmx(j,ij)
            if (k .eq. iion(k)) then
              if (ij .eq. 1 .and. j .eq. 1) write(6,*)'ion = ', k
c              a(j,ij) = ene(j,ij) - e_fermi
              orbx(1,j,ij) = orbx(1,j,ij) + totPnlmx(j,ij)
              if (orb_s) then
                orbx(2,j,ij) = orbx(2,j,ij) + sPnlmx(j,ij)
              endif
              if (orb_px) then
                orbx(3,j,ij) = orbx(3,j,ij) + pxPnlmx(j,ij)
              endif
              if (orb_py) then
                orbx(4,j,ij) = orbx(4,j,ij) + pyPnlmx(j,ij)
              endif
              if (orb_pz) then
                orbx(5,j,ij) = orbx(5,j,ij) + pzPnlmx(j,ij)
              endif
              if (orb_dxy) then
                orbx(6,j,ij) = orbx(6, j,ij) + dxyPnlmx(j,ij)
              endif
              if (orb_dxz) then
                orbx(7,j,ij) = orbx(7, j,ij) + dyzPnlmx(j,ij)
              endif
              if (orb_dz2) then
                orbx(8,j,ij) = orbx(8, j,ij) + dz2Pnlmx(j,ij)
              endif
              if (orb_dxz) then
                orbx(9,j,ij) = orbx(9, j,ij) + dxzPnlmx(j,ij)
              endif
              if (orb_dx2) then
                orbx(10,j,ij) = orbx(10,j,ij) + dx2Pnlmx(j,ij)
              endif
            endif
          enddo ! k
          read(10,*)
          do k = 1, ni
c Este eh o projetado Pnlmx|Ylm>
            read(10,*) naqui, sPnlmy(j,ij),
     :         pyPnlmy (j,ij), pzPnlmy (j,ij), pxPnlmy (j,ij),
     :         dxyPnlmy(j,ij), dyzPnlmy(j,ij), dz2Pnlmy(j,ij),
     :         dxzPnlmy(j,ij), dx2Pnlmy(j,ij), totPnlmy(j,ij)
            if (k .eq. iion(k)) then
              if (ij .eq. 1 .and. j .eq. 1) write(6,*)'ion = ', k
c              a(j,ij) = ene(j,ij) - e_fermi
              orby(1,j,ij) = orby(1,j,ij) + totPnlmy(j,ij)
              if (orb_s) then
                orby(2,j,ij) = orby(2,j,ij) + sPnlmy(j,ij)
              endif
              if (orb_px) then
                orby(3,j,ij) = orby(3,j,ij) + pxPnlmy(j,ij)
              endif
              if (orb_py) then
                orby(4,j,ij) = orby(4,j,ij) + pyPnlmy(j,ij)
              endif
              if (orb_pz) then
                orby(5,j,ij) = orby(5,j,ij) + pzPnlmy(j,ij)
              endif
              if (orb_dxy) then
                orby(6,j,ij) = orby(6, j,ij) + dxyPnlmy(j,ij)
              endif
              if (orb_dxz) then
                orby(7,j,ij) = orby(7, j,ij) + dyzPnlmy(j,ij)
              endif
              if (orb_dz2) then
                orby(8,j,ij) = orby(8, j,ij) + dz2Pnlmy(j,ij)
              endif
              if (orb_dxz) then
                orby(9,j,ij) = orby(9, j,ij) + dxzPnlmy(j,ij)
              endif
              if (orb_dx2) then
                orby(10,j,ij) = orby(10,j,ij) + dx2Pnlmy(j,ij)
              endif
            endif
          enddo ! k
          read(10,*)
          do k = 1, ni
c Este eh o projetado Pnlmx|Ylm>
            read(10,*) naqui, sPnlmz(j,ij),
     :         pyPnlmz (j,ij), pzPnlmz (j,ij), pxPnlmz (j,ij),
     :         dxyPnlmz(j,ij), dyzPnlmz(j,ij), dz2Pnlmz(j,ij),
     :         dxzPnlmz(j,ij), dx2Pnlmz(j,ij), totPnlmz(j,ij)
            if (k .eq. iion(k)) then
              if (ij .eq. 1 .and. j .eq. 1) write(6,*)'ion = ', k
c              a(j,ij) = ene(j,ij) - e_fermi
              orbz(1,j,ij) = orbz(1,j,ij) + totPnlmz(j,ij)
              if (orb_s) then
                orbz(2,j,ij) = orbz(2,j,ij) + sPnlmz(j,ij)
              endif
              if (orb_px) then
                orbz(3,j,ij) = orbz(3,j,ij) + pxPnlmz(j,ij)
              endif
              if (orb_py) then
                orbz(4,j,ij) = orbz(4,j,ij) + pyPnlmz(j,ij)
              endif
              if (orb_pz) then
                orbz(5,j,ij) = orbz(5,j,ij) + pzPnlmz(j,ij)
              endif
              if (orb_dxy) then
                orbz(6,j,ij) = orbz(6, j,ij) + dxyPnlmz(j,ij)
              endif
              if (orb_dxz) then
                orbz(7,j,ij) = orbz(7, j,ij) + dyzPnlmz(j,ij)
              endif
              if (orb_dz2) then
                orbz(8,j,ij) = orbz(8, j,ij) + dz2Pnlmz(j,ij)
              endif
              if (orb_dxz) then
                orbz(9,j,ij) = orbz(9, j,ij) + dxzPnlmz(j,ij)
              endif
              if (orb_dx2) then
                orbz(10,j,ij) = orbz(10,j,ij) + dx2Pnlmz(j,ij)
              endif
            endif
          enddo ! k
          read(10,*)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          read(10,*)
        enddo ! j
      enddo ! ij

c     printing output
      maxnorm = 0.

      maxcs   = 0.
      maxcpx  = 0.
      maxcpy  = 0.
      maxcpz  = 0.
      maxcdxy = 0.
      maxcdyz = 0.
      maxcdz2 = 0.
      maxcdxz = 0.
      maxcdx2 = 0.


      do ij = 1, nk
        do j = 1, nb
            if(orb(1,j,ij) .gt. maxtot) maxtot = orb(1,j,ij)
              tartot = perc*maxtot
          if (orb_s) then
            if (orb(2,j,ij) .gt. maxcs) maxcs = orb(2,j,ij) 
            targs = perc*maxcs
          endif 
          if (orb_px) then
            if (orb(3,j,ij) .gt. maxcpx) maxcpx = orb(3,j,ij) 
            targpx = perc*maxcpx
          endif 
          if (orb_py) then
            if (orb(4,j,ij) .gt. maxcpy) maxcpy = orb(4,j,ij) 
            targpy = perc*maxcpy
          endif 
          if (orb_pz) then
            if (orb(5,j,ij) .gt. maxcpz) maxcpz = orb(5,j,ij) 
            targpz = perc*maxcpz
          endif 
          if (orb_dxy) then
            if (orb(6,j,ij) .gt. maxcdxy) maxcdxy = orb(6,j,ij) 
            targdxy = perc*maxcdxy
          endif 
          if (orb_dyz) then
            if (orb(7,j,ij) .gt. maxcdyz) maxcdyz = orb(7,j,ij) 
            targdyz = perc*maxcdyz
          endif 
          if (orb_dz2) then
            if (orb(8,j,ij) .gt. maxcdz2) maxcdz2 = orb(8,j,ij) 
            targdz2 = perc*maxcdz2
          endif 
          if (orb_dxz) then
            if (orb(9,j,ij) .gt. maxcdxz) maxcdxz = orb(9,j,ij) 
            targdxz = perc*maxcdxz
          endif 
          if (orb_dx2) then
            if (orb(10,j,ij) .gt. maxcdx2) maxcdx2 = orb(10,j,ij) 
            targdxy = perc*maxcdx2
          endif 
          if (norm(j,ij) .gt. maxnorm) maxnorm = norm(j,ij)
            targnorm = perc*maxnorm
        enddo 
      enddo

      if (orb_s)   print*, 'Maximum = ', maxcs,   targs
      if (orb_px)  print*, 'Maximum = ', maxcpx,  targpx
      if (orb_py)  print*, 'Maximum = ', maxcpy,  targpy
      if (orb_pz)  print*, 'Maximum = ', maxcpz,  targpz
      if (orb_dxy) print*, 'Maximum = ', maxcdxy, targdxy
      if (orb_dyz) print*, 'Maximum = ', maxcdyz, targdyz
      if (orb_dz2) print*, 'Maximum = ', maxcdz2, targdz2
      if (orb_dxz) print*, 'Maximum = ', maxcdxz, targdxz
      if (orb_dx2) print*, 'Maximum = ', maxcdx2, targdx2

      do ij = 1, nk
        do j = 36,36 ! nb
c         if (norm(j,ij) .ne. 0 ) then
          write(13,*) XNada(ij),YNada(ij),orbx(1,j,ij),
     :    orby(1,j,ij),a(j,ij)
c         endif
       enddo
      enddo

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      ENDDO ! SUD

 100  format(a15, i4, a20, i4, a19, i4)
 101  format(a17, f14.8, a19)
 102  format(a18, 3f11.8)
 103  format(f10.7, 4X, i4)
 104  format(a4, i4, a9, f14.8, a19)
 111  format(f10.7, 2x, f15.7, 2x, f10.7)

      end
