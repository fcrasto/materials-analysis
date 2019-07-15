c************************************************************
c Last update 02/2015 Felipe Crasto de Lima
c************************************************************
	program band_vasp
	real, dimension(900000) :: x, Eneup, Enedo, x2, y2, z2, a
	real, dimension(900000) :: XNada, YNada, ZNada
	real  bx1,by1,bz1,bx2,by2,bz2,bx3,by3,bz3
	real  kp(3)
	real AA(3,3), XXNada, YYNada, ZZNada
	character(20) u,v,up,dow,ban
	real f
	integer b, nada, pk, nb, pn, soc
c*************************************************************
	u='EIGENVAL'
	v='inter'
	up='band-up.dat'
	 dow='band-down.dat'
	ban ='band.dat'
c*************************************************************
c read input
c*************************************************************
	open(42, file = 'inp.dat')
	read(42,*) soc
	read(42,*) spin
	read(42,*) f
	do m = 1, 3
	read(42,*) (AA(m,n), n = 1, 3)
	enddo
	close(42)
c*************************************************************
c arquivos de entrada
c*************************************************************
	open(10,file=u,status='old')
	open(11,file=v,status='unknown')
c*************************************************************
	do i=1,5,1
	read(10,*)
	end do
	read(10,*) nada, pk, nb
	pn=pk*nb
	read(10,*)
c*************************************************************
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                   SOC
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
	if(soc .eq. 1) then
	open(12,file=ban,status='unknown')
		do j=1,pk,1
		read(10,*) kp(1), kp(2), kp(3)
		XNada(j) = kp(1)*AA(1,1) + kp(2)*AA(2,1) + kp(3)*AA(3,1)
		YNada(j) = kp(1)*AA(1,2) + kp(2)*AA(2,2) + kp(3)*AA(3,2)
		ZNada(j) = kp(1)*AA(1,3) + kp(2)*AA(2,3) + kp(3)*AA(3,3)
		if (j .eq. 1) then
		x(j) = 0
		else
		XXNada = XNada(j-1) - XNada(j)
		YYNada = YNada(j-1) - YNada(j)
		ZZNada = ZNada(j-1) - ZNada(j)
		x(j) = sqrt(XXNada**2 + YYNada**2 + ZZNada**2)
		x(j) = x(j-1) + x(j)
		endif
		do i=1,nb,1
		read(10,*) a(i), Eneup(i)
		Eneup(i)=Eneup(i)-f
		write(11,*) x(j), Eneup(i) 
		end do
		if(j < pk)then
		read(10,*)
		end if
		end do
		close(10)
		close(11)
		open(11, file=v, status='old')
		do l=1,pn,1
		read(11,*) x2(l), y2(l)
		end do
		do k=1,nb,1
		do l=1,pk,1
		b=k+(l-1)*nb
		write(12,*) x2(b), y2(b)
		end do
		write(12,'(/)')
		end do
		close(10)
		close(11)
		close(12)
		stop
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                No-SOC
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
	elseif(soc .eq. 0) then
		if(spin .eq. 2) then
		open(12, file=up, status='unknown')
		open(13, file=dow, status='unknown')
			do j=1,pk,1
			read(10,*) kp(1), kp(2), kp(3)
			XNada(j) = kp(1)*AA(1,1) + kp(2)*AA(2,1) + kp(3)*AA(3,1)
			YNada(j) = kp(1)*AA(1,2) + kp(2)*AA(2,2) + kp(3)*AA(3,2)
			ZNada(j) = kp(1)*AA(1,3) + kp(2)*AA(2,3) + kp(3)*AA(3,3)
			if (j .eq. 1) then
			x(j) = 0
			else
			XXNada = XNada(j-1) - XNada(j)
			YYNada = YNada(j-1) - YNada(j)
			ZZNada = ZNada(j-1) - ZNada(j)
			x(j) = sqrt(XXNada**2 + YYNada**2 + ZZNada**2)
			x(j) = x(j-1) + x(j)
			endif
			do i=1,nb,1
			read(10,*) a(i), Eneup(i), Enedo(i)
			Eneup(i)=Eneup(i)-f
			Enedo(i)=Enedo(i)-f
			write(11,*) x(j), Eneup(i), Enedo(i)
			end do
			if(j < pk)then
			read(10,*)
			end if
			end do
			close(10)
			close(11)
			open(11, file=v, status='old')
			do l=1,pn,1
			read(11,*) x2(l), y2(l), z2(l)
			end do
			do k=1,nb,1
			do l=1,pk,1
			b=k+(l-1)*nb
			write(12,*) x2(b), y2(b)
			write(13,*) x2(b), z2(b)
			end do
			write(12,'(/)')
			write(13,'(/)') 
			end do
			close(10)
			close(11)
			close(12)
			close(13)
			stop
		elseif(spin .eq. 1) then
		open(12, file=ban, status='unknown')
			do j=1,pk,1
			read(10,*) kp(1), kp(2), kp(3)
			XNada(j) = kp(1)*AA(1,1) + kp(2)*AA(2,1) + kp(3)*AA(3,1)
			YNada(j) = kp(1)*AA(1,2) + kp(2)*AA(2,2) + kp(3)*AA(3,2)
			ZNada(j) = kp(1)*AA(1,3) + kp(2)*AA(2,3) + kp(3)*AA(3,3)
			if (j .eq. 1) then
			x(j) = 0
			else
			XXNada = XNada(j-1) - XNada(j)
			YYNada = YNada(j-1) - YNada(j)
			ZZNada = ZNada(j-1) - ZNada(j)
			x(j) = sqrt(XXNada**2 + YYNada**2 + ZZNada**2)
			x(j) = x(j-1) + x(j)
			endif
			do i=1,nb,1
			read(10,*) a(i), Eneup(i) 
			Eneup(i)=Eneup(i)-f
			write(11,*) x(j), Eneup(i) 
			 end do
			if(j < pk)then
			read(10,*)
			end if
			end do
			close(10)
			close(11)
			open(11, file=v, status='old')
			do l=1,pn,1
			read(11,*) x2(l), y2(l)
			end do
			do k=1,nb,1
			do l=1,pk,1
			b=k+(l-1)*nb
			write(12,*) x2(b), y2(b)
			end do
			write(12,'(/)')
			end do
			close(10)
			close(11)
			close(12)
			stop
		else
		write(*,*) 'Erro spin 1 ou 2'
		endif
	else
		write(*,*) 'Erro SOC 0 ou 1'
		stop
	endif
	end
