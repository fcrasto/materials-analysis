cccccccccccccccccccccccccccccccccccc
c   Programa para determinar SQS   c
c   em liga A_{x}B_{1-x}  ; x=a/b  c
c       Last update 02/2018        c
c      Felipe Crasto de Lima       c
cccccccccccccccccccccccccccccccccccc
	program SQS
	real, dimension(500000) :: x1, x2, x3
	real, dimension(500000) :: tipo
	real :: a11, a12, a13
	real :: a21, a22, a23
	real :: a31, a32, a33
	real :: cor1, cor2, cor3, cor 
	real :: conv1, conv2, conv3, conv
	real :: a, b  ! proporcao x= a/b
	real :: dis, dis1, dis2, dis3
	real :: ct1, ct2, ct3
	integer n0atm, ncell, na, nb
	integer nx, ny, nz    ! tamanho da célula SQS
	integer mx
	real :: nbond1, nbond2, nbond3
	character(20) w,v,p
c*************************************************************
c		Entrada
c*************************************************************
	w='POSCAR'
	v='POSCAR.out'
	p='test'
	conv = 1.
	conv1 = 1.
	conv2 = 1.
	conv3 = 1.
	write(*,*) "cccccccccccccccccccccccccccccccccccc"
	write(*,*) "c   Programa para determinar SQS   c"
	write(*,*) "c   em liga A_{x}B_{1-x}  ; x=a/b  c"
	write(*,*) "cccccccccccccccccccccccccccccccccccc"
	write(*,*) "Entre com a multiplicidade da célula SQS nxm: n m"
	 read(*,*) nx, ny, nz
	write(*,*) "Proporção da liga x=a/b: a b"
	read(*,*) a, b
	write(*,*) "Convergência abs(cor-cor_{inf})<=c: c1 c2 c3"
 	read(*,*) ct1, ct2, ct3
	
c	DO WHILE((conv1.gt.ct1).OR.(conv2.gt.ct2).OR.(conv3.gt.ct3))
 	open(12,file=w,status='old')
c leitura de dados
c	write(*,*) "Lendo dados das coordenadas..."
	read(12,*)
	read(12,*)
	read(12,*) a11, a12, a13
	read(12,*) a21, a22, a23
	read(12,*) a31, a32, a33
	read(12,*)
	read(12,*) n0atm
	read(12,*)
	do n=1,n0atm
		read(12,*) x1(n), x2(n), x3(n)
	enddo
	close(12)
c Passando Para coordenadas cartesianas
	do n=1,n0atm
		x1(n) = x1(n)*a11 + x2(n)*a21 + x3(n)*a31
		x2(n) = x1(n)*a12 + x2(n)*a22 + x3(n)*a32
		x3(n) = x1(n)*a13 + x2(n)*a23 + x3(n)*a33
c		write(*,*) x1(n), x2(n), x3(n)
	enddo
c contruindo célula SQS
	ncell = n0atm*nx*ny*nz
	mx = 0
	do n=1,nx
		do m=1,ny
			do l=1,n0atm
				x1(l+mx*n0atm) = x1(l)+((n-1)*a11 + (m-1)*a21)
				x2(l+mx*n0atm) = x2(l)+((n-1)*a12 + (m-1)*a22)
				x3(l+mx*n0atm) = x3(l)+((n-1)*a13 + (m-1)*a23)
			enddo
		mx = mx+1
		enddo
	enddo
	DO WHILE((conv1.gt.ct1).OR.(conv2.gt.ct2).OR.(conv3.gt.ct3))
c Sorteando cada sítio de acordo com a proporção x=a/b
c	write(*,*) "Sorteio dos Sítios..."
	do n=1,ncell
		call random_number(u)
		j = 1 + floor(b*u)
  		if(j .le. a)then
			tipo(n) = 1
		else
			tipo(n) = -1
		endif
c	write(*,*) x1(n), x2(n), x3(n), tipo(n)
	enddo
c contruindo imagens periódicas 
c	write(*,*) "Construindo imagens periódicas..."
	do l=1,ncell
		x1(l+1*ncell)= x1(l)+(nx*a11)
		x1(l+2*ncell)= x1(l)-(nx*a11)
		x1(l+3*ncell)= x1(l)+(ny*a21)
		x1(l+4*ncell)= x1(l)-(ny*a21)
		x1(l+5*ncell)= x1(l)+(nz*n31)
		x1(l+6*ncell)= x1(l)-(nz*n31) 
		x1(l+7*ncell)= x1(l)+(nx*a11 + ny*a21)
		x1(l+8*ncell)= x1(l)-(nx*a11 + ny*a21)
		x1(l+9*ncell)= x1(l)+(nx*a11 - ny*a21)
		x1(l+10*ncell)= x1(l)-(nx*a11 - ny*a21)

		x1(l+11*ncell)= x1(l)+(nx*a11 + nz*a31)
		x1(l+12*ncell)= x1(l)-(nx*a11 + nz*a31)
		x1(l+13*ncell)= x1(l)+(ny*a21 + nz*a31)
		x1(l+14*ncell)= x1(l)-(ny*a21 + nz*a31)
		
		x1(l+15*ncell)= x1(l)+(nx*a11 - nz*a31)
		x1(l+16*ncell)= x1(l)-(nx*a11 - nz*a31)
		x1(l+17*ncell)= x1(l)+(ny*a21 - nz*a31)
		x1(l+18*ncell)= x1(l)-(ny*a21 - nz*a31)

		x1(l+19*ncell)= x1(l)+(nx*a11 + ny*a21 + nz*a31)
		x1(l+20*ncell)= x1(l)-(nx*a11 + ny*a21 + nz*a31)
		x1(l+21*ncell)= x1(l)+(nx*a11 - ny*a21 + nz*a31)
		x1(l+22*ncell)= x1(l)-(nx*a11 - ny*a21 + nz*a31)
		
		x1(l+23*ncell)= x1(l)+(nx*a11 + ny*a21 - nz*a31)
		x1(l+24*ncell)= x1(l)-(nx*a11 + ny*a21 - nz*a31)
		x1(l+25*ncell)= x1(l)+(nx*a11 - ny*a21 - nz*a31)
		x1(l+26*ncell)= x1(l)-(nx*a11 - ny*a21 - nz*a31)
cccccc		
		x2(l+1*ncell)= x2(l)+(nx*a12)
		x2(l+2*ncell)= x2(l)-(nx*a12)
		x2(l+3*ncell)= x2(l)+(ny*a22)
		x2(l+4*ncell)= x2(l)-(ny*a22)
		x2(l+5*ncell)= x2(l)+(nz*a32)
		x2(l+6*ncell)= x2(l)-(nz*a32)
		x2(l+7*ncell)= x2(l)+(nx*a12 + ny*a22)
		x2(l+8*ncell)= x2(l)-(nx*a12 + ny*a22)
		x2(l+9*ncell)= x2(l)+(nx*a12 - ny*a22)
		x2(l+10*ncell)= x2(l)-(nx*a12 - ny*a22)

		x2(l+11*ncell)= x2(l)+(nx*a12 + nz*a32)
		x2(l+12*ncell)= x2(l)-(nx*a12 + nz*a32)
		x2(l+13*ncell)= x2(l)+(ny*a22 + nz*a32)
		x2(l+14*ncell)= x2(l)-(ny*a22 + nz*a32)
		
		x2(l+15*ncell)= x2(l)+(nx*a12 - nz*a32)
		x2(l+16*ncell)= x2(l)-(nx*a12 - nz*a32)
		x2(l+17*ncell)= x2(l)+(ny*a22 - nz*a32)
		x2(l+18*ncell)= x2(l)-(ny*a22 - nz*a32)

		x2(l+19*ncell)= x2(l)+(nx*a12 + ny*a22 + nz*a32)
		x2(l+20*ncell)= x2(l)-(nx*a12 + ny*a22 + nz*a32)
		x2(l+21*ncell)= x2(l)+(nx*a12 - ny*a22 + nz*a32)
		x2(l+22*ncell)= x2(l)-(nx*a12 - ny*a22 + nz*a32)
		
		x2(l+23*ncell)= x2(l)+(nx*a12 + ny*a22 - nz*a32)
		x2(l+24*ncell)= x2(l)-(nx*a12 + ny*a22 - nz*a32)
		x2(l+25*ncell)= x2(l)+(nx*a12 - ny*a22 - nz*a32)
		x2(l+26*ncell)= x2(l)-(nx*a12 - ny*a22 - nz*a32)
cccccc	
		x3(l+1*ncell)= x3(l)+(nx*a13)
		x3(l+2*ncell)= x3(l)-(nx*a13)
		x3(l+3*ncell)= x3(l)+(ny*a23)
		x3(l+4*ncell)= x3(l)-(ny*a23)
		x3(l+5*ncell)= x3(l)+(nz*a33)
		x3(l+6*ncell)= x3(l)-(nz*a33)
		x3(l+7*ncell)= x3(l)+(nx*a13 + ny*a23)
		x3(l+8*ncell)= x3(l)-(nx*a13 + ny*a23)
		x3(l+9*ncell)= x3(l)+(nx*a13 - ny*a23)
		x3(l+10*ncell)= x3(l)-(nx*a13 - ny*a23)

		x3(l+11*ncell)= x3(l)+(nx*a13 + nz*a33)
		x3(l+12*ncell)= x3(l)-(nx*a13 + nz*a33)
		x3(l+13*ncell)= x3(l)+(ny*a23 + nz*a33)
		x3(l+14*ncell)= x3(l)-(ny*a23 + nz*a33)
		
		x3(l+15*ncell)= x3(l)+(nx*a13 - nz*a33)
		x3(l+16*ncell)= x3(l)-(nx*a13 - nz*a33)
		x3(l+17*ncell)= x3(l)+(ny*a23 - nz*a33)
		x3(l+18*ncell)= x3(l)-(ny*a23 - nz*a33)

		x3(l+19*ncell)= x3(l)+(nx*a13 + ny*a23 + nz*a33)
		x3(l+20*ncell)= x3(l)-(nx*a13 + ny*a23 + nz*a33)
		x3(l+21*ncell)= x3(l)+(nx*a13 - ny*a23 + nz*a33)
		x3(l+22*ncell)= x3(l)-(nx*a13 - ny*a23 + nz*a33)
		
		x3(l+23*ncell)= x3(l)+(nx*a13 + ny*a23 - nz*a33)
		x3(l+24*ncell)= x3(l)-(nx*a13 + ny*a23 - nz*a33)
		x3(l+25*ncell)= x3(l)+(nx*a13 - ny*a23 - nz*a33)
		x3(l+26*ncell)= x3(l)-(nx*a13 - ny*a23 - nz*a33)
cccccc
		tipo(l+1*ncell) = tipo(l)
		tipo(l+2*ncell) = tipo(l)
		tipo(l+3*ncell) = tipo(l)
		tipo(l+4*ncell) = tipo(l)
		tipo(l+5*ncell) = tipo(l)
		tipo(l+6*ncell) = tipo(l)
		tipo(l+7*ncell) = tipo(l)
		tipo(l+8*ncell) = tipo(l)
		
		tipo(l+9*ncell) = tipo(l)
		tipo(l+10*ncell) = tipo(l)
		tipo(l+11*ncell) = tipo(l)
		tipo(l+12*ncell) = tipo(l)
		tipo(l+13*ncell) = tipo(l)
		tipo(l+14*ncell) = tipo(l)
		tipo(l+15*ncell) = tipo(l)
		tipo(l+16*ncell) = tipo(l)
		tipo(l+17*ncell) = tipo(l)
		tipo(l+18*ncell) = tipo(l)
		tipo(l+19*ncell) = tipo(l)
		tipo(l+20*ncell) = tipo(l)
		tipo(l+21*ncell) = tipo(l)
		tipo(l+22*ncell) = tipo(l)
		tipo(l+23*ncell) = tipo(l)
		tipo(l+24*ncell) = tipo(l)
		tipo(l+25*ncell) = tipo(l)
		tipo(l+26*ncell) = tipo(l)
	enddo
c Determinar distancia entre 1º, 2º, 3º vizinhos
c	write(*,*) "Determinando distancia entre vizinhos..."
	dis = 0.
	dis1 = 1000.
	dis2 = 1000.
	dis3 = 1000.
	do l=1,3
		do n=1,ncell
			do m=n+1,27*ncell
   		 dis = sqrt((x1(n)-x1(m))**2+(x2(n)-x2(m))**2+(x3(n)-x3(m))**2)
			 if(dis .le. dis1+0.1)then
				  dis1 = dis
			elseif(dis .le. dis2+0.1)then
				  dis2 = dis
			elseif(dis .le. dis3+0.1)then
				  dis3 = dis
			endif
			enddo
		enddo
	enddo
c	write(*,*) "Distancia 1º vizinhos:", dis1
c	write(*,*) "Distancia 2º vizinhos:", dis2
c	write(*,*) "Distancia 3º vizinhos:", dis3
c    Calculo da correlação
 	dis = 0.
	cor = 0.
	cor1 = 0.
	cor2 = 0.
	cor3 = 0.
	nbond1 = 0.
	nbond2 = 0.
	nbond3 = 0.
c	dis1 = dis1+0.1
c	dis2 = dis2+0.1
c	dis3 = dis3+0.1
c	write(*,*) "Calculo da corelação..."
	do n=1,ncell
		do m=1,27*ncell 
			 dis = sqrt((x1(n)-x1(m))**2+(x2(n)-x2(m))**2+(x3(n)-x3(m))**2)
			if(m .ne. n)then
 			if(dis .le. dis1+0.1)then
  				nbond1 = nbond1 + 1
  				cor1 = cor1 + tipo(n)*tipo(m) 
			elseif(dis .le. dis2+0.1)then
  				nbond2 = nbond2 + 1
  				cor2 = cor2 + tipo(n)*tipo(m)
			elseif(dis .le. dis3+0.1)then
  				nbond3 = nbond3 + 1
  				cor3 = cor3 + tipo(n)*tipo(m)
			endif
			endif
		enddo
	enddo
 	cor = (cor1/nbond1) + (cor2/nbond2) + (cor3/nbond3)
c	write(*,*) "correlation m=1, = ", cor1, nbond1
c	write(*,*) "correlation m=2, = ", cor2, nbond2
c	write(*,*) "correlation m=3, = ", cor3, nbond3
c	write(*,*) "total correlation", cor
	conv1=abs((cor1/nbond1) - ((2*(a/b) - 1)**2))
	conv2=abs((cor2/nbond2) - ((2*(a/b) - 1)**2))
	conv3=abs((cor3/nbond3) - ((2*(a/b) - 1)**2))
	write(*,*) "convergence parameter m=1: ", conv1
	write(*,*) "convergence parameter m=2: ", conv2
	write(*,*) "convergence parameter m=3: ", conv3
c	conv = conv1
c	if(conv2 .gt. conv)then
c	conv = conv2
c	endif
c	if(conv3 .gt. conv)then
c	conv = conv3
c	endif
	
	ENDDO
	write(*,*) "convergencia alcançada.", conv1, conv2, conv3
	write(*,*) "correlation m=1, = ", cor1, nbond1
	write(*,*) "correlation m=2, = ", cor2, nbond2
	write(*,*) "correlation m=3, = ", cor3, nbond3
	write(*,*) "total correlation", cor
c	write(*,*) "dis 1 2 3", dis1, dis2, dis3
c Writing SQS cell
	na=0
	nb=0
	do n=1,ncell
	if(tipo(n) .gt. 0)then
	na = na+1
	else
	nb= nb+1
	endif
	enddo
	 open(13,file=v,status='unknown')
	write(13,*) "A B   !", conv1, conv2, conv3
	write(13,*) "1.000000000"
	write(13,*) nx*a11, nx*a12, nx*a13
	write(13,*) ny*a21, ny*a22, ny*a23
	write(13,*) a31, a32, a33
	write(13,*) "A   B"
	write(13,*) na, nb
	write(13,*) "Cartesian"
	do n=1,ncell
		if(tipo(n) .eq. 1.0)then
		write(13,*) x1(n), x2(n), x3(n), tipo(n)
		endif
	enddo
	do n=1,ncell
		if(tipo(n) .eq. -1.0)then
		write(13,*) x1(n), x2(n), x3(n), tipo(n)
		endif
	enddo
	write(*,*) "Distancia 1º vizinhos:", dis1
	write(*,*) "Distancia 2º vizinhos:", dis2
	write(*,*) "Distancia 3º vizinhos:", dis3
	close(13)
	stop
	end
