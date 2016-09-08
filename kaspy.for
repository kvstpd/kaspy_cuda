c***********************************************************************
C     calculation of the wind response
C     result - 4 files: right, left, top and bottom interpolated 
c     velocities and pressure.
c     
C ---------------------------------------------------------
      PROGRAM POM3

      INCLUDE 'comblk.for'    

 


      COMMON/F_WIND/kx,ky,kt,kxu,kyu,ktu,kxv,kyv,ktv,
     1  XKI,XKA,YKI,YKA,TKI,TKA,XKUI,XKUA,YKUI,YKUA,TKUI,TKUA,
     2  XKVI,XKVA,YKVI,YKVA,TKVI,TKVA


c      COMMON/F_STATS/SEL,SSEL,SFA,SSFA,SFEL,SFAR,SSFAR,SFELR,
c     8     su,sv,ssu,ssv,ssuv,ssue,ssve

c      real*4
c     5     SEL(IM,JM),SSEL(IM,JM),SFA(IM,JM),SSFA(IM,JM),SFEL(IM,JM),
c     6     SFAR(IM,JM),SSFAR(IM,JM),SFELR(IM,JM),
c     8     su(im,jm),sv(im,jm),ssu(im,jm),ssv(im,jm),
c     9     ssuv(im,jm),ssue(im,jm),ssve(im,jm)
      


      CHARACTER*20 NAME,nameh




       arrays_marker = 3.1415926535897932384626433832795010
       arrays_end_marker = 0.9876543211234567890


		icycler = -1

       call cycler_create(icycler, dht, arrays_marker, kx)

		if (icycler.lt.0) then
			write(6,*) "cycler creation failed!"
			stop
		end if


      nameh='kaspi_1_5mR.grd'//CHAR(0)

		CALL TIDEGEN_C(nameh)



c c     SEL=0
c      SSEL=0
c      SFA=0
c      SSFA=0
c      SFEL=0
c      SFAR=0
c      SSFAR=0
c      SFELR=0
      
c      su=0
c      sv=0
c      ssu=0
c      ssv=0
c      ssuv=0
c      sue=0
c      sve=0



		call cycler_load(icycler)


        if (icycler.lt.0) then
			write(6,*) "cycler initialization failed!"
			stop
		end if



       call cycler_wsurf(icycler)

c		call cycler_get_data_back(icycler)


      call cycler_destroy(icycler)

      STOP
      
      END 



      SUBROUTINE READGRD(NX,NY,NDX,Z,NAME)
	CHARACTER*20 NAME
	real z(NDX,*)
	OPEN(1,FILE=NAME)
	READ(1,*)
	READ(1,*)NX,NY
	READ(1,*)XMI,XMA
	READ(1,*)YMI,YMA
	READ(1,*)
	READ(1,8202)((z(i,j),i=1,nx),j=1,ny)
	CLOSE(1)
8202   FORMAT(10G13.6)
	RETURN
	END



      SUBROUTINE WRITEGRD(NX,NY,NDX,Z,XMI,XMA,YMI,YMA,NAME, namef)
	CHARACTER*20 NAME,namef
	real*4 z(NDX,*),ZMI,ZMA

	real*4 fdata(322,442)




c	write(6,*) "writing to file: ", NAME

	OPEN(1,FILE=NAME)
	WRITE(1,200)
	WRITE(1,201)NX,NY
	WRITE(1,202)XMI,XMA
	WRITE(1,202)YMI,YMA
	ZMAX=Z(1,1)
	ZMIN=Z(1,1)
	DO J=1,NY
	DO I=1,NX
	  ZMAX=MAX(ZMAX,Z(I,J))
	  ZMIN=MIN(ZMIN,Z(I,J))
	END DO
	END DO  
	  
	WRITE(1,202)ZMIN,ZMAX
	DO J=1,NY
	  WRITE(1,202)(z(i,j),i=1,nx)
	END DO 
	CLOSE(1)


c	namef = 'ssel_f.grd'

	call READGRD(NX,NY,NDX,fdata,NAMEF)

	dmax = 0.0
	dcurr = 0.0

	DO J=1,NY
	DO I=1,NX
	  dcurr = abs(Z(I,J) - fdata(i,j))

		dmax=MAX(dmax,dcurr)
	END DO
	END DO

	write(6,*) "namef is: ", NAMEF
	write(6,*) "dmax is: ", dmax


200   FORMAT('DSAA')
201   FORMAT(2I6)
202   FORMAT(10G13.6)
	RETURN
	END


      