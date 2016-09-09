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



		call cycler_load(icycler)


        if (icycler.lt.0) then
			write(6,*) "cycler initialization failed!"
			stop
		end if



       call cycler_wsurf(icycler)

      call cycler_destroy(icycler)

      STOP
      
      END 






      SUBROUTINE WRITEGRD(NX,NY,NDX,Z,XMI,XMA,YMI,YMA,NAME, namef)
	CHARACTER*20 NAME,namef
	real*4 z(NDX,*),ZMI,ZMA


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





200   FORMAT('DSAA')
201   FORMAT(2I6)
202   FORMAT(10G13.6)
	RETURN
	END


      