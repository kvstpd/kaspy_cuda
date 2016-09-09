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








		icycler = -1

       call cycler_create(icycler, dht, arrays_marker, kx)

		if (icycler.lt.0) then
			write(6,*) "cycler creation failed!"
			stop
		end if



		call cycler_load(icycler)


        if (icycler.lt.0) then
			write(6,*) "cycler initialization failed!"
			stop
		end if



       call cycler_wsurf(icycler)

      call cycler_destroy(icycler)

      STOP
      
      END 





      