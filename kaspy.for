      PROGRAM POM3



		icycler = -1

       call cycler_create(icycler)

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





      