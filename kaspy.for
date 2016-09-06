c***********************************************************************
C     calculation of the wind response
C     result - 4 files: right, left, top and bottom interpolated 
c     velocities and pressure.
c     
C ---------------------------------------------------------
      PROGRAM POM3
C************************************************************************
C           TEST WITH SIMPLE FLAT BOTTOMED BASIN                        *
C                                                                       *
C THIS IS A VERSION OF THE THREE DIMENSIONAL, TIME DEPENDENT,           *
C PRIMITIVE EQUATION, OCEAN MODEL DEVELOPED BY ALLAN BLUMBERG AND ME    *
C WITH SUBSEQUENT CONTRIBUTIONS BY LEO OEY AND OTHERS. TWO REFERENCES   *
C ARE:                                                                  *
C                                                                       *
C     BLUMBERG, A.F. AND G.L. MELLOR, DIAGNOSTIC AND PROGNOSTIC         *
C     NUMERICAL CIRCULATION STUDIES OF THE SOUTH ATLANTIC BIGHT         *
C     J. GEOPHYS. RES. 88, 4579-4592, 1983.                             *
C                                                                       *
C     BLUMBERG, A.F. AND G.L. MELLOR, A DESCRIPTION OF A THREE          *
C     COASTAL OCEAN CIRCULATION MODEL, THREE DIMENSIONAL SHELF          *
C     MODELS, COASTAL AND ESTUARINE SCIENCES, 5, N.HEAPS, ED.,          *
C     AMERICAN GEOPHYSICAL UNION, 1987.                                 *
C                                                                       *
C IN SUBROUTINE PROFQ THE MODEL MAKES USE OF THE TURBULENCE CLOSURE     *
C SUB-MODEL MOST RECENTLY DESCRIBED IN:                                 *
C                                                                       *
C     MELLOR, G.L. AND T. YAMADA, DEVELOPMENT OF A TURBULENCE CLOSURE   *
C     MODEL FOR GEOPHYSICAL FLUID PROBLEMS, REV. GEOPHYS. SPACE PHYS.,  *
C     851-879, 1982.                                                    *
C                                                                       *
C                                          GEORGE MELLOR, 03/87         *
C                                                                       *
C IN THIS VERSION THE HORIZONTAL COORDINATES ARE GENERALIZED CURVI-     *
C LINEAR ORTHOGANAL COORDINATES AND THE MODEL HAS BEEN CONFIGURED       *
C TO RUN IN HALF PRECISION (32 BITS ) ON THE CYBER 205.                 *
C TO REDUCE ROUNDOFF ERROR, SALINTIES AND TEMPERATURES ARE REDUCED      *
C BY 35. AND 10. RESPECTIVELY. 64 BIT PRECISION IS USED LOCALLY         *
C IN SUBROUTINES DENS (WHEREIN THE T AND S OFFSET HAS BEEN TAKEN        *
C INTO ACCOUNT) AND PROFT.                                              *
C                                          GM, 09/87                    *
C                                                                       *
C THIS VERSION HAS BEEN CONVERTED TO STANDARD FORTRAN 77 BY STEVE       *
C BRENNER WITH SOME ADDITIONAL CHANGES BY ME IN PROFQ WHICH IS NOW      *
C SOMEWHAT SIMPLER. CODE RESTORED TO SINGLE PRECISION.                  *
C A USER'S GUIDE IS AVAILABLE:                                          *
C                                                                       *
C    MELLOR, G.L., A USER'S GUIDE FOR A THREE-DIMENSIONAL,              *
C    NUMERICAL OCEAN MODEL. PRINCETON UNIVERSITY REPORT, 1990.          *
C                                          GM, 04/90                    *
C                                                                       *
C Previous to this date the code was called pmod.f and all changes      *
C were recorded in pmod.change. Henceforth, the code will be called     * 
C pom.f and changes will be recorded in pom.change.                     *
C                                          GM, 11/92                    * 
C                                                                       *
C THE LAST CODE CHANGE, as recorded in pom.change, was on Feb. 7, 2001  *
C                                                                       *
C                                                                       *
C  Copyright 1996 The Trustees of Princeton University                  *
C                                                                       *
C This program is free software; you can redistribute it and/or         *
C modify it under the terms of the GNU General Public License           *
C as published by the Free Software Foundation, either version 2        *  
C of the License, or (at your option) any later version.                *
C                                                                       *
C This program is distributed in the hope that it will be useful,       *
C but WITHOUT ANY WARRANTY; without even the implied warranty of        *
C MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
C GNU General Public License for more details.                          *
C                                                                       *
C A copy of the GNU General Public License is available at              *
C        http://www.gnu.org/copyleft/gpl.html#SEC3                      *
C or by writing to The Free Software Foundation, Inc., 59 Temple Place  *
C Suite 330, Boston, MA  02111-1307, USA.                               *
C                                                                       *
C************************************************************************
C
      
c        !1 YEAR 1979
c        ! ALL TOGETHER 8760 hours


      INCLUDE 'comblk.for'    

 

      COMMON/F_FLOATS/ff_marker,
     3     fxf(im,jm),fyf(im,jm),fxb(im,jm),fyb(im,jm),
     4     FF(IM,JM),FB(IM,JM),
     7     ffu(im,jm),ffv(im,jm),fbu(im,jm),fbv(im,jm),
     8     ff_end

      COMMON/F_WIND/kx,ky,kt,kxu,kyu,ktu,kxv,kyv,ktv,
     1  XKI,XKA,YKI,YKA,TKI,TKA,XKUI,XKUA,YKUI,YKUA,TKUI,TKUA,
     2  XKVI,XKVA,YKVI,YKVA,TKVI,TKVA

      real*8
     5     SEL(IM,JM),SSEL(IM,JM),SFA(IM,JM),SSFA(IM,JM),SFEL(IM,JM),
     6     SFAR(IM,JM),SSFAR(IM,JM),SFELR(IM,JM),
     8     su(im,jm),sv(im,jm),ssu(im,jm),ssv(im,jm),
     9     ssuv(im,jm),ssue(im,jm),ssve(im,jm)
      

       real*8  ff_marker,ff_end



      DATA PI/3.141592654/,SMALL/1.E-10/

      DATA ri/0.01745329252/,GEE/9.807/
      CHARACTER*20 NAME,namep,nameu,namev,nameh
      integer ijloc(2)




       arrays_marker = 3.1415926535897932384626433832795010
       arrays_end_marker = 0.9876543211234567890

       ff_end = 100.001

		icycler = -1

       call cycler_create(icycler, dht, arrays_marker,
     1  ff_marker, kx)

		if (icycler.lt.0) then
			write(6,*) "cycler creation failed!"
			stop
		end if



C     READ IN GRID DATA AND INITIAL AND LATERAL BOUNDARY CONDITIONS
C------------------------------------------------------------------------
      nameh='kaspi_1_5mR.grd'//CHAR(0)

		CALL TIDEGEN_C(nameh)



C A very empirical specification of the bottom roughness
C          parameter follows
c      DO 45 J=1,JM
c      DO 45 I=1,IM
      
c      Z0B=.01
c      CBCMIN=.0025E0
c      IF (fsm(i,j).gt.0.5)THEN
c      CBC(I,J)=MAX(CBCMIN,.16E0/LOG(0.5*H(I,J)
c     1        /Z0B)**2)
c        ELSE
c        CBC(I,J)=0
c        END IF
c45    CONTINUE
C   Evaluate external CFL time step
c      DO 81 J=1,JM
c      DO 81 I=1,IM
c      D(I,J)=H(I,J)+EL(I,J)
c   81 CONTINUE
   
      FXF=0
      FYF=0
      FXB=0
      FYB=0
      FF=100
      FB=100
      FA=0
      SEL=0
      SSEL=0
      SFA=0
      SSFA=0
      SFEL=0
      SFAR=0
      SSFAR=0
      SFELR=0
      
      su=0
      sv=0
      ssu=0
      ssv=0
      ssuv=0
      sue=0
      sve=0



		call cycler_load(icycler)


        if (icycler.lt.0) then
			write(6,*) "cycler initialization failed!"
			stop
		end if



c        iold=0
c
c		DO 9000 IINT=1,IEND
c
c      timeh=(iint-1)*dti/3600.0
c      ihour_s=600/dti
c
c	iwrite=mod(iint,ihour_s)
c
c      if(iwrite.eq.1) then

c		 call cycler_find_elves(icycler)
c		write(6,1117) 't=',timeh,'h','Sea level=',elfmax,elfmin,'m'
c
c
c1117   FORMAT(a3,f12.4,1x,a1,5x,a10,f8.4,1x,f8.4,1x,a1,a1)
c c     end if
c*************************************************

c      itimeh=int(timeh)



c     	if(itimeh.gt.iold) then  !! writing to file and compute statiatics,
c STATISTICS WAS HERE
c		call cycler_get_data_back(icycler)
c		write(77,'(101f10.3)') timeh,(el(stx(kk),sty(kk)),kk=1,nstation)
c	   end if


       call cycler_wsurf(icycler)


c       iold=itimeh


C
C--------------------------------------------------------------------
C           BEGIN PRINT SECTION
C--------------------------------------------------------------------
c      ktime=timeh
c      if (ktime.gt.ktime0) then
c      print velocities
c        i0=191
c	  j0=71	 
c        write(9,*)
c        write(9,950) time
c	  do k=1,kb-1
c	    write(9,952) ((uf(i,j,k),i=2,imm1),j=2,jmm1)  
c        end do
c	  write(9,*)
c	  do k=1,kb-1
c	    write(9,952) ((vf(i,j,k),i=2,imm1),j=2,jmm1)  
c        end do
c        write(10,956)time,(uf(i0,j0,k),vf(i0,j0,k),k=1,kb-1)
c	  ktime0=ktime
c      end if
 950  format(f9.5,' Time')
 951  format('Left')
 953  format('Right')
 954  format('Bottom')
 955  format('Top')
 956  format(1x,13f8.4)

 952  format(1x,10f8.4)



c c7000 CONTINUE
C--------------------------------------------------------------------
C             END PRINT SECTION
C--------------------------------------------------------------------
c c9000                     CONTINUE


c       call cycler_get_data_back(icycler)

         
C***********************************************************************
C                                                                      *
C                  END NUMERICAL INTEGRATION                           *
C                                                                      *
C***********************************************************************
c      close(2)
c	close(3)



c      pause ' end'
844   format(80f7.2)

c      do k=1,kb-1
c      do j=2,jmm1
c		WRITE(2,946) (100*uf(i,j,k),i=2,IMm1)
c      end do
c	end do
c      do k=1,kb-1
c	do j=2,jmm1	  
c	    write(3,946) (100*vf(i,j,k),i=2,imm1)
c      end do
c	end do

      do j=2,jmm1
		WRITE(22,946) (100*EL(i,j),i=2,IMm1)
      end do
c	do j=2,jmm1	  
c	    write(23,946) (100*vaf(i,j),i=2,imm1)
c      end do
C     STATISTICS
      EPS=1.0E-10

c      NAME='SSFA.GRD'
c      CALL WRITEGRD(IM,JM,IM,SSFA,10.0,31.0,53.0,66.0,NAME)
cs       NAME='SSEL.GRD'
cs      CALL WRITEGRD(IM,JM,IM,SSEL,xmi,xma,ymi,yma,NAME)
c      NAME='SFEL.GRD'
c      CALL WRITEGRD(IM,JM,IM,SFEL,10.0,31.0,53.0,66.0,NAME)

c      NAME='ALFA1.GRD'
c      CALL WRITEGRD(IM,JM,IM,SFA,10.0,31.0,53.0,66.0,NAME)
c      NAME='ALFA2.GRD'
c      CALL WRITEGRD(IM,JM,IM,SEL,10.0,31.0,53.0,66.0,NAME)
      
cs      NAME='SSFAR.GRD'
cs      CALL WRITEGRD(IM,JM,IM,SSFAR,xmi,xma,ymi,yma,NAME)
cs      NAME='SFELR.GRD'
cs      CALL WRITEGRD(IM,JM,IM,SFELR,xmi,xma,ymi,yma,NAME)
c      NAME='ALFAR.GRD'
c      CALL WRITEGRD(IM,JM,IM,SFAR,10.0,31.0,53.0,66.0,NAME)
cs      NAME='ssu.GRD'
cs      CALL WRITEGRD(IM,JM,IM,ssu,xmi,xma,ymi,yma,NAME)
cs      NAME='ssv.GRD'
cs      CALL WRITEGRD(IM,JM,IM,ssv,xmi,xma,ymi,yma,NAME)
cs      NAME='ssuv.GRD'
cs      CALL WRITEGRD(IM,JM,IM,ssuv,xmi,xma,ymi,yma,NAME)
cs      NAME='ssue.GRD'
cs      CALL WRITEGRD(IM,JM,IM,ssue,xmi,xma,ymi,yma,NAME)
cs      NAME='ssve.GRD'
cs      CALL WRITEGRD(IM,JM,IM,ssve,xmi,xma,ymi,yma,NAME)
      
944   format('DSAA')
945   format(2i6)
946   format(10f7.2)
c      end if
c      end do
c	end do
	 
      close(2)
      close(3) 
        
      close(22)
      close(23) 

      call cycler_destroy(icycler)

      STOP
      
      END 







      SUBROUTINE WRITEGRD(NX,NY,NDX,Z,XMI,XMA,YMI,YMA,NAME)
	CHARACTER*20 NAME
	real*8 z(NDX,*),ZMI,ZMA
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


      