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
      LOGICAL SEAMT
 

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
      
      REAL, ALLOCATABLE :: PRESS(:,:,:),uwd(:,:,:),vwd(:,:,:)
c      REAL, ALLOCATABLE :: PRESS0(:,:),uwd0(:,:),vwd0(:,:)

       real*8  ff_marker,ff_end





C      PUBLIC PRESS
C      COMMON/wind/uwd,vwd
     

C--------------------------------------------------------------------
      REAL ISPI,ISP2I
      DATA PI/3.141592654/,SMALL/1.E-10/
      character crr
      integer icr,stx(100),sty(100)
      equivalence (icr,crr)
      data icr/13/
      DATA ri/0.01745329252/,GEE/9.807/
      CHARACTER*20 NAME,namep,nameu,namev,nameh
      integer ijloc(2)


	  call cycler_read_ini()




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


c 		iasize = kx*ky*kt

		

c		save_z_(int * nx,int * nsize,float * z, char * name)


C
C     READ IN GRID DATA AND INITIAL AND LATERAL BOUNDARY CONDITIONS
C------------------------------------------------------------------------
      nameh='kaspi_1_5mR.grd'//CHAR(0)
      CALL TIDEGEN(nameh)
C
C********************************************************************
C  Note that lateral thermodynamic boundary conditions are often set equal
C  to the initial conditions and are held constant thereafter. Users can
C  of course create varable boundary conditions.
C********************************************************************

C     END READING DATA      
      
	open(111,file='s66.txt'//CHAR(0))
	do i=1,100
	read(111,*,end=167) stx(i),sty(i)
	end do
167	nstation=i-1


        hours=(NH6-1)*dht
c	hours=720
        write (*,*) 'duration=',hours

c     read wind file
	 
      GRAV=9.806
      TBIAS=0.
      SBIAS=0.
      IINT=0
      TIMEh=0.
      ittt=0
C
      SEAMT=.FALSE.

c      ampl=0.0
c	awind=-22.3/3.6 
      time0=1.0
c
	     
      DAYS=5. 
      PRTD1=5.   
      ISPADV=10
      SMOTH=0.10
      HORCON=1.0
C  N.B.: TPRNI=0. yields zero  horizontal diffusivity !!
      TPRNI=0   
      UMOL=2.E-5
      MODE=2
      NREAD=0
C
      ISWTCH=100000
      IPRTD2=4
ccc	ISPLIT=19
c	DTE=21.0526315789

	ISPLIT=1
c	DTE=20.0

      DTI=DTE*ISPLIT 
      DTI=DTE*FLOAT(ISPLIT)
      DTE2=DTE*2
      DTI2=DTI*2

      IEND=HOURS*3600/DTI+2

      write(6,*)
      write(6,*)
      write(6,*)
      write(6,*)
      write(6,*) 'Duration (hours)',hours
c******
c      Read(*,*) hours
	


	tprni=1.0
      IEND=HOURS*3600/DTI+2
      IPRINT=PRTD1*3600/100
      ISWTCH=ISWTCH*3600
C
      ISPI=1.E0/FLOAT(ISPLIT)
      ISP2I=1.E0/(2.E0*FLOAT(ISPLIT))
C----------------------------------------------------------------------
C             ESTABLISH PROBLEM CHARACTERISTICS
C          ****** ALL UNITS IN M.K.S. SYSTEM ******
C      F,BLANK AND B REFERS TO FORWARD,CENTRAL AND BACKWARD TIME LEVELS.
C----------------------------------------------------------------------
      IP=IM
      ISK=2
      JSK=2

C
C  Additional initial conditions

      itime6_old=0


      RAMP=1.
C A very empirical specification of the bottom roughness
C          parameter follows
      DO 45 J=1,JM
      DO 45 I=1,IM
      
      Z0B=.01
      CBCMIN=.0025E0
      IF (fsm(i,j).gt.0.5)THEN
      CBC(I,J)=MAX(CBCMIN,.16E0/LOG(0.5*H(I,J)
     1        /Z0B)**2)
        ELSE
        CBC(I,J)=0
        END IF
45    CONTINUE
C   Evaluate external CFL time step
      DO 81 J=1,JM
      DO 81 I=1,IM
      D(I,J)=H(I,J)+EL(I,J)
   81 CONTINUE
   
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
      
c      ro_ratio=1.29/1020.0
      
      NSTAT=0
      BTIM=1
      FTIM=0
      
      AFSM=0
      DO J=1,JM
      DO I=1,IM
        AFSM=AFSM+FSM(I,J)
      END DO
      END DO 
      
      open(77,file='wnd_press_90d.txt')
c!      open(88,file='force1.txt') ALSO STATISTICS

C
c      TIME=TIMEI
C***********************************************************************
C                                                                      *
C                  BEGIN NUMERICAL INTEGRATION                         *
C                                                                      *
C***********************************************************************
C

		call cycler_load(icycler)


c		write(6,*) icycler

        if (icycler.lt.0) then
			write(6,*) "cycler initialization failed!"
			stop
		end if



c call display(surf,im,jm,im,jm,-1.0,1.0,0)
        iold=0

                     DO 9000 IINT=1,IEND

      timeh=(iint-1)*dti/3600.0
      ihour_s=600/dti

	iwrite=mod(iint,ihour_s)

      if(iwrite.eq.1) then

		 call cycler_find_elves(icycler)
		write(6,1117) 't=',timeh,'h','Sea level=',elfmax,elfmin,'m'


1117   FORMAT(a3,f12.4,1x,a1,5x,a10,f8.4,1x,f8.4,1x,a1,a1)
      end if
c*************************************************

      itimeh=int(timeh)



     	if(itimeh.gt.iold) then  !! writing to file and compute statiatics,
c STATISTICS WAS HERE
		call cycler_get_data_back(icycler)
		write(77,'(101f10.3)') timeh,(el(stx(kk),sty(kk)),kk=1,nstation)
	   end if


       call cycler_wsurf(icycler)


       iold=itimeh


C
C--------------------------------------------------------------------
C           BEGIN PRINT SECTION
C--------------------------------------------------------------------
      ktime=timeh 
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



 7000 CONTINUE
C--------------------------------------------------------------------
C             END PRINT SECTION
C--------------------------------------------------------------------
 9000                     CONTINUE


       call cycler_get_data_back(icycler)

         
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



C               
C                  
      SUBROUTINE SLPMIN(H,IM,JM,FSM,SL)
      DIMENSION H(IM,JM),FSM(IM,JM)
      DIMENSION SL(IM,JM)
C     This subroutine limits the maximum difference in H divided
C     by twice the average of H of adjacent cells.
C     The maximum possible value is unity.
C     
      SLMIN=0.2
C
      DO 100 LOOP=1,10
C    sweep right
      DO 3 J=2,JM-1
      DO 1 I=2,IM-1
      IF(FSM(I,J).EQ.0..OR.FSM(I+1,J).EQ.0.) GOTO 1
      SL(I,J)=ABS(H(I+1,J)-H(I,J))/(H(I,J)+H(I+1,J))
      IF(SL(I,J).LT.SLMIN) GOTO 1
      DH=0.5*(SL(I,J)-SLMIN)*(H(I,J)+H(I+1,J))
      SN=-1.
      IF(H(I+1,J).GT.H(I,J)) SN=1.
      H(I+1,J)=H(I+1,J)-SN*DH
      H(I,J)=H(I,J)+SN*DH
   1  CONTINUE
C    sweep left
      DO 2 I=IM-1,2,-1
      IF(FSM(I,J).EQ.0..OR.FSM(I+1,J).EQ.0.) GOTO 2
      SL(I,J)=ABS(H(I+1,J)-H(I,J))/(H(I,J)+H(I+1,J))
      IF(SL(I,J).LT.SLMIN) GOTO 2
      DH=0.5*(SL(I,J)-SLMIN)*(H(I,J)+H(I+1,J))
      SN=-1.
      IF(H(I+1,J).GT.H(I,J)) SN=1.
      H(I+1,J)=H(I+1,J)-SN*DH
      H(I,J)=H(I,J)+SN*DH
   2  CONTINUE
   3  CONTINUE
C   sweep up     
      DO 13 I=2,IM-1
      DO 11 J=2,JM-1
      IF(FSM(I,J).EQ.0..OR.FSM(I,J+1).EQ.0.) GOTO 11
      SL(I,J)=ABS(H(I,J+1)-H(I,J))/(H(I,J)+H(I,J+1))
      IF(SL(I,J).LT.SLMIN) GOTO 11
      DH=0.5*(SL(I,J)-SLMIN)*(H(I,J)+H(I,J+1))
      SN=-1.
      IF(H(I,J+1).GT.H(I,J)) SN=1.
      H(I,J+1)=H(I,J+1)-SN*DH
      H(I,J)=H(I,J)+SN*DH
   11 CONTINUE
C   sweep down
      DO 12 J=JM-1,2,-1
      IF(FSM(I,J).EQ.0..OR.FSM(I,J+1).EQ.0.) GOTO 12 
      SL(I,J)=ABS(H(I,J+1)-H(I,J))/(H(I,J)+H(I,J+1))
      IF(SL(I,J).LT.SLMIN) GOTO 12
      DH=0.5*(SL(I,J)-SLMIN)*(H(I,J)+H(I,J+1))
      SN=-1.
      IF(H(I,J+1).GT.H(I,J)) SN=1.
      H(I,J+1)=H(I,J+1)-SN*DH
      H(I,J)=H(I,J)+SN*DH
  12  CONTINUE
  13  CONTINUE
C
  100 CONTINUE 
      RETURN
      END
C
      SUBROUTINE tidegen(nameh)
       INCLUDE 'comblk.for'
      real sl(im,jm)
	DATA PI/3.141592654/,SMALL/1.E-10/DG/111111.1/
	character*20 nameh
      grav=9.81 
c      
      RI=PI/180.0
c        OPEN(1,FILE='Baltic_F.grd')
        OPEN(1,FILE=nameh)
	i251=imm1

c        OPEN(7,FILE='amp.txt')
c        OPEN(8,FILE='faz.txt')
        READ(1,*)
        READ(1,*)Nx,Ny
        READ(1,*)along1,along2
        READ(1,*)alat1,alat2
        READ(1,*)
        do j=2,jmm1 
	  READ(1,*) (H(I,J),I=2,IMm1)
	  end do
        CLOSE(1)
      dlat=(alat2-alat1)/(ny-1.0)
      dlong=(along2-along1)/(nx-1.0)
c     surrounded bourder
      do i=2,imm1
          h(i,jm)=h(i,jmm1)
          h(i,1)=h(i,2)
	end do
      do j=1,jm
		h(im,j)=h(imm1,j)
		h(1,j)=h(2,j)
      end do

      xmi=along1-dlong
      xma=along2+dlong

      ymi=alat1-dlat
      yma=alat2+dlat
      
      hmax=0.0
      kmin=0
	kmax=0
	hmin=4.0        
c	  tide_l0=0.5
        DO J=1,JM
        DO I=1,IM
          if (h(i,j).le.hmin.and.h(i,j).gt.0.5) then
	      kmin=kmin+1
	      h(i,j)=hmin
	    else
            if (h(i,j).le.0.5) then
              h(i,j)=1.5
            end if 
		end if 	 
        end do
	  end do
      
        hmax=0.0
      DO 12 J=1,JM
        yy=cos(ri*(alat1+dlat*(j-0.5)))
        DY(J)=dlat*dg
        DX(J)=dlong*dg*yy
        ART(J)=DX(J)*DY(J)
        ARU(j)=ART(J)
        ARV(j)=ART(J)
        COR(J)=pi/12.0/3600/yy
  12  CONTINUE
      DTMAX=10000
      do i=1,im
      do j=1,jm
      hmax=max(h(i,j),hmax)
      end do
      end do

      do j=1,jm
      DTM=DX(J)/SQRT(Hmax*GRAV)
      DTMAX=MIN(DTM,DTMAX)
      end do

      dte=0.5*dtmax   !!!! dte determination
      n=3600/dte+1
      dte=real(3600.0)/real(n)
C
C
      IP=IM
      ISK=2
      JSK=2
C     CALL PRXY('  TOPOGRAPHY 1 ',TIME, H,IP,ISK,JM,JSK,10.)
C
      DO J=1,JM
      DO I=1,IM
      FSM(I,J)=0.
      DUM(I,J)=0.
      DVM(I,J)=0.
      IF(H(I,J).GT.1.6) FSM(I,J)=1.
      ENDDO
      ENDDO
      DO J=2,JM
      DO I=2,IM
      DUM(I,J)=FSM(I,J)*FSM(I-1,J)                                     
      DVM(I,J)=FSM(I,J)*FSM(I,J-1)                                     
      ENDDO
      ENDDO


c      call SLPMIN(H,IM,JM,FSM,SL)



	   
C
      DO 51 J=1,JM
      DO 51 I=1,IM
      UAB(I,J)=0.0   
      VAB(I,J)=0.0
      ELB(I,J)=0.0           
      AAM2D=2
   51 CONTINUE
C------------------------------------------------
C  Set lateral boundary conditions for use in S.R.BCOND. 
C  In the seamount problem the east and west BCs are open 
C  whereas the south and north BCs are closed through specification
C  of the masks, FSM, DUM and DVM.
C------------------------------------------------
C
      RETURN
      END


      SUBROUTINE GAUSS(IX,S,AM,V)
c      use dfport
      v1=rand(0)
      v2=rand(0)
c      V1=URAND(IX)
c      V2=URAND(IX)
      V=S*SQRT(ABS((-2*ALOG(V1+1.E-20))))*COS(6.28*V2)+AM
      RETURN
      END




      SUBROUTINE READGRD(NX,NY,NDX,Z,XMI,XMA,YMI,YMA,NAME)
	CHARACTER*20 NAME
	real z(NDX,*)
	OPEN(1,FILE=NAME)
	READ(1,*)
	READ(1,*)NX,NY
	READ(1,*)XMI,XMA
	READ(1,*)YMI,YMA
	READ(1,*)
	READ(1,*)((z(i,j),i=1,nx),j=1,ny) 
	CLOSE(1)
	RETURN
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

      SUBROUTINE FFREADDIMGR3(NX,NY,NZ,NAME)
	CHARACTER*20 NAME
      INTEGER NX,NY,NZ

      OPEN(1,FILE=NAME)
	READ(1,*)
	READ(1,*)NX,NY,NZ
	CLOSE(1)
	RETURN
	END
      
      SUBROUTINE FFREADGR3(NX,NY,NZ,XMI,XMA,YMI,YMA,ZMI,ZMA,NAME,Z)
	REAL Z(NX,NY,NZ),XMI,XMA,YMI,YMA,ZMI,ZMA
      INTEGER NX,NY,NZ,I,J,K
      CHARACTER*20 NAME

      OPEN(1,FILE=NAME)
	READ(1,*)
	READ(1,*)NX,NY,NZ
      READ(1,*) XMI,XMA
      READ(1,*) YMI,YMA
      READ(1,*) ZMI,ZMA
	READ(1,*)(((Z(I,J,K),I=1,NX),J=1,NY),K=1,NZ)

	CLOSE(1)
	RETURN
	END
      