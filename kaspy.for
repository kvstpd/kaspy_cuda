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
      DIMENSION 
     1     ADVUA(IM,JM),ADVVA(IM,JM),
     2     surf(im,jm)
      COMMON/F_FLOATS/ff_marker,
     3     fxf(im,jm),fyf(im,jm),fxb(im,jm),fyb(im,jm),
     4     FF(IM,JM),FB(IM,JM),
     7     ffu(im,jm),ffv(im,jm),fbu(im,jm),fbv(im,jm),
     8     ff_end
      real*8
     5     SEL(IM,JM),SSEL(IM,JM),SFA(IM,JM),SSFA(IM,JM),SFEL(IM,JM),
     6     SFAR(IM,JM),SSFAR(IM,JM),SFELR(IM,JM),
     8     su(im,jm),sv(im,jm),ssu(im,jm),ssv(im,jm),
     9     ssuv(im,jm),ssue(im,jm),ssve(im,jm)
      
      REAL, ALLOCATABLE :: PRESS(:,:,:),uwd(:,:,:),vwd(:,:,:)
      REAL, ALLOCATABLE :: PRESS0(:,:),uwd0(:,:),vwd0(:,:)

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
      real*8 time_
      DATA ri/0.01745329252/,GEE/9.807/
      CHARACTER*20 NAME,namep,nameu,namev,nameh
      integer ijloc(2)


      
       ittest = 99


       write (*,*) 'tttt=',ittest

       ittest = test_me(ittest)

        vars_marker = 0.89898
                
        arrays_marker = 3.141592653589793238462643383279501010101010
        arrays_end_marker = 0.9876543211234567890

        ff_end = 100.001

	icycler =cycler_init(vars_marker, arrays_marker,ff_marker)
	

c	call cycler_destroy(icycler)

c	stop




      namep='pr1979.GR3'
      CALL READDIMGR3(KX,KY,KT,namep)
      ALLOCATE (PRESS(KX,KY,KT),PRESS0(KX,KY))
      CALL READGR3(KX,KY,KT,XKI,XKA,YKI,YKA,TKI,TKA,namep,PRESS)
ccc      NH6=KT   1 DURATION !!!
      NH6=100    !
      PRESS=PRESS/1000
   
      dht=(tka-tki)/(kt-1)
      
C     READ wind DATA
      nameu='u1979.GR3'
      CALL READDIMGR3(KXU,KYU,KTU,nameu)
      ALLOCATE (UWD(KXU,KYU,KTU),UWD0(KXU,KYU))
      CALL READGR3
     1 (KXU,KYU,KTU,XKUI,XKUA,YKUI,YKUA,TKUI,TKUA,nameu,UWD)
C  ÂÅÒÅÐ Â 0 !!!!!!
C      UWD=0
      
      namev='v1979.GR3'
      CALL READDIMGR3(KXV,KYV,KTV,namev)
      ALLOCATE (VWD(KXV,KYV,KTV),VWD0(KXV,KYV))
      CALL READGR3
     1 (KXV,KYV,KTV,XKVI,XKVA,YKVI,YKVA,TKVI,TKVA,namev,VWD)
      VWD=0
C
C     READ IN GRID DATA AND INITIAL AND LATERAL BOUNDARY CONDITIONS
C------------------------------------------------------------------------
      nameh='kaspi_1_5mR.grd'
      CALL TIDEGEN(nameh)
C
C********************************************************************
C  Note that lateral thermodynamic boundary conditions are often set equal
C  to the initial conditions and are held constant thereafter. Users can
C  of course create varable boundary conditions.
C********************************************************************

C     END READING DATA      
      
	open(111,file='s66.txt')
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
	
	call cycler_time()

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
      
      ro_ratio=1.29/1020.0
      
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
      open(88,file='force1.txt')

C
c      TIME=TIMEI
C***********************************************************************
C                                                                      *
C                  BEGIN NUMERICAL INTEGRATION                         *
C                                                                      *
C***********************************************************************
C

cc!        display = display_init_f()

cc!         display_data_width = im
cc!         display_data_height = jm


cc!        display_res =  display_load_data_f(display, 
cc!     & display_data_width, display_data_height, 
cc!     & surf, fsm, display_scale)


c call display(surf,im,jm,im,jm,-1.0,1.0,0)
        iold=0

                     DO 9000 IINT=1,IEND
C

      timeh=(iint-1)*dti/3600.0
      ihour_s=600/dti
c      set uwind and current,
c	 if time less than time0 acceleration otherwise constant current and wind
c      if (time.le.1) then	    	
c	   uwind=xwind*time
c	   vwind=ywind*time
c	   tide_l=0.5*ampl*(time/time0)**2 
c	else
c	   mtime=time
c	   uwind=xwind
c	   vwind=ywind
c         tide_l=ampl*(time/time0-0.5)
c	end if
	iwrite=mod(iint,ihour_s)

      if(iwrite.eq.1) then

      elfmax=0.0
	elfmin=0.0
	DO J=1,JM
      DO I=1,IM
        if (elfmax.lt.elf(i,j)) elfmax=elf(i,j)
        if (elfmin.gt.elf(i,j)) elfmin=elf(i,j)
      end do
c      call display_show_f(display)
	end do




        time_=timeh
        do i=1,im
        do j=1,jm
          if(fsm(i,j).gt.0) then
c            surf(i,j)=10*uf(i,j,1)
            surf(i,j)=elf(i,j)
          else
            surf(i,j)=-1.0
          end if
        end do
c        call display_show_f(display)
        end do
c        call display(surf,im,jm,im,jm,-1.0,1.0,0)


 

cc!        call display_update_f(display)

cc!        call display_show_f(display)
    
        write(6,1117) 't=',timeh,'h','Sea level=',elfmax,elfmin,'m'


1117   FORMAT(a3,f12.4,1x,a1,5x,a10,f8.4,1x,f8.4,1x,a1,a1)
      end if
c*************************************************
C------------------------------------------------------------------------
C
C------------------------------------------------------------------------

C********** BEGIN EXTERNAL MODE ***************************************
ccc                      DO 8000 IEXT=1,ISPLIT
c      TIMEH=(DTI*FLOAT(IINT-1)+(iext-1)*dte)/3600.
      itimeh=int(timeh)


CCCC      goto 8787   ! ÎÁÕÎÄ ÑÒÀÒÈÑÒÈÊÈ
     	if(itimeh.gt.iold) then  !! writing to file and compute statiatics, 
c                              !!    here every hour


	write(77,'(101f10.3)') timeh,(elf(stx(kk),sty(kk)),kk=1,nstation)
C     AFSM - NUMBER OF NON-ZERO POINTS IN THE FSM ARRAY
C==============================================CALCULATING STATISTICS===================================
        SSPRE=0
        DO J=1,JM
        DO I=1,IM
          FA=(btim*FB(i,j)+ftim*FF(i,j)-100)/10.0
          SSPRE=SSPRE+FA*FSM(I,J)/AFSM
          SFA(I,J)=SFA(I,J)+FA
          SSFA(I,J)=SSFA(I,J)+FA**2
          SEL(I,J)=SEL(I,J)+EL(I,J)
          SSEL(I,J)=SSEL(I,J)+EL(I,J)**2
          SFEL(I,J)=SFEL(I,J)+FA*EL(I,J)
          
c	    uw=0.9*(btim*fbu(i,j)+ftim*ffu(i,j))  
c	    vw=0.9*(btim*fbv(i,j)+ftim*ffv(i,j))

	    uw=(btim*fbu(i,j)+ftim*ffu(i,j))  
	    vw=(btim*fbv(i,j)+ftim*ffv(i,j))
	      
          su(i,j)=su(i,j)+uw
          sv(i,j)=sv(i,j)+vw
          ssu(i,j)=ssu(i,j)+uw**2
          ssv(i,j)=ssv(i,j)+vw**2
          ssuv(i,j)=ssuv(i,j)+uw*vw
          ssue(i,j)=ssue(i,j)+uw*el(i,j)
          ssve(i,j)=ssve(i,j)+vw*el(i,j)
        END DO
c        call display_show_f(display)
        END DO 
        DO J=1,JM
        DO I=1,IM
          FAR=(btim*FB(i,j)+ftim*FF(i,j)-100)/10.0-SSPRE
          SFAR(I,J)=SFAR(I,J)+FAR
          SSFAR(I,J)=SSFAR(I,J)+FAR**2
          SFELR(I,J)=SFELR(I,J)+FAR*EL(I,J)
        END DO
c        call display_show_f(display)
        END DO 
      
        NSTAT=NSTAT+1
	end if
c------------------------------------------------end if statistics------------------------------------
ccccc    8787  continue     ! Êîíåö îáõîäà ñòàòèñòèêè 
      timeh6=timeh/dht+1
      itime6=timeh6
      ftim=0
      ftim=(timeh6-itime6)
      btim=1.0-ftim
      
c================================================renew forcing fields ==================================
      if (itime6.gt.itime6_old) then
        itime6_old=itime6
        fxb=fxf
        fyb=fyf
        FB=FF
        fbu=ffu
        fbv=ffv
        press0(:,:)=press(:,:,itime6)
        call getnewpressureVAR(kx,ky,XKI,XKA,YKI,YKA,PRESS0,
     1 FF,fxf,fyf)                
        uwd0(:,:)=uwd(:,:,itime6)
        call getnewwindVAR(kxu,kyu,XKUI,XKUA,YKUI,YKUA,uwd0,ffu)
        vwd0(:,:)=vwd(:,:,itime6)
        call getnewwindVAR(kxv,kyv,XKVI,XKVA,YKVI,YKVA,vwd0,ffv)
        write(88,*) fyf(1,1),fyf(im,jm),fxf(1,1),fxf(im,jm)
      end if  
c----------------------------------------------------end of renew forcing-------------------------------


c====================================================computing wusurf(i,j), wvsurf/(i,j)=================
      DO J=2,JMM1
      DO I=2,IMM1
	uw=(btim*fbu(i,j)+ftim*ffu(i,j))
	vw=(btim*fbv(i,j)+ftim*ffv(i,j))  
	speed=sqrt(uw**2+vw**2) !******************************************************
!      speed=0
	windc=1.0e-3*(0.8+speed*0.065)*ro_ratio*speed  
      WUSURF(I,J)=-windc*uw
     1 	*.25*(DUM(I,J+1)+DUM(I+1,J)+DUM(I-1,J)+DUM(I,J-1))+
     2  0.5*(d(i,j)+d(i-1,j))*(btim*FxB(i,j)+ftim*FxF(i,j))
      WVSURF(I,J)=-windc*vw
     1 	*.25*(DVM(I,J+1)+DVM(I+1,J)+DVM(I-1,J)+DVM(I,J-1))+
     2  0.5*(d(i,j)+d(i,j-1))*(btim*FyB(i,j)+ftim*FyF(i,j))
      end do
      end do

c     write(6,'('' IEXT,TIME ='',I5,F9.2)') IEXT,TIME      

c======================405 cycle and 410 cycle: computing continuity equation
      DO 405 J=2,JM
      DO 405 I=2,IM
      FLUXUA(I,J)=.25E0*(D(I,J)+D(I-1,J))*(DY(j)+DY(j))*UA(I,J)
 405  FLUXVA(I,J)=.25E0*(D(I,J)+D(I,J-1))*(DX(j)+DX(j-1))*VA(I,J)
C
      DO 410 J=2,JMM1
      DO 410 I=2,IMM1
  410 ELF(I,J)=ELB(I,J)
     1    -DTE2*(FLUXUA(I+1,J)-FLUXUA(I,J)+FLUXVA(I,J+1)-FLUXVA(I,J))
     2                    /ART(J)
C
c--------------------end of continuity equation----------------------------------

      CALL BCOND(1) !!boundary condition: elevevation
C
      IF(MOD(IINT,ISPADV).EQ.0) CALL ADVAVE(ADVUA,ADVVA,MODE)
C  Note that ALPHA = 0. is perfectly acceptable. The value, ALPHA = .225
C  permits a longer time step.
c=======================main momentum update: cycles 420 and 430=================
      ALPHA=0.225     
      DO 420 J=2,JMM1
      DO 420 I=2,IM
      UAF1=ADVUA(I,J)   
     1    -.25*(COR(j)*D(I,J)*(VA(I,J+1)+VA(I,J))
     2              +cor(j)*D(I-1,J)*(VA(I-1,J+1)+VA(I-1,J)) )
     3         +.5E0*GRAV*dy(j)/aru(j)*(D(I,J)+D(I-1,J))
     4             *( (1.E0-2.E0*ALPHA)*(EL(I,J)-EL(I-1,J))
     4            +ALPHA*(ELB(I,J)-ELB(I-1,J)+ELF(I,J)-ELF(I-1,J)) )
     6      +WUSURF(I,J)-WUBOT(I,J)   
      UAF(I,J)=
     1         ((H(I,J)+ELB(I,J)+H(I-1,J)+ELB(I-1,J))*UAB(I,J)
     2                -4.E0*DTE*UAF1)
     3        /(H(I,J)+ELF(I,J)+H(I-1,J)+ELF(I-1,J))
420    CONTINUE
      DO 430 J=2,JM
      DO 430 I=2,IMM1
      VAF1=ADVVA(I,J)  
     1    +.25*(  COR(J)*D(I,J)*(UA(I+1,J)+UA(I,J))
     2               +COR(J-1)*D(I,J-1)*(UA(I+1,J-1)+UA(I,J-1)) )
     3         +.5E0*GRAV*DX(j)/arv(j)*(D(I,J)+D(I,J-1))
     4             *( (1.E0-2.E0*ALPHA)*(EL(I,J)-EL(I,J-1))
     4            +ALPHA*(ELB(I,J)-ELB(I,J-1)+ELF(I,J)-ELF(I,J-1)) )
     6    + WVSURF(I,J)-WVBOT(I,J)   

      VAF(I,J)=
     1        ((H(I,J)+ELB(I,J)+H(I,J-1)+ELB(I,J-1))*VAB(I,J)
     2              -4.E0*DTE*VAF1)
     3       /(H(I,J)+ELF(I,J)+H(I,J-1)+ELF(I,J-1))

c      call display_show_f(display)
  430 CONTINUE
      CALL BCOND(2)  ! boundary conditions for velocities
C
C
C  TEST FOR CFL VIOLATION. IF SO, PRINT AND STOP
C
      VMAXL=100.
      VAMAX=0.       
      UAMAX=0.       
      tps=sqrt(uaf**2+uvf**2)
      uamax=maxval(tps)
      ijloc=maxloc(tps)
      imax=ijloc(1)
      jmax=ijloc(2)
c      DO 442 J=1,JM
c      DO 442 I=1,IM
c      IF(ABS(VAF(I,J)).GE.VAMAX) THEN
c        VAMAX=ABS(VAF(I,J))   
c	  IMAX=I
c	  JMAX=J
c	end if
c      IF(ABS(UAF(I,J)).GE.UAMAX) THEN
c        UAMAX=ABS(UAF(I,J))   
c	  IMAX=I
c	  JMAX=J
c      ENDIF
c  442 CONTINUE
c      IF(VAMAX.GT.VMAXL) GO TO 9001
      IF(UAMAX.GT.VMAXL) GO TO 9001
C    
C       APPLY FILTER TO REMOVE TIME SPLIT. RESET TIME SEQUENCE.
      DO 445 J=1,JM
      DO 445 I=1,IM
      UA(I,J)=UA(I,J)+.5E0*SMOTH*(UAB(I,J)-2.E0*UA(I,J)+UAF(I,J))
      VA(I,J)=VA(I,J)+.5E0*SMOTH*(VAB(I,J)-2.E0*VA(I,J)+VAF(I,J))
      EL(I,J)=EL(I,J)+.5E0*SMOTH*(ELB(I,J)-2.E0*EL(I,J)+ELF(I,J))
      ELB(I,J)=EL(I,J)
      EL(I,J)=ELF(I,J)
      D(I,J)=H(I,J)+EL(I,J)
      UAB(I,J)=UA(I,J)
      UA(I,J)=UAF(I,J)
      VAB(I,J)=VA(I,J)
      VA(I,J)=VAF(I,J)
  445 CONTINUE
C

c      print external boundary conditions

c      writing boundary condition

C
C
       iold=itimeh

cc!       call display_show_f(display)

cccc 8000                    CONTINUE
      IF(VAMAX.GT.VMAXL) GO TO 9001
C---------------------------------------------------------------------
C          END EXTERNAL (2-D) MODE CALCULATION
C     AND CONTINUE WITH INTERNAL (3-D) MODE CALCULATION
C---------------------------------------------------------------------
      IF(IINT.EQ.1) GO TO 8200
      IF(MODE.EQ.2) GO TO 8200
C
 8200 CONTINUE
C
C
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
 9001 CONTINUE
      IF(VAMAX.GT.VMAXL) THEN
      write(*,*) imax,jmax
      write(*,*) 'vamax>vmax!!!'
c      STOP
      ENDIF
      IF(UAMAX.GT.VMAXL) THEN
      write(*,*) imax,jmax
       write(*,*) 'uamax>umax!!!'
c      STOP
      ENDIF
 7000 CONTINUE
C--------------------------------------------------------------------
C             END PRINT SECTION
C--------------------------------------------------------------------
 9000                     CONTINUE


cc!        call display_destroy_f(display)
         
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
      DO J=1,JM
      DO I=1,IM
        SFA(I,J)=SFA(I,J)/NSTAT      
        SEL(I,J)=SEL(I,J)/NSTAT      
        SFAR(I,J)=SFAR(I,J)/NSTAT      
        SSFA(I,J)=(SSFA(I,J)/NSTAT-SFA(I,J)**2)*1.0E4
        SSEL(I,J)=(SSEL(I,J)/NSTAT-SEL(I,J)**2)*1.0E4
        SFEL(I,J)=(SFEL(I,J)/NSTAT-SFA(I,J)*SEL(I,J))*1.0E4
        SSFAR(I,J)=(SSFAR(I,J)/NSTAT-SFAR(I,J)**2)*1.0E4
        SFELR(I,J)=(SFELR(I,J)/NSTAT-SFAR(I,J)*SEL(I,J))*1.0E4
        SFA(I,J)=SFEL(I,J)/(SSFA(I,J)+EPS)
c        SEL(I,J)=SQRT(SSEL(I,J)/(SSFA(I,J)+EPS))
        SFAR(I,J)=SFELR(I,J)/(SSFAR(I,J)+EPS)
        
        su(i,j)=su(i,j)/nstat
        sv(i,j)=sv(i,j)/nstat
        ssu(i,j)=(ssu(i,j)/NSTAT-su(i,j)**2)
        ssv(i,j)=(ssv(i,j)/NSTAT-sv(i,j)**2)
        ssuv(i,j)=(ssuv(i,j)/NSTAT-sv(i,j)*su(i,j))
        ssue(i,j)=(ssue(i,j)/NSTAT-su(i,j)*sel(i,j))*100
        ssve(i,j)=(ssve(i,j)/NSTAT-sv(i,j)*sel(i,j))*100                
      END DO
      END DO
c      NAME='SSFA.GRD'
c      CALL WRITEGRD(IM,JM,IM,SSFA,10.0,31.0,53.0,66.0,NAME)
      NAME='SSEL.GRD'
      CALL WRITEGRD(IM,JM,IM,SSEL,xmi,xma,ymi,yma,NAME)
c      NAME='SFEL.GRD'
c      CALL WRITEGRD(IM,JM,IM,SFEL,10.0,31.0,53.0,66.0,NAME)

c      NAME='ALFA1.GRD'
c      CALL WRITEGRD(IM,JM,IM,SFA,10.0,31.0,53.0,66.0,NAME)
c      NAME='ALFA2.GRD'
c      CALL WRITEGRD(IM,JM,IM,SEL,10.0,31.0,53.0,66.0,NAME)
      
      NAME='SSFAR.GRD'
      CALL WRITEGRD(IM,JM,IM,SSFAR,xmi,xma,ymi,yma,NAME)
      NAME='SFELR.GRD'
      CALL WRITEGRD(IM,JM,IM,SFELR,xmi,xma,ymi,yma,NAME)
c      NAME='ALFAR.GRD'
c      CALL WRITEGRD(IM,JM,IM,SFAR,10.0,31.0,53.0,66.0,NAME)
      NAME='ssu.GRD'
      CALL WRITEGRD(IM,JM,IM,ssu,xmi,xma,ymi,yma,NAME)
      NAME='ssv.GRD'
      CALL WRITEGRD(IM,JM,IM,ssv,xmi,xma,ymi,yma,NAME)
      NAME='ssuv.GRD'
      CALL WRITEGRD(IM,JM,IM,ssuv,xmi,xma,ymi,yma,NAME)
      NAME='ssue.GRD'
      CALL WRITEGRD(IM,JM,IM,ssue,xmi,xma,ymi,yma,NAME)
      NAME='ssve.GRD'
      CALL WRITEGRD(IM,JM,IM,ssve,xmi,xma,ymi,yma,NAME)
      
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
      SUBROUTINE ADVAVE(ADVUA,ADVVA,MODE)
       INCLUDE 'comblk.for'
      DIMENSION ADVUA(IM,JM),ADVVA(IM,JM)
C---------------------------------------------------------------------
C      CALCULATE U-ADVECTION & DIFFUSION
C---------------------------------------------------------------------
c      cycles: 7 cycles (460 470, 480, 860, 870, 880, 102) plus initialization with zeros
c      of the arrays advua, fluxua, advva, fluxva
c      writing to fluxua(i,j): 460 cycle, 870 cycle
c      writing to fluxva(i,j): 470 cycle, 860 cycle
c      writing to tps(i,j):  470 cycle
c      writing to advua(i,j) 480 cycle
c      writing to advva(i,j): 880 cycle
c      writing to wubot(i,j): 102 cycle
c      writing to wvbot(i,j): 102 cycle
c
c      read data from next arrays:
c      D(i,j)- full depth (h+elevation)
c      ua(i,j) - u velocity (central time layer)
c      uab(i,j) - u velocity (old time layer)
c      va(i,j) - v velocity (central time layer)
c      vab(i,j) - v velocity (old time layer)
c      cbc(i,j)  - empirical parameter of bottom roughness 
c      dxx(j),dy(j) - spatial steps (depends on j only)
c      aru(j), arv(j) - area u, area v (depends on j only)
c      aam2d  - horizondal friction (mm/s, here is a constant==2 

C
C-------- ADVECTIVE FLUXES INITIALIZATION-------------------------------------
C
      ADVUA=0
      FLUXUA=0

C----------- ADD VISCOUS FLUXES ---------------------------------
      DO 460 J=2,JM
      DO 460 I=2,IMM1
 460  FLUXUA(I,J)=DY(J)*(.125E0*((D(I+1,J)+D(I,J))*UA(I+1,J)
     1                 +(D(I,J)+D(I-1,J))*UA(I,J))
     2                  *(UA(I+1,J)+UA(I,J))
     3         -D(I,J)*2.E0*AAM2D*(UAB(I+1,J)-UAB(I,J))/DX(j))
      DO  470 J=2,JM
      DO  470 I=2,IM

      TPS(I,J)=(D(I,J)+D(I-1,J)+D(I,J-1)+D(I-1,J-1))
     1            *AAM2D
     2            *((UAB(I,J)-UAB(I,J-1))
     3                /(4*DY(j))
     4                 +(VAB(I,J)-VAB(I-1,J))
     5                /(4*DX(j)) )
      FLUXVA(I,J)=(.125E0*((D(I,J)+D(I,J-1))*VA(I,J)
     1                 +(D(I-1,J)+D(I-1,J-1))*VA(I-1,J))
     2                    *(UA(I,J)+UA(I,J-1))
     3 -TPS(I,J))*DX(j)
 470  CONTINUE
C----------------------------------------------------------------
      DO  480 J=2,JMM1
      DO  480 I=2,IMM1
 480  ADVUA(I,J)=(FLUXUA(I,J)-FLUXUA(I-1,J)
     1           +FLUXVA(I,J+1)-FLUXVA(I,J))/aru(j)
C----------------------------------------------------------------
C       CALCULATE V-ADVECTION & DIFFUSION
C----------------------------------------------------------------
C
      ADVVA=0.0
      FLUXVA=0.0
C
C------- ADD VISCOUS FLUXES -----------------------------------
      DO  860 J=2,JMM1
      DO  860 I=2,IM
 860  FLUXVA(I,J)=DX(J)*(.125E0*((D(I,J+1)+D(I,J))
     1       *VA(I,J+1)+(D(I,J)+D(I,J-1))*VA(I,J))
     2      *(VA(I,J+1)+VA(I,J))
     1        -D(I,J)*2.E0*AAM2D*(VAB(I,J+1)-VAB(I,J))/DY(j))
      DO  870 J=2,JM
      DO  870 I=2,IM
      FLUXUA(I,J)=(.125E0*((D(I,J)+D(I-1,J))*UA(I,J)
     1         +(D(I,J-1)+D(I-1,J-1))*UA(I,J-1))*
     2                        (VA(I-1,J)+VA(I,J))
     3  -TPS(I,J))*DY(j)
870   CONTINUE     
C---------------------------------------------------------------
      DO  880 J=2,JMM1
      DO  880 I=2,IMM1
 880  ADVVA(I,J)=(FLUXUA(I+1,J)-FLUXUA(I,J)
     1          +FLUXVA(I,J)-FLUXVA(I,J-1))/arv(j)
C
C---------------------------------------------------------------
      IF(MODE.NE.2) GO TO 5000
      DO 102 J=2,JMM1
      DO 102 I=2,IMM1
      WUBOT(I,J)=-0.5E0*(CBC(I,J)+CBC(I-1,J))
     1     *SQRT(UAB(I,J)**2+(.25E0*(VAB(I,J)
     2     +VAB(I,J+1)+VAB(I-1,J)+VAB(I-1,J+1)))**2)*UAB(I,J)
      WVBOT(I,J)=-0.5E0*(CBC(I,J)+CBC(I,J-1))
     1     *SQRT((.25E0*(UAB(I,J)+UAB(I+1,J)
     2     +UAB(I,J-1)+UAB(I+1,J-1)))**2+VAB(I,J)**2)*VAB(I,J)
  102 CONTINUE
  120 CONTINUE  
 5000 CONTINUE
C 
      RETURN
      END
C               
      SUBROUTINE BCOND(IDX)
       INCLUDE 'comblk.for'   
C
C  Closed boundary conditions are automatically enabled through
C  specification of the masks, DUM, DVM and FSM in which case open
C  boundary condition, included below, will be overwritten.
C
C
C*******************************************************************************
C                            POM (C-Grid)   
C                          ================
C
C    The diagram below and some of the text was provided by D.-S. Ko. 
C    It is for the case where U and V are the primary boundary conditions
C    together with T and S (co-located with E). E = EL itself is rather
C    unimportant and is substituted from an adjacent interior point.
C    
C    Inside ....... indicate the interior (non-boundary) grid points.
C    In general only those variables in the interior are computed and
C    variables at open boundary have to be specified.
C    All interpolations are centered in space except those at lateral
C    open boundary where an upstream scheme is usually used.
C
C    Horizontal locations of E(EL), T and S (etc.) are coincident.
C    NU = Not Used is indicated by *. However, for attractive output
C    adjacent interior values may be filled in at these points.   
C
C    I have been asked many times what kind of B.C. along the side wall
C    POM uses from people not acquainted with sigma coordinates. Although
C    the issue is not important as it might be for z-level grids, a direct
C    answer is "half slip" which, of course, is between free slip and
C    non-slip B.C.s.
C
C-------------------------------------------------------------------------------
C     |                               N O R T H
C     |
C     |    1     2     3           I-1   I    I+1         IM-2  IM-1   IM 
C-----+-------------------------------------------------------------------------
C     |   NU BC BC                                                    BC BC
C     |    v  v  v                                                     v  v
C     |
C     |BC> U* E  U  E  U  E  .  .  U  E  U  E  U  E  .  .  U  E  U  E  U  E  <BC
C     |    |     |     |           |     |     |           |     |     |      
C  JM |BC> +--V--+--V--+--V--   .  +--V--+--V--+--V--   .  +--V--+--V--+--V- <BC
C     |    |     | ....|...........|.....|.....|...........|.....|.... |      
C     |    U* E  U :E  U  E  .  .  U  E  U  E  U  E  .  .  U  E  U  E: U  E 
C     |    |     | :   |           |     |     |           |     |   : |      
C JM-1|    +--V--+--V--+--V--   .  +--V--+--V--+--V--   .  +--V--+--V--+--V-
C     |    |     | :   |           |     |     |           |     |   : |
C     |    U* E  U :E  U  E  .  .  U  E  U  E  U  E  .  .  U  E  U  E: U  E
C     |    |     | :   |           |     |     |           |     |   : |
C JM-2|    +--V--+--V--+--V--   .  +--V--+--V--+--V--   .  +--V--+--V--+--V-
C     |            :                                                 :          
C W   |       .  . :.  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .: .  .    E
C E   |            :                    Interior                     :         A
C S   |       .    :.     .     .     .     .     .     .     .     .:    .    S
C T   |            :                                                 :         T
C     |    |     | :   |           |     |     |           |     |   : |
C     |    U* E  U :E  U  E  .  .  U  E  U  E  U  E  .  .  U  E  U  E: U  E
C     |    |     | :   |           |     |     |           |     |   : |
C   3 |    +--V--+--V--+--V--   .  +--V--+--V--+--V--   .  +--V--+--V--+--V-
C     |    |     | :   |           |     |     |           |     |   : |      
C     |    U* E  U :E  U  E  .  .  U  E  U  E  U  E  .  .  U  E  U  E: U  E
C     |    |     | ....|...........|.....|.....|...........|.....|.... |      
C   2 |BC> +--V--+--V--+--V--   .  +--V--+--V--+--V--   .  +--V--+--V--+--V- <BC
C     |    |     |     |           |     |     |           |     |     |      
C     |BC> U* E  U  E  U  E  .  .  U  E  U  E  U  E  .  .  U  E  U  E  U  E  <BC
C     |    |     |     |           |     |     |           |     |     |      --
C   1 |NU> +--V*-+--V*-+--V*-   .  +--V*-+--V*-+--V*-   .  +--V*-+--V*-+--V* <NU
C     |
C     |    ^  ^  ^                                                     ^  ^
C     |   NU BC BC                                                    BC BC
C-----+-------------------------------------------------------------------------
C     |    1     2     3           I-1   I     I+1         IM-2  IM-1  IM
C     |
C     |                                S O U T H
C-------------------------------------------------------------------------------
C
C
C
	real*8 time_,timei
C
C
      GO TO (10,20,30), IDX
C
C-----------------------------------------------------------------------
C                   EXTERNAL B.C.'S
C  In this example the governing boundary conditions are a radiation   
C  condition on UAF in the east and in the west and VAF in the north
C  and south. The tangential velocities are set to zero on both boundaries.
C  These are only one set of possibilities and may not represent a choice
C  which yields the most physically realistic result.
C-----------------------------------------------------------------------
C
C
 10   CONTINUE
C---------ELEVATION --------------------
C  A4 BC     
c      tide_l=0.5
	DO 120 J=1,JM
	elf(2,j)=tide_l
	elf(imm1,j)=tide_l
      ELF(1,J)=ELF(2,J)
      ELF(IM,J)=elf(imm1,j)
120   CONTINUE
c   A2 set
      DO 130 I=1,IM
c        ELF(I,1)=TIDE_L
        ELF(I,1)=elf(i,2)
c      ELF(I,JM)=tide_l
      elf(i,jm)=elf(i,jmm1)
130   CONTINUE
C
      DO 140 J=1,JM
      DO 140 I=1,IM
 140  ELF(I,J)=ELF(I,J)*FSM(I,J)
      RETURN
C
C
 20   CONTINUE
C---------- VELOCITY --------------        
C          A4 set
      DO 210 J=2,JMM1
C            EAST
      if(dum(im,j).gt.0.5) then
      GAE=DTE*SQRT(GRAV*H(IM,j))/DX(j)
	UAF(IM,J)=GAE*UA(IMM1,j)+(1.-GAE)*UA(IM,j)
      else
	uaf(im,j)=0.0
	end if   
c	VAF(I,JM)=0.0
c      UAF(IM,J)=UAF(IMM1,J)
      VAF(IM,J)=0.0
C            WEST
      if(dum(2,j).gt.0.5) then  
      GAE=DTE*SQRT(GRAV*H(2,j))/DX(j)
	UAF(2,J)=GAE*UA(3,j)+(1.-GAE)*UA(2,j)
	else
	uaf(2,j)=0.0
	end if
c      UAF(2,J)=UAF(3,J)
      UAF(1,J)=UAF(2,J)
      VAF(1,J)=0.0

 210  CONTINUE
C-----------------------------------
c      A4 set
      DO 220 I=2,IMM1
C            NORTH
      if (dvm(i,jm).gt.0.5) then 
      GAE=DTE*SQRT(GRAV*H(I,JM))/DY(JM)
	VAF(I,JM)=GAE*VA(I,JMM1)+(1.-GAE)*VA(I,JM)
	else
	vaf(i,jm)=0.0
	end if
	UAF(I,JM)=0.0
C            SOUTH
      if (dvm(i,2).gt.0.5) then 
      GAE=DTE*SQRT(GRAV*H(I,2))/DY(1)
	VAF(I,2)=GAE*VA(I,3)+(1.-GAE)*VA(I,2)
	else
	vaf(i,2)=0.0
	end if
	vaf(i,1)=vaf(i,1)
	UAF(I,1)=0.0
 220  CONTINUE
C---------------------------
C
      DO 240 J=1,JM
      DO 240 I=1,IM
      UAF(I,J)=UAF(I,J)*DUM(I,J)
 240  VAF(I,J)=VAF(I,J)*DVM(I,J)
      RETURN
C
C
 30   CONTINUE
      RETURN
      END
C
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

      subroutine getnewpressureVAR(kx,ky,XKI,XKA,YKI,YKA,PRESS0,
     1 FF,fxf,fyf)
      INCLUDE 'comblk.for'
      real fxf(im,jm),fyf(im,jm),PRESS0(KX,KY),FF(IM,JM)
      CALL GETPRESScube(KX,KY,KX,PRESS0,XKI,XKA,YKI,YKA,
     1                     IM,JM,IM,FF,FXF,FYF,XMI,XMA,YMI,YMA)    
      return
      end

      subroutine getnewwindVAR(kxw,kyw,XKI,XKA,YKI,YKA,uw0,ffu)
      INCLUDE 'comblk.for'
      real uw0(KXw,KYw),ffu(IM,JM)
      CALL GETcube(KXw,KYw,KXw,uw0,XKI,XKA,YKI,YKA,
     1                     IM,JM,IM,FFU,XMI,XMA,YMI,YMA)
      return
      end

	  		           
	subroutine getCUBE(kx,ky,kd,pk,xki,xka,yki,yka,
     1                     nx,ny,nd,P,xmi,xma,ymi,yma)
c     interpolate pressure field to get derivatives.
c     inputs - pk(kx,ky)
c     inputs - xki,xka,yki,yka
c     inputs - xmi,xma,ymi,yma
c     inputs - Kx,Ky,Nx,Ny
c------------------------------------------------------------
c      output - px(Nx,Ny),Py(Nx,Ny)
	real pk(kd,*),P(ND,*)
	real PKK(50,50),C(4,4,50,50)


	dky=(yka-yki)/(ky-1)
	dkx=(xka-xki)/(kx-1)
      
	dy=(yma-ymi)/(ny-1)
	dx=(xma-xmi)/(nx-1)

C     SURROUNDING
	DO J=2,KY+1       
        DO I=2,KX+1
          PKK(I,J)=PK(I-1,J-1)
        END DO
      END DO
	DO J=2,KY+1       
        PKK(1,J)=2*PKK(2,J)-PKK(3,J)
        PKK(KX+2,J)=2*PKK(KX+1,J)-PKK(KX,J)
      END DO
	DO I=1,KX+2       
        PKK(I,1)=2*PKK(I,2)-PKK(I,3)
        PKK(I,KY+2)=2*PKK(I,KY+1)-PKK(I,KY)
      END DO
	CALL GETBICUBIC(KX+2,KY+2,50,PKK,C)

	do j=1,Ny
	  y=ymi+(j-1)*dy
	  j0=(y-yki)/dky+1
	  if (j0<1) j0=1
	  if (j0>ky-1) j0=ky-1
	  u=(y-(yki+(j0-1)*dky))/dky
	  
	  do i=1,Nx
	    x=xmi+(i-1)*dx
	    i0=(x-xki)/dkx+1
	    if (i0<1) I0=1
	    if (i0>kx-1) i0=kx-1
	    t=(x-(xki+(i0-1)*dkx))/dkx
          ay=0.
          DO K=4,1,-1
            ay=t*ay+((c(K,4,i0,j0)*u+c(k,3,i0,j0))*u+c(K,2,i0,j0))*u+
     1		  c(K,1,i0,j0)
          END DO
          p(i,j)=ay
	  end do
      END DO


	return
	end

	subroutine getpressCUBE(kx,ky,kd,pk,xki,xka,yki,yka,
     1                     nx,ny,nd,P,px,py,xmi,xma,ymi,yma)
c     interpolate pressure field to get derivatives.
c     inputs - pk(kx,ky)
c     inputs - xki,xka,yki,yka
c     inputs - xmi,xma,ymi,yma
c     inputs - Kx,Ky,Nx,Ny
c------------------------------------------------------------
c      output - px(Nx,Ny),Py(Nx,Ny)
	real px(nd,*),py(nd,*),pk(kd,*),P(ND,*)
	real PKK(50,50),C(4,4,50,50)

      c1=3.1415926/180.0      
      c2=111111.0

	dky=(yka-yki)/(ky-1)
	dkx=(xka-xki)/(kx-1)
      
	dy=(yma-ymi)/(ny-1)
	dx=(xma-xmi)/(nx-1)

C     SURROUNDING
	DO J=2,KY+1       
        DO I=2,KX+1
          PKK(I,J)=PK(I-1,J-1)
        END DO
      END DO
	DO J=2,KY+1       
        PKK(1,J)=2*PKK(2,J)-PKK(3,J)
        PKK(KX+2,J)=2*PKK(KX+1,J)-PKK(KX,J)
      END DO
	DO I=1,KX+2       
        PKK(I,1)=2*PKK(I,2)-PKK(I,3)
        PKK(I,KY+2)=2*PKK(I,KY+1)-PKK(I,KY)
      END DO
	CALL GETBICUBIC(KX+2,KY+2,50,PKK,C)



	do j=1,Ny
	  y=ymi+(j-1)*dy
	  j0=(y-yki)/dky+1
	  if (j0<1) j0=1
	  if (j0>ky-1) j0=ky-1
	  u=(y-(yki+(j0-1)*dky))/dky
	  
	  do i=1,Nx
	    x=xmi+(i-1)*dx
	    i0=(x-xki)/dkx+1
	    if (i0<1) I0=1
	    if (i0>kx-1) i0=kx-1
	    t=(x-(xki+(i0-1)*dkx))/dkx
          ay=0.
          a2=0.
          a1=0.
          DO K=4,1,-1
            ay=t*ay+((c(K,4,i0,j0)*u+c(k,3,i0,j0))*u+c(K,2,i0,j0))*u+
     1		  c(K,1,i0,j0)
            a2=t*a2+(3.*c(K,4,i0,j0)*u+2.*c(K,3,i0,j0))*u+c(K,2,i0,j0)
            a1=u*a1+(3.*c(4,K,i0,j0)*t+2.*c(3,K,i0,j0))*t+c(2,K,i0,j0)
          END DO
          a1=a1/dkx/c2/cos(c1*y)
          a2=a2/dky/c2

          p(i,j)=ay
	    px(i,j)=a1
	    py(i,j)=a2

	  end do
      END DO


	return
	end



	SUBROUTINE GETBICUBIC(NX,NY,ND,Z,C)
	REAL Z(ND,*), C(4,4,ND,*),Y(4),Y1(4),Y2(4),Y12(4),CC(4,4)
      D1=1
	D2=1
C     CALCULATE DERIVATIVE

      DO J=2,NY-2
	  DO I=2,NX-2
	     Y(1)=Z(I,J)
           Y(2)=Z(I+1,J)
		 Y(3)=Z(I+1,J+1)
	     Y(4)=Z(I,J+1)

		 Y1(1)=0.5*(Z(I+1,J)-Z(I-1,J)) 
		 Y1(4)=0.5*(Z(I+1,J+1)-Z(I-1,J+1)) 
		 Y1(2)=0.5*(Z(I+2,J)  -Z(I,J)) 
		 Y1(3)=0.5*(Z(I+2,J+1)-Z(I,J+1)) 
		 
		 Y2(1)=0.5*(Z(I,J+1)  -Z(I,J-1)) 
		 Y2(2)=0.5*(Z(I+1,J+1)-Z(I+1,J-1)) 
		 Y2(3)=0.5*(Z(I+1,J+2)-Z(I+1,J)) 
		 Y2(4)=0.5*(Z(I,J+2)-Z(I,J))
		 
		 Y12(1)=0.25*(Z(I+1,J+1)-Z(I+1,J-1)-Z(I-1,J+1)+Z(I-1,J-1)) 
		 Y12(2)=0.25*(Z(I+2,J+1)-Z(I+2,J-1)-Z(I,J+1)+Z(I,J-1)) 
		 Y12(3)=0.25*(Z(I+2,J+2)-Z(I+2,J)-Z(I,J+2)+Z(I,J)) 
		 Y12(4)=0.25*(Z(I+1,J+2)-Z(I+1,J)-Z(I-1,J+2)+Z(I-1,J)) 
           
		 CALL bcucof(y,y1,y2,y12,d1,d2,CC)
	     DO K=1,4
	     DO L=1,4
		   C(K,L,I-1,J-1)=CC(K,L)
	     END DO
	     END DO
        
	  END DO
	END DO  		  



      RETURN
	END

      SUBROUTINE bcucof(y,y1,y2,y12,d1,d2,c)
      REAL d1,d2,c(4,4),y(4),y1(4),y12(4),y2(4)
      INTEGER i,j,k,l
      REAL d1d2,xx,cl(16),wt(16,16),x(16)
      SAVE wt
      DATA wt/1,0,-3,2,4*0,-3,0,9,-6,2,0,-6,4,8*0,3,0,-9,6,-2,0,6,-4,10*
     *0,9,-6,2*0,-6,4,2*0,3,-2,6*0,-9,6,2*0,6,-4,4*0,1,0,-3,2,-2,0,6,-4,
     *1,0,-3,2,8*0,-1,0,3,-2,1,0,-3,2,10*0,-3,2,2*0,3,-2,6*0,3,-2,2*0,
     *-6,4,2*0,3,-2,0,1,-2,1,5*0,-3,6,-3,0,2,-4,2,9*0,3,-6,3,0,-2,4,-2,
     *10*0,-3,3,2*0,2,-2,2*0,-1,1,6*0,3,-3,2*0,-2,2,5*0,1,-2,1,0,-2,4,
     *-2,0,1,-2,1,9*0,-1,2,-1,0,1,-2,1,10*0,1,-1,2*0,-1,1,6*0,-1,1,2*0,
     *2,-2,2*0,-1,1/
      d1d2=d1*d2
      do 11 i=1,4
        x(i)=y(i)
        x(i+4)=y1(i)*d1
        x(i+8)=y2(i)*d2
        x(i+12)=y12(i)*d1d2
11    continue
      do 13 i=1,16
        xx=0.
        do 12 k=1,16
          xx=xx+wt(i,k)*x(k)
12      continue
        cl(i)=xx
13    continue
      l=0
      do 15 i=1,4
        do 14 j=1,4
          l=l+1
          c(i,j)=cl(l)
14      continue
15    continue
      return
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

      SUBROUTINE READDIMGR3(NX,NY,NZ,NAME)
	CHARACTER*20 NAME
      INTEGER NX,NY,NZ
	
      OPEN(1,FILE=NAME)
	READ(1,*)
	READ(1,*)NX,NY,NZ
	CLOSE(1)
	RETURN
	END 
      
      SUBROUTINE READGR3(NX,NY,NZ,XMI,XMA,YMI,YMA,ZMI,ZMA,NAME,Z)
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
      