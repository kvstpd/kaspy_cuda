      PARAMETER (IM=322,JM=442)
      PARAMETER (IMM1=IM-1,JMM1=JM-1)
      PARAMETER (IMM2=IM-2,JMM2=JM-2)
      PARAMETER (LIJ=IM*JM)
      COMMON/BLKCON/
     1          dht,GRAV,timeh,timeh6,tide_l,tide_l0,
     2          IINT,IPRINT,DTE,DTI,TPRNI,UMOL,xmi,xma,ymi,yma
C---------------- 2-D ARRAYS --------------------------------------
      COMMON/BLK2D/arrays_marker,H(IM,JM),DX(JM),DY(JM),D(IM,JM),
     1     ART(JM),ARU(JM),ARV(JM),CBC(IM,JM),
     3     DUM(IM,JM),DVM(IM,JM),FSM(IM,JM),COR(JM),
     4     WUSURF(IM,JM),WVSURF(IM,JM),WUBOT(IM,JM),WVBOT(IM,JM),
     5     TPS(IM,JM),AAM2D,padding,
     6     UAF(IM,JM),UA(IM,JM),UAB(IM,JM),VAF(IM,JM),VA(IM,JM),
     7     VAB(IM,JM),ELF(IM,JM),EL(IM,JM),ELB(IM,JM),
     8     FLUXUA(IM,JM),FLUXVA(IM,JM),arrays_end_marker


       real*4 dht,padding
       real*8 arrays_marker,timeh,timeh6,tide_l,TIDE_L0,dti,dte
       real*8 arrays_end_marker
       