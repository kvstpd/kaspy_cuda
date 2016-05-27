#ifndef	_FORTRAN_VARS_H_
#define	_FORTRAN_VARS_H_

//PARAMETER (IM=322,JM=442)
//PARAMETER (IMM1=IM-1,JMM1=JM-1)
//PARAMETER (IMM2=IM-2,JMM2=JM-2)
//PARAMETER (LIJ=IM*JM)

#define	F_DATA_WIDTH 322
#define	F_DATA_HEIGHT 442

#define	F_DATA_SIZE F_DATA_WIDTH * F_DATA_HEIGHT

#pragma pack(8)


//COMMON/BLKCON/
//1          GRAV,padding,timeh,timeh6,tide_l,tide_l0,
//2          IINT,IPRINT,DTE,DTI,TPRNI,UMOL,xmi,xma,ymi,yma
//real*4 padding
//real*8 timeh,timeh6,tide_l,TIDE_L0,dti,dte

typedef struct{
    float dht;
    float grav;
    double timeh;
    double timeh6;
    double tide_l;
    double tide_l0;
    int iint;
    int iprint;
    double dte;
    double dti;
    float tprni;
    float umol;
    float xmi;
    float xma;
    float ymi;
    float yma;
} fortran_common_vars;


/*COMMON/BLK2D/arrays_marker,H(IM,JM),DX(JM),DY(JM),D(IM,JM),
 1     ART(JM),ARU(JM),ARV(JM),CBC(IM,JM),
 3     DUM(IM,JM),DVM(IM,JM),FSM(IM,JM),COR(JM),
 4     WUSURF(IM,JM),WVSURF(IM,JM),WUBOT(IM,JM),WVBOT(IM,JM),
 5     TPS(IM,JM),AAM2D,padding,
 6     UAF(IM,JM),UA(IM,JM),UAB(IM,JM),VAF(IM,JM),VA(IM,JM),
 7     VAB(IM,JM),ELF(IM,JM),EL(IM,JM),ELB(IM,JM),
 8     FLUXUA(IM,JM),FLUXVA(IM,JM)*/


/*extern "C" {
 extern fortran_common_block blkcon_;
 }*/

/// #pragma pack(8)

typedef struct{
    double marker;
    float h[F_DATA_HEIGHT][F_DATA_WIDTH];
    float dx[F_DATA_HEIGHT];
    float dy[F_DATA_HEIGHT];
    float d[F_DATA_HEIGHT][F_DATA_WIDTH];
    float art[F_DATA_HEIGHT];
    float aru[F_DATA_HEIGHT];
    float arv[F_DATA_HEIGHT];
    float cbc[F_DATA_HEIGHT][F_DATA_WIDTH];
    float dum[F_DATA_HEIGHT][F_DATA_WIDTH];
    float dvm[F_DATA_HEIGHT][F_DATA_WIDTH];
    float fsm[F_DATA_HEIGHT][F_DATA_WIDTH];
    float cor[F_DATA_HEIGHT];
    float wusurf[F_DATA_HEIGHT][F_DATA_WIDTH];
    float wvsurf[F_DATA_HEIGHT][F_DATA_WIDTH];
    float wubot[F_DATA_HEIGHT][F_DATA_WIDTH];
    float wvbot[F_DATA_HEIGHT][F_DATA_WIDTH];
    float tps[F_DATA_HEIGHT][F_DATA_WIDTH];
    float aam2d;
    float padding;
    float uaf[F_DATA_HEIGHT][F_DATA_WIDTH];
    float ua[F_DATA_HEIGHT][F_DATA_WIDTH];
    float uab[F_DATA_HEIGHT][F_DATA_WIDTH];
    float vaf[F_DATA_HEIGHT][F_DATA_WIDTH];
    float va[F_DATA_HEIGHT][F_DATA_WIDTH];
    float vab[F_DATA_HEIGHT][F_DATA_WIDTH];
    float elf[F_DATA_HEIGHT][F_DATA_WIDTH];
    float el[F_DATA_HEIGHT][F_DATA_WIDTH];
    float elb[F_DATA_HEIGHT][F_DATA_WIDTH];
    float fluxua[F_DATA_HEIGHT][F_DATA_WIDTH];
    float fluxva[F_DATA_HEIGHT][F_DATA_WIDTH];
	float advua[F_DATA_HEIGHT][F_DATA_WIDTH];
	float advva[F_DATA_HEIGHT][F_DATA_WIDTH];
	float surf[F_DATA_HEIGHT][F_DATA_WIDTH];

	
    double end_marker;
} fortran_common_arrays;



/*COMMON/F_FLOATS/
3     ff_marker,fxf(im,jm),fyf(im,jm),fxb(im,jm),fyb(im,jm),
4     FF(IM,JM),FB(IM,JM),
7     ffu(im,jm),ffv(im,jm),fbu(im,jm),fbv(im,jm),ff_end*/


typedef struct{
    double marker;
    float fxf[F_DATA_HEIGHT][F_DATA_WIDTH];
    float fyf[F_DATA_HEIGHT][F_DATA_WIDTH];
    float fxb[F_DATA_HEIGHT][F_DATA_WIDTH];
    float fyb[F_DATA_HEIGHT][F_DATA_WIDTH];
    float ff[F_DATA_HEIGHT][F_DATA_WIDTH];
    float fb[F_DATA_HEIGHT][F_DATA_WIDTH];
    float ffu[F_DATA_HEIGHT][F_DATA_WIDTH];
    float ffv[F_DATA_HEIGHT][F_DATA_WIDTH];
    float fbu[F_DATA_HEIGHT][F_DATA_WIDTH];
    float fbv[F_DATA_HEIGHT][F_DATA_WIDTH];
    double end_marker;
} fortran_ffloats;


//       COMMON/F_WIND/kx,ky,kt,kxu,kyu,ktu,kxv,kyv,ktv,
//1  XKI,XKA,YKI,YKA,XKUI,XKUA,YKUI,YKUA,
//2  XKVI,XKVA,YKVI,YKVA

typedef struct{
    int kx;
    int ky;
    int kt;
    
    int kxu;
    int kyu;
    int ktu;
    
    int kxv;
    int kyv;
    int ktv;   
    
    float xki;
    float xka;
    float yki;
    float yka;
    
    float xkui;
    float xkua;
    float ykui;
    float ykua;
    
    float xkvi;
    float xkva;
    float ykvi;
    float ykva;
    
} fortran_wind_data;


#endif /* _FORTRAN_VARS_H_ */

