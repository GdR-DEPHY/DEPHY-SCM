!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ###############################################################################
SUBROUTINE COUPLING_TSZ0_n (DTCO, UG, U, USS, IM, DTZ, NDST, SLT,BLOWSNW, HPROGRAM, HCOUPLING,   &
                            PTSTEP, KYEAR, KMONTH, KDAY, PTIME, KI, KSV, KSW, PTSUN,&
                            PZENITH, PZENITH2, PAZIM, PZREF, PUREF, PZS, PU, PV,    &
                            PQA, PTA, PRHOA, PSV, PCO2, HSV, PRAIN, PSNOW, PLW,     &
                            PDIR_SW, PSCA_SW, PSW_BANDS, PPS, PPA, PSFTQ, PSFTH,    &
                            PSFTS, PSFCO2, PSFU, PSFV, PTRAD, PDIR_ALB, PSCA_ALB,   &
                            PEMIS, PTSURF, PZ0, PZ0H, PQSURF, PPEW_A_COEF,          &
                            PPEW_B_COEF, PPET_A_COEF, PPEQ_A_COEF, PPET_B_COEF,     &
                            PPEQ_B_COEF, HTEST      )  
!     ###############################################################################
!
!!****  *COUPLING_TSZ0_n * - Call of fluxes from vegetation scheme ISBA but 
!!        without temporal evolution of the soil/vegetation.
!!
!!    PURPOSE
!!    -------
!
!!**  METHOD
!!    ------
!!
!!    REFERENCE
!!    ---------
!!      
!!
!!    AUTHOR
!!    ------
!!     V. Masson 
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2004
!!      Modified    09/2012 : J. Escobar , SIZE(PTA) not allowed without-interface , replace by KI
!!      B. Decharme 04/2013 new coupling variables
!!      P. LeMoigne 12/2014 bug in "implicit" coefficients 
!!------------------------------------------------------------------
!
USE MODD_ISBA_n, ONLY : ISBA_P_t, ISBA_PE_t
USE MODD_SURFEX_n, ONLY : ISBA_MODEL_t
!
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_SURF_ATM_GRID_n, ONLY : SURF_ATM_GRID_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
USE MODD_SSO_n, ONLY : SSO_t
USE MODD_DATA_TSZ0_n, ONLY : DATA_TSZ0_t
USE MODD_DATA_ISBA_n, ONLY : DATA_ISBA_t
USE MODD_DST_n, ONLY : DST_NP_t
USE MODD_SLT_n, ONLY : SLT_t
USE MODD_BLOWSNW_n, ONLY : BLOWSNW_t
USE MODD_BLANK, ONLY : LDUMMY1
!
!
USE MODD_SURF_PAR, ONLY : XUNDEF
USE MODD_CSTS,   ONLY : XP00, XRD, XCPD, XG, XKARMAN
!
USE MODI_TSZ0
USE MODI_COUPLING_ISBA_OROGRAPHY_n
USE MODI_WIND_THRESHOLD

! 
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
TYPE(ISBA_MODEL_t), INTENT(INOUT) :: IM
!
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
TYPE(SURF_ATM_GRID_t), INTENT(INOUT) :: UG
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
TYPE(SSO_t), INTENT(INOUT) :: USS
TYPE(DATA_TSZ0_t), INTENT(INOUT) :: DTZ
TYPE(DST_NP_t), INTENT(INOUT) :: NDST
TYPE(SLT_t), INTENT(INOUT) :: SLT
TYPE(BLOWSNW_t), INTENT(INOUT) :: BLOWSNW
!
 CHARACTER(LEN=6),    INTENT(IN)  :: HPROGRAM  ! program calling surf. schemes
 CHARACTER(LEN=1),    INTENT(IN)  :: HCOUPLING ! type of coupling
                                              ! 'E' : explicit
                                              ! 'I' : implicit
INTEGER,             INTENT(IN)  :: KYEAR     ! current year (UTC)
INTEGER,             INTENT(IN)  :: KMONTH    ! current month (UTC)
INTEGER,             INTENT(IN)  :: KDAY      ! current day (UTC)
REAL,                INTENT(IN)  :: PTIME     ! current time since midnight (UTC, s)
INTEGER,             INTENT(IN)  :: KI        ! number of points
INTEGER,             INTENT(IN)  :: KSV       ! number of scalars
INTEGER,             INTENT(IN)  :: KSW       ! number of short-wave spectral bands
REAL, DIMENSION(KI), INTENT(IN)  :: PTSUN     ! solar time                    (s from midnight)
REAL,                INTENT(IN)  :: PTSTEP    ! atmospheric time-step                 (s)
REAL, DIMENSION(KI), INTENT(IN)  :: PZREF     ! height of T,q forcing                 (m)
REAL, DIMENSION(KI), INTENT(IN)  :: PUREF     ! height of wind forcing                (m)
!
REAL, DIMENSION(KI), INTENT(IN)  :: PTA       ! air temperature forcing               (K)
REAL, DIMENSION(KI), INTENT(IN)  :: PQA       ! air humidity forcing                  (kg/m3)
REAL, DIMENSION(KI), INTENT(IN)  :: PRHOA     ! air density                           (kg/m3)
REAL, DIMENSION(KI,KSV),INTENT(IN) :: PSV     ! scalar variables
!                                             ! chemistry:   first char. in HSV: '#'  (molecule/m3)
!                                             !
 CHARACTER(LEN=6), DIMENSION(KSV),INTENT(IN):: HSV  ! name of all scalar variables
REAL, DIMENSION(KI), INTENT(IN)  :: PU        ! zonal wind                            (m/s)
REAL, DIMENSION(KI), INTENT(IN)  :: PV        ! meridian wind                         (m/s)
REAL, DIMENSION(KI,KSW),INTENT(IN) :: PDIR_SW ! direct  solar radiation (on horizontal surf.)
!                                             !                                       (W/m2)
REAL, DIMENSION(KI,KSW),INTENT(IN) :: PSCA_SW ! diffuse solar radiation (on horizontal surf.)
!                                             !                                       (W/m2)
REAL, DIMENSION(KSW),INTENT(IN)  :: PSW_BANDS ! mean wavelength of each shortwave band (m)
REAL, DIMENSION(KI), INTENT(IN)  :: PZENITH   ! zenithal angle at t      (radian from the vertical)
REAL, DIMENSION(KI), INTENT(IN)  :: PZENITH2  ! zenithal angle at t+1    (radian from the vertical)
REAL, DIMENSION(KI), INTENT(IN)  :: PAZIM     ! azimuthal angle      (radian from North, clockwise)
REAL, DIMENSION(KI), INTENT(IN)  :: PLW       ! longwave radiation (on horizontal surf.)
!                                             !                                       (W/m2)
REAL, DIMENSION(KI), INTENT(IN)  :: PPS       ! pressure at atmospheric model surface (Pa)
REAL, DIMENSION(KI), INTENT(IN)  :: PPA       ! pressure at forcing level             (Pa)
REAL, DIMENSION(KI), INTENT(IN)  :: PZS       ! atmospheric model orography           (m)
REAL, DIMENSION(KI), INTENT(IN)  :: PCO2      ! CO2 concentration in the air          (kg/m3)
REAL, DIMENSION(KI), INTENT(IN)  :: PSNOW     ! snow precipitation                    (kg/m2/s)
REAL, DIMENSION(KI), INTENT(IN)  :: PRAIN     ! liquid precipitation                  (kg/m2/s)
!
!
REAL, DIMENSION(KI), INTENT(OUT) :: PSFTH     ! flux of heat                          (W/m2)
REAL, DIMENSION(KI), INTENT(OUT) :: PSFTQ     ! flux of water vapor                   (kg/m2/s)
REAL, DIMENSION(KI), INTENT(OUT) :: PSFU      ! zonal momentum flux                   (Pa)
REAL, DIMENSION(KI), INTENT(OUT) :: PSFV      ! meridian momentum flux                (Pa)
REAL, DIMENSION(KI), INTENT(OUT) :: PSFCO2    ! flux of CO2                           (m/s*kg_CO2/kg_air)
REAL, DIMENSION(KI,KSV),INTENT(OUT):: PSFTS   ! flux of scalar var.                   (kg/m2/s)
!
REAL, DIMENSION(KI), INTENT(OUT) :: PTRAD     ! radiative temperature                 (K)
REAL, DIMENSION(KI,KSW),INTENT(OUT):: PDIR_ALB! direct albedo for each spectral band  (-)
REAL, DIMENSION(KI,KSW),INTENT(OUT):: PSCA_ALB! diffuse albedo for each spectral band (-)
REAL, DIMENSION(KI), INTENT(OUT) :: PEMIS     ! emissivity                            (-)
!
REAL, DIMENSION(KI), INTENT(OUT) :: PTSURF    ! surface effective temperature         (K)
REAL, DIMENSION(KI), INTENT(OUT) :: PZ0       ! roughness length for momentum         (m)
REAL, DIMENSION(KI), INTENT(OUT) :: PZ0H      ! roughness length for heat             (m)
REAL, DIMENSION(KI), INTENT(OUT) :: PQSURF    ! specific humidity at surface          (kg/kg)
!
REAL, DIMENSION(KI), INTENT(IN) :: PPEW_A_COEF! implicit coefficients
REAL, DIMENSION(KI), INTENT(IN) :: PPEW_B_COEF! needed if HCOUPLING='I'
REAL, DIMENSION(KI), INTENT(IN) :: PPET_A_COEF
REAL, DIMENSION(KI), INTENT(IN) :: PPEQ_A_COEF
REAL, DIMENSION(KI), INTENT(IN) :: PPET_B_COEF
REAL, DIMENSION(KI), INTENT(IN) :: PPEQ_B_COEF
 CHARACTER(LEN=2),    INTENT(IN) :: HTEST ! must be equal to 'OK'
!
!*      0.2    declarations of local variables
!
!
TYPE(ISBA_P_t), POINTER :: PK
TYPE(ISBA_PE_t), POINTER :: PEK
!
REAL, DIMENSION(KI,IM%O%NGROUND_LAYER,IM%O%NPATCH) :: ZTG   ! soil temperature
REAL, DIMENSION(KI,IM%O%NGROUND_LAYER,IM%O%NPATCH) :: ZWG   ! soil water content
REAL, DIMENSION(KI,IM%O%NGROUND_LAYER,IM%O%NPATCH) :: ZWGI  ! soil ice content
REAL, DIMENSION(KI,IM%O%NPATCH) :: ZWR   ! interception reservoir
REAL, DIMENSION(KI,IM%O%NPATCH) :: ZRESA ! aerodynamical resistance
REAL, DIMENSION(KI,IM%NPE%AL(1)%TSNOW%NLAYER,IM%O%NPATCH) :: ZWSNOW! snow reservoir
REAL, DIMENSION(KI,IM%NPE%AL(1)%TSNOW%NLAYER,IM%O%NPATCH) :: ZRHOSN! snow density
REAL, DIMENSION(KI,IM%NPE%AL(1)%TSNOW%NLAYER,IM%O%NPATCH) :: ZHEASN! snow heat content
REAL, DIMENSION(KI,IM%O%NPATCH) :: ZALBSN! snow albedo
REAL, DIMENSION(KI,IM%O%NPATCH) :: ZEMISN! snow emissivity
!
REAL, DIMENSION(KI)     :: ZPEW_A_COEF ! implicit coefficients
REAL, DIMENSION(KI)     :: ZPEW_B_COEF ! needed if HCOUPLING='I'
REAL, DIMENSION(KI)     :: ZPET_A_COEF
REAL, DIMENSION(KI)     :: ZPEQ_A_COEF
REAL, DIMENSION(KI)     :: ZPET_B_COEF
REAL, DIMENSION(KI)     :: ZPEQ_B_COEF
INTEGER :: JP
REAL, DIMENSION(KI) :: ZUSTAR,ZL,ZTSTAR
REAL :: ZBM,ZBH
INTEGER  :: ZI
REAL, DIMENSION(KI)     :: ZVMOD! wind modulus with minimum threshold
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('COUPLING_TSZ0_N',0,ZHOOK_HANDLE)
!
!*      1.     Specified evolution of ISBA prognostic variables
!              ------------------------------------------------
!
DO JP = 1,IM%O%NPATCH
  CALL TSZ0(DTZ, PTIME, PTSTEP, IM%NK%AL(JP), IM%NPE%AL(JP))
ENDDO
!
!
!*      2.     Saves the prognostic variables
!              ------------------------------
!
DO JP = 1,IM%O%NPATCH
  PK => IM%NP%AL(JP)
  PEK => IM%NPE%AL(JP)
  !
  ZTG  (1:PK%NSIZE_P,:,JP) = PEK%XTG        (:,:)
  ZWG  (1:PK%NSIZE_P,:,JP) = PEK%XWG        (:,:)
  ZWGI (1:PK%NSIZE_P,:,JP) = PEK%XWGI       (:,:)
  ZWR  (1:PK%NSIZE_P,JP)   = PEK%XWR        (:)
  ZRESA(1:PK%NSIZE_P,JP)   = PEK%XRESA      (:)
  ZWSNOW(1:PK%NSIZE_P,:,JP)= PEK%TSNOW%WSNOW(:,:)
  ZRHOSN(1:PK%NSIZE_P,:,JP)= PEK%TSNOW%RHO  (:,:)
  ZALBSN(1:PK%NSIZE_P,JP)  = PEK%TSNOW%ALB  (:)
  IF (PEK%TSNOW%SCHEME=='3-L' .OR. PEK%TSNOW%SCHEME=='CRO') THEN
    ZHEASN(1:PK%NSIZE_P,:,JP)= PEK%TSNOW%HEAT (:,:)
    ZEMISN(1:PK%NSIZE_P,JP)  = PEK%TSNOW%EMIS (:)
  END IF
ENDDO
!
!
!*      3.     Call to surface scheme
!              ----------------------
!
 CALL COUPLING_ISBA_OROGRAPHY_n(DTCO, UG, U, USS, IM%SB, IM%NAG, IM%CHI, IM%NCHI, IM%MGN, IM%MSF,IM%DTV, IM%ID, &
                                IM%NGB, IM%GB, IM%ISS, IM%NISS, IM%G, IM%NG, IM%O, IM%S, IM%K, IM%NK, &
                                IM%NP, IM%NPE, NDST, SLT,BLOWSNW, HPROGRAM, 'E', 0.001, KYEAR,   &
                                KMONTH, KDAY, PTIME,  KI, KSV, KSW, PTSUN, PZENITH,       &
                                PZENITH2, PAZIM, PZREF, PUREF, PZS, PU, PV, PQA, PTA,     &
                                PRHOA, PSV, PCO2, HSV, PRAIN, PSNOW, PLW, PDIR_SW,        &
                                PSCA_SW, PSW_BANDS, PPS, PPA, PSFTQ, PSFTH, PSFTS, PSFCO2,&
                                PSFU, PSFV, PTRAD, PDIR_ALB, PSCA_ALB, PEMIS, PTSURF, PZ0,&
                                PZ0H, PQSURF, PPEW_A_COEF, PPEW_B_COEF, PPET_A_COEF,      &
                                PPEQ_A_COEF, PPET_B_COEF, PPEQ_B_COEF, 'OK'  )  
!

!ATTENTION : variation de Z0 dans le temps non implemente
IF(LDUMMY1) THEN
 !Coefficient de Beare et al. et Cuxart et al. (2006)
 ZBM=4.8
 ZBH=7.8
 ZL(:)=-9999
 ZVMOD= WIND_THRESHOLD(SQRT(PU(:)**2. + PV(:)**2.),PUREF)
!
 DO ZI=1,50
   ZUSTAR(:) = ZVMOD(:)/(LOG(PUREF(:)/IM%DTV%XPAR_Z0(1,1,1))/XKARMAN &
              +ZBM/(ZL(:)*XKARMAN)*(PUREF(:)-IM%DTV%XPAR_Z0(1,1,1)))

   ZTSTAR(:) = (PTA(:)*((XP00/PPA)**(XRD/XCPD))-ZTG(:,1,1)*((XP00/PPS)**(XRD/XCPD))) &
         /(LOG(PZREF(:)/IM%DTV%XPAR_Z0(1,1,1))/XKARMAN+ZBH/(ZL(:)*XKARMAN)*(PZREF(:)-IM%DTV%XPAR_Z0(1,1,1)))

   WHERE(ABS(ZTSTAR)<1.E-10)
     ZTSTAR(:) = 1.E-10
   END WHERE
  
   WHERE(ABS(ZUSTAR)<1.E-10)
     ZUSTAR(:) = 1.E-10
   END WHERE

   ZL(:)=(ZUSTAR(:)**2)/(XG/(PTA(:)*((XP00/PPA)**(XRD/XCPD)))*XKARMAN*ZTSTAR(:))
 END DO
!
 PSFU(:)  = -PU(:)/ZVMOD(:)*(ZUSTAR(:)**2)*PRHOA
 PSFV(:)  = -PV(:)/ZVMOD(:)*(ZUSTAR(:)**2)*PRHOA
 PSFTH(:) = -ZTSTAR(:)*ZUSTAR(:)*PRHOA*XCPD
END IF
!
!*      4.     Removes temporal evolution of ISBA variables
!              --------------------------------------------
!
!
DO JP = 1,IM%O%NPATCH
  PK => IM%NP%AL(JP)
  PEK => IM%NPE%AL(JP)
  !
  PEK%XTG        (:,:) = ZTG  (1:PK%NSIZE_P,:,JP)
  PEK%XWG        (:,:) = ZWG  (1:PK%NSIZE_P,:,JP)
  PEK%XWGI       (:,:) = ZWGI (1:PK%NSIZE_P,:,JP)
  PEK%XWR        (:)   = ZWR  (1:PK%NSIZE_P,JP)  
  PEK%XRESA      (:)   = ZRESA(1:PK%NSIZE_P,JP) 
  PEK%TSNOW%WSNOW(:,:) = ZWSNOW(1:PK%NSIZE_P,:,JP)
  PEK%TSNOW%RHO  (:,:) = ZRHOSN(1:PK%NSIZE_P,:,JP)
  PEK%TSNOW%ALB  (:)   = ZALBSN(1:PK%NSIZE_P,JP)
  IF (PEK%TSNOW%SCHEME=='3-L' .OR. PEK%TSNOW%SCHEME=='CRO') THEN
    PEK%TSNOW%HEAT (:,:) = ZHEASN(1:PK%NSIZE_P,:,JP)
    PEK%TSNOW%EMIS (:)   = ZEMISN(1:PK%NSIZE_P,JP) 
  END IF
ENDDO
!
IF (LHOOK) CALL DR_HOOK('COUPLING_TSZ0_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE COUPLING_TSZ0_n
