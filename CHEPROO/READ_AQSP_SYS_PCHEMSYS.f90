SUBROUTINE READ_AQSP_SYS_PCHEMSYS &
  (ITEMP,IACT,ICONV,NAQT,NPRI,TC2,LABEL,NAAQT,MXTOT,MXLABEL,IseRROR)
  
implicit double precision (a-h,o-z)  
implicit integer (i-n)
LOGICAL, INTENT(OUT) :: IseRROR
!-------------------------------------------------------------------------
!
!>   $Description:
!
!>   $Arguments:
!
 
 
 
!-------------------------------------------------------------------------
!
!>   $Pre-cond:
!
!>   $Post-cond:
!
!>   $License:
!
!-------------------------------------------------------------------------
 
!-------------------------------------------------------------------------
!
!>   $code
!
  
!-------------------------------------------------------------------------
!
!>   $Description:
!
!>   $Arguments:
!
CHARACTER(LEN=*) LABEL(MXLABEL)
CHARACTER(LEN=*) NAAQT(MXTOT) 
 
!-------------------------------------------------------------------------
!
!>   $Pre-cond:
!
!>   $Post-cond:
!
!>   $License:
!
!-------------------------------------------------------------------------
 
!-------------------------------------------------------------------------
!
!>   $code
!
 
IseRROR=.FALse. 
 
  
 
!%      read control parameters of the chemical system
 
  READ (1,'(A)',ERR=9001) LABEL(1)      !title of the problem
 
!%      read aqueous species of the system
  READ (1,'(A)',ERR=9002) LABEL(2)      !---------------------------
 
  READ (1,'(A)',ERR=9002) LABEL(3)      !aqueous species
 
  READ (1,'(f10.0,2i5)',ERR=9002) TC2,IACT,ICONV        !temp. to calculate
 
  IF (TC2.EQ.25.0D0) THEN
    ITEMP=0
  ELse
    ITEMP=1
  ENDIF
 
 
        READ (1,'(A)',ERR=9003) LABEL(4)      !aqueous primary species
            
        I=0
130     CONTINUE
        I=I+1
        READ (1,*,ERR=9003) NAAQT(I)
        IF (NAAQT(I).EQ.'*') THEN
          NPRI=I-1
          GO TO 120
        ENDIF
        GO TO 130
120     CONTINUE

!%       READ AQUEOUS COMPLEXES
        READ (1,'(A)',ERR=9004) LABEL(5)             !aqueous complexes
           
        J=NPRI
132     CONTINUE
        J=J+1
        READ (1,*,ERR=9004) NAAQT(J)
        IF (NAAQT(J).EQ.'*') THEN
          NAQT=J-1
          GO TO 122
        END IF
        GO TO 132
122     CONTINUE

        RETURN

9001    WRITE(*,*) 'ERROR READING TITLE'
        IseRROR=.TRUE.
        RETURN
9002    WRITE(*,*) 'ERROR READING TEMPERATURE AND THERMODYNAMIC MODEL'
        IseRROR=.TRUE.
        RETURN
9003    WRITE(*,*) 'ERROR READING PRIMARY AQUEOUS species'
        IseRROR=.TRUE.
        RETURN
9004    WRITE(*,*) 'ERROR READING AQUEOUS COMPLEXES'
        IseRROR=.TRUE.
        RETURN
 
END SUBROUTINE