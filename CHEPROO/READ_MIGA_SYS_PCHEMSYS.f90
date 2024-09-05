SUBROUTINE READ_MIGA_SYS_PCHEMSYS &
  (NGAS,NMIN, NMINEQ, NMINKIN, NAQT, NAAQT, &
    IMEQ, LABELM, NAGAS,NAMIN,MXLABEL,MXSP,iserror)
 
implicit double precision (a-h,o-z)  
implicit integer (i-n)
logical, intent(out)   :: iserror
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


 
   CHARACTER(LEN=*) LABELM(MXLABEL)
   CHARACTER(LEN=*) NAMIN(MXSP),NAGAS(MXSP), NAAQT(MXSP)
   DIMENSION IMEQ(MXSP)
 
iserror=.false. 
!FBP ini
!%        read minerals
         READ (1,'(A)',ERR=9001) LABELM(1)             !minerals
              
         M=0
         NMINEQ=0
         NMINKIN=0
140      CONTINUE
         M=M+1
         READ (1,*,ERR=9002) NAMIN(M),IMEQ(M)
         IF (NAMIN(M).EQ.'*') GOTO 150
         IF (IMEQ(M).EQ.0) THEN 
            NMINEQ=NMINEQ+1
         ELse
           NMINKIN=NMINKIN+1
         ENDIF
         GO TO 140
150      CONTINUE
         NMIN=NMINEQ+NMINKIN


!%       read kinetic data from kinetics.dat file 
!%        IF (NMINKIN.GT.0) CALL DATAKIN &
!%        (NAQT,  NMIN, NMINKIN, NAAQT, IMEQ, NAMIN,       &
!%           EA,  NKIN, THRESH, DINS,   FORA,    RK, PCAT, &
!%         NCAT, NACAT, filebase, kin_filename)  !FBP added kin_filename

!%       read gases
        READ (1,'(A)',ERR=9003) LABELM(3) 		!gases
        K=0
170     CONTINUE
        K=K+1
        READ (1,*,ERR=9004) NAGAS(K)

        IF (NAGAS(K).EQ.'*') THEN
          NGAS=K-1

          GO TO 190
        END IF
        GO TO 170
190     CONTINUE

        RETURN

9001    WRITE(*,*) 'ERROR READING HEADING OF MINERALS OF THE SYSTEM'
        iserror=.true.
        return
9002    WRITE(*,*) 'ERROR READING NAMES OF THE MINERALS OF THE SYSTEM'
        iserror=.true.
        return
9003    WRITE(*,*) 'ERROR READING HEADING OF GAseS OF THE SYSTEM'
        iserror=.true.
        return
9004    WRITE(*,*) 'ERROR READING NAMES OF THE GAseS OF THE SYSTEM'
        iserror=.true.
        return
 
END SUBROUTINE