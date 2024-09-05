SUBROUTINE READ_ADS_SYS_PCHEMSYS(NADS,NMOD,NPRI,NAPRI,NAADS, &
                           NAMOD,NADSMOD,MXDIM,MXLABEL,iserror) 

IMPLICIT NONE
!******************************************************************
!>  INTERNAL VARIABLES:SCALARS
!>    NAADSDUM 
!>  EXTERNAL VARIABLES:SCALARS
!>    NADS         number of total adsorbed species
!>    NSTT         number of total primary species of adsorption
!>    NMOD         number of adsorption models                           
!>  EXTERNAL VARIABLES:ARRAYS
!>    NADSMOD         adsorption model name
!>    NST          number of primary species of adsorption for each
!>                 model
!>    NADSMOD      number of adsorbed species for each model
!>    NAADS        name of total adsorbed species
!>    LABELS
!>  EXTERNAL VARIABLES:WORKSPACE
!>    NAST         name of primary species for each model
!>    NAADSMOD     name of adsorbed species for each model
!********************************************************************          
        INTEGER MXLABEL,MXDIM,NADS,NMOD,NPRI
        CHARACTER(LEN=*) NAADS(MXDIM),NAMOD(MXDIM),NAPRI(MXDIM)
        CHARACTER(LEN=20)JADSDUM(MXDIM),LABELS(MXLABEL),NAADSDUM
        INTEGER NADSMOD(MXDIM)
        INTEGER I,K,J  
        LOGICAL ISPRI,iserror
        
        iserror=.false.  

!>       read surface complexes
        READ (1,'(A)',ERR=9001) LABELS(1)              !surface complexes
   
        NMOD=0
        NADS=0
        I=0 
 
210     READ (1,*,ERR=9002) NAADSDUM   
        IF (NAADSDUM.EQ.'*') GO TO 250
215     I=I+1
        NMOD=I
        JADSDUM(I)= NAADSDUM
        K=0
        J=0

!>       Read kinetics model
        IF (NAADSDUM.EQ.'EX+KIN') THEN
           JADSDUM(I)='EX'
        ENDIF
        IF (NAADSDUM.EQ.'TL+KIN') THEN
           JADSDUM(I)='TL'

        ENDIF
        IF (NAADSDUM.EQ.'DL+KIN') THEN
           JADSDUM(I)='DL'
          
	  ENDIF
        IF (NAADSDUM.EQ.'CC+KIN') THEN
           JADSDUM(I)='CC'
	   
        ENDIF
        IF (NAADSDUM.EQ.'EQ+KIN') THEN
           JADSDUM(I)='EQ'
	   
        ENDIF
        IF (NAADSDUM.EQ.'LA+KIN') THEN
           JADSDUM(I)='LA'
	    
        ENDIF
        IF (NAADSDUM.EQ.'FR+KIN') THEN
           JADSDUM(I)='FR'
           
	  ENDIF
        IF (NAADSDUM.EQ.'KD+KIN') THEN
           JADSDUM(I)='KD'
          
	  ENDIF

      
       NAMOD(I)=JADSDUM(I)
      
!>       Read adsorbate names + kinetic constand + kinetic order
       
       K=0
       ISPRI=.TRUE.
220    READ (1,*,ERR=9003) NAADSDUM
	   	
		
		  IF (NAADSDUM.EQ.'EX'.OR.NAADSDUM.EQ.'TL'.OR. &
           NAADSDUM.EQ.'DL'.OR.NAADSDUM.EQ.'CC'.OR. &
           NAADSDUM.EQ.'FR'.OR.NAADSDUM.EQ.'LA'.OR. &
           NAADSDUM.EQ.'EQ'.OR.NAADSDUM.EQ.'KD'.OR. &
           NAADSDUM.EQ.'EX+KIN'.OR.NAADSDUM.EQ.'TL+KIN'.OR. &
           NAADSDUM.EQ.'DL+KIN'.OR.NAADSDUM.EQ.'CC+KIN'.OR. &
           NAADSDUM.EQ.'FR+KIN'.OR.NAADSDUM.EQ.'LA+KIN'.OR. &
           NAADSDUM.EQ.'EQ+KIN'.OR.NAADSDUM.EQ.'KD+KIN') THEN

          GO TO 215
            
		  END IF
     
          IF (NAADSDUM .EQ. '*') GOTO 250
	      IF (ISPRI) THEN 
	       ISPRI=.FALse.
	       seLECT CAse (NAMOD(I))
	       CAse('EX','TL','DL','CC')
	        NPRI=NPRI+1
		    NAPRI(NPRI)=NAADSDUM
		   END seLECT 
	      ELse
	       NADS=NADS+1
	       K=K+1
		   NAADS(NADS)=NAADSDUM   
		   NADSMOD(I)=K 
		  END IF 
 
          GO TO 220

        
 
 
250    CONTINUE
       IF (NADS==0.AND.NMOD/=0) GOTO 9003


       RETURN

9001   WRITE(*,*) 'ERROR READING HEADING SURF. COMPL. OF THE SYSTEM'
       iserror=.true.
       return

9002   WRITE(*,*) 'ERROR READING MODEL OF ADSORPTION OF THE SYSTEM'
       iserror=.true.
       return

9003   WRITE(*,*) 'ERROR READING NAMES OF SURF, COMPL. OF THE SYSTEM'
       iserror=.true.
       return

END SUBROUTINE 