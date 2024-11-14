      module m_constants_cheproo
**************************************************************
**************************************************************
*
*    Constants used in CHEPROO 
*
**************************************************************
**************************************************************
	real*8, parameter              ::
     . avogadro=6.022142d23,     ! Avogadro number
     . r0=1.0d-8,                ! minimun radio start the nucleation in minerals 
     . kgwmol=0.01801534d0,      ! Kilogram of water per mol of water 
     . secondyear=31536000.0d0,  ! Seconds per year
     . secondday=86400.0d0,      ! Seconds per day
     . secondhour=3600.0d0,      ! Seconds per hour 
     . rgas=8.3145d0,            ! Universal constant of gases J/mol K 
     . faraday=96490.0d0,        ! Faraday constant 
     . epsiz=8.854d-1,          ! 
     . cwater=55.50837d0,        ! Water concentration per water kilogram 
     . mpaatm=0.0980665d0,       ! Megapascal per atmosphere
     . melectron=9.109534d-31,   ! Electron mass [kg]
     . kgrgr=1.0d-3,             ! Kgr per gr
     . m3cm3=1.0d-6 ,            ! m3 per cm3 
     . pi=3.14159265359d0        ! Pi number  
**************************************************************
**************************************************************
**************************************************************
**************************************************************
* CODES USED IN CHEPROO
* PHASE CODES

      integer,parameter                       ::
     . ip_liquid_phase  =101,
     . ip_gas_phase     =102
     
* OTHER CODES

      integer,parameter   ::
     . ip_nameLenght=20

     
     
    	end module m_constants_cheproo