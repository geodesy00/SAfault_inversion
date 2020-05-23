	PROGRAM SABurchfault
C
C Program to perform simulated annealing to estimate parameters of the Burch fault in Wichita Uplift region 
The program performs an estimation based on EGM2008 "observations" of the
C field curvature magnitude,Gc=sqrt((Gyy-Gxx)^2+(2Gxy)^2) along a profile that is
C orthogonal to the general strike direction of the fault.
C It is assumed that the observation profile is at zero elevation (Z = 0). 
C The parameters to be estimated include fault location, dip angle, density
C contrasts, and depths, any of which may be fixed (not estimated).
C The simulated annealing uses dynamic step size changes for each of the parameters,
C as well as hard boundaries for each the parameters.
C Also a Markov chain is generated for each parameter independently. 

	
	USE MSIMSL

	IMPLICIT REAL*8 (A-H,O-Z)

C NOBS is the number of observation; M is the maximum number of parameters
	PARAMETER (NOBS=58,M=6)	

C OBS is the vector of observations, X is the vector of observation locations
C along the profile, GC is a vector for the computed model values
	DIMENSION OBS(NOBS),X(NOBS),GC(NOBS)

C PARM is the vector of parameters, 
C XL  and XU are the defined lower and upper bounds of the parameters; PARM1 is a work vector 
C The parameters are ordered as: location, dip angle, 2 depths, 2 density  contrasts 
 	DIMENSION PARM(M),PARMOPT(M),XL(M),XU(M),PARM1(M)

C DEL holds the parameter increments; DEL0 is a work vector; MJ is a vector of counters	
	DIMENSION DEL(M),DEL0(M),MJ(M)

C EPS1, PHI1 are used to control the termination of the annealing simulation	
	DIMENSION EPS1(4),PHI1(3)
C

C Define logical flags
	LOGICAL KFLAG,MFLAG
C
C Define COMMON areas
C MODEL is a common area for model constants, containing lower boundaries of density
C contrasts, 2 times Newton's gravitational constant times conversion to Eotvos
	COMMON /MODEL/Z11,Z12,G2E
C
C MET is a common area for the Metropolis algorithm, containing current optimum cost
C function value, parameters for the Markov chain, flag for output
	COMMON /MET/PHIOPT,L,N,MFLAG
C
C ISA is a common area for the Metropolis algorithm, containing the initial seed
C for the random number generator.
	COMMON /ISA/ISEED
C
C Define some constants (G is Newton's gravitational constant in mks units)
	DATA PI/3.141592653589793D0/,G/6.674D-11/,DUMMY/0.D0/
C
C Initialize seed for random number generator
	ISEED=54321
	CALL RNSET(ISEED)
C
C Define some global parameters
	DTR=PI/180.D0
C
C Define some fixed parameters of the model
C Two times Newton's gravitational constant time conversion to Eotvos
	G2E=2.D0*G*1.D9
C
C Depths (negative elevations) of lower boundaries of density contrasts
	Z11=-8300.D0
	Z12=-13000.D0
C
C Define initial parameters and bounds
C locations [m]
C	XL(1)=0.D0; PARM(1)=60000.D0; XU(1)=100000.D0  
	XL(1)=0.D0; PARM(1)=60000.D0; XU(1)=200000.D0
C
C dip angles [deg -> radian]
C	XL(2)=90.D0*DTR; PARM(2)=150.D0*DTR; XU(2)=170.D0*DTR
	XL(2)=10.D0*DTR; PARM(2)=150.D0*DTR; XU(2)=170.D0*DTR
C
C depths of upper boundaries [m]

	XL(4)=-13000.D0; PARM(4)=-2000.D0; XU(4)=-500.D0
	XL(3)=-8300.D0; PARM(3)=-3500.D0; XU(3)=-500.D0

C
C density contrasts [kg/m^3]
	
	XL(5)=50.D0; PARM(5)=200.D0; XU(5)=1000.D0
	XL(6)=50.D0; PARM(6)=400.D0; XU(6)=1000.D0

C
C Change initial values of parameters to be estimated
	FACTOR=0.8D0
C	FACTOR=0.01D0
C	FACTOR=0.D0
	IF (FACTOR .NE. 0.D0) THEN
	  PARM(1)=XL(1)+FACTOR*(XU(1)-XL(1))
	  PARM(2)=XL(2)+FACTOR*(XU(2)-XL(2))
	  PARM(3)=XL(3)+FACTOR*(XU(3)-XL(3))
	  PARM(4)=XL(4)+FACTOR*(XU(4)-XL(4))

	END IF
C
C Define initial step sizes of parameters
C (set any to zero that should not be estimated)
C location [meters]
	DEL(1)=100.D0
C
C dip angle [degrees -> radians]
	DEL(2)=1.D0*DTR
C
C depths (upper boundaries) [meters]
	DEL(3)=100.D0
	DEL(4)=100.D0
C	DEL(3)=0.D0
C	DEL(4)=0.D0
C
C density contrasts [kg/m^3]
	DEL(5)=0.D0
	DEL(6)=0.D0

C
C Define initial temperature
	T=10000.D0
C
C Define number of elements in subset of Markov chain
	L=20
C
C Define number of Markov chain subsets per temperature
	N=10
C
C Define maximum number of temperature reductions
	K0=500
C
C Define termination threshold
	EPS0=1.D-5
C
C Define ratio of temperature reduction
	RHOT=0.80D0

C
C Set KFLAG to TRUE if output should be writen to a file, otherwise set to FALSE
C	KFLAG=.TRUE.
	KFLAG=.FALSE.
C
C Set MFLAG to TRUE if intermediate Markov-chain candidates should be written to a file,
C otherwise set it to FALSE
C	MFLAG=.TRUE.
	MFLAG=.FALSE.
C
C Read observations
	OPEN(UNIT=10,FILE='ProfileAAB_EGM08_1080.dat',STATUS='OLD')
	I=0
	DO WHILE(.NOT.EOF(10))
	  I=I+1
	  READ(10,*) X(I),OBS(I)
	END DO
	IF (I .NE. NOBS) THEN
	  PRINT*,'Parameter NOBS is incorrect'
	  STOP
	END IF
C
C
C
C Compute gradients for initial parameters
	CALL FAULTGRADB(PARM,NOBS,X,GC)
C

C
C Compute cost function for initial parameters and set to current optimum
	PHI=COST(NOBS,GC,OBS)
	PHIOPT=PHI
C
C Initialize cost function trend
	DO I=1,3
	  PHI1(I)=PHI
	END DO
C
C Initialize optimum parameter values
	DO I=1,M
	  PARMOPT(I)=PARM(I)
	END DO
C
C Initialize loop counter and termination parameter
	K=0
	EPS=1.D0
C
C Open file for output
C	IF (KFLAG .OR. MFLAG)
C	&    OPEN(20,FILE=DIRINB//DIRIN1//FILOUT,STATUS='UNKNOWN')
C
C Write initial values
C	WRITE(6,10) K,T,PARM(1),PARM(2),PARM(3),PHI,EPS,
C	*            DEL(1),DEL(2),DEL(3)
C	IF (KFLAG) WRITE(20,20) K,T,PARM(1),PARM(2),PARM(3),PHI,EPS,
C	*             DEL(1),DEL(2),DEL(3)
C	IF (MFLAG .AND. .NOT.KFLAG) WRITE(20,30) PARM(1),PARM(2),PARM(3)
C
C Master loop for simulated annealing	estimation
	DO WHILE(EPS .GE. EPS0)
C
C Increment loop counter
	  K=K+1
C
C Obtain estimates based on Metropolis algorithm at given temperature
	  CALL METROPOLIS(T,M,NOBS,X,OBS,GC,
	&                  DEL,DEL0,PARM,PARM1,PARMOPT,MJ,XL,XU,PHI)
C
C Compute change in cost function trend
	  EPS1(1)=DABS(PHI-PHIOPT)
	  EPS1(2)=DABS(PHI-PHI1(3))
	  EPS1(3)=DABS(PHI-PHI1(2))
	  EPS1(4)=DABS(PHI-PHI1(1))
C
	  EPS=MAX(EPS1(1),EPS1(2),EPS1(3),EPS1(4))
C
C Update cost function trend
	  PHI1(1)=PHI1(2)
	  PHI1(2)=PHI1(3)
	  PHI1(3)=PHI
	  PHI=PHIOPT	! PHI now refers to optimum parameters out of Metro. alg.
C
C Initialize parameter values for next iteration of Metropolis algorithm
	  DO KK=1,M
	    PARM(KK)=PARMOPT(KK)
	  END DO
C
	  WRITE(6,10) K,T,PHI,EPS,
	&              PARM(1),PARM(2)/DTR,PARM(3),PARM(5),
	&              DUMMY,DUMMY,PARM(4),PARM(6)
   10	  FORMAT(1X,I3,1P3E12.4,/,(1P4E11.3))
C	  IF (KFLAG) WRITE(20,20) K,T,PARM(1),PARM(2),PARM(3),PHI,EPS,
C	*                DEL(1),DEL(2),DEL(3)
C   20	  FORMAT(I5,E14.5,3F12.5,2E13.6,3E13.5)
C	  IF (MFLAG .AND. .NOT.KFLAG) WRITE(20,30) PARM(1),PARM(2),PARM(3)
C   30	  FORMAT(19X,3F12.5)
C
C Reduce temperature
	  T=RHOT*T
C
	  IF (K .GT. K0) STOP
	END DO
C
	STOP
	END


	SUBROUTINE FAULTGRADB(PARM,N,X,GC)
C
C Subprogram to calculate the gradients due to a modeled fault in the
C Wichita Uplift region.  The fault is assumed to be dip-slip fault plane orthogonal
C to the observation profile, which is assumed to be in the X-direction.
C The parameters may include any or all of the following: location of fault [m],
C dip angle of fault [deg->rad], 2 upper boundary depths of density contrasts [m], and 
C 2 density contrasts [kg/m^3]. Depths are negative elevations.
C The gradient is the field curvature magnitude, given in this case by |Gzz| in units of [E].
C
C Input:
C   PARM = parameters (as described above)
C   N = number of computation points
C   X = vector of x-coordinates of points of computation [m]
C
C Output:
C   GC = gradient |Gzz| at points of computation [E]
C
	IMPLICIT REAL*8 (A-H,O-Z)
	DIMENSION X(1),GC(1),PARM(1)
	COMMON /MODEL/Z11,Z12,G2E
C
C Compute sine and cosine of dip angle
	SC=DSIN(PARM(2))
	CC=DCOS(PARM(2))
C
C Main computations for each observation point
	DO J=1,N
	  XJ=X(J)
C
C Compute gradient contributions
C GZZ1 (left of fault)
	  GZZ1=PARM(5)*GZZ(XJ,PARM(1),Z11,PARM(3),SC,CC)
C
C GZZ2 (right of fault)
	  GZZ2=-PARM(6)*GZZ(XJ,PARM(1),Z12,PARM(4),SC,CC)
C
C Combine and include units
	  GC(J)=G2E*DABS(GZZ1+GZZ2)
	END DO
C
	RETURN
	END

	REAL*8 FUNCTION THETA(X,X0,ZP,SC,CC)
	IMPLICIT REAL*8 (A-H,O-Z)
C
	THETA=DATAN2(X-X0-ZP*CC/SC,-ZP)
C
	RETURN
	END

	REAL*8 FUNCTION EE(X,X0,ZP,SC,CC)
	IMPLICIT REAL*8 (A-H,O-Z)
C
	EE=(X-X0)*SC-ZP*CC
C
	RETURN
	END

	REAL*8 FUNCTION FF(X,X0,ZP,SC,CC)
	IMPLICIT REAL*8 (A-H,O-Z)
C
	FF=(X-X0)*SC*CC-ZP
C
	RETURN
	END

	REAL*8 FUNCTION GG(X,X0,ZP,SC,CC)
	IMPLICIT REAL*8 (A-H,O-Z)
C
	GG=(X-X0-ZP*CC/SC)**2+ZP*ZP
C
	RETURN
	END

	REAL*8 FUNCTION AA(X,X0,ZP1,ZP2,SC,CC)
	IMPLICIT REAL*8 (A-H,O-Z)
C
	T1=THETA(X,X0,ZP1,SC,CC)*FF(X,X0,ZP1,SC,CC)
	T2=THETA(X,X0,ZP2,SC,CC)*FF(X,X0,ZP2,SC,CC)
	T3=DLOG(GG(X,X0,ZP1,SC,CC)/GG(X,X0,ZP2,SC,CC))
C
	AA=(T2-T1)/((X-X0)*SC*SC)-0.5D0*T3
C
	RETURN
	END

	REAL*8 FUNCTION DD(X,X0,ZP1,ZP2,SC,CC)
	IMPLICIT REAL*8 (A-H,O-Z)
C
	T1=THETA(X,X0,ZP1,SC,CC)*EE(X,X0,ZP1,SC,CC)
	T2=THETA(X,X0,ZP2,SC,CC)*EE(X,X0,ZP2,SC,CC)
C
	DD=(T2-T1)/((X-X0)**2*SC**3)
C
	RETURN
	END

	REAL*8 FUNCTION GZZ(X,X0,ZP1,ZP2,SC,CC)
	IMPLICIT REAL*8 (A-H,O-Z)
C
	T1=AA(X,X0,ZP1,ZP2,SC,CC)*SC*CC
	T2=(X-X0)*SC*SC*DD(X,X0,ZP1,ZP2,SC,CC)
C
	GZZ=T1-T2
C
	RETURN
	END
	
	SUBROUTINE METROPOLIS(T,M,NOBS,X,OBS,GC,
	&                      DEL,DEL0,PARM,PARM1,PARMOPT,MJ,XL,XU,PHI)
C
C Compute a Markov chain of random variables based on the Metropolis algorithm.
C
C COMMON area MET:
C   PHIOPT = current optimum cost function value
C   L = number of steps in subset of Markov chain with constant step size
C   N = number of subsets of Markov chain wtih constant step size
C   MFLAG = flag to indicate if intermediate Markov-chain candidates should be output
C
C Input:
C   T = temperature of system
C   M = number of parameters
C   NOBS = number of observation points
C   X = x-coordinates of observation points [m]
C   OBS = observations (field curvature magnitude)
C   DEL = incremental steps for parameters (will change on return; i.e. it is also output)
C   PARM = initial parameter estimates (will change on return; i.e. it is also output)
C   XL = vector of lower bounds allowed for PARM
C   XU = vector of upper bounds allowed for PARM
C   PHI = cost function for initial parameters (will change on return; i.e. it is also output)
C
C Work arrays
C   GC = dimension NOBS, used to compute the model values
C   DEL0 = dimension M, used for step sizes
C   PARM1 = dimension M, used as candidate for new parameter values
C   MJ = dimension M, used to record no. of acceptances
C
C Output:
C   PARM = final parameters estimates for this Markov chain (not necessarily optimum)
C   PARMOPT = optimum parameter values for Markov chain
C   PHI = final cost function value (NOT optimum for this Markov chain, based on LAST parameter values)
C   PHIOPT(in common /MET/) = cost function based on optimum parameters for this Markov chain
C   DEL = final step length vector
C   
	IMPLICIT REAL*8 (A-H,O-Z)
	DIMENSION PARM1(M),PARM(M),PARMOPT(M),DEL(M),DEL0(M),XL(M),XU(M)
	DIMENSION X(NOBS),OBS(NOBS),GC(NOBS)
	DIMENSION MJ(M)
	LOGICAL MFLAG
	COMMON /MET/PHIOPT,L,N,MFLAG
	COMMON /ISA/ISEED
C
C Define constant for step change algorithm
	DATA CC/2.D0/
C
C Save initial parameters and step sizes
	DO I=1,M
	  PARM1(I)=PARM(I)
	  DEL0(I)=DEL(I)
	END DO
C
C Loop to generate a Markov chain
	DO J=1,N
C
C Initialize acceptance counters
	  DO I=1,M
	    MJ(I)=0
	  END DO
C
	  DO IL=1,L
C
C Select a new candidate for each parameter using random number (2U-1 is between -1 and 1)
	    DO K=1,M
	      PARM1(K)=PARM(K)+DEL(K)*(2.D0*DRNUNF()-1.D0)
C
C Check if new candidate is within allowed bounds; if not assign random value inside bounds
	      IF (PARM1(K) .LT. XL(K) .OR. PARM1(K) .GT. XU(K))
	&          PARM1(K)=XL(K)+DRNUNF()*(XU(K)-XL(K))
C
C Compute model with these parameter values
	      CALL FAULTGRADB(PARM1,NOBS,X,GC)
C
C Compute cost function with this candidate parameter vector
	      PHI1=COST(NOBS,GC,OBS)
C
C Decide if new cost is consistent with temperature. Select with probability
C if the cost is greater than previous one; otherwise, accept with certainty
C (note CHECK is always <0)
	      DIF=-(PHI1-PHI)/T
		  CHECK=DLOG(DRNUNF())
	      IF (CHECK .LE. DIF) THEN
C
C Accept new parameter and reset cost function to corresponding new value
		    PARM(K)=PARM1(K)
	        PHI=PHI1
C
C Increment acceptance counter for this parameter
	        MJ(K)=MJ(K)+1
C
C Check if new cost is less than current optimum
	        IF (PHI .LT. PHIOPT) THEN
C
C If so, reset current optimum cost and all parameters
	          PHIOPT=PHI
	          PARMOPT(1:M)=PARM(1:M)
	        END IF
C
C Reject new candidate with probability - re-set to old value
	      ELSE
	        PARM1(K)=PARM(K)
	      END IF
C
C Otherwise, do not change the parameter value and do not reset cost function.
	    END DO
	  END DO
	  IF (MFLAG) WRITE(20,20) PARM(1),PARM(2),PARM(3)
   20	  FORMAT(19X,3F12.5)
C
C Change step length for next subset of Markov chain
	  DO I=1,M
	    XMJ=MJ(I)
	    IF (XMJ .LT. 0.4D0*L) THEN		! too few acceptances => decrease del
	      SCALE=1.D0/(1.D0+CC*(0.4D0-XMJ/L)/0.4D0)
	    ELSE IF (XMJ .GT. 0.6D0*L) THEN	! too many acceptances => increase del
	      SCALE=1.D0+CC*(XMJ/L-0.6D0)/0.4D0
	    ELSE
	      SCALE=1.
	    END IF
C
	    DEL(I)=DEL(I)*SCALE
	  END DO
C
C Check if step sizes are too large
	  DO I=1,M
	    IF (DEL(I) .GT. XU(I)-XL(I)) DEL(I)=XU(I)-XL(I)
	  END DO
C
	END DO
C
      RETURN
	END	
	  


	REAL*8 FUNCTION COST(N,GC,OBS)
C
C Define the cost function
C
C Input:
C   N = number of observation points
C   GC = modeled gradients [E]
C   OBS = observed gradients [E]
C
	IMPLICIT REAL*8 (A-H,O-Z)
	DIMENSION GC(1),OBS(1)
C									
C Compute cost function value
	SUM=0.D0
	DO J=1,N
	  SUM=SUM+(GC(J)-OBS(J))**2
	END DO
C
	COST=SUM
C
	RETURN
	END
