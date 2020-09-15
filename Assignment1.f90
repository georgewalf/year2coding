!----------------------------------------
!PHY2063 Assignment 1 - Monte Carlo Method
!URN 6457454
!November 20th 2017
!----------------------------------------

PROGRAM Assignment1
	IMPLICIT none
	INTEGER :: i, j
	INTEGER, DIMENSION (1:1) :: imax
	INTEGER, PARAMETER :: m=1050, n=10, taumax=100, taumin = 1
	REAL :: prior(0:m), prob(0:m), tau(0:m), decay(1:n), normconst(1:1), dtau=0.1

	WRITE(6,*) "This program will estimate the probability function for the lifetime"
	WRITE(6,*) "of an isotope, after a set of 10 decay times have been measured."
	WRITE(6,*) ""

	OPEN(unit=20, file='bayes_data.txt')
	READ(20,*) decay 																								!reads data for measured decays from bayes_data.txt data file
	
	DO i = 0,m 
		tau(i) = REAL(i) * dtau 																			!time counter goes up by 0.1s increments	
		IF (tau(i) <= taumax .AND. tau(i) >= taumin) THEN
			prior(i) = (1.0/(taumax - taumin))													!sets probability for tau to be between 1s and 100s			
		ELSE 
			prior(i) = 0																					 		  !sets all other probability = 0			
		END IF		
	prob(i) = prior(i) 																							!sets probability to prior	
	END DO

	DO j = 1,10 																										!outer do loop to ensure all values of decay are updated for all values in array
		DO i = 0,m 																										!inner do loop updates prob values in array
			IF (tau(i) == 0) THEN 																			
				prob(i) = 0
			ELSE
				prob(i) = (prob(i) * EXP((-decay(j)/tau(i)))/tau(i)) 			!calculates values of probabilities				
			END IF
		END DO
	END DO

	imax = MAXLOC(prob)																							
	normconst = (SUM(prob) - (0.5*(prob(0)+prob(imax)))) * dtau		 	!normalises probability density function
	prob = prob/normconst(1)
	
	OPEN(unit=30, file='results.dat')
	DO i = 0,m
	WRITE(30,'(5x,f5.2,5x,f7.5)') tau(i), prob(i)
	END DO
	
	WRITE(6,*) "The results can be plotted with GNUplot with the data that has been output to results.dat"
	WRITE(6,*) " "
	WRITE(6,'(a,f5.2)') " The most probable value of tau is: ", (iMax * DTau)
		
	CLOSE(20)
	CLOSE(30)

END PROGRAM






