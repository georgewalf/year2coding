!------------------------------------------
!Assignment 2 - Numerical Solutions of PDEs
!URN 6457454
!January 7th 2018
!------------------------------------------
PROGRAM hotcopperrod

	IMPLICIT none     
	INTEGER :: x, t
	INTEGER, PARAMETER :: n = 300, TimeMin = 0, TimeMax = 10, DisplacementMin = -30, DisplacementMax = 30
	REAL :: Temperature(-n:n) = 0.0, Displacement(-n:n) = 0.0, TemperatureDerivative(-n:n) = 0.0
	REAL, PARAMETER :: DiffusionConstant = 1.1, DisplacementStep = 0.1, TimeStep = 0.0001		
	
	WRITE(6,*) ""
	WRITE(6,*) "This program will calculate and produce a profile of a diffused temperature distribution along a copper wire."
	WRITE(6,*) ""
	WRITE(6,*) "The boundary conditions are as follows:"
	WRITE(6,*) ""
 	WRITE(6,*) " -From region of",DisplacementMin,"cm to ",DisplacementMax,'cm'
  	WRITE(6,'(a,f4.1,a)') " -With a grid spacing of ",DisplacementStep,"cm"
	WRITE(6,*) " -An evolution of temperature from ", TimeMin, "s to ", TimeMax, "s" 
	WRITE(6,'(a,f8.4,a)') " with a time step of ", TimeStep, "s" 
	WRITE(6,*) ""
	WRITE(6,*) "The initial boundary conditions are as follows:"
	WRITE(6,*) "At 0 degrees Celsius: x < -1cm"
	WRITE(6,*) "At 0 degrees Celsius: x > 1cm"
	WRITE(6,*) "At 20 degrees Celsius: -1cm <= x <= 1cm"
	WRITE(6,*) ""
	WRITE(6,'(a,f3.1,a)') "A Diffusion Constant of ", DiffusionConstant, " was used - defined" 
	WRITE(6,*) "for copper at room temperature."


	DO x= -10, 10
		Temperature(x) = 20.0				
				!sets the region -1cm <= x <= 1cm to T = 20
	END DO	

	
	DO t = TimeMin, INT(TimeMax/TimeStep)-1 		
				!integrating the PDE in small time steps
		DO x = -n+1, n-1 				
				!loop runs through -n+1 to n-1 to prevent error involving array accessing elements that aren't there
		
			TemperatureDerivative(x) = DiffusionConstant *((Temperature(x+1) - 2.0*Temperature(x) + Temperature(x-1))/DisplacementStep**2) 
				!finds the time derivative of temperature, using the formula given.
						
			Temperature(x) = Temperature(x) + (TimeStep*TemperatureDerivative(x)) 	
				!then integrates, calculating temperature values across the wire at next time step.
		END DO
	END DO
	
	OPEN(20,FILE='hotrod.dat')
		
	DO x = -n, n
		Displacement(x) = REAL(x) * DisplacementStep
		WRITE(20,*) Displacement(x), Temperature(x)    
				!writes temperature profile to dat file
	END DO

	CLOSE(20)
	
	WRITE(6,*) ""
	WRITE(6,*) "The solution has been exported to HOTROD.dat"
	
END PROGRAM





















