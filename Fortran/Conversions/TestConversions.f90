PROGRAM TestConversions

	USE TemperatureConversions

	REAL :: TempF1, TempF2, TempF3, TempF4
	REAL :: TempK1, TempK2

	TempF1 = 100.0     ! Fahrenheit
	TempF2 = -40.0	   ! Fahrenheit

	CALL FahrenheitToKelvin(TempF1, TempK1)

	CALL FahrenheitToKelvin(TempF2, TempK2)

	CALL KelvinToFahrenheit(TempK1, TempF3)
	
	CALL KelvinToFahrenheit(TempK2, TempF4)

	PRINT*, TempF1, TempF3
	PRINT*, TempF2, TempF4

END PROGRAM TestConversions
