MODULE TemperatureConversions

	IMPLICIT NONE

	! Specification Part

	REAL, PARAMETER :: CtoK = 273.15
	REAL, PARAMETER :: KtoC = -273.15
	REAL, PARAMETER :: FtoC = 0.5555
	REAL, PARAMETER :: CtoF = 1.8

CONTAINS

	! Internal Functions

	SUBROUTINE FahrenheitToCelsius(F,C)

		IMPLICIT NONE
		
		REAL, INTENT(IN)  :: F
		REAL, INTENT(OUT) :: C

		C = FtoC * (F - 32.0)

	END SUBROUTINE FahrenheitToCelsius

	SUBROUTINE CelsiusToFahrenheit(C,F)

		IMPLICIT NONE

		REAL, INTENT(IN)  :: C
		REAL, INTENT(OUT) :: F

		F = C * CtoF + 32.0

	END SUBROUTINE CelsiusToFahrenheit

	SUBROUTINE CelsiusToKelvin(C,K)

		IMPLICIT NONE

		REAL, INTENT(IN)  :: C
		REAL, INTENT(OUT) :: K

		K = C + CtoK

	END SUBROUTINE CelsiusToKelvin

	SUBROUTINE KelvinToCelsius(K,C)

		IMPLICIT NONE

		REAL, INTENT(IN)  :: K
		REAL, INTENT(OUT) :: C

		C = K + KtoC

	END SUBROUTINE KelvinToCelsius

	SUBROUTINE FahrenheitToKelvin(F,K)

		IMPLICIT NONE

		REAL, INTENT(IN)  :: F
		REAL, INTENT(OUT) :: K
		REAL              :: C

		CALL FahrenheitToCelsius(F,C)
		CALL CelsiusToKelvin(C,K)

	END SUBROUTINE FahrenheitToKelvin

	SUBROUTINE KelvinToFahrenheit(K,F)

		IMPLICIT NONE

		REAL, INTENT(IN)  :: K
		REAL, INTENT(OUT) :: F
		REAL              :: C

		CALL KelvinToCelsius(K,C)
		CALL CelsiusToFahrenheit(C,F)

	END SUBROUTINE KelvinToFahrenheit

END MODULE TemperatureConversions
