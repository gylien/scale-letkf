        !COMPILER-GENERATED INTERFACE MODULE: Thu Nov 21 01:39:18 2019
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE RADAR_ECHO__genmod
          INTERFACE 
            SUBROUTINE RADAR_ECHO(NLN,DHGT,KEXTTOT,SBACKTOT,LAMBDA,     &
     &KRADAR,Z)
              INTEGER(KIND=4), INTENT(IN) :: NLN
              REAL(KIND=8), INTENT(IN) :: DHGT(1:NLN)
              REAL(KIND=8), INTENT(IN) :: KEXTTOT(NLN)
              REAL(KIND=8), INTENT(IN) :: SBACKTOT(NLN)
              REAL(KIND=8), INTENT(IN) :: LAMBDA
              REAL(KIND=8), INTENT(IN) :: KRADAR
              REAL(KIND=8), INTENT(OUT) :: Z(NLN)
            END SUBROUTINE RADAR_ECHO
          END INTERFACE 
        END MODULE RADAR_ECHO__genmod
