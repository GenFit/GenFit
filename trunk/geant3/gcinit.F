      SUBROUTINE G3CINIT
      PARAMETER (NWGEANT=5000000)
      COMMON / GCBANK / GMEM(NWGEANT)
*
* Initialise geant system with a default dimension of NWGEANT
*
      CALL G3ZEBRA(NWGEANT)
      CALL G3INIT()
      CALL G3ZINIT()
*
      END
