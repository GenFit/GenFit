#include "geant321/pilot.h"
      subroutine g3tmany(level1)
#include "geant321/gcvol1.inc"
#include "geant321/gcvolu.inc"
#include "geant321/gcvdma.inc"
#include "geant321/gctrak.inc"
      integer level1

      if(mycoun.gt.1.and.nfmany.gt.0.and.step.ge.safety)then
        nlevel=manyle(nfmany)
        do 99 i=1,nlevel
          names(i)=manyna(nfmany,i)
          number(i)=manynu(nfmany,i)
 99     continue
        call glvolu(nlevel,names,number,ier)
        if(ier.ne.0)print *,'Fatal error in GLVOLU'
        ingoto=0
      endif
      if (level1.eq.1) then
            NLEVL1 = NLEVEL
            DO 20   I=1,NLEVEL
               NAMES1(I) = NAMES(I)
               NUMBR1(I) = NUMBER(I)
               LVOLU1(I) = LVOLUM(I)
   20       CONTINUE
      endif
      end
