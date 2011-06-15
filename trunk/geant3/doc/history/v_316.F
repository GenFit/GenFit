*
* $Id: v_316.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: v_316.F,v $
* Revision 1.1.1.1  2002/06/16 15:18:35  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:15  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:19:48  cernlib
* Geant
*
*
#include "sys/CERNLIB_machine.h"
#include "pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.18  by  S.Giani
*-- Author :
*
*     Correction cradle for version 3.15 of the GEANT Programme
*     and history of the evolution of 3.16.
*
*               Important !!!
*               =============
*
*     Only corrections flagged by "Introduced in the 3.15 correction
*     cradle" are applied to 3.15. All the other comments refer
*     to modifications which will appear in the 3.16 version of
*     the program.
*
********************************************************************
*
*===> 06/12/93
*
*           Release of the 3.11 correction cradle.
*
*===> 02/12/93
*
*           Correction in GBRELE and GRANGI. Thanks to F.Carminati.
*           Introduced in the 3.15 correction cradle.
*
*===> 29/11/93
*
*           Release of the 3.10/00 correction cradle.
*
*===> 21/10/93
*
*          Correction in FKDECA and HADEVV to permit the correct
*          treatment of charmed particles decay. Thanks to
*          A.Ferrari. Introduced in the 3.15 correction cradle.
*
*          Correction in FKDRES to avoid a loop. Thanks to
*          B.Slater and A.Ferrari. Introduced in the 3.15 correction
*          cradle.
*
*===> 20/10/93
*
*          Correction in GDRAW to update current pointer to IN
*          bank. Correction in GGPPAR to compute the internal
*          parameters of a TRAP with negative parameters. Thanks
*          to S.Banerjee. Introduced in the 3.15 correction cradle.
*
*===> 19/10/93
*
*           Correction in GTNEUT, GTHADR, to allow hadronic decay
*           in the vacuum. Thanks to R.Jones. Introduced in the
*           3.15 correction cradle.
*
*===> 01/10/93
*
*           Correction in GNSLWD to avoid numerical problem. Thanks
*           to A.Parri. Introduced in the 3.15 correction cradle.
*
*===> 15/09/93
*
*           Release of the 3.09/00 correction cradle.
*
*===> 13/09/93
*
*           Correction in GNOCTU, thanks to M.Nessi. Introduced in
*           the 3.15 correction cradle.
*
*===> 03/09/93
*
*           Correction in EVENTV to protect square root of a negative
*           number. Thanks to B.Slater and A.Ferrari. Introduced in
*           the 3.15 correction cradle.
*
*===> 11/08/93
*
*           New tratment for thin layer, mainly gas, of dE/dx
*           introduced. Thanks to V.Grishin, S.Kotelnikov,
*           K.Lassila-Perini and P.Nevski.
*
*           All bit manipulation functions replaced by Military
*           Standard functions. Thanks to F.Carminati.
*
*===> 19/07/93
*
*           Correction in FLUFIN in order to avoid 0 mechanism
*           number when IHADR=2. Thanks to P.Bloch. Introduced
*           in the 3.15 correction cradle.
*
*           Correction in GSROTM, NROTM updated when the bank is
*           pushed. Thanks to M.Maire. Introduced in the 3.15
*           correction cradle.
*
*           Correction in ERTRGO. Thanks to E.Nagy. Introduced
*           in the 3.15 correction cradle.
*
*===> 07/06/93
*
*           Release of version 3.08 of the correction cradle.
*
*           New version of the FLUKA interface introduced. New
*           model for energy straggling introduced (Ermilova&
*           Checkin) and synchrotron radiation. Thanks to K.
*           Lassila-Perini, A.Ferrari&P.Sala et al. and I.
*           Gavrilenko. Please see the documentation for more
*           information.
*
*           New routine GFANG introduced. It finds the angles
*           Theta and Phi which are used to rotate a vector in
*           the master system of the other. Thanks to F.Carminati
*
*           New rotine GVROT introduced. It will rotate in full
*           double precision a vector in the reference system of
*           another one. Thanks to F.Carminati
*
*           Corrections in GRIN/GFIN to calculate correctly the
*           number of rotation matrixes. Thanks to M.Maire.
*           Introduced in the 3.15 correction cradle.
*
*           Correction in FLUFIN to pass antiprotons to FLUKA when
*           they are 'stopping'. Thanks to A.Ferrari. Introduced in
*           the 3.15 correction cradle.
*
*
*===> 27/05/93
*
*           Corrections in FLUFIN to pass stopping antiprotons,
*           antineutrons and pions to FLUKA for annihilation. The
*           other hadrons are either decayed or (neutral particles)
*           passed to GHEISHA/GHSTOP.
*
*===> 04/05/93
*
*           Correction in FDNOPT to properly initialise an
*           array. Thanks to K.Lassila-Perini. Introduced in
*           the 3.15 correction cradle.
*
*===> 21/04/93
*
*           New routine from Laszlo Urban for the tracking of
*           ions, GTHION. Ions now have to be given particle
*           type 8. Thanks to L.Urban.
*
*===> 18/04/93
*
*           Addition in GPPERP, the sphere has been added. Thanks
*           to R.Jones.
*
*           Old Hadrin/Nucrin interface deleted. HADR 3 flag has
*           no meaning any more. Warning message introduced in
*           GPHYSI, when IHADR=3 then reset it to 1 and use GHEISHA.
*           Thanks to F.Carminati.
*
*           Tatina code eliminated from GEANT. Thanks to F.Carminati.
*
*           Small changes in GSCHIT to avoid integer/floating
*           conversions and to speed-up the code. Thanks to
*           F.Carminati.
*
*===> 07/04/93
*
*           Correction in the GPOISS routine to avoid overflows.
*           Thanks to S.Qian. Introduced in the 3.15 correction
*           cradle.
*
*===> 06/04/93
*
*           Correction in GNSPHR to avoid floating point overflow.
*           Thanks to S.Ljungfelt. Introduced in the 3.15 correction
*           cradle.
*
*===> 22/03/93
*
*           Correction in GEAMON to extend the field for CPU
*           time. Thanks to R.Brun. Introduced in the 3.15
*           correction cradle.
*
*===> 14/03/93
*
*           Release of version 3.07 of the correction cradle.
*
*           Corrections in GZINIT, GPHYSI and GDSHOW. The
*           working space is now allocated via calls to
*           GWORK. GWORK has been modified and it only
*           increases the working space, it never reduces it.
*           Thanks to M.Corden. Introduced in the 3.15
*           correction cradle.
*
*===> 10/03/93
*
*           Correction in GFIN/GFOUT to avoid spurious error
*           message. Thanks to M.Maire.
*
*===> 03/03/93
*
*           Correction in GLTRAC to improve the precision of
*           the calculation of the kinetic energy. Thanks to
*           C.Zeitnitz. Introduced in the 3.15 correction
*           cradle.
*
*           Correction in GTRACK to avoid output overflow in
*           format statement. Thanks to M.Corden. Introduced
*           in the 3.15 correction cradle.
*
*           Correction in GDXYZ to obtain the correct line color
*           for a track. Thanks to M.Corden. Introduced in the
*           3.15 correction cradle.
*
*===> 02/03/93
*
*           New treatment of the photoelectric effect introduced.
*           Thanks to J.Chwastowski. Please see documentation for
*           more information.
*
*===> 01/03/93
*
*           Changes in GKGEOM to avoid warnings with the new
*           version of KUIP. Thanks to A.Nathaniel. Introduced
*           in the 3.15 correction cradle.
*
*           Changes in GSDVN2 to support division of spheres
*           in PHI with initial offset. Thanks to F.Carminati.
*
*           Changes in GMEDIV and GNSPHR to correct precision
*           problems. Thanks to F.Carminati.
*
*           Introduction of the Cerenkov effect. Thanks to
*           R.Jones. Please see documentation for more
*           details.
*
*===> 12/02/93
*
*           Changes in GSDVN to support division in PHI (axis
*           N. 3) of the sphere. Thanks to F.Carminati.
*
*           Changes in GDRAWS to partially support the drawing
*           of parts of sphere. Thanks to F.Carminati.
*
*===> 09/02/93
*
*           Release of version 3.06 of the correction cradle.
*
*===> 28/01/93
*
*           Correction in GXINT for compatibility with the KUIP
*           installed with version 93a of the CERN Program Library.
*           Thanks to R.Jones. Introduced in the 3.15 correction
*           cradle.
*
*===> 27/01/93
*
*           Correction in GMUNU. Only GHEISHA used to handle the
*           interaction irrespective of the value of IHADR. Thanks
*           to K.Lassila-Perini and G.Battistoni. Introduced in the
*           3.15 correction cradle.
*
*           Corrections in GGPPAR to handle a tube with negative
*           parameters positioned into a tubes in the general case.
*           Thanks to F.Carminati.
*
*===> 19/01/93
*
*           Correction in FLUFIN to make sure that just neutrons
*           are passed to GHEISHA. Thanks to K.Lassila. Introduced
*           in the 3.15 correction cradle.
*
*===> 21/12/92
*
*           Option Q(uiet) introduced in GRIN/GROUT. Thanks to
*           M.Maire.
*
*           New version of GRKUTA protected against 0 fields and
*           improved for precision. Introduced in the 3.15 correction
*           cradle. Thanks to W.Manner.
*
*===> 15/12/92
*
*           IQ(JRUNG+29) set to the event sequence number of the
*           last generated event in GTRIGC. Thanks to M.Maire.
*
*===> 11/12/92
*
*           Call to LUPART replaced with LU1ENT in GLUDKY. LUPART
*           was an obsolete routine deleted from JETSET 7.3
*           library. Thanks to R.Gomes.
*
*===> 23/11/92
*
*           Release of the version 3.05/00 of the correction
*           cradle.
*
*           Corrections in GRANGI to account for the use of the
*           ERAN card. Thanks to P.Gumplinger. Introduced in the
*           3.15 correction cradle.
*
*===> 18/11/92
*
*           Correction in GDRAWS for very thin CTUBes. Thanks to
*           M.Guidal. Introduced in the 3.15 correction cradle.
*
*           New version of GRKUTA introduced. It will limit the
*           number of iterations and cut steps to prevent looping.
*           Thanks to V.Perevotchikov. Introduced in the 3.15
*           correction cradle.
*
*===> 11/11/92
*
*           Correction in GGDSPE to allow repeated division of
*           tubes along R. Thanks to M.Maire. Introduced in the
*           3.15 correction cradle.
*
*===> 11/11/92
*
*           Correction in GNPCON when the particle is slightly
*           outside of the last plane in Z. It could happen for
*           precision problems. Thank to C.Georgiopoulos.
*           Introduced in the 3.15 correction cradle.
*
*           Quiet option introduced in GFIN and GFOUT to suppress
*           printing. Thanks to A.Eisner.
*
*           Correction in the header of GFIN, GFOUT. Thanks to
*           A.Eisner.
*
*           Correction to GSNEXT. It now checks the values of the
*           input list. Repeated values or values which are not
*           within the allowed limits are ignored and a warning
*           message is printed. Thanks to F.Luehring.
*
*===> 04/11/92
*
*           Modification in GSVOLU. A warning is printed if the same
*           volume name is repeated twice. Thanks to V.Perevoztchikov.
*
*===> 23/10/92
*
*           Print statement in GSKING replaced by a WRITE(CHMAIL
*           Thanks to F.Carminati.
*
*===> 19/10/92
*
*           Correction in GDRAY to increase the precision of
*           the angle calculation. Thanks to F.Malek.
*           Introduced in the 3.15 correction cradle.
*
*===> 19/10/92
*
*           Release of version 3.04 of the correction cradle.
*
*           Corrections in FEREVV and HADEVV. Thanks to A.Ferrari.
*           Introduced in the 3.15 correction cradle.
*
*===> 16/10/92
*
*           Release of version 3.03 of the correction cradle.
*
*===> 13/10/92
*
*           Correction in GFIN and GFOUT to handle correctly
*           ZEBRA error conditions. Thanks to H.Uijterwaal.
*           Introduced in the 3.15 correction cradle.
*
*           Correction in GDECAY to improve precision of rotation
*           of secondaries. Thanks to F.Carminati. Introduced in
*           the 3.15 correction cradle.
*
*           Correction in FEREVV to avoid square root of negative
*           number. Thanks to A.Ferrari. Introduced in the 3.15
*           correction cradle.
*
*           Correction in GDRAW to avoid error when traversing the
*           tree. Thanks to F.Carminati. Included in the 3.15
*           correction cradle.
*
*===> 09/10/92
*
*           Correction in GXPICT, added definition of a character
*           variable. Thanks to F.Carminati. Introduced in the 3.15
*           correction cradle.
*
*===> 08/10/92
*
*           Corrections in ERPRDOC, EUFILV and ERTEST to change
*           Hollerith in characters for GEANE. Thanks to V.Innocente.
*           Introduced in the 3.15 correction cradle.
*
*===> 24/09/92
*
*           Correction in NUCREL to avoid square root of negative
*           numbers. Thanks to A.Ferrari. Introduced in the 3.15
*           correction cradle.
*
*===> 23/09/92
*
*           Correction in GMEDIA. INFR not set going up the tree to
*           avoid problem with user calls of GMEDIA during tracking.
*           Thanks to S.Banerjee. Introduced in the 3.15 correction
*           cradle.
*
*===> 21/09/92
*
*           Correction in GRANGI to avoid interpolation of DE/DX
*           beyond 10GeV. Thanks to V.Innocente. Introduced in the
*           3.15 correction cradle.
*
*===> 16/09/92
*
*           Correction in GBRELA. Auxiliary integration point for
*           muon dedx table was wrong. Thanks to V.Innocente.
*           Introduced in the 3.15 correction cradle.
*
*===> 15/09/92
*
*           Release of version 3.02 of the correction cradle.
*
*           Correction in GMEDIA, GTMEDI: JVIN not properly
*           initialised. Thanks to F.Carminati. Introduced in
*           the 3.15 correction cradle.
*
*           Release of version 3.01 of the correction cradle.
*
*           Correction to GSDVN for a mistype. Thanks to H.Plothow.
*           Introduced in the 3.15 correction cradle.
*
*           Small mod to GSNGTR for the CRAY. Function DSIGN changed
*           to SIGN. Thanks to F.Carminati. Introduced in the 3.15
*           correction cradle.
*
*           Release of version 3.00 of the correction cradle.
*
*           Correction in GDECAY. When IDCAY=2 the routine does not
*           produce secondaries. Thanks to R.Jones. Introduced in
*           the 3.15 correction cradle.
*
*           Corrections in GFIN, GFOUT. Thanks to B.Cole. Introduced
*           in the 3.15 correction cradle.
*
*===> 14/09/92
*
*           Correction in GMEDIV, GMEDIA, GINVOL, GTMEDI to
*           improve precision and prevent blocking of volumes.
*           Thanks to F.Carminati, R.Jones. Introduced in the
*           3.15 correction cradle.
*
*===> 03/09/92
*
*           GDLADY and GDMISS routines introduced. Thanks to
*           H.Plothow.
*
*           Correction in GXDRAW for the coordinates of the
*           drawing of the GIRL. Thanks to F.Carminati.
*           Introduced in the 3.15 correction cradle.
*
*           Correction in GPHSIG, typing mistake in a formula.
*           Thanks to M.Fanti. Introduced in the 3.15 correction
*           cradle.
*
*===> 02/09/92
*
*           Corrections in GSDVN, GGPPAR, GGDPAR, GMEDIV, GCENT,
*           GSDVN2, GSDVT, GSDVT2 to allow the division of a TRAP
*           along y. Thanks to F.Carminati. Introduced in the 3.15
*           correction cradle.
*
*           Corrections in GNPGO1 for precision problems. Thanks
*           to K.Lassila. Introduced in the 3.15 correction
*           cradle.
*
*           Corrections in GGDPAR, GSDVN for the z division of a
*           Polygone. Thanks to F.Carminati. Introduced in the
*           3.15 correction cradle.
*
*===> 01/09/92
*
*           Correction in GSDVT for the division of a parallepiped.
*           Thanks to Y.Fisiak. Introduced in the 3.15 correction
*           cradle.
*
*           Correction in GRKUTA. Helix is used if angle of
*           curvature is too large. In this case step is anyway
*           smaller than STMIN which is under user control.
*           Thanks to H.Plothow and V.Perevotchikov. Introduced
*           in the 3.15 correction cradle.
*
*===> 26/08/92
*
*           Correction in GKRZ and CDF for the command RZ/FILE.
*           Record lenght fixed at 1024 in GRFILE cannot be
*           Changed by the user. Thanks to H.Plothow. Introduced
*           in the 3.15 correction cradle.
*
*           Corrections in GTMEDI, GMEDIA and GINVOL for the case
*           of 'MANY' volumes. Uncleared blocking words found in
*           the tree. Thanks to H.Plothow. Introduced in the 3.15
*           correction cradle.
*
*===> 25/08/92
*
*           Correction in GNSPHR. The PHI limits were used
*           incorrectly in calculating the distance to the
*           sphere. This was slowing down calculation but still
*           giving correct results. Thanks to J.Stepanek.
*           Introduced in the 3.15 correction cradle.
*
*===> 21/08/92
*
*           Correction in GGORDQ. Incorrect counting of volumes
*           in case of phi ordering and phi division spanning
*           across 0 degrees. Thanks to L.Roberts and F.Carminati.
*           Introduced in the 3.15 correction cradle.
*
*===> 18/08/92
*
*           Version 2.00 of the correction cradle released
*
*
*           Correction in GDRELP and GDRELM. The shell corrections
*           where underestimated due to a wrong factor introduced
*           translating the original formula from ev to Gev. Thanks
*           to A.Ferrari and K.Lassila-Perini. Introduced in the
*           3.15 correction cradle.
*
*===> 17/08/92
*
*           Correction in GCENT in case of general trapezoid sliced
*           in z in a local development structure.
*           Thanks to S.Chekhtman. Introduced in the 3.15 correction
*           cradle.
*
*           Correction in GHEISH. Time of flight correctly initialised
*           for particles with number > 48. Thanks to A.Eisner.
*           Introduced in the 3.15 correction cradle.
*
*===> 14/08/92
*
*           Correction in GGPPAR for boxes with negative parameters
*           in traps. The current code supports correctly boxes with
*           one, two or three negative parameters positioned in a trap
*           not necessarily in the centre. The trap must not be twisted
*           and the ALP1 and ALP2 angles (PAR(7) and PAR(11)) must be
*           equal. This facility must be used with care because the
*           algorithm first maximises Z, then Y and then X. The result
*           may be not quite what is expected in all cases, and the
*           user is invited to check carefully the result. Thanks to
*           H.Plothow. Introduced in the 3.15 correction cradle.
*
*===> 12/08/92
*
*           Correction in GINVOL, GTMEDI and GMEDIA to handle correctly
*           the GSNEAR case when LIST(1)=1 and LIST(2)=0, that is when
*           all the volume are separated inside a mother. Optimisation
*           of the routine in case a large number of objects are
*           positioned at the same level and GUNEAR is used. Thanks to
*           R.Jones, M.Lefebvre and M.Pieri. Introduced in the 3.15
*           correction cradle.
*
*===> 07/08/92
*
*           Correction in GTNEXT. When the currents step brings
*           the particle on the surface of the next volume, signal
*           IGNEXT=1 to GTxxxx to 'push' into the next volume.
*           This is a genuine bug and some volumes may be lost
*           due to this. Thanks to M.Maire. Introduced in the 3.15
*           correction cradle.
*
*===> 02/08/92
*
*           Corrections to the precision in GNOPG1, GNCONE, GNOCON.
*           Thanks to R.Jones. Introduced in the 3.15 correction
*           cradle.
*
*===> 31/07/92
*
*           Correction in the printout of GTRACK when max number
*           of steps is reached. Event number and random numbers are
*           printed out as well. Thanks to M.Maire. Introduced in the
*           3.15 correction cradle.
*
*===> 29/07/92
*
*           Correction in PIMABS for HPUX only. If variable NT is
*           saved then the routine compiles correctly, otherwise
*           not. Thanks to F.Carminati. Introduced in the 3.15
*           correction cradle.
*
*           Correction in NIZL. All local variable are saved now.
*           Thanks to F.Carminati. Introduced in the 3.15 correction
*           cradle.
*
*===> 29/07/92
*
*           Correction in CDF and in GKDRAW. Minimum angle for the
*           drawing is 0 and not 1. Thanks to H.Plothow. Introduced
*           in the 3.15 correction cradle.
*
*           Minor update of GXFZ, a character length has been set
*           to 4 to keep KUIP happy. Thanks to F.Carminati. Introduced
*           in the 3.15 correction cradle.
*
*===> 28/07/92
*
*           Corrections in GXFZ. Now the array LUNITS is properly
*           updated when a file is closed.
*
*           Correction in GBRSGE, GBRSGM. Save statement introduced.
*           Thanks to L.Roberts. Introduced in the 3.15 correction
*           cradle.
*
*           Correction in GINVOL for many volumes. Thanks to S.Giani
*           and R.Jones. Introduced in the 3.15 correction cradle.
*
*           Correction in GMCOUL to prevent a division by 0. Thanks
*           to M.Maire. Introduced in the 3.15 correction cradle.
*
*===> 03/06/92
*
*           Correction in GSNGTR. Double precision variables removed
*           in case of 64 bits machines. Thanks to E.Boehm. Introduced
*           in the 3.15 correction cradle.
*
*           Correction in the SROTM command. The lower limit for
*           PHI3 is really 0. and not 1. Thanks to E.Boehm. Introduced
*           in the 3.15 correction cradle.
*
*===> 01/06/92
*
*           New version of GRKUTA. In fact this is the 3.14 version
*           with the multiple steps algorithm of 3.15 introduced.
*           If the bending angle in one step is too big, then the
*           integration is done in several steps rather than in
*           one. The IMPLICIT REAL*8 in the cradle is to avoid
*           problems with the CRAY in the absence of PATCHY delayed
*           control cards. Thanks to V.Perevotchikov. Introduced in
*           the 3.15 correction cradle.
*
*===> 25/05/92
*
*           Sequence GCONSP changed for 32 bits machines. All
*           constants are suffixed with D0 to avoid truncations.
*           Thanks to F.Carminati. Introduced in the 3.15 correction
*           cradle.
*
*===> 25/05/92
*
*           Version 1.00 of the correction cradle released
*
*===> 21/05/92
*
*          Correction in GHFNUC: TNNUC variable was undefined.
*          Thanks to B.Lund-Jensen.
*          Introduced in the 3.15 correction cradle.
*
*===> 08/05/92
*
*          Correction in TWOCLU to avoid infinite looping in some
*          pathological cases. Thanks to L.Roberts.
*          Introduced in the 3.15 correction cradle.
*
*===> 06/05/92
*
*          Correction in GNPGO1. Constants changed to double
*          precision to avoid problems on the IBM. Thanks to
*          D.Ward.
*          Introduced in the 3.15 correction cradle.
*
*===> 28/04/92
*
*          Correction in GXINT to obtain the correct time in
*          the monitoring routine when running interactively
*          on IBM/VM. Thanks to F.Carminati.
*          Introduced in the 3.15 correction cradle.
*
*===> 16/04/92
*
*          Correction in GMULTS. A particle going exactly along
*          the negative direction of the Z axis could be rotated
*          of 180 degrees after multiple scattering. Thanks to
*          P.Gabriel.
*          Introduced in the 3.15 correction cradle.
*
*===> 15/04/92
*
*          Correction in GXDRAW. The arguments of DHEAD were
*          fetched with the wrong order. Thanks to F.Videbaek.
*          Introduced in the 3.15 correction cradle.
*
*          Correction in GPHYSI. Variable NUMOLD is now reset to
*          zero. This may affect people running interactively and
*          using the command PHYSI. Thanks to F.Videbaek.
*          Introduced in the 3.15 correction cradle.
*
