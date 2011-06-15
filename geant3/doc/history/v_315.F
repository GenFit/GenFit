*
* $Id: v_315.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: v_315.F,v $
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
*
*        Maintenance and support for the GEANT MonteCarlo system
*        are performed by the Simulation Section of the CERN
*        Computing and Networks division. For any problem please
*        contact:
*
*                        Federico Carminati
*                        CERN-CN
*                        1211 Geneva 23
*                        Tel:        +41 22 767.4959
*                        Telefax:    +41 22 767.7155
*                        E-mail:
*                        BITnet/EARN: FCA@CERNVM
*                        DECnet     : VXCERN::FCA (Node 22.190)
*                        Internet   : fca@cernvm.cern.ch
*
*        An electronic GEANT discussion list has been set up on CERNVM.
*        To be included in the mailing list it is enough to send the
*        following message:
*
*       TELL LISTSERV AT CERNVM SUBSCRIBE LGEANT Firstname Lastname Inst
*
*        People who do not dispose of a TELL (or equivalent) bitnet
*        command should send a mail to LISTSERV@CERNVM containing the
*        only line:
*
*                   SUBSCRIBE LGEANT Firstname Lastname Inst
*
*================================= History
*
*===> 20/NOV/1990
*
*        Current plans set the release of GEANT 3.15 around the third
*        quarter of 1991.
*
*===> 27/NOV/1990
*
*        Correction in GTRACK to set correctly the value of PREC
*        to match machine precision. Thanks to Stephan Egli and
*        M.Maire.
*        Introduced in the 314 correction cradle.
*
*===> 28/NOV/1990
*
*        Correction in GNTUBE to avoid square root of a negative
*        number. Thanks to M.Maire
*        Introduced in the 314 correction cradle.
*
*        Eliminate TOFG=0 from GINIT, already done in GTRIGI
*        Thanks to F.Carminati
*        Introduced in the 314 correction cradle.
*
*===> 29/NOV/1990
*
*        UNIX and BSLASH flags correctly set inside *GEANT.
*        Thanks to A.Nathaniel
*        GPMATE now prints also the components in case of
*        mixtures. If NUMB is <0 then materials from 1 to
*        -NUMB are printed. Thanks to Stephan Egli.
*
*        Corrections in GRGET (obsolete), GGET, GRIN to reset
*        correctly the number in GCNUM. Thanks to Stephan Egli
*        Introduced in the 314 correction cradle.
*
*===> 30/NOV/1990
*
*        New version of GNOELT and GNELTU, thanks to Yoshihisa Iga
*        and R.Nierhaus
*        Introduced in the 314 correction cradle.
*
*        Correction in GTRAK, all mechanisms turned off in the
*        vacuum. Thanks to S.Egli
*        Introduced in the 314 correction cradle.
*
*
*===> 03/DEC/1990
*
*        New routine GXPICK for the interactive version (Command
*        DRAWING/PICK) to pick a point in a detector and return
*        volume/material characteristics. Thanks to S.Egli.
*
*===> 04/DEC/1990
*
*        MAXNOD from 2000 to 8000 in GDTR0. Suggestion by S.Egli.
*        Introduced in the 314 correction cradle.
*
*        Bug corrected in GGCLOS, do not touch ISEARC in case
*        of user search (ISEARC.GT.0). Thanks to D.Ward.
*        Introduced in the 314 correction cradle.
*
*===> 05/DEC/1990
*
*        Correction in TWOB, get correctly Tmin and Tmax and
*        their distribution. Thanks to Marko Mikuz.
*        Introduced in the 314 correction cradle.
*
*===> 06/DEC/1990
*
*        Correction in TWOCLU and GENXPT to improve energy
*        conservation. Thanks to S.Egli.
*        Introduced in the 314 correction cradle.
*
*===> 11/DEC/1990
*
*        Generic function LOG introduced in GHFRAK. Generic
*        functions MAX and MIN introduced in several routines
*        in GGEOM (HYPE) and CGPACK. Generic SIN, COS, ATAN and
*        SQRT introduced in CGPACK. Thanks to Federico Carminati.
*        Introduced in the 314 correction cradle.
*
*===> 13/DEC/1990
*
*        GTAU routine to be obsoleted.
*
*        Supporting link for the LSCAN data structure mispelled
*        in several places. Corrected now. Thanks to Isabel Josa.
*        Introduced in the 314 correction cradle.
*
*        Unused variable IPOT eliminated from GHFFER.
*
*        Changes in the interactive part for ATC GKS which does
*        not use FORTRAN unit for metafiles. (routines GXINT and
*        GXGCON). Thanks to W.Koellner
*
*===> 14/DEC/1990
*
*        REAL*8 changed to DOUBLE PRECISION in CGPACK.
*        Life of the D+/- mesons was 100 times too short.
*        Thanks to Cheng He Sheng.
*        Introduced in the 314 correction cradle.
*
*===> 19/DEC/1990
*
*        Modification in GDRAW and GDSPEC to correct a bug.
*        Modifications in GLVOLU not to print a warning message
*        caused by drawing.
*        Introduced in the 314 correction cradle.
*
*===> 20/DEC/1990
*
*        Correction in GRAYLI affecting Rayleigh effect in
*        mixtures.
*        Introduced in the 314 correction cradle.
*
*        Recode of GRAYLI and GHESIG to be safer with respect
*        to compiler optimization.
*        Thanks to Rafi Yaari.
*
*===> 8/JAN/1991
*
*        Corrections in GRGET, GGET, GRIN to set correctly
*        NPART, NTMED, NMATE, NVOLUM. Thanks to Lee Roberts.
*        Introduced in the 314 correction cradle.
*
*        New meaning for ISVOL. If <0 the material is not a
*        a detector but the tracking parameters are the same
*        as in the case ISVOL>0 (sensitive medium).
*        Thanks to Michel Maire.
*        Introduced in the 314 correction cradle.
*
*===> 9/JAN/1991
*
*        TMAXFD not recomputed with IGAUTO=1 if the value given
*        by the user >= 0. Default value reset to 20. instead of
*        60. degrees. Thanks to W.Mueller and M.Maire
*        Introduced in the 314 correction cradle.
*
*===> 10/JAN/1991
*
*        Corrections in GPFIS, GPFISI and GPHYSI for the fotofission
*        and fotoabsorbtion (Giant resonance) of photons on Uranium
*        and plutonium. Thanks to Harm Fesefeld.
*        Introduced in the 314 correction cradle.
*
*        Correction in GPHYSI in case ILOSS=2 and more than one
*        tracking medium with the same material number then a
*        spurious error message was printed. Thanks to Kevin Sparks
*        Introduced in the 314 correction cradle.
*
*        Corrections in GENXPT, TWOCLU, TWOB, CINEMA and EXNU.
*        Thanks to Harm Fesefeldt.
*        Introduced in the 314 correction cradle.
*
*===> 11/JAN/1991
*
*        Correction in GPPART in case of user words. Thanks to Walter
*        Mueller.
*        Introduced in the 314 correction cradle.
*
*        Correction in GSTMED to protect for EPSIL <=0. Thanks to
*        Michel Lefebvre.
*        Introduced in the 314 correction cradle.
*
*===> 14/JAN/1991
*
*        Correction in GRGET to read the banks in the correct
*        division. Thanks to M.Maire.
*        Introduced in the 314 correction cradle.
*
*        Improvement of GRIN, GROUT. New global keywords KINE and
*        TRIG can now be used. Thanks to Federico Carminati
*        Introduced in the 314 correction cradle.
*
*===> 17/JAN/1991
*
*        Correction in GSVOLU, GSDVN, GSDVN2, GSDVT, GSDVT2 not to
*        use uninitialized LINATT, set by GDINIT. This should take
*        care of the problem of empty drawings. Thanks to Simone Giani.
*        Introduced in the 314 correction cradle.
*
*        Correction in GPKINE to print correctly user words and vertex
*        information. Thanks to Kevin Sparks.
*        Introduced in the 314 correction cradle.
*
*        Correction in GPVERT to print correctly user words and vertex
*        information. Thanks to Soren G. Frederiksen.
*        Introduced in the 314 correction cradle.
*
*===> 18/JAN/1991
*
*        Max. GCKING stack size parametrized. The actual value of the
*        stack is not changed (100), but can be redefined in the
*        sequence GCKMAX. Important for LHC studies.
*        Thanks to Harm Fesefeldt.
*
*        Correction in GPRINT to call correctly GPDIGI and GPSETS.
*        Thanks to Michel Maire.
*        Introduced in the 314 correction cradle.
*
*===> 22/JAN/1991
*
*        New facilities of DZDOC introduced in the interactive menu.
*        Patchy flag (DZDOC) to assure backward compatibility.
*        Thanks to O.Schaile
*
*===> 31/JAN/1991
*
*        New facility using the JUMPT package to call user routines.
*        Patchy flag (USRJMP) to assure bacward compatibility.
*        Thanks to F.Carminati
*
*===> 05/FEB/1991
*
*        Correction in GINIT to print the level of the correction
*        cradle applied. Thanks to F.Carminati.
*        Introduced in the 314 correction cradle.
*
*        Correction Cradle 3.14/06 stamped and released.
*
*===> 06/FEB/1991
*
*        Introduction of the AIX370 flag in GEANX and GEANT.
*        Thanks to Roger Howard.
*
*===> 13/FEB/1991
*
*        Protection in GTRACK, GFTRAC and GLTRAC, do not update
*        pointers to cross section banks in the vacuum.
*        Thanks to C.Fuglesang, S.Banerjee, M.Maire.
*        Introduced in the 314 correction cradle.
*
*???>    Backward incompatibility
*
*        GFDETH returns NAMESH as a character*4 array and
*        GFDETD returns NAMESD as a character*4 array.
*        Thanks to M.Maire.
*
*
*===> 14/FEB/1991
*
*        Mods in GPMATE, GPPART, GPTMED, GSDK, GSMATE, GSMIXT, GSPART,
*        GSTMED, GSROTM, GPROTM
*        to warn against replacement of existing objects. Thanks to
*        Steve O'Neale.
*        Introduced in the 314 correction cradle.
*
*        Force load of user routines in GINIT.
*        Introduced in the 314 correction cradle.
*
*===> 15/FEB/1991
*
*        GSDVN to notify correctly when mother volume does not exists.
*        Thanks to Victor Perevoztchikov.
*        Introduced in the 314 correction cradle.
*
*===> 20/FEB/1991
*
*        Correction in GHEINI for the exponent range to be compatible
*        with IEEE machines. Original limits were for IBM/370.
*        Routine LIMDAT eliminated. Thanks to Hans-Jochen Trost.
*        Introduced in the 314 correction cradle.
*
*???>    Backward incompatibility
*
*        Parameter NAMATE returned from GFMATE now a CHARACTER*20
*        variable. Thanks to Michel Maire.
*
*        Monitoring introduced, flag MONITOR. Thanks to F.Carminati.
*        Mods in GXINT, GINIT, GLAST. Routine GEAMON introduced.
*        Introduced in the 314 correction cradle.
*
*===> 22/FEB/1991
*
*        GWORK in GINIT was (DUMMY) called without arguments.
*        GUDTIM was (dummy) called as a routine and it is a
*        function. Thanks to F.Carminati.
*        Introduced in the 314 correction cradle.
*
*        GUVIEW to specify 3rd argument as a character. Thanks to
*        Federico Carminati.
*        Introduced in the 314 correction cradle.
*
*        Call NUCRIN with the right number of arguments (the last
*        was not used anyway). Thanks to Federico Carminati.
*        Introduced in the 314 correction cradle.
*
*===> 04/MAR/1991
*
*        Correction in GMUSIG to avoid division by zero. Thanks to
*        M.Verzocchi.
*        Introduced in the 314 correction cradle.
*
*
*===> 20/MAR/1991
*
*        Corrections in GPRELM, muon e+e- soft radiation. Correction in
*        GDRELM for the ionization energy loss for muons. Thanks to
*        A.Nathaniel. Introduced in the 314 correction cradle.
*
*===> 21/MAR/1991
*
*        Correction in GPHYSI to initialise the JMULOF bank even in
*        the case of vacuum. Electron tracking in vacuum was not taking
*        into account the TMAXFD parameter. Correction in GRANGI to set
*        the range to BIG in case of 0 energy loss tabulated in JLOSS.
*        Correction in GPHYSI to switch off all the mechanisms in the
*        vacuum but decay. The corresponding code in GTRACK has been
*        removed. The previous corrections for IUPD in GLTRAC, GFTRAC
*        and GTRACK have been removed as well. Thanks to D.Ward, R.Brun,
*        F.Carminati
*        Introduced in the 314 correction cradle.
*
*        Corrections in GRANGI and in GTHADR to cure small steps in
*        tracking. Thanks to M.Maire.
*        Introduced in the 314 correction cradle.
*
*===> 25/MAR/1991
*
*???>    Backward incompatibility
*
*        GFPART returns NAPART as a CHARACTER*20 array. Thanks to R.Rui
*        Introduced in the 314 correction cradle.
*
*===> 26/MAR/1991
*
*        GFTMAT will return the correct values for the vacuum, except
*        for hadrons. Thanks to P.Gumplinger.
*        Introduced in the 314 correction cradle.
*
*        Correction in GPCXYZ when the number of mechanisms is greater
*        than 6. Thanks to Y.Iga.
*        Introduced in the 314 correction cradle.
*
*
*===> 27/MAR/1991
*
*        GDXYZ corrected to use HIGZ generic line types. Thanks to
*        M.Maire.
*        Introduced in the 314 correction cradle.
*
*        Correction in the CDF for RZ/IN and RZ/OUT. Thanks to
*        F.Carminati.
*        Introduced in the 314 correction cradle.
*
*===> 04/APR/1991
*
*        For reasons of speed GKS-type lines are used now in
*        GDCXYZ, GDXYZ. Thanks to F.Carminati.
*        Introduced in the 314 correction cradle.
*
*===> 18/APR/1991
*
*        Correction of the logic for TMAXFD. Now is:
*
*            AUTO 1: 0<TMAXFD<20  --> accepted
*                      TMAXFD>20  --> set to 20
*                      TMAXFD<=0  --> set to 20
*
*            AUTO 0:   TMAXFD>0   --> accepted
*                      TMAXFD<=0  --> set to 20
*
*        Thanks to Gerry Lynch.
*        Introduced in the 314 correction cradle.
*
*        Correction in GMULOF. STMIN=BIG for the vacuum. Thanks
*        to F.Carminati.
*        Introduced in the 314 correction cradle.
*
*        Correction in GRFILE, GRIN called with IDVERS=0 and not
*        999. Thanks to F.Carminati.
*        Introduced in the 314 correction cradle.
*
*===> 22/APR/1991
*
*        Corrections in GPRELA, GMULOF, GPHYSI and GTNEXT. Now
*        STMIN is set to 0. for vacuum independly of the AUTO
*        flag. Thanks to M.Maire.
*        Introduced in the 314 correction cradle.
*
*        Format correction in GSDVN. Thanks to Federico Carminati.
*        Introduced in the 314 correction cradle.
*
*===> 23/APR/1991
*
*        Correction in GINIT to initialize LIN via IQTTIN if
*        different from 0. Thanks to F.Carminati
*        Introduced in the 314 correction cradle.
*
*===> 24/APR/1991
*
*        New routines GPLMAT, GPRMAT, GPGKIN from Michel Maire.
*
*===> 25/APR/1991
*
*        geant.metafile is not opened any more by default.
*        The METAFILE command now accepts the metafile name
*        for packages like ATC-GKS which do not use FORTRAN.
*        Thanks to W.Koellner.
*
*===> 01/MAY/1991
*
*        Problem corrected in GHEISH for neutron cascade. A neutron
*        undergoing an hadronic reaction could generate two recoil
*        protons. Thanks to Henk den Bok.
*        Introduced in the 314 correction cradle.
*
*        Correction in GNEXT/GTNEXT. STMIN not used anymore for
*        MANY volumes. Thanks to M.Maire.
*        Introduced in the 314 correction cradle.
*
*===> 02/MAY/1991
*
*        Call to NVETIM routine removed from GHEISHA. Thanks to
*        F.Carminati.
*
*===> 07/MAY/1991
*
*        RNDM interactive command added to set the random number
*        generator seeds. Thanks to F.Carminati.
*
*===> 08/MAY/1991
*
*        New routine GDTRAK introduced to draw a track and delete
*        it from the JXYZ buffer optionally. Thanks to R.Brun,
*        F.Carminati.
*
*        Important correction in GTGAMA, GTNEUT to avoid problems
*        due to machine precision. When a push is made to cross a
*        volume boundary the condition that the change in coordinate
*        be not negligible is imposed. Thanks to F.Carminati, R.Brun.
*        Introduced in the 314 correction cradle.
*
*===> 09/MAY/1991
*
*        Correction in GMGAUS due to G.Lynch. The sigma of the multiple
*        scattering in the gaussian approssimation was too small.
*        Introduced in the 314 correction cradle.
*
*===> 14/MAY/1991
*
*        Corrections in GTELEC, GTHADR and GTMUON to introduce the
*        same boundary correction than in GTNEUT, GTGAMA. Thanks to
*        P.Gumplinger.
*        Introduced in the 314 correction cradle.
*
*        New menu in GXINT to handle FZ files. Thanks to F.Carminati
*
*===> 17/MAY/1991
*
*        Correction in GHFHDN to reset to 0. TOFD. Thanks to
*        R.Rui.
*        Introduced in the 314 correction cradle.
*
*===> 31/MAY/1991
*
*        Correction in GTRACK to allow PREC to decrease and avoid
*        problems when changing particle and position but not medium.
*        Thanks to R.Spiwoks.
*        Introduced in the 314 correction cradle.
*
*===>  4/JUN/1991
*
*        New version of the hidden line removal, thanks to S.Giani.
*        New menu CVOL introduced to cut volumes and new SHIFT
*        command to shift volumes for drawing.
*
*        Correction in GHEISH. The ENP parameters are reset to 0.
*        at the beginning of every event. Thanks to Henk den Bok.
*        Introduced in the 314 correction cradle.
*
*===>  5/JUN/1991
*
*        Correction introduced in GPHYSI. JTM recalculated after
*        bank push. Thanks to K.Sparks.
*        Introduced in the 314 correction cradle.
*
*===> 11/JUN/1991
*
*        Correction introduced in GTELEC, GTMUON, GTHADR. When the
*        energy loss is below machine precision, it is recalculated
*        via the De/Dx table. Thanks to H.Fesefeld and Ralf Spiwoks.
*        Introduced in the 314 correction cradle.
*
*        Booking of the bank IMAT-4 delayed till GPHYSI. Will save some
*        words in case of materials which are not inserted in a tracking
*        medium. Thanks to F.Carminati.
*
*???>    Backward incompatibility
*
*        The argument DMAXMS in the calling sequence of GSTMED, unused
*        in 314, will be called STEMAX and will assume the meaning of
*        maximum step allowed for a particle in the given material.
*        This will only affect people running with AUTO 0. In case of
*        automatic computation of the tracking media parameters, this
*        number will be set to BIG (=1.E10).
*        Thanks to F.Carminati, M.Maire.
*
*===> 12/JUN/1991
*
*        Corrections in GMUNU to allow the usage of a variable number of
*        energy bins. Corrections in PHASP and GMUSIG for very high
*        energies. Thanks to Hans-Jochen Trost.
*        Introduced in the 314 correction cradle.
*
*        Correction in GFLTHE to avoid division by 0. Thanks to
*        Roy Bossingham.
*        Introduced in the 314 correction cradle.
*
*===> 13/JUN/1991
*
*        Correction in GINPGO and GINPCO to check correctly the PHI
*        limits. Thanks to R.Bossingham.
*        Introduced in the 314 correction cradle.
*
*===> 19/JUN/1991
*
*        Correction in GMUNU to calculate correctly the angle of the
*        incoming muon. Thanks to H-J.Trost.
*        Introduced in the 314 correction cradle.
*
*===> 20/JUN/1991
*
*        Corrections in GLANDZ for high energy to avoid numerical
*        problems. Thanks to H-J.Trost.
*        Introduced in the 314 correction cradle.
*
*===> 26/JUN/1991
*
*        All the routines and the commons of the HADRIN/NUCRIN
*        package have been renamed. All routines begin now with
*        GHF and all commons with GCF. This to avoid interactions
*        with the interface with FLUKA. Thanks to F.Carminati
*
*===> 28/JUN/1991
*
*        New algorithm for pushing a view bank. Now push of the
*        max between MORGS and the 25% of the size of the view
*        bank. Thanks to S.Giani.
*
*===> 17/JUL/1991
*
*        Common HIATT of HIGZ removed from GXCONT and workstation type
*        retrieved in GXDZ from Workstation ID. Thanks to O.Couet.
*
*        Modification in GMUSIG to protect against possible division
*        by zero. Thanks to M.Sarris.
*
*        Correction in GDSHOW to recalculate COSPSI and SINPSI. Thanks
*        to M.Verzocchi.
*        Introduced in the 314 correction cradle.
*
*===> 18/JUL/1991
*
*        Correction in GINVOL. When tracking in magnetic
*        field they could return the wrong volume. Thanks to D.Greiner.
*        Introduced in the 314 correction cradle.
*
*===> 24/JUL/1991
*
*        Corrections in GRIN to allow tracking after reading a data
*        structure. Thanks to R.Brun.
*        Introduced in the 314 correction cradle.
*
*        Modifications in GXINT to allow running with the MOTIF
*        user interface. Thanks to R.Brun.
*
*===> 31/JUL/1991
*
*        INT=0 suppressed in GHEISH after nuclear fission.
*        Thanks to Kati Lassila.
*
*===> 01/AUG/1991
*
*        Protection introduced in GMUSIG. Thanks to Jochen Trost.
*        Introduced in the 314 correction cradle.
*
*===> 02/AUG/1991
*
*        Correction in GHSTOP to fix the calculation of the time of
*        flight for stopping particles. Thanks to L.Roberts.
*        Introduced in the 314 correction cradle.
*
*        Correction in GNPGON not to calculate safety for the inner
*        radius when this is 0. Thanks to Andrei Nomerotsky.
*        Introduced in the 314 correction cradle.
*
*        Correction in GNOTRP to return the correct SNXT. Thanks to
*        V.Innocente.
*        Introduced in the 314 correction cradle.
*
*        Correction in GHSTOP not to discard hadrons at rest, they
*        may decay. Thanks to H.Fesefeld.
*        Introduced in the 314 correction cradle.
*
*        Modifications in GHSTOP and GHEISH to handle user defined
*        particles. Thanks to P.Gumplinger.
*        Introduced in the 314 correction cradle.
*
*===> 12/AUG/1991
*
*        Call to GUINTI added in GINTRI to define user commands. Thanks
*        to V.Vercesi.
*        Introduced in the 314 correction cradle.
*
*===> 21/AUG/1991
*
*        Modification in GTMED to print a warning when FIELDM .EQ. 0
*        and IFIELD .NE. 0. Thanks to Federico Carminati
*
*        Modification in GTHADR to set CFLD = BIG in case FIELDM=0.
*        Thanks to V.Vercesi.
*        Introduced in the 314 correction cradle.
*
*===> 27/AUG/1991
*
*        Call to C dummy routine to initialise X11 introduced in
*        GXINT for IBM/VM. Thanks to M.Marquina.
*
*===> 19/SEP/1991
*
*        Rayleigh effect now the default. IRAYL is set to 1 by default
*        in GINIT. Thanks to M.Maire
*
*        Changes in the multiple scattering routines. GMOLI1 and GMOLI2
*        suppressed and GMOLIO used instead. Corresponding changes in
*        GMULTS and GMULOF. The code of GMOLS has been put in line inside
*        GMOLIE and GMOLS is obsolete. GPOISS and GMCOUL modified for
*        performance reasons.
*        Gaussian scattering is now generated according to the Rossi
*        gaussian formula and there is no logarithmic term in the
*        sigma of the gaussian. This gives a distribution where the
*        tails are underestimated, but which is consistent over many
*        steps.
*        Moliere and single Coulomb scattering are chosen according to
*        the value of Omega in the Moliere formula.
*        The new meaning of IMULS is the following:
*
*              IMULS = 0       No multiple scattering
*              IMULS = 1,2     Moliere or single Coulomb scattering
*              IMULS = 3       Gaussian scattering with Rossi formula
*
*        Thanks to G.Lynch.
*
*===> 14/OCT/1991
*
*         New algorithm for calculation of energy loss. The stopping
*         range is now calculated with a parabolical interpolation
*         instead than with a linear one. This gives a dE/dx curve
*         which is a piecewise linear function and not a step function.
*         Thanks to F.Carminati.
*
*         Modification of GNOPG1 taking care of the case of a particle
*         which is very near to the surface and it may appear inside due
*         to machine precision. This could cause the volume to be skipped.
*         Thanks to R.Nierhaus.
*
*         Routine GNPGON rewritten for the calculation of SNXT. The
*         previous routine was returning wrong results. Thanks to
*         R.Nierhaus.
*
*===> 28/OCT/1991
*
*         Modifications in the routines GFLPHI and GFCOOR. The ordering
*         for TUBS in phi was wrong. Volumes could become invisible.
*         Thanks to F.Carminati.
*
*         New routine GFVERT introduced. Retrieves the parameters of a
*         vertex. Thanks to F.Carminati.
*
*===> 29/OCT/1991
*
*         Modification in NUCREC to zero the whole of the PV array
*         to prevent the use of uninitiated variables. Thanks to
*         F.Ranjard. Introduced in the 3.14 correction cradel.
*
*===>  1/NOV/1991
*
*         Useless code commented out in GTRAIN. Thanks to F.Carminati.
*
*         GNCONE gone to double precision. Thanks to J.Toth.
*
*===>  5/NOV/1991
*
*         Variable USERW undefined in the RESULT common now set to
*         UPWGHT from GCTRAK common. Thanks to F.Ranjard.
*
*         Type declarations for GCSTAK completed. Thanks to M.Battle.
*
*===>  5/DEC/1991
*
*         New calculation of the range table. Simpson integration rule
*         used. Modifications in GRANGI, GCOEFF.
*
*===> 10/DEC/1991
*
*         Protection introduced against the reading of a pre-315 data
*         structure. The value of STEMAX is set to BIG in this case.
*         Modification in GPHYSI.
*
*         Better handling of version numbers in I/O operations.
*         Modifications in GPHYSI and GRIN.
*
*===> 08/JAN/1992
*
*         Correction in GLANDZ to avoid gaussian distribution for
*         very thin layers. The version of GEANT 3.13 has been
*         used for this. Thanks to F.Carminati and M.Maire.
*
*         Corrections in GTNEXT, GINVOL and GTMEDI to protect against
*         wrong values of INGOTO when using MANY volumes. Thanks to
*         R.Brun.
*
*===> 15/JAN/1992
*
*???>    Backward incompatibility
*
*         The random numbers seeds are stored at the end of every event
*         in the JRUNG bank at locations 19 and 20. If the JRUNG data
*         structure is read in and the data card RNDM or the interactive
*         command RNDM has not been issued (NRNDM(1), NRNDM(2) <> 0),
*         and if the words 19/20 are not 0, then the random number
*         generator is restarted with these seeds. Thanks to F.Carminati
*
*         The RNDM command now reads the values of the seeds into
*         NRNDM(1) and NRNDM(2) in common GCFLAG. The values 0 0
*         can be used for the random seeds. These values will not
*         alter the current status of the random number generator,
*         but zeroing the variables  NRNDM(1) and NRNDM(2) will
*         allow them to be reinitialised with the values stored
*         in a data structure read from disk. Thanks to F.Carminati
*
*         Routine GREND now needs an integer as input. This is to be
*         consistent with GRFILE.
*
*===> 27/JAN/1992
*
*         The following changes made to the names of FLUKA routines
*         to avoid clashes with ISAJET. Thanks to L.Roberts.
*
*        Type          Original name            New name
*        Routine       DECAY                    FKDECA
*        Routine       FLAVOR                   FKFLAV
*        Routine       SIGINT                   FKSIGI
*
*===> 31/JAN/1992
*
*        Corrections in GMEDIA in case of many volumes to avoid
*        program crash. Thanks to R.Jones.
*
*        Change of logics in GRUN. Now if NEVENT is <= 0, no event
*        is processed. Thanks to B.Lockman.
*
*        Common GCFLAX put in the GCFLAG sequence with the BATCH and
*        NOLOG variable from GXINT. Thanks to B.Cole.
*
*        Problem corrected in GDECAY. The mass of the particle was
*        altered. If the next particle was the same, the mass was
*        not reset to its correct value. Thanks to S.Tonse.
*
*        TIMINT is now really the time left after initiatlization
*        as is specified in the documentation and not the time
*        USED for initialization as it was till now.
*        Thanks to V.Ivanov.
*
*        Corrections in GNPGON. Improvement of code safety under
*        optimization in GNOGO1, GNPGON, GNPCON. Thanks to Y.Iga.
*
*===>  4/FEB/1992
*
*        Corrections in GMULOF and GTHADR. SFIELD takes precedence
*        on STMIN in case IFIELD=1 to make sure that the Runge-Kutta
*        approssimation still works. Thanks to R.Hawkings.
*
*===> 10/FEB/1992
*
*         Correction in GDRAY to improve the precision in the calculation
*         of the angle. Thanks to F.Carminati & P.Lubrano.
*
*         Change in GFKINE. The variable TOFG is not updated any more.
*         Thanks to F.Carminati
*
*         New GMEDIA, GTMEDI, GINVOL to take care of the problems with
*         MANY volumes. Thanks to R.Jones.
*
*===> 14/FEB/1992
*
*         Bug corrected in GNPGON when the particle was exactly on the
*         wall of the last Z section. Thanks to V.Palichik.
*
*===> 17/FEB/1992
*
*         Updated routines GLUND, GLUNDI and new sequences LUDAT1, LUDAT3,
*         LUJETS. The data cards MSTE, KTYP, PMAS, PWID, IDB have been
*         removed because either obsolete or not aplicable. The new code
*         runs with JETSET 7.3 upward. Thanks to F.Carminati, T.Sjostrand.
*
*         Update of GLUDKY to work with 7.3. Thanks to F.Carminati.
*
*===> 23/FEB/1992
*
*         New subroutine GFIN to handle sequential input. Routine GGET
*         has been maintained for backward compatibility.
*         Thanks to F.Carminati.
*
*===> 01/MAR/1992
*
*         New Runge-Kutta integration routine for the tracking in
*         magnetic field. If the result of the stepping is not accurate,
*         the step is divided in 2 parts and the integration repeated
*         and so on. Thanks to V.Perevotchikov.
*
*===> 05/MAR/1992
*
*         Correction in GLTRAC. When the particle fetched has an entry
*         in JKINE the correct vertex number is calculated. Thanks to
*         Y.Foka.
*
*===> 06/MAR/1992
*
*         Correction in GTELEC. The calculated range may be slightly
*         larger than the maximum allowed range due to precision
*         problems and this was leading to very small negative steps.
*         Thanks to R.Brun.
*
*===> 08/MAR/1992
*
*         New GRKUTA, GHELIX, GHELX3 from V.Perevotchikov.
*
*===> 16/MAR/1992
*
*          ENERGY renamed to FKENER in FLUKA. Thanks to F.Carminati
*
*===> 17/MAR/1992
*
*          VERTEX->FKVERT, ZEROIN->FKZERO, ERROR->FKERRO in FLUKA.
*          Thanks to F.Carminati
*
*          Corrections in PBANH. Thanks to M.Sasaki.
*
*===> 19/MAR/1992
*
*          Corrections in GFIN, GFOUT and GRIN, version 0 is now allowed
*          for a data structure in I/O. Thanks to B.Cole.
*
*          Modification in GBREME to improve precision for small angles.
*          Thanks to F.Carminati.
*
*          Inlining of rotation routines in GFTRAC, GINVOL, GTMEDI,
*          GTNEXT, GMEPOS, GMEDIA. Thanks to D.Kryn, F.Carminati.
*
*          New routines GDLENS, GDPLST and GDPRTR in the drawing package.
*          New version of GDTREE. Thanks to S.Giani.
*
*          Streamlining of GDTOM, GINROT, GITRAN, GMTOD, GRMTD, GRMUL,
*          GROT, GTRMUL and GTRNSF. Thanks to F.Carminati.
*
*===> 23/MAR/1992
*
*???>    Backward incompatibility
*           Tracking medium name in GFTMED changed in CHARACTER variable.
*           Thanks to V.Perevotchikov.
*
*           IMPULS renamed to FKIMPU
*           DRES   renamed to FKDRES
*           ERUP   renamed to FKERUP in FLUKA. Thanks to F.Carminati.
*
*===> 24/MAR/1992
*
*       COSI entry point renamed to FKCOSI
*       POLI             renamed to FKPOLI in FLUKA. Thanks to F.Carminati
*
*       Changes in GCOMP to increase the precision of the rotation
*       of the scattered photons and electrons in the reference
*       frame of the incoming particle. Thanks to F.Carminati
*
*       New version of the routines GNOPG1 and GNPGO1. Speed up of
*       a factor two or more achieved. Thanks to F.Carminati,
*       M.Roethlisberger.
*
*===> 29/MAR/1992
*
*       Improvement of the routines GMEDIA, GTMEDI and GINVOL. A better
*       use is made of the variable INGOTO. Thanks to R.Jones, F.Carminati.
*
*===> 01/APR/1992
*
*       Modification in GTRACK. If a particle tries for more than 5 times
*       consecutively to exit a volume, the precision used for tracking is
*       multiplied by 5 and so on every fifth attempt. Thanks to
*       F.Carminati.
*
*       Modification in FLUFIN to normalise in double precision the
*       direction cosines given to FLUKA. Thanks to A.Ferrari
*
*       Modification in GLANDZ to avoid peaks in energy loss for very
*       light materials. Thanks to F.Carminati
*
*       New routine GETVER to crack the title sequence and to return
*       the correct version number. Thanks to F.Carminati and M.Maire.
*
