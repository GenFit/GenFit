*
* $Id: v_314.F 220 2007-11-19 16:08:06Z rdm $
*
* $Log: v_314.F,v $
* Revision 1.1.1.1  2002/06/16 15:18:35  hristov
* Separate distribution  of Geant3
*
* Revision 1.1.1.1  1999/05/18 15:55:15  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:19:47  cernlib
* Geant
*
*
#include "sys/CERNLIB_machine.h"
#include "pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.18  by  S.Giani
*-- Author :
C
C
C           ************************************************
C           *                                              *
C           *         G E A N T  version  3.14.16          *
C           *                                              *
C           *           Official  R e l e a s e            *
C           *       ==============================         *
C           *                                              *
C           ************************************************
C
C
C
*    The version 3.14 of GEANT is being released. Apart from a few
*    features, reported below, the new version is compatible with the
*    previous version 3.13. Substantial developments have taken place,
*    in particular in the physics and tracking areas.
*
*    In addition to the GEANT team
*     (R.Brun:CERN/CN-AS,  F.Bruyant:CERN/ECP-PI and M.Maire:LAPP),
*    many people have contributed to the new version, especially
*     - a large number of users of the previous version who reported
*       their experience, found bugs and suggested several improvements,
*     - the guinea-pigs of 3.14 who gave us essential comments during
*       our attempt to automatize computation of the tracking parameters.
*                 ** A.Givernaud (UA1), F.Nessi, V.Vercesi(UA2/LHC)
*                 ** The LHC proto collaborations
*     - H.Fesefelt(Aachen) has provided a new version of the GHEISHA
*       package with help from N.Van Eijndhoven (CERN/CN-AS)
*     - H.J Trost from ANL has reported problems in the muon-nuclear
*       interactions routines and provided the relevant corrections.
*     - P.Pedroni (Pavia) has implemented a new interface called
*       by GHEISHA for low-energy hadronic processes (See NUCRIN below).
*     - the contributors to the geometry package:
*            for the introduction of new shapes
*                ELTU by A.Solano (ZEUS)
*                HYPE by M.Corden (ALEPH and SSC)
*                CTUB by A.McPherson (CEBAF)
*            for systematic investigation of problems with the old shapes
*                R.Nierhaus (CERN/CN-AS)
*     - S.Egli (H1) has proposed an automatic optimisation for the
*       geometry at initialisation time.
*     - G.Lynch from Berkeley has investigated the multiple-scattering
*       various strategies and proposed a new algorithm.
*     - A.Rotondi and P.Montagna have proposed a new technique for the
*       fast generation of Vavilov distribution.
*       K.S.Koelbig (CERN/CN-AS) has implemented a new routine GVAVIV
*       based on their work.
*     - J.Salt (CERN/CN-AS) has implemented the graphics interface to the
*       CG package from Serpukhov with the help of E.Chernaev.
*     - The IBM team (C.Guerin, M.Roelisberger) have investigated how
*       to speed-up the program. Their work has been coordinated
*       by F.Carminati (CERN/CN-AS).
*
*     GEANT3.14 documentation
*     =======================
*      The printing of a new manual describing the new version is
*      scheduled for the end of this year. The CERN Program Library
*      will not distribute anymore the old document (version 3.11).
*
*
*     Important notice to GEANT users
*     ===============================
*
*      Following the reorganisation of the CERN research divisions
*      in July 90, R.Brun is now in charge of the Application Software
*      group (AS) in the CN division and F.Bruyant is in charge of the
*      Production support and computers Infrastructure group (PI) in the
*      ECP division. They will nevertheless continue, with M.Maire, to be
*      actively involved in the development of GEANT. Federico Carminati
*      is the coordinator of a simulation software unit in the CN-AS group.
*      Users are strongly recommended to address their questions,etc
*      directly to him (email: FCA@CERNVM.CERN.CH on BITNET). In particular,
*      feedback from users making comparisons with real data will be most
*      appreciated.
*
*      The simulation team is preparing the ground for the next version
*      of GEANT. In view of the proposed new accelerators, the following
*      items are considered with high priority:
*        - Parametrization techniques. A survey of the various methods
*          used in current experiments has been done and the implementation
*          of a new algorithm based on the GSCAN geometry + GFLASH (H1)
*          is in progress.
*        - Improvements in the geometry package. General shape definition,
*          surface based algorithms.
*        - Detector data structure and data base. Interface with CAD systems.
*        - Parallelism (at event level and below).
*
*
*......................................................................
*
* *** Compatibility with version 3.13
*     ===============================
*     - If COMMON blocks GCMULO or/and GCJLOC were included in the
*       user code, the new GEANT sequences GCMULO and GCJLOC must be
*       inserted and the code recompiled.
*     - Initialisation data structures saved with the previous versions
*       cannot be read by the new version, because the binning for
*       the cross-sections and energy loss tables has been changed.
*     - The GCPHYS variables SOLOSS,STLOSS,SOMULS,STMULS are no more
*       defined (see comments below)
*     - The GCTRAK variable IDECAD is replaced by IGAUTO
*
*
* *** MAIN CHANGES IN THE TRACKING PROCEDURES
*     =======================================
*
*    The tracking control routines GTGAMA,GTELEC,GTHADR,GTNEUT,GTMUON
*    have been largely rewritten to reflect the changes to the energy-loss
*    and multiple scattering processes.
*
*
* *** The ENERGY RANGE of the cross section and energy loss tables can
*     be fixed by the user with the new data card :
*             'ERANG'   EKMIN  EKMAX    NKBIN
*     which defines nkbin bins from Ekmin to Ekmax in a logarithmic scale.
*     The default is, as before, 90 bins from 10 Kev to 10 Tev but in
*     logarithmic scale. NKBIN must be 50<NKBIN<200.
*
*     WARNING 1 : as a consequence the common GCMULO has been changed.
*     User applications referencing GCMULO must be recompiled.
*
*     WARNING 2 : as a consequence the ZEBRA data structure JMATE
*     contains now cross-sections and energy loss tables based on the new
*     energy range. Structures saved by previous versions of GEANT cannot
*     be used by the new version.
*
*     WARNING 3 : changing the 'ERANG' data card arguments requires to
*     recreate the initialisation data structures.
*
*     WARNING 4 : for Hadrons, GEANT tabulates the energy-loss tables
*     for a proton-kinetic energy equivalent. As a consequence, the value
*     given for EKMAX must be at least 7 times the maximum kinetic energy
*     of particles to be tracked (ratio proton/pion mass).
*     For LHC/SSC simulations, we recommend ERANG 1.E-5 1.E+5 100
*     The recommended value for NKBIN is 10*LOG10(EKMAX/EKBIN).
*
*
* *** The MEAN STOPPING RANGE of a charged particle (STOPMX) is computed and
*     tabulated as a function of the kinetic energy (routine GRANGI).
*     There are tables for electrons, muons, protons. The others particles are
*     scaled from the proton table.
*     As for the other tabulated quantities, the Stopping Range can be displayed
*     with the routines GPLMAT and/or GPRMAT (keyword 'RANG').
*
*     During the tracking, the MEAN ENERGY LOSS over a given step (DEMEAN) is
*     computed from the difference of the stopping range before and after the
*     transport of the particle.
*
*
* *** ENERGY LOSS STRAGGLING : when ILOSS=2 (default) the Landau/Vavilov
*     fluctuations are applied to the mean energy loss over the step.
*     The VAVILOV distribution has been improved so that the fluctuations are
*     less dependent on the step size.
*
*     ILOSS=1 is now equivalent to ILOSS=3 i.e. restricted fluctuations are
*     applied together with delta-ray production (see routine GLANDZ).
*     The user who wish to inhibit the energy fluctuations must set ILOSS=4
*     The latter option has been kept for debug purpose only.
*
*
* *** MULTIPLE SCATTERING.  The calculation of the step size (SMULS) is new :
*     Bethe criterion taken into account, limitation of step size at low
*     energy (see routine GMULOF).
*
*     New routine GMULTS : depending on the step size the Moliere distribution,
*     or its gaussian approximation, is used.
*
*     The routines GMOLI, GMOLIE (ex Gmol), GMOLS have been slightly modified.
*
*     In the routine GMGAUS (ex Gmuls), we use a new sigma (G. Lynch LBL-28165).
*     The lateral displacement has been deleted. As a consequence the tracking
*     parameter DMAXMS is not used anymore for the control of the step size.
*     DMAXMS.LE.0 forces no multiple scattering at all in the medium.
*
*     New routine GMCOUL : for very small steps single Coulomb scatters are
*     generated instead of Moliere distribution.
*
*     By default (IMULS=1) GEANT will select automatically the algorithm to use
*     (GMCOUL, GMOLIE or GMGAUS).
*
*
* *** NEW ROUTINE GMULOF. The step for energy loss (SLOSS), the step for
*     multiple scattering (SMULS), and the step  for the curvature in a magnetic
*     field (SFIELD) can be precomputed and tabulated, at least for electrons
*     and muons. The routine GMULOF, called from GPHYSI, tabulates:
*                  SMULOF = MIN ( SMULS , SLOSS , SFIELD )
*     which is the effective step due to the "continuous processes".
*     SMULOF can be displayed with the routines GPLMAT/GPRMAT, keyword 'STEP'.
*
*
* *** Multiple scattering and energy loss computations are systematically
*     applied at each step during tracking. Therefore the GCPHYS variables
*     SOLOSS, STLOSS, SOMULS, STMULS are OBSOLETE. However, for backward
*     compatibility with 3.13, STLOSS is set equal to STEP.
*
*
* *** TRACKING PARAMETERS STMIN and DEEMAX.
*
*     The meaning of STMIN is the following : at low energy (below 1 Mev) the
*     multiple scattering condition can induce a very small step SMULOF.
*     On the other hand the Stopping Range is also small. Therefore, if the
*     Stopping Range is smaller than STMIN, the constraint SMULOF is ignored.
*     The exact condition is :
*               IF (SMULOF.LE.STMIN)   SMULOF = MIN ( STOPMX, STMIN )
*     STMIN is no more than an accelerator factor for the stopping particles.
*
*     The DEEMAX parameter remains the main tracking parameter. It governs the
*     precision of the tracking by limitating SMULOF.
*
*
* *** AUTOMATIC COMPUTATION of STMIN and DEEMAX.
*
*     By default Geant3.14 overwrites STMIN and DEEMAX. The STMIN default value
*     correspond to a stopping range of 200 Kev above CUTELE.
*     The default for DEEMAX follows the following algorithm:
*       - For non-sensitive volumes (ISVOL=0) DEEMAX is set to 0.25 for
*         materials with a radiation length x0<2cm
*         and DEEMAX=0.25-0.2/sqrt(x0) for other materials.
*       - For sensitive volumes (ISVOL.NE.0) DEEMAX=0.2/sqrt(x0)
*     These values have been tuned empirically on a variety of setups.
*
*     A new data card 'AUTO' has been implemented :
*     Setting 'AUTO' 1 is equivalent to NO data card, i.e. automatic computation
*     Setting 'AUTO' 0 : NO automatic computation, EXCEPT if STMIN and/or DEEMAX
*     has been given a negative value by the user.
*
*     WARNING : The default values above have been found reasonable for any kind
*     of medium. The unexperienced user is invited to start with automatic
*     computation. Please check the actual parameters by calling GPRINT ('TMED'
*     after the call to GPHYSI.
*
*     We STRONGLY recommend to always run in AUTO mode (default). The AUTO
*     mode makes GEANT a predictive tool if all parameters are automatically
*     computed by the system as opposed to tuning data and Monte Carlo
*     via the tracking parameters.
*
* *** GPHYSI
*     The routine GPHYSI has been improved to include additional protections
*     and to take into account the changes due to the new energy-range tables.
*     It must be noted that is mandatory to always call GPHYSI at the end
*     of the initialisation phase even when initialisation data structures
*     are read from a file.
*
* *** GPIONS
*     This new routine (a complement to GPART) defines a subset of the most
*     common "stable" elements in the nature.
*     GPIONS can be called at the initialisation stage after GPART. It creates
*     particles with GEANT identifiers 61 to 112.
*     GPIONS has been written for the heavy ions experiments in view of an
*     interface with the program FRITIOF which is being developed by
*     P.Gorodetzky (Strasbourg).
*
* *** GDEBUG
*     This new user callable routine from GUSTEP will take the following actions
*     if the flag IDEBUG=1:
*
*     IF(IDEBUG.NE.0) THEN
*       IF((ISWIT(2).EQ.1).OR.(ISWIT(3).EQ.1)) CALL GSXYZ ! store point in JXYZ
*       IF (ISWIT(2).EQ.2) CALL GPCXYZ ! step by step printed debug
*       IF (ISWIT(1).EQ.2) CALL GPGKIN ! list of particles generated during step
*       IF (ISWIT(2).EQ.3) THEN
*          IF(ISWIT(4).EQ.3.AND.CHARGE.EQ.0.)RETURN
*          CALL GDCXYZ    ! interactive drawing of trajectories
*        ENDIF
*     ENDIF
*
* *** New user callable routine GBIRK to take into account BIRK's factors.
*     ====================================================================
*    This new routine is callable from the GUSTEP routine. EDEP=DESTEP
*    if no Birk factors have been given via GSTPAR.
*      SUBROUTINE GBIRK(EDEP)
*
*     *** apply BIRK's saturation law to energy deposition ***
*     *** only organic scintillators implemented in this version MODEL=1
*
*     Note : the material is assumed ideal, which means that impurities
*            and aging effects are not taken into account
*
*     algorithm : edep = destep / (1. + RKB*dedx + C*(dedx)**2)
*
*     the values of RKB and C can be entered via :
*
*     call gstpar(imate,'BIRK1',value) to set the model (value= 1. or 2.)
*     call gstpar(imate,'BIRK2',value) to set RKB
*     call gstpar(imate,'BIRK3',value) to set C
*
*     the basic units of the coefficient are g/(Mev*cm**2)
*     because the de/dx is obtained in Mev/cm
*
*     exp. values from NIM 80 (1970) 239-244 :
*
*     RKB = 0.013  g/mev*cm**2  and  C = 9.6e-6  g**2/(Mev**2)(cm**4)
*
*
*
*
* *** GEANH news
*
*
*     a serious bug in GMUNUI has been reported by
*     Hans-Jochen Trost from ANL:
*       : convert millibarns to cm**2 by a factor 10**-27
*         and obtain therefore NA * 1 millibarn/cm**2 a
*         factor 10 smaller than in the code. Hans has verified
*         that the original code gives to high an energy loss
*         due to muon-nuclear interactions on iron by a fac-
*         tor of about 5-10, comparing to W.Lohmann et al.,
*         CERN 85-03 whose predictions have recently been
*         verified by CCFR data up to 1.2 TeV.
*     Hans has also introduced protections in GMUNU.
*
*    A new version of the GHEISHA package has been introduced by
*     H.Fesefelt (contact him for details).
*
*     The TATINA package is considered obsolete. Default routines
*     GUPHAD/GUHADR call GHEISHA. We are planning to remove the code
*     of TATINA in the coming versions of GEANT.
*
*
* *** NUCRIN: A new model for hadronic showers at low energies
*     ========================================================
*
*
*    NUCRIN (see ref.1) simulates hadron-nucleus ( A > 3) inelastic
*    interactions from a few MeV/c up to about 4.5 GeV/c laboratory
*    momentum of the incoming particle.
*    NUCRIN is automatically called by the GHEISHA routine GHEISH
*    when the flag IHADR=3 (set by data card HADR).
*    It is assumed that these reactions  are the superimposition of three
*    basic processes:
*
*    (a) inelastic collision of the projectile hadron (allowed particles are :
*        p,pbar,n,nbar,pi0,pi+,pi-,k+,k-,k0,k0bar,lamda0,lamda0bar,sigma+,
*        sigma-,sigma0) with a target nucleon in the nucleus.
*        This interaction is simulated,taking into account of the nucleon
*        Fermi momentum,  using HADRIN (see ref.2) program.
*        The corresponding physical model is based on the experimental evidence
*        that, in the selected momentum range,the inelastic cross section shows
*        the typical threshold and resonance behaviour of meson production:
*        the primary hadron-nucleon system is excited to an isobaric state
*        which then decays into hadrons or other resonances.
*
*        If the interacting nucleus is hydrogen,HADRIN can also be used
*        in a separate way to simulate hadron-proton reaction.
*
*    (b) induced intranuclear cascade with resulting proton and neutron
*        emission;
*
*    (c) nuclear evaporation and deexcitation from residual nucleus. At the
*        output the total energy available for these processes is given as
*        "excitation energy".
*
*    The mean excitation and cascade energies and the average multiplicities
*    of cascade particles are parametrized, according to experimental
*    distributions.
*    In each event their value are sampled from gaussian distribution:
*    if they fall in the permitted kinematical region, energy and types of
*    cascade nucleons are calculated and the remaining energy is assigned to
*    the incoming particle.
*    For hadron-nucleon interactions all relevant kinematic variables are
*    Lorentz-trasformed into the target nucleon rest system. If interaction
*    kinetic energy is greater than the total available collision energy
*    a new Fermi momentun is sampled, otherwise an event is generated with
*    HADRIN,in which decays modes of 107 particle and resonances into about
*    450 different channels are tabulated and outgoing particle directions
*    and momenta are chosen to reproduced experimental momentum transfer
*    distributions.
*    Final state particles kinematical variables are transformed back into
*    laboratory system; reaction and sampled event energies are again compared:
*    if their difference is negative, energy is not conserved and generation
*    has to be started once more with a new Fermi momentum sampling or if it
*    is,on the contrary positive, particle momenta and energies are corrected
*    to reach conservation.
*    The sampled events conserve the energy, the momentum, the electric and
*    baryonic charge  and the strangeness.
*    NUCRIN and HADRIN are initialised, by default, before event generation,
*    with a call to subroutines HADDEN and CHANWX which estabilish internal
*    weight tables and decay channels.
*
*    -----------------------------------------------------------------------
*    (1) K. Hanssgen, J. Ranft , Comp. Phys. Comm. 39, 53 (1986)
*    (2) K. Hanssgen, J. Ranft , Comp. Phys. Comm. 39, 37 (1986)
*
*
* *** GEOMETRY PACKAGE: New shapes and many improvements
*     ==================================================
*
*    Automatic optimisation of the geometry structure:
*    A new data card OPTI has been introduced (S.Egli H1).
*      OPTI -1 : disable optimisation
*      OPTI  0 : only volumes GSORDered are optimised (as in 3.13)
*      OPTI  1 : volumes GSORDered are optimised along the axis
*                specified. All the other volumes are automatically
*                optimised along the best axis (1 to 7).
*      OPTI  2 : All volumes are optimised along the best axis.
*                Volumes for which GSORD was called are also optimised.
*      The default value for OPTI is 0.
*      In case OPTI >0, the result of the optimisation is printed.
*      The automatic optimisation is done at initialisation time by
*      a new routine GGORDQ called by GGCLOS.
*
*
*    Most of the geometry routines have been revisited and consolidated.
*    The following new shapes are available.
*
*      'ELTU'    is a cylinder with an elliptical section.
*                It has three parameters: the ellipse semi-axis in X,
*                the ellipse semi-axis in Y and the half length in Z.
*                Given the equation of the conical curve:
*                     X**2/A**2 + Y**2/B**2 = 1,
*                describing the volume,then:       PAR(1) = A
*                                                  PAR(2) = B
*                                                  PAR(3) = DZ
*                ELTU is not divisible.
*
*      'HYPE'    is a hyperbolic tube, ie the inner and outer surfaces
*                are hyperboloids, as would be formed by a system of
*                cylindrical wires which were then rotated
*                tangentially about their centres.  The 4 parameters
*                are the inner and outer radii, the half length in z,
*                and the "stereo angle" theta in degrees, such that
*                the hyperbolic surfaces are given by
*                r**2 = (z*tan(theta))**2 + (r at z=0)**2
*
*      'CTUB'   (for cut tube) is a TUBS whose end planes are not
*               perpendicular to the z axis. It has 11 parameters :
*               the 5 of the TUBS shape plus the components of the normal
*               to the end plane at the lower z (LXL,LYL,LZL) and
*               those at the higher z (LXH,LYH,LZH).
*               DZ means the half length in z for x = y = 0
*
*
*     Bugs have been fixed in the routines GNPCON, GNPGON and GNOPGO.
*
*     The new version of GNOTRP requires an extended parameter array.
*     In addition to the 11 specified parameters (of which 4 are
*     modified in subroutine GSVOLU or GSPOSP), the coefficients of
*     the implicit normalized plane equation for the 6 surfaces
*     of the hexahedron are stored.
*
*     IMPORTANT NOTE concerning the TRAP shape
*     ========================================
*     The Geant documentation describes the Geant shape TRAP as follows:
*
*          TRAP is a general trapezoid, i.e. one for which the
*          faces perpendicular to z are trapezia and their
*          centres are not at the same x, y. It has 11
*          parameters: Dz the half length in z, Th & Phi the
*          polar angles from the centre of the face at z=-Dz
*          to that at z=+Dz, H1 the half length in y  at
*          z=-Dz, LB1 the half length in x at z=-Dz and y=low
*          edge, LH1 the half length in x at z=-Dz and y=
*          high edge, Th1 the angle w.r.t. the y axis from
*          the centre of the low y  edge to the centre of the
*          high y edge, and H2, LB2, LH2, Th2 the
*          corresponding quantities to the 1s but at z=+Dz.
*
*     This seems to describe a general hexahedron with 3 constraints:
*     2 constraints follow from the fact that two faces are "trapezia".
*     (twice 2 edges parallel).
*     The 3rd constaint is that two faces are parallel, namely the "trapezia"
*     faces are both perpendicular to the Z-axis.
*     We will assume that shape TRAP is a hexahedron with 3
*     constraints and direct our attention to the degrees of freedom
*     of such a shape.
*     The shape has 8 vertices and therefore 24 coordinates.
*     If we first consider a volume with 8 vertices and 6 surfaces,
*     but do not make the assumption that the surfaces are plane,
*     we see that this shape has 18 degrees of freedom.
*     We loose 3 coordinates because of the translational invariance
*     of the shape and 3 coordinates because of its rotational
*     invariance.
*     If we now assume that the shape is a hexahedron, that is
*     bounded by plane surfaces, we have 6 constraints, one for each
*     surface, and our hexahedron has 12 degrees of freedom.
*     Considering the 3 constraints mentioned in the beginning,
*     we conclude that our shape trapezohedron has 9 degrees of freedom.
*     It is however described by 11 parameters.
*     Therefore we must either drop our assumption that shape TRAP is
*     a hexahedron, that is bounded by parallel surfaces, or we must
*     request that the user specifies the 11 parameters with certain
*     constraints.
*
*     To check that the user respected the constraints, we check the
*     coplanarity of the faces during the specification phase of shape TRAP.
*     We know the vertex coordinates, and we have the indices to the
*     vertices for each face.
*     Assuming that a face is tetrahedron, we compute its volume.
*     We divide by the surface of the base triangle, and get a measure
*     for the coplanarity of the face, which is actually a distance.
*     A warning message is printed in case of no-coplanar faces.
*
*    GSORD problem
*    =============
*    A bug has been found in GTNEXT (and alike) which is induced by a bug of
*    logic in GSORD/GGORD. User calls to GSORD, with ordering axis 4 (Rxy)
*    or 5 (Rxyz), may cause problems when the ordered contents are such that
*    one can jump from a given content to another one without crossing a
*    content which, along the given axis, occupies a position in between the
*    start and the end contents : e.g. coaxial TUBES with different Z-lengths
*    should not in general be ordered by GSORD along the axis IAX=4 (Rxy).
*    However, part of the information provided by such calls can still be
*    used, in the static context of GTMEDI for instance. Therefore, the
*    following convention has been introduced: If the user is sure that the
*    contents are positioned in such a way that the anomaly mentioned above
*    cannot occur, the call to GSORD can be modified by using IAX=14 (instead
*    of 4), or 15 (instead of 5), in which case the ordering techniques will
*    also be used in the dynamic context of GTNEXT. In case of doubt, the
*    user has better to keep the old code, with IAX=4 or 5.
*
* *** DRAWING PACKAGE: Interface to the CG package
*     ============================================
*
*    An interface to the CG (Combinatorial Geometry) package written
*    at Serpukhov by E.Cernaev et al is now available. The new package
*    is automatically called if the option 'HIDE' is selected. eg.
*
*         CALL GDOPT('HIDE','ON')  in a Fortran program
*         DOPT HIDE ON/OFF  in the interactive version
*
*    This new facility includes a hidden line and surface algorithm
*    which permits nice 3-D views of a detector.
*    In the frequent case of hermetic 4 PI detectors, a facility
*    to remove a box (The Cutting BOX) is also implemented. The Cutting BOX
*    specifies a region of the detector which must be Cut to see inside.
*    A new interactive command CBOX is available to specify the box limits.
*
*    The CG system is part of the GEANG file (See Patches CGCDES,CGPACK)
*    To activate the CG package, +USE,CG,*GEANG.
*
*
* *** The SCAN geometry
*     =================
*   A new Patch,GSCAN has been introduced on a provisional basis in the GEANG
*   file. See discussion about SCAN below.
*
* *** GUPARA: Parametrization interface
*     =================================
*   A new FFREAD data cards PCUT can be used to set parametrization cuts.
*   The first argument of the PCUT card is a integer flag which turns or
*   on off the parametrization mechanism. If the parametrization is turned on
*   and a particle falls below one of the 5 cuts specified by the PCUT card
*   (similar in kind to the cuts specified by the card CUTS), then the
*   routine GUPARA is called and tracking of the particle is abandoned. This
*   mechanism is provided for applying parametrization schemes which
*   replace a particle by a parametrized shower when it falls below
*   a certain threshold.
*
*.............................................................................
*
*
* *** GXINT: Interactive version
*     ==========================
*
*      New menu FORTRAN: CALL,FILE,CLOSE,FORTRAN
*      New menu HISTOGRAM: FILE,LIST,PLOT,DELETE,LEGO,HRIN,HROUT,PUT,GET
*                          ZONE,SET,OPTION,NULL
*      New menu PICTURE: FILE,LIST,DELETE,SCRATCH,PLOT,RENAME,IZOUT,IZIN,IGSET
*      New menu SCAN: PHI,THETA,SLIST,VERTEX,SFACTORS,LSCAN,HSCAN
*      New menu PHYSICS: ANNI,BREM,COMP,DCAY,DRAY,HADR,LOSS,MULS,MUNU,PAIR,
*                        PFIS,PHOT,RAYL,CUTS,PHYSI
*
*      The menus FORTRAN,HISTOGRAM and PICTURE are subsets of the similar
*      menus in the PAW system.
*
*      The menu PHYSICS gives the possibility to modify the run conditions.
*      In case physics conditions are changed (LOSS,DRAY,MULS,CUTS) it is
*      mandatory to call the command PHYSI (which calls GPHYSI) to recompute
*      the cross-section and/or energy loss tables.
*
*      The menu FORTRAN is similar to the FORTRAN menu of PAW. It contains
*      in addition a new very important command FORTRAN which gives the
*      possibility to describe the geometry (UGEOM) in a Fortran routine
*      that can be edited interactively with the local editor and also
*      executed interactively under the control of the Fortran interpreter
*      COMIS.
*
* ***  GEANT >FORTRAN  FNAME
*
*      The routines in the file FNAME will be compiled by COMIS.
*      If routines with names: UGEOM,GUKINE,GUOUT,UGLAST are found,
*      then they will be automatically called by GXINT instead of
*      the routines with the same names compiled with the standard
*      Fortran compiler and linked with the application.
*      The user callable routines from the GEANT library as well as
*      routines from PACKLIB (HBOOK,HPLOT,HIGZ,ZEBRA) may be called
*      from these user routines. All GEANT common blocks may be
*      referenced.
*      In case where the routine UGEOM is called several times,
*      it is important to DROP all the initialisation data structures
*      JVOLUM,JMATE,JTMED,etc already in memory by using the routine GIDROP.
*
*       Example of an interactive session where the routine UGEOM is modified:
*
*         GEANT > Edit ugeom.for
*         GEANT > Fortran ugeom.for
*         GEANT > Call GIDROP
*         GEANT > Call UGEOM
*         GEANT > Dtree
*         GEANT > Edit ugeom.for
*         GEANT > Fortran ugeom.for
*         GEANT > Call GIDROP
*         GEANT > Call UGEOM
*         GEANT > Dtree
*
*      If FNAME='-', calls to user routines is reset and standard
*      routines called instead.
*
* *** Interface to CG
*    The command DOPT has a new option HIDE (DOPT HIDE ON/OFF)
*    Type DOPT without parameters to get the list of all currently
*    available options.
*    When this option is ON, the subsequent graphics commands DCUT/DRAW
*    will invoke the CG system for hidden line/surface removal.
*    This algorithm requires a lot of memory and time. It is recommended
*    to set the visibility attributes (SATT SEEN 0) for many of the
*    volumes in case the command aborts for lake of memory or time.
*    This option can also be used with the view banks mechanism (DOPEN)
*
*    New command CBOX to specify the boundaries of the cutting box.
*
* *** PLMAT
*
*   The existing command PLMAT offers the new possibility to plot
*    various physics parameters (cross-sections, energy-loss tables,etc)
*   in graphics format (via HPLOT) if MECAN=ALLG
*   The Keywords 'STEP' or 'RANG' may also be specified to produce
*   an alphanumeric output of the step-size and energy-range tables.
*
*        PLMAT  IMATE IPART MECAN [ IDM ]
*
*       IMATE      I 'Material number'
*       IPART      I 'Particle number'
*       MECAN      C 'Mechanism'
*       IDM        I 'ID mode option' D=0
*
*
* *** New menu SCAN. The SCAN geometry
*     ================================
*
*    This new menu contains various commands for an interactive interface
*    to the SCAN geometry
*    The SCAN geometry algorithm has been designed as a tool to improve
*    the tracking speed. This new facility still requires substantial
*    developments in view of the new parametrisation algorithms which
*    are developed in collaboration between the CN/AS group and the
*    LEP/HERA/LHC/SSC and other interested groups.
*    The SCAN facility is being introduced in the version 3.14 on a trial
*    basis to familiarise potential users with the concept.
*    The SCAN geometry data structure JSCAN is automatically generated
*    either by calling the GSCAN routine in the PATCH,GSCAN of GEANG
*    or interactively by using the commands in the new menu SCAN.
*    Starting from the normal geometry data structure created by GSVOLU,
*    GSPOS,GSDVN,etc, the detector may be divided into a simpler geometry
*    structure (theta,phi) or (eta,phi).
*    Geantinos are tracked starting from a VERTEX position through
*    the NPHI,NTETA divisions. For each division, the SCAN procedure
*    will insert into the JSCAN data structure the following information
*    for every main detector component specified in the SLIST command
*    in the spherical R direction:
*      Total number of radiation lengths up to entry in each R
*      Total number of absorption lenghts
*      Detector identifier
*    When the interactive command TRIGGER is entered, the number of
*    Geantinos specified as parameter will be tracked. In case the
*    data structure JSCAN is not empty, the program will automatically
*    start with the first PHI,TETA division not yet filled. As the
*    number of Geantinos to be tracked can be very large (depending
*    on the granularity) this gives the possibility to fill the JSCAN
*    data structures in several passes.
*
*
*     New menu SCAN: PHI,TETA,SLIST,VERTEX,SFACTORS,LSCAN,HSCAN
*
* ==>   /SCAN/PHI  NPHI [ PHIMIN PHIMAX ]
*
*        NPHI       I 'Number of PHI divisions' D=90
*        PHIMIN     R 'Minimum PHI in degrees' D=0
*        PHIMAX     R 'Maximum PHI in degrees' D=360
*
*        To specify number of divisions along PHI.
*
*
* ==>   /SCAN/TETA  NTETA TETMIN TETMAX [ DIVTYP ]
*
*        NTETA      I 'Number of TETA divisions' D=90
*        TETMIN     R 'Minimum value of TETA' D=0
*        TETMAX     R 'Maximum value of TETA' D=180
*        DIVTYP     I 'Type of TETA division' D=1 R=1:2
*
*        To specify number of divisions along TETA.
*        If DIVTYP=1 divisions in degrees following the THETA angle.
*        If DIVTYP=2 divisions in pseudo-rapidity ETA.
*
*
* ==>   /SCAN/SLIST  LIST
*
*        LIST       C 'List of master volumes'
*
*        Only boundary crossings of volumes given in LIST
*        will be seen in the SCAN geometry.
*
*
* ==>   /SCAN/VERTEX  VX VY VZ
*
*        VX         R 'Scan vertex origin' D=0
*        VY         R 'Scan vertex origin' D=0
*        VZ         R 'Scan vertex origin' D=0
*
*        All Geantinos tracked will start from position VX,VY,VZ.
*
*
* ==>   /SCAN/SFACTORS  FACTX0 FACTL FACTR
*
*        FACTX0     R 'Scale factor for SX0' D=100
*        FACTL      R 'Scale factor for SL' D=1000
*        FACTR      R 'Scale factor for R' D=100
*
*        Set scale factors for SX0,SL and R. The given scale factors must be
*        such that:
*
*          SX0*FACTX0 < 2**15-1 (32767)
*          SL*FACTL   < 2**10-1 (1023)
*          SR*FACTR   < 2**17-1 (131071)
*
*
* ==>   /SCAN/LSCAN  ID [ VOLUME CHOPT ]
*
*        ID         I 'Lego plot identifier' D=2000
*        VOLUME     C 'Volume name' D='XXXX'
*        CHOPT      C 'List of options' D='OPX' R=' ,O,P,I,X,L'
*
*        Generates and plot a table of physics quantities such as
*        the total number of radiation lengths or interaction lengths
*        in function of the SCAN parameters TETA,PHI.
*
*          CHOPT='O' table is generated at Exit  of VOLUME.
*          CHOPT='I' table is generated at Entry of VOLUME.
*          CHOPT='X' radiation lengths
*          CHOPT='L' Interaction lengths
*          CHOPT='P' Plot the table
*
*        If VOLUME='XXXX' Mother volume is used.
*
*
* ==>   /SCAN/HSCAN  IDPHI [ VOLUME CHOPT ]
*
*        IDPHI      I 'Histogram/phi identifier' D=1000
*        VOLUME     C 'Volume name' D='XXXX'
*        CHOPT      C 'List of options' D='OPX' R=' ,O,P,I,X,L'
*
*        Generates and plot an histogram of physics quantities such as
*        the total number of radiation lengths or interaction lengths
*        in function of the SCAN parameter TETA for a given value of PHI.
*
*          CHOPT='O' histogram is generated at Exit  of VOLUME.
*          CHOPT='I' histogram is generated at Entry of VOLUME.
*          CHOPT='X' radiation lengths
*          CHOPT='L' Interaction lengths
*          CHOPT='P' Plot the histogram
*
*        If VOLUME='XXXX' Mother volume is used.
*        The histogram identifier IDPHI is used to also identify which
*        PHI division to plot. IPHI=MOD(IDPHI,1000).
*        If IPHI=0, then all PHI divisions are generated (not plotted)
*        with histogram identifiers IDPHI+PHI division number.
*
*
*
*
*
*
* ***   New commands FILE,REND,MDIR,CDIR,IN,OUT in the RZ menu.
*       =======================================================
*
*
* ==>   RZ/FILE  LUN FNAME [ CHOPT ]
*
*        LUN        I 'Logical unit number'
*        FNAME      C 'File name'
*        CHOPT      C 'Options' D=' ' R=' ,U,N,I,O'
*
*        Open a GEANT/RZ file. Call GRFILE (See below).
*
*         CHOPT=' ' readonly mode
*         CHOPT='U' update mode
*         CHOPT='N' create new file
*         CHOPT='I' Read all structures from existing file
*         CHOPT='O' Write all structures on file
*
*
* ==>   RZ/OUT  OBJECT [ IDVERS ]
*
*        OBJECT     C 'Structure name'
*        IDVERS     I 'Version number' D=1
*
*        Write data structure identified by OBJECT,IDVERS to RZ file.
*        Call GROUT (See below)
*
*          MATE write JMATE structure
*          TMED write JTMED structure
*          VOLU write JVOLUM structure
*          ROTM write JROTM structure
*          SETS write JSET  structure
*          PART write JPART structure
*          SCAN write JSCAN structure
*          *    write all structures
*
* ==>   RZ/IN  OBJECT [ IDVERS ]
*
*        OBJECT     C 'Structure name'
*        IDVERS     I 'Version number' D=1
*
*        Read data structure identified by OBJECT,IDVERS into memory.
*        Call GRIN (See below)
*
*          MATE read JMATE structure
*          TMED read JTMED structure
*          VOLU read JVOLUM structure
*          ROTM read JROTM structure
*          SETS read JSET  structure
*          PART read JPART structure
*          SCAN read JSCAN structure
*          *    read all structures
*
*
*
*
* *** New routines for direct access I/O in the GIOPA package
*     =======================================================
*
* ==>    SUBROUTINE GRFILE(LUN,CHFILE,CHOPT)
*.
*.           Routine to open a GEANT/RZ data base.
*.
*.           LUN logical unit number associated to the file
*.
*.           CHFILE RZ file name
*.
*.           CHOPT is a character string which may be
*.              'N' To create a new file
*.              'U' to open an existing file for update
*.              ' ' to open an existing file for read only
*.              'Q' The initial allocation (default 1000 records)
*.                  is given in IQUEST(10)
*.              'I' Read all data structures from file to memory
*.              'O' Write all data structures from memory to file
*.
*.           Note:
*.             If options 'I' or 'O' all data structures are read or
*.                written from/to file and the file is closed.
*.             See routine GRMDIR to create subdirectories
*.             See routines GROUT,GRIN to write,read objects
*.
*.
*.
*. ==>    SUBROUTINE GROUT(CHOBJ,IDVERS,CHOPT)
*.
*.           Routine to write GEANT object(s) in the RZ file
*.             at the Current Working Directory (See RZCDIR)
*.           Input is taken from the data structures in memory
*.               (VOLU,ROTM,TMED,MATE,SETS,PART,SCAN)
*.
*.           CHOBJ  The type of object to be written:
*.                  MATE write JMATE structure
*.                  TMED write JTMED structure
*.                  VOLU write JVOLUM structure
*.                  ROTM write JROTM structure
*.                  SETS write JSET  structure
*.                  PART write JPART structure
*.                  SCAN write LSCAN structure
*.                  INIT write all initialisation structures
*.
*.           IDVERS is a positive integer which specifies the version
*.               number of the object(s).
*.
*.           CHOPT List of options (none for the time being)
*.
*.        Note that if the cross-sections and energy loss tables
*.           are available in the data structure JMATE, then they are
*.           saved on the data base.
*.
*.
*.        The data structures saved by this routine can be retrieved
*.        with the routine GRIN.
*.
*.        Before calling this routine a RZ data base must have been
*.        created using GRFILE.
*.        The data base must be closed with RZEND.
*.          Ex: if LUN=1 CALL RZEND('LUN1')
*.
*.        The RZ data base can be transported between different
*.        machines in using the ZEBRA RZ utility RZTOFZ.
*.
*.        The interactive version of GEANT provides facilities
*.        to interactively update, create and display objects.
*.
*.          Example.
*.
*.          CALL GRFILE(1,'Geometry.dat','N')
*.          CALL GROUT('VOLU',1,' ')
*.          CALL GROUT('MATE',1,' ')
*.          CALL GROUT('TMED',1,' ')
*.          CALL GROUT('ROTM',1,' ')
*.          CALL GROUT('PART',1,' ')
*.          CALL GROUT('SCAN',1,' ')
*.          CALL GROUT('SETS',1,' ')
*.
*.          The same result can be achieved by:
*.          CALL GRFILE(1,'Geometry.dat','NO')
*.
*.
*.
*. ==>    SUBROUTINE GRIN(CHOBJ,IDVERS,CHOPT)
*.
*.           Routine to read GEANT object(s) fromin the RZ file
*.             at the Current Working Directory (See RZCDIR)
*.           The data structures from disk are read in memory
*.               (VOLU,ROTM,TMED,MATE,SETS,PART,SCAN)
*.
*.           CHOBJ  The type of object to be read:
*.                  MATE read JMATE structure
*.                  TMED read JTMED structure
*.                  VOLU read JVOLUM structure
*.                  ROTM read JROTM structure
*.                  SETS read JSET  structure
*.                  PART read JPART structure
*.                  SCAN read LSCAN structure
*.                  INIT read all initialisation structures
*.
*.           IDVERS is a positive integer which specifies the version
*.               number of the object(s).
*.
*.           CHOPT List of options (none for the time being)
*.
*.
*.        The RZ data base has been created via GRFILE/GROUT
*.
*.
*.          Example.
*.
*.          CALL GRFILE(1,'Geometry.dat',' ')
*.          CALL GRIN ('VOLU',1,' ')
*.          CALL GRIN ('MATE',1,' ')
*.          CALL GRIN ('TMED',1,' ')
*.          CALL GRIN ('ROTM',1,' ')
*.          CALL GRIN ('PART',1,' ')
*.          CALL GRIN ('SCAN',1,' ')
*.          CALL GRIN ('SETS',1,' ')
*.
*.          The same result can be achieved by:
*.          CALL GRFILE(1,'Geometry.dat','I')
*.
*.
*.
* ==>     SUBROUTINE GRMDIR(CHDIR,CHOPT)
*.
*.
*.           Routine to create a subdirectory
*.
*.           CHDIR Subdirectory name
*.
*.           CHOPT is a character string which may be
*.              ' ' To create a subdirectory
*.              'S' To create a subdirectory and set the new
*.                  Current Directory to this directory.
*.
*.
 
