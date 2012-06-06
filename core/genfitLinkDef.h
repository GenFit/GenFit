/* Copyright 2008-2010, Technische Universitaet Muenchen,
   Authors: Christian Hoeppner & Sebastian Neubert

   This file is part of GENFIT.

   GENFIT is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published
   by the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   GENFIT is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with GENFIT.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifdef __CINT__


#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;


#pragma link C++ class GFAbsTrackRep+;
#pragma link C++ class GFDafHit+;
#pragma link C++ class GFPlanarHitPolicy+;
#pragma link C++ class GFPseudoSpacepointWireHitPolicy+;
#pragma link C++ class GFSpacepointHitPolicy+;
#pragma link C++ class GFWirepointHitPolicy+;
#pragma link C++ class GFWireHitPolicy+;
#pragma link C++ class GFRecoHitIfc<GFPlanarHitPolicy>+;
#pragma link C++ class GFRecoHitIfc<GFPseudoSpacepointWireHitPolicy>+;
#pragma link C++ class GFRecoHitIfc<GFSpacepointHitPolicy>+;
#pragma link C++ class GFRecoHitIfc<GFWirepointHitPolicy>+;
#pragma link C++ class GFRecoHitIfc<GFWireHitPolicy>+;
#pragma link C++ class GFAbsRecoHit+;
#pragma link C++ class GFTrackCand+;
#pragma link C++ class GFTrack+;
#pragma link C++ class GFDetPlane+;
#pragma link C++ class GFBookkeeping-;
#pragma link C++ class GFAbsFinitePlane+;
#pragma link C++ class GFRectFinitePlane+;
#pragma link C++ class GFFieldManager+;

#endif
