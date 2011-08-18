/* Copyright 2008-2009, Technische Universitaet Muenchen,
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
// Description:
//      Function to extrapolate two tracks to common vertex
/**
 *  @author Christian H&ouml;ppner (Technische Universit&auml;t M&uuml;nchen, original author)
 *  @author Sebastian Neubert  (Technische Universit&auml;t M&uuml;nchen, original author)
 *  @author Sverre Doerheim  (Technische Universit&auml;t M&uuml;nchen, author)
 *
 */
#ifndef EXTRAPOLATETOPCA_H
#define EXTRAPOLATETOPCA_H

class GFTrack;
class GFAbsTrackRep; 
class TVector3;
//!Calculates poca between two tracks, changes the state of the track!
TVector3 trackProximity(GFTrack* trk1, GFTrack* trk2);  
//!Calculates poca between two track reps, changes the state of the trackrep!
TVector3 trackProximity(GFAbsTrackRep* rep1, GFAbsTrackRep* rep2);  

#endif


