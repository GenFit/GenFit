#ifndef G3NODE_H
#define G3NODE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: G3Node.h 220 2007-11-19 16:08:06Z rdm $ */

#include "TNode.h"

class G3Node : public TNode 
{
public:
    G3Node(){}
    G3Node(const char* name, const char* title, const char* shapename,
	    Double_t x = 0, Double_t y = 0, Double_t z = 0, const char* matrixname="",
	    Option_t* option="");

    G3Node(const char* name, const char* title, TShape* shape,
	    Double_t x = 0, Double_t y = 0, Double_t z = 0, TRotMatrix* matrix = 0,
	    Option_t* option="");
    G3Node(const G3Node &node, G3Node* parent);
    
    virtual ~G3Node(){}

    virtual void SetDivision(Int_t ndiv, Int_t axis, Float_t start, Float_t step);
    virtual void ExpandDivisions();
    virtual Int_t   Axis()   const {return fAxis;}
    virtual Int_t   Ndiv()   const {return fNDivision;}
    virtual Float_t Step()   const {return fStep;}
    virtual Float_t StartC() const {return fStartC;}
    virtual void    AddSons(TList* list);
    virtual void    AddSon(G3Node* node);
    
	    
	    
private:
    Int_t   fAxis;         // division axis
    Int_t   fNDivision;    // number of divisions
    Float_t fStep;         // number of steps
    Float_t fStartC;       // start coordinate

    G3Node &operator=(const G3Node &) {return *this;}

    ClassDef(G3Node,1) // G3 Node for the G3 GUI 
};

#endif








