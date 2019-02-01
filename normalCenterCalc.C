/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    normalCenterCalc.C

Group
    grpPostProcessingUtilities

Description
    Calculates the center of the patch, and the normal vector to the surface passing through that center.

Usage
	normalCenterCalc nameOfThePatch

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "argList.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    #include "addRegionOption.H"
    argList::validArgs.append("patchName");
    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createNamedMesh.H"

    const word patchName = args[1];
    //Info<< "arg[1]" << args[1] << endl;

    runTime.setTime(0,0); //(timeDirs[timeI], timeI);

    mesh.readUpdate();

    const label patchI = mesh.boundaryMesh().findPatchID(patchName);
	if (patchI < 0)
	{
		FatalError
			<< "Unable to find patch " << patchName << nl
			<< exit(FatalError);
	}
		
	// number of cells of the chosen patch:
	Info<< "Number of cells of the chosen patch: ";
	Info<< mesh.Cf().boundaryField()[patchI].size() << endl;
	
	vector centerPatch = gSum(mesh.Cf().boundaryField()[patchI])/mesh.Cf().boundaryField()[patchI].size() ;
	//Info<< "Center of the chosen patch: ";
	//Info<< centerPatch << endl;
	
	// back up original center of the patch 
        vector centerPatchOrg = centerPatch ;

    // Computation of the normal vector to the patch. 
    // Let's take the 3 center faces (from the mesh.Cf().boundaryField()[patchI].size() list) that are the most near to centerPatch.
    // then we will compute the cross product using those 3 center faces to 
    // 1) compute a new estimation of the "center" of the patch (if the surface if planar this point will coincide with centerPatch, but if the surface if curved this new point will be a good approximation of the projection of the centerPatch onto the surface)
    // 2) compute the cross product using those 3 points, and we get the normal vector to the patch passing through the new "center"  
        
    double first, second, third;
    double MAX_RANGE = 1e30; 
    first = second = third = MAX_RANGE;   

    // to save the 3 points nearest to the centerPatch
    vector centerP1, centerP2, centerP3; 
    centerP1 = centerP2 = centerP3 = {MAX_RANGE, MAX_RANGE, MAX_RANGE};
    
    double distanceSquare = 0.0; 
    
    for (int i = 0; i < mesh.Cf().boundaryField()[patchI].size() ; i ++)
    {   
		//Info<< "currentFace center " << mesh.Cf().boundaryField()[patchI][i][1] << endl;
		
		// square of the distance
		// (no necessary to make the sqrt, we just need to find the minimum)
		distanceSquare = sqr(centerPatch[0]-mesh.Cf().boundaryField()[patchI][i][0]) 
		               + sqr(centerPatch[1]-mesh.Cf().boundaryField()[patchI][i][1]) 
		               + sqr(centerPatch[2]-mesh.Cf().boundaryField()[patchI][i][2]);
		//Info<< "centerPatch[0] " << centerPatch[0] << endl;  
		//Info<< " " << endl;
		//Info<< "cycle " << i << endl;
		//Info<< "distance is: " << distanceSquare << endl;    
		
        // If current element is smaller than first 
        // then update both first and second 
        if (distanceSquare < first)
        {
			third = second;
            centerP3 = centerP2;            
            
            second = first;
            centerP2 = centerP1;
            
            first = distanceSquare;
            centerP1 = mesh.Cf().boundaryField()[patchI][i];
            //Info<< "center P1 is: " << centerP1 << endl;
            //Info<< "center P2 is: " << centerP2 << endl;
            //Info<< "center P3 is: " << centerP3 << endl;            
            //Info<< "first is: " << first << endl;
            //Info<< "distanceSquare is: " << distanceSquare << endl;
        }
 
        // If distanceSquare is in between first and second 
        // then update second and third.
        else if ( (distanceSquare < second) && (distanceSquare != first) )
        {    
			third = second;
            centerP3 = centerP2; 
            
            second = distanceSquare;
            centerP2 = mesh.Cf().boundaryField()[patchI][i];
            //Info<< "center P2 is: " << centerP2 << endl;
        }   
        

        // If distanceSquare is in between third and fourth 
        // then update third
        else if ( (distanceSquare < third) && (distanceSquare != first) && (distanceSquare != second) )
        {    
            third = distanceSquare;
            centerP3 = mesh.Cf().boundaryField()[patchI][i];
            //Info<< "center P3 is: " << centerP3 << endl;
        }           
         
		            
	}	 
    
    //Info<< "center P1 is: " << centerP1 << endl;
    //Info<< "center P2 is: " << centerP2 << endl;
    //Info<< "center P3 is: " << centerP3 << endl;    
      
    // Compute cross product to get the normal vector to the patch 
    vector normalVector = ((centerP2 - centerP3) ^ (centerP1 - centerP3));
    
    // back up normalVector value before normalization 
    vector normalVectorOrg = normalVector;
    normalVector = normalVector/mag(normalVector);
     
    // new center of the patch computed from averaging coordinates of P1, P2, and P3
    centerPatch = ( centerP1 + centerP2 + centerP3 ) / 3;
    
    // Performing test to know which is the convex and concave side of the patch. 
    // We compute the dot product between the "normalVector" and the vector "centerPatch-centerPatchOrg". 
    // If the angle between "normalVector" and "centerPatch-centerPatchOrg" is greater than pi/2, 
    // it means the "normalVector" is pointing towards the "inside" concave side of the patch, and therefore we must invert the direction
    // This assumption is correct only if we suppose surface where the outside see a CONVEX shape (like car, trains) it would not be correct for a cup.
    vector concaveDir = centerPatch - centerPatchOrg;
    // normalizing
    concaveDir = concaveDir/mag(concaveDir);
    
    // dot product (i.e. inner product)
    double dotProduct  = normalVector & concaveDir;
    
    // angle between normalVector and concaveDir
    double alphaAngle = Foam::acos(double(dotProduct));  // vector have already been normalized
    
    //Info<< "alphaAngle: " << alphaAngle << endl;
    //Info<< "alphaAngle: " << mag(alphaAngle) << endl;
    //Info<< "pi/2: " << degToRad(90) << endl;
    
    
    if ( mag(alphaAngle) > degToRad(90) )   // if angleAngle is > pi/2 then take the opposite direction.
	{
	    Info<< "Vector is pointing towards the inside of the model. Taking the opposite direction. " << endl;
            normalVector = -normalVector;
	    normalVectorOrg = -normalVectorOrg;
	}    
    
    
    Info<< " " << endl;
    Info<< "The normal vector to the patch is: " << normalVectorOrg << endl;
    // normalizing the output
    Info<< "Normalizing the vector: " << normalVector << endl;
    Info<< " " << endl;    
    
    Info<< "Center of the chosen patch: ";
	Info<< centerPatch << endl;

	Info<< endl;
    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
