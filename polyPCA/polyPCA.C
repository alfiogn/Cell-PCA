/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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
    cellPCA

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::tensor computePCA
(
    const polyMesh& mesh,
    const label& ci
)
{
    const cellList& cellFaces = mesh.cells();
    const faceList& facePoints = mesh.faces();
    const vectorField& faces = mesh.faceCentres();
    const vectorField& points = mesh.points();

    // PCA
    const labelList& facesI = cellFaces[ci];

    tensor pd(Zero);

    const vector& ctr = mesh.cellCentres()[ci];

    scalar surfArea = 0.0;

    forAll(facesI, fj)
    {
        surfArea += mesh.magFaceAreas()[facesI[fj]];

        const vector& cfr = faces[facesI[fj]];
        const labelList& pointsJ = facePoints[facesI[fj]];

        const vector ef(cfr - ctr);
        const tensor outEf
        (
            ef[0]*ef[0], ef[0]*ef[1], ef[0]*ef[2],
            ef[1]*ef[0], ef[1]*ef[1], ef[1]*ef[2],
            ef[2]*ef[0], ef[2]*ef[1], ef[2]*ef[2]
        );

        pd += outEf*scalar(pointsJ.size() - 1)/3.0;

        forAll(pointsJ, pj)
        {
            const vector& cfj = points[pointsJ[pj]];
            const vector ej(cfj - ctr);
            const tensor outEj
            (
                ej[0]*ej[0], ej[0]*ej[1], ej[0]*ej[2],
                ej[1]*ej[0], ej[1]*ej[1], ej[1]*ej[2],
                ej[2]*ej[0], ej[2]*ej[1], ej[2]*ej[2]
            );

            pd += outEj*2.0/3.0;
        }
    }

    pd /= surfArea;

    const vector evals = Foam::eigenValues(pd);
    const tensor evecs = Foam::eigenVectors(pd, evals);

    // Compute the system of reference weighted by evals
    pd.z() = evecs.z()*Foam::sqrt(evals.z());
    pd.y() = evecs.y()*Foam::sqrt(evals.y());
    pd.x() = evecs.x()*Foam::sqrt(evals.x());

    return pd;
}


void writeVTK
(
    const label& ci,
    const vector& ctr,
    const tensor& pd
)
{
    OFstream file("cellPCA_" + Foam::name(ci) + ".vtk");

    file<< "# vtk DataFile Version 2.0" << nl
        << "cell PCA\nASCII\nDATASET UNSTRUCTURED_GRID" << nl;

    file<< "POINTS 6 float" << nl;

    file<< ctr[0] - pd.x()[0]/2 << " "
        << ctr[1] - pd.x()[1]/2 << " "
        << ctr[2] - pd.x()[2]/2 << nl
        << ctr[0] + pd.x()[0]/2 << " "
        << ctr[1] + pd.x()[1]/2 << " "
        << ctr[2] + pd.x()[2]/2 << nl
        << ctr[0] - pd.y()[0]/2 << " "
        << ctr[1] - pd.y()[1]/2 << " "
        << ctr[2] - pd.y()[2]/2 << nl
        << ctr[0] + pd.y()[0]/2 << " "
        << ctr[1] + pd.y()[1]/2 << " "
        << ctr[2] + pd.z()[2]/2 << nl
        << ctr[0] - pd.z()[0]/2 << " "
        << ctr[1] - pd.z()[1]/2 << " "
        << ctr[2] - pd.z()[2]/2 << nl
        << ctr[0] + pd.z()[0]/2 << " "
        << ctr[1] + pd.z()[1]/2 << " "
        << ctr[2] + pd.z()[2]/2 << nl;

    file<< nl << "CELLS 3 9" << nl
        << "2 0 1" << nl
        << "2 2 3" << nl
        << "2 4 5" << nl;

    file<< nl << "CELL_TYPES 3" << nl
        << "3" << nl
        << "3" << nl
        << "3" << nl;

    file<< nl << "CELL_DATA 3" << nl
        << "FIELD attributes 1" << nl
        << "Eigenvalues 1 3 float" << nl;

    vector evals(mag(pd.x()), mag(pd.y()), mag(pd.z()));

    for (int i = 2; i > -1; i--)
    {
        file<< evals[i] << nl;
    }
}


int main(int argc, char *argv[])
{
    argList::addBoolOption("writeVTK", "write VTK files of principal directions");
    argList::addOption("cell", "label", "perform for a particular cell");
    argList::addOption("cellSet", "word", "perform for a cell set (TODO)");

    argList::addNote
    (
        "Tool to compute principal directions of mesh cells"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedMesh.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    if (args.optionFound("cell"))
    {
        label ci = args.optionRead<label>("cell");

        Pout<< "Computing principal directions of cell " << ci << endl;

        tensor pd = computePCA(mesh, ci);

        if (args.optionFound("writeVTK"))
        {
            const vector& ctr = mesh.cellCentres()[ci];
            writeVTK(ci, ctr, pd);
        }
    }
    else
    {
        forAll(mesh.cells(), ci)
        {
            Pout<< "Computing principal directions of cell " << ci << endl;

            tensor pd = computePCA(mesh, ci);

            if (args.optionFound("writeVTK"))
            {
                const vector& ctr = mesh.cellCentres()[ci];
                writeVTK(ci, ctr, pd);
            }
        }
    }


    Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
