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
#include "SortableList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::tensor computePCA
(
    const polyMesh& mesh,
    const label& ci
)
{
    const faceList& facePoints = mesh.faces();
    const vectorField& faces = mesh.faceCentres();
    const vectorField& points = mesh.points();

    // PCA
    const labelList& facesI = mesh.cells()[ci];

    symmTensor pca(Zero);

    const vector& ctr = mesh.cellCentres()[ci];

    //label nPts = 0;
    scalar surfArea = 0.0;

    forAll(facesI, fj)
    {
        //nPts++;
        surfArea += mesh.magFaceAreas()[facesI[fj]];

        const vector& cfr = faces[facesI[fj]];
        const labelList& pointsJ = facePoints[facesI[fj]];

        const vector ef(cfr - ctr);
        symmTensor outEf
        (
            ef[0]*ef[0], ef[0]*ef[1], ef[0]*ef[2],
                         ef[1]*ef[1], ef[1]*ef[2],
                                      ef[2]*ef[2]
        );
        outEf *= (scalar(pointsJ.size()) - 1.0)/3.0;

        pca += outEf;

        forAll(pointsJ, pj)
        {
            //nPts++;

            const vector& cfj = points[pointsJ[pj]];
            const vector ej(cfj - ctr);
            symmTensor outEj
            (
                ej[0]*ej[0], ej[0]*ej[1], ej[0]*ej[2],
                             ej[1]*ej[1], ej[1]*ej[2],
                                          ej[2]*ej[2]
            );
            outEj *= 2.0/3.0;

            pca += outEj;
        }
    }

    //pca /= scalar(nPts - 1);
    pca /= surfArea;

    const vector evals = Foam::eigenValues(pca);
    const tensor evecs = Foam::eigenVectors(pca, evals);

    // Compute the system of reference weighted by evals
    return tensor
    (
        evecs.z()*Foam::sqrt(evals.z()),
        evecs.y()*Foam::sqrt(evals.y()),
        evecs.x()*Foam::sqrt(evals.x())
    );
}


Foam::tensor computePDedges
(
    const polyMesh& mesh,
    const label& ci
)
{
    const edgeList& edges = mesh.edges();
    const vectorField& points = mesh.points();

    // PD sorting edges by length
    const labelList& edgesI = mesh.cellEdges()[ci];

    List<vector> edgesV(edgesI.size());
    SortableList<scalar> magEdges(edgesI.size());

    forAll(edgesI, ej)
    {
        edgesV[ej] =
            points[edges[edgesI[ej]][1]]
          - points[edges[edgesI[ej]][0]];
        magEdges[ej] = mag(edgesV[ej]);
    }

    magEdges.sort();

    vector v0 = edgesV[magEdges.indices().last()];
    vector v2 = edgesV[magEdges.indices()[0]];

    SortableList<scalar> minCos(magEdges.size() - 2);
    for (int ej = 1; ej < (magEdges.size() - 1); ej++)
    {
        minCos[ej - 1] = (edgesV[ej] & v0) + (edgesV[ej] & v2);
    }

    minCos.sort();
    const label& idx = magEdges.indices()[minCos.indices()[0] + 1];

    // Gram-schmidt process to find orthogonal base
    vector u1 =
        edgesV[idx] - (edgesV[idx] & v0)/magSqr(v0)*v0;
    vector u2 = v2 - (v2 & v0)/magSqr(v0)*v0 - (v2 & u1)/magSqr(u1)*u1;

    return tensor(v0, u1, u2);
    //return tensor(v0, edgesV[idx], v2);
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

    file<< ctr[0] - pd.xx()/2 << " "
        << ctr[1] - pd.xy()/2 << " "
        << ctr[2] - pd.xz()/2 << nl
        << ctr[0] + pd.xx()/2 << " "
        << ctr[1] + pd.xy()/2 << " "
        << ctr[2] + pd.xz()/2 << nl
        << ctr[0] - pd.yx()/2 << " "
        << ctr[1] - pd.yy()/2 << " "
        << ctr[2] - pd.yz()/2 << nl
        << ctr[0] + pd.yx()/2 << " "
        << ctr[1] + pd.yy()/2 << " "
        << ctr[2] + pd.yz()/2 << nl
        << ctr[0] - pd.zx()/2 << " "
        << ctr[1] - pd.zy()/2 << " "
        << ctr[2] - pd.zz()/2 << nl
        << ctr[0] + pd.zx()/2 << " "
        << ctr[1] + pd.zy()/2 << " "
        << ctr[2] + pd.zz()/2 << nl;

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


void writeVTK
(
    const fvMesh& mesh,
    const tensorField& pd
)
{
    word parRun =
        Pstream::parRun() ?
        "_proc" + Foam::name(Pstream::myProcNo()) : "";

    OFstream file("cellPCA" + parRun + ".vtk");

    file<< "# vtk DataFile Version 2.0" << nl
        << "cell PCA\nASCII\nDATASET UNSTRUCTURED_GRID" << nl;

    label nc = mesh.nCells();

    file<< "POINTS " << 6*nc << " float" << nl;

    forAll(mesh.cells(), ci)
    {
        const vector& ctr = mesh.cellCentres()[ci];
        file<< ctr[0] - pd[ci].xx()/2 << " "
            << ctr[1] - pd[ci].xy()/2 << " "
            << ctr[2] - pd[ci].xz()/2 << nl
            << ctr[0] + pd[ci].xx()/2 << " "
            << ctr[1] + pd[ci].xy()/2 << " "
            << ctr[2] + pd[ci].xz()/2 << nl
            << ctr[0] - pd[ci].yx()/2 << " "
            << ctr[1] - pd[ci].yy()/2 << " "
            << ctr[2] - pd[ci].yz()/2 << nl
            << ctr[0] + pd[ci].yx()/2 << " "
            << ctr[1] + pd[ci].yy()/2 << " "
            << ctr[2] + pd[ci].yz()/2 << nl
            << ctr[0] - pd[ci].zx()/2 << " "
            << ctr[1] - pd[ci].zy()/2 << " "
            << ctr[2] - pd[ci].zz()/2 << nl
            << ctr[0] + pd[ci].zx()/2 << " "
            << ctr[1] + pd[ci].zy()/2 << " "
            << ctr[2] + pd[ci].zz()/2 << nl;
    }

    file<< nl << "CELLS " << 3*nc << " " << 9*nc << nl;
    forAll(mesh.cells(), ci)
    {
        file<< "2 " <<  6*ci      << " " << (6*ci + 1) << nl
            << "2 " << (6*ci + 2) << " " << (6*ci + 3) << nl
            << "2 " << (6*ci + 4) << " " << (6*ci + 5) << nl;
    }

    file<< nl << "CELL_TYPES " << 3*nc << nl;
    forAll(mesh.cells(), ci)
    {
        file<< "3" << nl
            << "3" << nl
            << "3" << nl;
    }

    file<< nl << "CELL_DATA " << 3*nc << nl
        << "FIELD attributes 1" << nl
        << "Eigenvalues 1 " << 3*nc << " float" << nl;
    forAll(mesh.cells(), ci)
    {
        vector evals(mag(pd[ci].x()), mag(pd[ci].y()), mag(pd[ci].z()));

        for (int i = 2; i > -1; i--)
        {
            file<< evals[i] << nl;
        }
    }

}


int main(int argc, char *argv[])
{
    argList::addBoolOption("writeVTK", "write VTK files of principal directions");
    argList::addBoolOption("usePCA", "use PCA to find principal directions");
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

    bool usePCA = args.optionFound("usePCA");

    if (args.optionFound("cell"))
    {
        label ci = args.optionRead<label>("cell");

        Pout<< "Computing principal directions of cell " << ci << endl;

        tensor pd = usePCA ? computePCA(mesh, ci) : computePDedges(mesh, ci);
        Pout<<pd<<nl;

        if (args.optionFound("writeVTK"))
        {
            const vector& ctr = mesh.cellCentres()[ci];
            writeVTK(ci, ctr, pd);
        }
    }
    else
    {
        tensorField pd(mesh.nCells());
        forAll(mesh.cells(), ci)
        {
            Pout<< "Computing principal directions of cell " << ci << endl;

            pd[ci] = usePCA ? computePCA(mesh, ci) : computePDedges(mesh, ci);
        }

        if (args.optionFound("writeVTK"))
        {
            writeVTK(mesh, pd);
        }
    }


    Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
