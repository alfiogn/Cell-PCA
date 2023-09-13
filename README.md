# Poly-PCA - Compute principal directions of polyhedra

An OpenFOAM based application to compute principal directions (PDs) of a mesh cell.

The PCA is computed using surface points of a cell, using as correlation matrix.
Let first ${\textbf{p}\_i}$, $i=1,...,N$, be the cell nodes and
let ${\textbf{f}\_i}$, $i=1,...,N\_f$ be the face barycentres.
Let $K$ be a mesh cell and let **o** be its barycentre, then the correlation matrix C reads
$$
\text{C} \simeq \cfrac{1}{\left\vert\partial K\right\vert}
\int_{\partial K} (\mathbf{x} - \textbf{o})(\mathbf{x} - \textbf{o})^\intercal \text{d} \mathbf{x},
$$
where the integral is computed with a quadrature formula using cell nodes and face barycentres.


## Installation
Once OpenFOAM is sourced, go into the `polyPCA` folder and type `wmake`.
Then type `polyPCA -help` to see the application usage.

Tested on OpenFOAM-10.


## Usage
`example-of10`: usage of `polyPCA` on a one-cell mesh.


## TODO
* Input for general polyhedra (not only OpenFOAM mesh cells)
* Pass from surface integration to volume integration for correlation matrix
* Compute PDs on a cell set or a cell zone

