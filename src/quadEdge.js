/**
 * Guibas & Stolfi quad-edge primitives for the divide-and-conquer Euclidean
 * Delaunay triangulation. Adapted from the O(n log n) algorithm in
 * Guibas & Stolfi (1985), "Primitives for the Manipulation of General
 * Subdivisions and the Computation of Voronoi Diagrams"; the port follows the
 * MIT-licensed reference implementation by Philippe Legault
 * (github.com/Bathlamos/delaunay-triangulation).
 *
 * Each undirected edge is four directed records (q0,q1,q2,q3) forming a rot
 * cycle: q0 is orig->dest, q2 is dest->orig, q1/q3 are the dual (Voronoi) side.
 */

function makeEdge(orig, dest){
    const q0 = {onext: null, rot: null, orig};
    const q1 = {onext: null, rot: null, orig: null};
    const q2 = {onext: null, rot: null, orig: dest};
    const q3 = {onext: null, rot: null, orig: null};

    q0.onext = q0; q2.onext = q2;
    q1.onext = q3; q3.onext = q1;

    q0.rot = q1; q1.rot = q2; q2.rot = q3; q3.rot = q0;
    return q0;
}

function sym(e){
    return e.rot.rot;
}

function dest(e){
    return sym(e).orig;
}

function rotSym(e){
    return e.rot.rot.rot;
}

function oprev(e){
    return e.rot.onext.rot;
}

function lnext(e){
    return rotSym(e).onext.rot;
}

function rprev(e){
    return sym(e).onext;
}

function splice(a, b){
    const alpha = a.onext.rot;
    const beta = b.onext.rot;
    const t2 = a.onext;
    const t3 = beta.onext;
    const t4 = alpha.onext;

    a.onext = b.onext;
    b.onext = t2;
    alpha.onext = t3;
    beta.onext = t4;
}

function connect(a, b){
    const q = makeEdge(dest(a), b.orig);
    splice(q, lnext(a));
    splice(sym(q), b);
    return q;
}

function deleteEdge(q){
    splice(q, oprev(q));
    splice(sym(q), oprev(sym(q)));
}

export {makeEdge, sym, dest, rotSym, oprev, lnext, rprev, splice, connect, deleteEdge};
