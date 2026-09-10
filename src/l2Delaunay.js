/**
 * Divide-and-conquer Euclidean Delaunay triangulation (Guibas & Stolfi, 1985).
 * This is the L2 counterpart of divideConquer.js: it splits the x-sorted sites,
 * recurses, then merges the two hulls via the lower common tangent and the
 * "rising bubble" InCircle walk. Only the Delaunay neighbour graph is consumed
 * downstream (see l2Voronoi.js); coordinates are normalised to [0,1] so the
 * floating-point predicates are well conditioned.
 */

import {makeEdge, sym, dest, lnext, oprev, rprev, splice, connect, deleteEdge} from './quadEdge.js';

function sign(x){
    return x > 0 ? 1 : (x < 0 ? -1 : 0);
}

function orient(a, b, c){
    return sign((b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x));
}

function ccw(a, b, c){
    return orient(a, b, c) > 0;
}

function leftOf(x, e){
    return orient(x, e.orig, dest(e)) > 0;
}

function rightOf(x, e){
    return orient(x, dest(e), e.orig) > 0;
}

function valid(e, basel){
    return rightOf(dest(e), basel);
}

/**
 * d lies inside the circumcircle of the CCW triangle a,b,c.
 */
function inCircle(a, b, c, d){
    if(samePoint(a, d) || samePoint(b, d) || samePoint(c, d)){
        return false;
    }

    const sa = a.x * a.x + a.y * a.y;
    const sb = b.x * b.x + b.y * b.y;
    const sc = c.x * c.x + c.y * c.y;
    const sd = d.x * d.x + d.y * d.y;

    const d1 = sc - sd;
    const d2 = c.y - d.y;
    const d3 = c.y * sd - sc * d.y;
    const d4 = c.x - d.x;
    const d5 = c.x * sd - sc * d.x;
    const d6 = c.x * d.y - c.y * d.x;

    const det = a.x * (b.y * d1 - sb * d2 + d3)
        - a.y * (b.x * d1 - sb * d4 + d5)
        + sa * (b.x * d2 - b.y * d4 + d6)
        - b.x * d3 + b.y * d5 - sb * d6;

    return det > 1e-12;
}

function samePoint(a, b){
    return a.x === b.x && a.y === b.y;
}

function delaunay(s){
    let a, b, c, t;

    if(s.length === 2){
        a = makeEdge(s[0], s[1]);
        return {le: a, re: sym(a)};
    }

    if(s.length === 3){
        a = makeEdge(s[0], s[1]);
        b = makeEdge(s[1], s[2]);
        splice(sym(a), b);

        if(ccw(s[0], s[1], s[2])){
            connect(b, a);
            return {le: a, re: sym(b)};
        }
        if(ccw(s[0], s[2], s[1])){
            c = connect(b, a);
            return {le: sym(c), re: c};
        }
        return {le: a, re: sym(b)};
    }

    const half = Math.ceil(s.length / 2);
    const left = delaunay(s.slice(0, half));
    const right = delaunay(s.slice(half));

    let ldo = left.le;
    let ldi = left.re;
    let rdi = right.le;
    let rdo = right.re;

    for(;;){
        if(leftOf(rdi.orig, ldi)){
            ldi = lnext(ldi);
        }
        else if(rightOf(ldi.orig, rdi)){
            rdi = rprev(rdi);
        }
        else{
            break;
        }
    }

    let basel = connect(sym(rdi), ldi);
    if(ldi.orig === ldo.orig){
        ldo = sym(basel);
    }
    if(rdi.orig === rdo.orig){
        rdo = basel;
    }

    for(;;){
        let lcand = sym(basel).onext;
        if(valid(lcand, basel)){
            while(inCircle(dest(basel), basel.orig, dest(lcand), dest(lcand.onext))){
                t = lcand.onext;
                deleteEdge(lcand);
                lcand = t;
            }
        }

        let rcand = oprev(basel);
        if(valid(rcand, basel)){
            while(inCircle(dest(basel), basel.orig, dest(rcand), dest(oprev(rcand)))){
                t = oprev(rcand);
                deleteEdge(rcand);
                rcand = t;
            }
        }

        if(!valid(lcand, basel) && !valid(rcand, basel)){
            break;
        }

        if(!valid(lcand, basel) || (valid(rcand, basel) && inCircle(dest(lcand), lcand.orig, rcand.orig, dest(rcand)))){
            basel = connect(rcand, sym(basel));
        }
        else{
            basel = connect(sym(basel), sym(lcand));
        }
    }

    return {le: ldo, re: rdo};
}

function collectNeighbors(le){
    const adjacency = new Map();
    const visited = new Set();
    const stack = [le];

    function link(a, b){
        if(!adjacency.has(a.i)){
            adjacency.set(a.i, new Set());
        }
        adjacency.get(a.i).add(b.i);
    }

    while(stack.length > 0){
        const e = stack.pop();
        if(visited.has(e)){
            continue;
        }
        visited.add(e);
        link(e.orig, dest(e));
        stack.push(sym(e), lnext(e), oprev(e));
    }

    return adjacency;
}

/**
 * Build the Delaunay neighbour graph of the given sites.
 *
 * @param {Array<Array<number>>} sites
 * @returns {Map<number, Set<number>>} site index -> set of neighbour indices
 */
export function delaunayNeighbors(sites){
    if(sites.length < 2){
        return new Map();
    }

    const maxAbs = Math.max(1, ...sites.map(p => Math.max(Math.abs(p[0]), Math.abs(p[1]))));
    const points = sites.map((p, i) => ({x: p[0] / maxAbs, y: p[1] / maxAbs, i}));

    points.sort((a, b) => a.x === b.x ? a.y - b.y : a.x - b.x);

    const {le} = delaunay(points);
    return collectNeighbors(le);
}
