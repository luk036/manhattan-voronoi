/**
 * Euclidean (L2) Voronoi cells. Each site's cell is the intersection of the
 * canvas rectangle with the perpendicular-bisector half-planes against its
 * Delaunay neighbours (the neighbour graph is produced by the divide-and-conquer
 * triangulation in l2Delaunay.js). Results go through the shared cell.js DTO,
 * so L1 and L2 consumers see the same shape.
 */

import {delaunayNeighbors} from './l2Delaunay.js';
import {toCell} from './cell.js';

function dedupePoints(points){
    const seen = new Set();
    const out = [];
    for(const p of points){
        const key = `${p[0]},${p[1]}`;
        if(!seen.has(key)){
            seen.add(key);
            out.push(p.slice());
        }
    }
    return out;
}

function signedArea(polygon){
    let sum = 0;
    for(let i = 0, j = polygon.length - 1; i < polygon.length; j = i++){
        sum += polygon[j][0] * polygon[i][1] - polygon[i][0] * polygon[j][1];
    }
    return sum / 2;
}

function clipHalfPlane(polygon, site, neighbor){
    const out = [];
    const f = p => (p[0] - site[0]) * (p[0] - site[0]) + (p[1] - site[1]) * (p[1] - site[1])
        - (p[0] - neighbor[0]) * (p[0] - neighbor[0]) - (p[1] - neighbor[1]) * (p[1] - neighbor[1]);

    for(let i = 0; i < polygon.length; i++){
        const A = polygon[i];
        const B = polygon[(i + 1) % polygon.length];
        const fa = f(A);
        const fb = f(B);

        if(fa <= 0){
            out.push(A);
        }
        if((fa < 0 && fb > 0) || (fa > 0 && fb < 0)){
            const ratio = fa / (fa - fb);
            out.push([A[0] + ratio * (B[0] - A[0]), A[1] + ratio * (B[1] - A[1])]);
        }
    }
    return out;
}

function normalizePolygon(polygon){
    const eps = 1e-9;
    const out = [];

    for(const p of polygon){
        const last = out[out.length - 1];
        if(!last || Math.abs(p[0] - last[0]) > eps || Math.abs(p[1] - last[1]) > eps){
            out.push(p);
        }
    }

    if(out.length > 1){
        const first = out[0];
        const last = out[out.length - 1];
        if(Math.abs(first[0] - last[0]) <= eps && Math.abs(first[1] - last[1]) <= eps){
            out.pop();
        }
    }

    if(out.length >= 3 && signedArea(out) < 0){
        out.reverse();
    }
    return out.length >= 3 ? out : [];
}

function clipCell(site, neighbors, width, height){
    let polygon = [[0, 0], [width, 0], [width, height], [0, height]];

    for(const neighbor of neighbors){
        polygon = clipHalfPlane(polygon, site, neighbor);
        if(polygon.length < 3){
            return [];
        }
    }
    return normalizePolygon(polygon);
}

/**
 * Generate an L2 (Euclidean) Voronoi diagram.
 *
 * @param {array} sitePoints
 * @param {number} width
 * @param {number} height
 * @param {boolean} nudgeData
 * @returns {Array<Site>}
 */
export function generateL2Voronoi(sitePoints, width, height, nudgeData = true){
    const sites = nudgeData ? dedupePoints(sitePoints) : sitePoints.map(p => p.slice());

    if(sites.length === 0){
        return [];
    }

    if(sites.length === 1){
        const polygonPoints = [[0, 0], [width, 0], [width, height], [0, height]];
        return [toCell({site: sites[0], polygonPoints, neighbors: [], bisectors: []})];
    }

    const adjacency = delaunayNeighbors(sites);

    return sites.map((site, index) => {
        const neighborPoints = (adjacency.has(index) ? [...adjacency.get(index)] : []).map(j => sites[j]);
        const polygonPoints = clipCell(site, neighborPoints, width, height);
        const cell = toCell({
            site,
            polygonPoints,
            neighbors: neighborPoints.map(p => p.slice()),
            bisectors: []
        });
        if(polygonPoints.length === 0){
            cell.d = '';
        }
        return cell;
    });
}
