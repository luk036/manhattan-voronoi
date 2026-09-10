'use strict';

/**
 * Brute-force L2 (Euclidean) nearest-site oracle for validating the
 * divide-and-conquer Euclidean Voronoi output. Test support only.
 */

function l2Distance(a, b){
    const dx = a[0] - b[0];
    const dy = a[1] - b[1];
    return Math.sqrt(dx * dx + dy * dy);
}

function nearestSiteIndex(sites, point){
    let best = 0;
    let bestDist = Infinity;
    for(let i = 0; i < sites.length; i++){
        const d = l2Distance(sites[i], point);
        if(d < bestDist){
            bestDist = d;
            best = i;
        }
    }
    return best;
}

function pointInPolygon(polygon, point){
    let inside = false;
    for(let i = 0, j = polygon.length - 1; i < polygon.length; j = i++){
        const xi = polygon[i][0], yi = polygon[i][1];
        const xj = polygon[j][0], yj = polygon[j][1];
        const crosses = (yi > point[1]) !== (yj > point[1]);
        if(crosses && point[0] < ((xj - xi) * (point[1] - yi)) / (yj - yi) + xi){
            inside = !inside;
        }
    }
    return inside;
}

function pointSegmentDistance(p, a, b){
    const vx = b[0] - a[0];
    const vy = b[1] - a[1];
    const wx = p[0] - a[0];
    const wy = p[1] - a[1];
    const len2 = vx * vx + vy * vy;
    let t = len2 === 0 ? 0 : (wx * vx + wy * vy) / len2;
    t = Math.max(0, Math.min(1, t));
    const dx = p[0] - (a[0] + t * vx);
    const dy = p[1] - (a[1] + t * vy);
    return Math.sqrt(dx * dx + dy * dy);
}

function distanceToBoundary(polygon, point){
    let best = Infinity;
    for(let i = 0, j = polygon.length - 1; i < polygon.length; j = i++){
        best = Math.min(best, pointSegmentDistance(point, polygon[j], polygon[i]));
    }
    return best;
}

/**
 * Every sampled grid point must lie in the cell of the site nearest to it
 * under Euclidean distance. Samples within `boundaryTolerance` of a cell
 * boundary are skipped, because point-in-polygon and rounding legitimately
 * disagree exactly on a shared edge.
 *
 * @param {number} [options.step=5] grid sampling stride
 * @param {number} [options.boundaryTolerance=0.75]
 * @returns {{checked:number, skipped:number, violations:Array<object>}}
 */
function verifyL2Cells(sites, cells, width, height, options = {}){
    const step = options.step || 5;
    const boundaryTolerance = options.boundaryTolerance ?? 0.75;
    const byKey = new Map(cells.map(c => [`${c.site[0]},${c.site[1]}`, c]));
    const violations = [];
    let checked = 0;
    let skipped = 0;

    for(let x = 0; x <= width; x += step){
        for(let y = 0; y <= height; y += step){
            const point = [x, y];
            const expected = sites[nearestSiteIndex(sites, point)];
            const cell = byKey.get(`${expected[0]},${expected[1]}`);
            if(!cell || !cell.polygonPoints || cell.polygonPoints.length < 3){
                violations.push({point, expected, got: null, reason: 'missing cell'});
                continue;
            }
            if(pointInPolygon(cell.polygonPoints, point)){
                checked++;
            }
            else if(distanceToBoundary(cell.polygonPoints, point) <= boundaryTolerance){
                skipped++;
            }
            else{
                violations.push({point, expected, got: cell.site, reason: 'point outside its nearest site cell'});
            }
        }
    }

    return {checked, skipped, violations};
}

module.exports = {l2Distance, nearestSiteIndex, pointInPolygon, distanceToBoundary, verifyL2Cells};
