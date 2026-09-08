/**
 * Pure point/segment geometry shared across the Voronoi pipeline.
 */

/**
 * Find L1 distance
 *
 * @param {Array} P1 - [x,y]
 * @param {Array} P2 - [x,y]
 * @returns {number}
 */
export function distance(P1, P2){
    return Math.abs(P1[0] - P2[0]) + Math.abs(P1[1] - P2[1]);
}

/**
 * Determine if two points are the same point
 *
 * @param {Array} P1 - [x,y]
 * @param {Array} P2 - [x,y]
 */
export function samePoint(P1, P2){
    return P1[0] === P2[0] && P1[1] === P2[1];
}

/**
 * Angle of vector P2-P1 in [0, 2*PI)
 *
 * @param {Array} P1 - [x,y]
 * @param {Array} P2 - [x,y]
 */
export function angle(P1, P2){
    let result = Math.atan2(P2[1] - P1[1], P2[0] - P1[0]);

    if(result < 0){
        result = Math.PI + Math.PI + result;
    }

    return result;
}

/**
 * Find intersection of two line segments, if it exists.
 * Returns null when the segments do not intersect.
 *
 * @param {LineSegment} L1 - [[x,y],[x,y]]
 * @param {LineSegment} L2 - [[x,y],[x,y]]
 * @returns {Array|null}
 */
export function segmentIntersection(L1, L2){

    var ua, ub, denom = (L2[1][1] - L2[0][1])*(L1[1][0] - L1[0][0]) - (L2[1][0] - L2[0][0])*(L1[1][1] - L1[0][1]);

    if (denom == 0) {
        return null;
    }
    ua = ((L2[1][0] - L2[0][0])*(L1[0][1] - L2[0][1]) - (L2[1][1] - L2[0][1])*(L1[0][0] - L2[0][0]))/denom;
    ub = ((L1[1][0] - L1[0][0])*(L1[0][1] - L2[0][1]) - (L1[1][1] - L1[0][1])*(L1[0][0] - L2[0][0]))/denom;

    if(
        !(ua >= 0 && ua <= 1 &&
        ub >= 0 && ub <= 1)
    ){
        return null;
    }

    return [
        L1[0][0] + ua*(L1[1][0] - L1[0][0]),
        L1[0][1] + ua*(L1[1][1] - L1[0][1])
    ];

}
