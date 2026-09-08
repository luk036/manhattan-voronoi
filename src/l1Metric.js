/**
 * The L1 metric: bisector construction and the upward-test used while walking
 * merge lines. The rest of the pipeline consumes these through a curried
 * `findBisector(siteA, siteB)` function, so this module is the seam where a
 * different metric (e.g. L-infinity) would slot in.
 */

import {samePoint} from './geometry.js';
import {createBisector} from './bisector.js';

/**
 * Curry a bisector factory with fixed canvas dimensions.
 *
 * @param {function} callback
 * @param {number} width
 * @param {number} height
 * @return {function}
 */
export function curryFindBisector(callback, width, height){
    return function(P1, P2){
        return callback(P1, P2, width, height);
    }
}

/**
 * Generate an L1 bisector between two sites.
 *
 * @param {Site} P1
 * @param {Site} P2
 * @param {number} width
 * @param {number} height
 * @returns {Bisector}
 */
export function findL1Bisector(P1, P2, width, height){

    let xDistance = P1.site[0] - P2.site[0];
    let yDistance = P1.site[1] - P2.site[1];

    let midpoint = [
        (P1.site[0] + P2.site[0]) / 2,
        (P1.site[1] + P2.site[1]) / 2
    ];

    let vertexes = [];
    let up = null;

    if(samePoint(P1.site,P2.site)){
        throw new Error(`Duplicate point: Points ${JSON.stringify(P1)} and ${JSON.stringify(P2)} are duplicates. please remove one`);
    }

    if(Math.abs(xDistance) === 0){
        vertexes = [
            [0, midpoint[1]],
            [width, midpoint[1]]
        ];

        return {sites:[P1, P2], up:false, points:vertexes, intersections:[], compound:false};
    }

    if(Math.abs(yDistance) === 0){
        vertexes = [
            [midpoint[0], 0],
            [midpoint[0], height]
        ];

        return {sites:[P1, P2], up:true, points:vertexes, intersections:[], compound:false};
    }

    let slope = yDistance/xDistance > 0 ? -1 : 1;
    let intercept = midpoint[1] - midpoint[0] * slope;

    if(Math.abs(xDistance) > Math.abs(yDistance)){
        vertexes = [
            [(P1.site[1] - intercept) / slope, P1.site[1]],
            [(P2.site[1] - intercept) / slope, P2.site[1]]
        ];

        up = true;
    }
    else if(Math.abs(xDistance) < Math.abs(yDistance)){
        vertexes = [
            [P1.site[0] , (P1.site[0] * slope) + intercept ],
            [P2.site[0] , (P2.site[0] * slope) + intercept ]
        ];

        up = false;
    }
    else { // |dx| === |dy| (square bisector)
        if (slope === 1){
            vertexes = [
                [P1.site[1] - intercept, P1.site[1]],
                [P2.site[1] - intercept, P2.site[1]]
            ];

            up = true;
        }
        else { // slope === -1
            vertexes = [
                [P1.site[0] , -P1.site[0] + intercept ],
                [P2.site[0] , -P2.site[0] + intercept ]
            ];

            up = false;
        }
    }

    let bisector = createBisector([P1, P2], up);

    if(up){
        const sortedVerts = vertexes.sort((a,b) => a[1] - b[1]);

        bisector.points = [
            [sortedVerts[0][0], 0],
            ...sortedVerts,
            [sortedVerts[1][0], height]
        ].sort((a,b) => a[1] - b[1]);

    }
    else{
        const sortedVerts = vertexes.sort((a,b) => a[0] - b[0]);

        bisector.points = [
            [0,sortedVerts[0][1]],
            ...sortedVerts,
            [width,sortedVerts[1][1]]
        ].sort((a,b) => a[0] - b[0]);
    }

    return bisector;
}

/**
 * Check whether the bisector between a hop pair travels upward relative to the
 * merge walk direction. Metric-specific geometry on raw site coordinates.
 *
 * @param {Site} hopTo
 * @param {Site} hopFrom
 * @param {Site} site
 * @param {boolean} goUp
 * @returns {boolean}
 */
export function isNewBisectorUpward(hopTo, hopFrom, site, goUp){

    if(hopTo.site[0] - site.site[0] === 0){
        return site.site[1] > hopTo.site[1];
    }

    let slope = (hopTo.site[1] - site.site[1])/(hopTo.site[0] - site.site[0]);
    let intercept = hopTo.site[1] - (slope * hopTo.site[0]);

    if(Math.abs(slope) === Infinity){
        return site.site[1] > hopTo.site[1];
    }

    let isAboveLine = hopFrom.site[1] > (slope * hopFrom.site[0]) + intercept;

    return isAboveLine;
}
