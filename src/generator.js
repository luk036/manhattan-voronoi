/**
 * Public generator entry point. generateL1Voronoi composes the pipeline
 * stages (preprocess -> sort/init -> divide & conquer -> polygonize). The
 * naive brute-force generator lives in tests/ as an oracle.
 */

import {recursiveSplit} from './divideConquer.js';
import {polygonizeSite} from './polygonizer.js';
import {toCell} from './cell.js';
import {cleanData} from './preprocess.js';
import {createL1Metric} from './l1Metric.js';
import {createSite} from './bisector.js';

function compareByXY(a, b){
    if(a[0] !== b[0]){
        return a[0] - b[0];
    }
    return a[1] - b[1];
}

/**
 * Generate an L1 Voronoi diagram using Lee & Wong's algorithm.
 *
 * @param {array} sitePoints
 * @param {number} width
 * @param {number} height
 * @param {boolean} nudgeData
 * @returns {Array<Site>}
 */
export function generateL1Voronoi(sitePoints, width, height, nudgeData = true){

    let workingPoints = sitePoints;
    if(nudgeData){
        // nudge a copy so the caller's input is never modified
        workingPoints = cleanData(sitePoints.map(p => p.slice()));
    }

    let sites = workingPoints.slice().sort(compareByXY).map(createSite);

    const metric = createL1Metric(width, height);
    const graph = recursiveSplit(sites, metric);

    const cells = graph.map(site => polygonizeSite(site, metric));

    return cells.map(toCell);
}
