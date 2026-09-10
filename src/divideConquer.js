/**
 * The divide step of the algorithm: recursively split a sorted site list into
 * singletons, building base bisectors for pairs and walking merge lines to
 * combine sibling halves.
 */

import {linkBisectorToSites, clearOutOrphans} from './bisector.js';
import {walkMergeLine, determineStartingBisector} from './mergeLine.js';

/**
 * Recursively split and merge sets of points.
 *
 * @param {Array<Site>} splitArray
 * @param {Metric} metric
 * @returns {Array<Site>}
 */
export function recursiveSplit(splitArray, metric){

    if(splitArray.length > 2){
        let splitPoint = (splitArray.length - splitArray.length % 2) / 2

        let L = recursiveSplit(splitArray.slice(0, splitPoint), metric);
        let R = recursiveSplit(splitArray.slice(splitPoint), metric);

        // Order the right-hand sites by distance from the merge start point.
        R.sort((a,b) => metric.distance(L[L.length - 1].site, a.site) - metric.distance(L[L.length - 1].site, b.site));

        let startingInfo = determineStartingBisector(L[L.length - 1], R[0], metric);

        let initialBisector = startingInfo.startingBisector;
        let initialR = startingInfo.nearestNeighbor;
        let initialL = startingInfo.w;

        let upStrokeArray = walkMergeLine(initialR, initialL, initialBisector, [metric.width, metric.height], true, metric);
        let downStrokeArray = walkMergeLine(initialR, initialL, initialBisector, [0, 0], false, metric);

        let mergeArray = [initialBisector, ...upStrokeArray, ...downStrokeArray];

        mergeArray.forEach(bisector => {
            bisector.mergeLine = splitArray.length;
            bisector.sites[0].bisectors = clearOutOrphans(bisector.sites[0], bisector.sites[1], metric);
            bisector.sites[1].bisectors = clearOutOrphans(bisector.sites[1], bisector.sites[0], metric);
            linkBisectorToSites(bisector);
        });

        return [...L, ...R];
    }
    else if(splitArray.length === 2){
        let bisector = metric.bisector(...splitArray);
        linkBisectorToSites(bisector);
        return splitArray;
    }
    else{
        return splitArray;
    }
}
