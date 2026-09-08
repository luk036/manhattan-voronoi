/**
 * The divide step of the algorithm: recursively split a sorted site list into
 * singletons, building base bisectors for pairs and walking merge lines to
 * combine sibling halves.
 */

import {distance} from './geometry.js';
import {linkBisectorToSites, clearOutOrphans} from './bisector.js';
import {walkMergeLine, determineStartingBisector} from './mergeLine.js';

/**
 * Recursively split and merge sets of points.
 *
 * @param {Array<Site>} splitArray
 * @param {function} findBisector
 * @param {number} width
 * @param {number} height
 * @returns {Array<Site>}
 */
export function recursiveSplit(splitArray, findBisector, width, height){

    if(splitArray.length > 2){
        let splitPoint = (splitArray.length - splitArray.length % 2) / 2

        let L = recursiveSplit(splitArray.slice(0, splitPoint), findBisector, width, height);
        let R = recursiveSplit(splitArray.slice(splitPoint), findBisector, width, height);

        // Order the right-hand sites by distance from the merge start point.
        R.sort((a,b) => distance(L[L.length - 1].site, a.site) - distance(L[L.length - 1].site, b.site));

        let startingInfo = determineStartingBisector(L[L.length - 1], R[0], width, null, findBisector);

        let initialBisector = startingInfo.startingBisector;
        let initialR = startingInfo.nearestNeighbor;
        let initialL = startingInfo.w;

        let upStrokeArray = walkMergeLine(initialR, initialL, initialBisector, [width, height], true, null, [], findBisector);
        let downStrokeArray = walkMergeLine(initialR, initialL, initialBisector, [0, 0], false, null, [], findBisector);

        let mergeArray = [initialBisector, ...upStrokeArray, ...downStrokeArray];

        mergeArray.forEach(bisector => {
            bisector.mergeLine = splitArray.length;
            bisector.sites[0].bisectors = clearOutOrphans(bisector.sites[0], bisector.sites[1]);
            bisector.sites[1].bisectors = clearOutOrphans(bisector.sites[1], bisector.sites[0]);
            linkBisectorToSites(bisector);
        });

        return [...L, ...R];
    }
    else if(splitArray.length === 2){
        let bisector = findBisector(...splitArray);
        linkBisectorToSites(bisector);
        return splitArray;
    }
    else{
        return splitArray;
    }
}
