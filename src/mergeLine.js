/**
 * The merge-line traversal of the Lee & Wong divide-and-conquer merge. Walks a
 * bisector between two merged halves, hopping across the nearest intersecting
 * bisector on either side until the merge line exits the canvas (or an orphaned
 * bisector needs to be unwound).
 */

import {angle, samePoint} from './geometry.js';
import {bisectorIntersection, trimBisector, findHopTo, isBisectorTrapped, getExtremePoint, removeBisector} from './bisector.js';

/**
 * Determine which border is crossed first: "right", "left", or null when both
 * are equidistant.
 *
 * @param {{bisector: Bisector, point: Array}} cropR
 * @param {{bisector: Bisector, point: Array}} cropL
 * @param {Array} currentCropPoint - [x,y]
 * @returns {string|null}
 */
export function determineFirstBorderCross(cropR, cropL, currentCropPoint){
    if(Math.abs(cropR.point[1] - currentCropPoint[1]) === Math.abs(cropL.point[1] - currentCropPoint[1])){
        return null;
    }
    else{
        return Math.abs(cropR.point[1] - currentCropPoint[1]) < Math.abs(cropL.point[1] - currentCropPoint[1]) ? "right" : "left";
    }
}

/**
 * Candidates on one side of the current merge line: the side-site bisectors
 * that intersect the merge line, sorted by the angle of their far site.
 *
 * @param {Site} sideSite
 * @param {Site} otherSite
 * @param {Bisector} currentBisector
 * @param {Array} currentCropPoint - [x,y]
 * @param {Bisector} crossedBorder
 * @param {boolean} goUp
 * @param {Metric} metric
 * @param {boolean} ascending - sort direction (left side sorts descending)
 * @returns {Array<{bisector: Bisector, point: Array}>}
 */
function cropCandidates(sideSite, otherSite, currentBisector, currentCropPoint, crossedBorder, goUp, metric, ascending){
    return sideSite.bisectors
        .map(e => {return {bisector: e, point: bisectorIntersection(currentBisector, e)}})
        .filter(e => {
            let hopTo = findHopTo(e.bisector, sideSite);
            return e.point && (goUp === metric.isUpward(hopTo, sideSite, otherSite, goUp)) && (!samePoint(e.point, currentCropPoint) || e.bisector !== crossedBorder);
        })
        .sort((a, b) => {
            let angleA = angle(sideSite.site, findHopTo(a.bisector, sideSite).site);
            let angleB = angle(sideSite.site, findHopTo(b.bisector, sideSite).site);
            return ascending ? angleA - angleB : angleB - angleA;
        })
        .filter((e, i, candidates) => {
            let hopTo = findHopTo(e.bisector, sideSite);
            let newMergeLine = metric.bisector(otherSite, hopTo);
            trimBisector(newMergeLine, e.bisector, e.point, metric);
            return candidates.every(d => !isBisectorTrapped(findHopTo(d.bisector, sideSite), newMergeLine, metric) || findHopTo(d.bisector, sideSite) === hopTo);
        });
}

function noCrop(goUp){
    return {bisector: null, point: goUp ? [Infinity, Infinity] : [-Infinity, -Infinity]};
}

/**
 * Walk the merge line between two partially-merged diagrams, recording each
 * bisector segment of the merge line in mergeArray.
 *
 * @param {Site} currentR
 * @param {Site} currentL
 * @param {Bisector} currentBisector
 * @param {Array} currentCropPoint - [x,y]
 * @param {boolean} goUp
 * @param {Metric} metric
 * @param {Bisector} crossedBorder
 * @param {Array} mergeArray
 * @returns {Array<Bisector>}
 */
export function walkMergeLine(currentR, currentL, currentBisector, currentCropPoint, goUp, metric, crossedBorder = null, mergeArray = []){

    while(true){

        if(
            !currentBisector.sites.every(e => e === currentR || e === currentL)
        ){
            currentBisector = metric.bisector(currentR, currentL);
            trimBisector(currentBisector, crossedBorder, currentCropPoint, metric);
            mergeArray.push(currentBisector);
        }

        let cropLArray = cropCandidates(currentL, currentR, currentBisector, currentCropPoint, crossedBorder, goUp, metric, false);
        let cropRArray = cropCandidates(currentR, currentL, currentBisector, currentCropPoint, crossedBorder, goUp, metric, true);

        let cropL = cropLArray.length > 0 ? cropLArray[0] : noCrop(goUp);
        let cropR = cropRArray.length > 0 ? cropRArray[0] : noCrop(goUp);

        // If no intersection, the merge line is finished.
        if(
            !cropL.bisector && !cropR.bisector
        ){
            // Check for orphaned bisectors on either side.
            let leftOrphan = checkForOrphans(currentR, currentL, goUp, metric);
            let rightOrphan = checkForOrphans(currentL, currentR, goUp, metric);

            if(leftOrphan){
                removeBisector(leftOrphan);
                let hopTo = findHopTo(leftOrphan, currentL);
                currentR = findCorrectW(currentR, hopTo, metric);
                let newMergeBisector = metric.bisector(hopTo, currentR);
                mergeArray.push(newMergeBisector);
                currentBisector = newMergeBisector;
                currentL = hopTo;
                continue;
            }
            else if(rightOrphan){
                removeBisector(rightOrphan);
                let hopTo = findHopTo(rightOrphan, currentR);
                currentL = findCorrectW(currentL, hopTo, metric);
                let newMergeBisector = metric.bisector(hopTo, currentL);
                mergeArray.push(newMergeBisector);
                currentBisector = newMergeBisector;
                currentR = hopTo;
                continue;
            }

            return mergeArray;
        }

        // Cross the nearest intersecting bisector (or both when equidistant).
        if(determineFirstBorderCross(cropR, cropL, currentCropPoint) === "right"){
            trimBisector(cropR.bisector, currentBisector, cropR.point, metric);
            trimBisector(currentBisector, cropR.bisector, cropR.point, metric);
            currentBisector.intersections.push(cropR.point);
            crossedBorder = cropR.bisector;
            currentR = findHopTo(cropR.bisector, currentR);
            currentCropPoint = cropR.point;
        }
        else if(determineFirstBorderCross(cropR, cropL, currentCropPoint) === "left"){
            trimBisector(cropL.bisector, currentBisector, cropL.point, metric);
            trimBisector(currentBisector, cropL.bisector, cropL.point, metric);
            currentBisector.intersections.push(cropL.point);
            crossedBorder = cropL.bisector;
            currentL = findHopTo(cropL.bisector, currentL);
            currentCropPoint = cropL.point;
        }
        else{
            if(cropR.bisector){
                trimBisector(cropR.bisector, currentBisector, cropR.point, metric);
                trimBisector(currentBisector, cropR.bisector, cropR.point, metric);
                currentBisector.intersections.push(cropR.point);
                crossedBorder = cropR.bisector;
                currentR = findHopTo(cropR.bisector, currentR);
                currentCropPoint = cropR.point;
            }
            if(cropL.bisector){
                trimBisector(cropL.bisector, currentBisector, cropL.point, metric);
                trimBisector(currentBisector, cropL.bisector, cropL.point, metric);
                currentBisector.intersections.push(cropL.point);
                crossedBorder = cropL.bisector;
                currentL = findHopTo(cropL.bisector, currentL);
                currentCropPoint = cropL.point;
            }
        }
    }
}

/**
 * Determine the starting bisector for the merge process.
 *
 * @param {Site} w - starting site
 * @param {Site} nearestNeighbor
 * @param {Metric} metric
 * @param {Array} lastIntersect - [x,y]
 * @returns {{startingBisector: Bisector, w: Site, nearestNeighbor: Site, startingIntersection: Array}}
 */
export function determineStartingBisector(w, nearestNeighbor, metric, lastIntersect = null){

    let z = [metric.width, w.site[1]];

    if(!lastIntersect){
        lastIntersect = w.site;
    }

    let zline = {points: [w.site, z]};

    let intersection = nearestNeighbor.bisectors.map(bisector => {
        return {point: bisectorIntersection(zline, bisector), bisector: bisector}
    }).find(intersection => intersection.point);

    if(intersection && metric.distance(w.site, intersection.point) > metric.distance(nearestNeighbor.site, intersection.point)){
        var startingBisector = metric.bisector(w, nearestNeighbor);
        return {
            startingBisector: startingBisector,
            w: w,
            nearestNeighbor: nearestNeighbor,
            startingIntersection: intersection.point ? intersection.point : w.site
        };
    }
    else if(intersection && metric.distance(w.site, intersection.point) < metric.distance(nearestNeighbor.site, intersection.point) && intersection.point[0] > lastIntersect[0]){
        let nextR = findHopTo(intersection.bisector, nearestNeighbor);
        return determineStartingBisector(w, nextR, metric, intersection.point);
    }
    else{
        w = findCorrectW(w, nearestNeighbor, metric);

        let startingBisector = metric.bisector(w, nearestNeighbor);

        return {
            startingBisector: startingBisector,
            w: w,
            nearestNeighbor: nearestNeighbor,
            startingIntersection: intersection ? intersection.point : w.site
        };
    }
}

/**
 * Ensure the starting point does not produce a trapped bisector.
 *
 * @param {Site} w
 * @param {Site} nearestNeighbor
 * @param {Metric} metric
 * @returns {Site}
 */
export function findCorrectW(w, nearestNeighbor, metric){

    var startingBisector = metric.bisector(w, nearestNeighbor);

    let wTrap = w.bisectors.map(e => {
        let hopTo = findHopTo(e, w);
        return {hopTo: hopTo, isTrapped: isBisectorTrapped(hopTo, startingBisector, metric)}
    })
    .filter(e => e.isTrapped)
    .sort((a,b) => metric.distance(a.hopTo.site, nearestNeighbor.site) - metric.distance(b.hopTo.site, nearestNeighbor.site))[0];

    if(wTrap){
        return findCorrectW(wTrap.hopTo, nearestNeighbor, metric);
    }
    else{
        return w;
    }
}

/**
 * Recursively find an orphaned bisector on the trapped side of a merge.
 *
 * @param {Site} trapper
 * @param {Site} trapped
 * @param {boolean} goUp
 * @param {Metric} metric
 * @returns {Bisector|null}
 */
export function checkForOrphans(trapper, trapped, goUp, metric){

    let orphan = trapped.bisectors.filter(bisector => {
        let hopTo = findHopTo(bisector, trapped);
        return goUp === hopTo.site[1] < trapped.site[1] && isBisectorTrapped(trapper, bisector, metric);
    }).sort((a,b) => {

        let hopToA = findHopTo(a, trapped);
        let hopToB = findHopTo(b, trapped);

        let mergeLineA = metric.bisector(hopToA, trapper);
        let mergeLineB = metric.bisector(hopToB, trapper);

        let extremeA = getExtremePoint(mergeLineA, goUp);
        let extremeB = getExtremePoint(mergeLineB, goUp);

        return goUp ? extremeB - extremeA : extremeA - extremeB;
    })[0];

    return orphan ? orphan : null;
}
