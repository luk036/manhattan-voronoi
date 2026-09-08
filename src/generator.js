/**
 * Public generator entry points. generateL1Voronoi composes the pipeline
 * stages (preprocess -> sort/init -> divide & conquer -> polygonize); the
 * exported functions are the library facade.
 */

import {recursiveSplit} from './divideConquer.js';
import {polygonizeSite} from './polygonizer.js';
import {cleanData} from './preprocess.js';
import {curryFindBisector, findL1Bisector} from './l1Metric.js';
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

    const findBisector = curryFindBisector(findL1Bisector, width, height);
    const graph = recursiveSplit(sites, findBisector, width, height);

    return graph.map(site => polygonizeSite(site, width, height));
}

/**
 * Generate Voronoi points via a basic, naive algorithm. Takes any distance
 * callback.
 *
 * @param {array} points
 * @param {number} width
 * @param {number} height
 * @param {function} distanceCallback
 * @returns {Array<Array<number>>}
 */
export function generateVoronoiPoints(points, width, height, distanceCallback){

    let colors = points.map(e =>{ return {point:e, color: new Array(3).fill(0).map(d => Math.ceil(Math.random() * 255))}})

    let imageData = new Array(width * height).fill(0).map((point, index) => {
        let coordinate = [index % height , Math.ceil(index / height)];
        let closest = colors.reduce((c,e) => {

            if(Array.isArray(c)){
                return c.every(d => distanceCallback(d.point, coordinate) < distanceCallback(e.point, coordinate) ) ? c : e;
            }
            else if(distanceCallback(c.point, coordinate) === distanceCallback(e.point, coordinate)){
                return [c,e];
            }
            else{
                return distanceCallback(c.point, coordinate) < distanceCallback(e.point, coordinate) ? c : e;
            }

        }, {point:[Infinity,Infinity]});

        return Array.isArray(closest) ? [0,0,0] : closest.color;
    });

    return imageData;
}
