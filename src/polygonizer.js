/**
 * Post-processing stage: for each site, chain its bisectors into an ordered
 * polygon, then enrich the site with the polygon points and its neighbor list.
 */

import {angle, samePoint} from './geometry.js';
import {isPointOnEdge, arePointsOnSameEdge, bisectorIntersection, findHopTo} from './bisector.js';

/**
 * Chain the bisectors of a site into a closed polygon point list.
 *
 * @param {Site} site
 * @param {Metric} metric
 * @returns {Array<Array<number>>}
 */
export function chainBisectorPoints(site, metric){
    if(site.bisectors.length === 0){
        // a site with no bisectors owns the whole canvas
        return [
            [0, 0],
            [metric.width, 0],
            [metric.width, metric.height],
            [0, metric.height]
        ];
    }

    return site.bisectors.reduce((total, bisector, index, bisectors)=>{

        if(index === 0){

            // start from a bisector that touches the canvas edge, if any
            let startBisector = bisectors.find(e => {
                return e.points.some(e => isPointOnEdge(e, metric.width, metric.height));
            }) || bisector;

            let startingPoints = startBisector.points;

            if(isPointOnEdge(startingPoints[startingPoints.length - 1], metric.width, metric.height)){
                startingPoints = startingPoints.reverse();
            }

            return {
                points: startingPoints,
                used: [startBisector]
            };
        }
        else{
            let last = total.points[total.points.length - 1];

            let nextBisector = bisectors.filter(e => total.used.every(d => e !== d)).reduce((c,e) => {

                let eDistance = metric.distance(last, e.points[0]) < metric.distance(last, e.points[e.points.length - 1]) ? metric.distance(last, e.points[0]) : metric.distance(last, e.points[e.points.length - 1]);
                let cDistance = metric.distance(last, c.points[0]) < metric.distance(last, c.points[c.points.length - 1]) ? metric.distance(last, c.points[0]) : metric.distance(last, c.points[c.points.length - 1]);

                return eDistance < cDistance ? e : c;
            },{points:[[Infinity, Infinity]]});

            let nextPoints = nextBisector.points;

            if(samePoint(nextPoints[nextPoints.length - 1], last)){
                nextPoints = nextPoints.reverse();
            }

            return {
                points: [...total.points, ...nextPoints],
                used: [...total.used, nextBisector]
            };
        }
    }, {}).points;
}

/**
 * Fill in canvas corners when the polygon opens onto two different edges.
 *
 * @param {Site} site
 * @param {Array<Array<number>>} polygonPoints
 * @param {Metric} metric
 * @returns {Array<Array<number>>}
 */
function appendOpenEdgeCorners(site, polygonPoints, metric){
    const corners = [
        [0, 0],
        [metric.width, 0],
        [metric.width, metric.height],
        [0, metric.height]
    ];

    if(
        isPointOnEdge(polygonPoints[0], metric.width, metric.height) &&
        isPointOnEdge(polygonPoints[polygonPoints.length - 1], metric.width, metric.height) &&
        !arePointsOnSameEdge(polygonPoints[0], polygonPoints[polygonPoints.length - 1], metric.width, metric.height)
    ){
        const filteredCorners = corners.filter(e => {
            return site.bisectors.every(d => !bisectorIntersection({points: [e, site.site]}, d));
        });
        return [...polygonPoints, ...filteredCorners];
    }
    return polygonPoints;
}

/**
 * Polygonize and enrich a site in place: adds polygonPoints and neighbors.
 *
 * @param {Site} site
 * @param {Metric} metric
 * @returns {Site}
 */
export function polygonizeSite(site, metric){

    let polygonPoints = chainBisectorPoints(site, metric);

    polygonPoints = appendOpenEdgeCorners(site, polygonPoints, metric);

    site.polygonPoints = polygonPoints.sort((a,b)=> angle(site.site, a) - angle(site.site, b));

    site.neighbors = site.bisectors.map(e => findHopTo(e, site).site);

    return site;
}
