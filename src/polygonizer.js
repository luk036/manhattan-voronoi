/**
 * Post-processing stage: for each site, chain its bisectors into an ordered
 * polygon, then enrich the site with the polygon points, its SVG path and its
 * neighbor list. The SVG path helper is the presentation adapter for the
 * computational core.
 */

import {distance, angle, samePoint} from './geometry.js';
import {isPointOnEdge, arePointsOnSameEdge, bisectorIntersection, findHopTo} from './bisector.js';

/**
 * Chain the bisectors of a site into a closed polygon point list.
 *
 * @param {Site} site
 * @param {number} width
 * @param {number} height
 * @returns {Array<Array<number>>}
 */
export function chainBisectorPoints(site, width, height){
    if(site.bisectors.length === 0){
        // a site with no bisectors owns the whole canvas
        return [
            [0, 0],
            [width, 0],
            [width, height],
            [0, height]
        ];
    }

    return site.bisectors.reduce((total, bisector, index, bisectors)=>{

        if(index === 0){

            // start from a bisector that touches the canvas edge, if any
            let startBisector = bisectors.find(e => {
                return e.points.some(e => isPointOnEdge(e, width, height));
            }) || bisector;

            let startingPoints = startBisector.points;

            if(isPointOnEdge(startingPoints[startingPoints.length - 1], width, height)){
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

                let eDistance = distance(last, e.points[0]) < distance(last, e.points[e.points.length - 1]) ? distance(last, e.points[0]) : distance(last, e.points[e.points.length - 1]);
                let cDistance = distance(last, c.points[0]) < distance(last, c.points[c.points.length - 1]) ? distance(last, c.points[0]) : distance(last, c.points[c.points.length - 1]);

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
 * @param {number} width
 * @param {number} height
 * @returns {Array<Array<number>>}
 */
function appendOpenEdgeCorners(site, polygonPoints, width, height){
    const corners = [
        [0, 0],
        [width, 0],
        [width, height],
        [0, height]
    ];

    if(
        isPointOnEdge(polygonPoints[0], width, height) &&
        isPointOnEdge(polygonPoints[polygonPoints.length - 1], width, height) &&
        !arePointsOnSameEdge(polygonPoints[0], polygonPoints[polygonPoints.length - 1], width, height)
    ){
        const filteredCorners = corners.filter(e => {
            return site.bisectors.every(d => !bisectorIntersection({points: [e, site.site]}, d));
        });
        return [...polygonPoints, ...filteredCorners];
    }
    return polygonPoints;
}

/**
 * Render an SVG path string for a polygon point list.
 *
 * @param {Array<Array<number>>} polygonPoints
 * @returns {string}
 */
export function toSVGPath(polygonPoints){
    return `M ${ polygonPoints.map(e => e.join(" ")).join(" L")} Z`;
}

/**
 * Polygonize and enrich a site in place: adds polygonPoints, d, neighbors.
 *
 * @param {Site} site
 * @param {number} width
 * @param {number} height
 * @returns {Site}
 */
export function polygonizeSite(site, width, height){

    let polygonPoints = chainBisectorPoints(site, width, height);

    polygonPoints = appendOpenEdgeCorners(site, polygonPoints, width, height);

    site.polygonPoints = polygonPoints.sort((a,b)=> angle(site.site, a) - angle(site.site, b));

    site.d = toSVGPath(site.polygonPoints);

    site.neighbors = site.bisectors.map(e => findHopTo(e, site).site);

    return site;
}
