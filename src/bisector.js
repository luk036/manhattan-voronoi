/**
 * Site and bisector graph helpers: object factories plus the predicates and
 * mutation helpers that keep the site <-> bisector graph consistent.
 */

import {distance, samePoint, segmentIntersection} from './geometry.js';

/**
 * Create a site object holding a raw coordinate.
 *
 * @param {Array} point - [x,y]
 * @returns {Site}
 */
export function createSite(point){
    return {site: point, bisectors: []};
}

/**
 * Create an empty bisector between two sites.
 *
 * @param {Array<Site>} sites
 * @param {boolean} up
 * @returns {Bisector}
 */
export function createBisector(sites, up){
    return {sites: sites, up: up, points: [], intersections: [], compound: false};
}

/**
 * Register a bisector with both of its sites.
 *
 * @param {Bisector} bisector
 * @returns {Bisector}
 */
export function linkBisectorToSites(bisector){
    bisector.sites.forEach(site => site.bisectors.push(bisector));
    return bisector;
}

/**
 * Unregister a bisector from both of its sites.
 *
 * @param {Bisector} bisector
 */
export function removeBisector(bisector){
    bisector.sites.forEach(site => {
        site.bisectors = site.bisectors.filter(e => e !== bisector);
    });
}

/**
 * Find the other site across a bisector.
 *
 * @param {Bisector} bisector
 * @param {Site} hopFrom
 * @returns {Site}
 */
export function findHopTo(bisector, hopFrom){
    return bisector.sites.find(e => e !== hopFrom);
}

/**
 * Check if a point lies on the canvas boundary.
 *
 * @param {Array} point - [x,y]
 * @param {number} width
 * @param {number} height
 * @returns {boolean}
 */
export function isPointOnEdge(point, width, height){
    return point[0] === 0 ||
           point[0] === width ||
           point[1] === 0 ||
           point[1] === height;
}

/**
 * Check if two points lie on the same canvas boundary.
 *
 * @param {Array} P1 - [x,y]
 * @param {Array} P2 - [x,y]
 * @param {number} width
 * @param {number} height
 * @returns {boolean}
 */
export function arePointsOnSameEdge(P1, P2, width, height){
    return (P1[0] === P2[0] && P1[0] === 0)     ||
           (P1[0] === P2[0] && P1[0] === width) ||
           (P1[1] === P2[1] && P1[1] === 0)     ||
           (P1[1] === P2[1] && P1[1] === height);
}

/**
 * Determine if a bisector is trapped inside a site's polygon.
 * Trapped means every point of the bisector is at least as close to the
 * trap site as to either of the bisector's own sites.
 *
 * @param {Site} trapPoint
 * @param {Bisector} bisector
 * @returns {boolean}
 */
export function isBisectorTrapped(trapPoint, bisector){
    return bisector.points.every(point => distance(trapPoint.site, point) <= distance(bisector.sites[0].site, point) && distance(trapPoint.site, point) <= distance(bisector.sites[1].site, point));
}

/**
 * Find the highest (goUp) or lowest point of a bisector.
 *
 * @param {Bisector} bisector
 * @param {boolean} goUp
 * @returns {number}
 */
export function getExtremePoint(bisector, goUp){
    return bisector.points.reduce((c,e)=>{
        return goUp ? Math.max(e[1],c) : Math.min(e[1],c);
    }, goUp ? -Infinity : Infinity);
}

/**
 * Trim a bisector at a particular point, discarding the points lying inside
 * the other site's polygon. Mutates target.points.
 *
 * @param {Bisector} target
 * @param {Bisector} intersector
 * @param {Array} intersection - [x,y]
 */
export function trimBisector(target, intersector, intersection){

    let polygonSite = intersector.sites.find(e => target.sites.find(d => d === e) === undefined);

    let newPoints = target.points.filter(e => {
        return distance(e, target.sites[0].site) < distance(e, polygonSite.site) && distance(e, target.sites[1].site) < distance(e, polygonSite.site);
    });

    newPoints.push(intersection);

    target.points = newPoints.sort((a,b) => {
        if(target.up){
            return a[1] - b[1];
        }
        else{
            return a[0] - b[0];
        }
    });

}

/**
 * Find the intersection of two bisectors, if it exists.
 * Returns null when the bisectors do not intersect.
 *
 * @param {Bisector} B1
 * @param {Bisector} B2
 * @returns {Array|null}
 */
export function bisectorIntersection(B1, B2){
    if(B1 === B2){
        return null;
    }
    for(let i = 0; i < B1.points.length - 1; i++){
        for(let j = 0; j < B2.points.length - 1; j++){
            let intersect = segmentIntersection([B1.points[i], B1.points[i+1]], [B2.points[j], B2.points[j+1]]);

            if(intersect){
                return intersect;
            }
        }
    }

    return null;
}

/**
 * Filter out the bisectors of a site that would be trapped by another site.
 *
 * @param {Site} orphanage
 * @param {Site} trapPoint
 * @returns {Array<Bisector>}
 */
export function clearOutOrphans(orphanage, trapPoint){
    return orphanage.bisectors.filter(bisector => !isBisectorTrapped(trapPoint, bisector));
}
