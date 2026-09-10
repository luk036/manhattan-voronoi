/**
 * Presentation boundary: map the internal site/bisector graph onto plain result
 * objects. Consumers and adapters (SVG, JSON, Canvas) see these DTOs, so the
 * live, mutable graph stays an implementation detail.
 */

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
 * Map an internal bisector onto a plain result object.
 *
 * @param {Bisector} bisector
 * @returns {object}
 */
export function toBisectorDTO(bisector){
    return {
        sites: bisector.sites.map(site => ({site: site.site.slice()})),
        up: bisector.up,
        points: bisector.points.map(point => point.slice()),
        intersections: bisector.intersections.map(point => point.slice()),
        compound: bisector.compound,
        mergeLine: bisector.mergeLine === undefined ? null : bisector.mergeLine
    };
}

/**
 * Map an internal site onto its public cell, applying the SVG adapter.
 *
 * @param {Site} site
 * @returns {object}
 */
export function toCell(site){
    return {
        site: site.site.slice(),
        polygonPoints: site.polygonPoints.map(point => point.slice()),
        d: toSVGPath(site.polygonPoints),
        neighbors: site.neighbors.map(neighbor => neighbor.slice()),
        bisectors: site.bisectors.map(toBisectorDTO)
    };
}
