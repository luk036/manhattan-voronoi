'use strict';

/**
 * Brute-force Voronoi oracle. Assigns every pixel of a width x height grid to
 * its nearest site under the supplied distance callback. This is
 * O(width * height * sites) and is kept only as a test reference; the shipped
 * library uses the Lee & Wong divide-and-conquer generator instead.
 */

/**
 * Generate a per-pixel color grid by brute force.
 *
 * @param {Array<Array<number>>} points - sites in the form [x,y]
 * @param {number} width
 * @param {number} height
 * @param {function} distanceCallback - (pointA, pointB) => number
 * @returns {Array<Array<number>>} one [r,g,b] per pixel
 */
function generateVoronoiPoints(points, width, height, distanceCallback){

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

module.exports = {generateVoronoiPoints};
