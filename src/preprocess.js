/**
 * Input preprocessing: nudges degenerate point configurations (duplicate
 * points, points on a square) that the algorithm cannot handle exactly.
 */

/**
 * Nudge points to eliminate square bisectors and duplicate configurations.
 * Mutates and returns the input array; callers pass a copy when the original
 * must be preserved.
 *
 * @param {Array<Array<number>>} data - points in the form [x,y]
 * @returns {Array<Array<number>>}
 */
export function cleanData(data){
    data.forEach((e,i)=> {
        data.forEach((d,j) => {
            if(
                i !== j &&
                Math.abs(d[0] - e[0]) === Math.abs(d[1] - e[1])
            ){
                d[0] = d[0] + 1e-10*d[1];
                d[1] = d[1] + 2e-10*d[0];
            }
        });
    });
    return data;
}
