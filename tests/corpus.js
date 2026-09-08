'use strict';

/**
 * Deterministic corpus of Voronoi test cases, shared by the JS golden-master
 * tests and the JS <-> Python differential harness.
 *
 * Randomness is a fixed-seed LCG so both languages generate identical cases.
 */

function makeRng(seed) {
    let s = seed >>> 0;
    return function () {
        s = (s * 1664525 + 1013904223) >>> 0;
        return s / 4294967296;
    };
}

function pickInt(rng, lo, hi) {
    return lo + Math.floor(rng() * (hi - lo + 1));
}

/**
 * Generate a case with `n` unique integer points in [0,w] x [0,h].
 * @param {number} seed
 * @param {number} n
 * @param {number} w
 * @param {number} h
 * @param {Array<Array<number>>} forcedSites - extra fixed sites appended first
 * @returns {{sites: Array<Array<number>>, w: number, h: number}}
 */
function genCase(seed, n, w, h, forcedSites) {
    const rng = makeRng(seed);
    const seen = new Set();
    const sites = [];
    function add(x, y) {
        const k = x + ',' + y;
        if (!seen.has(k)) {
            seen.add(k);
            sites.push([x, y]);
        }
    }
    (forcedSites || []).forEach(function (p) { add(p[0], p[1]); });
    let guard = 0;
    while (sites.length < n && guard++ < 100000) {
        add(pickInt(rng, 0, w), pickInt(rng, 0, h));
    }
    return {sites: sites, w: w, h: h};
}

/**
 * Structured cases that exercise known edge conditions plus seeded random cases.
 * @returns {Array<{sites: Array<Array<number>>, w: number, h: number}>}
 */
function corpus() {
    const cases = [];

    cases.push(genCase(1, 2, 100, 100));
    cases.push(genCase(2, 3, 100, 100));
    cases.push(genCase(3, 4, 100, 100));
    cases.push(genCase(4, 5, 100, 100));

    cases.push(genCase(5, 8, 300, 100));
    cases.push(genCase(6, 12, 80, 240));

    cases.push(genCase(7, 6, 100, 100, [[20, 20], [20, 40], [20, 60]]));
    cases.push(genCase(8, 6, 100, 100, [[20, 20], [40, 20], [60, 20], [80, 20]]));

    cases.push(genCase(9, 5, 100, 100, [[10, 10], [30, 30], [50, 50]]));
    cases.push(genCase(10, 6, 100, 100, [[10, 90], [30, 70], [50, 50]]));

    cases.push(genCase(11, 16, 200, 200));
    cases.push(genCase(12, 32, 200, 200));
    cases.push(genCase(13, 64, 400, 400));

    const sizes = [2, 3, 5, 8, 13, 21, 34, 55, 89];
    const canvases = [[400, 400], [640, 480], [300, 300]];
    let seed = 100;
    sizes.forEach(function (n) {
        canvases.forEach(function (c) {
            cases.push(genCase(seed++, n, c[0], c[1]));
        });
    });

    cases.push(genCase(500, 12, 100, 100,
        [[0, 0], [0, 100], [100, 0], [100, 100], [0, 50], [50, 0], [100, 50], [50, 100]]));

    cases.push({
        sites: [[100, 100], [101, 100], [100, 101], [101, 101], [1000, 1000]],
        w: 2000, h: 2000
    });

    return cases;
}

module.exports = {corpus: corpus, genCase: genCase};
