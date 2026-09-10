'use strict';

const test = require('node:test');
const assert = require('node:assert');
const voronoi = require('../dist/voronoi.js');

test('single site owns the whole canvas', function () {
    const result = voronoi.generateL1Voronoi([[10, 10]], 100, 100, false);
    assert.strictEqual(result.length, 1);
    assert.strictEqual(result[0].polygonPoints.length, 4);
    assert.deepStrictEqual(result[0].neighbors, []);
});

test('does not mutate the caller input array', function () {
    const input = [[30, 30], [10, 10], [20, 20]];
    const snapshot = JSON.stringify(input);
    voronoi.generateL1Voronoi(input, 100, 100, false);
    assert.strictEqual(JSON.stringify(input), snapshot);
});

test('nudgeData nudges diagonal alignments when enabled', function () {
    const aligned = [[10, 10], [30, 30]];
    const nudged = voronoi.generateL1Voronoi(aligned, 100, 100, true);
    const coords = nudged.map(function (s) { return s.site; });
    assert.deepStrictEqual(coords, [[10, 10], [30.000000003, 30.000000006]]);
});

test('nudgeData disabled leaves coordinates untouched', function () {
    const aligned = [[10, 10], [30, 30]];
    const raw = voronoi.generateL1Voronoi(aligned, 100, 100, false);
    assert.deepStrictEqual(raw.map(function (s) { return s.site; }), aligned);
});

test('two sites produce reciprocal neighbor relationship', function () {
    const result = voronoi.generateL1Voronoi([[10, 10], [90, 90]], 100, 100, false);
    const byKey = {};
    result.forEach(function (s) { byKey[s.site.join(',')] = s; });
    assert.deepStrictEqual(byKey['10,10'].neighbors, [[90, 90]]);
    assert.deepStrictEqual(byKey['90,90'].neighbors, [[10, 10]]);
});

test('returned cell exposes the documented public shape', function () {
    const result = voronoi.generateL1Voronoi(
        [[10, 10], [90, 90], [50, 20]], 100, 100, false);
    result.forEach(function (s) {
        assert.ok(Array.isArray(s.site) && s.site.length === 2, 'site is [x,y]');
        assert.ok(Array.isArray(s.bisectors), 'bisectors is an array');
        assert.ok(Array.isArray(s.polygonPoints), 'polygonPoints is an array');
        assert.strictEqual(typeof s.d, 'string', 'd is an SVG path string');
        assert.ok(/^M /.test(s.d), 'd starts with an SVG moveto');
        assert.ok(Array.isArray(s.neighbors), 'neighbors is an array');
    });
});

test('reported neighbors never include the site itself', function () {
    const result = voronoi.generateL1Voronoi(
        [[23, 37], [50, 71], [1, 2], [80, 20], [60, 90]], 100, 100, false);
    result.forEach(function (s) {
        s.neighbors.forEach(function (n) {
            assert.notDeepStrictEqual(n, s.site);
        });
    });
});
