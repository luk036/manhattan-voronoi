'use strict';

const test = require('node:test');
const assert = require('node:assert');
const voronoi = require('../dist/voronoi.js');
const {verifyL2Cells} = require('./l2Oracle.js');

function cellArea(cell){
    const p = cell.polygonPoints;
    let sum = 0;
    for(let i = 0, j = p.length - 1; i < p.length; j = i++){
        sum += p[j][0] * p[i][1] - p[i][0] * p[j][1];
    }
    return Math.abs(sum) / 2;
}

const corpora = [
    {name: 'two sites', sites: [[20, 50], [80, 50]], width: 100, height: 100},
    {name: 'three sites', sites: [[10, 10], [90, 90], [50, 20]], width: 100, height: 100},
    {name: 'five sites', sites: [[23, 37], [50, 71], [1, 2], [80, 20], [60, 90]], width: 100, height: 100},
    {name: 'nine-site grid', sites: [[10, 10], [10, 50], [10, 90], [50, 10], [50, 50], [50, 90], [90, 10], [90, 50], [90, 90]], width: 100, height: 100},
    {name: 'non-square canvas', sites: [[30, 20], [200, 40], [120, 150], [260, 90], [60, 120]], width: 300, height: 200}
];

corpora.forEach(function (corpus) {
    test('L2 matches the brute-force oracle: ' + corpus.name, function () {
        const cells = voronoi.generateL2Voronoi(corpus.sites, corpus.width, corpus.height, true);
        assert.strictEqual(cells.length, corpus.sites.length);
        const report = verifyL2Cells(corpus.sites, cells, corpus.width, corpus.height, {step: 5, boundaryTolerance: 1});
        assert.deepStrictEqual(report.violations, []);
        assert.ok(report.checked > 0, 'some samples were verified');
    });
});

test('L2 single site owns the whole canvas', function () {
    const result = voronoi.generateL2Voronoi([[10, 10]], 100, 100, false);
    assert.strictEqual(result.length, 1);
    assert.strictEqual(result[0].polygonPoints.length, 4);
    assert.deepStrictEqual(result[0].neighbors, []);
});

test('L2 does not mutate the caller input array', function () {
    const input = [[30, 30], [10, 10], [20, 20]];
    const snapshot = JSON.stringify(input);
    voronoi.generateL2Voronoi(input, 100, 100, true);
    assert.strictEqual(JSON.stringify(input), snapshot);
});

test('L2 cells expose the documented public shape', function () {
    const result = voronoi.generateL2Voronoi([[10, 10], [90, 90], [50, 20]], 100, 100, false);
    result.forEach(function (s) {
        assert.ok(Array.isArray(s.site) && s.site.length === 2, 'site is [x,y]');
        assert.ok(Array.isArray(s.bisectors), 'bisectors is an array');
        assert.ok(Array.isArray(s.polygonPoints), 'polygonPoints is an array');
        assert.strictEqual(typeof s.d, 'string', 'd is an SVG path string');
        assert.ok(/^M /.test(s.d), 'd starts with an SVG moveto');
        assert.ok(Array.isArray(s.neighbors), 'neighbors is an array');
    });
});

test('L2 cells partition the canvas area', function () {
    const sites = [[23, 37], [50, 71], [1, 2], [80, 20], [60, 90]];
    const cells = voronoi.generateL2Voronoi(sites, 100, 100, true);
    const total = cells.reduce(function (acc, cell) { return acc + cellArea(cell); }, 0);
    assert.ok(Math.abs(total - 100 * 100) < 1e-6, 'total area was ' + total);
});

test('L2 neighbors are reciprocal', function () {
    const sites = [[23, 37], [50, 71], [1, 2], [80, 20], [60, 90]];
    const cells = voronoi.generateL2Voronoi(sites, 100, 100, true);
    const byKey = new Map(cells.map(function (c) { return [c.site.join(','), c]; }));
    cells.forEach(function (c) {
        c.neighbors.forEach(function (n) {
            const other = byKey.get(n.join(','));
            assert.ok(other, 'neighbor cell exists');
            assert.ok(other.neighbors.some(function (x) { return x.join(',') === c.site.join(','); }), 'reciprocal');
        });
    });
});

test('L2 handles cocircular sites', function () {
    const sites = [[30, 50], [50, 30], [70, 50], [50, 70]];
    const cells = voronoi.generateL2Voronoi(sites, 100, 100, true);
    const report = verifyL2Cells(sites, cells, 100, 100, {step: 5, boundaryTolerance: 1});
    assert.deepStrictEqual(report.violations, []);
});

test('L2 handles collinear sites', function () {
    const sites = [[10, 50], [30, 50], [50, 50], [70, 50], [90, 50]];
    const cells = voronoi.generateL2Voronoi(sites, 100, 100, true);
    const report = verifyL2Cells(sites, cells, 100, 100, {step: 5, boundaryTolerance: 1});
    assert.deepStrictEqual(report.violations, []);
});

test('L2 collapses duplicate sites', function () {
    const cells = voronoi.generateL2Voronoi([[10, 10], [10, 10], [90, 90]], 100, 100, true);
    assert.strictEqual(cells.length, 2);
});
