'use strict';

/**
 * Golden-master characterization tool. Records exact JS library output for the
 * corpus into tests/fixtures/js-golden.json, and (in compare mode) fails if the
 * rebuilt dist output deviates from the recorded snapshot.
 */

const fs = require('fs');
const path = require('path');
const assert = require('assert');
const voronoi = require('../dist/voronoi.js');
const corpusModule = require('./corpus.js');

const FIXTURE = path.join(__dirname, 'fixtures', 'js-golden.json');

function runCase(c) {
    const sites = c.sites.map(function (p) { return p.slice(); });
    return voronoi.generateL1Voronoi(sites, c.w, c.h, false);
}

function serializeResult(sites) {
    return sites.map(function (site) {
        return {
            site: site.site,
            polygonPoints: site.polygonPoints,
            d: site.d,
            neighbors: site.neighbors,
            bisectors: site.bisectors.map(function (b) {
                return {
                    sites: b.sites.map(function (s) { return s.site; }),
                    up: b.up,
                    points: b.points,
                    intersections: b.intersections,
                    compound: b.compound,
                    mergeLine: b.mergeLine === undefined ? null : b.mergeLine
                };
            })
        };
    });
}

function buildSnapshot() {
    return corpusModule.corpus().map(function (c) {
        return {
            w: c.w,
            h: c.h,
            sites: serializeResult(runCase(c))
        };
    });
}

function record() {
    fs.mkdirSync(path.dirname(FIXTURE), {recursive: true});
    fs.writeFileSync(FIXTURE, JSON.stringify(buildSnapshot(), null, 1));
    console.log('recorded golden snapshot: ' + FIXTURE);
}

function compare() {
    const expected = JSON.parse(fs.readFileSync(FIXTURE, 'utf8'));
    const actual = buildSnapshot();
    assert.strictEqual(actual.length, expected.length, 'corpus size changed; re-record golden');
    let diffCount = 0;
    actual.forEach(function (snap, i) {
        try {
            assert.deepStrictEqual(snap, expected[i]);
        } catch (e) {
            diffCount++;
            if (diffCount <= 3) {
                console.error('golden mismatch in corpus case ' + i + ':\n' + e.message.slice(0, 2000));
            }
        }
    });
    if (diffCount > 0) {
        throw new Error(diffCount + ' corpus case(s) differ from golden snapshot');
    }
    console.log('golden: ' + actual.length + ' cases match snapshot');
}

module.exports = {record: record, compare: compare, serializeResult: serializeResult};

if (require.main === module) {
    const mode = process.argv[2];
    if (mode === '--record') {
        record();
    } else {
        compare();
    }
}
