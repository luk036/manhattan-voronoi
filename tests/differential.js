'use strict';

/**
 * Differential harness: JS library (dist) vs Python port (py_ai), run over the
 * shared corpus with nudge disabled. JS and Python currently diverge on some
 * cases (pre-existing port differences in the crop-trap filter). This tool
 * fails only when a case that previously AGREED starts diverging, or vice versa.
 *
 * Usage:
 *   node tests/differential.js               compare against recorded baseline
 *   node tests/differential.js --refresh     re-record the divergence baseline
 */

const fs = require('fs');
const os = require('os');
const path = require('path');
const cp = require('child_process');
const voronoi = require('../dist/voronoi.js');
const corpusModule = require('./corpus.js');

const BASELINE = path.join(__dirname, 'fixtures', 'diff-baseline.json');

function summarizeJs(c) {
    const sites = voronoi.generateL1Voronoi(
        c.sites.map(function (p) { return p.slice(); }), c.w, c.h, false);
    return sites.map(function (site) {
        return {
            site: site.site,
            neighbors: site.neighbors
                .map(function (n) { return [n[0], n[1]]; })
                .sort(function (a, b) { return a[0] - b[0] || a[1] - b[1]; }),
            polyCount: site.polygonPoints ? site.polygonPoints.length : null
        };
    });
}

function mapBySite(sites) {
    const m = {};
    sites.forEach(function (s) {
        m[s.site.join(',')] = {neighbors: s.neighbors, polyCount: s.polyCount};
    });
    return m;
}

function run() {
    const cases = corpusModule.corpus();
    const workDir = fs.mkdtempSync(path.join(os.tmpdir(), 'voronoi-diff-'));
    const corpusPath = path.join(workDir, 'corpus.json');
    const pyOutPath = path.join(workDir, 'py_out.json');

    fs.writeFileSync(corpusPath, JSON.stringify(cases));

    const pyScript = path.join(__dirname, 'differential_py.py');
    const proc = cp.spawnSync('python', [pyScript, corpusPath, pyOutPath],
        {encoding: 'utf8', maxBuffer: 64 * 1024 * 1024});
    if (proc.status !== 0) {
        console.error(proc.stdout);
        console.error(proc.stderr);
        fs.rmSync(workDir, {recursive: true, force: true});
        throw new Error('python differential side failed with status ' + proc.status);
    }

    const pyOut = JSON.parse(fs.readFileSync(pyOutPath, 'utf8'));
    fs.rmSync(workDir, {recursive: true, force: true});
    const jsOut = cases.map(summarizeJs);

    const diverged = [];
    cases.forEach(function (c, i) {
        const jsMap = mapBySite(jsOut[i]);
        const pyMap = mapBySite(pyOut[i]);
        const keys = new Set([].concat(Object.keys(jsMap), Object.keys(pyMap)));
        let caseDiverged = false;
        keys.forEach(function (key) {
            const a = jsMap[key];
            const b = pyMap[key];
            const na = JSON.stringify(a && a.neighbors);
            const nb = JSON.stringify(b && b.neighbors);
            const ca = a && a.polyCount;
            const cb = b && b.polyCount;
            if (na !== nb || ca !== cb) {
                caseDiverged = true;
            }
        });
        if (caseDiverged) {
            diverged.push(i);
        }
    });

    return {total: cases.length, diverged: diverged};
}

function loadBaseline() {
    return JSON.parse(fs.readFileSync(BASELINE, 'utf8'));
}

if (require.main === module) {
    const result = run();
    const mode = process.argv[2];

    if (mode === '--refresh') {
        fs.mkdirSync(path.dirname(BASELINE), {recursive: true});
        fs.writeFileSync(BASELINE, JSON.stringify(result.diverged));
        console.log('differential baseline refreshed: ' + JSON.stringify(result.diverged));
        return;
    }

    const baseline = loadBaseline();
    const newDivergences = result.diverged.filter(function (i) { return baseline.indexOf(i) < 0; });
    const healed = baseline.filter(function (i) { return result.diverged.indexOf(i) < 0; });

    console.log('differential: ' + result.total + ' cases, ' + result.diverged.length +
        ' diverging (baseline ' + baseline.length + '), healed ' + healed.length);
    if (newDivergences.length > 0) {
        console.error('REGRESSION: cases newly diverging from python port: ' +
            JSON.stringify(newDivergences));
        process.exit(1);
    }
    console.log('no new divergences');
}

module.exports = {run: run, summarizeJs: summarizeJs, loadBaseline: loadBaseline};
