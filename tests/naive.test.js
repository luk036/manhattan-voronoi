'use strict';

const test = require('node:test');
const assert = require('node:assert');
const {generateVoronoiPoints} = require('./naiveOracle.js');

function l1(a, b){
    return Math.abs(a[0] - b[0]) + Math.abs(a[1] - b[1]);
}

test('naive oracle assigns a color to every pixel', function () {
    const result = generateVoronoiPoints([[0, 0], [9, 0]], 10, 10, l1);
    assert.strictEqual(result.length, 100);
    result.forEach(function (color) {
        assert.ok(Array.isArray(color) && color.length === 3);
        color.forEach(function (channel) {
            assert.ok(Number.isInteger(channel) && channel >= 0 && channel <= 255);
        });
    });
});
