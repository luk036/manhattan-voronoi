'use strict';

const test = require('node:test');
const golden = require('./golden.js');

test('JS output matches recorded golden snapshot', function () {
    golden.compare();
});
