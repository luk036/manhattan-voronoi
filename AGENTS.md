# AGENTS.md - Developer Guidelines for manhattan-voronoi

## Project Overview

A JavaScript library for generating L1 (Manhattan distance) Voronoi diagrams using Lee and Wong's algorithm. Exports ES6 modules from `src/voronoi.js`.

---

## Build, Test & Development Commands

### Install Dependencies
```bash
npm install
```

### Build Commands

| Task | Command | Description |
|------|---------|-------------|
| Build library | `npx gulp build-library` | Compiles `src/**/*.js` with Babel to `dist/` |
| Build demo | `npx gulp build` | Bundles `main.js` with Browserify to `build/build.js` |
| Full build | `npx gulp` or `npx gulp default` | Runs both library + demo build |
| Watch mode | `npx gulp watch` | Rebuild demo on file changes |
| Watch library | `npx gulp watch-library` | Rebuild library on file changes |

### Testing
```bash
# Golden-master tests (42 corpus cases) + smoke tests. Rebuilds the library first.
npm test

# JS <-> Python differential harness (py_ai port). Fails only on NEW divergences
# from the recorded baseline.
npm run test:diff

# Re-record the golden snapshot after an intentional behavior change.
npm run test:golden:record

# Re-record the differential divergence baseline after a deliberate change.
npm run test:diff:refresh

# Python port's own pytest suite
cd py_ai && python -m pytest -q
```

### Running the Demo
```bash
# Serve files locally - use any static server, e.g.:
npx serve .
# Then open http://localhost:3000 in browser
```

---

## Code Style Guidelines

### Language & Syntax
- **Language**: Vanilla JavaScript (ES6+)
- **Module System**: ES6 modules (`export`/`import`)
- **No TypeScript** - plain JS only

### Formatting
- Use 4 spaces for indentation (match existing codebase)
- No enforced formatter - keep consistent with surrounding code
- Max line length: ~100 characters (soft guideline)

### Naming Conventions

| Type | Convention | Examples |
|------|------------|-----------|
| Functions | camelCase | `generateL1Voronoi`, `cleanData` |
| Variables | camelCase | `sitePoints`, `mergeArray` |
| Constants | camelCase | No uppercase constants in this codebase |
| File names | kebab-case | `voronoi.js`, `main.js`, `gulpfile.js` |
| SVG attributes | lowercase | `setAttribute("points", ...)` |

### Imports/Exports
- Use named exports only:
  ```javascript
  export {generateL1Voronoi, cleanData};
  ```
- Import with destructuring:
  ```javascript
  import {generateL1Voronoi, cleanData} from "./src/voronoi.js";
  ```

### JSDoc Comments
- Use JSDoc for public API functions (like existing code):
  ```javascript
  /**
   * Generate an L1 Voronoi diagram
   * @param {array} sitePoints
   * @param {number} width
   * @param {number} height
   * @param {boolean} nudgeData
   * @returns {Array<Site>}
   */
  ```
- Keep existing comment style for consistency

### Error Handling
- Use descriptive `Error` messages:
  ```javascript
  throw new Error(`Square bisector: are points on a square.`);
  throw new Error(`Duplicate point: Points ${JSON.stringify(P1)} and ${JSON.stringify(P2)} are duplicates.`);
  ```
- Fail fast with clear error messages

### Data Structures

| Pattern | Usage | Example |
|---------|-------|---------|
| Points | Array `[x, y]` | `[100, 200]` |
| Colors | Array `[r, g, b]` | `[255, 128, 0]` |
| Sites | Object `{site: [x,y], bisectors: []}` | `{site: [10,20], bisectors: []}` |
| Bisectors | Object `{sites:[P1,P2], up:bool, points:[], intersections:[], compound:bool}` | `{sites:[site1,site2], up:true, points:[[0,50],[100,50]], intersections:[], compound:false}` |

### Array/Object Patterns
- Prefer `forEach`, `map`, `filter`, `reduce` over `for` loops
- Use `.sort()` with comparator functions
- Prefer concise arrow functions where readability isn't harmed

### Performance Considerations
- This is a computational geometry library - algorithm efficiency matters
- Reuse arrays where possible (avoid `new Array()` in hot paths)
- The codebase mutates some arrays in-place - document such behavior

---

## Directory Structure

```
manhattan-voronoi/
├── src/
│   ├── voronoi.js          # Library facade (public re-exports)
│   ├── generator.js        # Pipeline composition + naive generator
│   ├── geometry.js         # Pure point/segment geometry predicates
│   ├── bisector.js         # Site/Bisector factories + graph predicates
│   ├── l1Metric.js         # L1 bisector construction (metric seam)
│   ├── mergeLine.js        # Merge-line walker (Lee & Wong merge step)
│   ├── divideConquer.js    # Recursive split step
│   ├── polygonizer.js      # Bisector chaining (polygon construction)
│   ├── cell.js             # Result DTO + SVG path adapter (presentation boundary)
│   └── preprocess.js       # Input nudging (cleanData)
├── dist/                   # Built library (Babel output, mirrors src/)
├── build/
│   └── build.js            # Bundled demo
├── py_ai/                  # Python port, mirrors the src/ module map
│   ├── voronoi.py          # Facade (public re-exports)
│   ├── generator.py        # generate_l1_voronoi + naive generator
│   ├── geometry.py         # angle, distance, same_point, segment_intersection
│   ├── bisector.py         # Site/Bisector factories + graph predicates
│   ├── l1_metric.py        # L1 metric strategy (create_l1_metric)
│   ├── merge_line.py       # Merge-line walker
│   ├── divide_conquer.py   # Recursive split step
│   ├── polygonizer.py      # Bisector chaining
│   ├── cell.py             # Result DTO + SVG path adapter
│   ├── preprocess.py       # Input nudging (clean_data)
│   └── naive_oracle.py     # Brute-force oracle (test support)
├── tests/
│   ├── corpus.js           # Deterministic shared test corpus (42 cases)
│   ├── golden.js           # Golden-master snapshot tool (record/compare)
│   ├── golden.test.js      # node:test runner for golden snapshot
│   ├── smoke.test.js       # Behavior smoke tests (node:test)
│   ├── naiveOracle.js      # Brute-force oracle (moved out of the library)
│   ├── naive.test.js       # Oracle smoke test (node:test)
│   ├── differential.js     # JS <-> Python differential harness
│   └── differential_py.py  # Python side of the harness
├── main.js                 # Demo entry point
├── index.html              # Demo HTML
├── styles.css              # Demo styles
├── gulpfile.js             # Build configuration
├── package.json            # Project config
└── README.md               # Documentation
```

---

## Key Constraints

1. **No TypeScript** - plain JavaScript only
2. **No linting configured** - manually check code quality
3. **Tests live in `tests/`** (node:test + golden snapshot + Python differential)
4. **Build outputs to `dist/` and `build/`** - don't edit built files directly
5. **Compatibility** - ES6+ (no IE11 support needed)

---

## Common Development Tasks

| Task | Approach |
|------|---------|
| Add new function | Add to the matching module in `src/`, re-export from `src/voronoi.js` if public, rebuild library |
| Fix bug | Edit the owning module in `src/`, run `npm test` + `npm run test:diff` |
| Change output behavior | Deliberately re-record: `npm run test:golden:record` (and `npm run test:diff:refresh` if the Python port is affected) |
| Add a test case | Extend the shared corpus in `tests/corpus.js`, then `npm run test:golden:record` |
| Add demo feature | Edit `main.js` or `index.html` |
| Change build process | Edit `gulpfile.js` |
| Add dependencies | `npm install --save <package>`, update this file |

---

## Notes for AI Agents

- This is a math-heavy computational geometry library
- The algorithm (Lee and Wong) has specific edge cases: duplicate points, points on squares
- The code has intentional nudging behavior controlled by `nudgeData` parameter
- When modifying the algorithm, verify with the demo (various `?points=N` values)
- The demo renders bisectors with color coding by merge level - useful for debugging