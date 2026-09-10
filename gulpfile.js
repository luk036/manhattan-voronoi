const fs = require('fs');
const gulp = require('gulp');
const watch = require('gulp-watch');
const sourcemaps = require('gulp-sourcemaps');
const babel = require('gulp-babel');
const concat = require('gulp-concat');
const browserify = require('browserify');
const babelify = require('babelify');
const watchify = require('watchify');
const babelPreset = require('babel-preset-es2015');

gulp.task('build-library', () => {
    return gulp.src(['src/**/*.js'])
        .pipe(sourcemaps.init())
        .pipe(babel({
            presets: ['es2015']
        }))
        .pipe(sourcemaps.write('.'))
        .pipe(gulp.dest('dist'));
});

function bundle(entry, output, watch) {
  console.log("rebuilding " + entry + "...");
  var bundler = browserify(entry, { debug: true });
  bundler.transform(babelify.configure({presets: ["es2015"]}));
  if (watch) { bundler = watchify(bundler); }

  function rebundle() {
    fs.mkdirSync('./build', { recursive: true });
    return bundler.bundle()
      .on('error', function(err) { console.error(err); this.emit('end'); })
      .pipe(fs.createWriteStream(output));
  }

  if (watch) {
    bundler.on('update', function() {
      console.log('-> bundling...');
      rebundle();
    });
  }

  return rebundle();
}

function compile(watch) {
  return bundle('./main.js', './build/build.js', watch);
}

function watchSource() {
  return watch(['src/**/*.js','main.js'], function(){

      compile(true);
  });
};

gulp.task('build', function() { return compile(); });

gulp.task('build-l2', function() { return bundle('./l2main.js', './build/l2.js', false); });

gulp.task('watch',function(){
    return gulp.watch(['src/**/*.js','main.js','l2main.js'], gulp.series('build', 'build-l2'));
});

gulp.task('watch-library',function(){
    return gulp.watch(['src/**/*.js'], gulp.series('build-library'));
});

gulp.task('default', gulp.series('build-library', 'build', 'build-l2'));
