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

function compile(watch) {
  console.log("rebuilding...");
  var bundler = browserify('./main.js', { debug: true });
  bundler.transform(babelify.configure({presets: ["es2015"]}));
  if (watch) { bundler = watchify(bundler); }

  function rebundle() {
    fs.mkdirSync('./build', { recursive: true });
    return bundler.bundle()
      .on('error', function(err) { console.error(err); this.emit('end'); })
      .pipe(fs.createWriteStream('./build/build.js'));
  }

  if (watch) {
    bundler.on('update', function() {
      console.log('-> bundling...');
      rebundle();
    });
  }

  return rebundle();
}

function watchSource() {
  return watch(['src/**/*.js','main.js'], function(){

      compile(true);
  });
};

gulp.task('build', function() { return compile(); });

gulp.task('watch',function(){
    return gulp.watch(['src/**/*.js','main.js'], gulp.series('build'));
});

gulp.task('watch-library',function(){
    return gulp.watch(['src/**/*.js'], gulp.series('build-library'));
});

gulp.task('default', gulp.series('build-library', 'build'));
