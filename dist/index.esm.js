import { createFunction } from 'agora-graph';

var defaultSource = Math.random;

var randomInt = (function sourceRandomInt(source) {
  function randomInt(min, max) {
    if (arguments.length < 2) max = min, min = 0;
    min = Math.floor(min);
    max = Math.floor(max) - min;
    return function() {
      return Math.floor(source() * max + min);
    };
  }

  randomInt.source = sourceRandomInt;

  return randomInt;
})(defaultSource);

// https://en.wikipedia.org/wiki/Linear_congruential_generator#Parameters_in_common_use
const mul = 0x19660D;
const inc = 0x3C6EF35F;
const eps = 1 / 0x100000000;

function lcg(seed = Math.random()) {
  let state = (0 <= seed && seed < 1 ? seed / eps : Math.abs(seed)) | 0;
  return () => (state = mul * state + inc | 0, eps * (state >>> 0));
}

/**
 * Initialize the RNG state using the given seed.
 */
function rk_seed(seed) {
  return {
    lcg: lcg(seed)
  };
}

/**
 * Returns a random unsigned long between 0 and ULONG_MAX inclusive
 */
function rk_interval(max, state) {
  return randomInt.source(state.lcg)(max + 1)();
}

function layout(pos, sizes) {
  return {
    pos: pos,
    sizes: sizes,
    N: pos.length
  };
}
function term(i, j, d, w, o) {
  return {
    i: i,
    j: j,
    d: d,
    w: w,
    o: o
  };
}
var params = {
  K: 4,
  ALPHA: 2,
  MINIMUM_MOVEMENT: 1e-6,
  MAX_ITER: 30,
  MAX_PASSES: 100,
  SCALE_STEP: 0.1,
  distance_fn: FORBID_delta,
  delta: 0.03,
  eps: 0.01,
  seed: 0,
  PRIME: false
};

/** https://github.com/LoannGio/FORBID/blob/78e1003e017ca7b366c1c3cb427a897109d32dc9/src/FORBID.cpp#L13 */
function maxOfLeft(p1, s1, p2, s2) {
  /* double maxOfLeft(Coord const &p1, const Size &s1, Coord const &p2, const Size &s2)
     {
         double left1 = p1.x() - s1.width() / 2;
         double left2 = p2.x() - s2.width() / 2;
         return left1 > left2 ? left1 : left2;
     } */

  var left1 = p1.x - s1.width / 2;
  var left2 = p2.x - s2.width / 2;
  return left1 > left2 ? left1 : left2;
}
/** https://github.com/LoannGio/FORBID/blob/78e1003e017ca7b366c1c3cb427a897109d32dc9/src/FORBID.cpp#L20 */
function minOfRight(p1, s1, p2, s2) {
  /* double minOfRight(Coord const &p1, const Size &s1, Coord const &p2, const Size &s2)
     {
         double right1 = p1.x() + s1.width() / 2;
         double right2 = p2.x() + s2.width() / 2;
         return right1 < right2 ? right1 : right2;
     } */

  var right1 = p1.x + s1.width / 2;
  var right2 = p2.x + s2.width / 2;
  return right1 < right2 ? right1 : right2;
}

/** https://github.com/LoannGio/FORBID/blob/78e1003e017ca7b366c1c3cb427a897109d32dc9/src/FORBID.cpp#L27 */
function minOfTop(p1, s1, p2, s2) {
  /* double minOfTop(Coord const &p1, const Size &s1, Coord const &p2, const Size &s2)
     {
         double top1 = p1.y() + s1.height() / 2;
         double top2 = p2.y() + s2.height() / 2;
         return top1 < top2 ? top1 : top2;
     } */

  var top1 = p1.y + s1.height / 2;
  var top2 = p2.y + s2.height / 2;
  return top1 < top2 ? top1 : top2;
}

/** https://github.com/LoannGio/FORBID/blob/78e1003e017ca7b366c1c3cb427a897109d32dc9/src/FORBID.cpp#L34 */
function maxOfBot(p1, s1, p2, s2) {
  /* double maxOfBot(Coord const &p1, const Size &s1, Coord const &p2, const Size &s2)
     {
         double bot1 = p1.y() - s1.height() / 2;
         double bot2 = p2.y() - s2.height() / 2;
         return bot1 > bot2 ? bot1 : bot2;
     } */

  var bot1 = p1.y - s1.height / 2;
  var bot2 = p2.y - s2.height / 2;
  return bot1 > bot2 ? bot1 : bot2;
}

/** https://github.com/LoannGio/FORBID/blob/78e1003e017ca7b366c1c3cb427a897109d32dc9/src/FORBID.cpp#L41 */
function nodeRectanglesIntersection(p1, s1, p2, s2) {
  /* tuple<double, double> nodeRectanglesIntersection(Coord const &p1, const Size &s1, Coord const &p2, const Size &s2)
     {
         double min_of_right = minOfRight(p1, s1, p2, s2);
         double max_of_left = maxOfLeft(p1, s1, p2, s2);
         double min_of_top = minOfTop(p1, s1, p2, s2);
         double max_of_bot = maxOfBot(p1, s1, p2, s2);
     
         double intersection_width = min_of_right - max_of_left;
         double intersection_height = min_of_top - max_of_bot;
         return make_tuple(intersection_width, intersection_height);
     } */

  var min_of_right = minOfRight(p1, s1, p2, s2);
  var max_of_left = maxOfLeft(p1, s1, p2, s2);
  var min_of_top = minOfTop(p1, s1, p2, s2);
  var max_of_bot = maxOfBot(p1, s1, p2, s2);
  var intersection_width = min_of_right - max_of_left;
  var intersection_height = min_of_top - max_of_bot;
  return [intersection_width, intersection_height];
}

/** https://github.com/LoannGio/FORBID/blob/78e1003e017ca7b366c1c3cb427a897109d32dc9/src/FORBID.cpp#L53 */
function sort_indexes(v) {
  /* vector<size_t> sort_indexes(const vector<double> &v)
     {
     
         vector<size_t> idx(v.size());
         iota(idx.begin(), idx.end(), 0);
     
         stable_sort(idx.begin(), idx.end(),
                     [&v](size_t i1, size_t i2)
                     { return v[i1] < v[i2]; });
     
         return idx;
     } */

  var idx = Uint32Array.from(v, function (_v, i) {
    return i;
  });
  idx.sort(function (i1, i2) {
    return v[i1] - v[i2];
  });
  return idx;
}

/** https://github.com/LoannGio/FORBID/blob/78e1003e017ca7b366c1c3cb427a897109d32dc9/src/FORBID.cpp#L66 */
function sortNodesByX(X, S, sortedX, sortedS) {
  /* vector<size_t> sortNodesByX(vector<double> &X, vector<double> &S, vector<double> &sortedX, vector<double> &sortedS)
     {
         int N = X.size() / 2;
     
         vector<double> xs(N);
         for (int i = 0; i < N; ++i)
             xs[i] = X[i * 2];
         vector<size_t> idx = sort_indexes(xs);
         for (int i = 0; i < N; ++i)
         {
             sortedX[i * 2] = X[2 * idx[i]];
             sortedX[i * 2 + 1] = X[2 * idx[i] + 1];
         
             sortedS[i * 2] = S[2 * idx[i]];
             sortedS[i * 2 + 1] = S[2 * idx[i] + 1];
         }
         return idx;
     } */
  var N = X.length / 2;
  var xs = new Array(N);
  for (var i = 0; i < N; i++) xs[i] = X[i * 2];
  var idx = sort_indexes(xs);
  for (var _i = 0; _i < N; _i++) {
    sortedX[_i * 2] = X[2 * idx[_i]];
    sortedX[_i * 2 + 1] = X[2 * idx[_i] + 1];
    sortedS[_i * 2] = S[2 * idx[_i]];
    sortedS[_i * 2 + 1] = S[2 * idx[_i] + 1];
  }
  return idx;
}

/** https://github.com/LoannGio/FORBID/blob/78e1003e017ca7b366c1c3cb427a897109d32dc9/src/FORBID.cpp#L85
 * @use agora-graph#overlap
 * @see overlap
 */
function overlapCheck(p1, p2, s1, s2) {
  /* bool overlapCheck(Coord const &p1, Coord const &p2, const Size &s1, const Size &s2)
     {
         auto [i_w, i_h] = nodeRectanglesIntersection(p1, s1, p2, s2);
         if (i_w > 0 && i_h > 0)
             return true;
         return false;
     } */

  var rects = nodeRectanglesIntersection(p1, s1, p2, s2);
  var i_w = rects[0];
  var i_h = rects[1];
  if (i_w > 0 && i_h > 0) {
    return true;
  }
  return false;
}

/**
 *  https://github.com/LoannGio/FORBID/blob/78e1003e017ca7b366c1c3cb427a897109d32dc9/src/FORBID.cpp#L93
 * @see hasOverlap
 */
function scanLineOverlapCheck(X, S) {
  /* bool scanLineOverlapCheck(vector<double> &X, vector<double> &S)
     {
         int N = X.size() / 2;
         vector<double> sortedX(2 * N);
         vector<double> sortedS(2 * N);
         sortNodesByX(X, S, sortedX, sortedS);
     
         Coord p1, p2;
         Size s1, s2;
         double left_j, right_i;
         for (int i = 0; i < N; ++i)
         {
             p1.set(
                 sortedX[i * 2],
                 sortedX[i * 2 + 1]);
             s1.set(
                 sortedS[i * 2],
                 sortedS[i * 2 + 1]);
             
             right_i = p1.x() + s1.width() / 2;
             for (int j = i + 1; j < N; ++j)
             {
                 p2.set(sortedX[j * 2], sortedX[j * 2 + 1]);
                 s2.set(sortedS[j * 2], sortedS[j * 2 + 1]);
             
                 left_j = p2.x() - s2.width() / 2;
                 if (overlapCheck(p1, p2, s1, s2))
                     return true;
                 else if (left_j > right_i)
                     break;
             }
         }
         return false;
     } */

  var N = X.length / 2;
  var sortedX = new Array(2 * N);
  var sortedS = new Array(2 * N);
  sortNodesByX(X, S, sortedX, sortedS);
  for (var i = 0; i < N; i++) {
    var p1 = {
      x: sortedX[i * 2],
      y: sortedX[i * 2 + 1]
    };
    var s1 = {
      width: sortedS[i * 2],
      height: sortedS[i * 2 + 1]
    };
    var right_i = p1.x + s1.width / 2;
    for (var j = i + 1; j < N; j++) {
      var p2 = {
        x: sortedX[j * 2],
        y: sortedX[j * 2 + 1]
      };
      var s2 = {
        width: sortedS[j * 2],
        height: sortedS[j * 2 + 1]
      };
      var left_j = p2.x - s2.width / 2;
      if (overlapCheck(p1, p2, s1, s2)) return true;else if (left_j > right_i) break;
    }
  }
  return false;
}

/**
 *  https://github.com/LoannGio/FORBID/blob/78e1003e017ca7b366c1c3cb427a897109d32dc9/src/FORBID.cpp#L128 */
function getAllOverlaps(X, S) {
  /* vector<tuple<int, int>> getAllOverlaps(vector<double> &X, vector<double> &S)
     {
     
         int N = X.size() / 2;
     
         vector<double> sortedX(X);
         vector<double> sortedS(S);
         vector<size_t> mapping = sortNodesByX(X, S, sortedX, sortedS);
     
         Coord p1, p2;
         Size s1, s2;
         double left_j, right_i;
         vector<tuple<int, int>> overlaps;
         for (int i = 0; i < N; ++i)
         {
             p1.set(sortedX[i * 2], sortedX[i * 2 + 1]);
             s1.set(sortedS[i * 2], sortedS[i * 2 + 1]);
             right_i = p1.x() + s1.width() / 2;
             for (int j = i + 1; j < N; ++j)
             {
                 p2.set(sortedX[j * 2], sortedX[j * 2 + 1]);
                 s2.set(sortedS[j * 2], sortedS[j * 2 + 1]);
                 left_j = p2.x() - s2.width() / 2;
                 if (overlapCheck(p1, p2, s1, s2))
                     overlaps.emplace_back(make_tuple(mapping[i], mapping[j]));
                 else if (left_j > right_i)
                     break;
             }
         }
         return overlaps;
     } */

  var N = X.length / 2;
  var sortedX = new Array(2 * N);
  var sortedS = new Array(2 * N);
  var mapping = sortNodesByX(X, S, sortedX, sortedS);
  var overlaps = [];
  for (var i = 0; i < N; i++) {
    var p1 = {
      x: sortedX[i * 2],
      y: sortedX[i * 2 + 1]
    };
    var s1 = {
      width: sortedS[i * 2],
      height: sortedS[i * 2 + 1]
    };
    var right_i = p1.x + s1.width / 2;
    for (var j = i + 1; j < N; j++) {
      var p2 = {
        x: sortedX[j * 2],
        y: sortedX[j * 2 + 1]
      };
      var s2 = {
        width: sortedS[j * 2],
        height: sortedS[j * 2 + 1]
      };
      var left_j = p2.x - s2.width / 2;
      if (overlapCheck(p1, p2, s1, s2)) overlaps.push([mapping[i], mapping[j]]);else if (left_j > right_i) break;
    }
  }
  return overlaps;
}

/** https://github.com/LoannGio/FORBID/blob/78e1003e017ca7b366c1c3cb427a897109d32dc9/src/FORBID.cpp#L160 */
function FORBID_delta(ci, cj, si, sj, intersec_width, intersec_height) {
  /* double FORBID_delta(Coord const &ci, Coord const &cj, const Size &si, const Size &sj, double intersec_width, double intersec_height)
     {
         return sqrt(pow(si.width() / 2 + sj.width() / 2, 2.) + pow(si.height() / 2 + sj.height() / 2, 2.));
     } */
  return Math.sqrt(Math.pow(si.width / 2 + sj.width / 2, 2) + Math.pow(si.height / 2 + sj.height / 2, 2));
}

/** https://github.com/LoannGio/FORBID/blob/78e1003e017ca7b366c1c3cb427a897109d32dc9/src/FORBID.cpp#L180 */
function eucl(x1, y1, x2, y2) {
  /* double eucl(double x1, double y1, double x2, double y2)
  {
    return sqrt(pow(x1 - x2, 2.) + pow(y1 - y2, 2.));
  } */

  return Math.sqrt(Math.pow(x1 - x2, 2) + Math.pow(y1 - y2, 2));
}

/** https://github.com/LoannGio/FORBID/blob/78e1003e017ca7b366c1c3cb427a897109d32dc9/src/FORBID.cpp#L185 */
function euclCoord(p1, p2) {
  /* double eucl(Coord const &p1, Coord const &p2)
     {
         return sqrt(pow(p1.x() - p2.x(), 2.) + pow(p1.y() - p2.y(), 2.));
     } */

  return Math.sqrt(Math.pow(p1.x - p2.x, 2) + Math.pow(p1.y - p2.y, 2));
}

/** https://github.com/LoannGio/FORBID/blob/78e1003e017ca7b366c1c3cb427a897109d32dc9/src/FORBID.cpp#L190 */
function vecNorm2D(vec_x, vec_y) {
  /* double vecNorm2D(double vec_x, double vec_y)
     {
         return sqrt(vec_x * vec_x + vec_y * vec_y);
     } */
  return Math.sqrt(vec_x * vec_x + vec_y * vec_y);
}

/** https://github.com/LoannGio/FORBID/blob/78e1003e017ca7b366c1c3cb427a897109d32dc9/src/FORBID.cpp#L195 */
function maxScaleRatio(X, S) {
  /* double maxScaleRatio(vector<double> &X, vector<double> &S)
     {
         double padding = 1e-4;
         double maxRatio = 1.;
         double optimalDist, actualDist, ratio, unoverlapRatio;
         double actualX, actualY, desiredWidth, desiredHeight, widthRatio, heightRatio;
         auto overlaps = getAllOverlaps(X, S);
         for (unsigned int i = 0; i < overlaps.size(); ++i)
         {
             auto [u, v] = overlaps[i];
             actualDist = eucl(X[u * 2], X[u * 2 + 1], X[v * 2], X[v * 2 + 1]);
         
             actualX = X[u * 2] - X[v * 2];
             actualY = X[u * 2 + 1] - X[v * 2 + 1];
             desiredWidth = (S[u * 2] + S[v * 2]) / 2 + padding;
             desiredHeight = (S[u * 2 + 1] + S[v * 2 + 1]) / 2 + padding;
         
             widthRatio = desiredWidth / actualX;
             heightRatio = desiredHeight / actualY;
         
             unoverlapRatio = min(abs(widthRatio), abs(heightRatio));
             actualX *= unoverlapRatio;
             actualY *= unoverlapRatio;
         
             optimalDist = vecNorm2D(actualX, actualY);
             ratio = optimalDist / actualDist;
         
             // ratio = max(
             // 1.,
             // min((S[u*2] / 2 + S[v*2] / 2 + padding) / abs(X[u*2] - X[v*2]),
             //     (S[u*2+1] / 2 + S[v*2+1] / 2 + padding) / abs(X[u*2+1] - X[v*2+1])));
         
             maxRatio = max(maxRatio, ratio);
         }
         return maxRatio;
     } */

  var padding = 1e-4;
  var maxRatio = 1;
  var overlaps = getAllOverlaps(X, S);
  for (var i = 0; i < overlaps.length; ++i) {
    var uv = overlaps[i];
    var u = uv[0];
    var v = uv[1];
    var actualDist = eucl(X[u * 2], X[u * 2 + 1], X[v * 2], X[v * 2 + 1]);
    var actualX = X[u * 2] - X[v * 2];
    var actualY = X[u * 2 + 1] - X[v * 2 + 1];
    var desiredWidth = (S[u * 2] + S[v * 2]) / 2 + padding;
    var desiredHeight = (S[u * 2 + 1] + S[v * 2 + 1]) / 2 + padding;
    var widthRatio = desiredWidth / actualX;
    var heightRatio = desiredHeight / actualY;
    var unoverlapRatio = Math.min(Math.abs(widthRatio), Math.abs(heightRatio));
    actualX *= unoverlapRatio;
    actualY *= unoverlapRatio;
    var optimalDist = vecNorm2D(actualX, actualY);
    var ratio = optimalDist / actualDist;
    maxRatio = Math.max(maxRatio, ratio);
  }
  return maxRatio;
}

/** https://github.com/LoannGio/FORBID/blob/78e1003e017ca7b366c1c3cb427a897109d32dc9/src/FORBID.cpp#L232 */
function scaleLayout(X, scaleFactor) {
  /* void scaleLayout(vector<double> &X, double scaleFactor)
     {
         int N = X.size() / 2;
         for (int i = 0; i < N; ++i)
         {
             X[i * 2] *= scaleFactor;
             X[i * 2 + 1] *= scaleFactor;
         }
     } */
  var N = X.length / 2;
  for (var i = 0; i < N; ++i) {
    X[i * 2] *= scaleFactor;
    X[i * 2 + 1] *= scaleFactor;
  }
}

/** https://github.com/LoannGio/FORBID/blob/78e1003e017ca7b366c1c3cb427a897109d32dc9/src/FORBID.cpp#L242 */
function isCurrentScaleSolvable(X, S) {
  /* bool isCurrentScaleSolvable(vector<double> &X, vector<double> &S)
     {
         int N = X.size() / 2;
         double areas_sum = 0;
         double min_x = DBL_MAX;
         double min_y = DBL_MAX;
         double max_x = -DBL_MAX;
         double max_y = -DBL_MAX;
     
         double left, right, top, bot;
         for (int i = 0; i < N; ++i)
         {
             areas_sum += S[i * 2] * S[i * 2 + 1];
             left = X[i * 2] - S[i * 2] / 2;
             right = X[i * 2] + S[i * 2] / 2;
             top = X[i * 2 + 1] + S[i * 2 + 1] / 2;
             bot = X[i * 2 + 1] - S[i * 2 + 1] / 2;
             if (left < min_x)
                 min_x = left;
             if (right > max_x)
                 max_x = right;
             if (bot < min_y)
                 min_y = bot;
             if (top > max_y)
                 max_y = top;
         }
         double bb_area = (max_x - min_x) * (max_y - min_y);
       
         return bb_area >= areas_sum;
     } */

  var N = X.length / 2;
  var areas_sum = 0;
  var min_x = Number.MAX_VALUE;
  var min_y = Number.MAX_VALUE;
  var max_x = -Number.MAX_VALUE;
  var max_y = -Number.MAX_VALUE;
  for (var i = 0; i < N; ++i) {
    areas_sum += S[i * 2] * S[i * 2 + 1];
    var left = X[i * 2] - S[i * 2] / 2;
    var right = X[i * 2] + S[i * 2] / 2;
    var top = X[i * 2 + 1] + S[i * 2 + 1] / 2;
    var bot = X[i * 2 + 1] - S[i * 2 + 1] / 2;
    if (left < min_x) min_x = left;
    if (right > max_x) max_x = right;
    if (bot < min_y) min_y = bot;
    if (top > max_y) max_y = top;
  }
  var bb_area = (max_x - min_x) * (max_y - min_y);
  return bb_area >= areas_sum;
}

/** https://github.com/LoannGio/FORBID/blob/78e1003e017ca7b366c1c3cb427a897109d32dc9/src/FORBID.cpp#L289 */
function loadParams(defaultParams, params) {
  var _params$K, _params$ALPHA, _params$MINIMUM_MOVEM, _params$MAX_ITER, _params$MAX_PASSES, _params$SCALE_STEP, _params$distance_fn, _params$delta, _params$eps, _params$seed, _params$PRIME;
  /* void loadParams(Parametrizer &params, char *argv[])
     {
         params.ALPHA = atof(argv[2]);
         params.K = atof(argv[3]);
         params.MINIMUM_MOVEMENT = atof(argv[4]);
         params.MAX_ITER = atoi(argv[5]);
         params.MAX_PASSES = atoi(argv[6]);
         params.SCALE_STEP = atof(argv[7]);
         params.PRIME = atoi(argv[8]);
     } */

  return {
    K: (_params$K = params.K) !== null && _params$K !== void 0 ? _params$K : defaultParams.K,
    ALPHA: (_params$ALPHA = params.ALPHA) !== null && _params$ALPHA !== void 0 ? _params$ALPHA : defaultParams.ALPHA,
    MINIMUM_MOVEMENT: (_params$MINIMUM_MOVEM = params.MINIMUM_MOVEMENT) !== null && _params$MINIMUM_MOVEM !== void 0 ? _params$MINIMUM_MOVEM : defaultParams.MINIMUM_MOVEMENT,
    MAX_ITER: (_params$MAX_ITER = params.MAX_ITER) !== null && _params$MAX_ITER !== void 0 ? _params$MAX_ITER : defaultParams.MAX_ITER,
    MAX_PASSES: (_params$MAX_PASSES = params.MAX_PASSES) !== null && _params$MAX_PASSES !== void 0 ? _params$MAX_PASSES : defaultParams.MAX_PASSES,
    SCALE_STEP: (_params$SCALE_STEP = params.SCALE_STEP) !== null && _params$SCALE_STEP !== void 0 ? _params$SCALE_STEP : defaultParams.SCALE_STEP,
    distance_fn: (_params$distance_fn = params.distance_fn) !== null && _params$distance_fn !== void 0 ? _params$distance_fn : defaultParams.distance_fn,
    delta: (_params$delta = params.delta) !== null && _params$delta !== void 0 ? _params$delta : defaultParams.delta,
    eps: (_params$eps = params.eps) !== null && _params$eps !== void 0 ? _params$eps : defaultParams.eps,
    seed: (_params$seed = params.seed) !== null && _params$seed !== void 0 ? _params$seed : defaultParams.seed,
    PRIME: (_params$PRIME = params.PRIME) !== null && _params$PRIME !== void 0 ? _params$PRIME : defaultParams.PRIME
  };
}

/** https://github.com/LoannGio/FORBID/blob/78e1003e017ca7b366c1c3cb427a897109d32dc9/src/FORBID.cpp#L300 */
function initLayout(g) {
  /* vector<double> initLayout(Layout const &g)
     {
         vector<double> X(2 * g.N());
         for (unsigned int i = 0; i < g.N(); i++)
         {
             X[i * 2] = g.pos[i].x();
             X[i * 2 + 1] = g.pos[i].y();
         }
         return X;
     } */
  var X = new Array(2 * g.N);
  for (var i = 0; i < g.N; i++) {
    X[i * 2] = g.pos[i].x;
    X[i * 2 + 1] = g.pos[i].y;
  }
  return X;
}

/** https://github.com/LoannGio/FORBID/blob/78e1003e017ca7b366c1c3cb427a897109d32dc9/src/FORBID.cpp#L311 */
function initSizes(g) {
  /* vector<double> initSizes(Layout const &g)
     {
         vector<double> S(2 * g.N());
         for (unsigned int i = 0; i < g.N(); i++)
         {
             S[i * 2] = g.sizes[i].width();
             S[i * 2 + 1] = g.sizes[i].height();
         }
         return S;
     } */
  var S = new Array(2 * g.N);
  for (var i = 0; i < g.N; i++) {
    S[i * 2] = g.sizes[i].width;
    S[i * 2 + 1] = g.sizes[i].height;
  }
  return S;
}

/** https://github.com/LoannGio/FORBID/blob/78e1003e017ca7b366c1c3cb427a897109d32dc9/src/FORBID.cpp#L322 */
function copyLayout(ins, out) {
  /* void copyLayout(vector<double> const &in, vector<double> &out)
     {
         auto N = out.size() / 2;
         for (unsigned int i = 0; i < 2 * N; ++i)
             out[i] = in[i];
     } */
  var N = out.length / 2;
  for (var i = 0; i < 2 * N; ++i) out[i] = ins[i];
}

/** https://github.com/LoannGio/FORBID/blob/78e1003e017ca7b366c1c3cb427a897109d32dc9/src/FORBID.cpp#L356 */
function schedule(terms, t_max, eps) {
  /* vector<double> schedule(const vector<term> &terms, int t_max, double eps)
     {
         double w_min = terms[0].w, w_max = terms[0].w;
         for (unsigned i = 1; i < terms.size(); i++)
         {
             const double &w = terms[i].w;
             if (w < w_min)
                 w_min = w;
             if (w > w_max)
                 w_max = w;
         }
         double eta_max = 1.0 / w_min;
         double eta_min = eps / w_max;
         double lambda = log(eta_max / eta_min) / ((double)t_max - 1);
       
         // initialize step sizes
         vector<double> etas;
         etas.reserve(t_max);
         for (int t = 0; t < t_max; t++)
             etas.push_back(eta_max * exp(-lambda * t));
       
         return etas;
     } */
  var w_min = terms[0].w,
    w_max = terms[0].w;
  for (var i = 1; i < terms.length; i++) {
    var w = terms[i].w;
    if (w < w_min) w_min = w;
    if (w > w_max) w_max = w;
  }
  var eta_max = 1.0 / w_min;
  var eta_min = eps / w_max;
  var lambda = Math.log(eta_max / eta_min) / (t_max - 1);

  // initialize step sizes
  var etas = [];
  for (var t = 0; t < t_max; t++) etas.push(eta_max * Math.exp(-lambda * t));
  return etas;
}

/** https://github.com/LoannGio/FORBID/blob/78e1003e017ca7b366c1c3cb427a897109d32dc9/src/FORBID.cpp#L380 */
function fisheryates_shuffle(terms, rstate) {
  /* void fisheryates_shuffle(vector<term> &terms, rk_state &rstate)
     {
         int n = terms.size();
         for (unsigned i = n - 1; i >= 1; i--)
         {
             unsigned j = rk_interval(i, &rstate);
             term temp = terms[i];
             terms[i] = terms[j];
             terms[j] = temp;
         }
     } */
  var n = terms.length;
  for (var i = n - 1; i >= 1; i--) {
    var j = rk_interval(i, rstate);
    var temp = terms[i];
    terms[i] = terms[j];
    terms[j] = temp;
  }
}

/** https://github.com/LoannGio/FORBID/blob/78e1003e017ca7b366c1c3cb427a897109d32dc9/src/FORBID.cpp#L392 */
function layoutToTerms(X, init_X, S, params) {
  /* vector<term> layoutToTerms(vector<double> &X, vector<double> &init_X, vector<double> &S, Parametrizer &params)
     {
         int N = X.size() / 2;
         vector<term> terms;
         double d_ij, w_ij; //, d_eucl;
         double xi, xj, yi, yj;
         Coord ci, cj;
         Coord init_ci, init_cj;
         Size si, sj;
         bool overlap;
         for (int i = 0; i < N; i++)
         {
             ci.set(X[i * 2], X[i * 2 + 1]);
             init_ci.set(init_X[i * 2], init_X[i * 2 + 1]);
             si.set(S[i * 2], S[i * 2 + 1]);
         
             for (int j = i + 1; j < N; j++)
             {
                 if (i != j)
                 {
                     cj.set(X[j * 2], X[j * 2 + 1]);
                     init_cj.set(init_X[j * 2], init_X[j * 2 + 1]);
                     sj.set(S[j * 2], S[j * 2 + 1]);
                     auto [i_w, i_h] = nodeRectanglesIntersection(ci, si, cj, sj);
                 
                     if (i_w > 0 && i_h > 0) // there is overlap
                     {
                         d_ij = (*params.distance_fn)(init_ci, init_cj, si, sj, i_w, i_h);
                         w_ij = pow(d_ij, params.K * params.ALPHA);
                         overlap = true;
                     }
                     else
                     {
                         d_ij = eucl(ci, cj);
                         w_ij = pow(d_ij, params.ALPHA);
                         overlap = false;
                     }
                     terms.push_back(term(i, j, d_ij, w_ij, overlap));
                 }
             }
         }
         return terms;
     } */

  var N = X.length / 2;
  var terms = [];
  var d_ij, w_ij; //, d_eucl;
  var overlap;
  for (var i = 0; i < N; i++) {
    var _ci = {
      x: X[i * 2],
      y: X[i * 2 + 1]
    };
    var init_ci = {
      x: init_X[i * 2],
      y: init_X[i * 2 + 1]
    };
    var _si = {
      width: S[i * 2],
      height: S[i * 2 + 1]
    };
    for (var j = i + 1; j < N; j++) {
      if (i != j) {
        var _cj = {
          x: X[j * 2],
          y: X[j * 2 + 1]
        };
        var init_cj = {
          x: init_X[j * 2],
          y: init_X[j * 2 + 1]
        };
        var _sj = {
          width: S[j * 2],
          height: S[j * 2 + 1]
        };
        var rects = nodeRectanglesIntersection(_ci, _si, _cj, _sj);
        var i_w = rects[0];
        var i_h = rects[1];
        if (i_w > 0 && i_h > 0) {
          // there is overlap
          d_ij = params.distance_fn(init_ci, init_cj, _si, _sj, i_w, i_h);
          w_ij = Math.pow(d_ij, params.K * params.ALPHA);
          overlap = true;
        } else {
          d_ij = euclCoord(_ci, _cj);
          w_ij = Math.pow(d_ij, params.ALPHA);
          overlap = false;
        }
        terms.push(term(i, j, d_ij, w_ij, overlap));
      }
    }
  }
  return terms;
}

/** https://github.com/LoannGio/FORBID/blob/78e1003e017ca7b366c1c3cb427a897109d32dc9/src/FORBID.cpp#L437 */
function OPTIMIZATION_PASS(X, terms, init_X, S, etas, params) {
  /* // S_GD2 optim algorithm, adapted from https://github.com/jxz12/s_gd2/blob/master/cpp/s_gd2/layout.cpp
     void OPTIMIZATION_PASS(vector<double> &X, vector<term> &terms, vector<double> &init_X, vector<double> &S, const vector<double> &etas, Parametrizer &params)
     {
         rk_state rstate;
         rk_seed(params.seed, &rstate);
         double mvt_sum;
         unsigned int i_eta;
         for (i_eta = 0; i_eta < etas.size(); i_eta++)
         {
             const double eta = etas[i_eta];
         
             unsigned n_terms = terms.size();
             if (n_terms == 0)
                 return;
             fisheryates_shuffle(terms, rstate);
         
             mvt_sum = 0;
             for (unsigned i_term = 0; i_term < n_terms; i_term++)
             {
                 const term &t = terms[i_term];
                 const int &i = t.i, &j = t.j;
                 const double &w_ij = t.w;
                 const double &d_ij = t.d;
                 if (true || t.o)
                 {
                     double mu = eta * w_ij;
                     if (mu > 1)
                         mu = 1;
                     double dx = X[i * 2] - X[j * 2];
                     double dy = X[i * 2 + 1] - X[j * 2 + 1];
                     double mag = sqrt(dx * dx + dy * dy);
                 
                     double r = 0;
                     if (mag != 0)
                         r = (mu * (mag - d_ij)) / (2 * mag);
                     double r_x = r * dx;
                     double r_y = r * dy;
                     mvt_sum += abs(r_x) + abs(r_y);
                     X[i * 2] -= r_x;
                     X[i * 2 + 1] -= r_y;
                     X[j * 2] += r_x;
                     X[j * 2 + 1] += r_y;
                 }
             }
             if (mvt_sum < params.MINIMUM_MOVEMENT)
             {
                 return;
             }
             terms = layoutToTerms(X, init_X, S, params);
         }
         return;
     } */

  var rstate = rk_seed(params.seed);
  for (var i_eta = 0; i_eta < etas.length; i_eta++) {
    var eta = etas[i_eta];
    var n_terms = terms.length;
    if (n_terms === 0) return;
    fisheryates_shuffle(terms, rstate);
    var mvt_sum = 0;
    for (var i_term = 0; i_term < n_terms; i_term++) {
      var t = terms[i_term];
      var i = t.i,
        j = t.j;
      var w_ij = t.w;
      var d_ij = t.d;
      {
        var mu = eta * w_ij;
        if (mu > 1) mu = 1;
        var dx = X[i * 2] - X[j * 2];
        var dy = X[i * 2 + 1] - X[j * 2 + 1];
        var mag = Math.sqrt(dx * dx + dy * dy);
        var r = 0;
        if (mag !== 0) r = mu * (mag - d_ij) / (2 * mag);
        var r_x = r * dx;
        var r_y = r * dy;
        mvt_sum += Math.abs(r_x) + Math.abs(r_y);
        X[i * 2] -= r_x;
        X[i * 2 + 1] -= r_y;
        X[j * 2] += r_x;
        X[j * 2 + 1] += r_y;
      }
    }
    if (mvt_sum < params.MINIMUM_MOVEMENT) return;
    terms = layoutToTerms(X, init_X, S, params);
  }
}
/** https://github.com/LoannGio/FORBID/blob/78e1003e017ca7b366c1c3cb427a897109d32dc9/src/FORBID.cpp#L356 */
function passInOptim(X, initX, S, params) {
  /* void passInOptim(vector<double> &X, vector<double> &initX, vector<double> &S, Parametrizer &params)
     {
         vector<term> orig_terms = layoutToTerms(X, initX, S, params);
         vector<double> etas = schedule(orig_terms, params.MAX_ITER, params.eps);
         OPTIMIZATION_PASS(X, orig_terms, initX, S, etas, params);
     } */
  var orig_terms = layoutToTerms(X, initX, S, params);
  var etas = schedule(orig_terms, params.MAX_ITER, params.eps);
  OPTIMIZATION_PASS(X, orig_terms, initX, S, etas, params);
}

/** https://github.com/LoannGio/FORBID/blob/78e1003e017ca7b366c1c3cb427a897109d32dc9/src/FORBID.cpp#L496 */
function main(g, options) {
  /* int main(int argc, char *argv[])
     {
         if (argc != 9)
         {
             cerr << "Wrong parameters" << endl;
             return EXIT_FAILURE;
         }
         string input_layout_path = argv[1];
         Layout g = loadOriginalDistance(input_layout_path);
       
         Parametrizer params;
         loadParams(params, argv);
       
         auto init_X = initLayout(g);
         vector<double> scaled_init_X(g.N() * 2);
         auto S = initSizes(g);
         vector<double> X(g.N() * 2);
         copyLayout(init_X, scaled_init_X);
         copyLayout(init_X, X);
       
         int n_passes = 0;
         double scaleFactor;
       
         double maxRatio = maxScaleRatio(X, S);
         double upperScale = maxRatio;
         double lowerScale = 1.;
         double oldScale = 1.;
         double curScale = 1.;
       
         auto start = std::chrono::steady_clock::now();
         bool stop = !scanLineOverlapCheck(X, S); // do not enter if there is no overlap
       
         // if possible, try to solve the problem in current scale
         if (!stop && isCurrentScaleSolvable(X, S))
         {
             passInOptim(X, init_X, S, params);
             stop = !scanLineOverlapCheck(X, S);
         }
       
         while (!stop)
         {
             curScale = (upperScale + lowerScale) / 2;
             scaleFactor = curScale / oldScale;
             oldScale = curScale;
             scaleLayout(scaled_init_X, scaleFactor);
             if (params.PRIME)
                 copyLayout(scaled_init_X, X);
             else
                 scaleLayout(X, scaleFactor);
         
             passInOptim(X, scaled_init_X, S, params);
         
             if (scanLineOverlapCheck(X, S))
                 lowerScale = curScale;
             else // no overlap
             {
                 if (upperScale - lowerScale < params.SCALE_STEP)
                     stop = true;
                 else
                     upperScale = curScale;
             }
           
             ++n_passes;
             if (n_passes >= params.MAX_PASSES)
                 break;
         }
         auto elapsed = std::chrono::duration<double, std::milli>     (std::chrono::steady_clock::now() - start).count();
         string savepath = save(X, S, input_layout_path, elapsed, params.PRIME);
       
         cout << "DONE; saved to " << savepath << endl;
         return EXIT_SUCCESS;
     }*/

  var init_X = initLayout(g);
  var scaled_init_X = new Array(g.N * 2);
  var S = initSizes(g);
  var X = new Array(g.N * 2);
  copyLayout(init_X, scaled_init_X);
  copyLayout(init_X, X);
  var n_passes = 0;
  var maxRatio = maxScaleRatio(X, S);
  var upperScale = maxRatio;
  var lowerScale = 1;
  var oldScale = 1;
  var hasOverlap = scanLineOverlapCheck(X, S); // do not enter if there is no overlap

  // if possible, try to solve the problem in current scale
  if (hasOverlap && isCurrentScaleSolvable(X, S)) {
    passInOptim(X, init_X, S, params);
    hasOverlap = scanLineOverlapCheck(X, S);
  }
  while (hasOverlap) {
    var curScale = (upperScale + lowerScale) / 2;
    var scaleFactor = curScale / oldScale;
    oldScale = curScale;
    scaleLayout(scaled_init_X, scaleFactor);
    if (params.PRIME) copyLayout(scaled_init_X, X);else scaleLayout(X, scaleFactor);
    passInOptim(X, scaled_init_X, S, params);
    if (scanLineOverlapCheck(X, S)) lowerScale = curScale;
    // no overlap
    else {
      if (upperScale - lowerScale < params.SCALE_STEP) hasOverlap = false;else upperScale = curScale;
    }
    ++n_passes;
    if (n_passes >= params.MAX_PASSES) break;
  }
  return [X, S];
}

/**
 * Implementation of FORBID algorithm
 *
 * FORBID: Fast Overlap Removal By Stochastic Gradient Descent for Graph Drawing,
 * https://github.com/LoannGio/FORBID
 */

/**
 * Executes the FORBID algorithm for this graph
 *
 * @param {Graph} graph the graph to update
 * @param {object} [options] options to pass to the algorith
 * @param {number} options.padding padding to add between nodes
 *
 * @returns {Result} the updated graph
 */
var forbid = createFunction(function (graph) {
  var options = arguments.length > 1 && arguments[1] !== undefined ? arguments[1] : {};
  /* https://github.com/LoannGio/FORBID/blob/78e1003e017ca7b366c1c3cb427a897109d32dc9/src/FORBID.cpp#L496 */

  var g = layout(graph.nodes, graph.nodes);
  var result = main(g, loadParams(params, options));
  var X = result[0];
  var S = result[1];
  return {
    graph: {
      nodes: graph.nodes.map(function (v, i) {
        return {
          x: X[i * 2],
          y: X[i * 2 + 1],
          width: S[i * 2],
          height: S[i * 2 + 1],
          index: v.index,
          label: v.label,
          meta: v.meta,
          up: v.up
        };
      }),
      edges: graph.edges
    }
  };
});
var FORBIDAlgorithm = {
  name: 'FORBID',
  algorithm: forbid
};

export { FORBIDAlgorithm, FORBIDAlgorithm as default, forbid };
