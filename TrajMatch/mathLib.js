var mathLib = {}
let pi = 3.141592654
var erf = require('math-erf')
var seedrandom = require('seedrandom')
var rng = seedrandom('43553')

mathLib.pnorm = function (x, mu = 0, sd = 1, lower_tail = true, give_log = false) {
  if (sd < 0) {
    return NaN
  }
  let ans = 1 / 2 * (1 + erf((x - mu) / sd / Math.sqrt(2)))
  if (!lower_tail) {
    ans = 1 - ans
  }
  if (give_log) {
    ans = Math.log(ans)
  }
  return ans
}

mathLib.rnorm = function (mu = 0, sd = 1) {
  var val = Math.sqrt(-2 * Math.log(rng())) * Math.cos(2 * pi * rng())
  return val * sd + mu
}

mathLib.dpois = function (x, lambda) {
  let ans, total = 0
  if (isNaN(x) || isNaN(lambda) || lambda < 0) {
    return NaN
  }
  if (!Number.isInteger(x)) {
    return 0
  }
  if (x < 0 || !isFinite(x)) {
    return 0
  }
  x = Math.round(x)
  ans = -lambda + x * Math.log(lambda)
  for (let i = 1; i <= x; i++) {
    total += Math.log(i)
  }
  let logAns = ans - total
  return Math.exp(logAns)
}
mathLib.interpolator = function (points) {
  var first, n = points.length - 1,
    interpolated,
    leftExtrapolated,
    rightExtrapolated;

  if (points.length === 0) {
    return function () {
      return 0
    }
  }

  if (points.length === 1) {
    return function () {
      return points[0][1]
    }
  }

  points = points.sort(function (a, b) {
    return a[0] - b[0]
  })
  first = points[0]

  leftExtrapolated = function (x) {
    var a = points[0], b = points[1];
    return a[1] + (x - a[0]) * (b[1] - a[1]) / (b[0] - a[0])
  }

  interpolated = function (x, a, b) {
    return a[1] + (x - a[0]) * (b[1] - a[1]) / (b[0] - a[0])
  }

  rightExtrapolated = function (x) {
    var a = points[n - 1], b = points[n];
    return b[1] + (x - b[0]) * (b[1] - a[1]) / (b[0] - a[0])
  }

  return function (x) {
    var i
    if (x <= first[0]) {
      return leftExtrapolated(x)
    }
    for (i = 0; i < n; i += 1) {
      if (x > points[i][0] && x <= points[i + 1][0]) {
        return interpolated(x, points[i], points[i + 1])
      }
    }
    return rightExtrapolated(x);
  }
}


sum = function (array) {
  var sum = []  
  for(i = 0; i < array[0].length; i++){
    var s= 0
    for (j = 0; j < array.length; j++) {
       s += array[j][i] 
    }
    sum.push(s)
  }
  return sum
}

sp = function (scalar, array) {
  var sum = []
  for(i = 0; i < array.length; i++){
   sum.push(scalar * array[i]);
  }
  return sum
}
abs = function (array) {
  var sum = 0
  for(i = 0; i < array.length; i++){
   sum += Math.pow(Math.abs(array[i]), 2)
  }
  return Math.sqrt(sum)
}


mathLib.odeMethod = function (method, func, N, t, h, params, pop, birthrate, tol = 1e4) {
  let tempArray
  let k1, k2, k3, k4, k5, k6, y, z, s
  let a, b, b2, c, d, out
  
  switch (method) {
  case 'euler':
    tempArray = func(t, N, params, pop, birthrate)
    return [sum ([N,sp(h,tempArray)]), h]
    
  case 'rk4':
    c = [0, 1/3, 2/3, 1]
    a21 = 1/3 ,a31 = -1/3 ,a32 = 1, a41 = 1 ,a42 = -1 , a43 = 1
    b = [0, 1/8, 3/8, 3/8, 1/8]   
    k1 = func(t          , N, params, pop, birthrate)
    k2 = func(t + c[2] * h , sum([N , sp(h * a21 , k1)]), params, pop, birthrate)
    k3 = func(t + c[3] * h , sum([N , sp(h * a31 , k1), sp(h * a32, k2)]), params, pop, birthrate)
    k4 = func(t + c[4] * h , sum([N , sp(h * a41, k1),  sp(h *a42, k2), sp(h *a43, k3)]), params, pop, birthrate)
    return [sum ([N, sp (h *  b[1] , sum ([k2 , k3])) ,sp(h * b[2] ,sum ([k1 , k2]))]), h]
    
  case 'rkf45':
    c = [0,0, 1/4, 3/8, 12/13, 1, 1/2]
    a = [[0, 0, 0, 0, 0],
             [0,0, 0, 0, 0, 0],
             [0,1/4, 0, 0, 0, 0],
             [0,3/32, 9/32, 0, 0, 0],
             [0,1932/2197, -7200/2197, 7296/2197, 0, 0],
             [0,439/216, -8, 3680/513, -845/4104, 0],
             [0,-8/27, 2, -3544/2565, 1859/4104, -11/40]]      
    b = [0,25/216, 0, 1408/2565, 2197/4104, -1/5, 0]
    b2 = [0,16/135,   0,  6656/12825,   28561/56430,  -9/50,  2/55]
    k1 = func(t          , N, params, pop, birthrate)
    k2 = func(t + c[2] * h , sum([N , sp(h * a[2][1] , k1)]), params, pop, birthrate)
    k3 = func(t + c[3] * h , sum([N , sp(h * a[3][1] , k1), sp(h * a[3][2], k2)]), params, pop, birthrate)
    k4 = func(t + c[4] * h , sum([N , sp(h * a[4][1], k1),  sp(h *a[4][2], k2), sp(h *a[4][3], k3)]), params, pop, birthrate)
    k5 = func(t + c[5] * h , sum([N , sp(h * a[5][1], k1), sp(h *a[5][2], k2), sp(h *a[5][3], k3), sp(h *a[5][4], k4)]), params, pop, birthrate)
    k6 = func(t + c[6] * h , sum([N , sp(h * a[6][1], k1), sp(h *a[6][2], k2), sp(h *a[6][3], k3), sp(h *a[6][4], k4), sp(h *a[6][5], k5)]), params, pop, birthrate)
    y = sum ([N, sp (h *  b[1], k1), sp (h * b[2], k2) ,sp(h * b[3], k3), sp (h *  b[4], k4), sp (h * b[5], k5) ,sp(h * b[6], k6)])
    z = sum ([N, sp (h *  b2[1], k1), sp (h * b2[2], k2) ,sp(h * b2[3], k3), sp (h *  b2[4], k4), sp (h * b2[5], k5) ,sp(h * b2[6], k6)])
    s = abs(sum([z,sp(-1,y)])) ?  Math.sqrt(Math.sqrt((tol * h)/ (2 * abs(sum([z,sp(-1,y)]))))) : h;
    // console.log(s)
    return[z, s] 
    
  case 'rk45dp7':
    c = [0, 0, 1/5, 3/10, 4/5, 8/9, 1, 1]
    a = [[0,0, 0, 0, 0, 0, 0],
            [0,0, 0, 0, 0, 0, 0],
            [0, 1/5, 0, 0, 0, 0, 0],
            [0, 3/40, 9/40, 0, 0, 0, 0],
            [0, 44/45, -56/15, 32/9, 0, 0, 0],
            [0, 19372/6561, -25360/2187, 64448/6561, -212/729, 0, 0],
            [0, 9017/3168, -355/33, 46732/5247, 49/176, -5103/18656, 0],
            [0, 35/384, 0, 500/1113, 125/192, -2187/6784, 11/84]]      
    b = [0, 5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40]
    b2 = [0, 35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0]
    d = [-12715105075/11282082432, 0, 87487479700/32700410799, -10690763975/1880347072, 701980252875/199316789632, -1453857185/822651844, 69997945/29380423]  
    k1 = func(t          , N, params, pop, birthrate)
    k2 = func(t + c[2] * h , sum([N , sp(h * a[2][1] , k1)]), params, pop, birthrate)
    k3 = func(t + c[3] * h , sum([N , sp(h * a[3][1] , k1), sp(h * a[3][2], k2)]), params, pop, birthrate)
    k4 = func(t + c[4] * h , sum([N , sp(h * a[4][1], k1),  sp(h *a[4][2], k2), sp(h *a[4][3], k3)]), params, pop, birthrate)
    k5 = func(t + c[5] * h , sum([N , sp(h * a[5][1], k1), sp(h *a[5][2], k2), sp(h *a[5][3], k3), sp(h *a[5][4], k4)]), params, pop, birthrate)
    k6 = func(t + c[6] * h , sum([N , sp(h * a[6][1], k1), sp(h *a[6][2], k2), sp(h *a[6][3], k3), sp(h *a[6][4], k4), sp(h *a[6][5], k5)]), params, pop, birthrate)
    y = sum ([N, sp (h *  b[1], k1), sp (h * b[2], k2) ,sp(h * b[3], k3), sp (h *  b[4], k4), sp (h * b[5], k5) ,sp(h * b[6], k6)])
    z = sum ([N, sp (h *  b2[1], k1), sp (h * b2[2], k2) ,sp(h * b2[3], k3), sp (h *  b2[4], k4), sp (h * b2[5], k5) ,sp(h * b2[6], k6)])
    s = Math.sqrt(Math.sqrt((tol * h)/ (2 * abs(sum([z,sp(-1,y)])))));
    // console.log(h, (2 * abs(sum([z,sp(-1,y)]))))
    return [z, s]
    
  }
  // return out  
}

module.exports = mathLib



