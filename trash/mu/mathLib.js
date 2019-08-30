
var mathLib = {}
let exp = 2.718281828
let pi = 3.141592654
var erf = require('math-erf')
var seedrandom = require('seedrandom')
var rng = seedrandom('915909831')

//* Set the seed for rnorm-In R:RNGkind("Mersenne-Twister",normal.kind="Ahrens-Dieter");set.seed(1234) 
//or normal.kind="Box-Muller" or "BuggyKindermanRamage" or "Inversion" or "KindermanRamage"
const libR = require('lib-r-math.js')
var {
  Normal,
  rng: {
    LecuyerCMRG,
    normal: { BoxMuller }
  }
} = libR
const ad = new LecuyerCMRG(1234)
const { rnorm } = Normal(new BoxMuller(ad))
// console.log(rnorm(10))

//* Set the seed for rgamma-In R:RNGkind("L'Ecuyer-CMRG", normal.kind="Box-Muller"); set.seed(1234)
const {
  Gamma
  // rng: {
  //   LecuyerCMRG,
  //   normal: { BoxMuller }
  // }
} = libR
const lc = new LecuyerCMRG(1234)
const { rgamma } = Gamma(new BoxMuller(lc))
// console.log(rgamma(1, 1, 0.1))

//* Set the seed for rbinom-In R: RNGkind("Knuth-TAOCP-2002");set.seed(1234)
const {
  Binomial,
  rng: { KnuthTAOCP2002 }
} = libR
const kn = new KnuthTAOCP2002(1234)
const { rbinom } = Binomial(kn)
// console.log(rbinom(2,40,.5))
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

mathLib.factorial = function (intValue) {
  var i, nextNumber, carret, result
  if (intValue === 0) {
    return '1'
  }
  if (!intValue) {
    return ''
  }
  result = intValue.toString().split('').reverse().map(Number)
  while (--intValue) {
    i = carret = 0
    while ((nextNumber = result[i++]) !== undefined || carret) {
      carret = (nextNumber || 0) * intValue + carret
      result[i - 1] = carret % 10
      carret = parseInt(carret / 10)
    }
  }
  return result.reverse().join('')
}

mathLib.reulermultinom = function (m = 1, size, rateAdd, dt, transAdd, rate, trans) {
  // console.log(size)
  var p = 0
  var j, k
  if ((size < 0) || (dt < 0) || (Math.floor(size + 0.5) !== size)) {
    for (k = 0; k < m; k++) trans[k + transAdd] = NaN
    return 0
  }
  for (k = 0; k < m; k++) {
    if (rate[k + rateAdd] < 0.0) {
      for (j = 0; j < m; j++) trans[j + transAdd] = NaN
      return 0
    }
    p += rate[k + rateAdd]// total event rate
  }
  if (p > 0) {
    size = rbinom(1, size, 1 - Math.exp(-p * dt)) // total number of events
    if (!(isFinite(size)))
      throw 'result of binomial draw is not finite.'
    m -= 1
    for (k = 0; k < m; k++) {
      if (rate[k + rateAdd] > p) p = rate[k + rateAdd]
      trans[k + transAdd] = ((size > 0) && (p > 0)) ? rbinom(1, size, rate[k + rateAdd] / p) : 0
      if (!(isFinite(size) && isFinite(p) && isFinite(rate[k + rateAdd]) && isFinite(trans[k + transAdd]))) {
        throw 'result of binomial draw is not finite.'
      }
      size -= trans[k + transAdd]
      p -= rate[k + rateAdd]
    }
    trans[m + transAdd] = size
  } else {
    for (k = 0; k < m; k++) trans[k + transAdd] = 0
  }
}

mathLib.logit = function (p) {
  return Math.log(p / (1 - p))
}
mathLib.expit = function (x) {
  return 1 / (1 + Math.exp(-x))
}

mathLib.rgammawn = function (sigma, dt) {
  var sigmasq
  sigmasq = Math.pow(sigma, 2)
  return (sigmasq > 0) ? rgamma(1, dt / sigmasq, sigmasq) : dt
}

mathLib.fromLogBarycentric = function (xN, n) {
  var sum = 0, xt = []
  for (let i = 0; i < n; i++) {
    xt.push(Math.exp(xN[i]))
    sum += xt[i]
  }
  for (let i = 0; i < n; i++) {
    xN[i] = xt[i] / sum
  }
  return xN
}

mathLib.toLogBarycentric = function (xN, n) {
  var sum = 0
  for (let i = 0; i < n; i++) {
    sum += xN[i]
  }
  for (let i = 0; i < n; i++) {
    xN[i] = Math.log(xN[i] / sum)
  }
  return xN
}

mathLib.rnorm = function (mu = 0, sd = 1) {
  var val = Math.sqrt(-2.0 * Math.log(rng())) * Math.cos(2.0 * pi * rng())
  return val * sd + mu
}
mathLib.matrix = function (Nrows) {
  var ary = []
  for (var i = 0; i < Nrows; i++) {
    ary[i] = []
  }
  return ary
}
module.exports = mathLib;

