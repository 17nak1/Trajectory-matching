var mathLib = {}
let pi = 3.141592654
var erf = require('math-erf')
var seedrandom = require('seedrandom')
var rng = seedrandom('915909831')
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
  var val = Math.sqrt(-2.0 * Math.log(rng())) * Math.cos(2.0 * pi * rng())
  return val * sd + mu
}

// mathLib.rnorm = function (mu = 0, sd = 1) {
//   if (sd < 0) {
//     return NaN
//   }
//   let arr, val, a = 1 / Math.sqrt(2 * pi )
//   val = rng()
//   arr = a * Math.exp(- Math.pow(val,2)/ 2 )
//   return (arr * sd + mu)  // Why the fraction makes cases similar?
// }
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

module.exports = mathLib
