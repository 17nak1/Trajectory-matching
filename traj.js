var mathLib = require('./mathLib')
var linearInterpolator = require('linear-interpolator/node_main')
let fmin = require('fmin')

//* main function****************************************************************
function traj_match (dataCovar, dataCases, params, indx) {
  var estim = []
  var place = []
  var delT = 0.03832991102// = 2/52
  for (let i = 0; i < params.length; i++) {
    params[i] = Number(params[i])
    if (indx[i] === 1) {
      place.push(i)
      estim.push(Math.log(params[i]))
    }
  }
  var d1 = []// read time and population from 1st data and make interpolation function
  var d2 = []// read time and birthrate from 1st data and make interpolation function
  for (let i = 0; i < dataCovar.length - 1; i++) {
    d1.push([Number(dataCovar[i][0]), Number(dataCovar[i][1])])
    d2.push([Number(dataCovar[i][0]), Number(dataCovar[i][2])])
  }
  var interpolPop = linearInterpolator(d1)
  var interpolBirth = linearInterpolator(d2)

  //* ODE function and ...
  function poly (params, time, N) {
    var va = 0, seas, dy = []
    var R0 = params[0], amplitude = params[1], gamma = params[2], mu = params[3], sigma = params[4] 
    var beta0 = R0 * (gamma + mu) * (sigma + mu) / sigma
    var S = N[0], E = N[1], R = N[2], I = N[3]
    var pop = interpolPop(time)
    var birthrate = interpolBirth(time)
    var tt = (time - Math.floor(time)) * 365.25
    if ((tt >= 7 && tt <= 100) || (tt >= 115 && tt <= 199) || (tt >= 252 && tt <= 300) || (tt >= 308 && tt <= 356)) {
      seas = 1 + amplitude * 0.2411 / 0.7589
    } else {
      seas = 1 - amplitude
    }
    var Beta = beta0 * seas / pop
    dy[0] = birthrate * (1 - va) - Beta * S * I - mu * S
    dy[1] = Beta * S * I - (sigma + mu) * E
    dy[2] = gamma * I - mu * R + birthrate * va
    dy[3] = sigma * E - (gamma + mu) * I
    dy[4] = gamma * I
    return dy
  }

  // ODE solver
  function EulersMethod (params, covarData, delT) {
    var rho = params[5], psi = params[6], t0 = params[11], tdata = params[12]
    var steps = 800, arr2, arr = [], pop = interpolPop(t0),
      m = pop / (params[7] + params[8] + params[9] + params[10]),
      S = Math.round(m * params[7]),
      E = Math.round(m * params[8]),
      R = Math.round(m * params[9]),
      I = Math.round(m * params[10]),
      H = 0,
      N = [S, E, R, I, H]
    for (let k = t0; k <= Number(dataCases[dataCases.length - 2][0]) + delT / 3; k += delT) {
      N[4] = 0
      if (k <= tdata && k > tdata - delT) {
        k = tdata
      }
      for (let stp = 0; stp < steps; stp++) { // steps in each time interval
        arr2 = poly(params, k + stp / steps * delT, N)
        N = N.map((a, i) => a + arr2[i] * 1 / steps * delT)
      }
      H = Math.round(N[4])
      var mn = rho * H
      var v = mn * (1.0 - rho + psi * psi * mn)
      var tol = 1.0e-18
      var cases = mathLib.rnorm(mn, Math.sqrt(v) + tol)
      if (cases > 0) {
        cases = Math.round(cases)
      } else {
        cases = 0
      }
      if (k > tdata - 2 * delT) {
        if (k <= tdata - delT) {
          k = tdata - delT
        }
        arr.push(H)
      // arr.push([k + delT , N[0], N[1], N[2], N[3], H])
      }
    }
    return arr
  }
  //* *********************************************************************
  function logLik (estim) {
    var lik, loglik = 0
    for (let i = 0; i < estim.length; i++) {
      params[place[i]] = Math.exp(estim[i])
    }
    var rho = params[5], psi = params[6]
    var simCases = EulersMethod(params, dataCovar, delT)
    for (let i = 0; i < simCases.length; i++) {
      var mn = rho * simCases[i]
      var v = mn * (1.0 - rho + psi * psi * mn)
      var tol = 1.0e-18
      var modelCases = Number(dataCases[i][1])
      if (modelCases > 0.0) {
        lik = mathLib.pnorm(modelCases + 0.5, mn, Math.sqrt(v) + tol, 1, 0) - mathLib.pnorm(modelCases - 0.5, mn, Math.sqrt(v) + tol, 1, 0) + tol
      } else {
        lik = mathLib.pnorm((modelCases + 0.5, mn, Math.sqrt(v) + tol)) + tol
      }
      loglik = loglik + Math.log(lik)
    }
    // console.log(loglik, estim)
    return [-(loglik).toFixed(6)]
  }
  //* Optimizer
  var idx = 0
  var solution = fmin.nelderMead(logLik, estim)
  for (let j = 0; j < params.length; j++) {
    if (indx[j] === 1) {
      params[j] = Math.exp(solution.x[idx])
      idx++
    }
  }
  activeDownload()
  // params.splice(11, 2)
  return[params, -solution.fx]
}

module.exports = {
  traj_match
}
