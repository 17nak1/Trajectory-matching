/**
 *  @file       traj.js        This function attempts to match trajectories of a model's deterministic skeleton to data.
 *                             Trajectory matching is equivalent to maximum likelihood estimatedation under the assumption 
 *                             that process noise is entirely absent, i.e., that all stochasticity is measurement error.
 *                             Accordingly, this method uses only the skeleton and dmeasure components of a POMP model.
 *
 *  @author     Nazila Akhavan
 *  @date       Jan 2019
 */

let linearInterpolator = require('linear-interpolator/node_main')
let snippet = require('./modelSnippet.js')
let mathLib = require('./mathLib')
let fmin    = require('fmin')

/**
  * interpolBirth : Linear interpolation for birth rates.
  * interpolPopulation   : Linear interpoplation for population.
  * dataCovar     : Matrix created in "creadataStartTimeSet.R". 
  * dataCases     : Matrix created in "creadataStartTimeSet.R". 
  * params        : Array of parameters and initial states i.e. [R0, amplitude, gamma, mu, sigma, rho, psi, S_0, E_0, I_0, R_0]
  * times         : Array of 2 values, t0 and start time in dataCases
  * index         : Array with value 1 at the place of the parameters who are going to be estimated.
  * place         : Array include the index of parameters who are going to be estimated.
  * estimated     : Array of initial value of the parameters who are going to be estimated.
  * deltaT        : It is considered biweekly (2/52).
  * states        : Empty array for calculated states.
  * solution      : Include the result of  the optimizer function using Nelder Mead method.
  * data1         : Read time and population from dataCovar for the interpolation function
  * data2         : Read time and birthrate from dataCovar for the interpolation function
 */
function traj_match (dataCovar, dataCases, params, times, index) {
  var deltaT = 0.03832991102
  var tempIndex = 0
  var estimated = []
  var place = []
  var states = []
  var data1 = []
  var data2 = []
  var solution
   //*Change the initial values of estimating parameters(with index one) to the log or logit scale.
  // From those amplitude and rho (parameters number 1 and 5) are in logit scale and the rest are in log scale
  for (let i = 0; i < params.length; i++) {
    params[i] = Number(params[i])
    if (index[i] === 1) {
      place.push(i)
      if ((i === 1) || (i === 5)) {
        estimated.push(Math.log(params[i] / (1 - params[i]))) //logit scale
      } else {
        estimated.push(Math.log(params[i])) //log scale
      }
    }
  }
  for (let i = 0; i < dataCovar.length - 1; i++) {
    data1.push([Number(dataCovar[i][0]), Number(dataCovar[i][1])])
    data2.push([Number(dataCovar[i][0]), Number(dataCovar[i][2])])
  }
  var interpolPopulation = linearInterpolator(data1)
  var interpolBirth = linearInterpolator(data2)

  //* Optimizer function using Nelder Mead method
  solution = fmin.nelderMead(logLik, estimated)
  for (let j = 0;j < params.length; j++) {
    if (index[j] === 1) { // Using exp and expit to get back to the regular scale.
      if ((j === 1) || (j === 5)){
        params[j] = 1/ (1 + Math.exp(-solution.x[tempIndex]))
      } else {
        params[j] = Math.exp(solution.x[tempIndex])
      }
      tempIndex++
    }
  }
  activateDownload()
  //* calculate log likelihood
  function logLik (estimated) {
    var likvalue = 0
    var loglik = 0
    var rho 
    var psi 
    for (let i = 0; i < estimated.length; i++) {
      if ((place[i] === 1) || (place[i] === 5)) { //Change to the exp scale and let optimizer to search all real numbers.
        params[place[i]] = 1 / (1 + Math.exp(-estimated[i]))
      } else {
        params[place[i]] = Math.exp(estimated[i])
      }
    }
    rho = params[5]
    psi = params[6]
    var simH = EulersMethod(interpolPopulation, interpolBirth, params, times, deltaT)
    for (let i = 0; i < simH.length; i++) {
      likvalue = snippet.dmeasure(rho, psi, simH[i], dataCases[i][1])
      loglik = loglik + Math.log(likvalue)
    }
    return [-(loglik).toFixed(6)]
  }
  
  return[params, -solution.fx]
}

//* ODE solver
function EulersMethod (interpolPopulation, interpolBirth, params, times, deltaT) {
  var steps = 1000 // Total number of steps in the Euler method
  var t0 = times[0]
  var dataStartTime = times[1]
  var datadataEndTime = times[2]
  var rho = params[5]
  var psi = params[6]
  var states = []
  var arr = []
  var pop 
  var birthrate
  var tempArray
  var N = snippet.initz(interpolPopulation(t0), params[7], params[8], params[9], params[10])
  for (let k = t0; k <= datadataEndTime + deltaT / 3; k += deltaT) {
    N[4] = 0
    if (k <= dataStartTime && k > dataStartTime - deltaT) {
      k = dataStartTime
    }
    for (let stp = 0; stp < steps; stp++) { 
      var pop = interpolPopulation(k + stp / steps * deltaT)
      var birthrate = interpolBirth(k + stp / steps * deltaT)
      tempArray = snippet.skeleton(params, k + stp / steps * deltaT, N, pop, birthrate)
      N = N.map((a, i) => a + tempArray[i] * 1 / steps * deltaT)
    }
    H = Math.round(N[4])
    var simCases = snippet.rmeasure(N[4], rho, psi)
    if (k > dataStartTime - 2 * deltaT) {
      if (k <= dataStartTime - deltaT) {
        k = dataStartTime - deltaT
      }
      arr.push(H)
      states.push([k , N[0], N[1], N[3], N[2], H]) // Generate the part we have as traj@states. 
    }
  }
  return arr
}
module.exports = {
  traj_match
}
