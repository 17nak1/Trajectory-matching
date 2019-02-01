let mathLib = require('./mathLib')
let snippet = require('./modelSnippet.js')
let linearInterpolator = require('linear-interpolator/node_main')
let fmin = require('fmin')

//* main function
//* dataCovar and dataCases created in "creatDataSet.R " 
//* params all parameters and initials
//* times = [t0, data start time]
//* indx = an array with value 1 at the place of the parameters who are going to estimated
function traj_match (dataCovar, dataCases, params, times, indx) {
  var estim = []
  var place = []
  var states = []
  var delT = 0.03832991102//Biweekly =  2/52
  for (let i = 0; i < params.length; i++) {
    params[i] = Number(params[i])
    if (indx[i] === 1) {//***change the initial values of estimating parameters to the log scale
      place.push(i)
      if ((i === 1) || (i === 5)) {
        estim.push(Math.log(params[i] / (1 - params[i])))
      } else {
        estim.push(Math.log(params[i]))
      }
    }
  }
  var d1 = []// read time and population from dataCovar for the interpolation function
  var d2 = []// read time and birthrate from dataCovar for the interpolation function
  for (let i = 0; i < dataCovar.length - 1; i++) {
    d1.push([Number(dataCovar[i][0]), Number(dataCovar[i][1])])
    d2.push([Number(dataCovar[i][0]), Number(dataCovar[i][2])])
  }
  var interpolPop = linearInterpolator(d1)
  var interpolBirth = linearInterpolator(d2)

  //* ODE solver
  function EulersMethod (params, times, delT) {
    var rho = params[5], psi = params[6], t0 = times[0], tdata = times[1]
    var steps = 1000// total number of steps in the Euler method
    var arr2, arr = [], pop = interpolPop(t0)
    states = []
    var N = snippet.initz(interpolPop(t0), params[7], params[8], params[9], params[10])
    for (let k = t0; k <= Number(dataCases[dataCases.length - 2][0]) + delT / 3; k += delT) {
      N[4] = 0
      if (k <= tdata && k > tdata - delT) {
        k = tdata
      }
      for (let stp = 0; stp < steps; stp++) { // steps in each time interval
        var pop = interpolPop(k + stp / steps * delT)
        var birthrate = interpolBirth(k + stp / steps * delT)
        arr2 = snippet.skeleton(params, k + stp / steps * delT, N, pop, birthrate)
        N = N.map((a, i) => a + arr2[i] * 1 / steps * delT)
      }
      H = Math.round(N[4])
      var simCases = snippet.rmeasure(N[4], rho, psi)
      if (k > tdata - 2 * delT) {
        if (k <= tdata - delT) {
          k = tdata - delT
        }
        arr.push(H)
        states.push([k , N[0], N[1], N[3], N[2], H])// generate the part we have as traj@states
      }
    }
    return arr
  }
  //* calculate log likelihood
  function logLik (estim) {
    var loglik = 0
    for (let i = 0; i < estim.length; i++) {
      if ((place[i] === 1) || (place[i] === 5)) {// change to the exp scale and let optimizer to search all real numbers
        params[place[i]] = 1 / (1 + Math.exp(-estim[i]))
      } else {
        params[place[i]] = Math.exp(estim[i])
      }
    }
    var rho = params[5], psi = params[6]
    var simH = EulersMethod(params, times, delT)
    for (let i = 0; i < simH.length; i++) {
      var likld = snippet.dmeasure(rho, psi, simH[i], dataCases[i][1])
      loglik = loglik + Math.log(likld)
    }
    return [-(loglik).toFixed(6)]
  }
  //* Optimizer function using Nelder Mead method
  var idx = 0
  var solution = fmin.nelderMead(logLik, estim)
  for (let j = 0;j < params.length; j++) {
    if (indx[j] === 1) {//*** now use exp() to get back to the regular scale
      if ((j === 1) || (j === 5)){
        params[j] = 1/ (1 + Math.exp(-solution.x[idx]))
      } else {
        params[j] = Math.exp(solution.x[idx])
      }
      idx++
    }
  }
  activeDownload()
  // const createCsvWriter = require('csv-writer').createArrayCsvWriter//download the states values
  // const csvWriter = createCsvWriter({
  //   header: ['time', 'S', 'E', 'I', 'R', 'H'],
  //   path: './trajStates.csv'
  // }) 
  // csvWriter.writeRecords(states)
  //   .then(() => {
  // })
  return[params, -solution.fx]
}

module.exports = {
  traj_match
}
