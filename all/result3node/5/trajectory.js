/**
 *  @file       traj.js        This function attempts to match trajectories of a model's deterministic skeleton to data.
 *                             Trajectory matching is equivalent to maximum likelihood estimatedation under the assumption 
 *                             that process noise is entirely absent, i.e., that all stochasticity is measurement error.
 *                             Accordingly, this method uses only the skeleton and dmeasure components of a POMP model.
 *
 *  @author     Nazila Akhavan
 *  @date       Jan 2019
 */

let snippet = require('./modelSnippet.js')
let mathLib = require('./mathLib')
let fmin    = require('fmin')
let fs = require('fs')

var LondonBidata, LondonCovar
var param = fs.readFileSync('./paramSet_psi.csv').toString()
var params = []
var lines = param.split('\n')
for (let i = 1; i < lines.length - 1; i++) {
  params.push(lines[i].split(','))
}
let date = new Date()
var index =  new Array(12)
index[0] = 1; index[1] = 1; 
index[3] = 1; index[5] = 1

//* 1st data set
var London_covar = fs.readFileSync('./London_covar.csv').toString()
var dataCovar = []
var lines = London_covar.split('\n')
for (let i = 1; i < lines.length; i++) {
  dataCovar.push(lines[i].split(','))
}
//* 2nd data set
dataCases = []
var London_BiData = fs.readFileSync('./London_BiData.csv').toString()
var lines = London_BiData.split('\n')
for (let i = 1; i < lines.length; i++) {
  dataCases.push(lines[i].split(','))
}
let times = [1940, Number(dataCases[0][0]), Number(dataCases[dataCases.length - 2][0])]
var d1 = [] // read time and population from 1st data and make interpolation function
var d2 = [] // read time and birthrate from 1st data and make interpolation function
for (let i = 0; i < dataCovar.length - 1; i++) {
  d1.push([Number(dataCovar[i][0]), Number(dataCovar[i][1])])
  d2.push([Number(dataCovar[i][0]), Number(dataCovar[i][2])])
}
var interpolPopulation = mathLib.interpolator(d1)
var interpolBirth = mathLib.interpolator(d2)
var res = []
let len = params[0].length
for (let j = 0; j < params.length; j++) {
  for (let i = 0; i < params[j].length; i++) {
      params[j][i] = Number(params[j][i])
    }
   
  st = params[j][len - 4] + params[j][len - 3] + params[j][len - 2] + params[j][len - 1];
  for (let i = 4; i > 0; i--) {
    params[j][len - i] = params[j][len - i] /st
  }
}

function traj_match (interpolPopulation, interpolBirth, dataCases, params, times, index) {
  let deltaT = 0.03832991102
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

    var simH = integrate(interpolPopulation, interpolBirth, params, times, deltaT)
    for (let i = 0; i < simH.length; i++) {
      likvalue = snippet.dmeasure(rho, psi, simH[i], dataCases[i][1])
      loglik = loglik + Math.log(likvalue)
    }
    return [-(loglik).toFixed(6)]
  }
  console.log([...params, -solution.fx])
  return[...params, -solution.fx]
}

//* ODE solver
function integrate (interpolPopulation, interpolBirth, params, times, deltaT) {
  var steps = 200 // Total number of steps in the Euler method
  var t0 = times[0]
  var dataStartTime = times[1]
  var dataEndTime = times[2]
  var rho = params[5]
  var psi = params[6]
  // var states = []
  var arr = []
  var pop 
  var birthrate
  var timetemp
  var Npre

  var N = snippet.initz(interpolPopulation(t0), params[7], params[8], params[9], params[10])
  var  k = t0 , count
  var flag = 0
  dt = deltaT
  while ( flag === 0 ) {
    Npre = N
    for (let stp = 0; stp < steps; stp++) { 
      pop = interpolPopulation(k + stp / steps * dt)
      birthrate = interpolBirth(k + stp / steps * dt)
      N45 = mathLib.odeMethod('rkf45', snippet.skeleton, N, k + stp / steps * dt, 1 / steps * dt, params, pop, birthrate)
      N = N45[0]
    }
    timetemp = k
    k += dt
    if (k > dataStartTime) {
      k = timetemp
      dt = dataStartTime - timetemp 
      N = Npre
    }
    if (k >= dataStartTime) {  
      k = timetemp + dt
      flag = 1
      arr.push(N[4])
      // states.push([k , N[0], N[1], N[3], N[2], N[4]])
    }
  }

  count = 0
  while (k < dataEndTime) {
    
    if (Number(dataCases[count + 1][0]) !== "undefined") {
      dt = Number(dataCases[count + 1][0]) - Number(dataCases[count][0])
    } else {
      dt = deltaT
    }
    
    N[4] = 0
    for (let stp = 0; stp < steps; stp++) { 
      pop = interpolPopulation(k + stp / steps * dt)
      birthrate = interpolBirth(k + stp / steps * dt)
      N45 = mathLib.odeMethod('rkf45', snippet.skeleton, N, k + stp / steps * dt, 1 / steps * dt, params, pop, birthrate)
      N = N45[0]
      H = N[4]
    }
    k += dt
    count++
    arr.push(H)
    // states.push([k , N[0], N[1], N[3], N[2], H])
  }
  return arr
}
for ( let id = 0; id < params.length; id++) {
  
  res.push(traj_match (interpolPopulation, interpolBirth, dataCases, params[id], times, index))
}

const createCsvWriter = require('csv-writer').createArrayCsvWriter;
  const csvWriter = createCsvWriter({
    header: ['R0', 'amplitude', 'gamma', 'mu', 'sigma', 'rho', 'psi', 'S_0', 'E_0', 'R_0', 'I_0', 'LogLik'],
    path: './resPsi.csv'
  })   
  csvWriter.writeRecords(res)
    .then(() => {
    console.log('...Done')
  })