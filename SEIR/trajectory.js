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
const paramObject = {
  R0Index : 0,
  AMPLITUDE : 1,
  GAMMA : 2,
  MU : 3,
  SIGMA : 4,
  RHO : 5,
  PSI : 6
}

var LondonBidata, LondonCovar
var param = fs.readFileSync('./ParamSet_DeterministicSEIR_run1.csv').toString()
var params = []
var lines = param.split('\n')
for (let i = 1; i < lines.length - 1; i++) {
  params.push(lines[i].split(','))
}
let date = new Date()
var index =  new Array(12)
index[0] = 1; 
index[1] = 1; 
index[3] = 1; index[5] = 1; index[6] = 1

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
  // From those amplitude and rho are in logit scale and the rest are in log scale
  for (let i = 0; i < params.length; i++) {
    params[i] = Number(params[i])
    if (index[i] === 1 ) {
      place.push(i)
      if ((i === paramObject.AMPLITUDE) || (i === paramObject.RHO)) {
        estimated.push(Math.log(params[i] / (1 - params[i]))) //logit scale
      } else {
        estimated.push(Math.log(params[i])) //log scale
      }
    }
  }
  // estimated = snippet.toEst(params)

  //* Optimizer function using Nelder Mead method
  solution = fmin.nelderMead(logLik, estimated)
  
  //* calculate log likelihood
  function logLik (estimated) {
    var likvalue = 0
    var loglik = 0
    var rho 
    var psi 
    for (let i = 0; i < estimated.length; i++) {
      if ((place[i] === paramObject.AMPLITUDE) || (place[i] === paramObject.RHO)) { //Change to the exp scale and let optimizer to search all real numbers.
        params[place[i]] = 1 / (1 + Math.exp(-estimated[i]))
      } else {
        params[place[i]] = Math.exp(estimated[i])
      }
    }
    rho = params[5]
    psi = params[6]

    var simH = integrate(interpolPopulation, interpolBirth, params, times, deltaT)
    for (let i = 0; i < simH.length; i++) {
      likvalue = snippet.dmeasure(rho, psi, simH[i], dataCases[i][1], 1)
      loglik = loglik + likvalue
    }
    ;console.log(params, loglik)
    return [-(loglik).toFixed(6)]
  }
  return[params, -solution.fx]
}

//* ODE solver
function integrate (interpolPopulation, interpolBirth, params, times, deltaT) {
  var steps = 200 // Total number of steps in the each interval.
  var t0 = times[0]
  var dataStartTime = times[1]
  var dataEndTime = times[2]
  var rho = params[5]
  var psi = params[6]
  var arr = []
  var pop 
  var birthrate
  var timetemp
  var Npre

  var N = snippet.initz(interpolPopulation(t0), params[7], params[8], params[9], params[10])
  var k = t0 , count
  var flag = 0
  dt = deltaT
  while ( flag === 0 ) {
    Npre = N
    for (let stp = 0; stp < steps; stp++) { 
      pop = interpolPopulation(k + stp / steps * dt)
      birthrate = interpolBirth(k + stp / steps * dt)
      N = mathLib.odeMethod('rkf45', snippet.skeleton, N, k + stp / steps * dt, 1 / steps * dt, params, pop, birthrate)
    }
    timetemp = k
    k += dt
    if (k > dataStartTime) {
      k = timetemp;
      dt = dataStartTime - timetemp ;
      N = Npre
    }
    if (k >= dataStartTime) {  
      k = timetemp + dt
      flag = 1
      arr.push(N[4])
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
      N = mathLib.odeMethod('rkf45', snippet.skeleton, N, k + stp / steps * dt, 1 / steps * dt, params, pop, birthrate)
      H = N[4]
    }
    k += dt
    count++
    arr.push(H)
  }
  return arr
}

for ( let id = 0; id < params.length; id++) {
  
  res.push(traj_match (interpolPopulation, interpolBirth, dataCases, params[id], times, index))
}

const createCsvWriter = require('csv-writer').createArrayCsvWriter;
  const csvWriter = createCsvWriter({
    header: ['R0', 'amplitude', 'gamma', 'mu', 'sigma', 'rho', 'psi', 'S_0', 'E_0', 'R_0', 'I_0', 'LogLik'],
    path: './sobolres1.csv'
  })   
  csvWriter.writeRecords(res)
    .then(() => {
    console.log('...Done')
  })