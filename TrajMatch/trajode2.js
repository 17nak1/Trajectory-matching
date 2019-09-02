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
let ode2 = require('./ode2.js')

var LondonBidata, LondonCovar
var param = fs.readFileSync('./ParamSet_DeterministicSEIR_run26.csv').toString()
var params = []
var lines = param.split('\n')
for (let i = 1; i < lines.length; i++) {
  params.push(lines[i].split(','))
}
var index =  new Array(12)
// indx[1] = 1; 
// indx[3] = 1; indx[5] = 1; indx[6] = 1

//* 1st data set
var London_covar = fs.readFileSync('./London_covar.csv').toString()
var dataCovar = []
var lines = London_covar.split('\n')
for (let i = 1; i < lines.length ; i++) {
  dataCovar.push(lines[i].split(','))
}
//* 2nd data set
dataCases = []
var London_BiData = fs.readFileSync('./London_BiData.csv').toString()
var lines = London_BiData.split('\n')
for (let i = 1; i < lines.length - 1; i++) {
  dataCases.push(lines[i].split(','))
}
let times = []
for (let i = 0; i < dataCases.length; i++) {
  times.push(Number(dataCases[i][0]))//;console.log(times[0])
}

var d1 = []// read time and population from 1st data and make interpolation function
var d2 = []// read time and birthrate from 1st data and make interpolation function
for (let i = 0; i < dataCovar.length - 1; i++) {
  d1.push([Number(dataCovar[i][0]), Number(dataCovar[i][1])])
  d2.push([Number(dataCovar[i][0]), Number(dataCovar[i][2])])
}
var interpolPopulation = mathLib.interpolator(d1)
var interpolBirth = mathLib.interpolator(d2)
params = params[0]
for(i = 0; i < params.length; i++) {
  params[i] = Number(params[i])
}
var t0 = 1940, deltaT = 2/52
var rho = params[5] ,psi = params[6], arr = [] ,steps = 1000

var N = snippet.initz(interpolPopulation(t0), params[7], params[8], params[9], params[10])
for (let k = t0; k <= times[0]; k += deltaT) {//
  for (let stp = 0; stp < steps; stp++) { 
    var pop = interpolPopulation(k + stp / steps * deltaT)
    var birthrate = interpolBirth(k + stp / steps * deltaT)
    var N45 = mathLib.odeMethod('rkf45', snippet.skeleton, N, k + stp / steps * deltaT, 1 / steps * deltaT, params, pop, birthrate)
    N = N45[0]//;console.log(N)
  }
  arr.push([k, ...N])
}
N = [1.030473e+05, 1.223305e+02, 1.294772e+06, 8.176444e+01, 2.049497e+05 ]
// console.log(N)
var N45 = ode2.ode(times,  N, snippet.skeleton, params, interpolPopulation, interpolBirth) 
//       N = N45[0]
//       console.log(N)

// function traj_match (interpolPopulation, interpolBirth, dataCases, params, times, index) {
//   var deltaT = 0.03832991102
//   var tempIndex = 0
//   var estimated = []
//   var place = []
//   var states = []
//   var data1 = []
//   var data2 = []
//   var solution
//   var t0 = 1940
//    //*Change the initial values of estimating parameters(with index one) to the log or logit scale.
//   // From those amplitude and rho (parameters number 1 and 5) are in logit scale and the rest are in log scale
//   for (let i = 0; i < params.length; i++) {
//     params[i] = Number(params[i])
//     if (index[i] === 1) {
//       place.push(i)
//       if ((i === 1) || (i === 5)) {
//         estimated.push(Math.log(params[i] / (1 - params[i]))) //logit scale
//       } else {
//         estimated.push(Math.log(params[i])) //log scale
//       }
//     }
//   }

//     var likvalue = 0
//     var loglik = 0
//     var rho 
//     var psi 
//     var larr =[]
//     for (let i = 0; i < estimated.length; i++) {
//       if ((place[i] === 1) || (place[i] === 5)) { //Change to the exp scale and let optimizer to search all real numbers.
//         params[place[i]] = 1 / (1 + Math.exp(-estimated[i]))
//       } else {
//         params[place[i]] = Math.exp(estimated[i])
//       }
//     }
//     rho = params[5]
//     psi = params[6]
//     var simH = integrate(interpolPopulation, interpolBirth, params, times, t0, deltaT)
// }

// //* ODE solver
// function integrate (interpolPopulation, interpolBirth, params, times, t0, deltaT) {
//   var rho = params[5] ,psi = params[6], arr = [] ,steps = 1000
//   var N = snippet.initz(interpolPopulation(t0), params[7], params[8], params[9], params[10])
//   for (let k = t0; k <= times[0]; k += deltaT) {
//     for (let stp = 0; stp < steps; stp++) { 
//       var pop = interpolPopulation(k + stp / steps * deltaT)
//       var birthrate = interpolBirth(k + stp / steps * deltaT)
//       var N45 = mathLib.odeMethod('euler', snippet.skeleton, N, params, k + stp / steps * deltaT, 1 / steps * deltaT, pop, birthrate)
//       N = N45[0]
//     }
//     arr.push([k,...N])
//   }
  const createCsvWriter = require('csv-writer').createArrayCsvWriter;
  const csvWriter = createCsvWriter({
    header: ['times' ,'S', 'E', 'R', 'I', 'H'],
    path: './res.csv'
  }) 
  csvWriter.writeRecords(arr)
    .then(() => {
    console.log('...Done')
  })
//   return arr
// }

// console.log(traj_match (interpolPopulation, interpolBirth, dataCases, params, times, index) )
//    