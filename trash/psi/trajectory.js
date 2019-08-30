fs = require('fs')
let fmin = require ('fmin')
var mathLib  = require('./mathLib')
var linearInterpolator = require('linear-interpolator/node_main')

var delT = 0.03832991102// = 2/52
var LondonBidata, LondonCovar
// params = [26.7306094868902,  0.273464921798168, 73.05, 0.043283234551231, 45.66, 0.513095687309477, 0.966454196355743, 0.008210807, 0.001633519, 0.990155663, 1.10E-08,1940,  1944]
// LondonData: time, cases
// LondonCovar: time, population(time), birthrate(time -4)
var param = fs.readFileSync('./paramSet_psi.csv').toString()
var paramss = []
var lines = param.split('\n')
for (let i = 1; i < lines.length; i++) {
  paramss.push(lines[i].split(','))
}
if(lines[0].split(',')[9] !== "R_0"){
  for (i = 0; i < paramss.length; i++) {
    var tem = paramss[i][9]
  paramss[i][9] = paramss[i][10]
  paramss[i][10] = tem
  }
}
var paramSet = [paramss][0]
var indx =  new Array(12)
indx[0] = 1; indx[1] = 1; 
indx[3] = 1; indx[5] = 1

//* 1st data set
var London_covar = fs.readFileSync('./London_covar.csv').toString()
var LondonCovar = []
var lines = London_covar.split('\n')
for (let i = 1; i < lines.length; i++) {
  LondonCovar.push(lines[i].split(','))
}
var dataCovar = [LondonCovar][0]
//* 2nd data set
LondonBidata = []
var London_BiData = fs.readFileSync('./London_BiData.csv').toString()
var lines = London_BiData.split('\n')
for (let i = 1; i < lines.length; i++) {
  LondonBidata.push(lines[i].split(','))
}
var dataCases = [LondonBidata][0]
//* main function****************************************************************

var d1 = []// read time and population from 1st data and make interpolation function
var d2 = []// read time and birthrate from 1st data and make interpolation function
for (let i = 0; i < dataCovar.length - 1; i++) {
  d1.push([Number(dataCovar[i][0]), Number(dataCovar[i][1])])
  d2.push([Number(dataCovar[i][0]), Number(dataCovar[i][2])])
}
var interpolPop = linearInterpolator(d1)
var interpolBirth = linearInterpolator(d2)

//* ODE function and ...
var res = []
for (let id = 0; id < paramSet.length - 1; id++) {
  var params = []
  for (let j = 0; j < paramSet[0].length; j++) {
    params.push(Number(paramSet[id][j]))
  }
var [R0, amplitude, gamma, mu, sigma, rho, psi, S_0, E_0, R_0, I_0, t0, t1] = params
var estim = []
var place = []
for (let i = 0; i < params.length - 1; i++) {
  if (indx[i] == 1) {
    place.push(i)
    if ((i === 1) || (i === 5)) {
      estim.push(Math.log(params[i] / (1 - params[i])))
    } else {
      estim.push(Math.log(params[i]))
    }
    
  }
}
function poly (params, t, N) {
  var va
  if (t < 1968)
    va = 0
  else if (t >= 1968 && t <= 1969)
    va = 0.33
  else if (t >= 1969 && t <= 1970)
    va = 0.46
  else if (t >= 1970 && t <= 1971)
    va = 0.51
  else if (t >= 1971 && t <= 1972)
    va = 0.53
  else if (t >= 1972 && t <= 1973)
    va = 0.52
  else if (t >= 1973 && t <= 1974)
    va = 0.46
  else if (t >= 1974 && t <= 1975)
    va = 0.46
  else if (t >= 1975 && t <= 1976)
    va = 0.48
  else if (t >= 1976 && t <= 1977)
    va = 0.48
  else if (t >= 1977 && t <= 1978)
    va = 0.51
  else if (t >= 1978 && t <= 1979)
    va = 0.53;
  else if (t >= 1979 && t <= 1980)
    va = 0.55;
  else if (t >= 1980 && t <= 1981)
    va = 0.58;
  else if (t >= 1981 && t <= 1982)
    va = 0.60
  else if (t >= 1982 && t <= 1983)
    va = 0.63
  else if (t >= 1983 && t <= 1984)
    va = 0.68
  else if (t >= 1984 && t <= 1985)
    va = 0.71
  else if (t >= 1985 && t <= 1988)
    va = 0.76
  else if (t >= 1988 && t <= 1989)
    va = 0.814
  else if (t >= 1989 && t <= 1990)
    va = 0.9488
  else if (t >= 1990 && t <= 1991)
    va = 0.9818
  else if (t >= 1991 && t <= 1992)
    va = 0.90
  else if (t >= 1992 && t <= 1993)
    va = 0.92
  else if (t >= 1993 && t <= 1994)
    va = 0.91
  else if (t >= 1994 && t <= 1995)
    va = 0.91
  else if (t >= 1995 && t <= 1996)
    va = 0.92
  else if (t >= 1996 && t <= 1997)
    va = 0.92
  else if (t >= 1997 && t <= 1998)
    va = 0.91
  else if (t >= 1998 && t <= 1999)
    va = 0.88
  else if (t >= 1999 && t <= 2000)
    va = 0.88
  else if (t >= 2000 && t <= 2001)
    va = 0.87
  else if (t >= 2001 && t <= 2002)
    va = 0.84
  else if (t >= 2002 && t <= 2003)
    va = 0.82
  else if (t >= 2003 && t <= 2004)
    va = 0.80
  else if (t >= 2004 && t <= 2005)
    va = 0.81
  else if (t >= 2005 && t <= 2006)
    va = 0.84
  else if (t >= 2006 && t <= 2007)
    va = 0.85
  else if (t >= 2007 && t <= 2008)
    va = 0.85
  else if (t >= 2008 && t <= 2009)
    va = 0.85
  else if (t >= 2009 && t <= 2010)
    va = 0.88
  else
    va = 0.89
  var seas, dy = []
  var R0 = params[0], amplitude = params[1], gamma = params[2], mu = params[3], sigma = params[4] 
  var beta0 = R0 * (gamma + mu) * (sigma + mu) / sigma
  var S = N[0], E = N[1], R = N[2],I = N[3]
  var pop = interpolPop(t)
  var birthrate = interpolBirth(t)

  var tt = (t - Math.floor(t)) * 365.25
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
function EulersMethod (params, delT) {
  var rho = params[5], psi = params[6], t0 = 1940, tdata = 1944
  var steps = 1000, arr2, arr = [], pop = interpolPop(t0),
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
      if (k <= tdata - delT ) {
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
    if ((place[i] === 1) || (place[i] === 5)) {
      params[place[i]] = 1 / (1 + Math.exp(-estim[i]))
    } else {
      params[place[i]] = Math.exp(estim[i])
    }
  }
  var rho = params[5], psi = params[6]
  var simCases = EulersMethod(params, delT)
  for (let i = 0; i < simCases.length ; i++) {
    var mn = rho * simCases[i]
    var v = mn * (1.0 - rho + psi * psi * mn)
    var tol = 1.0e-18
    var modelCases = Number(dataCases[i][1])
    if(!isNaN(modelCases)){
      if (modelCases > 0.0) {
        lik = mathLib.pnorm(modelCases + 0.5, mn, Math.sqrt(v) + tol, 1, 0) - mathLib.pnorm(modelCases - 0.5, mn, Math.sqrt(v) + tol, 1, 0) + tol
      } else {
        lik = mathLib.pnorm((modelCases + 0.5, mn, Math.sqrt(v) + tol)) + tol
      }
    } else {
    lik = 1
    }
    loglik = loglik + Math.log(lik)
  }
  // console.log(loglik, estim)
  return [-(loglik).toFixed(6)]
}
//* Optimizer
// logLik(estim)
var idx = 0
var solution = fmin.nelderMead(logLik, estim)
// console.log(solution.x)
for (let j = 0;j < params.length; j++) {
  if (indx[j] === 1) {
    if ((j === 1) || (j === 5)){
      params[j] = 1/ (1 + Math.exp(-solution.x[idx]))
    } else {
      params[j] = Math.exp(solution.x[idx])
    }
    idx++
  }
}
  console.log(params,- solution.fx)
  res.push([...params,- solution.fx])
}

const createCsvWriter = require('csv-writer').createArrayCsvWriter;
const csvWriter = createCsvWriter({
  header: ['R0', 'amplitude', 'gamma', 'mu', 'sigma', 'rho', 'psi', 'S_0', 'E_0', 'R_0', 'I_0', 'LogLik'],
  path: './DeterministicSEIR_run1.csv'
})
 
 
csvWriter.writeRecords(res)
  .then(() => {
  console.log('...Done')
})