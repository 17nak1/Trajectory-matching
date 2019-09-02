let snippet = require('./modelSnippet.js')
let mathLib = require('./mathLib')
let fmin    = require('fmin')
let fs = require('fs')


var LondonBidata, LondonCovar
var param = fs.readFileSync('./ParamSet_DeterministicSEIR_run28.csv').toString()
var params = []
var lines = param.split('\n')
for (let i = 1; i < lines.length - 1; i++) {
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
  times.push(Number(dataCases[i][0]))
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
for ( i = 0; i < params.length; i++) {
  params[i] = Number(params[i])
}
t0 = 1940, deltaT = 2/52
var rho = params[5] ,psi = params[6], arr = [] ,steps = 1000
// N = [1.030473e+05, 1.223305e+02, 1.294772e+06, 8.176444e+01, 2.049497e+05 ]
  N = [87531.4660, 560.6357, 1546060.3460, 402.3514, 265004.7898]

sum = function (array) {
  let sum = []  
  for(i = 0; i < array[0].length; i++){
    let s= 0
    for (j = 0; j < array.length; j++) {
       s += array[j][i] 
    }
    sum.push(s)
  }
  return sum
}

sp = function (scalar, array) {
  let sum = []
  for(i = 0; i < array.length; i++){
   sum.push(scalar * array[i]);
  }
  return sum
}

abs = function (array) {
  let sum = 0
  for(i = 0; i < array.length; i++){
   sum += Math.pow(Math.abs(array[i]), 2)
  }
  return Math.sqrt(sum)
}
maxerr = function(y0, y1, y2, atol, rtol, n) {
  var serr = 0, scal, delta;
  for (let i = 0; i < n; i++) {
    /* y2 is used to estimate next y-value */
    scal  = atol[i] + Math.max(Math.abs(y0[i]), Math.abs(y2[i])) * rtol[i]
    delta = Math.abs(y2[i] - y1[i])
    if (scal > 0) serr += Math.pow(delta/scal, 2);
  }
  return(Math.sqrt(serr/n)); /* Euclidean norm */
}




ode = function (times, y, f2, params, intPop, intBirth) {
  var t = times[0] , t1, R, hdefault = .1 * (times[1] - t)
  var k1 = [], k2 = [], k3 = [], k4 = [], k5 = [], k6 = []
  var res = [], y_rk4 = [], y_rk5 = [], difmax = 0, count, flag = 0, accept = 0, yaccept, h, tOld, qerr = 4
  let facmin = 0.2, fac = 0.9 , facmax = 10
  let beta = 0; // 1 / qerr
  let alpha = 1 / qerr - 0.75 * beta
  let errOld = 0
  let hmax = hdefault; hmin = 0
  let n = y.length
  let atol = new Array(n).fill(1e-16)
  let rtol = new Array(n).fill(1e-16)
  let err = 0
  let dtnew
  res.push([t, ...y])
  for (let i = 1; i < times.length ; i++) {//times.length - 1
    t1 = times[i]
    h= hdefault
    flag = 0
  
    while ( flag === 0) {
      count = 0
            
      while(t < t1 && count < 2) {
        if ( (t1 - t ) < 1e-10) {
          t = t1 - h
          break
        }
        count++
        k1 = sp(h , f2( t, y, params, intPop(t), intBirth(t)))
        k2 = sp(h , f2( t + h, sum([y,sp(1/4, k1)]), params, intPop(t + h/4), intBirth(t + h/4) ))
        k3 = sp(h , f2( t + 3*h/8, sum([y, sp(3/32,k1),sp(9/32, k2)]), params, intPop(t + 3*h/8), intBirth(t + 3*h/8) )) 
        k4 = sp(h , f2( t + 12*h/13, sum([y, sp(1932/2197,k1),sp(-7200/2197, k2),sp(7296/2197, k3)]), params, intPop(t + 12*h/13), intBirth(t + 12*h/13))) 
        k5 = sp(h , f2( t + h, sum([y, sp(439/216, k1),sp(-8,k2), sp(3680/513, k3), sp(-845/4104, k4)]), params, intPop(t + h), intBirth(t + h)))
        k6 = sp(h , f2( t + h/2 ,sum([y, sp(-8/27, k1),sp(2, k2), sp(-3544/2565, k3), sp(1859/4104,k4), sp(-11/40, k5)]), params, intPop(t + h/2), intBirth(t + h/2)))
        y_rk4 = sum([y, sp(25/216, k1), sp(1408/2565, k3), sp(2197/4104, k4), sp(-1/5, k5)])                  // compute RK4
        y_rk5 = sum([y, sp(16/135, k1), sp(6656/12825, k3) ,sp(28561/56430, k4), sp( -9/50, k5) ,sp( 2/55, k6)]);       // compute RK5
        
        if ( accept == 1) {
          // accept = 0
          break
        }
        
        err = maxerr(y, y_rk4, y_rk5, atol, rtol, y.length)
        dtnew = h;
        if (err == 0) {  /* use max scale if all tolerances are zero */
          dtnew  = Math.min(dt * 10, hmax);
          errold = Math.max(err, 1e-4); /* 1e-4 taken from Press et al. */
          accept = TRUE;
        } else if (err < 1) {//console.log("break", count)
          if (accept) {
            dtnew = Math.min(hmax, dt * Math.min(safe * Math.pow(err, -alpha) * Math.pow(errold, beta), maxscale));
          }
          errold = Math.max(err, 1e-4); /* 1e-4 taken from Press et al. */
          accept = TRUE;
          break
        } else if (err > 1) {
          nreject++;    /* count total number of rejected steps */
          accept = FALSE;
          dtnew = dt * Math.max(safe * Math.pow(err, -alpha), minscale);
        }
        h = dtnew
      } // while loopin
      // console.log("end", t, t1,h, count)
      yaccept = y
      y = y_rk5
      // console.log(y, y_rk5 == yaccept)
      // console.log(y)
      tOld = t
      t += h
      // console.log(t,t1)
      if (t > t1 ) {
        t = t - h
        h = t1 - t
        y= yaccept //;console.log("t>t1",t == t1)
        accept = 1
      } 
      if (t == t1) flag = 1;//console.log("t",t == t1)
    } // while loop

    res.push([t, ...y]); console.log( y)
    if ( t >= Number(dataCases[0][0])) {
      y[4] = 0
    }
  } // for loop
  return res
} // function

var arr = ode(times,  N, snippet.skeleton, params, interpolPopulation, interpolBirth) 

const createCsvWriter = require('csv-writer').createArrayCsvWriter;
  const csvWriter = createCsvWriter({
    header: ['time' ,'S', 'E', 'R', 'I', 'H'],
    path: './trajexam.csv'
  }) 
  csvWriter.writeRecords(arr)
    .then(() => {
    console.log('...Done')
  })