const sum = function (array) {
  var sum = []  
  for(i = 0; i < array[0].length; i++){
    var s= 0
    for (j = 0; j < array.length; j++) {
       s += array[j][i] 
    }
    sum.push(s)
  }
  return sum
}

const sp = function (scalar, array) {
  var sum = []
  for(i = 0; i < array.length; i++){
   sum.push(scalar * array[i]);
  }
  return sum
}
const abs = function (array) {
  var sum = 0
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

// initial values
let  y0 = [1, 0, 0.9], istep
const f2 = function(t, y) {
  var dy1, dy2, dy3
  dy1 = -2 * y[1] * y[2]
  dy2 = 1.25* y[0] * y[2]
  dy3 = -0.5* y[0] * y[1]
  return [dy1, dy2, dy3]
}

times = []
for (let i = 0; i <= 20; i = istep + .1) {
  istep = Number(i.toFixed(8))
  times.push(istep)
}

const ode = function (times, y, f2) {
  var t = times[0] , t1, R, hdefault = (times[1] - t)
  var k1 = [], k2 = [], k3 = [], k4 = [], k5 = [], k6 = []
  var res = [], y_rk4 = [], y_rk5 = [] , epsilon = 1e-11, difmax = 0, count, flag = 0, accept = 0, yaccept
  let facmin = 0.2, fac = 0.9 , facmax = 10
  let n = y.length
  let atol = new Array(n).fill(1e-16)
  let rtol = new Array(n).fill(1e-16)
  let err = 0
  for (let i = 1; i < times.length - 1 ; i++) {//times.length - 1
    t1 = times[i]
    h= hdefault
    while ( flag === 0) {
      count = 0
      // console.log(t, t1)
      while(t < t1 && count < 10) {
        count++
        k1 = sp(h , f2( t, y))
        k2 = sp(h , f2( t + h/4, sum([y,sp(1/4, k1)])))
        k3 = sp(h , f2( t + 3*h/8, sum([y, sp(3/32,k1),sp(9/32, k2)]))) 
        k4 = sp(h , f2( t + 12*h/13, sum([y, sp(1932/2197,k1),sp(-7200/2197, k2),sp(7296/2197, k3)]))) 
        k5 = sp(h , f2( t + h, sum([y, sp(439/216, k1),sp(-8,k2), sp(3680/513, k3), sp(-845/4104, k4)])))
        k6 = sp(h , f2( t + h/2 ,sum([y, sp(-8/27, k1),sp(2, k2), sp(-3544/2565, k3), sp(1859/4104,k4), sp(-11/40, k5)])))
        y_rk4 = sum([y, sp(25/216, k1), sp(1408/2565, k3), sp(2197/4104, k4), sp(-1/5, k5)])                       // compute RK4
        y_rk5 = sum([y, sp(16/135, k1), sp(6656/12825, k3) ,sp(28561/56430, k4), sp( -9/50, k5) ,sp( 2/55, k6)]);  // compute RK5
        if ( accept == 1) {
          accept = 0
          break
        }
        err = maxerr(y, y_rk4, y_rk5, atol, rtol, y.length)
        if(err <= 1) {
          break
        } else { 
          h = hdefault * Math.min(facmax, Math.max(facmin, fac * Math.pow(err, -0.2)))//;console.log("change", h )
        }
        // if ((t1 -t) < 1e-16) t = t1
      } // while loop
      yaccept = y
      y = y_rk5
      t += h
      if (t > t1) {
        t = t - h
        h = t1 - t
        y= yaccept
        accept = 1
      }
      if (t == t1) break
    }// while loop
    res.push([t, ...y]);console.log(t, y)
  } // for loop
  return res
} // function

    
var r = ode (times, y0, f2)
// // console.log(r)
const createCsv = require('csv-writer').createArrayCsvWriter;
const csv = createCsv({
    header: [],
    path: './r.csv'
  })
csv.writeRecords(r).then(() => {
    console.log('...Done')
  })
// maxerr = function(y0, y1, y2, atol, rtol, n) {
//   var serr = 0, scal, delta;
//   for (let i = 0; i < n; i++) {
//     /* y2 is used to estimate next y-value */
//     scal  = atol[i] + Math.max(Math.abs(y0[i]), Math.abs(y2[i])) * rtol[i]
//     delta = Math.abs(y2[i] - y1[i])
//     if (scal > 0) serr += Math.pow(delta/scal, 2);
//   }
//   return(sqrt(serr/n)); /* Euclidean norm */
// }

// /*====================================================================*/
//     /*      stepsize adjustment                                           */
//     /*====================================================================*/
//     //hmax, hmin  alpha, beta, errold, maxscale, minscale
// var hmax = //deltaT
// let safe = 0.9, maxscale = 10, minscale = 0.2 //rkauto.c
//     err = maxerr(y0, y1, y2, atol, rtol, neq);
//     hNew = h;
//     if (err == 0) {   use max scale if all tolerances are zero 
//       hNew = Math.min(h * 10, hmax);
//       errold = Math.max(err, 1e-4); /* 1e-4 taken from Press et al. */
//       accept = TRUE;
//     } else if (err < 1.0) {
//       /* increase step size only if last one was accepted */
//       if (accept) {
//         hNew = Math.min(hmax, h * Math.min(safe * pow(err, -alpha) * pow(errold, beta), maxscale));
//       }
//       errold = Math.max(err, 1e-4); /* 1e-4 taken from Press et al. */
//       accept = TRUE;
//     } else if (err > 1.0) {
//       nreject++;    /* count total number of rejected steps */
//       accept = FALSE;
//       hNew = h * Math.max(safe * pow(err, -alpha), minscale);
//     }

//     if (hNew < hmin) {
//       accept = TRUE;
//       // if (verbose) Rprintf("warning, h < Hmin\n");
//       istate[0] = -2;
//       hNew = hmin;
//     }
