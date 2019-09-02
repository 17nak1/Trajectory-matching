

ode2 = {}

ode2.ode = function (times, y, f2, params, intPop, intBirth) {
  var t = times[0] , t1, R, hdefault = .1 * (times[1] - t)
  var k1 = [], k2 = [], k3 = [], k4 = [], k5 = [], k6 = []
  var res = [], y_rk4 = [], y_rk5 = [], difmax = 0, count, flag = 0, accept = 0, yaccept, h
  let facmin = 0.2, fac = 0.9 , facmax = 10
  let n = y.length
  let atol = new Array(n).fill(1e-16)
  let rtol = new Array(n).fill(1e-16)
  let err = 0
  for (let i = 1; i < 2 ; i++) {//times.length - 1
    t1 = times[i]
    h= hdefault
    // console.log(t,t1, h, count)
  //   t = 1944.0383258648972
  //   y = [ 113880.05379259076,
  // 190.7672492366465,
  // 1285107.828471409,
  // 106.15831151665354,
  // 205913.22059070622 ]
  //   h = 0.0000040461227399646305
    console.log(accept , h)
    while ( flag === 0) {
      count = 0
      // h= hdefault
      // console.log(t,t1, h)
      
      while(t < t1 && count < 2) {
        if ( (t1 - t ) < 1e-10) {
          t = t1 - h
          break
        }
        count++
        k1 = ode2.sp(h , f2( t, y, params, intPop(t), intBirth(t)))
        k2 = ode2.sp(h , f2( t + h/4, ode2.sum([y,ode2.sp(1/4, k1)]), params, intPop(t + h/4), intBirth(t + h/4) ))
        k3 = ode2.sp(h , f2( t + 3*h/8, ode2.sum([y, ode2.sp(3/32,k1),ode2.sp(9/32, k2)]), params, intPop(t + 3*h/8), intBirth(t + 3*h/8) )) 
        k4 = ode2.sp(h , f2( t + 12*h/13, ode2.sum([y, ode2.sp(1932/2197,k1),ode2.sp(-7200/2197, k2),ode2.sp(7296/2197, k3)]), params, intPop(t + 12*h/13), intBirth(t + 12*h/13))) 
        k5 = ode2.sp(h , f2( t + h, ode2.sum([y, ode2.sp(439/216, k1),ode2.sp(-8,k2), ode2.sp(3680/513, k3), ode2.sp(-845/4104, k4)]), params, intPop(t + h), intBirth(t + h)))
        k6 = ode2.sp(h , f2( t + h/2 ,ode2.sum([y, ode2.sp(-8/27, k1),ode2.sp(2, k2), ode2.sp(-3544/2565, k3), ode2.sp(1859/4104,k4), ode2.sp(-11/40, k5)]), params, intPop(t + h/2), intBirth(t + h/2)))
        y_rk4 = ode2.sum([y, ode2.sp(25/216, k1), ode2.sp(1408/2565, k3), ode2.sp(2197/4104, k4), ode2.sp(-1/5, k5)])                  // compute RK4
        y_rk5 = ode2.sum([y, ode2.sp(16/135, k1), ode2.sp(6656/12825, k3) ,ode2.sp(28561/56430, k4), ode2.sp( -9/50, k5) ,ode2.sp( 2/55, k6)]);       // compute RK5
        
        if ( accept == 1) {
          // accept = 0
          break
        }
        
        err = ode2.maxerr(y, y_rk4, y_rk5, atol, rtol, y.length)
        if(err <= 1) {//console.log("break", count)
          break
        } else { 
          // console.log("change", h, count )
          h = h * Math.min(facmax, Math.max(facmin, fac * Math.pow(err, -0.2)))//;console.log("change", h, count )
        }
        // console.log("end", t, t1,h, count)
      } // while loopin
      // console.log("end", t, t1,h, count)
      yaccept = y
      y = y_rk5
      // console.log(y, y_rk5 == yaccept)
      // console.log(y)
      t += h
      // console.log(t,t1)
      if (t > t1) {
        t = t - h
        h = t1 - t
        y= yaccept//;console.log("t>t1",h)
        accept = 1
      }
      if (t == t1) break
    } // while loop

    res.push([t, ...y]); console.log( y)
  } // for loop
  return res
} // function

ode2.sum = function (array) {
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

ode2.sp = function (scalar, array) {
  let sum = []
  for(i = 0; i < array.length; i++){
   sum.push(scalar * array[i]);
  }
  return sum
}

ode2.abs = function (array) {
  let sum = 0
  for(i = 0; i < array.length; i++){
   sum += Math.pow(Math.abs(array[i]), 2)
  }
  return Math.sqrt(sum)
}
ode2.maxerr = function(y0, y1, y2, atol, rtol, n) {
  var serr = 0, scal, delta;
  for (let i = 0; i < n; i++) {
    /* y2 is used to estimate next y-value */
    scal  = atol[i] + Math.max(Math.abs(y0[i]), Math.abs(y2[i])) * rtol[i]
    delta = Math.abs(y2[i] - y1[i])
    if (scal > 0) serr += Math.pow(delta/scal, 2);
  }
  return(Math.sqrt(serr/n)); /* Euclidean norm */
}

module.exports = ode2
