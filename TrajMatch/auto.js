rkMethods = require('./rkmethods')
let snippet = require('./modelSnippet.js')

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
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// start from y = initial value
rk = function(y, times, func, parms, rtol = 1e-6, atol = 1e-6,
  tcrit = null, hmin = 0, hmax = null, hini = hmax,
  ynames = true, /*method = rkMethods("rk45dp7", ... ),*/ maxsteps = 5000,
  nout = 0) {

  if (hmax == null) {
    hmax = (times == null)? 0 : maxDiff(times) // checkInput in functions.R
  }
 
  if (hmax == 0) hmax = Number.MAX_SAFE_INTEGER
  n = y.length
  Nstates = y.length

  if (maxsteps < 0)       throw "maxsteps must be positive"
  if (!isFinite(maxsteps)) maxsteps = Number.MAX_SAFE_INTEGER - 1
  if (tcrit == null) {
    tcrit = Math.max(...times)
  }

  atol = new Array(Nstates).fill(atol)
  rtol = new Array(Nstates).fill(rtol)

  nsteps  = Math.min(Number.MAX_SAFE_INTEGER - 1,
                   Math.max(maxsteps * times.length,    // max. total number of steps
                       maxDiff(times) / hini + 1))   // but not less than required

  // Methods with variable step size
  if (hini == null) {
    hini = hmax
  }
  
  out = call_rkAuto (y, times,
    func, 0, parms, 0, 0,
    0, 0, atol,
    rtol, tcrit, 0,
    hmin, hmax, hini,
    0, 0, 0,
    nsteps, 0)
// console.log(out)
return out 
} 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
// start at initial and return the whole matrix of times and values
 call_rkAuto = function ( Xstart,  Times,  Func,  Initfunc,
   Parms,  eventfunc,  elist,  Nout,  Rho,
   Atol,  Rtol,  Tcrit,  Verbose,
   Hmin,  Hmax,  Hini,  Rpar,  Ipar,
   Method,  Maxsteps,  Flist) {
        /**  Initialization **/
  let nprot = 0;
  let tt = null, xs = null;

  let y, f, Fj, tmp, FF, rr;
  let R_yout;
  let y0, y1, y2, dy1, dy2, out, yout;
  let yt = []
  let errold = 0, t, dt, tmax;

  let R_FSAL, Alpha, Beta;
  let fsal = false       /* assume no FSAL */

  let interpolate = true
  let i = 0, j = 0, it = 0, it_tot = 0, it_ext = 0, nt = 0, neq = 0, it_rej = 0;
  let isForcing, isEvent;
   /*------------------------------------------------------------------------*/
  /* Processing of Arguments                                                */
  /*------------------------------------------------------------------------*/
  let rtol = [], atol = []

  for (let j = 0; j < Rtol.length; j++) rtol[j] = Rtol[j];
  for (let j = 0; j < Atol.length; j++) atol[j] = Atol[j];

  let  tcrit = Tcrit
  let  hmin  = Hmin
  let  hmax  = Hmax
  let  hini  = Hini
  let  maxsteps = Maxsteps
  // int  nout     = INTEGER(Nout)[0]; /* number of global outputs is func is in a DLL */

  let stage = 6//method.stage

  let R_A, R_B1, R_B2, R_C, R_D, R_densetype;
  let  A, bb1, bb2 = null, cc = null, dd = null;
  // R_A = method.a
  // R_B1 = method.b1 
  // R_B2 = method.b2 
  // R_C = method.c 
  // R_D = method.d
  /* dense output Cash-Karp: densetype = 2 */
  // let densetype = 0;
  // R_densetype = Method.densetype
  // if (R_densetype.length) densetype = R_densetype

  let  qerr = 4//Method.Qerr
  let  beta = 0;      /* 0.4/qerr; */
  let  alpha = 1 / qerr - 0.75 * beta;
  tt = Times
  nt = Times.length
  xs  = Xstart
  neq = Xstart.length

  y0  =  new Array(neq)
  y1  =  new Array(neq)
  y2  =  new Array(neq)
  dy1 =  new Array(neq)
  dy2 =  new Array(neq)
  f   =  new Array(neq)
  y   =  new Array(neq)
  Fj  =  new Array(neq)
  tmp =  new Array(neq)
  FF  =  new Array(stage).fill(Array(neq))
  rr  =  new Array(neq * 5)

  /* matrix for polynomial interpolation */
  let R_nknots;
  let nknots = 6;  /* 6 = 5th order polynomials by default*/
  let iknots = 0;  /* counter for knots buffer */
  let yknots;
  // PROTECT(R_nknots = getListElement(Method, "nknots")); nprot++;
  // if (length(R_nknots)) nknots = INTEGER(R_nknots)[0] + 1;
  // if (nknots < 2) {nknots = 1; interpolate = FALSE;}
  // if (densetype > 0) interpolate = true;
  // yknots = (double*) R_alloc((neq + 1) * (nknots + 1), sizeof(double));

  /* matrix for holding states and global outputs */
  let nout = 0//defined
  yout = new Array(nt).fill(Array(neq + nout + 1))

  /* attribute that stores state information, similar to lsoda */
  let istate = Array(22)
  istate[0] = 0; /* assume succesful return */
  for (i = 0; i < 22; i++) istate[i] = 0;
  /*------------------------------------------------------------------------*/
  /* Initialization of Integration Loop                                     */
  /*------------------------------------------------------------------------*/
  yout[0][0]   = tt[0];              /* initial time                 */
  // yknots[0] = tt[0];              /* for polynomial interpolation */
  for (i = 0; i < neq; i++) {
    y0[i]        = xs[i];         /* initial values               */
    yout[0][(i + 1)] = y0[i];   /* output array                 */
    // yknots[iknots + nknots * (i + 1)] = xs[i]; /* for polynomials */
  }
  iknots++;

  t = tt[0];
  tmax = Math.max(tt[nt - 1], tcrit);
  dt   = Math.min(hmax, hini);
  hmax = Math.min(hmax, tmax - t);

  /* Initialize work arrays (to be on the safe side, remove this later) */
  for (i = 0; i < neq; i++)  {
    y1[i] = 0;
    y2[i] = 0;
    Fj[i] = 0;
    for (j= 0; j < stage; j++)  { // FF is a matrix stage*neq
      FF[j][i] = 0; //FF[i + j * neq] = 0;
    }
  } 
  /*------------------------------------------------------------------------*/
  /* Main Loop                                                              */
  /*------------------------------------------------------------------------*/
  it     = 1; /* step counter; zero element is initial state   */
  it_ext = 0; /* counter for external time step (dense output) */
  it_tot = 0; /* total number of time steps                    */
  it_rej = 0;
  
  /* integrate separately between external time steps; do not interpolate */
  for (let j = 0; j < nt - 1; j++) {
    t = tt[j];
    tmax = Math.min(tt[j + 1], tcrit);
    dt = tmax - t;
    // if (isEvent) {
    //   updateevent(&t, y0, istate);
    // }
    // if (verbose) Rprintf("\n %d th time interval = %g ... %g", j, t, tmax);

    y22 = rk_auto(
      fsal, neq, stage, 0, 0, 0, 0, 0,
      0, maxsteps, nt,
      0, it, it_ext, it_tot, it_rej,
      0, 0,
      t,  tmax, hmin, hmax, alpha, beta,
      dt, errold,
      tt, y0, y1, y2, dy1, dy2, f, y, Fj, tmp, FF, rr, A,
      out, bb1, bb2, cc, dd, atol, rtol, yknots, yout,
      Func, Parms, Rho
    );

    /* in this mode, internal interpolation is skipped,
        so we can simply store the results at the end of each call */
    // yout[j + 1][0] = tmax;
    // for (i = 0; i < neq; i++) yout[j + 1][1 + i] = y2[i];
    yt.push(y22)
  }

    /*====================================================================*/
  /* call derivs again to get global outputs                            */
  /* j = -1 suppresses unnecessary internal copying                     */
  /*====================================================================*/
  // if (nout > 0) {
  //   for (int j = 0; j < nt; j++) {
  //     t = yout[j];
  //     for (i = 0; i < neq; i++) tmp[i] = yout[j + nt * (1 + i)];
  //     derivs(Func, t, tmp, Parms, Rho, FF, out, -1, neq, ipar, isDll, isForcing);
  //     for (i = 0; i < nout; i++) {
  //       yout[j + nt * (1 + neq + i)] = out[i];
  //     }
  //   }
  // }
    
  return(yt);
}  
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// use rk method
rk_auto = function(
       /* integers */
        fsal,  neq,  stage,
        isDll,  isForcing,  verbose,
        nknots,  erpolate,  densetype,  maxsteps,  nt,
       /*  pointers */
        _iknots,  _it,  _it_ext,  _it_tot,  _it_rej,
        istate,   ipar,
       /* double */
        t,  tmax,  hmin,  hmax,
        alpha,  beta,
       /*  pointers */
        _dt,  _errold,
       /* arrays */
        tt,  y0,  y1,  y2,  dy1,  dy2,
        f,  y,  Fj,  tmp,
        FF,  rr,  A,  out,
        bb1,  bb2,  cc,  dd,
        atol,  rtol,  yknots,  yout,
       /* SEXPs */
        Func,  Parms,  Rho)
{

  let i = 0, j = 0, j1 = 0, k = 0, accept = false, nreject = _it_rej, one = 1; 
  let DBL_EPSILON = 1e-9
  // int iknots = *_iknots, it = *_it, it_ext = *_it_ext, it_tot = *_it_tot;
  let err, dtnew, t_ext;
  let dt = _dt, errold = _errold;

  /* todo: make this user adjustable */
  let minscale = 0.2, maxscale = 10.0, safe = 0.9;
  

  /*------------------------------------------------------------------------*/
  /* Main Loop                                                              */
  /*------------------------------------------------------------------------*/
  do {
    y = y0
  //   /*  save former results of last step if the method allows this
  //      (first same as last)                                             */
  //   /* Karline: improve by saving "accepted" FF, use this when rejected */
  //   if (fsal && accept){
  //     j1 = 1;
  //     for (let i = 0; i < neq; i++) FF[0][i] = FF[stage - 1][i]
  //   } else {
  //     j1 = 0;
  //   }
  //   /******  Prepare Coefficients from Butcher table ******/
  //   for (let j = j1; j < stage; j++) {
  //     for(let i = 0; i < neq; i++) Fj[i] = 0;
  //       k = 0;
  //       while(k < j) {
  //         for(i = 0; i < neq; i++)
  //           Fj[i] = Fj[i] + A[k][j] * FF[k][i] * dt;
  //         k++;
  //       }
  //       for (int i = 0; i < neq; i++) {
  //         tmp[i] = Fj[i] + y0[i];
  //       }
  //       /******  Compute Derivatives ******/
  //       // /* pass option to avoid unnecessary copying in derivs */
  //       // derivs(Func, t + dt * cc[j], tmp, Parms, Rho, FF, out, j, neq,
  //              // ipar, isDll, isForcing);    
  //   }

    /*====================================================================*/
    /* Estimation of new values                                           */
    /*====================================================================*/

    // /* use BLAS wrapper with reduced error checking */
    // blas_matprod1(FF, neq, stage, bb1, stage, one, dy1);
    // blas_matprod1(FF, neq, stage, bb2, stage, one, dy2);

    // it_tot++; /* count total number of time steps */
    // for (i = 0; i < neq; i++) {
    //   y1[i] = y0[i] + dt * dy1[i];
    //   y2[i] = y0[i] + dt * dy2[i];
    // }
    console.log(t)
    k1 = sp(dt , func( t, y))
    k2 = sp(dt , func( t + dt/4, sum([y,sp(1/4, k1)])))
    k3 = sp(dt , func( t + 3*dt/8, sum([y, sp(3/32,k1),sp(9/32, k2)]))) 
    k4 = sp(dt , func( t + 12*dt/13, sum([y, sp(1932/2197,k1),sp(-7200/2197, k2),sp(7296/2197, k3)]))) 
    k5 = sp(dt , func( t + dt, sum([y, sp(439/216, k1),sp(-8,k2), sp(3680/513, k3), sp(-845/4104, k4)])))
    k6 = sp(dt , func( t + dt/2 ,sum([y, sp(-8/27, k1),sp(2, k2), sp(-3544/2565, k3), sp(1859/4104,k4), sp(-11/40, k5)])))
    y1= sum([y, sp(25/216, k1), sp(1408/2565, k3), sp(2197/4104, k4), sp(-1/5, k5)])                       // compute RK4
    y2 = sum([y, sp(16/135, k1), sp(6656/12825, k3) ,sp(28561/56430, k4), sp( -9/50, k5) ,sp( 2/55, k6)]);  // compute RK5
  
    /*====================================================================*/
    /*      stepsize adjustment                                           */
    /*====================================================================*/

    err = maxerr(y0, y1, y2, atol, rtol, neq)
    dtnew = dt;
    if (err == 0) {  /* use max scale if all tolerances are zero */
      dtnew  = Math.min(dt * 10, hmax);
      errold = Math.max(err, 1e-4); /* 1e-4 taken from Press et al. */
      accept = true;
    } else if (err < 1.0) {
      /* increase step size only if last one was accepted */
      if (accept)
        dtnew = Math.min(hmax, dt *
          Math.min(safe * Math.pow(err, -alpha) * Math.pow(errold, beta), maxscale));
      errold = Math.max(err, 1e-4); /* 1e-4 taken from Press et al. */
      accept = true;
    } else if (err > 1.0) {
      nreject++;    /* count total number of rejected steps */
      accept = FALSE;
      dtnew = dt * Math.max(safe * Math.pow(err, -alpha), minscale);
    }

    if (dtnew < hmin) {
      accept = true;
      // if (verbose) Rprintf("warning, h < Hmin\n");
      istate[0] = -2;
      dtnew = hmin;
    }
    /*====================================================================*/
    /*      Interpolation and Data Storage                                */
    /*====================================================================*/
    interpolate = 0
    if (accept) {
      if (interpolate) {
      /*--------------------------------------------------------------------*/
      /* case A1) "dense output type 1": built-in polynomial interpolation  */
      /* available for certain rk formulae, e.g. for rk45dp7                */
      /*--------------------------------------------------------------------*/
      if (densetype == 1) {
        denspar(FF, y0, y2, dt, dd, neq, stage, rr);
        t_ext = tt[it_ext];
        while (t_ext <= t + dt) {
          densout(rr, t, t_ext, dt, tmp, neq);
          /* store outputs */
          if (it_ext < nt) {
            yout[it_ext] = t_ext;
            for (i = 0; i < neq; i++)
              yout[it_ext + nt * (1 + i)] = tmp[i];
          }
          if(it_ext < nt-1) t_ext = tt[++it_ext]; else break;
        }

        /*--------------------------------------------------------------------*/
        /* case A2) dense output type 2: the Cash-Karp method                 */
        /*--------------------------------------------------------------------*/
      } else if (densetype == 2)  {   /* dense output method 2 = Cash-Karp */
        derivs(Func, t + dt, y2, Parms, Rho, dy2, out, 0, neq,
               ipar, isDll, isForcing);

        t_ext = tt[it_ext];

        while (t_ext <= t + dt) {
          densoutck(t, t_ext, dt, y0, FF, dy2, tmp, neq);
          /* store outputs */
          if (it_ext < nt) {
            yout[it_ext] = t_ext;
            for (i = 0; i < neq; i++)
              yout[it_ext + nt * (1 + i)] = tmp[i];
          }
          if(it_ext < nt-1) t_ext = tt[++it_ext]; else break;
       }
       /* FSAL (first same as last) for Cash-Karp */
       for (i = 0; i < neq; i++) FF[i + neq * (stage - 1)] = dy2[i] ;

        /*--------------------------------------------------------------------*/
        /* case B) Neville-Aitken-Interpolation for integrators               */
        /* without dense output                                               */
        /*--------------------------------------------------------------------*/
      } else {
          /* (1) collect number "nknots" of knots in advance */
          yknots[iknots] = t + dt;   /* time is first column */
          for (i = 0; i < neq; i++) yknots[iknots + nknots * (1 + i)] = y2[i];
          if (iknots < (nknots - 1)) {
            iknots++;
          } else {
      // /* (2) do polynomial interpolation */
      //       t_ext = tt[it_ext];
      //       while (t_ext <= t + dt) {
      //         neville(yknots, &yknots[nknots], t_ext, tmp, nknots, neq);
      //         /* (3) store outputs */
      //         if (it_ext < nt) {
      //           yout[it_ext] = t_ext;
      //           for (i = 0; i < neq; i++)
      //             yout[it_ext + nt * (1 + i)] = tmp[i];
      //         }
      //         if(it_ext < nt-1) t_ext = tt[++it_ext]; else break;
      //       }
      //       shiftBuffer(yknots, nknots, neq + 1);
          }
        }
      } else {
        /*--------------------------------------------------------------------*/
        /* Case C) no interpolation at all (for step to step integration);    */
        /*         results are stored after the call                          */
        /*--------------------------------------------------------------------*/
      }
      /*--------------------------------------------------------------------*/
      /* next time step                                                     */
      /*--------------------------------------------------------------------*/
      t = t + dt;
      // it++;
      for (i=0; i < neq; i++) y0[i] = y2[i];

    } /* else rejected time step */
    dt = Math.min(dtnew, tmax - t)
    // if (it_ext > nt) {
    //   Rprintf("error in RK solver rk_auto.c: output buffer overflow\n");
    //   break;
    // }
    // if (it_tot > maxsteps) {
    //   istate[0] = -1;
    //   warning("Number of time steps %i exceeded maxsteps at t = %g\n", it, t);
    //   break;
    // }
    /* tolerance to avoid rounding errors */
  } while (t < (tmax - 100.0 * DBL_EPSILON * dt)); /* end of rk main loop */
  return y2
  /* return reference values */
  // _iknots = iknots; _it = it; _it_ext = it_ext; _it_rej = nreject;
  // _it_tot = it_tot; _dt = dtnew; _errold = errold;
}


// helper
maxDiff = function (arr) {
  let d = 0
  for ( i = 1 ; i < arr.length; i++) {
    d1 = Math.abs(arr[i] - arr[i - 1])
    if ( d1 > d) d = d1
  }
  return d
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

////////////////////////////////////////
let  y = [1, 0, 0.9]
const func = function(t, y) {
  var dy1, dy2, dy3
  dy1 = -2 * y[1] * y[2]
  dy2 = 1.25* y[0] * y[2]
  dy3 = -0.5* y[0] * y[1]
  return [dy1, dy2, dy3]
}

let times = []
for (let i = 0; i <= 20; i = istep + .1) {
  istep = Number(i.toFixed(8))
  times.push(istep)
}
let hmax = .1
let parms = []
p = rk(y, times, func, parms,  1e-6, 1e-6,
  null,  0,  null,  hmax,
  true, 5000,
   0)


const createCsvWriter = require('csv-writer').createArrayCsvWriter;
  const csvWriter = createCsvWriter({
    header: [],
    path: './p.csv'
  })   
  csvWriter.writeRecords(p)
    .then(() => {
    console.log('...Done')
  })