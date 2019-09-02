//============================================================================
//Interface to a generalized code for solving explicit variable and fixed
//step ODE solvers of the Runge-Kutta family, see helpfile for details.
//============================================================================

rk = function(y, times, func, parms, rtol = 1e-6, atol = 1e-6,
  verbose = false, tcrit = null, hmin = 0, hmax = null, hini = hmax,
  ynames = true, method = "rk45f", maxsteps = 5000,
  dllname = null, initfunc = dllname, initpar = parms,
  rpar = null,  ipar = null, nout = 0, outnames = null, forcings = null,
  initforc = null, fcontrol = null, events = null) {

  

  // if (is.list(func)) {            // a list of compiled functions
  //     if (!is.null(initfunc) & "initfunc" %in% names(func))
  //        stop("If 'func' is a list that contains initfunc, argument 'initfunc' should be null")
  //     if (!is.null(dllname) & "dllname" %in% names(func))
  //        stop("If 'func' is a list that contains dllname, argument 'dllname' should be null")
  //     if (!is.null(initforc) & "initforc" %in% names(func))
  //        stop("If 'func' is a list that contains initforc, argument 'initforc' should be null")
  //     if (!is.null(events$func) & "eventfunc" %in% names(func))
  //        stop("If 'func' is a list that contains eventfunc, argument 'events$func' should be null")
  //     if ("eventfunc" %in% names(func)) {
  //        if (! is.null(events))
  //          events$func =func$eventfunc
  //        else
  //          events =list(func = func$eventfunc)
  //     }
  //    if (!is.null(func$initfunc)) initfunc =func$initfunc
  //    if (!is.null(func$dllname))  dllname =func$dllname
  //    if (!is.null(func$initforc)) initforc =func$initforc
  //    func =func$func
  // }
  method = rkMethod(method)
  varstep = method.varstep
  if (!varstep && (hmin != 0 | !is.null(hmax))) {
    throw "'hmin' and 'hmax' are ignored (fixed step Runge-Kutta method).\n"
  }

  if (hmax == null) {
    hmax = (times == null)? 0 : maxDiff(times) // checkInput in functions.R
  }
  if (hmax == 0) hmax = Number.MAX_SAFE_INTEGER
  n = y.length
  

  if (maxsteps < 0)       throw "maxsteps must be positive"
  if (!isFinite(maxsteps)) maxsteps = Number.MAX_SAFE_INTEGER - 1
  if (tcrit == null) {
    tcrit = Math.max(...times)
  }

  // ToDo: check for nonsense-combinations of densetype and d
  if (!isNaN(method.densetype)) {
    // make this an integer to avoid errors on the C level
    method,densetype = method.densetype
    if (!(method.densetype === 1 || method.densetype = 2)) {
      throw "Unknown value of densetype; set to null."
      method.densetype = null
    }
  } 

  // Checks and ajustments for Neville-Aitken interpolation
  // - starting from deSolve >= 1.7 this interpolation method
  //   is disabled by default.
  // - Dense output for special RK methods is enabled and
  //   all others adjust internal time steps to hit external time steps
  if (isNaN(method.nknots)) {
      method.nknots = 0
  } else {
    method.nknots = Math.ceil(method.nknots)
  }

  nknots <- method.nknots
  if (nknots > 8) {
      throw "Large number of nknots does not make sense."
  } else if (nknots < 2) {
    //method without or with disabled interpolation
    method.nknots = 0
  } else {
    trange = diff(Math.min(...times), Math.max(...times))
    //ensure that we have at least nknots + 2 data points; + 0.5 for safety)
    //to allow 3rd order polynomial interpolation
    //for methods without built-in dense output
    if ((isNaN(method.d) &&                        // has no "dense output"?
          isNaN(method.densetype) &&                // or no dense output type
      (hmax > 1 / (nknots + 2.5) * trange))) {      // or time steps too large?
      //in interpolation mode: automatic adjustment of step size arguments
      //to ensure the required minimum of knots
      hmax = 1 / (nknots + 2.5) * trange
      hini = hmax
      if (hmin < hini) {
        hmin = hini
      }
      console.log("\nNote: Method "+ method.ID + " needs intermediate steps for interpolation\n")
      console.log("hmax decreased to" + hmax + "\n")
    }
  }


    // Model as shared object (DLL)?
    Ynames = null //attr(y, "names");
    Initfunc = null
    Eventfunc = null
    events =checkevents(events, times, Ynames, dllname)
    if (! is.null(events$newTimes)) times =events$newTimes

    // dummy forcings
    flist    <-list(fmat = 0, tmat = 0, imat = 0, ModelForc = null)
    Nstates =length(y) //assume length of states is correct

    // function specified in a DLL or inline compiled
    if (is.character(func) | class(func) == "CFunc") {
      DLL =checkDLL(func, null, dllname,
                      initfunc, verbose, nout, outnames)

      Initfunc  =DLL$ModelInit
      Func      =DLL$Func
      Nglobal   =DLL$Nglobal
      Nmtot     =DLL$Nmtot
      Eventfunc =events$func

      if (! is.null(forcings))
        flist =checkforcings(forcings, times, dllname, initforc, verbose, fcontrol)

      if (is.null(ipar)) ipar =0
      if (is.null(rpar)) rpar =0
      // preparation for events in R if function is a DLL (added by KS)
      if (is.function(Eventfunc))
        rho =environment(Eventfunc)
      else
        rho =null

    } else {
      // parameter initialisation not needed if function is not a DLL
      initpar =null
      rho =environment(func)

      // func is overruled, either including ynames, or not
      // This allows to pass the "..." arguments and the parameters
      if(ynames) {
        Func   =function(time, state, parms){
          attr(state, "names") =Ynames
          func(time, state, parms, ...)}
        if (! is.null(events$Type))
          if (events$Type == 2)
            Eventfunc =function(time, state) {
              attr(state, "names") =Ynames
              events$func(time, state, parms, ...)
            }
      } else {                            // no ynames...
        Func   =function(time, state, parms)
          func(time, state, parms, ...)
        if (! is.null(events$Type))
          if (events$Type == 2)
            Eventfunc =function(time, state)
              events$func(time, state, parms, ...)
      }

      // Call func once to figure out whether and how many "global"
      // results it wants to return and some other safety checks
      FF =checkFuncEuler(Func, times, y, parms, rho, Nstates)
      Nglobal =FF$Nglobal
      Nmtot   =FF$Nmtot

      if (! is.null(events$Type))
        if (events$Type == 2) checkEventFunc(Eventfunc, times, y, rho)
    }

    // handle length of atol and rtol
    if (Nstates %% length(atol))
      warning("length of atol does not match number of states")
    if (Nstates %% length(rtol))
      warning("length of rtol does not match number of states")

    atol =rep(atol, length.out = Nstates)
    rtol =rep(rtol, length.out = Nstates)

    // Number of steps until the solver gives up
    // nsteps  =min(.Machine$integer.max -1, maxsteps * length(times))

    // total number of time steps is set to
    // average number per time step * number of time steps
    // but not less than required for the largest time step with given hini

    nsteps  =min(.Machine$integer.max - 1,
                   max(maxsteps * length(times),    // max. total number of steps
                       max(diff(times))/hini + 1)   // but not less than required
               )

    vrb =false //true forces some internal debugging output of the C code
    // Implicit methods
    on.exit(.C("unlock_solver"))
    implicit =method$implicit
    if (is.null(implicit)) implicit =0
    if (implicit) {
      if (is.null(hini)) hini =0
      out =.Call("call_rkImplicit", as.double(y), as.double(times),
        Func, Initfunc, parms, Eventfunc, events,
        as.integer(Nglobal), rho,
        as.double(tcrit), as.integer(vrb),
        as.double(hini), as.double(rpar), as.integer(ipar), method,
        as.integer(nsteps), flist)

    } else if (varstep) {// Methods with variable step size
      if (is.null(hini)) hini =hmax
      out =.Call("call_rkAuto", as.double(y), as.double(times),
        Func, Initfunc, parms, Eventfunc, events,
        as.integer(Nglobal), rho, as.double(atol),
        as.double(rtol), as.double(tcrit), as.integer(vrb),
        as.double(hmin), as.double(hmax), as.double(hini),
        as.double(rpar), as.integer(ipar), method,
        as.integer(nsteps), flist)
    } else { // Fixed step methods
      // hini = 0 for fixed step methods means
      // that steps in "times" are used as they are
      if (is.null(hini)) hini =0
      out =.Call("call_rkFixed", as.double(y), as.double(times),
        Func, Initfunc, parms, Eventfunc, events,
        as.integer(Nglobal), rho,
        as.double(tcrit), as.integer(vrb),
        as.double(hini), as.double(rpar), as.integer(ipar), method,
        as.integer(nsteps), flist)
    }

    // output cleanup
    out =saveOutrk(out, y, n, Nglobal, Nmtot,
                     iin = c(1, 12:15), iout = c(1:3, 13, 18))

    attr(out, "type") ="rk"
    if (verbose) diagnostics(out)
    return(out)
}
