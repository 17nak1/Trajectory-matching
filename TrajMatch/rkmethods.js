var rkMethods = {}
sum = function (array) {
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

sp = function (scalar, array) {
  var sum = []
  for(i = 0; i < array.length; i++){
   sum.push(scalar * array[i]);
  }
  return sum
}
abs = function (array) {
  var sum = 0
  for(i = 0; i < array.length; i++){
   sum += Math.pow(Math.abs(array[i]), 2)
  }
  return Math.sqrt(sum)
}

  var rkMethods = {}
sum = function (array) {
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

sp = function (scalar, array) {
  var sum = []
  for(i = 0; i < array.length; i++){
   sum.push(scalar * array[i]);
  }
  return sum
}
abs = function (array) {
  var sum = 0
  for(i = 0; i < array.length; i++){
   sum += Math.pow(Math.abs(array[i]), 2)
  }
  return Math.sqrt(sum)
}


rkMethods.method = function (method) {
  let A, b1, b2, c, d
  let stage , Qerr
  switch (method) {
  case "rk45f":
    varstep = TRUE,
    method.A  = [[0, 0, 0, 0, 0],
          [1/4, 0, 0, 0, 0],
          [3/32, 9/32, 0, 0, 0],
          [1932/2197, -7200/2197, 7296/2197, 0, 0],
          [439/216, -8, 3680/513, -845/4104, 0],
          [-8/27, 2, -3544/2565, 1859/4104, -11/40]]
                
    method.b1 = [25/216,   0,  1408/2565,  2197/4104,  -1/5,   0]
    method.b2 = [16/135,   0,  6656/12825,   28561/56430,  -9/50,  2/55]
    method.c  = [0,  1/4,  3/8,  12/13,  1,  1/2]
    method.stage = 6
    method.Qerr  = 4
    
}
  module.exports = rkMethods