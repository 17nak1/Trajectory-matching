/**
 *  @file      combineTables.js       This function clean, sort and combine tables in csv files.
 *                                    modelType : Name of the model, include: 'DeterministicSEIR', 'StochasticSEIR'
 *                                    runs : The array of parameters that wants to combine their tables, include; [R0Index, AMPLITUDE, MU, RHO, PSI]
 *  @return                           A sorted table based on the "LogLik" column.
 *  @reference                        Aaron A. King.
 *  @author    Nazila Akhavan
 *  @date      March 2019
 */

fs= require('fs')
combineTables = {}

combineTables.combine = function (modelType, runs) {
  let allSets = [] 
  var newSet
  var dataset
  for ( run = 0; run < runs.length; run++){
    var table = [], dataset = []
    var temp
    var data = fs.readFileSync('./' + modelType + runs[run] +'.csv').toString()
    var lines = data.split('\n')
    for (let i = 1; i < lines.length; i++) {
      temp = lines[i].split(',')
      if (temp[temp.length - 1] !== '"NaN"' && temp[temp.length - 1] < 0) {
        table.push(temp)
      }
    }

    table.sort(sortFunction)

    newSet = {}
    table.forEach(function(arr){
      newSet[arr.join("|")] = arr
    })
    dataset = Object.keys(newSet).map(function(k){
      return newSet[k]
    })
    allSets.push(...dataset) 
  }

  allSets.sort(sortFunction)
  let finalSet = [allSets[0]]
  let size = allSets[0].length - 1
  for (let i = 1; i < allSets.length; i++) {
    if(allSets[i - 1][0] !== allSets[i][0]) {
      finalSet.push(allSets[i])
    } else {
      if(allSets[i - 1][size] !== allSets[i][size]) {
        finalSet.push(allSets[i])
      }
    } 
  }
  
  const createCsvWriter = require('csv-writer').createArrayCsvWriter;
  const csvWriter = createCsvWriter({
    header: ['R0', 'amplitude', 'gamma', 'mu', 'sigma', 'rho', 'psi', 'S_0', 'E_0', 'I_0','R_0','LogLik' ],
    path: './allcombined.csv'
  })
  csvWriter.writeRecords(finalSet)
    .then(() => {
    console.log('...done')
  })
}


// Helper function
function sortFunction(a, b) {
  if (Number(a[a.length - 1]) === Number(b[a.length - 1])) {
    return 0
  }
  else {
    return (Number(a[a.length - 1]) < Number(b[a.length - 1])) ? 1 : -1;
  }
}
module.exports = combineTables