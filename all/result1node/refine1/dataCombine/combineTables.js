

fs= require('fs')

combineTables = function (modeltype, runs) {
  let allSets = [] 

  for ( run = 0; run < runs.length; run++){
    var table = [], dataset = []
    var temp
    var data = fs.readFileSync('./' + modeltype + runs[run] +'.csv').toString()
    var lines = data.split('\n')
    for (let i = 1; i < lines.length; i++) {
      temp = lines[i].split(',')
      if (temp[temp.length - 1] !== '"NaN"' && temp[temp.length - 1] < 0) {
        table.push(temp)
      }
    }

    table.sort(sortFunction)

    var newSet = {}
    table.forEach(function(arr){
      newSet[arr.join("|")] = arr
    })

    var dataset = Object.keys(newSet).map(function(k){
      return newSet[k]
    })

    console.log(run, dataset.length)
    allSets.push(...dataset) 
  }

  allSets.sort(sortFunction)
  
  const createCsvWriter = require('csv-writer').createArrayCsvWriter;
  const csvWriter = createCsvWriter({
    header: ['R0', 'amplitude', 'gamma', 'mu', 'sigma', 'rho', 'psi', 'S_0', 'E_0', 'R_0','I_0','LogLik' ],
    path: './all.csv'
  }) 
  csvWriter.writeRecords(allSets)
    .then(() => {
    console.log('...done')
  })
}
combineTables ('sobolres', ['1', '2','3','4','5'])
function sortFunction(a, b) {
  if (Number(a[a.length - 1]) === Number(b[a.length - 1])) {
    return 0
  }
  else {
    return (Number(a[a.length - 1]) < Number(b[a.length - 1])) ? 1 : -1;
  }
}