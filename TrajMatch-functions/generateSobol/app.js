/* 
*  User needs to define lower bounds and upper bounds for each parameter
*  and also the number of points to be generated. 
*/


let sobolSeq = require('./sobolSeq.js');

let lowerBounds = {R0:0, amplitude:0, gamma:73.05, mu:0, sigma:45.66, rho:0, psi:0, S_0:0, E_0:0, I_0:0, R_0:0}
let upperBounds = {R0:200, amplitude:1, gamma:73.05, mu:1, sigma:45.66, rho:1, psi:1, S_0:0.1, E_0:1e-4, I_0:1e-4, R_0:0.9}
let numberOfPoints = 2000

let sobolSet = sobolSeq.sobolDesign( lowerBounds,  upperBounds, numberOfPoints)

// Write the result as a csv file
var headerArray =[]
for(i = 0; i < Object.values(lowerBounds).length; ++i) {
  headerArray.push(Object.keys(upperBounds)[i])
}

const createCsvWriter = require('csv-writer').createArrayCsvWriter
const csvWriter = createCsvWriter({
  header: headerArray, 
  path: "./paramSetSobol.csv"
})  
csvWriter.writeRecords(sobolSet)
  .then(() => {
  console.log('...done')
})
