/* 
 *                                    User needs to provide values for 
 *                                    modelType : Name of the model, include: 'DeterministicSEIR', 'StochasticSEIR'
 *                                    runs : The array of parameters that wants to combine their tables, include; [R0Index, AMPLITUDE, MU, RHO, PSI]
 */

let combineTables = require('./combineTables.js')

combineTables.combine ('', ['all', 'R0','amp', 'mu', 'rho', 'psi'])