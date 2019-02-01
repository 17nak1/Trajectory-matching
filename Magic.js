/**
 *  @file       Magic.js  
 */

let traj_match = require('traj').traj_match
var dataCovar = [], dataCases = [], inputArr = [], init = [], res = []
var indx = new Array(12)

function start () {
  let req = new XMLHttpRequest()
  req.open('GET', './demo-compiled.js')
  req.onload = function () {
    code = req.responseText
  }
  var fileChooser = document.getElementById('file1-upload')
  fileChooser.onclick = function () {
    this.value = ''
  }
  document.getElementById('file1-upload').onchange = function () {
    document.getElementById('label-file1').innerHTML = 'Uploaded'
    document.getElementById('label-file1').style.backgroundColor = '#ffbf00'
    var file = this.files[0]
    var outArr = []
    dataCovar = []
    var reader = new FileReader()
    reader.onload = function () {
      var lines = this.result.split('\n')
      for (var line = 1; line < lines.length; line++) {
        outArr.push(lines[line].split(','))
      }
      dataCovar.push(outArr)
    }
    reader.readAsText(file)
  }
  var fileChooser = document.getElementById('file2-upload')
  fileChooser.onclick = function () {
    this.value = ''
  }
  document.getElementById('file2-upload').onchange = function () {
    document.getElementById('label-file2').innerHTML = 'Uploaded'
    document.getElementById('label-file2').style.backgroundColor = '#ffbf00'
    var file = this.files[0]
    var outArr = []
    dataCases = []
    var reader = new FileReader ()
    reader.onload = function () {
      var lines = this.result.split('\n')
      for (var line = 1; line < lines.length; line++) {
        outArr.push(lines[line].split(','))
      }
      dataCases.push(outArr)
    }
    reader.readAsText(file)
  }
  var fileChooser = document.getElementById('file3-upload')
  fileChooser.onclick = function () {
    this.value = ''
  }
  document.getElementById('file3-upload').onchange = function () {
    document.getElementById('label-file3').innerHTML = 'Uploaded'
    document.getElementById('label-file3').style.backgroundColor = '#ffbf00'
    var file = this.files[0]
    var outArr = []
    init = []
    var reader = new FileReader()
    reader.onload = function () {
      var lines = this.result.split('\n')
      for (var line = 1; line < lines.length; line++) {
        outArr.push(lines[line].split(','))
      }
      init.push(outArr)// Upload parameters
    }
    reader.readAsText(file)
  }
  let computeButton = document.querySelector('button#calc')
  let downloadButton = document.querySelector('button#download')
  downloadButton.style.display = 'none'
  computeButton.onclick = function () {
    inputArr = [], res = []
    computeButton.style.display = 'none'
    downloadButton.style.display = ''
    downloadButton.style.backgroundColor = '#D3D3D3'
    let setup1 = document.querySelector('#setup')
    let table = setup1.querySelector('table#table')
    let rows = table.querySelectorAll('tr')
    for (let p = 0; p < indx.length; p++) {
      if (indx[p] === 1) {
        document.querySelector('#setup').querySelector('table#table').querySelectorAll('tr')[p + 1].style.backgroundColor ='#2ed573'
      } else {
        rows[p + 1].bgColor = '#FFFFFF'
      }
    }
    for (let i = 1; i < rows.length; i++) {
      let row = rows[i]
      let cols = row.querySelectorAll('td')
      let cell = cols[cols.length - 1]
      var input = cell.querySelector('input#valInp').value;
      inputArr.push(Number(input))// Read parameters from the table
      if (init.length) {
        cell.querySelector('input#valInp').disabled = 'true'
        if (i !== 12) {
          cell.querySelector('input#valInp').value = ''
        } 
      }
    }
    let times = [inputArr[11], Number(dataCases[0][0][0]), Number(dataCases[0][dataCases[0].length - 2][0])]
    inputArr.pop()
    setTimeout(function () {
      if (init.length) {
        for (let i = 0; i < init[0].length - 1; i++) {
          var ans = traj_match(dataCovar[0], dataCases[0], init[0][i], times, indx)
          res.push(ans)
        }
        console.log(res)
      } else {
        var tem = inputArr[9]
        inputArr[9] = inputArr[10]
        inputArr[10] = tem
        res.push(traj_match(dataCovar[0], dataCases[0], inputArr, times, indx))
        console.log(res)
      }
      res.splice(0, 0, ['R0', 'amplitude', 'gamma', 'mu', 'sigma', 'rho', 'psi', 'S_0', 'E_0', 'R_0', 'I_0', 'LogLik'])
    }, 0)
  }
  downloadButton.onclick = function () {
    Csv()
  }
}
function activateRadioButton (radioBtn, i) {
  if (radioBtn.checked) {
    indx[i] = 1
  }
  document.querySelector('button#calc').style.display = ''
  document.querySelector('button#download').style.display = 'none'
}
function deactivateRadioButton (radioBtn, i) {
  if (radioBtn.checked) {
    indx[i] = 0
    document.querySelector('#setup').querySelector('table#table').querySelectorAll('tr')[i + 1].style.backgroundColor = '#FFFFFF'
  }
  document.querySelector('button#calc').style.display = ''
  document.querySelector('button#download').style.display = 'none'
}
function Csv () {
  var csv = ''
  res.forEach(function (row) {
    csv += row.join(',')
    csv += '\n'
  })
  var hiddenElement = document.createElement('a')

  hiddenElement.href = 'data:text/csv;charset=utf-8,' + encodeURI(csv)
  hiddenElement.setAttribute('download', 'result.csv')
  hiddenElement.click()
}
function activateDownload () {
  document.querySelector('button#download').disabled = false
  document.querySelector('button#download').style.backgroundColor = '#2ed573'
}
