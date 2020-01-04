function setTime()
*------------------------------
'q dims'
tline = sublin(result, 5)
tstring = subwrd(tline, 6)
tt = subwrd(tline, 9)
tth = (tt-1)*21600
rc = math_format('%06.0f', tth) 
ttf = rc

yr = substr(tstring, 9, 4)
mnf = substr(tstring, 6, 3)
dy = substr(tstring, 4, 2)
hr = substr(tstring, 1, 2)

if (mnf = 'JAN')
  mn = '01'
endif
if (mnf = 'FEB')
  mn = '02'
endif
if (mnf = 'MAR')
  mn = '03'
endif
if (mnf = 'APR')
  mn = '04'
endif
if (mnf = 'MAY')
  mn = '05'
endif
if (mnf = 'JUN')
  mn = '06'
endif
if (mnf = 'JUL')
  mn = '07'
endif
if (mnf = 'AUG')
  mn = '08'
endif
if (mnf = 'SEP')
  mn = '09'
endif
if (mnf = 'OCT')
  mn = '10'
endif
if (mnf = 'NOV')
  mn = '11'
endif
if (mnf = 'DEC')
  mn = '12'
endif

tfilename = 'f'ttf

'set t 1'
'q dims'
tline = sublin(result, 5)
trun = subwrd(tline, 6)
'set t 'tt

*------------------------------
return tt' 'tfilename' 'tstring' 'trun
