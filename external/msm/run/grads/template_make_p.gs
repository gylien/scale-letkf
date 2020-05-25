'sdfopen <--NCFILE_P-->'
'set gxout fwrite'
'set fwrite -be <--GRADSFILE_P-->'
'set x 1 241'
'set y 1 253'
t=1
while (t<=34)
'set t 't
z=1
while (z<=16)
'set z 'z
'd z'
z=z+1
endwhile

z=1
while (z<=16)
'set z 'z
'd w'
z=z+1
endwhile

z=1
while (z<=16)
'set z 'z
'd u'
z=z+1
endwhile

z=1
while (z<=16)
'set z 'z
'd v'
z=z+1
endwhile

z=1
while (z<=16)
'set z 'z
'd temp'
z=z+1
endwhile

z=1
while (z<=16)
'set z 'z
'd rh'
z=z+1
endwhile

t=t+1
endwhile

'disable fwrite'
'quit'
