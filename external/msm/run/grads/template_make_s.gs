'sdfopen <--NCFILE_S-->'
'set gxout fwrite'
'set fwrite -be <--GRADSFILE_S-->'
'set x 1 481'
'set y 1 505'
t=1
while (t<=34)
'set t 't
'd psea'
'd sp'
'd u'
'd v'
'd temp'
'd rh'
'd r1h'
t=t+1
endwhile
'disable fwrite'
'quit'
