program read_test

character*30::ctime

open(11,file='test.txt')
 read(11,'(A14,E11.3,E11.3,I)') ctime,rdata,rdata2,inum
! read(11,*) rdatam
 write(*,*) rdata
close(11)


stop
end program read_test
