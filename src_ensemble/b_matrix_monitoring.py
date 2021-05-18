#!/usr/bin/env python

no=3
n=12

fac=[]
for io in range(no):
    fac.append(2**(no-io))
print(fac)

nn=[]
for io in range(no):
    nn.append(n/fac(io))
print(nn)


# do io=1,no
#    allocate(delta(no,nn(io),nn(io))
#    delta(io,:,:)=0.0
#    do ib=1,nn(io)
#       delta(io,ib,:)=0.0
#       do id=1,nn(io)
#          delta(io,nn(io),ib)=1
#       write(*,'(a)') 'delta:',delta   
