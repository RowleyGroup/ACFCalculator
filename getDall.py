import os
import numpy
import sys
import math

lin=numpy.array([])

for i in range(1, 6):
    ds=[]
    d=[]
    x=[]
    p=[]
    temp=numpy.array([])
    for j in range(39, -1, -1):
        fh=open(os.path.expanduser('~/out_up'+str(i)+'/out_up_sim_'+str(j)+'_2.out'), 'r')
        lines=fh.readlines()
        fh.close()
        
        for l in lines:
	    f=l.split();
            if(f[0]=='#avg'):
                x.append(float(f[-1]))
            if(f[0]=='#D'):
                d.append(float(f[-2]))
            if(f[0]=='#Ds'):
                ds.append(float(f[-2]))
		p.append(int(-j))

    for k in range(0, 40):
        fh=open(os.path.expanduser('~/out_down'+str(i)+'/out_down_sim_'+str(k)+'_2.out'), 'r')
        lines=fh.readlines()
        fh.close()

        for l in lines:
            f=l.split();
            if(f[0]=='#avg'):
                x.append(float(f[-1]))
            if(f[0]=='#D'):
                d.append(float(f[-2]))
            if(f[0]=='#Ds'):
                ds.append(float(f[-2]))
		p.append(int(k))
    
    if len(lin)==0:
      lin = [x, d, d, ds, ds]
    else: 
      temp = numpy.dstack((lin, [x, d, d, ds, ds]))
      lin = temp

    f_out = open(os.path.expanduser('~/out_ds'+str(i)+'_3.out'), 'w')
    for l in range(0, len(ds)):
    	f_out.write(str(p[l])+'\t'+str(ds[l])+ '\n')
    f_out.close()

sys.stdout.write('{:16} {:15} {:20} {:15} {:20} \n'.format('x', 'Dacf', '+-','Dvacf', '+-'))

for i in range(0, len(lin[0])):	
     sys.stdout.write('{:16} {:15} {:20} {:15} {:20} \n'.format(str(numpy.average(lin[0][i])), str(numpy.average(lin[1][i])), str(numpy.std(lin[2][i])), str(numpy.average(lin[3][i])), str(numpy.std(lin[4][i]))))    
