import numpy as np 

from pylab import plot, savefig, xlabel,ylabel, loglog, semilogy,semilogx, show,gca,ylim,xlim,fill_between,axhline,legend

relhum_test=np.genfromtxt('clima_relhum_test.txt',unpack=True,skip_header=2,skip_footer=1)
basic_test=np.genfromtxt('clima_basic_test.txt',unpack=True,skip_header=2,skip_footer=1)
alt_relhum_test=np.genfromtxt('alt_relhum_test.txt',unpack=True,skip_header=1)

pres=relhum_test[1]
alt=relhum_test[2]

temper=relhum_test[3]




pres_basic=basic_test[1]
alt_basic=basic_test[2]
temper_basic=basic_test[3]

conv_basic=basic_test[4]

cloud=alt_relhum_test[2]

conv=relhum_test[4]
iconv=np.where(conv >= 1.)
iconv_basic=np.where(conv_basic >= 1.)


icloud=np.where(cloud >= 1.e-3)



print(icloud)
print()
print(iconv)
plot(temper,alt,color='blue',label='Cloudy')

for i in range(len(alt_basic[iconv_basic])):
        axhline(alt_basic[iconv_basic][i],-10,10,color='red')
for i in range(len(alt[iconv])):

    axhline(alt[iconv][i],-10,10,color='blue',ls='-.')



plot(temper_basic,alt_basic,color='red',label='Clear')
#gca().invert_yaxis()
xlabel('Temperature (K)')
ylabel('Altitude (km)')
legend(loc='best')
show()

plot(temper,alt,color='blue')

for i in range(len(alt_basic[iconv_basic])):
        axhline(alt_basic[iconv_basic][i],-10,10,color='red')
for i in range(len(alt[iconv])):

    axhline(alt[iconv][i],-10,10,color='blue',ls='-.')



plot(temper_basic,alt_basic,color='red')
#gca().invert_yaxis()
ylim(0,20)
xlabel('Temperature (K)')
ylabel('Altitude (km)')
show()

def dt(T1,T2,P1,P2):
    print( T2-T1)
    return (T2-T1)/(P2-P1)
    
    

deltat=np.zeros(len(temper))
for i in range(len(temper)-1):
    deltat[i]=deltat[i] + dt(temper[i],temper[i+1],alt[i],alt[1+i])
    print(i,i+1)



deltat_basic=np.zeros(len(temper_basic))
for k in range(len(temper_basic)-1):
    deltat_basic[k]=deltat_basic[k] + dt(temper_basic[k],temper_basic[k+1],alt_basic[k],alt_basic[1+k])

print(deltat_basic)
xdel=np.linspace(-12,12,len(pres[iconv]))
# semilogy(deltat,pres,color='blue')
# semilogy(deltat_basic,pres_basic,color='red')

# #plot(deltat,alt)
# ylim(1e-4,1.2)
# #fill_between(xdel,np.max(pres[iconv]),np.min(pres[iconv]),alpha=0.1)
# for i in range(len(pres_basic[iconv_basic])):

#     axhline(pres_basic[iconv_basic][i],-12,12,color='red')
# for i in range(len(pres[iconv])):

#     axhline(pres[iconv][i],-12,12,color='blue',ls='-.')


# #xlim(-1e4,1e4)
# gca().invert_yaxis()
# show()

plot(deltat,alt,color='blue',label='dt/dz cloudy')
plot(deltat_basic,alt,color='red',label='dt/zd clear')
for i in range(len(alt_basic[iconv_basic])):

    axhline(alt_basic[iconv_basic][i],-12,12,color='red')
for i in range(len(alt[iconv])):

    axhline(alt[iconv][i],-12,12,color='blue',ls='-.')
#for i in range(len(alt[icloud])):
#    axhline(alt[icloud][i],-12,12,color='red',ls='-.')

xlabel('dt/dz')
ylabel('alt (km)')
fill_between(xdel,np.max(alt[icloud]),np.min(alt[icloud]),alpha=0.1,color='black')
xlim(-12,12)
legend(loc='best')
show()


plot(deltat,alt,color='blue',label='dt/dz cloudy')
plot(deltat_basic,alt,color='red',label=' dt/dz clear')
for i in range(len(alt_basic[iconv_basic])):

    axhline(alt_basic[iconv_basic][i],-12,12,color='red')
for i in range(len(alt[iconv])):

    axhline(alt[iconv][i],-12,12,color='blue',ls='-.')
#for i in range(len(alt[icloud])):
#    axhline(alt[icloud][i],-12,12,color='red',ls='-.')

xlabel('dt/dz')
ylabel('alt (km)')
fill_between(xdel,np.max(alt[icloud]),np.min(alt[icloud]),alpha=0.1,color='black')
xlim(-12,12)
ylim(0,20)
legend(loc='best')
show()


