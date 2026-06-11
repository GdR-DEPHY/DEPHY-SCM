import math
from numpy import zeros,array
import pymeteo.skewt as skewt
import matplotlib.pyplot as plt

#z_trop   = 10000.0
z_trop   = 12000.0

th_trop  = 343.0
t_trop   = 213.0
th_sfc   = 300.0
prs_sfc  = 100000.0
qv_pbl   = 0.016
#qv_pbl   = 0.014
g      = 9.81
p00=100000.0
rp00   = 1.0/p00
XAVOGADRO   = 6.0221367E+23
XBOLTZ      = 1.380658E-23
XMD    = 28.9644E-3
rd= XAVOGADRO * XBOLTZ / XMD
cp = 7.* rd /2.
rovcp  = rd/cp
#print (cp)
rv     = 461.5
eps    = rd/rv
reps = rv/rd


rd     = 287.0
#cp     = 1015.0
cp     = 1005.7
cpv    = 1870.0
#eps    = rd/rv
#reps = rv/rd
rovcp  = rd/cp
eps    = rd/rv
reps = rv/rd

XMD    = 28.9644E-3
XMV    = 18.0153E-3
XRD    = XAVOGADRO * XBOLTZ / XMD
XRV    = XAVOGADRO * XBOLTZ / XMV
XCPD   = 7.* XRD /2.
XCPV   = 4.* XRV

#print (XCPD)
#print (XCPV)

XCL    = 4.218E+3
XCI    = 2.106E+3
XTT    = 273.16
XLVTT  = 2.5008E+6
XLSTT  = 2.8345E+6
XLMTT  = XLSTT - XLVTT
XESTT  = 611.14
XGAMW  = (XCL - XCPV) / XRV
XBETAW = (XLVTT/XRV) + (XGAMW * XTT)
XALPW  = math.log(XESTT) + (XBETAW /XTT) + (XGAMW *math.log(XTT))
XGAMI  = (XCI - XCPV) / XRV
XBETAI = (XLSTT/XRV) + (XGAMI * XTT)
XALPI  = math.log(XESTT) + (XBETAI /XTT) + (XGAMI *math.log(XTT))
ZEPS   = XMV / XMD


def rslf(pres, temp):
    esl=611.2 * math.exp( 17.67 * ( temp  - 273.15 ) / ( temp  - 29.65 ))
    esl = min( esl , pres*0.5 )
#    print (esl)
    #rv     = 461.5
    #XAVOGADRO   = 6.0221367E+23
    #XBOLTZ      = 1.380658E-23
    #XMD    = 28.9644E-3
    #rd= XAVOGADRO * XBOLTZ / XMD
#    print (rd)
    eps    = rd/rv
#    print (eps)
    rslf= eps*esl/(pres-esl)
    return rslf

def rfoes(pres,temp):
    """"""
    esl=math.exp(XALPW - XBETAW/temp - XGAMW*math.log(temp))
    esl = min( esl , pres*0.5 )
    rfoes= ZEPS * esl /(pres -esl)
    return rfoes

def rfoesi(pres,temp):
    esl=math.exp(XALPI - XBETAI/temp - XGAMI*math.log(temp))
    rfoesi= ZEPS * esl /(pres -esl)
    return rfoesi

def esat(pres,temp):
    esl=math.exp(XALPW - XBETAW/temp - XGAMW*math.log(temp))
    esat = min( esl , pres*0.5 )
    return esat

pi_sfc  = (prs_sfc/p00)**(rd/cp)
#T=th_sfc*pi_sfc
#esl=611.2 * math.exp( 17.67 * ( T  - 273.15 ) / ( T  - 29.65 ))
#p=prs_sfc
#esl = min( esl , p*0.5 )
#rv     = 461.5

qv_sfc  = rslf(prs_sfc,th_sfc*pi_sfc)
#qv_sfc  = rfoes(prs_sfc,th_sfc*pi_sfc)

#print(qv_sfc)
thv_sfc = th_sfc*(1.0+qv_sfc*reps)/(1.0+qv_sfc)


#print (qv_sfc,th_sfc, thv_sfc)
#stop

zh=zeros(4100)
rh0=zeros(4100)
th0=zeros(4100)

for i in range(0,4100):
    zh[i] = i * 5
    if (zh[i] <= z_trop) :
        th0[i]=th_sfc+(th_trop-th_sfc)*((zh[i]/z_trop)**1.25)
        rh0[i]=1.0-0.75*((zh[i]/z_trop)**1.25)
#        rh0[i]=1.0-0.70*((zh[i]/10000)**1.25)
#        rh0[i]=1.0-0.85*((zh[i]/10000)**1.25)
    else :
        th0[i]=th_trop*math.exp((g/(t_trop*cp))*(zh[i]-z_trop))
#        th0[i]=334.24*math.exp((g/(231*cp))*(zh[i]-10000))
#        rh0[i]=0.30
        rh0[i]=0.25
#        rh0[i]=0.15

thv0=zeros(4100)
pi0=zeros(4100)
prs0=zeros(4100)
qv0=zeros(4100)

pi0[0]=pi_sfc
for j in range(0,20):
    for i in range(0,4100):
        thv0[i]=th0[i]*(1.0+reps*qv0[i])/(1.0+qv0[i])
    pi0[1]=pi_sfc-g*zh[1]/(cp*0.5*(thv_sfc+thv0[1]))
    for i in range(2,4100):
        pi0[i]=pi0[i-1]-g*(zh[i]-zh[i-1])/(cp*0.5*(thv0[i]+thv0[i-1]))
    for i in range(0,4100):
        prs0[i]=p00*(pi0[i]**(cp/rd))
    for i in range(0,4100):
        qv0[i]=rh0[i]*rslf(prs0[i],th0[i]*pi0[i])
#        qv0[i]=rh0[i]*rfoes(prs0[i],th0[i]*pi0[i])
        if(qv0[i] > qv_pbl):
            qv0[i]=qv_pbl



rh03=zeros(4100)

rh02=zeros(4100)
qs0=zeros(4100)
lnexp=zeros(4100)
td0=zeros(4100)
u0=zeros(4100)
v0=zeros(4100)



for i in range(0,4100):
    rh02[i]=qv0[i]/(rslf(prs0[i],th0[i]*pi0[i]))
#    rh02[i]=qv0[i]/(rfoes(prs0[i],th0[i]*pi0[i]))
    qs0[i]=qv0[i]/(1+qv0[i])
    if (rh02[i]*math.exp((17.62*(th0[i]*pi0[i]-273.15))/(243.12 - 273.15 + th0[i]*pi0[i])) > 0):
        lnexp[i]= math.log(rh02[i]*math.exp((17.62*(th0[i]*pi0[i]-273.15))/(243.12 - 273.15 + th0[i]*pi0[i])))
        td0[i]= 243.12 * lnexp[i]/(17.62 - lnexp[i])
        rh03[i]=100*esat(prs0[i],td0[i]+273.15)/esat(prs0[i],th0[i]*pi0[i])
        u0[i]=1
        v0[i]=1
psurf = prs_sfc
tsurf = th_sfc*((psurf*rp00)**rovcp)

if( qv0[1] < qv_pbl ):
#qsurf = cgs1*qv0(1)+cgs2*qv0(2)+cgs3*qv0(3)
    print ('pbleme')
else:
    qsurf = qv_pbl
    print ('ok')

#for i in range(0,4100,10):
#    print (zh[i],prs0[i],thv0[i],th0[i],th0[i]*pi0[i],td0[i]+273.15,rh0[i],rh02[i],qv0[i],qs0[i])

#for i in range(0,4100,10):
#    print (zh[i],prs0[i],th0[i]*pi0[i],td0[i]+273.15,rh0[i],rh02[i],rh03[i],qv0[i])

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
#ax.plot(rh02,prs0)
#ax.plot(th0*pi0,prs0)
#ax.plot(td0+273.15,prs0)
ax.plot(th0*pi0,zh)
ax.plot(td0+273.15,zh)


#ax.set_ylim([100000, 2000])
plt.show()

fichier = open("dataRH25.txt", "w")
for i in range(0,4100,20):
    fichier.write(str(prs0[i])+'\t')
    fichier.write(str(th0[i]*pi0[i])+'\t')
    fichier.write(str(td0[i]+273.15)+'\n')

fichier.close()

skewt.plot(None, zh, th0, prs0, qv0, u0, v0, 'RS_ideal', title="00 UTC")
