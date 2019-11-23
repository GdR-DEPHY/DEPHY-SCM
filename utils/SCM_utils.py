import os, sys
#import time

#import netCDF4 as nc
import numpy as np
import math
from scipy import interpolate
import matplotlib.pyplot as plt

kappa = 2./7.
p0 = 100000.
g = 9.80665 #9.81
R = 287.0596736665907 #287 

def rt2qt(rt,units='kg kg-1'):
    if units == 'kg kg-1':
        return rt/(1+rt)
    elif units == 'g kg-1':
        return rt/(1.+rt/1000.)
    else:
        print 'units unknown:', units
        sys.exit()


def advrt2advqt(rt=None,advrt=None,rt_units='kg kg-1'):

    if rt is None:
        print "rt is missing"
        sys.exit()
    if advrt is None:
        print "advrt is missing"
        sys.exit()

    if rt_units == 'kg kg-1':
        return advrt/((1+rt)*(1+rt))
    elif rt_units == 'g kg-1':
        return advrt/((1.+rt/1000.)*(1.+rt/1000.))
    else:
        print 'units unknown for rt:', rt_units
        sys.exit()


def theta2t(p=None,theta=None,p0=p0):

    if theta is None:
        print "theta is missing"
        sys.exit()
    if p is None:
        print "p is missing"
        sys.exit()

    tmp = theta*(p/p0)**kappa
    return tmp

def t2theta(p,t,p0=p0):

    if t is None:
        print "t is missing"
        sys.exit()
    if p is None:
        print "p is missing"
        sys.exit()

    tmp = t*(p0/p)**kappa
    return tmp


def z2p(theta=None,z=None,ps=None,g=g,R=R,p0=p0):

    if theta is None:
        print "theta is missing"
        sys.exit()
    if z is None:
        print "z is missing"
        sys.exit()
    if ps is None:
        print "ps is missing"
        sys.exit()

    nlev, = theta.shape

    integ = 0.

    p = np.zeros((nlev,),dtype=np.float64)
    p[0] = ps

    for ilev in range(1,nlev):
      dz = z[ilev]- z[ilev-1]    
#      print 'dz =', dz
      integ = integ + (g/(R*theta[ilev-1])+g/(R*theta[ilev]))/2*dz
#      print 'integ =', integ
      tmp = ps**kappa-p0**kappa*kappa*integ
      p[ilev] = math.exp(math.log(tmp)/kappa)

    return p

def plot(x,y,x2=None,y2=None,xlim=None,ylim=None,xlabel=None,ylabel=None,title=None,rep_images=None,name=None,label="",label2=""):
    plt.plot(x,y,'k',label=label)
    if not(xlim is None): plt.xlim(xlim)
    if not(ylim is None): plt.ylim(ylim)
    if not(xlabel is None): plt.xlabel(xlabel)
    if not(ylabel is None): plt.ylabel(ylabel)
    if not(title is None): plt.title(title)

    if not(x2 is None) and not(y2 is None):
        plt.plot(x2,y2,'r',label=label2)
        plt.legend(loc='best')

    if rep_images is None:
      plt.savefig(name)
    else:
      plt.savefig(rep_images + name)
    plt.close()

def plot2D(x,y,z,xlim=None,ylim=None,xlabel=None,ylabel=None,title=None,rep_images=None,name=None):
    nt,nz = z.shape
    X = np.tile(x,(nz,1))
    Y = np.tile(y,(nt,1))
    plt.pcolormesh(X,np.transpose(Y),np.transpose(z))
    if not(xlim is None): plt.xlim(xlim)
    if not(ylim is None): plt.ylim(ylim)
    if not(xlabel is None): plt.xlabel(xlabel)
    if not(ylabel is None): plt.ylabel(ylabel)
    if not(title is None): plt.title(title)
    plt.colorbar()
    if rep_images is None:
      plt.savefig(name)
    else:
      plt.savefig(rep_images + name)
    plt.close()


class Axis:

    def __init__(self,axisid,data,name=None,units=None,**kwargs):

        self.id = axisid
        self.data = np.array(data,dtype=np.float64)
        self.length, = self.data.shape
        self.units = units
        self.name = name
        for x in kwargs.keys():
            self.__dict__[x] = kwargs[x]

    def info(self,data=False):
        print '-'*10, 'axis id:', self.id
        print '-'*10, 'name:', self.name
        print '-'*10, 'units:', self.units
        print '-'*10, 'length:', self.length
        for x in self.__dict__.keys():
            if not(x in ['id','name','units','length','data']):
                print '-'*10, '{0}: {1}'.format(x,self.__dict__[x])
        if data:
            print '-'*10, 'data:', self.data


    def write(self,filein):

        if not(self.id in filein.dimensions):
            dim = filein.createDimension(self.id, self.length)
            ax = filein.createVariable(self.id,"f8",(self.id,))
            ax[:] = self.data
            ax.long_name = self.name
            ax.units = self.units
            for x in self.__dict__.keys():
                if not(x in ['id','units','name','data','length']):
                    ax.setncattr(x,self.__dict__[x])


class Variable:

    def __init__(self,varid,data=None,units=None,name=None,level=None,time=None,lat=None,lon=None,axlist=None,axes=None):

        self.id = varid
        self.units = units
        self.name = name

        self.data = None
        if not(data is None): 
            self.data = np.array(data,dtype=np.float32)
            self.sh = len(self.data.shape)

        if axes is None:
            self.axes = []
            self.axlist = []

            self.time = time
            self.level = level
            self.lat = lat
            self.lon = lon

            for ax in ['time','level','lat','lon']:
                if not(self.__dict__[ax] is None):
                    self.axes.append(self.__dict__[ax])
                    self.axlist.append(self.__dict__[ax].id)
        else:
            self.axes = axes
            self.axlist = axlist

            self.time = None
            self.level = None
            self.lat = None
            self.lon = None

            for ax in axes:
                if ax.id == 't0' or ax.id[0:4] == 'time':
                    self.time = ax
                elif ax.id[0:3] == 'lev' or ax.id == 'nlev':
                    self.level = ax
                elif ax.id == 'lat':
                    self.lat = ax
                elif ax.id == 'lon':
                    self.lon = ax
                else:
                    print 'Axis unexpected:', ax.id
                    ax.info()
                    sys.exit()

        self.coord = " ".join(self.axlist)

    def copy(self,varid):

        tmp = Variable(varid,data=self.data,units=self.units,name=self.name,axes=self.axes,axlist=self.axlist)

        return tmp


    def info(self):
        print '-'*10, 'Variable:', self.id
        print '-'*10, 'Name:', self.name
        print '-'*10, 'Units:', self.units
        print '-'*10, 'Axes:', self.axlist
        print '-'*10, 'Coordinates:', self.coord
        

    def set_coordinates(self,*coord):
        self.axlist = tuple(coord)
        self.coord = " ".join(coord)

    def write(self,filein):

        for ax in self.axes:
            ax.write(filein)

        if not(self.data is None):
          tmp = filein.createVariable(self.id,"f4",self.axlist)
          tmp[:] = self.data
          tmp.long_name = self.name
          tmp.units = self.units
          tmp.coordinates = self.coord

    def plot(self,rep_images=None,var2=None,label="",label2=""):

        if not(self.time is None) and not(self.level is None):

            if self.time.length == 1:
                if var2 is None:
                    plot(self.data[0,:,0,0],self.level.data,xlabel='{0} ({1})'.format(self.id,self.units),ylabel='Altitude ({0})'.format(self.level.units),title='{0} ({1})'.format(self.name,self.time.name),rep_images=rep_images,name='{0}.png'.format(self.id))
                else:
                    plot(self.data[0,:,0,0],self.level.data,x2=var2.data[0,:,0,0],y2=var2.level.data,xlabel='{0} ({1})'.format(self.id,self.units),ylabel='Altitude ({0})'.format(self.level.units),title='{0} ({1})'.format(self.name,self.time.name),rep_images=rep_images,name='{0}.png'.format(self.id),label=label,label2=label2)
            else:
                plot2D(self.time.data,self.level.data,self.data[:,:,0,0],xlabel=self.time.units,ylabel='Altitude ({0})'.format(self.level.units),title=self.name,rep_images=rep_images,name='{0}.png'.format(self.id))

        elif not(self.time is None):

            if self.time.length > 1:
                if var2 is None:
                    plot(self.time.data,self.data[:,0,0],ylabel='{0} ({1})'.format(self.id,self.units),xlabel=self.time.units,title=self.name,rep_images=rep_images,name='{0}.png'.format(self.id)) 
                else:
                    plot(self.time.data,self.data[:,0,0],x2=var2.time.data,y2=var2.data[:,0,0],ylabel='{0} ({1})'.format(self.id,self.units),xlabel=self.time.units,title=self.name,rep_images=rep_images,name='{0}.png'.format(self.id),label=label,label2=label2)
            else:
                print 'no plot for variable', self.id

        elif not(self.level is None):

            if var2 is None:
                plot(self.data[:,0,0],self.level.data,xlabel='{0} ({1})'.format(self.id,self.units),ylabel='Altitude ({0})'.format(self.level.units),title=self.name,rep_images=rep_images,name='{0}.png'.format(self.id))
            else:
                plot(self.data[:,0,0],self.level.data,x2=var2.data[:,0,0],y2=var2.level.data,xlabel='{0} ({1})'.format(self.id,self.units),ylabel='Altitude ({0})'.format(self.level.units),title='{0} ({1})'.format(self.name,self.time.name),rep_images=rep_images,name='{0}.png'.format(self.id),label=label,label2=label2)

        else:
            print 'no plot for variable', self.id



def read(name,filein):

    tmp = filein[name]

    axlist = []
    axes = []
    for ax in tmp.dimensions:
        try:
            axes.append(Axis(ax,filein[ax][:],name=filein[ax].long_name,units=filein[ax].units))
        except AttributeError:
            axes.append(Axis(ax,filein[ax][:],units=filein[ax].units))

        axlist.append(ax)

    try:
        varout = Variable(name,data=tmp[:],name=tmp.long_name,units=tmp.units,axes=axes,axlist=axlist)
    except AttributeError:
        varout = Variable(name,data=tmp[:],units=tmp.units,axes=axes,axlist=axlist)

    return varout


def interpol(var,levout=None,timeout=None):

    if not(var.time is None) and not(var.level is None):

        ntin, nlevin = var.time.length, var.level.length    
        if not(levout is None):
            nlevout = levout.length
            tmp = np.zeros((ntin,nlevout,1,1),dtype=np.float64)
            for it in range(0,ntin):
                ff = interpolate.interp1d(var.level.data,var.data[it,:,0,0],bounds_error=False,fill_value="extrapolate")
                tmp[it,:,0,0] = ff(levout.data)
        else:
            tmp = var.data
            nlevout = nlevin
            levout = var.data.level

        if not(timeout is None) and ntin > 1:
            ntout = timeout.length
            tmp2 = np.zeros((ntout,nlevout,1,1),dtype=np.float64)
            for ilev in range(0,nlevout):
                ff = interpolate.interp1d(var.time.data,tmp[:,ilev,0,0],bounds_error=False,fill_value="extrapolate")
                tmp2[:,ilev,0,0] = ff(timeout.data)
        else:
            tmp2 = tmp
            ntout = ntin
            timeout = var.time

        varout = Variable(var.id,name=var.name,units=var.units,data=tmp2,time=timeout,level=levout,lat=var.lat,lon=var.lon)

    elif not(var.time is None):

        ntin = var.time.length
        if not(timeout is None) and ntin > 1:
            ntout = timeout.length
            tmp = np.zeros((ntout,1,1),dtype=np.float64)
            ff = interpolate.interp1d(var.time.data,var.data[:,0,0],bounds_error=False,fill_value="extrapolate")
            tmp[:,0,0] = ff(timeout.data)
        else:
            tmp = var.data
            ntout = ntin
            timeout = var.time

        varout = Variable(var.id,name=var.name,units=var.units,data=tmp,time=timeout,lat=var.lat,lon=var.lon)

    elif not(var.level is None):

        nlevin = var.level.length    
        if not(levout is None):
            nlevout = levout.length
            tmp = np.zeros((nlevout,1,1),dtype=np.float64)
            ff = interpolate.interp1d(var.level.data,var.data[:,0,0],bounds_error=False,fill_value="extrapolate")
            tmp[:,0,0] = ff(levout.data)
        else:
            tmp = var.data
            nlevout = nlevin
            levout = data.level

        varout = Variable(var.id,name=var.name,units=var.units,data=tmp,level=levout,lat=var.lat,lon=var.lon)

    else:
        print 'Weird, time and level are None for var=', var
        sys.exit()

    return varout
 

