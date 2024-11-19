#!/usr/bin/python
# -*-coding:utf-8 -*-

import math

"""
Calculs thermodynamiques

----------------------------
ATTENTION !

Unités utilisées :
- hPa  pour les pressions
- K    pour les températures
- %    pour les humidités relatives
- g/kg pour les rapports de mélange, humidités spécifiques, etc.
----------------------------

Création le 23/07/2013 par CNRM/GMEI/4M (DT) à partir de la lib ASPIC thermo.c
Modif. le 03/06/2016 : ajout calcul humidité absolue par CNRM/GMEI/4M (AB)

"""


#=========================================

# Constantes (SI)

T0 = 273.15 # K

Ra = 287.05 # Specific gas constant for dry air (J.kg-1.K-1)
Rw = 461.53    # Specific gas constant for water vapour (J.kg-1.K-1)

Cp = 1004.6 # Specific heat at constant pressure for dry air (J.kg-1.K-1)



#=========================================

def ew(t):
  """ Calcul de la pression de vapeur saturante par rapport à l'eau liquide ew (hPa)
  """
  global T0

  t -= T0 # conversion en °C
  ew =  6.112 * math.exp ( (17.67 * t) / (t + 243.12) )

  return ew


#=========================================

def ewi(t):
  """ Calcul de la pression de vapeur saturante par rapport à la glace ewi (hPa)
      à partir de la température t (K)
  """
  global T0

  t -= T0 # conversion en °C
  ewi =  6.112 * math.exp ( (22.46 * t) / (t + 272.62) )

  return ewi


#=========================================

def rm_from_hu(p,t,hu):
  """ Calcul du rapport de mélange (g/kg)
      à partir de la pression p (hPa)
      de la température t (K)
      et de l'humidité relative hu (%)
  """

  e = ew(t) * hu / 100.
  rm = 0.622 * 1000. * e / (p - e)

  return rm


#=========================================

def rm_from_td(p,td):
  """ Calcul du rapport de mélange (g/kg)
      à partir de la pression p (hPa)
      et de la température du point de rosée td (K)
  """

  e = ew(td)
  rm = 0.622 * 1000. * e / (p - e)

  return rm


#=========================================

def qv_from_hu(p,t,hu):
  """ Calcul de l'humidité spécifique q (g/kg)
      à partir de la pression p (hPa)
      de la température t (K)
      et de l'humidité relative hu (%)
  """

  rm = rm_from_hu(p,t,hu) / 1000.
  qv = rm / (1. + rm)
  qv *= 1000.

  return qv


#=========================================

def qv_from_rm(rm):
  """ Calcul de l'humidité spécifique qv (g/kg)
      à partir du rapport de mélange rm (g/kg)
  """

  rm /= 1000.
  qv = rm / (1. + rm)
  qv *= 1000.

  return qv


#=========================================

def td_from_rm(p,rm):
  """ Calcul de la température du point de rosée td (K)
      à partir de la pression p (hPa)
      et du rapport de mélange (g/kg)
  """
  global T0

  e = p * rm / ( 622. + rm ) # pression de vapeur = ew(td)
  loge= math.log( e / 6.112 )
  td = ( 243.5 * loge ) / ( 17.67 - loge) + T0

  return td


#=========================================

def td_from_hu(t,hu):
  """ Calcul de la température du point de rosée td (K)
      à partir de la température t (K)
      et de l'humidité relative hu (%)
  """
  global T0

  ewtd = hu * ew(t) / 100.
  loge= math.log( ewtd / 6.112)
  td = ( 243.5 * loge ) / ( 17.67 -  loge) + T0

  return td


#=========================================

def hu_from_td(t,td):
  """ Calcul de l'humidité relative hu (%)
      à partir de la température t (K)
      et de la température du point de rosée td (K)
  """
  global T0

  hu = 100. * ew(td) / ew(t)

  return hu


#=========================================

def tcond_from_rm(p,t,rm):
  """ Calcul de la température du point de condensation tcond (K)
      à partir de la pression (hPa)
      de la température t (K)
      et du rapport de mélange (g/kg)

  """

  rmels = rm_from_hu(p,t,100.)
  humi = rm / (622. + rm) / rmels * (622. + rmels)

  a = 1. / ( t - 55. )
  b = math.log(humi) / 2840.
  tcond = 1. / (a - b)  +  55.

  return tcond


#=========================================

def tcond_from_hu(t,hu):
  """ Calcul de la température du point de condensation tcond (K)
      à partir de la température t (K)
      et de l'humidité relative (%)
  """

  a = 1. / ( t - 55. )
  b = math.log(hu / 100.) / 2840.;
  tcond = 1. / (a - b)  +  55.

  return tcond


#=========================================

def pcond_from_hu(p,t,hu):
  """ Calcul de la pression du point de condensation pcond (hPa)
      à partir de la pression (hPa)
      de la température t (K)
      et de l'humidité relative (%)

  """

  e = ew(t) * hu / 100.;
  tc = tcond_from_hu(t,hu);
  pcond = p * ew(tc) / e

  return pcond


#=========================================

def theta_from_t(p,t):
  """ Calcul de la température potentielle theta (K)
      à partir de la pression p (hPa)
      et de la température t (K)
  """
  global Ra
  global Cp

  theta = t * math.pow(1000./p,Ra/Cp)

  return theta


#=========================================

def t_from_theta(p,theta):
  """ Calcul de la température t (K)
      à partir de la pression p (hPa)
      et de la température potentielle theta (K)
  """
  global Ra
  global Cp

  t = theta * math.pow(p/1000.,Ra/Cp)

  return t


#=========================================

def tv_from_hu(p,t,hu):
  """ Calcul de la température virtuelle tv (K)
      à partir de la pression p (hPa)
      de la température t (K)
      et de l'humidité relative hu (%)
  """
  rm = rm_from_hu(p,t,hu) / 1000.
  tv = t * ( 1. + 1.608 * rm ) / (1. + rm)

  return tv


#=========================================

def thetav_from_hu(p,t,hu):
  """ Calcul de la température potentielle virtuelle thetav (K)
      à partir de la pression p (hPa)
      de la température t (K)
      et de l'humidité relative hu (%)
  """
  theta = theta_from_t(p,t)
  qv = qv_from_hu(p,t,hu) / 1000.
  thetav = theta * ( 1. + 0.608 * qv )

  return thetav


#=========================================

def ha_from_hu(hu,t):
    """ Calcul de l'humidité absolue ha (g/m^3)
    à partir de la température t (K)
    et de l'humidité relative hu (%)
    """
    e = ew(t) * hu
    ha=e / (Rw * t) * 1000
    return ha

#=========================================


def ha_from_td(td,t):
    """Calcul de l'humidité absolue ha (g/m^3)
    à partir de la température t (K)
    et de la température du point de rosée td (K)
    """
    e = ew(td)
    ha=e / (Rw * t) * 1000 * 100
    return ha
