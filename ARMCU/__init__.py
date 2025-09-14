#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ARM-Cumulus case definition.

This module implements the ARM-Cumulus (ARMCU) case definition with two
subcases:

- "REF" : Reference case definition.
- "MESONH" : MESO-NH adapted definition.

The class builds initial profiles, forcings, geostrophic wind and surface
fluxes, and provides a method to produce an SCM-ready version of the case.

Original case description:
http://projects.knmi.nl/eurocs/ARM/case_ARM_html/
"""

from copy import deepcopy
import os

import netCDF4 as nc
import numpy as np

from dephycf.Case import Case
from dephycf import thermo
from dephycf import constants as cc

__all__ = ["ARMCU"]

_dirloc = _dirDEPHYCF = os.path.dirname(os.path.abspath(__file__))

class ARMCU(Case):
    """
    ARM-Cumulus case container.

    Parameters
    ----------
    subcase : {'REF', 'MESONH'}
        Subcase identifier.

    Raises
    ------
    ValueError
        If `subcase` is not recognized.
    """

    def __init__(self,subcase):
        if subcase not in ['REF','MESONH']:
           raise ValueError(f'Unknown subcase: {subcase}')

        # --- Step 1: Case metadata ---
        lat = 36.
        lon = -97.5
        startDate = "1997-06-21 11:30:00"
        endDate = "1997-06-22 02:00:00"
        title = "Forcing and initial conditions for ARM-Cumulus case - Original definition"
        authors = "R. Roehrig"
        reference = (
            "http://projects.knmi.nl/eurocs/ARM/case_ARM_html ; "
            "Brown et al. (2002, QJRMS)"
        )
        surfaceType='land'
        zorog = 314.

        if subcase == 'MESONH':
            lat = 36.6
            endDate = "1997-06-22 02:30:00"
            title = "Forcing and initial conditions for ARM-Cumulus case - MESO-NH definition"
            authors = "R. Roehrig, F. Couvreux"

        # Initialize parent Case
        super().__init__(
            f'ARMCU/{subcase}',
            lat=lat,
            lon=lon,
            startDate=startDate,
            endDate=endDate,
            surfaceType=surfaceType,
            zorog=zorog
        )

        # store subcase and metadata
        self._subcase = subcase
        self.set_title(title)
        self.set_reference(reference)
        self.set_author(authors)

        # --- Step 2: initial state
        self._build_initial_state()

        # --- Step 3: atmospheric forcing
        self._build_atm_forcing()

        # --- Step 4: surface forcing
        self._build_sfc_forcing()

        # Default outputs used by SCM()
        # Default new vertical grid, 10-m resolution from surface to 6000 m (above the surface)
        self._levout = np.arange(0.0, 20001.0, 10, dtype=np.float64)
        # Default time frequency for new temporal grid (in seconds)
        self._timefreq = 1800

    # --------------------------
    # Initialization builders
    # --------------------------
    def _build_initial_state(self):
        """Build and register initial profiles depending on the subcase."""

        # Surface pressure (Pa)
        ps = 97000.
        self.add_init_ps(ps)

        if self._subcase == 'REF':
            init_profiles = np.array(
                [         
            #      z (m) theta (K) rt (g kg-1) u (m s-1) v (m s-1)
                    0.0,   299.00,   15.20,      10.0,     0.0,
                   50.0,   301.50,   15.17,      10.0,     0.0,
                  350.0,   302.50,   14.98,      10.0,     0.0,
                  650.0,   303.53,   14.80,      10.0,     0.0,
                  700.0,   303.70,   14.70,      10.0,     0.0,
                 1300.0,   307.13,   13.50,      10.0,     0.0,
                 2500.0,   314.00,    3.00,      10.0,     0.0,
                 5500.0,   343.20,    3.00,      10.0,     0.0,
                ],
                dtype=np.float64,
            )

            z = init_profiles[0::5]
            theta = init_profiles[1::5]
            rt = init_profiles[2::5]/1000.0 # g/kg -> kg/kg
            zwind = z
            u = init_profiles[3::5]
            v = init_profiles[4::5]

            # Register
            self.add_init_theta(theta, lev=z, levtype='altitude')
            self.add_init_rt(rt, lev=z, levtype='altitude')
            self.add_init_wind(u=u, v=v, lev=zwind, levtype='altitude')

            # Turbulent Kinetic Energy (TKE)
            ztke = [0.0, 150.0, 5500.0]
            tke = np.zeros(len(ztke),dtype=np.float64)
            for iz, zz in enumerate(ztke):
                if zz < 150.0:
                  tke[iz] = 0.15 * (1.0 - zz / 150.0)
                else:
                  tke[iz] = 0.0
            self.add_init_tke(tke, lev=ztke, levtype='altitude')

        elif self._subcase == 'MESONH':
            init_profiles = np.array(
                [
            #     z (m), theta (K), rv (kg kg-1)
                    0.0,    299.00,   15.20e-3,
                   50.0,    301.50,   15.17e-3,
                  350.0,    302.50,   14.98e-3,
                  650.0,    303.53,   14.80e-3,
                  700.0,    303.70,   14.70e-3,
                 1300.0,    307.13,   13.50e-3,
                 2500.0,    314.00,    3.00e-3,
                 5500.0,    343.20,    3.00e-3,
                ],
                dtype=np.float64,
            )

            z = init_profiles[0::3]
            theta = init_profiles[1::3]
            rv = init_profiles[2::3]

            # wind profile: z, u, v
            init_profiles = np.array(
                [  
            #     z (m), u (m s-1), v (m s-1)
                    0.0,      10.0,     0.0,
                 5500.0,      10.0,     0.0,
                ],
                dtype=np.float64,
            )

            zwind = init_profiles[0::3]
            u = init_profiles[1::3]
            v = init_profiles[2::3]

            self.add_init_theta(theta, lev=z, levtype='altitude')
            self.add_init_rv(rv, lev=z, levtype='altitude')
            self.add_init_wind(u=u,v=v, lev=zwind, levtype='altitude')

        else:
            raise ValueError(f'Subcase {subcase} is not coded')

        # Keep wind information as attributes for use in forcing
        self._init_zwind = zwind
        self._init_u = u
        self._init_v = v

    # --------------------------
    # Atmospheric forcing builder
    # --------------------------
    def _build_atm_forcing(self):
        """Build geostrophic wind and advections."""

        # Constant Geostrophic wind across the simulation
        if self._subcase == 'REF':
            zforc = self._init_zwind
            ug = self._init_u
            vg = self._init_v
            timeF = None
        elif self._subcase == 'MESONH':
            zforc = np.array([0.0, 1000.0, 3000.0, 5000.0], dtype=np.float64)
            timeF = np.array([41400.0, 52200.0, 63000.0, 73800.0, 84600.0, 86400.0 + 9000.0], dtype=np.float64) - 41400.
            ntf = timeF.size
            nzf = zforc.size
            ug = np.zeros((ntf,nzf),dtype=np.float64) + 10.0
            vg = np.zeros((ntf,nzf),dtype=np.float64) + 0.0
        else:
            raise ValueError(f'Subcase {subcase} is not coded')
        
        # register geostrophic wind forcing
        self.add_geostrophic_wind(ug=ug, vg=vg, lev=zforc, levtype='altitude', time=timeF)

        # --- Advections and radiative tendency ---
        if self._subcase == 'REF':
            #       t (s), A_theta (K hour-1) R_theta (K hour-1) A_rt (g kg-1 hour-1)
            advF = np.array(
                [
            #       t (s), A_theta (K hour-1), R_theta (K hour-1), A_rt (g kg-1 hour-1)
                    41400,              0.000,             -0.125,               0.080,
                    52200,              0.000,              0.000,               0.020,
                    63000,              0.000,              0.000,              -0.040,
                    73800,             -0.080,              0.000,              -0.100,
                    84600,             -0.160,              0.000,              -0.160,
                    93600,             -0.160,             -0.100,              -0.300,
                ],
                dtype=np.float64,
            )

            timeF = advF[0::4] - 41400 # Forcing time axis should be counted from starting date
            ntf = timeF.size
            A_theta = advF[1::4]
            R_theta = advF[2::4]
            A_rt = advF[3::4]

            zforc = np.array([0.0, 1000.0, 3000.0, 5500.0], dtype=np.float64)
            nzf = zforc.size
            forc_theta = np.zeros((ntf,nzf), dtype=np.float64)
            forc_rt = np.zeros((ntf,nzf), dtype=np.float64)

            # 2000 (Brown et al. 2002, QJRMS) ou 3000 m (http://projects.knmi.nl/eurocs/ARM/case_ARM_html/ and used in Meso-NH)
            # We take 3000 m here.

            for it in range(ntf):
                for iz in range(nzf):
                    zlev = zforc[iz]
                    if zlev <= 1000.0:
                        forc_theta[it,iz] = A_theta[it] + R_theta[it]
                        forc_rt[it,iz] = A_rt[it]
                    elif zlev <= 3000.0 :
                        factor = 1.0 - (zlev - 1000.0) / 2000.0
                        forc_theta[it,iz] = (A_theta[it] + R_theta[it]) * factor
                        forc_rt[it,iz] = A_rt[it] * factor
                    else:
                        forc_theta[it,iz] = 0.
                        forc_rt[it,iz] = 0.

            # convert to per second and SI
            adv_theta = forc_theta / 3600.0 # K s-1
            adv_rt = forc_rt / 3600.0 / 1000.0 # kg kg-1 s-1

            # register advective forcings
            self.add_theta_advection(adv_theta, time=timeF, lev=zforc, levtype='altitude', include_rad=True)
            self.add_rt_advection(adv_rt, time=timeF, lev=zforc, levtype='altitude')


        elif self._subcase == 'MESONH':
            # Potential temperature advection (in K s-1)
            tmp_theta = np.array(
                [
                    41400.,       -3.47222E-05, -3.47222E-05, 0.,    0.,
                    52200.,        0.,           0.,          0.,    0.,
                    63000.,        0.,           0.,          0.,    0.,
                    73800.,       -2.22222E-05, -2.22222E-05, 0.,    0.,
                    84600.,       -4.44444E-05, -4.44444E-05, 0.,    0.,
                    86400.+9000., -7.77777E-05, -7.77777E-05, 0.,    0.,
                ],
                dtype=np.float64,
            )


            timeF = tmp_theta[0::5] - 41400.0 # Forcing time axis should be counted from starting date
            zadv = [0.0, 1000.0, 3000.0, 5000.0]
            ntf = timeF.size
            nzf = len(zadv)
            adv_theta = np.zeros((ntf,nzf), dtype=np.float64)
            for it in range(ntf):
                adv_theta[it,:] = tmp_theta[5 * it + 1 : 5 * it + 5]

            # Water vapor mixing ratio advection (in kg kg-1 s-1)
            tmp_rv = np.array(
                [
                    41400.,        2.22222E-08,  2.22222E-08, 0.,    0.,
                    52200.,        5.55555E-09,  5.55555E-09, 0.,    0.,
                    63000.,       -1.11111E-08, -1.11111E-08, 0.,    0.,
                    73800.,       -2.77778E-08, -2.77778E-08, 0.,    0.,
                    84600.,       -4.44444E-08, -4.44444E-08, 0.,    0.,
                    86400.+9000., -9.11111E-08, -9.11111E-08, 0.,    0.,
                ],
                dtype=np.float64,
            )

            adv_rv = np.zeros((ntf,nzf), dtype=np.float64)
            for it in range(ntf):
                adv_rv[it,:] = tmp_rv[5 * it + 1 : 5 * it + 5]

            # register advective forcings
            self.add_theta_advection(adv_theta, time=timeF, lev=zadv, levtype='altitude', include_rad=True)
            self.add_rv_advection(adv_rv, time=timeF, lev=zadv, levtype='altitude')
        
        else:
            raise ValueError(f'Subcase {subcase} is not coded')

    # --------------------------
    # Surface forcing builder
    # --------------------------
    def _build_sfc_forcing(self):
        """Build surface forcing."""

        if self._subcase == 'REF':
            # Surface Forcing
            #            t (s) H (W m-2) LE (W m-2)
            sfc_forc= np.array(
                [
                    41400,  -30,       5,\
                    55800,   90,     250,\
                    64800,  140,     450,\
                    68400,  140,     500,\
                    77400,  100,     420,\
                    86400,  -10,     180,\
                    93600,  -10,       0,
                ],
                dtype=np.float64,
            )

            time_sfc = sfc_forc[0::3] - 41400.0 # Forcing time axis should be counted from starting date
            sens = sfc_forc[1::3]
            lat = sfc_forc[2::3]

        elif self._subcase == 'MESONH':
            txt_path = os.path.join(_dirloc, "MESONH", "surface_flux_forcings.txt")
            time_sfc, evap, sens = np.genfromtxt(txt_path, dtype=float, skip_header=0, usecols=[0,1,2]).transpose()
            lat = evap * cc.Lv

        else:
            raise ValueError(f'Subcase {subcase} is not coded')

        # register surface forcings
        self.add_surface_fluxes(sens=sens, lat=lat, time=time_sfc, forc_wind='z0', z0=0.035)


    def SCM(self, levout=None, timefreq=None):

        if levout is None:
            levout = self._levout
        if timefreq is None:
            timefreq = self._timefreq

        tmin = 0.0
        tmax = float(self.tend - self.tstart)
        timeout = np.arange(tmin, tmax + timefreq, timefreq)

        case_scm = deepcopy(self)

        # Constantly extend winds up to hmax 
        hmax = 20000
        case_scm.extend_init_wind(height=hmax)

        # Extend moisture profile to zero above defined levels
        if self._subcase == 'REF':
            case_scm.extend_init_rt(rt=[0.0, 0.0], height=[5600.0, hmax])
        if self._subcase == 'MESONH':
            case_scm.extend_init_rv(rv=[0.0, 0.0], height=[5600.0, hmax])

        # ERA5-based temperature extension above 12 km
        # 12 km allows us to have a nicer and more stable extrapolated profile.
        era5_path = os.path.join(_dirloc, "aux", "ERA5", "ERA5_SGP_19970621000000-19970622230000.nc")
        with nc.Dataset(era5_path) as f:
            temp = np.average(np.squeeze(f["ta"][:,::-1]), axis=0)
            pa = f["plev"][::-1]
            theta = thermo.t2theta(p=pa, temp=temp)
            height = np.average(np.squeeze(f["zg"][:,::-1]), axis=0)
            mask = height >= 12000
            theta_sel = theta[mask]
            height_sel = height[mask]
            case_scm.extend_init_theta(theta=theta_sel, height=height_sel)

        # Add a surface temperature based on SGP observations
        tskin_path = os.path.join(_dirloc, "aux", "tskin", "tskin_SGP_C1_irt10m_19970621003000-19970622233000.nc")
        with nc.Dataset(tskin_path) as f:
            dates = nc.num2date(f["time"][:], units=f["time"].units)
            index = [ (d.day == 21 and d.hour >= 11) or (d.day == 22 and d.hour <= 2) for d in dates ]
            times = f["time"][index]
            tskin = f["tskin"][index]

            case_scm.add_surface_skin_temp(tskin, time=times - times[0])

        # convert to SCM grid
        case_scm = case_scm.convert2SCM(time=timeout, lev=levout, levtype='altitude')

        # update a few attributes
        case_scm.set_title("Forcing and initial conditions for ARM-Cumulus case - SCM-enabled version")

        return case_scm


# --------------------------
# Command-line usage
# --------------------------
if __name__ == "__main__":
    lverbose = True
    lplot = True
    lcompare = True
    lwrite = True

    case_def = ARMCU('REF')
    case_scm = case_def.SCM()

    if lverbose:
        case_def.info()
        case_scm.info()

    if lplot:
        case_def.plot(rep_images='./test_REF/driver_DEF/',timeunits='hours')
        case_scm.plot(rep_images='./test_REF/driver_SCM/',timeunits='hours')

    if lcompare:
        case_scm.plot_compare(case_def,
            rep_images='./test_REF/compare/',
            label1="SCM-enabled",
            label2="Original",
            timeunits='hours')

    if lwrite:
        case_def.write('ARMCU_REF_DEF_driver.nc')
        case_scm.write('ARMCU_REF_SCM_driver.nc')

    case_def = ARMCU('MESONH')
    case_scm = case_def.SCM()

    if lverbose:
        case_def.info()
        case_scm.info()

    if lplot:
        case_def.plot(rep_images='./test_MESONH/driver_DEF/',timeunits='hours')
        case_scm.plot(rep_images='./test_MESONH/driver_SCM/',timeunits='hours')

    if lcompare:
        case_scm.plot_compare(case_def,
            rep_images='./test_MESONH/compare/',
            label1="SCM-enabled",
            label2="Original",
            timeunits='hours')

    if lwrite:
        case_def.write('ARMCU_MESONH_DEF_driver.nc')
        case_scm.write('ARMCU_MESONH_SCM_driver.nc')
