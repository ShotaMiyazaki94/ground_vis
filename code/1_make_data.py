from oem import OrbitEphemerisMessage
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import cm
from astropy.coordinates import CartesianRepresentation, GCRS
from astropy.coordinates import (SkyCoord, Distance, Galactic, EarthLocation, AltAz, ICRS, TEME)
from astropy.coordinates import get_body_barycentric
import astropy.units as u
from astropy.time import Time
import seaborn as sns
sns.set(font="times",font_scale=1.6,style="ticks")
from scipy.interpolate import interp1d

#oem = OrbitEphemerisMessage.open("data/RST_032427.oem") 
oem = OrbitEphemerisMessage.open("../misasa2/v2/data/RST_Ephemeris_Files_From_NASA_Flight_Dynamics/RST_EPH_PRED_LONG_2027047_2032046_01.oem") 
pos_r  = []
time_r = []
vel_r  = []
for state in oem.states:
    time_r.append(state.epoch)
    pos_r.append(state.position)
    vel_r.append(state.velocity)
X,Y,Z    = np.array(pos_r).T
Vx,Vy,Vz = np.array(vel_r).T
time_r   = Time(time_r)
time_jd  = time_r.jd

X_m = interp1d(time_jd,X)
Y_m = interp1d(time_jd,Y)
Z_m = interp1d(time_jd,Z)
time_use = Time("2027-02-17 00:00:00") + u.year*np.arange(0,4.99,1.0/365.25/24./6) #10min bin
X_use    = X_m(time_use.jd) 
Y_use    = Y_m(time_use.jd) 
Z_use    = Z_m(time_use.jd)

orb_roman = SkyCoord(x=X_use,y=Y_use,z=Z_use, representation_type='cartesian', frame=GCRS(obstime=time_use), unit=u.km)
orb_roman  = orb_roman.transform_to("icrs")
orb_Earth  = SkyCoord(get_body_barycentric('earth', time_use), frame='icrs')

Misasa      = EarthLocation(lat=36.13*u.deg,  lon=138.35*u.deg, height=1612.75*u.m)
Whitesands  = EarthLocation(lat=32.507*u.deg, lon=-106.611*u.deg, height=1192*u.m)
NNS         = EarthLocation(lat=-(31+2/60+54/3600)*u.deg, lon=(116+11/60+28/3600)*u.deg, height=252*u.m)
coord_roman = SkyCoord(ra=orb_roman.ra, dec=orb_roman.dec, unit=u.deg, frame="icrs")

altaz_MDSS_roman = coord_roman.transform_to(AltAz(obstime=time_use,location=Misasa))
altaz_WS_roman   = coord_roman.transform_to(AltAz(obstime=time_use,location=Whitesands))
altaz_NNS_roman  = coord_roman.transform_to(AltAz(obstime=time_use,location=NNS))

jwst = pd.read_csv("../misasa2/v2/data/jwst_20220101-20250101_1hour_xyz.txt",comment="#",header=None,
                   names=["jd","calender","x","y","z","vx","vy","vz","w"])
time_j = Time(jwst.jd,format="jd")
orb_jwst = SkyCoord(x=jwst.x,y=jwst.y,z=jwst.z,representation_type='cartesian', frame='icrs',unit=u.km)
orb_jwst.representation_type = 'spherical'
coord_jwst = SkyCoord(ra=orb_jwst.ra,dec=orb_jwst.dec,unit=u.deg)
altaz_MDSS_jwst = coord_jwst.transform_to(AltAz(obstime=time_j,location=Misasa))
altaz_WS_jwst   = coord_jwst.transform_to(AltAz(obstime=time_j,location=Whitesands))
altaz_NNS_jwst  = coord_jwst.transform_to(AltAz(obstime=time_j,location=NNS))

roman_from_Earth = SkyCoord(x=orb_roman.cartesian.x - orb_Earth.cartesian.x,
                            y=orb_roman.cartesian.y - orb_Earth.cartesian.y,
                            z=orb_roman.cartesian.z - orb_Earth.cartesian.z, representation_type='cartesian', frame='icrs')
distance = np.linalg.norm(np.array([roman_from_Earth.x.value, roman_from_Earth.y.value, roman_from_Earth.z.value]).T,axis=1)
distance = distance - 6378.0
Roman_output = np.array([time_use.jd, altaz_MDSS_roman.alt.to_value(), distance])
np.save("data/Roman_time_alt_distance",Roman_output)

fig,ax=plt.subplots(2,1,figsize=(12,10),sharey=True)
fig.subplots_adjust(hspace=0.25)

ax[0].set_title("Roman")

ax[0].plot(time_use.decimalyear,altaz_NNS_roman.alt.to_value(),lw=1,label="New Norcia Station",alpha=0.5)
ax[0].plot(time_use.decimalyear,altaz_WS_roman.alt.to_value(),lw=1,label="White Sands Space Harbor",alpha=0.5)
ax[0].plot(time_use.decimalyear,altaz_MDSS_roman.alt.to_value(),lw=1,label="Misasa Deep Space Station",alpha=0.5)

ax[0].hlines(10,2028,2031,ls="--",color="k",lw=1)
ax[0].hlines(15,2028,2031,ls=":",color="k",lw=1)
ax[0].set_xlim(2027,2032.2)
ax[0].set_ylim(-90,90)
ax[0].legend(loc="upper left",fontsize=12)
ax[0].set_ylabel("elevation (deg)")

ax[1].set_title("JWST")
ax[1].plot(time_j.decimalyear,altaz_NNS_jwst.alt.to_value(),lw=1,label="New Norcia Station",alpha=0.5)
ax[1].plot(time_j.decimalyear,altaz_WS_jwst.alt.to_value(),lw=1,label="White Sands Space Harbor",alpha=0.5)
ax[1].plot(time_j.decimalyear,altaz_MDSS_jwst.alt.to_value(),lw=1,label="Misasa Deep Space Station",alpha=0.5)

ax[1].hlines(10,2022,2025,ls="--",color="k",lw=1)
ax[1].hlines(15,2022,2025,ls=":",color="k",lw=1)
ax[1].set_xlim(2022,2025)
ax[1].set_ylim(-90,90)
ax[1].legend(loc="upper left",fontsize=12)
ax[1].set_ylabel("elevation (deg)")
ax[1].set_xlabel("year (UT)")
ax[0].set_yticks(np.arange(-90,120,30))
plt.savefig("plots/elevation_roman_jwst.png", dpi=200, bbox_inches="tight")
plt.close()

fig,ax=plt.subplots(2,1,figsize=(10,8),sharey=False, sharex=True)
fig.subplots_adjust(hspace=0.1)

ax[0].set_title("Roman seen from the ground stations")
ax[0].plot(time_use.decimalyear,altaz_MDSS_roman.alt.to_value(),lw=0.7,label="Misasa Deep Space Station")
ax[0].plot(time_use.decimalyear,altaz_WS_roman.alt.to_value(),lw=0.5,label="White Sands Space Harbor")
ax[0].plot(time_use.decimalyear,altaz_NNS_roman.alt.to_value(),lw=0.3,label="New Norcia Station")
ax[0].hlines(10,2028,2031,ls="--",color="k",lw=2,label="10 deg")
ax[0].hlines(15,2028,2031,ls=":",color="k",lw=2,label="15 deg")
ax[0].set_xlim(2028,2029)
ax[0].set_ylim(-90,90)
ax[0].set_yticks(np.arange(-90,120,30))
ax[0].set_xticks(np.arange(2028,2029.1,2/12),["2028/1","2028/3","2028/5","2028/7","2028/9","2028/11","2029/1"])
ax[0].legend(bbox_to_anchor=(0.58,0.48),fontsize=12)
ax[0].set_ylabel("elevation (deg)")
ax[1].plot(time_use.decimalyear,distance,lw=2)
ax[1].set_xlabel("time (UT)")
ax[1].set_xticks(np.arange(2028,2029.1,2/12),["2028/1","2028/3","2028/5","2028/7","2028/9","2028/11","2029/1"])
#ax.legend(bbox_to_anchor=(1,0.6),fontsize=15)
ax[1].set_ylabel("distance (km)")
ax[1].grid(ls="--")
plt.savefig("plots/elevation_roman.png", dpi=200, bbox_inches="tight")
plt.close()




