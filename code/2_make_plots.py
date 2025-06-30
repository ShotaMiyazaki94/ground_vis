import itur
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(font="times",font_scale=1.5,style="ticks")
plt.rcParams['text.usetex'] = True
import pandas as pd
from scipy.interpolate import interp1d
from astropy.time import Time

el = np.linspace(5, 90, 100) # Vectorize the elevation angle
f = 26.675 * u.GHz            # Link frequency
p = 5                         # Unavailability (Values exceeded 1% of time)

lat_GS,lon_GS = 36.13,138.35  # Ground station coordinates (Misasa)
#D = 1 * u.m                  # Antenna diameters for Misasa
D = 54 * u.m                  # Antenna diameters for Misasa
Att_MDSS = itur.atmospheric_attenuation_slant_path(lat_GS, lon_GS, f, el, p, D)
Att_cont_MDSS = itur.atmospheric_attenuation_slant_path(lat_GS, lon_GS, f, el, p, D,return_contributions=True)

lat_GS,lon_GS =32.943242, -106.419531 # Ground station coordinates (White Sands)
D = 18 * u.m       # Antenna diameters for WS
Att_WS = itur.atmospheric_attenuation_slant_path(lat_GS, lon_GS, f, el, p, D)
Att_cont_WS = itur.atmospheric_attenuation_slant_path(lat_GS, lon_GS, f, el, p, D,return_contributions=True)

lat_GS,lon_GS =-(31+2/60+54/3600), (116+11/60+28/3600) # Ground station coordinates (NNS)
D = 35 * u.m       # Antenna diameters for NNS
Att_NNS = itur.atmospheric_attenuation_slant_path(lat_GS, lon_GS, f, el, p, D)
Att_cont_NNS = itur.atmospheric_attenuation_slant_path(lat_GS, lon_GS, f, el, p, D,return_contributions=True)

fig,ax = plt.subplots(1,3,figsize=(15,5),sharex=True,sharey=True)
fig.subplots_adjust(wspace=0.1)
labels=["gas","cloud","rain","scintillation","total"]
colors=["red","blue","orange","green","k"]
for i in range(5):
    ax[0].plot(el, Att_cont_MDSS[i].value,lw=3,label=labels[i]+" (MDSS)",color=colors[i],ls="-",alpha=0.8)
    ax[1].plot(el, Att_cont_WS[i].value,lw=3,label=labels[i]+" (WSSH)",color=colors[i],ls="-",alpha=0.8)
    ax[2].plot(el, Att_cont_NNS[i].value,lw=3,label=labels[i]+" (NNS)",color=colors[i],ls="-",alpha=0.8)
ax[0].set_title(r"MDSS (%s, D=54m, %s"%(f,100-p)+"\\%"+" availability)",fontsize=13)
ax[1].set_title(r"White Sands (%s, D=18m, %s"%(f,100-p)+"\\%"+" availability)",fontsize=13)
ax[2].set_title(r"New Norcia (%s, D=35m, %s"%(f,100-p)+"\\%"+" availability)",fontsize=13)
ax[0].legend(fontsize=14)
ax[1].legend(fontsize=14)
ax[2].legend(fontsize=14)
ax[0].set_xticks(np.arange(0,100,10))
ax[0].set_yticks(np.arange(0,25,1))
ax[0].set_ylabel("attenuation (dB)")
ax[0].set_xlabel("elevation (deg)")
ax[1].set_xlabel("elevation (deg)")
ax[2].set_xlabel("elevation (deg)")
ax[0].set_xlim(0,90)
ax[0].set_ylim(0,20)
ax[0].grid(ls=":")
ax[1].grid(ls=":")
ax[2].grid(ls=":")
plt.savefig("plots/attenuation_each_site.png", dpi=200, bbox_inches="tight")
plt.close()

plt.figure(figsize=(6,7))

labels=["gas","cloud","rain","scintillation","total"]
colors=["red","blue","orange","green","k"]
ls    =["-",":","--","-."]

for i in range(5):
    plt.plot(el, Att_cont_MDSS[i].value,lw=3,label=labels[i]+" (MDSS)",color=colors[i],ls="-",alpha=0.8)
    plt.plot(el, Att_cont_WS[i].value,lw=2.5,label=labels[i]+" (White Sands)",color=colors[i],ls="-.",alpha=0.8)
    plt.plot(el, Att_cont_NNS[i].value,lw=2,label=labels[i]+" (New Norcia)",color=colors[i],ls=":",alpha=1)
plt.legend(loc="upper left",bbox_to_anchor=(1.03,1))
plt.title(r"ITU-R P.618 (%s,D=%s, %s"%(f,D,100-p)+"\\%"+" availability)",fontsize=15)
plt.xlabel('elevation angle [deg]')
plt.ylabel('attenuation [dB]')
plt.grid(which='major', linestyle=':',alpha=.8)
plt.xticks(np.arange(0,100,10))
plt.yticks(np.arange(0,20,1))
#plt.yscale("log")
plt.xlim(0,90)
plt.savefig("plots/attenuation_each_site_comp.png", dpi=200, bbox_inches="tight")
plt.close()

Roman = np.load("data/Roman_time_alt_distance.npy")
time  = Time(Roman[0],format="jd") 
lat_GS,lon_GS = 36.13,138.35  # Ground station coordinates (Misasa)
D = 54 * u.m                  # Antenna diameters for Misasa
el_Roman = Roman[1]
el_Roman = np.where(el_Roman>15,el_Roman,0)
Att = itur.atmospheric_attenuation_slant_path(lat_GS, lon_GS, f, el_Roman, p, D)
Att2 = np.where(np.isinf(Att.value),np.nan,Att.value)

Att_m = interp1d(time.decimalyear,Att2)
time_plot = np.arange(2028,2029,1/365/24/60)
Att_plot  = Att_m(time_plot)

from matplotlib import colors
fig,ax=plt.subplots(figsize=(8,5))
plt.hlines(6.24,np.min(time.decimalyear),np.max(time.decimalyear),zorder=3,color="cyan",ls="-",lw=3)
plt.hlines(7.2,np.min(time.decimalyear),np.max(time.decimalyear),zorder=3,color="magenta",ls="--",lw=3)
plt.hist2d(time_plot,np.nan_to_num(Att_plot,nan=1e-5),bins=1000,
           norm=colors.LogNorm(vmin=1),cmap="coolwarm_r")
cb=plt.colorbar(pad=0.01)
cb.set_label("observable time (min)")
cb.set_ticks([1,10,100,1000],labels=["1","10","100","1000"])
#plt.scatter(time_plot,Att_plot,color="gray",alpha=0.01,s=20)
#plt.plot(time_plot,Att_plot)
plt.xlim(2028,2029)
plt.ylim(0.5,7.5)
plt.ylabel("attenuation (dB)")
plt.xlabel("year (UT)")
plt.title("%s, %s, %s"%(f,D,100-p)+"\%"+" availability, $\\theta>15$ deg")
plt.savefig("plots/attenuation_misasa_year.png", dpi=200, bbox_inches="tight")
plt.close()

GHz = 26.675
FSL_r  = 92.45 + 20*np.log10(Roman[2]) + 20*np.log10(GHz)  
plt.subplots(figsize=(8,6))
plt.hlines(246.1,np.min(time.decimalyear),np.max(time.decimalyear),zorder=3,color="red",ls="-",lw=3,label="ICD value (246.1 dB)")
plt.plot(time.decimalyear,FSL_r,label="Expected FSL from Roman orbit",lw=3)
plt.ylabel("Free Space Loss (dB)")
plt.xlabel("year (UT)")
plt.title("Expected Free Space Loss")
plt.xlim(2028,2032)
plt.grid(alpha=0.5)
plt.legend(loc="lower center")
plt.savefig("plots/free_space_loss_misasa_year.png", dpi=200, bbox_inches="tight")
plt.close()

plt.figure(figsize=(8,5))
plt.hlines(252.34,2028,2032,zorder=3,color="cyan",ls="-",lw=3)
plt.hlines(253.3,2028,2032,zorder=3,color="magenta",ls="--",lw=3)
plt.plot(time.decimalyear,FSL_r+Att2,alpha=0.8)
plt.xlabel("year (UT)")
plt.ylabel("attenuation + free space loss (dB)")
plt.title("%s, %s, %s"%(f,D,100-p)+"\%"+" availability, $\\theta>15$ deg")
plt.xlim(2028,2032)
plt.grid(alpha=0.5)
plt.legend(["ICD value (252.34dB, winter)","ICD value (253.3dB, all seasons)","expected value"],
           fontsize=14,loc="lower right")
plt.savefig("plots/att_FSL_misasa_year.png", dpi=200, bbox_inches="tight")
plt.close()

plt.figure(figsize=(8,6))
plt.title("%s, %s, %s"%(f,D,100-p)+"\%"+" availability, $\\theta>15$ deg")
plt.vlines(252.34,0,1,zorder=3,color="cyan",ls="-",lw=3)
plt.vlines(253.3,0,1,zorder=3,color="magenta",ls="--",lw=3)
plt.legend(["ICD (252.34dB, winter)","ICD (253.3dB, all seasons)"],loc="upper left")
plt.hist(FSL_r+Att2,cumulative=True,bins=30,lw=3,histtype="step",density=True)
plt.xticks(np.arange(244,255,1))
plt.yticks(np.arange(0,1.1,0.1),["0\%","10\%","20\%","30\%","40\%","50\%","60\%","70\%","80\%","90\%","100\%"])
plt.grid(ls=":",alpha=0.5)
plt.xlabel("Attenuation + Free Space Loss (dB)")
plt.ylabel("cumulative fraction of observing time")
plt.savefig("plots/att_FSL_misasa_year_cumu.png", dpi=200, bbox_inches="tight")
plt.close()

total_10deg = []
total_15deg = []
total_20deg = []
total_25deg = []
for i in range(365):
    use = (time.decimalyear>2028+i/365)&(time.decimalyear<2028+(i+1)/365)&(Roman[1]>10)
    total_10deg.append(np.sum(use)/6) # hour/day
    use = (time.decimalyear>2028+i/365)&(time.decimalyear<2028+(i+1)/365)&(Roman[1]>15)
    total_15deg.append(np.sum(use)/6) # hour/day
    use = (time.decimalyear>2028+i/365)&(time.decimalyear<2028+(i+1)/365)&(Roman[1]>20)
    total_20deg.append(np.sum(use)/6) # hour/day
    use = (time.decimalyear>2028+i/365)&(time.decimalyear<2028+(i+1)/365)&(Roman[1]>25)
    total_25deg.append(np.sum(use)/6) # hour/day
plt.figure(figsize=(7,5))
plt.plot(total_10deg,"-",label="$>10$ deg")
plt.plot(total_15deg,"-",label="$>15$ deg")
plt.plot(total_20deg,"-",label="$>20$ deg")
plt.plot(total_25deg,"-",label="$>25$ deg")
plt.legend()
plt.xlim(0,365)
plt.title("Roman orbit in 2028 seen from MDSS")
plt.yticks(np.arange(0,14,2))
plt.xticks(np.arange(0,365,30.41),np.arange(0,13,1))
plt.grid(ls=":")
plt.ylabel("observable time (hour/day)")
plt.xlabel("end of month")
plt.ylim(0,13)
plt.savefig("plots/obstime_misasa_year.png", dpi=200, bbox_inches="tight")
plt.close()

fig,ax=plt.subplots(2,1,figsize=(6,7),sharey=True,sharex=True)
fig.subplots_adjust(hspace=0.05)
months = ["Jan.","Feb.","Mar.","Apr.","May","Jun.",
          "Jul.","Aug.","Sep.","Oct.","Nov.","Dec."]
for i in range(12):
    use = (time.decimalyear>2028+i/12)&(time.decimalyear<2028+(i+1)/12)&(Roman[1]>15)
    ele_use = Roman[1][use]
    if i<6:
        ax[0].hist(ele_use,bins=np.arange(10,90,2),histtype="step",lw=2,
                   label="%s, %.f hours"%(months[i],len(ele_use)/6))
    else:
        ax[1].hist(ele_use,bins=np.arange(15,90,2),histtype="step",lw=2,
                   label="%s, %.f hours"%(months[i],len(ele_use)/6))
ax[0].set_yticks(np.arange(0,80,20)*6,np.arange(0,80,20))
ax[1].set_xticks(np.arange(0,100,10))
ax[0].legend(bbox_to_anchor=(1,1))
ax[1].legend(bbox_to_anchor=(1,1))
ax[0].set_title("Roman seen from MDSS")
ax[1].set_xlabel("elevation (deg)")
ax[0].set_ylabel("observable time (hours)")
ax[1].set_ylabel("observable time (hours)")
plt.savefig("plots/obstime_misasa_months.png", dpi=200, bbox_inches="tight")
plt.close()


