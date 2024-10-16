include,"/root/Yorick/STECKMAP/Pierre/POP/pop_paths.i"
include,"sfit.i"
//wavel=[4313.0, 5217.0]
mask=[[5880,5900]]
wavel=[5100.0, 6060.0] // MEGARA LR-V
fv1="/root/Yorick/shared_directory/input/spectrum_xpos_19_ypos_21_fiber_number_0839_NGC7025_LR-V.fits"
a1=convert_all(fv1,z0=0.016571,SNR0=100,nosav=0,noplot=0)
b=bRbasis3([0.63e8,1.78e10],basisfile="/root/Yorick/STECKMAP/Pierre/POP/BASE/MILES_SSP.yor",nbins=30,wavel=wavel,FWHMpix=7.90064096)
x=sfit(a1,b,kin=1,epar=3,nde=40,noskip=1,sav=1,pl=0,pr=1,RMASK=mask)
system, "mv /root/Yorick/shared_directory/input/*.res* /root/Yorick/shared_directory/output/"

include,"/root/Yorick/STECKMAP/Pierre/POP/pop_paths.i"
include,"sfit.i"
//wavel=[4313.0, 5217.0]
mask=[[5880,5900]]
wavel=[5100.0, 6060.0] // MEGARA LR-V
fv1="/root/Yorick/shared_directory/input/spectrum_xpos_19_ypos_22_fiber_number_0840_NGC7025_LR-V.fits"
a1=convert_all(fv1,z0=0.016571,SNR0=100,nosav=0,noplot=0)
x=sfit(a1,b,kin=1,epar=3,nde=40,noskip=1,sav=1,pl=0,pr=1,RMASK=mask)
system, "mv /root/Yorick/shared_directory/input/*.res* /root/Yorick/shared_directory/output/"