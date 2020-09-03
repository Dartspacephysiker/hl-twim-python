import hltwimm;
import numpy as np

# import importlib;
# importlib.reload(hltwim);

# ORIG
doy,ut=np.int32(300),np.float32(12.)
gdlat,gdlon=np.float32(78.15),np.float32(16.03)
kp=np.float32(6.)
# w,mw = hltwimm.hltwim(doy,ut,gdlat,gdlon,kp);

# NUMPY
doy,ut=np.asfortranarray(np.array([doy])),np.asfortranarray(np.array([ut]))
gdlat,gdlon=np.asfortranarray(np.array([gdlat])),np.asfortranarray(np.array([gdlon]))
kp=np.asfortranarray(np.array([kp]))
w,mw = np.asfortranarray(np.zeros(2)),np.asfortranarray(np.zeros(2))
print(doy,ut,gdlat,gdlon,kp,w,mw)
print(doy.dtype,ut.dtype,gdlat.dtype,gdlon.dtype,kp.dtype,w.dtype,mw.dtype)
w[:],mw[:] = hltwimm.hltwim(doy,ut,gdlat,gdlon,kp);
# print(w,mw)

# mlt = np.asfortranarray([12])
# mlat,mlon = gdlat,gdlon
# mmwind,mzwind = hltwimm.hltwimx(doy,mlt,mlat,kp,mlon);
# print(mmwind,mzwind)
