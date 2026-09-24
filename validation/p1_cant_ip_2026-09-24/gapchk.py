from setup import *
import json
for Om in (0.339343, 0.808929, 0.940032):
    for xm,ng,na in ((44.6,20,400),(36.0,20,400),(30.0,10,240)):
        t=time.time(); r=full_search(ip.fast,Om,xmax=xm,ngrid=ng,nax=na)
        rs=sorted(r,key=lambda z:(abs(z.imag),abs(z)))
        print(Om,xm,ng,na,len(r),[complex(round(z.real,2),round(z.imag,2)) for z in rs][:12],round(time.time()-t),flush=True)
