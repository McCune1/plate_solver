import json
from setup import *
from inv import inventory, argcount
r=inventory(0.35); print('n',len(r),[complex(round(z.real,3),round(z.imag,3)) for z in r],flush=True)
json.dump({'0.35':[[z.real,z.imag] for z in r]},open('inv_cache%s_0.35.json'%TAG,'w'))
for X,Y in ((12,22),(30,45)):
    w,mj=argcount(0.35,X,Y); print('argcount',X,Y,round(w,3),round(mj,3),'found',sum(1 for z in r if z.real<X and z.imag<Y),flush=True)
