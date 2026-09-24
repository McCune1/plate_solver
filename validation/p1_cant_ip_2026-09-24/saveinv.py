import json,sys
from inv import inventory
from setup import TAG
out={}
for Om in map(float,sys.argv[1:]):
    r=inventory(Om); out[str(Om)]=[[z.real,z.imag] for z in r]; print(Om,len(r),flush=True)
json.dump(out,open('inv_cache%s_%s.json'%(TAG,sys.argv[1]),'w'))
