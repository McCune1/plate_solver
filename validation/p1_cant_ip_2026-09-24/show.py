import json,sys
for fn in sys.argv[1:]:
    try: L=[json.loads(l) for l in open(fn)]
    except FileNotFoundError: print(fn,'none'); continue
    print('==',fn, len(L))
    cur=None; row=[]
    for d in L:
        if d['kind']=='inv': print(' inv', d['n_inv'], d.get('NQ'), d.get('MM'))
        elif d['kind']=='scan':
            if d['n']!=cur:
                if row: print('  n=%d cnt=%d nr=%s:'%(cur,cnt,nr),' '.join(row))
                cur=d['n']; row=[]
            cnt=d['cnt']; nr=d['nr']; row.append('%.2f'%d['s'])
        else:
            if row: print('  n=%d cnt=%d:'%(cur,cnt),' '.join(row)); row=[]
            print('  MIN n=%d cnt=%d Om=%.6f s=%.2f edge=%s t=%s'%(d['n'],d['cnt'],d['Om'],d['s'],d['grid_argmin_edge'],d['t']))
    if row: print('  n=%d cnt=%d (partial):'%(cur,cnt),' '.join(row))
