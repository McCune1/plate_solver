import json,sys
from collections import defaultdict
for fn in sys.argv[1:]:
    G=defaultdict(list); R=defaultdict(list)
    for l in open(fn):
        d=json.loads(l)
        (G if d['kind']=='g' else R)[d['n']].append(d)
    for n in sorted(set(G)|set(R)):
        g=sorted(G[n],key=lambda d:d['Om'])
        print('== %s n=%d cnt=%s pts=%d range %.4f-%.4f'%(fn,n,g[0]['cnt'] if g else '?',len(g),g[0]['Om'] if g else 0,g[-1]['Om'] if g else 0))
        print('   u sign:', ''.join('+' if d['u'][0]>0 else '-' for d in g))
        print('   c sign:', ''.join('+' if d['c'][0]>0 else '-' for d in g))
        print('   u smin:', ' '.join('%.1f'%d['u'][2] for d in g))
        for r in sorted(R[n],key=lambda d:d['Om']): print('   ROOT %s Om=%.6f smin=%.2f ld=%.2f'%(r['mat'],r['Om'],r['smin'],r['ld']))
