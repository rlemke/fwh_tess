
import json, statistics, collections
rows=json.load(open('./near40.json'))
free=[r for r in rows if not (r['toi'] or r['ctoi'] or r['planet']) and len(r['sectors'])>=3 and r['bp_rp'] and float(r['bp_rp'])>1.8]
free.sort(key=lambda r:r['d_pc'])
# dedupe tic (multiple gaia matches)
seen=set(); out=[]
for r in free:
    if r['tic'] in seen: continue
    seen.add(r['tic']); out.append(r)
pick=out[:300]
ns=[len(r['sectors']) for r in pick]
print('n',len(pick),'d max',pick[-1]['d_pc'],'sectors median',statistics.median(ns),'max',max(ns),'sum',sum(ns))
print('spans (yrs) median', statistics.median((max(r['sectors'])-min(r['sectors']))*27.4/365 for r in pick))
print(collections.Counter(min(n,20) for n in ns).most_common())
json.dump([r['tic'] for r in pick], open('./nearby300.json','w'))
json.dump(pick, open('./nearby300_meta.json','w'), indent=0)
