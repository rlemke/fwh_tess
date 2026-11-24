
import json, glob, re, csv, math, collections, sys
D='' + __import__('os').environ.get('TESS_DATA', 'data') + ''
by=collections.defaultdict(set); pos={}
for p in sorted(glob.glob(D+'/lc/manifest-s*.json')):
    s=int(re.search(r's(\d{4})',p).group(1))
    if s>96: continue
    for t in json.load(open(p))['targets']: by[int(t)].add(s)
for p in sorted(glob.glob('./tl/S*.txt')):
    for line in open(p):
        if line.startswith('#') or not line.strip(): continue
        f=line.split()
        try: pos.setdefault(int(f[0]), (float(f[4]), float(f[5]), float(f[3])))
        except (ValueError, IndexError): pass
print('2-min stars', len(by), 'with coords', sum(1 for t in by if t in pos))
cell=collections.defaultdict(list)
for t,(ra,dec,tm) in pos.items():
    if t in by: cell[(int(dec*60), int(ra*60*math.cos(math.radians(dec))) )].append((t,ra,dec,tm))
def near(ra,dec):
    k0=int(dec*60); best=None
    for dk in (-1,0,1):
        dd=k0+dk
        c=math.cos(math.radians(dec))
        k1=int(ra*60*c)
        for dr in (-1,0,1):
            for (t,r2,d2,tm) in cell.get((dd,k1+dr),[]):
                sep=math.hypot((r2-ra)*c,(d2-dec))*3600
                if sep<4 and (best is None or sep<best[1]): best=(t,sep,tm)
    return best
rows=[]
for r in csv.DictReader(open(sys.argv[1])):
    ra,dec,plx=float(r['ra']),float(r['dec']),float(r['parallax'])
    pmra=float(r['pmra'] or 0); pmdec=float(r['pmdec'] or 0)
    c=math.cos(math.radians(dec))
    ra0=ra-pmra*16/3.6e6/c; dec0=dec-pmdec*16/3.6e6
    m=near(ra0,dec0)
    if m: rows.append({'tic':m[0],'sep':round(m[1],2),'tmag':m[2],'d_pc':round(1000/plx,2),'bp_rp':r['bp_rp'],'g':r['phot_g_mean_mag'],'gaia':r['source_id'],'sectors':sorted(by[m[0]])})
print('gaia d<40pc', sum(1 for _ in open(sys.argv[1]))-1, 'matched to 2-min stars', len(rows))
def tics(path, col):
    out=set()
    with open(path) as f:
        lines=[l for l in f if not l.startswith('#')]
    for r in csv.DictReader(lines):
        v=r.get(col)
        try: out.add(int(float(v)))
        except: pass
    return out
toi=tics(D+'/tois-2026-11-23.csv','TIC')
ctoi=set()
for l in open(D+'/ctois-2026-11-23.csv'):
    v=l.split('|')[0].strip()
    if v.isdigit(): ctoi.add(int(v))
cp=tics(D+'/pscomppars-2026-11-23.csv','tic_id') if False else set()
import io
# pscomppars tic_id like 'TIC 123'
with open(D+'/pscomppars-2026-11-23.csv') as f:
    lines=[l for l in f if not l.startswith('#')]
for r in csv.DictReader(lines):
    v=(r.get('tic_id') or '').replace('TIC','').strip()
    if v.isdigit(): cp.add(int(v))
print('toi hosts',len(toi),'ctoi',len(ctoi),'confirmed',len(cp))
for r in rows:
    r['toi']=r['tic'] in toi; r['ctoi']=r['tic'] in ctoi; r['planet']=r['tic'] in cp
json.dump(rows, open(sys.argv[2],'w'))
free=[r for r in rows if not (r['toi'] or r['ctoi'] or r['planet'])]
for k in (1,2,3,4,6,10):
    sel=[r for r in free if len(r['sectors'])>=k]
    print('>=',k,'sectors, no TOI/CTOI/planet:', len(sel), ' d<20:', sum(1 for r in sel if r['d_pc']<20), ' M dwarfs (bp_rp>1.8):', sum(1 for r in sel if r['bp_rp'] and float(r['bp_rp'])>1.8))
import statistics
sel=[r for r in free if len(r['sectors'])>=3]
print('sectors per star (>=3): median', statistics.median(len(r['sectors']) for r in sel), 'total star-sectors', sum(len(r['sectors']) for r in sel))
