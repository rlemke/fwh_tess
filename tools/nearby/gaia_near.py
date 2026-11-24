
import urllib.request, urllib.parse, time, sys, re
base='https://gea.esac.esa.int/tap-server/tap/async'
q="SELECT source_id, ra, dec, parallax, pmra, pmdec, phot_g_mean_mag, bp_rp FROM gaiadr3.gaia_source WHERE parallax > 25 AND parallax_over_error > 10"
data=urllib.parse.urlencode({'REQUEST':'doQuery','LANG':'ADQL','FORMAT':'csv','QUERY':q,'PHASE':'RUN'}).encode()
class NoRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, *a, **k): return None
op=urllib.request.build_opener(NoRedirect)
try:
    r=op.open(urllib.request.Request(base,data=data)); loc=r.headers.get('Location')
except urllib.error.HTTPError as e:
    loc=e.headers.get('Location')
print('job', loc, flush=True)
while True:
    ph=urllib.request.urlopen(loc+'/phase',timeout=60).read().decode()
    print(ph, flush=True)
    if ph in ('COMPLETED','ERROR','ABORTED'): break
    time.sleep(15)
if ph=='COMPLETED':
    open(sys.argv[1],'wb').write(urllib.request.urlopen(loc+'/results/result',timeout=600).read())
    print('saved')
