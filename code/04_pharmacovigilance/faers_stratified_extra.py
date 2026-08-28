# -*- coding: utf-8 -*-
"""Extra robustness queries for age 40-64 and >=65 strata (dedup + masking)."""
import urllib.parse, urllib.request, json, time, math, os

BASE = "https://api.fda.gov/drug/event.json"
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "faers_stratified.json")

def count(search, tries=4):
    url = BASE + "?search=" + urllib.parse.quote(search) + "&limit=1"
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=40) as r:
                return json.load(r)["meta"]["results"]["total"]
        except urllib.error.HTTPError as e:
            if e.code == 404: return 0
            if i == tries-1: return f"ERR:{e}"
            time.sleep(6*(i+1))
        except Exception as e:
            if i == tries-1: return f"ERR:{e}"
            time.sleep(6*(i+1))

AF  = 'patient.reaction.reactionmeddrapt:"atrial fibrillation"'
SER = 'serious:1'
V2  = 'safetyreportversion:[2 TO 100]'
CAF = 'patient.drug.medicinalproduct:"caffeine"'
MASK_DRUGS = '"amiodarone" OR "dronedarone" OR "adenosine" OR "dofetilide" OR "ibutilide" OR "flecainide" OR "propafenone" OR "sotalol" OR "digoxin" OR "dobutamine" OR "dopamine" OR "epinephrine" OR "norepinephrine" OR "milrinone" OR "theophylline"'
NOTMASK = f'NOT patient.drug.medicinalproduct:({MASK_DRUGS})'
STRATA2 = {
 'age_40_64': 'patient.patientonsetage:[40 TO 64] AND patient.patientonsetageunit:801',
 'age_ge65':  'patient.patientonsetage:[65 TO 120] AND patient.patientonsetageunit:801',
}

out = json.load(open(OUT)) if os.path.exists(OUT) else {}
def q(key, s):
    if key in out and isinstance(out[key], int):
        print(key, out[key], "(cached)", flush=True); return
    v = count(s); out[key] = v
    json.dump(out, open(OUT, "w"), indent=1)
    print(key, v, flush=True); time.sleep(2)

for nm, s in STRATA2.items():
    for tag, extra in [('dedup', f'AND NOT {V2}'), ('nomask', f'AND {NOTMASK}')]:
        q(f'{nm}_{tag}_bg_all',   f'{SER} AND {s} {extra}')
        q(f'{nm}_{tag}_bg_af',    f'{SER} AND {AF} AND {s} {extra}')
        q(f'{nm}_{tag}_caf_total',f'{SER} AND {CAF} AND {s} {extra}')
        q(f'{nm}_{tag}_caf_af',   f'{SER} AND {CAF} AND {AF} AND {s} {extra}')

def ror(a, et, ba, bt):
    b = et - a; c = ba - a; d = (bt - ba) - b
    if min(a, b, c, d) <= 0: return None
    r = (a/b)/(c/d); se = math.sqrt(1/a + 1/b + 1/c + 1/d)
    return r, math.exp(math.log(r)-1.96*se), math.exp(math.log(r)+1.96*se)

for nm in STRATA2:
    for tag in ['dedup','nomask']:
        r = ror(out[f'{nm}_{tag}_caf_af'], out[f'{nm}_{tag}_caf_total'], out[f'{nm}_{tag}_bg_af'], out[f'{nm}_{tag}_bg_all'])
        print(f'{nm} {tag}: ROR {r[0]:.2f} ({r[1]:.2f}-{r[2]:.2f})')
print("DONE")
