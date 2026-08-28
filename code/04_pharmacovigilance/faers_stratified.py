# -*- coding: utf-8 -*-
"""FAERS post-hoc stratification: caffeine x AF disproportionality by sex and age group.
Queries openFDA drug/event endpoint (public, no key). Follows faers_rebuild.py conventions.
Output: faers_stratified.json + printed ROR table."""
import urllib.parse, urllib.request, json, time, math, os

BASE = "https://api.fda.gov/drug/event.json"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "faers_stratified.json")

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
V2  = 'safetyreportversion:[2 TO 100]'   # EXCLUDE follow-up versions (dedup: keep v1 only)
CAF = 'patient.drug.medicinalproduct:"caffeine"'
MASK_DRUGS = '"amiodarone" OR "dronedarone" OR "adenosine" OR "dofetilide" OR "ibutilide" OR "flecainide" OR "propafenone" OR "sotalol" OR "digoxin" OR "dobutamine" OR "dopamine" OR "epinephrine" OR "norepinephrine" OR "milrinone" OR "theophylline"'
NOTMASK = f'NOT patient.drug.medicinalproduct:({MASK_DRUGS})'

STRATA = {
    "male":    'patient.patientsex:1',
    "female":  'patient.patientsex:2',
    "age_lt40":   'patient.patientonsetage:[0 TO 39] AND patient.patientonsetageunit:801',
    "age_40_64":  'patient.patientonsetage:[40 TO 64] AND patient.patientonsetageunit:801',
    "age_ge65":   'patient.patientonsetage:[65 TO 120] AND patient.patientonsetageunit:801',
}

out = {}
if os.path.exists(OUT):
    out = json.load(open(OUT))

def q(key, s):
    if key in out and isinstance(out[key], int):
        print(key, out[key], "(cached)", flush=True); return
    v = count(s)
    out[key] = v
    json.dump(out, open(OUT, "w"), indent=1)
    print(key, v, flush=True)
    time.sleep(2)

for name, st in STRATA.items():
    q(f"{name}_bg_all", f"{SER} AND {st}")
    q(f"{name}_bg_af",  f"{SER} AND {AF} AND {st}")
    q(f"{name}_caf_total", f"{SER} AND {CAF} AND {st}")
    q(f"{name}_caf_af",    f"{SER} AND {CAF} AND {AF} AND {st}")

# robustness for the youngest stratum: dedup + masking
for tag, extra in [("dedup", f"AND NOT {V2}"), ("nomask", f"AND {NOTMASK}")]:
    s = STRATA["age_lt40"]
    q(f"age_lt40_{tag}_bg_all", f"{SER} AND {s} {extra}")
    q(f"age_lt40_{tag}_bg_af",  f"{SER} AND {AF} AND {s} {extra}")
    q(f"age_lt40_{tag}_caf_total", f"{SER} AND {CAF} AND {s} {extra}")
    q(f"age_lt40_{tag}_caf_af",    f"{SER} AND {CAF} AND {AF} AND {s} {extra}")

def ror(a, et, ba, bt):
    b = et - a; c = ba - a; d = (bt - ba) - b
    if min(a, b, c, d) <= 0: return None
    r = (a/b)/(c/d)
    se = math.sqrt(1/a + 1/b + 1/c + 1/d)
    return r, math.exp(math.log(r)-1.96*se), math.exp(math.log(r)+1.96*se)

print("\n=== ROR (caffeine vs stratum-matched background) ===")
rows = {}
for name in STRATA:
    try:
        r = ror(out[f"{name}_caf_af"], out[f"{name}_caf_total"], out[f"{name}_bg_af"], out[f"{name}_bg_all"])
        rows[name] = r
        print(f"{name}: a={out[f'{name}_caf_af']}/{out[f'{name}_caf_total']} | bg {out[f'{name}_bg_af']}/{out[f'{name}_bg_all']} | ROR {r[0]:.2f} ({r[1]:.2f}-{r[2]:.2f})")
    except Exception as e:
        print(name, "ERR", e)
for tag in ["dedup", "nomask"]:
    try:
        r = ror(out[f"age_lt40_{tag}_caf_af"], out[f"age_lt40_{tag}_caf_total"], out[f"age_lt40_{tag}_bg_af"], out[f"age_lt40_{tag}_bg_all"])
        print(f"age_lt40 {tag}: ROR {r[0]:.2f} ({r[1]:.2f}-{r[2]:.2f})")
    except Exception as e:
        print("age_lt40", tag, "ERR", e)
json.dump(out, open(OUT, "w"), indent=1)
print("SAVED", OUT)
