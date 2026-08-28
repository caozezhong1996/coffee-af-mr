# -*- coding: utf-8 -*-
"""FAERS control analysis: caffeine x malignant ventricular arrhythmia / cardiac arrest.
Discriminating test for the <40y AF signal: if caffeine disproportionality is also absent
for ventricular outcomes (the AHA-flagged high-dose hazard), the young-age AF signal more
likely reflects reporting structure; if present, dose-toxicity is corroborated.
Mirrors faers_stratified.py conventions (serious reports, same masking list, same ROR).
Output: faers_va_control.json + printed ROR table."""
import urllib.parse, urllib.request, json, time, math, os

BASE = "https://api.fda.gov/drug/event.json"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "faers_va_control.json")

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

SER = 'serious:1'
V2  = 'safetyreportversion:[2 TO 100]'
CAF = 'patient.drug.medicinalproduct:"caffeine"'
MASK_DRUGS = '"amiodarone" OR "dronedarone" OR "adenosine" OR "dofetilide" OR "ibutilide" OR "flecainide" OR "propafenone" OR "sotalol" OR "digoxin" OR "dobutamine" OR "dopamine" OR "epinephrine" OR "norepinephrine" OR "milrinone" OR "theophylline"'
NOTMASK = f'NOT patient.drug.medicinalproduct:({MASK_DRUGS})'

# composite malignant ventricular arrhythmia / arrest outcome (MedDRA PTs)
VA = ('patient.reaction.reactionmeddrapt:("ventricular tachycardia" OR "ventricular fibrillation" '
      'OR "ventricular arrhythmia" OR "torsade de pointes" OR "cardiac arrest" OR "sudden cardiac death")')
VT = 'patient.reaction.reactionmeddrapt:"ventricular tachycardia"'   # single-PT sensitivity

SCOPES = {
    "overall":   None,
    "age_lt40":  'patient.patientonsetage:[0 TO 39] AND patient.patientonsetageunit:801',
    "age_40_64": 'patient.patientonsetage:[40 TO 64] AND patient.patientonsetageunit:801',
    "age_ge65":  'patient.patientonsetage:[65 TO 120] AND patient.patientonsetageunit:801',
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

for name, st in SCOPES.items():
    cond = f" AND {st}" if st else ""
    q(f"{name}_bg_all", f"{SER}{cond}")
    q(f"{name}_bg_va",  f"{SER} AND {VA}{cond}")
    q(f"{name}_caf_total", f"{SER} AND {CAF}{cond}")
    q(f"{name}_caf_va",    f"{SER} AND {CAF} AND {VA}{cond}")

# single-PT VT sensitivity, overall + young
for name in ["overall", "age_lt40"]:
    st = SCOPES[name]; cond = f" AND {st}" if st else ""
    q(f"{name}_bg_vt",  f"{SER} AND {VT}{cond}")
    q(f"{name}_caf_vt", f"{SER} AND {CAF} AND {VT}{cond}")

# robustness for the youngest stratum: dedup + nomask (composite VA)
for tag, extra in [("dedup", f"AND NOT {V2}"), ("nomask", f"AND {NOTMASK}")]:
    s = SCOPES["age_lt40"]
    q(f"age_lt40_{tag}_bg_all", f"{SER} AND {s} {extra}")
    q(f"age_lt40_{tag}_bg_va",  f"{SER} AND {VA} AND {s} {extra}")
    q(f"age_lt40_{tag}_caf_total", f"{SER} AND {CAF} AND {s} {extra}")
    q(f"age_lt40_{tag}_caf_va",    f"{SER} AND {CAF} AND {VA} AND {s} {extra}")

def ror(a, et, ba, bt):
    b = et - a; c = ba - a; d = (bt - ba) - b
    if min(a, b, c, d) <= 0: return None
    r = (a/b)/(c/d)
    se = math.sqrt(1/a + 1/b + 1/c + 1/d)
    return r, math.exp(math.log(r)-1.96*se), math.exp(math.log(r)+1.96*se)

print("\n=== ROR caffeine -> VA composite ===")
for name in SCOPES:
    try:
        r = ror(out[f"{name}_caf_va"], out[f"{name}_caf_total"], out[f"{name}_bg_va"], out[f"{name}_bg_all"])
        print(f"{name}: a={out[f'{name}_caf_va']}/{out[f'{name}_caf_total']} | bg {out[f'{name}_bg_va']}/{out[f'{name}_bg_all']} | ROR {r[0]:.2f} ({r[1]:.2f}-{r[2]:.2f})")
    except Exception as e:
        print(name, "ERR", e)
print("\n=== ROR caffeine -> VT (single PT) ===")
for name in ["overall", "age_lt40"]:
    try:
        r = ror(out[f"{name}_caf_vt"], out[f"{name}_caf_total"], out[f"{name}_bg_vt"], out[f"{name}_bg_all"])
        print(f"{name}: a={out[f'{name}_caf_vt']}/{out[f'{name}_caf_total']} | bg {out[f'{name}_bg_vt']}/{out[f'{name}_bg_all']} | ROR {r[0]:.2f} ({r[1]:.2f}-{r[2]:.2f})")
    except Exception as e:
        print(name, "ERR", e)
print("\n=== <40 robustness (VA composite) ===")
for tag in ["dedup", "nomask"]:
    try:
        r = ror(out[f"age_lt40_{tag}_caf_va"], out[f"age_lt40_{tag}_caf_total"], out[f"age_lt40_{tag}_bg_va"], out[f"age_lt40_{tag}_bg_all"])
        print(f"age_lt40 {tag}: ROR {r[0]:.2f} ({r[1]:.2f}-{r[2]:.2f})")
    except Exception as e:
        print("age_lt40", tag, "ERR", e)
json.dump(out, open(OUT, "w"), indent=1)
print("SAVED", OUT)
