import math
import numpy as np
import matplotlib.pyplot as plt

import sys,pathlib,collections,importlib.abc; r=next(p for p in pathlib.Path(__file__).resolve().parents if p.name=="PSPL_Rocket_A");   m=collections.defaultdict(list); [m[f.stem.casefold()].append(f) for f in r.rglob("*.py") if f.stem.isidentifier() and f.name != "__init__.py" and not (set(f.relative_to(r).parts)&{".git",".venv","__pycache__","build","dist"})]; dup={k:v for k,v in m.items() if   len(v)>1}; sys.meta_path.insert(0,type("AmbiguousBareImportBlocker",(importlib.abc.MetaPathFinder,),{"find_spec":lambda   self,fullname,path=None,target=None: (_ for _ in ()).throw(ImportError(f"The import {fullname!r} could refer to any following packages: "+" ".join("\n\t" + str(p.relative_to(r)) for p in dup[fullname.casefold()])+f"\n\nSpecify which package it is by using the folder.\nFor example:\n\t'import SFD.{fullname}'")) if "." not in fullname and fullname.casefold() in dup else None})()); sys.path.insert(0,str(r)) if str(r) not   in sys.path else None; [sys.path.append(str(v[0].parent)) for k,v in m.items() if k not in dup and str(v[0].parent) not in sys.path]
import constants as c
import verifyrelief

Cv = 0.8
P1 = 300 * c.BAR2PA * c.PA2PSI
T = 300 * c.KELVIN2RANK
Gs = 0.967

qdot = Cv * 13.61*P1*math.sqrt(1/(T*Gs))
print(f"qdot: {qdot:.2f}")

scfm_history = []
pressure_history = np.linspace(3000, 5000, 100)
for P1 in pressure_history:
    qdot = Cv * 13.61*P1*math.sqrt(1/(T*Gs))
    scfm_history.append(qdot)

plt.plot(pressure_history, scfm_history)
plt.xlabel('Pressure (psia)')
plt.ylabel('Flow Rate (SCFM)')
plt.title('Reg Fail Flow Rate vs Inlet Pressure')
plt.axvline(300 * c.BAR2PA * c.PA2PSI, color='r', linestyle='--', label='COPV MAWP')
plt.axvline(3300, color='b', linestyle='--', label='Nominal Inlet Pressure')
plt.axhline(729.33, color='g', linestyle='--', label='Max Relief Flow Rate at set 475 psig')
plt.axhline(842.38, color='c', linestyle='--', label='Max Relief Flow Rate at set 550 psig')
plt.axhline(729.33 * 3, color='m', linestyle='--', label='3X Max Relief Flow Rate at set 475 psig')
plt.axhline(842.38 * 3, color='y', linestyle='--', label='3X Max Relief Flow Rate at set 550 psig')
plt.ylim(0, 3000)
plt.xlim(3000, 5000)
plt.legend()
plt.show()