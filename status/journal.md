# CryoLab Journal

## 2026-03-12

### Repo Reorganization
- Moved ~92 scripts into organized structure: hardware/{gui,daemons,campaigns,diagnostics}
- Created Series/ hierarchy: Camp1_BlackBox, Camp2_RodMagnet, Camp3_CalibratedCryo, Camp4_NbTiN_Meissner
- Training materials split into laws/ (inviolable) and maxims/ (conditional)
- Added REPO_ROOT to cryolab/__init__.py, extracted LakeshoreProxy to cryolab/hardware.py

### Camp4 NbTiN Meissner — 10G Tc Sweep
- **Goal:** Measure Meissner effect onset via ODMR contrast change across Tc
- **Range:** 11.0K → 15.0K in 0.1K steps (41 stops)
- **Field:** 10G (stepper position -1655000)
- **Status:** Running autonomously, stop 3+ of 41
- Stop 1 (11.0K): 4 pairs, avg contrast 4.4665%, σ=0.0188%
- Stop 2 (11.1K): 10 pairs, avg contrast 4.4726%, σ=0.0462%

### GitHub Setup
- Created PUBLIC repo: bogie970/cryolab
- Whitelist .gitignore: only references/ PDFs and status/ files allowed
- All secrets (Discord webhooks, bot tokens, Twilio creds) excluded
- Added Git & Secrets law to training/laws/safety_limits.md
