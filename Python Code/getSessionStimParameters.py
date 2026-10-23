# Parameter Range Extraction Across Sessions
# Extracts: maxTargetTimeMS, taskGaborSD, RFGaborSD, speed, catch trial percentage

from usefulFns import *
import numpy as np
import os

fileList = ['Meetz_221010', 'Meetz_221013', 'Meetz_221108', 'Meetz_221110',
            'Meetz_221115', 'Meetz_221117', 'Meetz_221124', 'Meetz_221128',
            'Meetz_221206', 'Meetz_221208', 'Meetz_221229', 'Meetz_230123',
            'Meetz_230126', 'Akshan_240927', 'Akshan_241016', 'Akshan_241017',
            'Akshan_241025']

# Save starting directory once — always return here after each session
home_dir = os.getcwd()

# ── Storage dicts (keyed by monkey) ──────────────────────────────────────────
results = {
    'Meetz':  {'maxTargetTimeMS': [], 'taskGaborSD': [], 'rfGaborSD': [],
               'speed': [], 'catchPct': []},
    'Akshan': {'maxTargetTimeMS': [], 'taskGaborSD': [], 'rfGaborSD': [],
               'speed': [], 'catchPct': []},
}

# ── Per-session collection ────────────────────────────────────────────────────
for fileIterator in fileList:
    print(f'\nProcessing {fileIterator}')
    monkeyName, seshDate = fileIterator.split('_')
    fileName = f'{monkeyName}_{seshDate}_MTNC_Spikes.mat'
    allTrials, header = loadMatFilePyMat(monkeyName, seshDate, fileName)

    # Reset cwd immediately after load — loadMatFilePyMat leaves us inside
    # the session folder; all subsequent sessions need to start from home_dir
    os.chdir(home_dir)

    monkey_key = monkeyName  # 'Meetz' or 'Akshan'

    # ── 1. Max target time (from header) ─────────────────────────────────────
    # header['maxTargetTimeMS'] comes back as {'data': 2000, 'timeMS': ...}
    maxTargetTimeMS = header['maxTargetTimeMS']['data']
    results[monkey_key]['maxTargetTimeMS'].append(maxTargetTimeMS)
    print(f'  maxTargetTimeMS : {maxTargetTimeMS} ms')

    # ── 2. Task Gabor SD — fixed per session, read from trial index 3 (0-based = trials{1,4}) ──
    ref_trial = allTrials[3]  # trials{1,4} in MATLAB → index 3 in Python
    taskGaborSD = ref_trial['taskGabor']['data']['sigmaDeg']
    results[monkey_key]['taskGaborSD'].append(taskGaborSD)
    print(f'  taskGaborSD     : {taskGaborSD}')

    # ── 3. RF Gabor SD — fixed per session, same reference trial ─────────────
    rfGaborSD = ref_trial['rfGabor']['data']['sigmaDeg']
    results[monkey_key]['rfGaborSD'].append(rfGaborSD)
    print(f'  rfGaborSD       : {rfGaborSD}')

    # ── 4. Speed = temporalFreqHz / spatialFreqCPD ───────────────────────────
    spatialFreqCPD = ref_trial['rfGabor']['data']['spatialFreqCPD']
    temporalFreqHz = ref_trial['rfGabor']['data']['temporalFreqHz']
    speed = temporalFreqHz / spatialFreqCPD          # deg/s
    results[monkey_key]['speed'].append(speed)
    print(f'  speed           : {speed:.2f} deg/s  '
          f'(tf={temporalFreqHz} Hz / sf={spatialFreqCPD} cpd)')

    # ── 5. Catch trial percentage ─────────────────────────────────────────────
    n_valid = 0
    n_catch = 0

    for trial in allTrials:
        trial_data = trial.get('trial', {}).get('data', {})

        # Skip instruct trials
        if trial_data.get('instructTrial', 1) == 1:
            continue

        n_valid += 1
        if trial_data.get('catchTrial', 0) == 1:
            n_catch += 1

    catch_pct = (n_catch / n_valid * 100) if n_valid > 0 else float('nan')
    results[monkey_key]['catchPct'].append(catch_pct)
    print(f'  catch trials    : {n_catch}/{n_valid} valid = {catch_pct:.1f}%')


# ── Summary ───────────────────────────────────────────────────────────────────
print('\n' + '='*60)
print('PARAMETER RANGES ACROSS SESSIONS')
print('='*60)

for monkey in ['Meetz', 'Akshan']:
    r = results[monkey]
    if not r['maxTargetTimeMS']:
        continue
    print(f'\n── Monkey {monkey} ──')

    def rng(arr, fmt='.2f'):
        a = np.array(arr)
        u = np.unique(a)
        if len(u) == 1:
            return f'{u[0]:{fmt}}  (constant)'
        return f'{np.min(a):{fmt}} – {np.max(a):{fmt}}  (unique: {sorted(u.tolist())})'

    print(f'  maxTargetTimeMS : {rng(r["maxTargetTimeMS"], ".0f")} ms')
    print(f'  taskGaborSD     : {rng(r["taskGaborSD"])} deg')
    print(f'  rfGaborSD       : {rng(r["rfGaborSD"])} deg')
    print(f'  speed           : {rng(r["speed"])} deg/s')
    print(f'  catch trial %%  : {rng(r["catchPct"])} %%')