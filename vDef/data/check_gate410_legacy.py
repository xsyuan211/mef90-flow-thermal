#!/usr/bin/env python3
"""Gate 4.10 -- root-level stale source MIRROR removal verdict.

Pure standard library (CI has no numpy / PyYAML).  The repository used to
carry a complete set of root-level Fortran mirrors plus top-level Makefile
mirrors.  Every build rule in the top Makefile, in the five live subdir
Makefiles and in GitHub Actions points at a SUBDIRECTORY source, so the
root copies never entered the compiler.  56 of the 66 same-named Fortran
mirrors were byte-identical dead weight; the other 10 (and one Makefile
mirror) had DIVERGED from the live sources -- the classic "edit the wrong
file" trap.  Gate 4.10 archives all 93 of them into legacy/root-mirror/
with git mv (nothing deleted) and hard-guards the invariant.

Modes
-----
  report : check_gate410_legacy.py
      INFORMATIONAL, always exits 0.

  gate   : check_gate410_legacy.py --gate
      HARD Gate 4.10, prints GATE410_RESULT: PASS|FAIL and exits 0 / 1.
      Missing / unreadable / tampered data ALWAYS fails.

  selftest : check_gate410_legacy.py --selftest
      Offline replay with a synthetic state: 1 positive + 9 negative
      controls.  No git, no PETSc, no CI.

Criteria
--------
  W1  no tracked .F90/.inc remains at repo root; any untracked ones must
      be gitignored local diagnostics
  W2  legacy/root-mirror/ holds exactly EXPECTED_ARCHIVE files and the
      MANIFEST.sha256 entry set matches it exactly
  W3  every manifest sha256 matches the archived file on disk
  W4  all 66 same-named Fortran mirrors have a live counterpart in the
      active source subdirectories
  W5  the five live Makefiles reference only subdirectory .F90 sources
      and define BASE_DIR as ${MEF90_DIR}/<subdir>
  W6  the top Makefile and build.yml use only -f ../../<subdir>/Makefile
      and never reference an archived mirror Makefile
  W7  no build entry point mentions legacy/ (archive cannot be compiled)
  W8  README.md documents the archive, the diverged files and the no-build
      contract
"""
import copy
import hashlib
import os
import re
import subprocess
import sys

EXPECTED_ARCHIVE = 93
EXPECTED_MIRRORS = 66
LIVE_SUBDIRS = ('MEF90', 'm_DefMech', 'm_HeatXfer', 'm_Flow', 'vDef')
ROOT_ONLY_ARCHIVED = ('vDefDefault_remote.F90', 'test20.F90')
DIVERGED = ('m_MEF90_DefMech.F90', 'm_MEF90_DefMechAssembly.F90',
            'm_MEF90_DefMechCtx.F90', 'm_MEF90_Flow.F90',
            'm_MEF90_FlowAssembly.F90', 'm_MEF90_FlowCtx.F90',
            'm_MEF90_HeatXferAssembly.F90', 'm_MEF90_HeatXferCtx.F90',
            'm_vDefDefault.F90', 'vDef.F90', 'Makefile_m_Flow')
README_KEYS = ('不参与构建', 'git mv', 'MANIFEST.sha256')
MANIFEST_LINE = re.compile(r'^([0-9a-f]{64})\s{2}(\S+)\s*$')
ROOT_SRC_REF = re.compile(r'\$\{MEF90_DIR\}/([A-Za-z0-9_]+\.F90)')
F_TARGET = re.compile(r'-f\s+(\S*Makefile\S*)')


def read_text(path):
    try:
        with open(path, 'r', encoding='utf-8', errors='replace') as fh:
            return fh.read()
    except Exception:                                  # noqa: BLE001
        return ''


def sha256(path):
    """Hash CRLF-normalised bytes so the verdict is stable on Windows
    (core.autocrlf=true) AND after CI's dos2unix pass over the archive.
    The manifest is generated with exactly the same normalisation."""
    try:
        with open(path, 'rb') as fh:
            data = fh.read().replace(b'\r\n', b'\n')
        return hashlib.sha256(data).hexdigest()
    except Exception:                                  # noqa: BLE001
        return ''


def git_lines(root, *args):
    try:
        p = subprocess.run(['git'] + list(args), cwd=root,
                           capture_output=True, text=True)
        return p.stdout.splitlines() if p.returncode == 0 else []
    except Exception:                                  # noqa: BLE001
        return []


def git_ignored(root, name):
    try:
        return subprocess.run(['git', 'check-ignore', '-q', '--', name],
                              cwd=root, capture_output=True).returncode == 0
    except Exception:                                  # noqa: BLE001
        return False


def collect(root):
    """Observe the real repository; every field feeds one W-criterion."""
    s = {}
    tracked = git_lines(root, 'ls-files')
    top_tracked = [p for p in tracked if '/' not in p]
    s['root_f90'] = sorted(p for p in top_tracked
                           if p.endswith(('.F90', '.inc')))
    s['tracked_count'] = len(tracked)

    disk = []
    try:
        disk = [f for f in os.listdir(root)
                if os.path.isfile(os.path.join(root, f))
                and f.endswith(('.F90', '.inc'))]
    except Exception:                                  # noqa: BLE001
        pass
    extra = sorted(set(disk) - set(top_tracked))
    s['extra'] = extra
    s['notignored'] = [f for f in extra if not git_ignored(root, f)]

    ld = os.path.join(root, 'legacy', 'root-mirror')
    try:
        files = sorted(f for f in os.listdir(ld)
                       if os.path.isfile(os.path.join(ld, f))
                       and f not in ('README.md', 'MANIFEST.sha256'))
    except Exception:                                  # noqa: BLE001
        files = []
    s['legacy_files'] = files

    entries = {}
    for line in read_text(os.path.join(ld, 'MANIFEST.sha256')).splitlines():
        m = MANIFEST_LINE.match(line)
        if m:
            entries[m.group(2)] = m.group(1)
    s['manifest_entries'] = entries
    s['manifest_extra'] = sorted(set(entries) - set(files))
    s['manifest_missing'] = sorted(set(files) - set(entries))
    s['manifest_hash_ok'] = bool(files) and all(
        name in entries and entries[name] == sha256(os.path.join(ld, name))
        for name in files)

    live = {}
    for p in tracked:
        if '/' in p and not p.startswith('legacy/') and p.endswith(('.F90', '.inc')):
            live.setdefault(os.path.basename(p), []).append(p)
    mirrors = [f for f in files
               if f.endswith(('.F90', '.inc'))
               and f not in ROOT_ONLY_ARCHIVED
               and not f.startswith('test_')]
    s['mirrors'] = mirrors
    s['missing_live'] = sorted(f for f in mirrors if f not in live)

    bad_refs, base_defs = [], {}
    for d in LIVE_SUBDIRS:
        mf = os.path.join(root, d, 'Makefile')
        txt = read_text(mf)
        if not txt:
            bad_refs.append(d + ':missing-Makefile')
            continue
        bad_refs += ['%s:root-src:%s' % (d, m.group(1))
                     for m in ROOT_SRC_REF.finditer(txt)]
        for m in re.finditer(r'^\s*BASE_DIR\s*=\s*(\S+)', txt, re.M):
            val = m.group(1)
            base_defs[d] = val
            if not val.startswith('${MEF90_DIR}/' + d):
                bad_refs.append('%s:BASE_DIR=%s' % (d, val))
        for m in re.finditer(r'^\s*BASE_DIR\s*=\s*\$\{MEF90_DIR\}\s*$', txt, re.M):
            bad_refs.append('%s:BASE_DIR-points-at-root' % d)
    s['bad_refs'] = bad_refs
    s['base_defs'] = base_defs

    top = read_text(os.path.join(root, 'Makefile'))
    ci = read_text(os.path.join(root, '.github', 'workflows', 'build.yml'))
    f_targets = sorted(set(F_TARGET.findall(top + '\n' + ci)))
    s['f_targets'] = f_targets
    # Any -f target must be a SUBDIRECTORY Makefile. Do not enumerate the
    # directories: the top Makefile has historical optional targets such as
    # ../../m_Elasticity/Makefile (no such dir today) that are still
    # structurally correct. Root-level mirror Makefiles are caught by
    # mirror_refs below, and legacy/ by W7.
    sub_re = re.compile(r'^\.\./\.\./[A-Za-z0-9_]+/Makefile$')
    s['bad_f'] = [t for t in f_targets if not sub_re.match(t)]
    mirror_names = [f for f in files if f.startswith('Makefile')]
    blob = top + '\n' + ci
    s['mirror_refs'] = sorted(
        n for n in mirror_names
        if re.search(r'(?:^|[\s/])' + re.escape(n) + r'(?:$|[\s"\'])', blob))

    build_entries = [os.path.join(root, 'Makefile'),
                     os.path.join(root, 'Makefile.include'),
                     os.path.join(root, '.github', 'workflows', 'build.yml')]
    for d in LIVE_SUBDIRS:
        build_entries += [os.path.join(root, d, 'Makefile'),
                          os.path.join(root, d, 'Makefile.include')]
    legacy_refs = []
    for p in build_entries:
        for i, line in enumerate(read_text(p).splitlines(), 1):
            if 'legacy/' in line or 'legacy\\' in line:
                legacy_refs.append('%s:%d' % (os.path.relpath(p, root), i))
    s['legacy_refs'] = legacy_refs

    rd = read_text(os.path.join(ld, 'README.md'))
    s['readme_ok'] = (bool(rd)
                      and all(n in rd for n in DIVERGED)
                      and all(k in rd for k in README_KEYS))
    return s


def evaluate(s):
    chk = {}
    chk['W1'] = (not s['root_f90'] and not s['notignored'],
                 'tracked root .F90/.inc=%d; untracked-not-ignored=%s'
                 % (len(s['root_f90']), s['notignored'] or 'none'))
    chk['W2'] = (len(s['legacy_files']) == EXPECTED_ARCHIVE
                 and len(s['manifest_entries']) == EXPECTED_ARCHIVE
                 and not s['manifest_extra'] and not s['manifest_missing'],
                 'archived=%d (want %d), manifest=%d, extra=%s, missing=%s'
                 % (len(s['legacy_files']), EXPECTED_ARCHIVE,
                    len(s['manifest_entries']), s['manifest_extra'] or 'none',
                    s['manifest_missing'] or 'none'))
    chk['W3'] = (s['manifest_hash_ok'] and not s['manifest_extra'],
                 'all %d archived sha256 match MANIFEST.sha256: %s'
                 % (len(s['legacy_files']), s['manifest_hash_ok']))
    chk['W4'] = (len(s['mirrors']) == EXPECTED_MIRRORS
                 and not s['missing_live'],
                 'same-named mirrors=%d (want %d), without live twin=%s'
                 % (len(s['mirrors']), EXPECTED_MIRRORS,
                    s['missing_live'] or 'none'))
    chk['W5'] = (not s['bad_refs'] and len(s['base_defs']) >= 4,
                 'live Makefile root-source refs / bad BASE_DIR=%s; '
                 'BASE_DIR defs=%d'
                 % (s['bad_refs'] or 'none', len(s['base_defs'])))
    chk['W6'] = (bool(s['f_targets']) and not s['bad_f']
                 and not s['mirror_refs'],
                 '-f targets=%d, invalid=%s, archived-mirror refs=%s'
                 % (len(s['f_targets']), s['bad_f'] or 'none',
                    s['mirror_refs'] or 'none'))
    chk['W7'] = (not s['legacy_refs'],
                 'build entries referencing legacy/ = %s'
                 % (s['legacy_refs'] or 'none'))
    chk['W8'] = (s['readme_ok'],
                 'README documents archive + %d diverged files + no-build contract: %s'
                 % (len(DIVERGED), s['readme_ok']))
    return chk


def print_report(chk):
    print('=== Gate 4.10 -- root legacy mirror removal verdict ===')
    for name in sorted(chk):
        ok, detail = chk[name]
        print('  %-4s %-4s %s' % (name, 'PASS' if ok else 'FAIL', detail))
    return [k for k, v in chk.items() if not v[0]]


# --------------------------------------------------------------------------
# Offline selftest -- synthetic state, no git / no CI
# --------------------------------------------------------------------------
def _good_state():
    mirrors = ['m%02d.F90' % i for i in range(EXPECTED_MIRRORS)]
    others = (['vDefDefault_remote.F90', 'test20.F90']
              + ['test_%02d.F90' % i for i in range(16)]
              + ['Makefile_x%02d' % i for i in range(9)])
    files = sorted(mirrors + others)
    assert len(files) == EXPECTED_ARCHIVE, len(files)
    return dict(
        root_f90=[], notignored=[], extra=['FlowCtx.F90'],
        legacy_files=files,
        manifest_entries={f: 'a' * 64 for f in files},
        manifest_extra=[], manifest_missing=[], manifest_hash_ok=True,
        mirrors=mirrors, missing_live=[],
        bad_refs=[], base_defs={d: '${MEF90_DIR}/' + d for d in
                                ('m_DefMech', 'm_HeatXfer', 'm_Flow', 'vDef')},
        f_targets=['../../MEF90/Makefile', '../../m_Flow/Makefile'],
        bad_f=[], mirror_refs=[], legacy_refs=[], readme_ok=True)


def selftest():
    print('=== Gate 4.10 offline selftest (synthetic state) ===')
    good = _good_state()

    def mutated(key, value):
        s = copy.deepcopy(good)
        s[key] = value
        return s

    cases = [
        ('positive (clean archive)', good, True),
        ('mirror copied back to root', mutated('root_f90', ['m_MEF90_Flow.F90']), False),
        ('untracked unignored .F90 at root', mutated('notignored', ['New.F90']), False),
        ('one archived file lost',
         mutated('legacy_files', good['legacy_files'][:-1]), False),
        ('manifest hash tampered', mutated('manifest_hash_ok', False), False),
        ('manifest has stray entry',
         mutated('manifest_extra', ['ghost.F90']), False),
        ('mirror lost its live twin',
         mutated('missing_live', ['m00.F90']), False),
        ('live Makefile points at root source',
         mutated('bad_refs', ['m_Flow:root-src:m_MEF90_Flow.F90']), False),
        ('build.yml uses archived Makefile',
         mutated('bad_f', ['../../Makefile_m_Flow']), False),
        ('build entry references legacy/',
         mutated('legacy_refs', ['Makefile:12']), False),
        ('README missing the diverged list',
         mutated('readme_ok', False), False),
    ]
    all_ok = True
    for label, state, want in cases:
        fails = [k for k, v in evaluate(state).items() if not v[0]]
        got = not fails
        ok = got == want
        all_ok &= ok
        print('  %-42s -> %-4s (want %-4s) %s  failing=%s'
              % (label, 'PASS' if got else 'FAIL',
                 'PASS' if want else 'FAIL',
                 'OK' if ok else 'SELFTEST ERROR', ','.join(fails) or 'none'))
    print('SELFTEST_RESULT: %s' % ('PASS' if all_ok else 'FAIL'))
    return 0 if all_ok else 1


def main(argv):
    if argv and argv[0] == '--selftest':
        return selftest()
    gate = bool(argv) and argv[0] == '--gate'
    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.normpath(os.path.join(here, '..', '..'))
    chk = evaluate(collect(root))
    fails = print_report(chk)
    print('  repo root: %s' % root)
    if gate:
        if fails or len(chk) < 8:
            print('GATE410_RESULT: FAIL (%s)'
                  % (','.join(fails) or 'incomplete criteria'))
            return 1
        print('GATE410_RESULT: PASS (W1 root clean, W2/W3 93 archived files '
              'sha256-verified, W4 66 live twins, W5-W7 build only uses '
              'subdirectory sources, W8 archive documented)')
        return 0
    print('VERDICT: %s (report mode, informational only)'
          % ('FAIL ' + ','.join(fails) if fails else 'PASS'))
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
