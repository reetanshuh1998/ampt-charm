#!/usr/bin/env python3
"""code_map.py OUT.tex : LaTeX tables of every code change, generated from the repository itself so that the line
numbers in docs/CHANGES.pdf are exact for build/ampt-src/ (= official AMPT + patches/*.patch + src_new/).
  Table 1: each hunk of the patch (stock file, line range in the patched file, lines added/removed, its first comment)
  Table 2: each routine of src_new/charm_mod.f and src_new/coal_new.f with its line range
  Table 3: each environment switch, the line that reads it, and its value in the reference configuration"""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = sys.argv[1]


def tex(s, brk=True):
    s = s.replace('\\', '/').replace('_', r'\_').replace('&', r'\&').replace('%', r'\%').replace('#', r'\#')
    s = s.replace('^', r'\^{}').replace('~', r'\~{}').replace('$', r'\$').replace('{', r'\{').replace('}', r'\}')
    return s.replace(r'\_', r'\_\allowbreak{}') if brk else s


L = []
# ------------------------------------------------------------------ table 1: patch hunks
hunks, cur = [], None
for line in open(ROOT / 'patches/ampt-v1.26t9b-v2.26t9b_charm.patch', errors='replace'):
    if line.startswith('+++ '):
        fn = line.split()[1].split('/', 1)[1]
    elif line.startswith('@@'):
        m = re.match(r'@@ -(\d+)(?:,(\d+))? \+(\d+)(?:,(\d+))? @@', line)
        cur = dict(file=fn, start=int(m.group(3)), n=int(m.group(4) or 1), add=0, rem=0, note='', first=None, last=None, k=0)
        hunks.append(cur)
    elif cur is not None and line[:1] in '+- ':
        if line[0] == '+':
            cur['add'] += 1
            if cur['first'] is None:
                cur['first'] = cur['start'] + cur['k']
            cur['last'] = cur['start'] + cur['k']
            t = line[1:].rstrip()
            if not cur['note'] and re.match(r'^[cC*!]', t) and len(t.strip()) > 3:
                cur['note'] = re.sub(r'^[cC*!]\s+', '', t).strip()
        elif line[0] == '-':
            cur['rem'] += 1
        if line[0] != '-':
            cur['k'] += 1
for h in hunks:
    if not h['note']:     # code-only hunk: the first added statement is shown below
        h['note'] = 'code:'
L += [r'{\small\begin{longtable}{@{}p{2.7cm}p{1.5cm}p{1.3cm}p{9.4cm}@{}}',
      r'\caption{Every hunk of \texttt{patches/ampt-v1.26t9b-v2.26t9b\_charm.patch}. Lines = the added lines in the patched file '
      r'\texttt{build/ampt-src/}; $+/-$ = lines added/removed. The note is the first comment of the hunk (tagged '
      r'\texttt{cdmeson} in the source).}\label{tab:hunks}\\',
      r'\toprule file & lines & $+/-$ & note \\ \midrule \endfirsthead', r'\toprule file & lines & $+/-$ & note \\ \midrule \endhead']
src = {}
for h in hunks:
    f = h['file']
    if f not in src:
        src[f] = (ROOT / 'build/ampt-src' / f).read_text(errors='replace').splitlines() if (ROOT / 'build/ampt-src' / f).exists() else []
    note = h['note']
    if note.startswith('code:'):
        a = h['first'] or h['start']
        note = 'code: ' + (src[f][a - 1].strip() if src[f] and a - 1 < len(src[f]) else '')
    rng = f'{h["first"]}--{h["last"]}' if h['first'] and h['last'] != h['first'] else f'{h["first"] or h["start"]}'
    L.append(rf'\texttt{{{tex(f)}}} & {rng} & +{h["add"]}/$-${h["rem"]} & {tex(note[:110])} \\')
L += [r'\bottomrule \end{longtable}}', '']
# ------------------------------------------------------------------ table 2: routines
L += [r'{\small\begin{longtable}{@{}p{1.9cm}p{2.0cm}p{1.6cm}p{9.4cm}@{}}',
      r'\caption{Routines of the new files (\texttt{src\_new/}), with their line ranges and the first line of their '
      r'header comment.}\label{tab:routines}\\', r'\toprule file & routine & lines & purpose \\ \midrule \endfirsthead',
      r'\toprule file & routine & lines & purpose \\ \midrule \endhead']
for f in ('charm_mod.f', 'coal_new.f'):
    lines = (ROOT / 'src_new' / f).read_text(errors='replace').splitlines()
    starts = [(i + 1, m.group(2)) for i, l in enumerate(lines)
              if (m := re.match(r'^\s+((?:double precision |real |integer )?(?:SUBROUTINE|FUNCTION|subroutine|function|BLOCK DATA))\s+(\w+)', l, re.I))]
    for k, (a, nm) in enumerate(starts):
        b = (starts[k + 1][0] - 1) if k + 1 < len(starts) else len(lines)
        com = next((re.sub(r'^[cC*]\s*', '', l).strip() for l in lines[a:a + 15]
                    if re.match(r'^[cC*]\s+\S', l) and not l.strip('c*= -').startswith('---')), '')
        L.append(rf'\texttt{{{tex(f, False)}}} & \texttt{{{tex(nm)}}} & {a}--{b} & {tex(com[:120])} \\')
L += [r'\bottomrule \end{longtable}}', '']
# ------------------------------------------------------------------ table 3: switches
ref = dict(re.findall(r'(AMPT_\w+)=([^\s"]+)', (ROOT / 'configs/reference_auau.env').read_text()))
refpp = dict(re.findall(r'(AMPT_\w+)=([^\s"]+)', (ROOT / 'configs/reference_pp.env').read_text()))
L += [r'{\small\begin{longtable}{@{}p{4.2cm}p{3.4cm}p{2.8cm}p{2.8cm}@{}}',
      r'\caption{Environment switches: the source line that reads each one, and its value in the reference Au+Au and p+p '
      r'configurations (\texttt{configs/reference\_*.env}). Unset = the stock behaviour.}\label{tab:switches}\\',
      r'\toprule switch & read at & reference Au+Au & reference p+p \\ \midrule \endfirsthead',
      r'\toprule switch & read at & reference Au+Au & reference p+p \\ \midrule \endhead']
seen = set()
for f in ('charm_mod.f', 'coal_new.f'):
    for i, l in enumerate((ROOT / 'src_new' / f).read_text(errors='replace').splitlines()):
        for v in re.findall(r"'(AMPT_\w+)'", l):
            if v not in seen:
                seen.add(v)
                L.append(rf'\texttt{{{tex(v)}}} & \texttt{{{tex(f)}:{i + 1}}} & {tex(ref.get(v, "unset"))} & {tex(refpp.get(v, "unset"))} \\')
L += [r'\bottomrule \end{longtable}}']
open(OUT, 'w').write('\n'.join(L) + '\n')
print(f'{len(hunks)} hunks, {len(seen)} switches -> {OUT}')
