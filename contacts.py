#!/usr/bin/env python3
"""
Residue-residue contacts using MDTraj. Contacts are computed using a cutoff and normalized 
by the allowed residue–residue pairs. Intra-chain contacts exclude bonded and 1–3 neighbor pairs.

Usage:
  python script.py -f traj.dcd -s top.pdb -b 0 -e 1000 -dt 1 -o AA AB BB -c 5.0
"""

import argparse, sys, os
import numpy as np, mdtraj as md
from tqdm import tqdm
from collections import defaultdict
from queue import SimpleQueue
from datetime import datetime
import concurrent.futures

def ga():
    p = argparse.ArgumentParser(
        description="Residue contacts with MDTraj.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("-f", "--trajectory", required=True, metavar="", help="Trajectory file (e.g., DCD)")
    p.add_argument("-s", "--topology", required=True, metavar="", help="Topology file (PDB/PSF)")
    p.add_argument("-b", "--begin", type=float, default=0.0, metavar="", help="Begin time (ps)")
    p.add_argument("-e", "--end", type=float, default=None, metavar="", help="End time (ps)")
    p.add_argument("-dt", "--skip", type=int, default=1, metavar="", help="Process every nth frame")
    p.add_argument("-o", "--output", nargs="*", metavar="", help="Chain pair specifiers (e.g., AA, AB)")
    p.add_argument("-c", "--cutoff", type=float, default=5.0, metavar="", help="Cutoff in Å")
    return p.parse_args()

def parse_pdb(f):
    info = []
    try:
        with open(f, "r") as fp:
            for l in fp:
                if l.startswith("ATOM") or l.startswith("HETATM"):
                    ch = l[21].strip()
                    seg = l[72:76].strip()
                    info.append((seg, ch))
    except Exception as e:
        sys.exit(f"Error reading {f}: {e}")
    return info

def build_cs(info):
    cs = defaultdict(lambda: defaultdict(list))
    for i, (seg, ch) in enumerate(info):
        cs[ch][seg].append(i)
    return cs

def precomp_ex_seg(cs, topo, info):
    seg = defaultdict(list)
    for i, (sid, _) in enumerate(info):
        seg[sid].append(i)
    ex = {}
    for sid, idxs in seg.items():
        es = {i: {i} for i in idxs}
        for b in topo.bonds:
            i1, i2 = b.atom1.index, b.atom2.index
            if i1 in es and i2 in es:
                es[i1].add(i2)
                es[i2].add(i1)
        for i in idxs:
            vis = {i}
            q = SimpleQueue()
            q.put((i, 0))
            while not q.empty():
                cur, d = q.get()
                if d >= 2: continue
                for nb in es[cur]:
                    if nb not in vis:
                        vis.add(nb)
                        es[i].add(nb)
                        q.put((nb, d+1))
        ex[sid] = es
    return ex, seg

def precomp_allowed(seg, ex):
    ap = {}
    for sid, idxs in seg.items():
        n = len(idxs)
        if n < 2:
            ap[sid] = (np.empty(0, int), np.empty(0, int))
            continue
        ii, jj = np.triu_indices(n, k=1)
        allow = []
        for k in range(len(ii)):
            gi = idxs[ii[k]]
            gj = idxs[jj[k]]
            allow.append(gj not in ex[sid].get(gi, set()))
        allow = np.array(allow)
        ap[sid] = (ii[allow], jj[allow])
    return ap

def cci(p1, p2, cut):
    d = np.linalg.norm(p1[:, None, :] - p2[None, :, :], axis=-1)
    nc = np.count_nonzero(d < cut)
    tot = p1.shape[0] * p2.shape[0]
    return nc, tot

def cci_intra(idxs, ap_idx, pos, cut):
    ii, jj = ap_idx
    if len(ii) == 0: return 0, 0
    d = np.linalg.norm(pos[ii] - pos[jj], axis=1)
    nc = np.count_nonzero(d < cut)
    tot = len(ii)
    return nc, tot

def proc_fr(fr, cs, ex, ap, cut):
    res = {}
    chs = sorted(cs.keys())
    for i, ch1 in enumerate(chs):
        for ch2 in chs[i:]:
            lab = f"{min(ch1,ch2)}-{max(ch1,ch2)}"
            cnt, tot = 0, 0
            if ch1 == ch2:
                sids = list(cs[ch1].keys())
                for sid in sids:
                    idxs = cs[ch1][sid]
                    pos_seg = fr[idxs]
                    n, t = cci_intra(idxs, ap[sid], pos_seg, cut)
                    cnt += n; tot += t
                for i1 in range(len(sids)):
                    for i2 in range(i1+1, len(sids)):
                        pos1 = fr[cs[ch1][sids[i1]]]
                        pos2 = fr[cs[ch1][sids[i2]]]
                        n, t = cci(pos1, pos2, cut)
                        cnt += n; tot += t
            else:
                idx1 = []; idx2 = []
                for lst in cs[ch1].values(): idx1.extend(lst)
                for lst in cs[ch2].values(): idx2.extend(lst)
                p1 = fr[idx1]; p2 = fr[idx2]
                n, t = cci(p1, p2, cut)
                cnt += n; tot += t
            res[lab] = cnt/tot if tot > 0 else 0
    return res

def bkup_file(f):
    if os.path.exists(f):
        base, ext = os.path.splitext(f)
        ts = datetime.now().strftime("%d-%m-%Y_%H-%M-%S")
        os.rename(f, f"{base}_{ts}{ext}")

def main():
    arg = ga()
    try:
        traj = md.load(arg.trajectory, top=arg.topology)
    except Exception as e:
        sys.exit(f"Error loading files: {e}")
    info = parse_pdb(arg.topology)
    if len(info) != traj.n_atoms:
        sys.exit("Error: Atom count mismatch between PDB and trajectory.")
    cs = build_cs(info)
    if arg.output:
        req = set()
        for spec in arg.output:
            if len(spec) != 2:
                sys.exit(f"Invalid spec '{spec}'; must be 2 chain IDs.")
            a, b = spec[0], spec[1]
            if a not in cs or b not in cs:
                sys.exit(f"Chain '{a if a not in cs else b}' not found.")
            req.add(f"{min(a,b)}-{max(a,b)}")
        pairs = sorted(req)
    else:
        chs = sorted(cs.keys())
        pairs = []
        for i, a in enumerate(chs):
            for b in chs[i:]:
                pairs.append(f"{a}-{b}")
    ex, seg = precomp_ex_seg(cs, traj.topology, info)
    ap = precomp_allowed(seg, ex)
    t = traj.time if traj.time is not None and len(traj.time)==traj.n_frames else np.arange(traj.n_frames, float)
    s_idx = np.searchsorted(t, arg.begin, side='left')
    e_idx = np.searchsorted(t, arg.end, side='right') if arg.end is not None else traj.n_frames
    f_idxs = list(range(s_idx, e_idx, arg.skip))
    
    res_list = []
    def worker(i):
        ti = t[i]
        fr = traj.xyz[i]
        r = proc_fr(fr, cs, ex, ap, arg.cutoff)
        return (ti, r)
    
    with concurrent.futures.ThreadPoolExecutor() as executor:
        for r in tqdm(executor.map(worker, f_idxs), total=len(f_idxs), desc="Processing frames"):
            res_list.append(r)
    
    out = {p: [] for p in pairs}
    for ti, r in res_list:
        for p in pairs:
            if p in r:
                out[p].append(f"{ti:.3f} {r[p]:.6f}\n")
    
    for p in pairs:
        fname = f"contacts_{p}.dat"
        bkup_file(fname)
        try:
            with open(fname, "w") as f:
                f.writelines(out[p])
        except Exception as e:
            sys.exit(f"Error writing to {fname}: {e}")

if __name__ == "__main__":
    main()
