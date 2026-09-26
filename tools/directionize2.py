"""
C003 有向化（改訂版）: パラメータ q による統一族
  各無向辺を確率 q で双方向化、1-q でランダム一方向化。
  補正: T が単一 SCC に収まるまで、一方向化した辺をランダム順に双方向化（逆辺追加）。
        補正辺は必ず元の無向辺の逆向きなので、有向辺集合 ⊆ 無向辺集合の向き付き版に収まる。
  最小プレフィックス長は二分探索で求める（辺追加で強連結性は壊れないので単調）。
  案A = q=1, 案B = q=0(+補正), 案C = q=0.5(+補正)
"""
import random, sys
from typing import List, Tuple
import networkx as nx
sys.path.insert(0, "/home/claude/scss")
from stp_parse import STPInstance

def _T_in_single_scc(n, dedges, terminals):
    G = nx.DiGraph(); G.add_nodes_from(range(n)); G.add_edges_from(dedges)
    t0 = terminals[0]
    for comp in nx.strongly_connected_components(G):
        if t0 in comp:
            return all(t in comp for t in terminals)
    return False

def directionize_q(inst: STPInstance, q: float, rng: random.Random):
    """戻り値: (dedges, n_bidir_by_q, n_bidir_by_fix, n_oneway_final)"""
    n = inst.n
    base = []      # 確定した有向辺
    oneway = []    # 一方向化した辺 (u->v)。補正候補。
    nb_q = 0
    for u, v in inst.edges:
        if rng.random() < q:
            base.append((u, v)); base.append((v, u)); nb_q += 1
        else:
            if rng.random() < 0.5: oneway.append((u, v))
            else: oneway.append((v, u))
    dedges = base + oneway
    if _T_in_single_scc(n, dedges, inst.terminals):
        return dedges, nb_q, 0, len(oneway)
    # 補正: oneway をシャッフルし、逆辺を先頭 k 本追加する最小 k を二分探索
    order = oneway[:]; rng.shuffle(order)
    rev = [(v, u) for u, v in order]
    lo, hi = 1, len(rev)
    if not _T_in_single_scc(n, dedges + rev, inst.terminals):
        # 全部双方向化しても不可 = 無向グラフで T が非連結（ベンチマークでは起きないはず）
        return dedges + rev, nb_q, len(rev), 0
    while lo < hi:
        mid = (lo + hi) // 2
        if _T_in_single_scc(n, dedges + rev[:mid], inst.terminals): hi = mid
        else: lo = mid + 1
    return dedges + rev[:lo], nb_q, lo, len(oneway) - lo
