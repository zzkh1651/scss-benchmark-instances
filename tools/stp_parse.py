"""
STP (SteinLib) format parser for SCSS benchmark instances.

学士論文の gen_er.py / 既存 stp_parse.py（共同研究者環境には未配置）とは独立に、
今回の有向化検討用として実装したもの。0-indexed に変換して返す。

想定フォーマット:
    SECTION Graph / Section Graph
    Nodes <n>
    Edges <m>
    E u v [weight]
    END / End
    SECTION Terminals / Section Terminals
    Terminals <p>
    T <t>
    END / End

大文字・小文字の揺れ（281件中 231件が "SECTION", 残りは "Section" 系）に対応するため
キーワード比較はすべて大文字化して行う。
"""
from dataclasses import dataclass
from typing import List, Tuple


@dataclass
class STPInstance:
    name: str
    n: int
    edges: List[Tuple[int, int]]  # 0-indexed, 無向辺として (u, v), u<v は保証しない
    terminals: List[int]          # 0-indexed


def parse_stp(path: str) -> STPInstance:
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()

    n = None
    m_declared = None
    edges = []
    terminals = []

    mode = None  # "graph" or "terminals" or None
    for raw in lines:
        line = raw.strip()
        if not line:
            continue
        upper = line.upper()

        if upper.startswith("SECTION GRAPH"):
            mode = "graph"
            continue
        if upper.startswith("SECTION TERMINALS"):
            mode = "terminals"
            continue
        if upper in ("END", "EOF") or upper.startswith("END"):
            # "SECTION" 内の END でモード解除。EOF はファイル末尾。
            if upper == "EOF":
                break
            mode = None
            continue

        if mode == "graph":
            if upper.startswith("NODES"):
                n = int(line.split()[1])
            elif upper.startswith("EDGES"):
                m_declared = int(line.split()[1])
            elif upper.startswith("E "):
                parts = line.split()
                u, v = int(parts[1]), int(parts[2])
                edges.append((u - 1, v - 1))  # 0-indexed
        elif mode == "terminals":
            if upper.startswith("T "):
                t = int(line.split()[1])
                terminals.append(t - 1)  # 0-indexed

    if n is None:
        raise ValueError(f"{path}: Nodes not found")
    if not terminals:
        raise ValueError(f"{path}: Terminals not found")

    name = path.split("/")[-1]
    return STPInstance(name=name, n=n, edges=edges, terminals=terminals)


if __name__ == "__main__":
    import sys
    inst = parse_stp(sys.argv[1])
    print(f"name={inst.name} n={inst.n} m={len(inst.edges)} p={len(inst.terminals)}")
