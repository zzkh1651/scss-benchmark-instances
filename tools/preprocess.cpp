// ============================================================
//  SSD Algorithm for SCSS Problem
//  + Lagrangian Relaxation (Subgradient Method) as lower bound
//  + LP Relaxation (CPLEX) for comparison
//
//  Problem:  p-Strongly Connected Steiner Subgraph (SCSS)
//  Algorithm: Branch-and-Bound (Algorithm 3 / 4) with
//             Lagrangian & LP lower bound comparison
// ============================================================

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <queue>
#include <stack>
#include <algorithm>
#include <numeric>
#include <chrono>
#include <iomanip>
#include <filesystem>
#include <tuple>
#include <cmath>
#include <limits>
#include <cassert>
#include <cstring>
#include <unordered_map>

// CPLEXライブラリがある場合はマクロで有効化
#ifdef USE_CPLEX
#include <ilcplex/ilocplex.h>
#endif

namespace fs = std::filesystem;

// ============================================================
//  Bitset  (dynamic, std::vector<bool>-based)
// ============================================================
using Bitset = std::vector<bool>;

inline void  bs_set  (Bitset& b, int i)              { b[i] = true; }
inline int   bs_count(const Bitset& b)               { int c=0; for(bool v:b) c+=v; return c; }
inline bool  bs_empty(const Bitset& b)               { for(bool v:b) if(v) return false; return true; }
inline bool  bs_subset(const Bitset& sub,const Bitset& sup){
    for(int i=0;i<(int)sub.size();++i) if(sub[i]&&!sup[i]) return false;
    return true; }
inline Bitset bs_and  (const Bitset& a,const Bitset& b){
    Bitset r(a.size()); for(int i=0;i<(int)a.size();++i) r[i]=a[i]&&b[i]; return r; }
inline Bitset bs_or   (const Bitset& a,const Bitset& b){
    Bitset r(a.size()); for(int i=0;i<(int)a.size();++i) r[i]=a[i]||b[i]; return r; }
inline Bitset bs_minus(const Bitset& a,const Bitset& b){
    Bitset r(a.size()); for(int i=0;i<(int)a.size();++i) r[i]=a[i]&&!b[i]; return r; }
inline bool  bs_equal (const Bitset& a,const Bitset& b){ return a==b; }

// ============================================================
//  Dominator Tree  (Lengauer-Tarjan, dynamic N)
// ============================================================
struct DominatorTree {
    int N;
    std::vector<std::vector<int>> g, rg;
    std::vector<int> idom, semi, par;
    std::vector<int> ord;

    struct UF {
        std::vector<int> par, best;
        const std::vector<int>& semi;
        UF(int n,const std::vector<int>& s):par(n),best(n),semi(s){
            std::iota(par.begin(),par.end(),0);
            std::iota(best.begin(),best.end(),0); }
        int find(int v){
            if(par[v]==v) return v;
            int r=find(par[v]);
            if(semi[best[v]]>semi[best[par[v]]]) best[v]=best[par[v]];
            return par[v]=r; }
        int eval(int v){ find(v); return best[v]; }
        void link(int p,int c){ par[c]=p; }
    };

    explicit DominatorTree(int n):N(n),g(n),rg(n),idom(n,-1),semi(n,-1),par(n,0){ ord.reserve(n); }
    void add_edge(int u,int v){ g[u].push_back(v); }

    void dfs(int v){
        semi[v]=(int)ord.size(); ord.push_back(v);
        for(int w:g[v]) if(semi[w]<0){ dfs(w); par[w]=v; } }

    void build(int root){
        for(int u=0;u<N;++u) for(int v:g[u]) rg[v].push_back(u);
        semi.assign(N,-1); ord.clear();
        dfs(root);
        UF uf(N,semi);
        std::vector<std::vector<int>> bucket(N);
        std::vector<int> U(N,0);
        for(int i=(int)ord.size()-1;i>=1;--i){
            int x=ord[i];
            for(int v:rg[x]){ if(semi[v]<0) continue; int u=uf.eval(v); if(semi[x]>semi[u]) semi[x]=semi[u]; }
            bucket[ord[semi[x]]].push_back(x);
            for(int v:bucket[par[x]]) U[v]=uf.eval(v);
            bucket[par[x]].clear();
            uf.link(par[x],x);
        }
        for(int i=1;i<(int)ord.size();++i){ int x=ord[i],u=U[x]; idom[x]=(semi[x]==semi[u])?semi[x]:idom[u]; }
        for(int i=1;i<(int)ord.size();++i){ int x=ord[i]; idom[x]=ord[idom[x]]; }
        idom[root]=root;
    }
    int operator[](int k) const { return idom[k]; }
};

// ============================================================
//  Graph utilities
// ============================================================
constexpr int INF = 1000000000;

void distbfs(const Bitset& src, const Bitset& S,
             const std::vector<std::vector<int>>& adj,
             std::vector<int>& dist)
{
    int N=(int)S.size(); dist.assign(N,INF);
    std::queue<int> q;
    for(int i=0;i<N;++i) if(src[i]){ dist[i]=0; q.push(i); }
    while(!q.empty()){ int u=q.front();q.pop();
        for(int v:adj[u]) if(S[v]&&dist[v]==INF){ dist[v]=dist[u]+1; q.push(v); } }
}

void reachbfs(const Bitset& src, const Bitset& S,
              const std::vector<std::vector<int>>& adj,
              int theta, int& counter, std::vector<int>& clist)
{
    counter++; if(counter>=INF){ counter=1; std::fill(clist.begin(),clist.end(),0); }
    std::queue<std::pair<int,int>> q;
    for(int i=0;i<(int)S.size();++i) if(src[i]){ clist[i]=counter; q.push({i,0}); }
    while(!q.empty()){ auto[u,d]=q.front();q.pop();
        for(int v:adj[u]) if(S[v]&&clist[v]!=counter){
            clist[v]=counter; if(d+2<theta) q.push({v,d+1}); } }
}

Bitset findSCC(const std::vector<std::vector<int>>& g,
               const std::vector<std::vector<int>>& rg, int start)
{
    int N=(int)g.size();
    std::stack<int> fin; std::vector<bool> vis(N,false);
    for(int s=0;s<N;++s){
        if(vis[s]) continue;
        std::stack<std::pair<int,int>> dfs;
        dfs.push({s,0});
        while(!dfs.empty()){
            auto&[v,idx]=dfs.top(); if(!vis[v]) vis[v]=true;
            if(idx<(int)g[v].size()){ int w=g[v][idx++]; if(!vis[w]) dfs.push({w,0}); }
            else{ fin.push(v); dfs.pop(); } } }
    std::fill(vis.begin(),vis.end(),false);
    Bitset result(N,false);
    while(!fin.empty()){
        int v=fin.top();fin.pop(); if(vis[v]) continue;
        Bitset scc(N,false); std::stack<int> s2;
        s2.push(v); vis[v]=true; scc[v]=true;
        while(!s2.empty()){ int u=s2.top();s2.pop();
            for(int w:rg[u]) if(!vis[w]){ vis[w]=true; scc[w]=true; s2.push(w); } }
        if(scc[start]){ result=scc; break; } }
    return result;
}

// ============================================================
//  Preprocess-only driver (C003 / C005 用)
//  入力: BB 形式 (V E T / u v ×E / terminal list)
//  出力: 1 行  "<file> V E T feasible |C| |T_hat| time_ms"
//  ロジックは bb_cur の solve() 冒頭 (findSCC + 支配木 T_hat) と同一
// ============================================================
int main(int argc,char* argv[]){
    if(argc<2){ std::cerr<<"usage: preprocess <bb_format_file>...\n"; return 1; }
    for(int a=1;a<argc;++a){
        std::ifstream in(argv[a]);
        if(!in.is_open()){ std::cout<<argv[a]<<" ERROR open\n"; continue; }
        int V,E,T; in>>V>>E>>T;
        std::vector<std::vector<int>> adj(V),rg(V);
        for(int i=0;i<E;++i){ int u,v; in>>u>>v; adj[u].push_back(v); rg[v].push_back(u); }
        in.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
        std::vector<int> terminal; std::string line;
        while(std::getline(in,line)) if(!line.empty()) break;
        { std::istringstream ss(line); int t; while(ss>>t) terminal.push_back(t); }
        auto start=std::chrono::steady_clock::now();
        Bitset bitC=findSCC(adj,rg,terminal[0]);
        bool feasible=true;
        for(int t:terminal) if(!bitC[t]){ feasible=false; break; }
        int Csz=bs_count(bitC), Tsz=-1;
        if(feasible){
            Bitset bitT(V,false);
            for(int t:terminal) bs_set(bitT,t);
            for(int u:terminal){
                DominatorTree u_DT(V);
                for(int x=0;x<V;++x) if(bitC[x]) for(int y:adj[x]) if(bitC[y]) u_DT.add_edge(x,y);
                u_DT.build(u);
                for(int v:terminal){ if(v==u) continue; int cur=v;
                    while(cur!=u&&cur!=-1&&cur!=u_DT[cur]){ cur=u_DT[cur]; if(cur!=-1) bs_set(bitT,cur); } }
            }
            Tsz=bs_count(bitT);
        }
        double ms=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-start).count();
        std::cout<<argv[a]<<" "<<V<<" "<<E<<" "<<terminal.size()<<" "<<(feasible?1:0)<<" "<<Csz<<" "<<Tsz<<" "<<ms<<"\n";
    }
    return 0;
}
