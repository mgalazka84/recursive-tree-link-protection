// Reproduction program for The value of topology information for link protection.
// Compile with a C++17 compiler. See README.md for the experimental modes.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <queue>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using Parents = std::vector<int>;
using Values = std::vector<double>;

struct Moments {
    int count = 0;
    double mean = 0, sum_sq = 0;
    void add(double x) {
        ++count;
        const double delta = x - mean;
        mean += delta / count;
        sum_sq += delta * (x - mean);
    }
    double se() const { return std::sqrt(sum_sq / (count - 1) / count); }
};

// Uniform 32-bit integer in [0,m). Multiplication with rejection avoids modulo bias.
uint32_t below(std::mt19937& rng, uint32_t m) {
    uint64_t product = uint64_t(rng()) * m;
    if (uint32_t(product) < m) {
        const uint32_t threshold = -m % m;
        while (uint32_t(product) < threshold) product = uint64_t(rng()) * m;
    }
    return uint32_t(product >> 32);
}

double unit(std::mt19937& rng) {
    const uint64_t high = rng() >> 6;
    const uint64_t low = rng() >> 6;
    const uint64_t bits = (high << 26) | low;
    return (double(bits) + 0.5) * 0x1p-52;
}

void random_tree(Parents& parent, std::mt19937& rng) {
    for (int i = 1; i < int(parent.size()); ++i) parent[i] = int(below(rng, i));
}

Values weights(const Parents& parent, double p) {
    Values w(parent.size(), 1);
    for (int i = int(w.size()) - 1; i > 0; --i) w[parent[i]] += p * w[i];
    return w;
}

struct Gains { Values age, map; };
Gains solve(const Parents& parent, double p, int budget) {
    const int n = int(parent.size());
    const double q = 1 - p;
    const Values w = weights(parent, p);
    Gains gains{Values(budget + 1), Values(budget + 1)};
    for (int j = 1; j <= budget; ++j) gains.age[j] = gains.age[j-1] + q*w[j];
    if (budget == 1) {
        double best = 0;
        for (int i = 1; i < n; ++i) if (parent[i] == 0) best = std::max(best, w[i]);
        gains.map[1] = q*best;
        return gains;
    }

    std::vector<int> size(n,1), first(n,-1), sibling(n,-1);
    for (int i = n-1; i > 0; --i) {
        size[parent[i]] += size[i];
        sibling[i] = first[parent[i]];
        first[parent[i]] = i;
    }
    std::vector<int> order, stack{0};
    order.reserve(n);
    while (!stack.empty()) {
        const int v = stack.back();
        stack.pop_back();
        order.push_back(v);
        for (int j = first[v]; j != -1; j = sibling[j]) stack.push_back(j);
    }
    Values prev(n+1), cur(n+1);
    for (int j = 1; j <= budget; ++j) {
        for (int i = n-1; i >= 1; --i) {
            const int v = order[i];
            cur[i] = std::max(cur[i+size[v]], w[v] + prev[i+1]);
        }
        gains.map[j] = q*cur[1];
        cur.swap(prev);
    }
    return gains;
}

double age_exact(int n, int k, double q) {
    double log_base = 0, log_prefix = 0;
    for (int i = 2; i <= n; ++i) {
        const double term = std::log1p(-q/i);
        log_base += term;
        if (i <= k+1) log_prefix += term;
    }
    return n * std::exp(log_base) * std::expm1(-log_prefix);
}

void finite_case(int n, double c, int samples, uint32_t seed, const std::vector<int>& budgets) {
    const double q = c/std::log(n), p = 1-q;
    std::mt19937 rng(seed);
    Parents parent(n);
    const int cap = budgets.back();
    std::vector<Moments> map(cap+1), age(cap+1);
    for (int s = 0; s < samples; ++s) {
        random_tree(parent,rng);
        const Gains g = solve(parent,p,cap);
        for (int k : budgets) { map[k].add(g.map[k]); age[k].add(g.age[k]); }
    }
    for (int k : budgets) {
        const double a = age_exact(n,k,q);
        std::cout << n << ',' << c << ',' << k << ',' << samples << ',' << seed << ','
                  << a << ',' << age[k].mean << ',' << map[k].mean << ',' << map[k].se()
                  << ',' << map[k].mean/a << ',' << map[k].se()/a << '\n';
    }
}

void finite_runs(bool pilot) {
    std::cout << "n,c,k,trees,seed,age_exact,age_sample,map_mean,map_se,ratio,ratio_se\n";
    if (pilot) { finite_case(10000,1,40000,20261109,{1}); return; }
    const std::vector<int> k{1,2,4,8,16,32,64};
    finite_case(128,1,8000,20261110,k);
    finite_case(1024,1,5000,20261111,k);
    finite_case(8192,1,1200,20261112,k);
    finite_case(1024,0.5,2500,20261113,k);
    finite_case(1024,2,2500,20261114,k);
}

void preferential_tree(Parents& parent, std::mt19937& rng, std::vector<int>& tickets) {
    // A site owns one ticket for itself and one for each child.
    tickets.clear();
    tickets.push_back(0);
    for (int child = 1; child < int(parent.size()); ++child) {
        const int chosen = tickets[below(rng, uint32_t(tickets.size()))];
        parent[child] = chosen;
        tickets.push_back(chosen);
        tickets.push_back(child);
    }
}

double expected_preferential_oldest(int n, double p) {
    double expected_weight = 1;
    for (int t = 2; t < n; ++t)
        expected_weight += ((1 + p) * expected_weight - 1) / (2*t - 1);
    return (1-p) * expected_weight;
}

void sensitivity_case(int n, int samples, bool preferential, uint32_t seed) {
    const Values failure{0.1, 0.2, 0.5, 1/std::log(n)};
    const std::vector<std::string> regimes{"fixed", "fixed", "fixed", "critical"};
    std::vector<Moments> full(failure.size()), oldest(failure.size());
    std::vector<Moments> partial(failure.size()), residual(failure.size());
    Parents parent(n);
    std::vector<int> sizes(n), branches, tickets;
    if (preferential) tickets.reserve(2*n-1);
    std::mt19937 rng(seed);

    for (int draw = 0; draw < samples; ++draw) {
        if (preferential) preferential_tree(parent,rng,tickets);
        else random_tree(parent,rng);
        std::fill(sizes.begin(), sizes.end(), 1);
        for (int j = n-1; j >= 1; --j) sizes[parent[j]] += sizes[j];
        branches.clear();
        for (int j = 1; j < n; ++j)
            if (parent[j] == 0) branches.push_back(j);
        // Increasing label order selects the oldest in a tie.
        const int largest = *std::max_element(branches.begin(), branches.end(),
                             [&](int i,int j){ return sizes[i] < sizes[j]; });
        for (size_t r = 0; r < failure.size(); ++r) {
            const double q = failure[r];
            const Values wt = weights(parent, 1-q);
            double maximum = 0;
            for (int j : branches) maximum = std::max(maximum,wt[j]);
            full[r].add(q*maximum);
            oldest[r].add(q*wt[1]);
            partial[r].add(q*wt[largest]);
            residual[r].add(q*(maximum-wt[largest]));
        }
    }
    for (size_t r = 0; r < failure.size(); ++r) {
        const double q = failure[r], p = 1-q;
        const double a = preferential ? expected_preferential_oldest(n,p) : age_exact(n,1,q);
        std::cout << (preferential ? "preferential" : "uniform") << ',' << regimes[r] << ','
                  << n << ',' << q << ',' << samples << ',' << seed << ',' << a << ','
                  << oldest[r].mean << ',' << oldest[r].se() << ',' << full[r].mean << ','
                  << full[r].se() << ',' << full[r].mean/a << ',' << full[r].se()/a << ','
                  << partial[r].mean << ',' << residual[r].mean << ',' << residual[r].se() << '\n';
    }
}

void sensitivity_runs() {
    std::cout << "growth,regime,n,q,trees,seed,age_exact,age_sample,age_se,map_mean,map_se,"
                 "ratio,ratio_se,size_policy_mean,map_minus_size,paired_gap_se\n";
    const std::vector<std::pair<int,int>> sizes{{128,15000},{1024,10000},{8192,5000},{65536,2500}};
    for (int model = 0; model <= 1; ++model) {
        for (size_t i = 0; i < sizes.size(); ++i)
            sensitivity_case(sizes[i].first, sizes[i].second, model!=0, 711901+20*model+i);
    }
}

void limiting_runs() {
    constexpr int samples = 20000, cap = 1024;
    constexpr uint32_t seed = 20261115;
    std::mt19937 rng(seed);
    std::vector<Moments> top(cap+1);
    for (int s = 0; s < samples; ++s) {
        // A true entry is a vertex, a false entry is unexplored sibling mass.
        std::priority_queue<std::pair<double,bool>> heap;
        heap.emplace(1.0,false);
        double sum = 0;
        int recorded = 0;
        while (recorded < cap) {
            const auto [mass, vertex] = heap.top();
            heap.pop();
            if (vertex) {
                sum += mass;
                ++recorded;
                if ((recorded & (recorded-1)) == 0) top[recorded].add(sum);
                if (recorded == cap) break;
            }
            const double left = mass * unit(rng);
            heap.emplace(left,true);
            heap.emplace(mass-left,false);
        }
    }
    std::cout << "k,trees,seed,age_constant,map_constant,map_se,ratio,ratio_se,gap,jensen_gap\n";
    double age = 0;
    for (int k = 1; k <= cap; ++k) {
        age += 1.0/(k+1);
        if ((k & (k-1)) != 0) continue;
        const auto& v = top[k];
        std::cout << k << ',' << samples << ',' << seed << ',' << age << ',' << v.mean << ','
                  << v.se() << ',' << v.mean/age << ',' << v.se()/age << ','
                  << v.mean-age << ',' << std::log1p(k)-v.mean << '\n';
    }
}

int bitcount(int x) {
    int k = 0;
    while (x) { x &= x-1; ++k; }
    return k;
}

double direct_reach(const Parents& par, int mask, double p) {
    Values reachable(par.size(),1);
    for (int v = 1; v < int(par.size()); ++v)
        reachable[v] = reachable[par[v]] * ((mask & (1 << (v-1))) ? 1 : p);
    return std::accumulate(reachable.begin(),reachable.end(),0.0);
}

Values tails(int n, int mask, double p) {
    Values law(n+1);
    law[1] = 1;
    for (int t=1; t<n; ++t) {
        Values next(n+1);
        const double r = (mask & (1 << (t-1))) ? 1 : p;
        for (int x=1; x<=t; ++x) {
            const double jump = r*x/t;
            next[x] += law[x]*(1-jump);
            next[x+1] += law[x]*jump;
        }
        next.swap(law);
    }
    for (int x=n-1; x>0; --x) law[x] += law[x+1];
    return law;
}

void verify() {
    const Values probabilities{0.25,0.5,0.8};
    long long trees = 0, a_sets = 0, comparisons = 0;
    double error = 0, violation = 0, exact_error = 0;
    for (int n=2; n<=8; ++n) {
        Parents par(n);
        std::vector<Values> totals(probabilities.size(),Values(n));
        int count = 0;
        std::function<void(int)> enumerate = [&](int i) {
            if (i<n) {
                for (int parent=0; parent<i; ++parent) { par[i]=parent; enumerate(i+1); }
                return;
            }
            ++trees; ++count;
            for (size_t idx=0; idx<probabilities.size(); ++idx) {
                const double p = probabilities[idx], base=direct_reach(par,0,p);
                const Gains g=solve(par,p,n-1);
                Values brute(n);
                for (int mask=0; mask<(1<<(n-1)); ++mask) {
                    ++a_sets;
                    const int k=bitcount(mask);
                    const double gain=direct_reach(par,mask,p)-base;
                    brute[k]=std::max(brute[k],gain);
                    if (mask == (1<<k)-1) error=std::max(error,std::abs(gain-g.age[k]));
                }
                const Gains single = solve(par,p,1);
                error=std::max(error,std::abs(single.map[1]-brute[1]));
                for (int k=0; k<n; ++k) {
                    error=std::max(error,std::abs(brute[k]-g.map[k]));
                    totals[idx][k]+=g.age[k];
                }
            }
        };
        enumerate(1);
        for (size_t idx=0; idx<probabilities.size(); ++idx) for (int k=0;k<n;++k)
            exact_error=std::max(exact_error,std::abs(totals[idx][k]/count-age_exact(n,k,1-probabilities[idx])));
    }
    for (int n=2; n<=10; ++n) for (double p : probabilities) {
        std::vector<Values> early(n);
        for (int k=0; k<n; ++k) early[k]=tails(n,(1<<k)-1,p);
        for (int mask=0; mask<(1<<(n-1)); ++mask) {
            const Values law=tails(n,mask,p);
            for (int t=1; t<=n; ++t) {
                ++comparisons;
                violation=std::max(violation,law[t]-early[bitcount(mask)][t]);
            }
        }
    }
    std::cout << "recursive_trees=" << trees << "\nprobability_tree_subset_evaluations=" << a_sets
              << "\ntail_comparisons=" << comparisons << "\nmax_direct_error=" << error
              << "\nmax_age_expectation_error=" << exact_error << "\nmax_tail_violation=" << violation << '\n';
    if (std::max({error,exact_error,violation})>1e-9) throw std::runtime_error("Verification failed");
}

int main(int argc, char** argv) {
    std::cout << std::setprecision(13);
    const std::string mode = argc==2 ? argv[1] : "";
    if (mode=="verify") verify();
    else if (mode=="one") finite_runs(true);
    else if (mode=="finite") finite_runs(false);
    else if (mode=="limit") limiting_runs();
    else if (mode=="sensitivity") sensitivity_runs();
    else { std::cerr << "Usage: experiments {verify|one|finite|limit|sensitivity}\n"; return 2; }
}
