#pragma once
#include "Scheduler.hpp"
class MQTTScheduler : private Scheduler
{
private:
	multiset<pair<ll, VLL_WRAPPER>> J;
public:
	MQTTScheduler()
	{
		J = multiset<pair<ll, VLL_WRAPPER>>();
	}
	void insert(ll due, ll process, ll id)
	{
		J.insert({ due, VLL_WRAPPER(process, id)});
	}
	void clear()
	{
		J.clear();
	}
	VLL_WRAPPER run()
	{
		vector<pair<pair<ll, ll>, pair<ll, multiset<VLL_WRAPPER>>>> X;
		multiset<VLL_WRAPPER> T;
		T.insert(VLL_WRAPPER(0LL, 0LL));
		ll P = 0LL, dsharp, count = 1LL, sum = 0LL, cur = 0LL, k = 0LL, Pi, size;
		bool next = false;
		for (auto i = J.begin(); i != J.end(); i++)
		{
			auto j = i;
			if (i == J.begin() || (*i).first != (*(--j)).first)
			{
				X.push_back({ {(*i).first, 0LL}, {-1LL, multiset<VLL_WRAPPER>()} });
				X[X.size() - 1LL].second.second.insert(VLL_WRAPPER(0, 0));
			}
			X[X.size() - 1LL].second.second.insert((*i).second);
			X[X.size() - 1LL].first.second += ll((*i).second);
			P += ll((*i).second);
		}
		dsharp = X.size();
		for (ll i = 0LL; i < dsharp; i++)
		{
			if (X[i].first.second > (ll)(pow(P, 0.6L)))
			{
				X[i].second.first = 0LL;
				sum = 0LL;
				if (next) count++;
				continue;
			}
			next = true;
			if ((X[i].second.first = count) && (sum = sum + X[i].first.second) > (ll)(pow(P, 0.6L)))
			{
				sum = X[i].first.second;
				X[i].second.first = ++count;
			}
		}
		vector<multiset<VLL_WRAPPER>> S(dsharp);
		for (ll i = 0LL; i < dsharp; i++)
		{
			S[i].insert(VLL_WRAPPER(0, 0));
			if (X[i].second.first == 0LL)
			{
				T = sumset(T, (S[i] = subsetsum(X[i].second.second)));
				for (auto j = T.begin(); j != T.end(); ++j) if (ll(*j) > X[i].first.first) j = --(T.erase(j));
				continue;
			}
			if (X[i].second.first != cur) cur = X[(k = i)].second.first;
			if (i == dsharp - 1LL || X[i].second.first != X[i + 1LL].second.first)
			{
				Pi = 0LL;
				for (ll j = k; j <= i; j++)
				{
					S[j] = subsetsum(X[j].second.second);
					Pi += X[j].first.second;
				}
				vector<vector<VLL_WRAPPER>> M(i - k + 1LL, vector<VLL_WRAPPER>(2LL * Pi + 1LL, VLL_WRAPPER(0LL, 0LL)));
				for (ll j = k; j <= i; j++)
				{
					for (ll l = 0LL; l <= 2LL * Pi; l++)
					{
						if (l == 0LL) M[j - k][l] = VLL_WRAPPER(X[j].first.first, 0);
						else if (S[j].find(VLL_WRAPPER(l, 0)) != S[j].end() && l <= X[j].first.first) M[j - k][l] = VLL_WRAPPER(X[j].first.first - l, 0) + (*(S[j].find(VLL_WRAPPER(l, 0)))).v;
						else M[j - k][l] = VLL_WRAPPER(N_INF, 0);
					}
				}
				size = i - k + 1LL;
				while (size > 1LL)
				{
					for (ll j = 1LL; j <= size / 2LL; j++)
					{
						M[2LL * j - 2LL][0LL] = M[2LL * j - 1LL][0LL];
						vector<VLL_WRAPPER> b = M[2LL * j - 2LL], a = M[2LL * j - 1LL];
						for (ll k = 0LL; k <= 2LL * Pi; k++) a[k] = a[k] + VLL_WRAPPER(k, 0);
						M[j - 1LL] = maxminskewedconvolution(a, b);
						M[j - 1LL].resize(2LL * Pi + 1);
						for (ll k = 0LL; k <= 2LL * Pi; k++) if ((M[j - 1LL][k] = M[j - 1LL][k] - k) < 0LL) M[j - 1LL][k] = VLL_WRAPPER(N_INF, 0);
					}
					for (ll j = size / 2LL + 1LL; 2LL * j - 1LL <= size; j++) M[j - 1LL] = M[2LL * j - 2LL];
					size = ll(ceil(size / 2.0L));
				}
				multiset<VLL_WRAPPER> Si, prefix = T, suffix = T;
				Si.insert(VLL_WRAPPER(0, 0));
				for (ll j = 0LL; j <= Pi; j++) if (M[0LL][j] >= 0LL) Si.insert(VLL_WRAPPER(j, 0) + M[0LL][j].v);
				prefix.erase(prefix.upper_bound(VLL_WRAPPER(max(X[k].first.first - Pi, 0LL), 0)), prefix.end());
				suffix.erase(suffix.begin(), suffix.lower_bound(VLL_WRAPPER(max(X[k].first.first - Pi, 0LL), 0)));
				if (X[k].first.first >= Pi) for (auto j : sumset(prefix, Si)) if (T.find(j) == T.end()) T.insert(j);
				vector<VLL_WRAPPER> Mprime(2LL * Pi + 1LL, VLL_WRAPPER(N_INF, 0));
				for (auto j : suffix) Mprime[ll(j) - max(X[k].first.first - Pi, 0LL)] = VLL_WRAPPER(0LL, 0LL);
				for (ll j = 0LL; j <= 2LL * Pi; j++) M[0LL][j] = M[0LL][j] - max(0LL, X[k].first.first - Pi);
				vector<VLL_WRAPPER> b = Mprime, a = M[0LL];
				for (ll j = 0LL; j <= 2LL * Pi; j++) a[j] = a[j]+ j;
				vector<VLL_WRAPPER> C = maxminskewedconvolution(a, b);
				C.resize(2LL * Pi + 1LL);
				for (ll j = 0LL; j <= 2LL * Pi; j++) if ((C[j] = C[j] - j) < 0LL) C[j] = VLL_WRAPPER(N_INF, 0);
				for (ll j = 0LL; j <= 2LL * Pi; j++) if (C[j] == 0LL) if (T.find(VLL_WRAPPER(j + max(X[k].first.first - Pi, 0LL), 0)) == T.end()) T.insert(VLL_WRAPPER(j + max(X[k].first.first - Pi, 0LL), 0) + C[j].v);
				for (auto j = T.begin(); j != T.end(); ++j) if (ll(*j) > X[i].first.first) j = --(T.erase(j));
			}
		}
		if (T.size() == 0LL)
		{
			return VLL_WRAPPER(-1LL, 0);
		}
		else
		{
			VLL_WRAPPER ans = *(--T.end());
			ans.sum = P - ans.sum;
			return ans;
		}
	}
	~MQTTScheduler()
	{
		J.clear();
	}
};