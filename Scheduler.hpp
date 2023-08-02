#pragma once
#include "MQTT.hpp"
class VLL_WRAPPER
{
public:
    ll sum = 0;
    vector<ll> v;
    VLL_WRAPPER()
    {
        sum = 0;
    }
    VLL_WRAPPER(ll x, ll id)
    {
        sum = x;
        if (x != 0)
        {
            v.push_back(id);
        }
    }
    VLL_WRAPPER(const VLL_WRAPPER& x)
    {
        sum = x.sum;
        v = x.v;
    }
    bool operator<(const VLL_WRAPPER& obj) const {
        return sum < obj.sum;
    }
    bool operator>(const VLL_WRAPPER& obj) const {
        return sum > obj.sum;
    }
    bool operator<=(const VLL_WRAPPER& obj) const {
        return sum <= obj.sum;
    }
    bool operator>=(const VLL_WRAPPER& obj) const {
        return sum >= obj.sum;
    }
    bool operator==(const VLL_WRAPPER& obj) const {
        return sum == obj.sum;
    }
    bool operator<(const ll compare) const {
        return sum < compare;
    }
    bool operator>(const ll compare) const {
        return sum > compare;
    }
    bool operator<=(const ll compare) const {
        return sum <= compare;
    }
    bool operator>=(const ll compare) const {
        return sum >= compare;
    }
    bool operator==(const ll compare) const {
        return sum == compare;
    }
    VLL_WRAPPER operator*(const ll compare) const {
        VLL_WRAPPER next(0, 0);
        next.sum = sum;
        next.sum *= compare;
        next.v.insert(next.v.end(), v.begin(), v.end());
        return next;
    }
    VLL_WRAPPER operator+(const ll compare) const {
        VLL_WRAPPER next(0, 0);
        next.sum = sum;
        next.sum += compare;
        next.v.insert(next.v.end(), v.begin(), v.end());
        return next;
    }
    VLL_WRAPPER operator-(const ll compare) const {
        VLL_WRAPPER next(0, 0);
        next.sum = sum;
        next.sum -= compare;
        next.v.insert(next.v.end(), v.begin(), v.end());
        return next;
    }
    VLL_WRAPPER operator/(const ll compare) const {
        VLL_WRAPPER next(0, 0);
        next.sum = sum;
        next.sum /= compare;
        next.v.insert(next.v.end(), v.begin(), v.end());
        return next;
    }
    explicit operator ll() const {
        return sum;
    }
    VLL_WRAPPER operator+(const VLL_WRAPPER& obj) const
    {
        VLL_WRAPPER next(0, 0);
        next.sum = sum + obj.sum;
        next.v.insert(next.v.end(), v.begin(), v.end());
        next.v.insert(next.v.end(), obj.v.begin(), obj.v.end());
        return next;
    }
    VLL_WRAPPER operator+(const vector<ll>& obj) const
    {
        VLL_WRAPPER next(0, 0);
        next.sum = sum;
        next.v.insert(next.v.end(), v.begin(), v.end());
        next.v.insert(next.v.end(), obj.begin(), obj.end());
        return next;
    }
    friend ostream& operator<<(ostream& os, const VLL_WRAPPER& val)
    {
        os << val.sum << " ( ";
        for (auto b : val.v)
        {
            if (b != 0)
                os << b << " ";
        }
        os << ")";
        return os;
    }
};
class Scheduler
{
public:
    virtual void fft(vector<cd>& a, bool invert)
    {
        ll n = a.size();
        for (ll i = 1LL, j = 0LL; i < n; i++)
        {
            ll bit = n >> 1LL;
            for (; j & bit; bit >>= 1LL) j ^= bit;
            j ^= bit;
            if (i < j) swap(a[i], a[j]);
        }
        for (ll len = 2LL; len <= n; len <<= 1LL)
        {
            double ang = 2LL * M_PI / len * (invert ? -1LL : 1LL);
            cd wlen(cos(ang), sin(ang));
            for (ll i = 0LL; i < n; i += len)
            {
                cd w(1LL);
                for (ll j = 0LL; j < len / 2LL; j++)
                {
                    cd u = a[i + j], v = a[i + j + len / 2LL] * w;
                    a[i + j] = u + v;
                    a[i + j + len / 2LL] = u - v;
                    w *= wlen;
                }
            }
        }
        if (invert) for (cd& x : a) x /= n;
    }
    /*
    virtual vector<ll> multiply(vector<ll>const& a, vector<ll> const& b)
    {
        vector<cd> fa(a.begin(), a.end()), fb(b.begin(), b.end());
        ll n = 1LL;
        while (n < ll(a.size() + b.size())) n <<= 1LL;
        fa.resize(n);
        fb.resize(n);
        fft(fa, false);
        fft(fb, false);
        for (ll i = 0LL; i < n; i++) fa[i] *= fb[i];
        fft(fa, true);
        vector<ll> result(n);
        for (ll i = 0LL; i < n; i++) result[i] = (ll)round(fa[i].real());
        return result;
    }
    */
    virtual vector<VLL_WRAPPER> multiply(vector<VLL_WRAPPER>const& a, vector<VLL_WRAPPER> const& b)
    {
        ll n = 1LL;
        while (n < ll(a.size() + b.size())) n <<= 1LL;
        vector<VLL_WRAPPER> result(n);
        for (ll i = 0; i < a.size(); i++)
        {
            for (ll j = 0; j < b.size(); j++)
            {
                result[i + j] = a[i] + b[j];
            }
        }
        return result;
    }
    virtual multiset<VLL_WRAPPER> sumset(multiset<VLL_WRAPPER> const& a, multiset<VLL_WRAPPER> const& b)
    {
        ll n = 1LL;
        while (n < ll(a.size() + b.size())) n <<= 1LL;
        multiset<VLL_WRAPPER> result;
        for (multiset<VLL_WRAPPER>::const_iterator i = a.begin(); i != a.end(); ++i)
        {
            for (multiset<VLL_WRAPPER>::const_iterator j = b.begin(); j != b.end(); ++j)
            {
                result.insert((*i) + (*j));
            }
        }
        return result;
    }
    virtual multiset<VLL_WRAPPER> subsetsum(multiset<VLL_WRAPPER>& X)
    {
        if (X.size() == 1LL)
        {
            X.insert(VLL_WRAPPER(0LL, 0LL));
            return X;
        }
        multiset<VLL_WRAPPER> a, b;
        multiset<VLL_WRAPPER>::iterator i = X.begin();
        multiset<VLL_WRAPPER>::iterator it = X.begin();
        for (ll j = 1LL; j < (ll)X.size() / 2LL + 1LL; j++) ++i;
        for (; it != i; ++it) a.insert(*(it));
        for (; it != X.end(); ++it) b.insert(*(it));
        return sumset(subsetsum(a), subsetsum(b));
    }
    virtual bool ckgreaterthanv(vector<VLL_WRAPPER> const& a, vector<VLL_WRAPPER> const& b, vector<VLL_WRAPPER> const& d, vector<pair<VLL_WRAPPER, ll>> const& a1, vector<pair<VLL_WRAPPER, ll>> const& b1, vector<vector<vector<VLL_WRAPPER>>> const& matrix, ll k, ll v)
    {
        assert(a.size() == b.size());
        ll n = a.size();
        ll i1, j1, i2, j2, i3, j3;
        for (i1 = 0LL; i1 < n; i1++) if (a1[i1].first > v) break;
        for (j1 = 0LL; j1 < n; j1++) if (b1[j1].first > v - ll(d[k])) break;
        if (i1 == n || j1 == n) return false;
        ll p = (ll)floor(pow(n, 2.0L / 3.0L));
        i2 = min((ll)(ceil((i1 + 1.0L) / p) * p), n) - 1LL;
        j2 = min((ll)(ceil((j1 + 1.0L) / p) * p), n) - 1LL;
        i3 = i1 / p;
        j3 = j1 / p;
        if (matrix[i3][j3][k] > 0LL) return true;
        for (ll i = i1; i < i2; i++) if (k - a1[i].second >= 0LL && k - a1[i].second < n && (b[k - a1[i].second]) / (4LL * n) >= (b1[j1].first) / (4LL * n)) return true;
        for (ll j = j1; j < j2; j++) if (k - b1[j].second >= 0LL && k - b1[j].second < n && (a[k - b1[j].second]) / (4LL * n) >= (a1[i1].first) / (4LL * n)) return true;
        return false;
    }
    virtual vector<VLL_WRAPPER> maxminskewedconvolution(vector<VLL_WRAPPER>& a, vector<VLL_WRAPPER>& b)
    {
        assert(a.size() == b.size());
        ll n = a.size();
        vector<pair<VLL_WRAPPER, ll>> a1(n), b1(n);
        vector<VLL_WRAPPER> a2(n), b2(n);
        vector<VLL_WRAPPER> d(2LL * n);
        for (ll i = 0LL; i < n; i++) a[i] = a[i] * 4LL * n + i;
        for (ll i = n; i < 2LL * n; i++) b[i - n] = b[i - n] * 4LL * n + i;
        for (ll i = 0LL; i < 2LL * n; i++) d[i] = VLL_WRAPPER(4LL * n * i, 0);
        for (ll i = 0LL; i < n; i++) a1[i] = { a[i], i };
        for (ll i = 0LL; i < n; i++) b1[i] = { b[i], i };
        ll p = (ll)floor(pow(n, 2.0L / 3.0L));
        sort(a1.begin(), a1.end());
        sort(b1.begin(), b1.end());
        vector<vector<vector<VLL_WRAPPER>>> matrix;
        ll i = 0LL;
        for (ll u = 0LL; u < n; u += p)
        {
            if (u == 0LL) continue;
            matrix.push_back({});
            for (ll w = 0; w < n; w += p)
            {
                if (w == 0LL) continue;
                for (ll i = 0LL; i < n; i++) a2[i] = (a[i] / (4LL * n) >= a1[u - 1LL].first / (4LL * n) ? (a[i] / (4LL * n)) : VLL_WRAPPER(0, 0));
                for (ll i = 0LL; i < n; i++) b2[i] = (b[i] / (4LL * n) >= b1[w - 1LL].first / (4LL * n) ? (b[i] / (4LL * n)) : VLL_WRAPPER(0, 0));
                matrix[i].push_back(multiply(a2, b2));
            }
            i++;
        }
        i = 0LL;
        for (ll u = 0LL; u < n; u += p)
        {
            ll w = n;
            if (u == 0LL) continue;
            for (ll i = 0LL; i < n; i++) a2[i] = (a[i] / (4LL * n) >= a1[u - 1LL].first / (4LL * n) ? (a[i] / (4LL * n)) : VLL_WRAPPER(0, 0));
            for (ll i = 0LL; i < n; i++) b2[i] = (b[i] / (4LL * n) >= b1[w - 1LL].first / (4LL * n) ? (b[i] / (4LL * n)) : VLL_WRAPPER(0, 0));
            matrix[i].push_back(multiply(a2, b2));
            i++;
        }
        matrix.push_back({});
        for (ll w = 0; w < n; w += p)
        {
            ll u = n;
            if (w == 0LL) continue;
            for (ll i = 0LL; i < n; i++) a2[i] = (a[i] / (4LL * n) >= a1[u - 1LL].first / (4LL * n) ? (a[i] / (4LL * n)) : VLL_WRAPPER(0, 0));
            for (ll i = 0LL; i < n; i++) b2[i] = (b[i] / (4LL * n) >= b1[w - 1LL].first / (4LL * n) ? (b[i] / (4LL * n)) : VLL_WRAPPER(0, 0));
            matrix[i].push_back(multiply(a2, b2));
        }
        ll u = n, w = n;
        for (ll i = 0LL; i < n; i++) a2[i] = (a[i] / (4LL * n) >= a1[u - 1LL].first / (4LL * n) ? (a[i] / (4LL * n)) : VLL_WRAPPER(0, 0));
        for (ll i = 0LL; i < n; i++) b2[i] = (b[i] / (4LL * n) >= b1[w - 1LL].first / (4LL * n) ? (b[i] / (4LL * n)) : VLL_WRAPPER(0, 0));
        matrix[i].push_back(multiply(a2, b2));
        vector<VLL_WRAPPER> c(2LL * n);
        for (ll k = 0LL; k < 2LL * n; k++)
        {
            vector<VLL_WRAPPER> options;
            ll l = 0LL, r = a1.size() - 1LL;
            while (l < r)
            {
                ll m = l + (r - l) / 2LL;
                if (ckgreaterthanv(a, b, d, a1, b1, matrix, k, ll(a1[m].first))) l = m + 1LL;
                else r = m;
            }
            VLL_WRAPPER candidateone = a1[l].first;
            l = 0LL;
            r = b1.size() - 1LL;
            while (l < r)
            {
                ll m = l + (r - l) / 2LL;
                if (ckgreaterthanv(a, b, d, a1, b1, matrix, k, ll(b1[m].first + d[k]))) l = m + 1LL;
                else r = m;
            }
            VLL_WRAPPER candidatetwo = b1[r].first + d[k];
            c[k] = min(candidateone, candidatetwo) / (4LL * n);
        }
        return c;
    }
};