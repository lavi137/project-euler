// euler.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include <stdio.h>
#include <iostream>
#include <algorithm>
#include <string>
#include <tuple>
#include <stack>
#include <queue>
#include <cmath>
#include <set>
#include <iomanip>
#include <functional> 
#include <chrono>
#include <fstream>
#include <bit>
#include <bitset>
#include <cstdint>
#include <numbers>
#include <limits>
#include <string_view>
#define Find(v,target) find(v,0,v.size()-1,target)
using namespace std;
typedef long long ll;
typedef vector<long long> vl;
typedef vector<vl> vvl;
typedef pair<ll, ll> pl;
typedef long double ld;
typedef ll(*func)(ll, ll, ll);//func is every function which gets (ll,ll,ll) and returns ll
typedef priority_queue<ll>pq;
ll string_to_int(string s)
{
    ll flag = 0;
    if (s[0] == '-')flag = 1;
    ll sum = 0;
    ll mult = 1;
    for (int i = s.size() - 1; i >= flag; i--)
    {
        sum += mult * (s[i] - '0');
        mult *= 10;
    }
    if (flag)return -sum;
    return sum;
}
inline ll power(ll a, ll b, ll limit = 1e18, func f = [](ll l1, ll l2, ll limit) {return (l1 * l2) % limit; }, ll d = 1)//can do other logrithmic operations
{
    if (a == 0 && b == 0)return d;
    ll ans = d, tmp;
    while (b != 0)
    {
        tmp = (b & 1) ? a : d;
        ans = f(ans, tmp, limit);
        a = f(a, a, limit);
        b /= 2;
    }
    return ans;
}
ll find(vector<ll>& v, ll s, ll e, ll target)//NEED TO CHECK IT WORKS AND THE USES OF IT IN OTHER PROGRAMS
{
    ll mid;
    while (s <= e)
    {
        mid = (s + e) / 2;
        if (v[mid] > target)e = mid - 1;
        else if (v[mid] < target)s = mid + 1;
        else if (v[mid] == target)return mid;
    }
    return -1;
}
ll log(ll a, ll base = (10))
{
    ll cnt = 0;
    while (a /= base)cnt++;
    return cnt + 1;
}
void sieve(vector<ll>& prime, vector<ll>& div, vector<bool>& p, ll n)//find primes
{
    for (ll i = 2; i < n; i++)
    {
        if (div[i] == 0)
        {
            div[i] = i; prime.push_back(i);
            p[i] = 1;
        }
        for (ll j = 0; j < prime.size() && prime[j] * i < n; j++)
        {
            if (prime[j] > div[i])break;
            div[prime[j] * i] = prime[j];
        }
    }
}
void nsieve(vl& p, ll n) {//faster than sieve but doesnt have div
    const ll S = sqrt(n);//can be changed but i think its quite optimal. should check it sometimes
    vector<ll> primes;
    ll nsqrt = sqrt(n);
    vector<char> is_prime(nsqrt + 2, 1);
    for (ll i = 2; i <= nsqrt; i++) {
        if (is_prime[i]) {
            primes.push_back(i);
            for (ll j = i * i; j <= nsqrt; j += i)
                is_prime[j] = 0;
        }
    }
    vector<char> block(S);
    for (ll k = 0; k * S <= n; k++) {//can probably start with k=1 because first loop is up.
        fill(block.begin(), block.end(), 1);
        ll start = k * S;
        for (ll p : primes) {
            ll start_idx = (start + p - 1) / p;
            ll j = max(start_idx, p) * p - start;
            for (; j < S; j += p)
                block[j] = 0;
        }
        if (k == 0)
            block[0] = block[1] = 0;
        for (ll i = 0; i < S && start + i <= n; i++) {
            if (block[i])
            {
                p.push_back(start + i);
            }
        }
    }
}
ll gcd(ll a, ll b)
{
    a = abs(a); b = abs(b);
    while (a)
    {
        b %= a;
        swap(a, b);
    }
    return b;
}
ll egcd(ll a, ll b, ll& x, ll& y) {//extended gcd
    x = 1, y = 0;
    ll x1 = 0, y1 = 1, a1 = a, b1 = b;
    while (b1) {
        ll q = a1 / b1;
        tie(x, x1) = make_tuple(x1, x - q * x1);
        tie(y, y1) = make_tuple(y1, y - q * y1);
        tie(a1, b1) = make_tuple(b1, a1 - q * b1);
    }
    return a1;
}
ll spin(ll a)//spins number digits 123->312
{
    ll n = 0, tmp = a, b;
    while (tmp /= 10)
    {
        n++;
    }
    tmp = a;
    b = tmp % 10;
    tmp /= 10;
    return b * power(10, n) + tmp;
}
bool palindrome(ll a, ll base = 10)
{
    ll tmp = a, cnt = 0;
    while (tmp > 0)
    {
        tmp /= base; cnt++;
    }
    for (ll i = 1; i <= cnt / 2; i++)
    {
        if (((a % power(base, i)) / power(base, i - 1)) != ((a % power(base, cnt - i + 1)) / power(base, cnt - i)))return 0;
    }
    return 1;
}
bool is_prime(ll p, ll iteration = 40)
{//last 3 in a for divisibillity issues in last 3 numbers in base (they arent prime)
    vl prime = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,193,407521,299210837 };
    vl base = { 2, 325, 9375, 28178, 450775, 9780504, 1795265022 };
    for (ll i : prime)
    {
        if (p == i)return true;
        if (p % i == 0)return false;
    }
    if (p < 2)
    {
        return false;
    }
    if (p != 2 && p % 2 == 0)
    {
        return false;
    }
    ll s = p - 1;
    while (s % 2 == 0)s /= 2;
    if (p < 3215031751)//you need here only 4 bases
    {
        for (ll i = 0; i < 4; i++)//these 7 bases cover everything until 2^64 so u find if prime in 7 iterations
        {
            ll a = prime[i], temp = s;
            ll mod = power(a, temp, p, [](ll l1, ll l2, ll m1) {return power(l1, l2, m1, [](ll x, ll y, ll z) {return (x + y) % z; }, 0); }, 1);
            while (temp != p - 1 && mod != 1 && mod != p - 1)
            {
                mod = power(mod, mod, p, [](ll x, ll y, ll z) {return (x + y) % z; }, 0);
                temp *= 2;
            }
            if (mod != p - 1 && temp % 2 == 0)
            {
                return false;
            }
        }
    }
    else if (p < 1e18)
    {
        for (ll i = 0; i < 7; i++)//these 7 bases cover everything until 2^64 so u find if prime in 7 iterations
        {
            ll a = base[i], temp = s;
            ll mod = power(a, temp, p, [](ll l1, ll l2, ll m1) {return power(l1, l2, m1, [](ll x, ll y, ll z) {return (x + y) % z; }, 0); }, 1);
            while (temp != p - 1 && mod != 1 && mod != p - 1)
            {
                mod = power(mod, mod, p, [](ll x, ll y, ll z) {return (x + y) % z; }, 0);
                temp *= 2;
            }
            if (mod != p - 1 && temp % 2 == 0)
            {
                return false;
            }
        }
    }
    else
    {
        for (ll i = 0; i < iteration; i++)
        {
            ll a = rand() % (p - 1) + 1, temp = s;
            ll mod = power(a, temp, p, [](ll l1, ll l2, ll m1) {return power(l1, l2, m1, [](ll x, ll y, ll z) {return (x + y) % z; }, 0); }, 1);
            while (temp != p - 1 && mod != 1 && mod != p - 1)
            {
                mod = power(mod, mod, p, [](ll x, ll y, ll z) {return (x + y) % z; }, 0);
                temp *= 2;
            }
            if (mod != p - 1 && temp % 2 == 0)
            {
                return false;
            }
        }
    }
    return true;
}
namespace divisors
{
    ll f(ll x, ll c, ll mod) {
        return (power(x, x, mod, [](ll l1, ll l2, ll m1) {return (l1 + l2) % m1; }, 0) + c) % mod;
    }
    ll brent(ll n, ll x0 = 2, ll c = 1) {
        ll x = x0;
        ll g = 1;
        ll q = 1;
        ll xs, y;
        ll m = 128;
        ll l = 1;
        while (g == 1) {
            y = x;
            for (ll i = 1; i < l; i++)
                x = f(x, c, n);
            ll k = 0;
            while (k < l && g == 1) {
                xs = x;
                for (ll i = 0; i < m && i < l - k; i++) {
                    x = f(x, c, n);
                    q = power(q, abs(y - x), n, [](ll l1, ll l2, ll m1) {return (l1 + l2) % m1; }, 0);
                }
                g = gcd(q, n);
                k += m;
            }
            l *= 2;
        }
        if (g == n) {
            do {
                xs = f(xs, c, n);
                g = gcd(abs(xs - y), n);
            } while (g == 1);
        }
        return g;
    }
    ll find_div(ll n)//i'll assume when calling this function that the number isnt prime
    {
        //if (is_prime(n))return n;
        vl prime = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97 };
        ll tmp;
        for (int x : prime)if (n % x == 0)return x;
        for (ll x0 = 2;; x0++)
        {
            tmp = brent(n, x0);
            if (tmp < n)return tmp;
        }
        return n;
    }
    vl find_primes(ll n)
    {
        if (n == 1)return {};
        if (is_prime(n))
        {
            return { n };
        }
        ll tmp = find_div(n);
        tmp = min(tmp, n / tmp);
        vl v1 = find_primes(tmp);
        tmp = n / tmp;
        for (ll p : v1)while (tmp % p == 0)tmp /= p;
        vl v2 = find_primes(tmp);
        for (ll x : v2)v1.push_back(x);
        return v1;
    }
    vl find_divisors(ll n)
    {
        vl p = find_primes(n);
        ll s = p.size(), tmp = n;
        vl d(s);
        vl div = { 1 };
        ll sz;
        for (ll i = 0; i < s; i++)
        {
            sz = div.size();
            while (tmp % p[i] == 0)
            {
                d[i]++; tmp /= p[i];
            }
            ll m = 1;
            for (ll j = 0; j < d[i]; j++)
            {
                m *= p[i];
                for (ll k = 0; k < sz; k++)
                {
                    div.push_back(m * div[k]);
                }
            }
        }
        return div;
    }
}
ll phi(ll n) {
    if (is_prime(n))return n - 1;
    ll ans = 1;
    vl p = divisors::find_primes(n);
    for (ll d : p)
    {
        ans *= (d - 1); n /= d;
        while (n % d == 0)
        {
            ans *= d; n /= d;
        }
    }
    return ans;
}
ll crt(ll p, ll q, ll a, ll b, bool p1 = 0, bool q1 = 0, ll p2 = 0, ll q2 = 0)//finds solution for number that is a mod p and b mod q
{// p1=1 is p prime,p1=2 is p power of p2 rn i only work with p1=0,1
    func f = [](ll x, ll y, ll z) {return (x + y) % z; };
    if (p1 == 1)
    {
        swap(p, q); swap(p1, q1);
    }
    ll n = power(p, phi(q) - 1, q, [](ll l1, ll l2, ll m1) {return power(l1, l2, m1, [](ll x, ll y, ll z) {return (x + y) % z; }, 0); }, 1);
    ll m = a + power(power(b - a + p * q, p, p * q, f, 0), n, p * q, f, 0);
    return m;
}
bool miillerTest(ll d, ll n)
{
    ll a = 2 + rand() % (n - 4);
    ll x = power(a, d, n, [](ll l1, ll l2, ll m1) {return power(l1, l2, m1, [](ll x, ll y, ll z) {return (x + y) % z; }, 0); }, 1);
    if (x == 1 || x == n - 1)
        return true;
    while (d != n - 1)
    {
        x = power(x, 2, n, [](ll l1, ll l2, ll m1) {return power(l1, l2, m1, [](ll x, ll y, ll z) {return (x + y) % z; }, 0); }, 1);
        d *= 2;
        if (x == 1) return false;
        if (x == n - 1)    return true;
    }
    return false;
}
bool isPrime(ll n, ll k)//outdated kept only for functions that used it before
{
    if (n <= 1 || n == 4)  return false;
    if (n <= 3) return true;
    ll d = n - 1;
    while (d % 2 == 0)
        d /= 2;
    for (ll i = 0; i < k; i++)
        if (!miillerTest(d, n))
            return false;

    return true;
}
ll stonelli(ll n, ll p)//finds r such that r^2=n mod p
{
    n %= p;
    ll nr, q, s = 0, m, c, t, r, b, tmp;
    if (p == 2)return n;
    if (p % 4 == 3)
    {
        ll r = power(n, (p + 1) / 4, p);
        if (r * r % p == n)return r;
        else return -1;
    }
    for (ll i = 1;; i++)if (power(i, (p - 1) / 2, p) != 1)
    {
        nr = i; break;
    }
    q = p - 1;
    while (q % 2 == 0)
    {
        s++; q >>= 1;
    }
    m = s, c = power(nr, q, p), t = power(n, q, p), r = power(n, (q + 1) / 2, p);
    while (1)
    {
        if (!t)return 0;
        if (t == 1)return r;
        tmp = t;
        for (ll i = 1; i <= m; i++)
        {
            if (i == m)return -1;
            tmp *= tmp; tmp %= p;
            if (tmp == 1)
            {
                b = power(c, power(2, m - i - 1, p - 1), p);
                m = i, c = b * b % p, t = t * b % p * b % p, r = r * b % p;
            }
        }
    }
    return -1;
}
namespace dsu
{
    void make_set(vl &parent,vl &size,ll v) {
        parent[v] = v;
        size[v] = 1;
    }
    ll find_set(ll v,vl &parent) {
        if (v == parent[v])
            return v;
        return parent[v] = find_set(parent[v], parent);
    }
    void union_sets(ll a, ll b,vl &parent,vl &size) {
        a = find_set(a,parent);
        b = find_set(b,parent);
        if (a != b) {
            if (size[a] < size[b])
                swap(a, b);
            parent[b] = a;
            size[a] += size[b];
        }
    }
}
/*struct rational
{
    pl r;
    rational(ll p, ll q)
    {
        r.first = p; r.second = q;
    }
    rational gcd()
    rational operator+(rational& o)
    {
        rational c = rational(r.first * o.second + r.second * o.first, r.second * o.second);
        return gcd(c.first,)
    }
    pl gcd(pl c)
    {
        ll a = abs(c.first), b = abs(c.second);
        if (a < b)swap(a, b);
        while (b > 0)
        {
            a %= b;
            if (a < b)swap(a, b);
        }
        return { c.first / a,c.second / a };
    }
    pl mult(pl a, pl b)
    {
        pl c = { a.first * b.first,a.second * b.second };
        if (c.second < 0)c = { -a.first * b.first,-a.second * b.second };
        return gcd(c);
    }
    pl div(pl a, pl b)
    {
        return mult(a, { b.second,b.first });
    }
    pl add(pl a, pl b)
    {
        pl c = { a.first * b.second + a.second * b.first,a.second * b.second };
        return gcd(c);
    }
    pl sub(pl a, pl b)
    {
        return add(a, { -b.first,b.second });
    }
    bool bigger(pl a, pl b)
    {
        if (a.first * b.second > a.second * b.first)return 1;
        return 0;
    }
}*/
namespace Q1//175 microseconds
{
    ll solve()
    {
        ll n = 1e3, a = 3, b = 5, sum = 0;
        for (ll i = 1; i < n; i++)
        {
            if (i % 3 == 0 || i % 5 == 0)sum += i;
        }
        return sum;
    }
}
namespace Q2//213 microseconds
{
    ll solve()
    {
        ll n = 4e6, sum = 0; vl f = { 1,1 };
        for (ll i = 2;; i++)
        {
            f.push_back(f[i - 1] + f[i - 2]);
            if (f[i] > n)break;
            if (f[i] % 2 == 0)sum += f[i];
        }
        return sum;
    }
}
namespace Q3//163 microseconds
{
    ll solve()
    {
        ll n = 600851475143;
        for (ll i = 2; n > 1; i++)
        {
            while (n % i == 0 && n != i)
            {
                n /= i;
            }
            if (n == i)return i;
        }
    }
}
namespace Q4//20,522 microseconds
{
    ll solve()
    {
        ll m = -1;
        for (ll i = 100; i < 1e3; i++)for (ll j = i; j < 1e3; j++)
        {
            if (palindrome(i * j))m = max(m, i * j);
        }
        return m;
    }
}
namespace Q5//231 microseconds
{
    ll solve()//lcm
    {
        return 16 * 9 * 5 * 7 * 11 * 13 * 17 * 19;
    }
}
namespace Q6//197 microseconds
{
    ll solve()//trivial
    {
        ll n = 100;
        return (n * (n + 1) / 2) * (n * (n + 1) / 2) - n * (n + 1) * (2 * n + 1) / 6;
    }
}
namespace Q7//7,092 microseconds
{
    ll solve()
    {
        ll n = 1e6;
        vl div(n), prime;
        vector<bool>q(n);
        sieve(prime, div, q, n);
        return prime[10000];
    }
}
namespace Q8//ill solve later idk how to take that input
{
    ll solve()
    {
        ll n = 1000, m = -1, cur;
        string s; cin >> s;
        for (ll i = 0; i < n - 13; i++)
        {
            cur = 1;
            for (ll j = 0; j < 13; j++)
            {
                cur *= s[i + j] - '0';
            }
            m = max(m, cur);
        }
        return m;
    }
}
namespace Q9//225 microseconds
{
    ll solve()//m(m+n)=500,m>n so m=20,n=5
    {
        ll m = 20, n = 5;
        return (m * m - n * n) * (m * m + n * n) * (2 * m * n);
    }
}
namespace Q10//13,981 microseconds
{
    ll solve()
    {
        ll n = 2e6, sum = 0;
        vl div(n), prime;
        vector<bool>q(n);
        sieve(prime, div, q, n);
        for (ll x : prime)sum += x;
        return sum;
    }
}
ll Q11(vvl& v)
{
    ll k = 1;
    ll n = v.size(), m = v[0].size();
    for (ll i = 0; i < n; i++)
    {
        for (ll j = 0; j < m; j++)
        {
            if (j + 3 < m)
            {
                k = max(k, v[i][j] * v[i][j + 1] * v[i][j + 2] * v[i][j + 3]);
            }
            if (i + 3 < n)
            {
                k = max(k, v[i][j] * v[i + 1][j] * v[i + 2][j] * v[i + 3][j]);
            }
            if (j + 3 < m && i + 3 < n)
            {
                k = max(k, v[i][j] * v[i + 1][j + 1] * v[i + 2][j + 2] * v[i + 3][j + 3]);
            }
            if (j + 3 < m && i - 3 >= 0)
            {
                k = max(k, v[i][j] * v[i - 1][j + 1] * v[i - 2][j + 2] * v[i - 3][j + 3]);
            }
        }
    }
    return k;
}
namespace Q12//942,545 microseconds
{
    ll solve()
    {
        ll n = 1e8, d, cur, tmp, m;
        vl div(n), prime;
        vector<bool>q(n);
        sieve(prime, div, q, n);
        for (ll i = 1;; i++)
        {
            cur = i * (i + 1) / 2, m = 1;
            if (cur >= n)return -1;
            tmp = 1;
            while (cur > 1)
            {
                tmp = div[cur];
                d = 1;
                while (div[cur] == tmp)
                {
                    d++; cur /= div[cur];
                }
                m *= d;
            }
            if (m > 500)return i * (i + 1) / 2;
        }
        return -1;
    }
}
ll Q13(vl& v)
{
    ll sum = 0;
    for (int i = 0; i < 100; i++)
    {
        sum += v[i];
    }
    ll x = 1e10;
    while (sum >= x)sum /= 10;
    return sum;
}
ll Q14(ll n)
{
    ll m = 0;
    ll p = -1;
    for (int i = 2; i < n; i++)
    {
        ll cnt = 0, tmp = i;
        while (tmp != 1)
        {
            if (tmp & 1)tmp = 3 * tmp + 1;
            else tmp /= 2;
            cnt++;
        }
        if (m < cnt)
        {
            m = cnt;
            p = i;
        }
    }
    return p;
}
ll Q15()
{
    return 137846528820;//it is 40 choose 20
}
namespace Q16
{
    //have in python
}
ll Q17()
{
    return 21124;//annoying to do
}
ll Q18(vvl& v)
{
    vvl ans = v;
    for (int i = 1; i < v.size(); i++)
    {
        for (int j = 0; j < v[i].size(); j++)
        {
            if (j == 0)ans[i][j] += ans[i - 1][j];
            else if (j == v[i].size() - 1)ans[i][j] += ans[i - 1][j - 1];
            else
            {
                ans[i][j] += max(ans[i - 1][j], ans[i - 1][j - 1]);
            }
        }
    }
    ll m = -1;
    for (int i = 0; i < v[v.size() - 1].size(); i++)
    {
        m = max(m, ans[v.size() - 1][i]);
    }
    return m;
}
ll Q19()
{
    return 171;//annoying
}
namespace Q20
{
    //have in python
}
ll Q21(ll n)
{
    vl ans(n + 1);
    for (int i = 1; i < n; i++)
    {
        for (int j = 1; j < i; j++)
        {
            if (i % j == 0)ans[i] += j;
        }
    }
    ll sum = 0;
    for (int i = 1; i < 10000; i++)
    {
        if (i == ans[i])continue;
        if (i == ans[ans[i]])sum += i;
    }
    return sum;
}
ll Q22(vector<string>& v)
{
    sort(v.begin(), v.end());
    ll sum = 0;
    for (int i = 0; i < v.size(); i++)
    {
        ll tmp = 0;
        for (int j = 0; j < v[i].size(); j++)
        {
            tmp += (v[i][j] - 'A' + 1);
        }
        sum += tmp * (i + 1);
    }
    return sum;
}
ll Q23()
{
    vl v;
    v.push_back(12);
    for (int i = 13; i <= 28123; i++)
    {
        ll sum = 0;
        for (int j = 1; j < i; j++)
        {
            if (i % j == 0)sum += j;
        }
        if (sum > i)v.push_back(i);
    }
    vector<bool>seen(100000);
    for (int i = 0; i < v.size(); i++)
    {
        for (int j = i; j < v.size(); j++)
        {
            seen[v[i] + v[j]] = true;
        }
    }
    ll sum = 0;
    for (int i = 1; i <= 28123; i++)
    {
        if (!seen[i])sum += i;
    }
    return sum;
}
ll Q24(ll n)
{
    vector<ll>ans;
    ll num = 1;
    ll tmp = 1;
    ll MAX = 1e6;
    while (ans.size() < n)
    {
        tmp = 1;
        for (int i = 2; i <= n - ans.size(); i++)
        {
            tmp *= (i - 1);
        }
        for (int j = 1; j < 100; j++)
        {
            num += tmp;
            if (num > MAX)
            {
                num -= tmp;
                ans.push_back(j);
                break;
            }
        }
    }
    return 1;
}
namespace Q25//151 microseconds
{
    ll solve()//closed form
    {
        double b = (1 + sqrt(5)) / 2;
        return (999 * log(10) + log(sqrt(5))) / log(b) + 1;
    }
}
ll Q26(ll n)
{
    ll d = 1;
    ll len = 0;
    ll MAX = 1e18;
    for (ll i = 2; i < n; i++)
    {
        if (MAX % i == 0)continue;
        else
        {
            ll tmp = 1;
            ll num = 10;
            ll j = i;
            while (j % 2 == 0)j /= 2;
            while (j % 5 == 0)j /= 5;
            while (num % j != 1)
            {
                tmp++; num *= 10;
                num %= j;
            }
            if (len < tmp)
            {
                d = i; len = tmp;
            }
        }
    }
    return d;
}
ll Q27(ll lim)
{
    ll sum = 0;
    const ll N = 5 * lim * lim;
    vector<ll> lp(N + 1);
    vector<ll> pr;
    for (ll i = 2; i <= N; ++i) {
        if (lp[i] == 0) {
            lp[i] = i;
            pr.push_back(i);
        }
        for (ll j = 0; j < (int)pr.size() && pr[j] <= lp[i] && i * pr[j] <= N; ++j) {
            lp[i * pr[j]] = pr[j];
        }
    }
    ll maxlen = 0;
    ll mult = 0;
    for (ll i = -lim + 1; i < lim; i++)
    {
        for (ll k = 0; pr[k] <= lim; k++)
        {
            ll j = pr[k];
            ll cnt = 0;
            while (find(pr, 0, pr.size() - 1, cnt * cnt + i * cnt + j) != -1)
            {
                cnt++;
            }
            if (cnt > maxlen)
            {
                maxlen = cnt;
                mult = i * j;
                ll c = 2;
            }
        }
    }
    return mult;
}
ll Q28(ll n)
{
    ll sum = 1;
    for (int i = 3; i <= n; i += 2)
    {
        sum += (4 * i * i - 6 * (i - 1));
    }
    return sum;
}
ll Q29(ll a, ll b)
{
    set<ll>s;
    for (ll i = 2; i <= a; i++)
    {
        for (ll j = 2; j <= b; j++)
        {
            s.insert(log(power(i, j, 1e9 + 7)));
        }
    }
    return s.size();
    //the answer is 9183 (mod problems cant give it)
}
ll Q30(ll n)
{
    ll sum = 0;
    for (ll i = 10; i < 400000; i++)
    {
        ll tmp = i, check = 0;
        while (tmp != 0)
        {
            check += power(tmp % 10, 5);
            tmp /= 10;
        }
        if (i == check)
        {
            sum += check;
        }
    }
    return sum;
}
ll Q31(vl& coin, ll x)
{
    ll n = coin.size();
    vector<ll>window(x + 1);
    sort(coin.begin(), coin.end());
    for (int i = 0; i < n; i++)
    {
        window[coin[i]]++;
        for (int j = coin[i] + 1; j <= x; j++)
        {
            window[j] += window[j - coin[i]];
        }
    }
    return window[x];
}
ll Q32()
{
    ll sum = 0;
    vector<ll>ans;
    for (ll i = 1; i <= 10000; i++)
    {
        bool flag = false;
        if (flag)continue;
        for (ll j = 1; j < i; j++)
        {
            if (flag)
            {
                j = i;
                continue;
            }
            if (i % j != 0)continue;
            ll k1 = i / j;
            ll k2 = j;
            ll k3 = i;
            vector<bool>seen(10, false);
            ll len = 0;
            while (k1 != 0)
            {
                len++;
                seen[k1 % 10] = true;
                k1 /= 10;
            }
            while (k3 != 0)
            {
                len++;
                seen[k3 % 10] = true;
                k3 /= 10;
            }
            while (k2 != 0)
            {
                len++;
                seen[k2 % 10] = true;
                k2 /= 10;
            }
            bool check = !seen[0];
            for (ll m = 1; m < 10; m++)
            {
                check = (check && seen[m]);
            }
            if (check && len == 9)
            {
                sum += i;
                ans.push_back(i);
                flag = true;
            }
        }
    }
    return sum;
}
ll Q33()
{
    ll nume = 1, deno = 1;
    for (ll i = 10; i < 100; i++)
    {
        for (ll j = i + 1; j < 100; j++)
        {
            if (j % 10 == 0 && i % 10 == 0)continue;
            if (i * (j / 10) == j * (i / 10) && (j % 10 == i % 10))
            {
                nume *= j; deno *= i;
            }
            else if (i * (j % 10) == j * (i / 10) && (j / 10 == i % 10))
            {
                nume *= j; deno *= i;
            }
            else if (i * (j % 10) == j * (i % 10) && (j / 10 == i / 10))
            {
                nume *= j; deno *= i;
            }
            else if (i * (j / 10) == j * (i % 10) && (j % 10 == i / 10))
            {
                nume *= j; deno *= i;
            }
        }
    }
    return 1;
}
ll Q34()
{
    int fact[10] = { 1,1,2,6,24,120,720,5040,40320,362880 };
    int sum = 0, cur, tmp;
    for (int i = 10; i < 7 * 362880; i++)
    {
        cur = 0; tmp = i;
        while (tmp > 0)
        {
            cur += fact[tmp % 10];
            tmp /= 10;
        }
        if (cur == i)sum += i;
    }
    return sum;
}
ll Q35()
{
    ll n = 1e6, tmp, sum = 0;
    vector <ll> prime;
    vector<ll>div(n);
    vector<bool>p(n);
    sieve(prime, div, p, n);
    for (int i = 2; i < n; i++)
    {
        int flag = 1; tmp = i;
        for (int j = 0; j < 6; j++)
        {
            if (!p[tmp])flag = 0;
            tmp = spin(tmp);
        }
        if (flag)
        {
            sum += 1;
        }
    }
    return sum;
}
ll Q36()
{
    int n = 1e6, sum = 0;

    for (int i = 1; i < n; i++)
    {
        if (palindrome(i, 10) && palindrome(i, 2))sum += i;
    }
    return sum;
}
ll Q37()
{
    ll n = 1e6, tmp, sum = -17, flag;//-17 for -2-3-5-7
    vector <ll> prime;
    vector<ll>div(n);
    vector<bool>p(n);
    sieve(prime, div, p, n);
    for (ll i = 0; i < prime.size(); i++)
    {
        flag = 1;
        tmp = prime[i];
        while (tmp > 0)
        {
            if (!p[tmp])
            {
                flag = 0;
                break;
            }
            tmp /= 10;
        }
        if (!flag)continue;
        tmp = prime[i];
        while (tmp > 0)
        {
            if (!p[tmp])
            {
                flag = 0;
                break;
            }
            tmp = tmp % power(10, log(tmp) - 1);
        }
        if (flag)
        {
            sum += prime[i];
        }
    }
    return sum;
}
ll Q38()
{
    ll n = 1e5, m = 0, tmp, cur = 0;
    for (ll i = 1; i < n; i++)
    {
        tmp = 0;
        vector<bool>pandigital(10, false);
        for (ll j = 1; log(tmp) < 9; j++)
        {
            tmp *= power(10, log(j * i));
            tmp += j * i;
        }
        cur = tmp;
        if (log(tmp) > 9)continue;
        while (tmp > 0)
        {
            pandigital[tmp % 10] = 1;
            tmp /= 10;
        }
        tmp = 1;
        for (ll j = 1; j < 10; j++)
        {
            if (!pandigital[j])
            {
                tmp = 0; break;
            }
        }
        if (tmp)
        {
            m = max(m, cur);
        }
    }
    return m;
}
ll Q39()
{
    vector<ll>ans(1001);
    ll x, y, z, m = 0, p;
    for (x = 1; x < 1000; x++)
    {
        for (y = 1; y < x; y++)
        {
            for (z = 1; z < y && (x + y + z) < 1001; z++)
            {
                if (x * x == y * y + z * z)ans[x + y + z]++;
            }
        }
    }
    for (ll i = 0; i < 1001; i++)
    {
        if (m < ans[i])
        {
            m = ans[i]; p = i;
        }
    }
    return p;
}
ll Q40()
{
    ll ans = 1, cur = 0, cnt = 1;
    for (ll i = 0; i < 7; i++)
    {
        while (cur < power(10, i))
        {
            cur += log(cnt);
            cnt++;
        }
        cnt--;
        ans *= (cnt / power(10, cur - power(10, i))) % 10;
        cnt++;
    }
    return ans;
}
ll Q41()
{
    ll n = 1e8, tmp;
    vector <ll> prime;
    vector<ll>div(n);
    vector<bool>p(n);
    sieve(prime, div, p, n);
    for (ll i = prime.size() - 1;; i--)
    {
        vector<bool>pandigital(10);
        tmp = prime[i];
        while (tmp > 0)
        {
            pandigital[tmp % 10] = 1;
            tmp /= 10;
        }
        tmp = 1;
        for (ll j = 1; j <= log(prime[i]); j++)
        {
            if (!pandigital[j])tmp = 0;
        }
        if (tmp)return prime[i];
    }
    return 1;
}
ll Q42()
{
    //cancer
    return 162;
}
namespace Q43//544,112 microseconds
{//took most of the code from my cses - creating substrings
    vector <string> pos;
    void fill(string s, int i, int n)
    {
        string d = "";
        if (i == n - 1)
        {
            pos.push_back(s);
            return;
        }
        fill(s, i + 1, n);
        for (int j = i + 1; j < n; j++)
        {
            if (s[i] == s[j])continue;
            d = "";
            for (int m = 0; m < n; m++)
            {
                if (m < i)d += s[m];
                if (m > i && m != j)d += s[m];
                if (m == i)
                {
                    d += s[j];
                    d += s[i];
                }
            }
            fill(d, i + 1, n);
        }
    }
    ll solve()
    {
        string s = "0123456789";
        vector <int>a(26);
        vl num;
        ll sum = 0;
        vl p = { 17,13,11,7,5,3,2 };
        for (int i = 0; i < s.size(); i++)
        {
            int b = s[i] - '0';
            a[b]++;
        }
        s = "";
        for (int i = 0; i < 26; i++)
        {
            for (int j = 0; j < a[i]; j++)
            {
                char y = i + '0';
                s += y;
            }
        }
        fill(s, 0, s.size());
        for (string str : pos)
        {
            if (str[0] == '0')continue;
            num.push_back(string_to_int(str));
        }
        for (ll n : num)
        {
            ll flag = 1;
            for (ll i = 0; i < 7; i++)
            {
                if (n / power(10, i) % 1000 % p[i] != 0)
                {
                    flag = 0; break;
                }
            }
            if (flag)sum += n;
        }
        return sum;;
    }
}
namespace Q44
{
    //in python
}
ll Q45()
{
    vector<ll>tri = { 0 }, pent = { 0 }, hex = { 0 };
    for (ll i = 1; i < 1e6; i++)
    {
        tri.push_back((i * i + i) / 2);
        pent.push_back((3 * i * i - i) / 2);
        hex.push_back(2 * i * i - i);
    }
    for (ll i = 286; i < 1e6; i++)
    {
        if (find(pent, 0, 1e6, tri[i]) != -1 && find(hex, 0, 1e6, tri[i]) != -1)return tri[i];
    }
}
namespace Q46//6,642 microseconds
{
    ll solve()
    {
        ll n = 1e6;
        vl div(n), prime;
        vector<bool>p(n);
        sieve(prime, div, p, n);
        for (ll i = 3; i < n; i += 2)
        {
            if (p[i])continue;
            ll flag = 1;
            for (ll j = sqrt(i / 2); j > 0; j--)
            {
                if (p[i - 2 * j * j])
                {
                    flag = 0; break;
                }
            }
            if (flag)return i;
        }
        return -1;
    }
}
namespace Q47//8,583 microseconds
{
    ll four(ll n, vl& div)
    {
        ll cnt = 0;
        while (n > 1)
        {
            ll tmp = div[n];
            cnt++;
            while (tmp == div[n])n /= tmp;
        }
        if (cnt > 3)return 1;
        return 0;
    }
    ll solve()
    {
        ll n = 1e6;
        vl div(n), prime;
        vector<bool>q(n);
        sieve(prime, div, q, n);
        for (ll i = 2; i < n - 3; i++)
        {
            if (i == 134043)
            {
                ll c = 2;
            }
            ll flag = 1;
            for (ll j = 0; j < 4; j++)
            {
                if (!four(i + j, div))
                {
                    flag = 0; break;
                }
            }
            if (flag)return i;
        }
        return -1;
    }
}
namespace Q48//2,109 microseconds
{
    ll solve()
    {
        ll n = 1e3, p = 1e10, s = 0;
        for (ll i = 1; i <= 1e3; i++)
        {
            s += power(i, i, p, [](ll l1, ll l2, ll m1) {return power(l1, l2, m1, [](ll x, ll y, ll z) {return (x + y) % z; }, 0); }, 1);
            s %= p;
        }
        return s;
    }
}
namespace Q49//425,482 microseconds
{
    ll solve()
    {
        ll n = 1e4;
        vl div(n), prime;
        vector<bool>p(n);
        sieve(prime, div, p, n);
        for (ll i = 1e3; i < n; i++)
        {
            if (!p[i])continue;
            for (ll j = 1; j <= (n - i) / 2; j++)
            {
                if (i == 1487 && j == 3330)continue;//the one we already know
                if (!p[i + j] || !p[i + 2 * j])continue;
                vl tmp = { i,i + j,i + 2 * j };
                vvl v(3);
                for (ll k = 0; k < 3; k++)
                {
                    for (ll m = 0; m < 4; m++)
                    {
                        v[k].push_back(tmp[k] % 10); tmp[k] /= 10;
                    }
                    sort(v[k].begin(), v[k].end());
                }
                if (v[0] != v[1] || v[0] != v[2])continue;
                return power(10, 8) * i + power(10, 4) * (i + j) + i + 2 * j;
            }
        }
        return -1;
    }
}
ll Q50()
{
    ll n = 1e6;
    vector<ll> prime;
    vector<ll>a(n);
    vector<bool>b(n);
    pl ans = { -1,-1 };
    sieve(prime, a, b, n);
    vector<ll>psum = { 0 };
    for (ll i : prime)
    {
        psum.push_back(psum[psum.size() - 1] + i);
    }
    for (ll i = 1; i < psum.size(); i++)
    {
        for (ll j = 0; j < i; j++)
        {
            if (find(prime, 0, prime.size() - 1, psum[i] - psum[j]) != -1)
            {
                if (ans.second < i - j)
                {
                    ans = { psum[i] - psum[j],i - j };
                }
            }
        }
    }
    return ans.first;
}
namespace Q51//382,807 microseconds
{
    ll solve()//kinda disgusting, also i assume the answer will be less than 1mil so only check replacing 3 digits
    {
        ll n = 1e6, tmp;
        vl p, div(n), seen(n);
        vector<bool>c(n);
        sieve(p, div, c, n);
        for (ll i = 1000; i < n; i++)
        {
            for (ll j1 = 0; j1 < log(i); j1++)
            {
                for (ll j2 = 0; j2 < j1; j2++)
                {
                    for (ll j3 = 0; j3 < j2; j3++)
                    {
                        ll flag = 1;
                        if (j1 == log(i) - 1)flag = 0;
                        ll cnt = 0;
                        tmp = i;
                        tmp -= (tmp / power(10, j1) % 10 * power(10, j1));
                        tmp -= (tmp / power(10, j2) % 10 * power(10, j2));
                        tmp -= (tmp / power(10, j3) % 10 * power(10, j3));
                        if (seen[tmp])continue;
                        for (ll x = 0; x < 10; x++)
                        {
                            ll a = tmp + x * (power(10, j1) + power(10, j2) + power(10, j3));
                            if (c[a])cnt++;
                            seen[a] = 1;
                        }
                        if (!flag && c[tmp])cnt--;
                        if (cnt > 7)
                        {
                            tmp = i;
                            tmp -= (tmp / power(10, j1) % 10 * power(10, j1));
                            tmp -= (tmp / power(10, j2) % 10 * power(10, j2));
                            tmp -= (tmp / power(10, j3) % 10 * power(10, j3));
                            for (ll x = 1; x < 10; x++)
                            {
                                ll a = tmp + x * (power(10, j1) + power(10, j2) + power(10, j3));
                                if (c[a])return a;
                            }
                        }
                    }
                }
            }
        }
        return 0;
    }
}
ll Q52()
{
    for (ll i = 1; i < 1e8; i++)
    {
        vector<vector<ll>>d;
        for (ll j = i; j < 7 * i; j += i)
        {
            vector<ll>v;
            ll tmp = j;
            while (tmp != 0)
            {
                v.push_back(tmp % 10);
                tmp /= 10;
            }
            sort(v.begin(), v.end());
            d.push_back(v);
        }
        ll flag = 1;
        for (ll i = 0; i < 5; i++)
        {
            for (ll j = 0; j < d[i].size(); j++)
            {
                if (d[i].size() != d[i + 1].size() || d[i][j] != d[i + 1][j])
                {
                    flag = 0;
                    j = 1000000;
                    i = 5;
                }
            }
        }
        if (flag)return i;
    }
}
namespace Q53
{
    //python
}
namespace Q54//15,025 microseconds
{
    pl convert(string a)
    {
        pl b;
        b.first = a[0] - '0'; b.second = a[1] - 'A';
        if (a[0] == 'T')b.first = 10;
        if (a[0] == 'J')b.first = 11;
        if (a[0] == 'Q')b.first = 12;
        if (a[0] == 'K')b.first = 13;
        if (a[0] == 'A')b.first = 14;
        return b;
    }
    bool royalflush(vector<pl>& v)
    {
        for (ll i = 0; i < 5; i++)
        {
            if (v[i].first != 10 + i)return 0;
        }
        for (ll i = 0; i < 4; i++)
        {
            if (v[i].second != v[i + 1].second)return 0;
        }
        return 1;
    }
    bool straightflush(vector<pl>& v)
    {
        for (ll i = 0; i < 4; i++)
        {
            if (v[i + 1].first - v[i].first != 1)return 0;
            if (v[i].second != v[i + 1].second)return 0;
        }
        return 1;
    }
    bool four(vector<pl>& v)
    {
        ll f1 = 1, f2 = 1;
        for (ll i = 0; i < 3; i++)
        {
            if (v[i].first != v[i + 1].first)f1 = 0;
        }
        for (ll i = 1; i < 4; i++)
        {
            if (v[i].first != v[i + 1].first)f2 = 0;
        }
        if (f1 || f2)return 1;
        return 0;
    }
    bool fullhouse(vector<pl>& v)
    {
        ll f1 = 1, f2 = 1;
        if (v[0].first != v[1].first)f1 = 0;
        for (ll i = 2; i < 4; i++)
        {
            if (v[i].first != v[i + 1].first)f1 = 0;
        }
        if (v[3].first != v[4].first)f2 = 0;
        for (ll i = 0; i < 2; i++)
        {
            if (v[i].first != v[i + 1].first)f2 = 0;
        }
        if (f1 || f2)return 1;
        return 0;
    }
    bool flush(vector<pl>& v)
    {
        for (ll i = 0; i < 4; i++)
        {
            if (v[i].second != v[i + 1].second)return 0;
        }
        return 1;
    }
    bool straight(vector<pl>& v)
    {
        for (ll i = 0; i < 4; i++)
        {
            if (v[i + 1].first - v[i].first != 1)return 0;
        }
        return 1;
    }
    bool three(vector<pl>& v)
    {
        vl p(3, 1);
        for (ll i = 0; i < 3; i++)
        {
            for (ll j = i; j < i + 2; j++)
            {
                if (v[j].first != v[j + 1].first)p[i] = 0;
            }
        }
        for (ll i = 0; i < 3; i++)if (p[i])return 1;
        return 0;
    }
    bool twopair(vector<pl>& v)
    {
        ll cnt = 0;
        for (ll i = 0; i < 4; i++)
        {
            if (v[i].first == v[i + 1].first)cnt++;
        }
        if (cnt > 1)return 1;
        return 0;
    }
    bool pair(vector<pl>& v)
    {
        for (ll i = 0; i < 4; i++)
        {
            if (v[i].first == v[i + 1].first)return 1;
        }
        return 0;
    }
    ll value(vector<pl>& v)
    {
        if (royalflush(v))return 10;
        if (straightflush(v))return 9;
        if (four(v))return 8;
        if (fullhouse(v))return 7;
        if (flush(v))return 6;
        if (straight(v))return 5;
        if (three(v))return 4;
        if (twopair(v))return 3;
        if (pair(v))return 2;
        return 1;
    }
    ll solve()//only had 1 thing to debug!!!!!!
    {
        ifstream F("C:\\Users\\user\\source\\repos\\computer science\\project euler\\Q54 input.txt");
        ll a, b, cnt = 0, s1 = 0, s2 = 0;
        vector<string>t1(5), t2(5);
        vector<pl>h1(5), h2(5);
        while (F >> t1[0])
        {
            for (ll i = 1; i < 5; i++)F >> t1[i];
            for (ll i = 0; i < 5; i++)F >> t2[i];
            for (ll i = 0; i < 5; i++)h1[i] = convert(t1[i]);
            for (ll i = 0; i < 5; i++)h2[i] = convert(t2[i]);
            sort(h1.begin(), h1.end()); sort(h2.begin(), h2.end());
            a = value(h1), b = value(h2);
            if (a > b)cnt++;
            if (a == b)
            {
                if (a == 9)
                {
                    if (h1[0].first > h2[0].first)cnt++;
                }
                else if (a == 8 || a == 7)
                {
                    if (h1[2].first > h2[2].first)cnt++;
                    if (h1[2].first == h2[2].first)
                    {
                        s1 = s2 = 0;
                        for (pl i : h1)s1 += i.first;
                        for (pl i : h2)s2 += i.first;
                        if (s1 > s2)cnt++;
                    }
                }
                else if (a == 6)
                {
                    for (ll i = 4; i >= 0; i--)
                    {
                        if (h1[i].first > h2[i].first)
                        {
                            cnt++; break;
                        }
                        if (h1[i].first < h2[i].first)break;
                    }
                }
                else if (a == 5)
                {
                    if (h1[0].first > h2[0].first)cnt++;
                }
                else if (a == 4)
                {
                    if (h1[2].first > h2[2].first)cnt++;
                    if (h1[2].first == h2[2].first)
                    {
                        vl p1, p2;
                        for (ll i = 0; i < 5; i++)if (h1[i].first != h1[2].first)p1.push_back(h1[i].first);
                        for (ll i = 0; i < 5; i++)if (h2[i].first != h2[2].first)p2.push_back(h2[i].first);
                        if (p1[0] < p1[1])swap(p1[0], p1[1]); if (p2[0] < p2[1])swap(p2[0], p2[1]);
                        if (p1[0] > p2[0])cnt++;
                        else if (p1[0] == p2[0])if (p1[1] > p2[1])cnt++;
                    }
                }
                else if (a == 3)
                {
                    s1 = s2 = 0;
                    vl p1, p2;
                    for (ll i = 0; i < 4; i++)
                    {
                        if (h1[i].first == h1[i + 1].first)p1.push_back(h1[i].first);
                    }
                    for (ll i = 0; i < 4; i++)
                    {
                        if (h2[i].first == h2[i + 1].first)p2.push_back(h2[i].first);
                    }
                    if (p1[0] < p1[1])swap(p1[0], p1[1]); if (p2[0] < p2[1])swap(p2[0], p2[1]);
                    if (p1[0] > p2[0])cnt++;
                    else if (p1[0] == p2[0])
                    {
                        if (p1[1] > p2[1])cnt++;
                        else if (p1[1] == p2[1])
                        {
                            for (pl i : h1)s1 += i.first;
                            for (pl i : h2)s2 += i.first;
                            if (s1 > s2)cnt++;
                        }
                    }
                }
                else if (a == 2)
                {
                    vl p1, p2;
                    for (ll i = 0; i < 4; i++)
                    {
                        if (h1[i].first == h1[i + 1].first)p1.push_back(h1[i].first);
                    }
                    for (ll i = 0; i < 4; i++)
                    {
                        if (h2[i].first == h2[i + 1].first)p2.push_back(h2[i].first);
                    }
                    if (p1[0] > p2[0])cnt++;
                    else if (p1[0] == p2[0])
                    {
                        vl p3, p4;
                        for (ll i = 0; i < 5; i++)
                        {
                            if (h1[i].first != p1[0])p3.push_back(h1[i].first);
                        }
                        for (ll i = 0; i < 5; i++)
                        {
                            if (h2[i].first != p2[0])p4.push_back(h2[i].first);
                        }
                        sort(p3.begin(), p3.end()); sort(p4.begin(), p4.end());
                        for (ll i = 2; i >= 0; i--)
                        {
                            if (p3[i] > p4[i])
                            {
                                cnt++; break;
                            }
                            if (p3[i] < p4[i])break;
                        }
                    }
                }
                else
                {
                    for (ll i = 4; i >= 0; i--)
                    {
                        if (h1[i].first > h2[i].first)
                        {
                            cnt++; break;
                        }
                        if (h1[i].first < h2[i].first)break;
                    }
                }
            }
        }
        return cnt;
    }
}
namespace Q55
{
    //in python
}
namespace Q56
{
    //in python
}
namespace Q57
{
    //in python
}
ll Q58()
{
    ll n = 1e9;
    ll cnt = 0;
    vector<ll>prime, div(n); vector<bool>p(n);
    sieve(prime, div, p, n);
    for (ll m = 3; m < sqrt(n); m += 2)
    {
        cnt += p[m * m - m + 1] + p[m * m - 2 * m + 2] + p[m * m - 3 * m + 3];
        if (10 * cnt < 2 * m - 1)return m;
    }
    return 1;
}
namespace Q60//48,298,309 microseconds. the code gives the write answer but im pretty sure it has bugs so should write it again
{//disgusting code
    ll inc(vector<pl>& v, vvl& p, vl& prime, ll sum)//inc runs over all possible ordered primes such that their sum isnt more than sum
    {
        if (v.size() * v[0].first > sum || v[0].first >= p.size())return 0;
        ll a = prime[v[0].first];
        for (ll i = 1; i < v.size(); i++)//immediately change value
        {
            if (a + (v.size() - i) * prime[v[i].first] > sum)
            {
                i--;
                ll flag = 0;
                v[i] = { i ? p[v[i - 1].first][v[i].second + 1] : v[0].first + 1,i ? v[i].second + 1 : 0 };
                for (ll j = i + 1; j < v.size(); j++)
                {
                    if (flag || p[v[j - 1].first].size() == 0)
                    {
                        v[j] = { 0,1e18 };
                        flag = 1;
                    }
                    else v[j] = { p[v[j - 1].first][0],0 };
                }
                if (flag)inc(v, p, prime, sum);
                return 1;
            }
            a += prime[v[i].first];
        }
        for (ll i = v.size() - 1; ~i; i--)
        {
            ll x = i ? p[v[i - 1].first].size() - 1 - v[i].second : 1;
            if (x > 0)
            {
                ll flag = 0;
                v[i] = { i ? p[v[i - 1].first][v[i].second + 1] : v[0].first + 1,i ? v[i].second + 1 : 0 };
                for (ll j = i + 1; j < v.size(); j++)
                {
                    if (flag || p[v[j - 1].first].size() == 0)
                    {
                        v[j] = { 0,1e18 };
                        flag = 1;
                    }
                    else v[j] = { p[v[j - 1].first][0],0 };
                }
                if (flag)inc(v, p, prime, sum);
                return 1;
            }
        }
    }
    bool conc(vl& v)//checks concetration 
    {
        ll a = 0; for (ll i : v)a += i;
        ll x;
        vl s(v.size());
        for (ll i = 0; i < v.size(); i++)s[i] = log(v[i], 10);
        for (ll i = 0; i < v.size(); i++)for (ll j = 0; j < v.size(); j++)
        {
            if (i == j || i == j + 1 || i == j - 1)continue;
            x = v[i] + power(10, s[i]) * v[j];
            if (!is_prime(x))
            {
                return 0;
            }
        }
        return 1;
    }
    ll solve()
    {
        ll k = 20000, s = 5, sum = 26033, tmp;//already found such sum from running the code earlier
        vl div(k), prime;
        vector<bool>q(k);
        vvl p(2);//without 2
        p[0].push_back(3); p[1].push_back(3);
        sieve(prime, div, q, k);
        for (ll i = 2; i < prime.size(); i++)
        {
            if (prime[i] % 3 == 1)p[0].push_back(prime[i]);
            else p[1].push_back(prime[i]);
        }
        for (ll i : {0, 1})
        {
            vector<pl> v(s);
            vvl pair(p[i].size());
            for (ll m = 0; m < pair.size(); m++)
            {
                for (ll j = m + 1; j < pair.size(); j++)
                {
                    if (!is_prime(p[i][m] + power(10, log(p[i][m])) * p[i][j]))continue;
                    if (!is_prime(p[i][j] + power(10, log(p[i][j])) * p[i][m]))continue;
                    pair[m].push_back(j);
                }
            }
            v[0] = { 0,0 };
            for (ll i = 1; i < s; i++)v[i] = { pair[v[i - 1].first][0],0 };
            v[s - 1].second = -1;
            while (inc(v, pair, p[i], sum))
            {
                vl c(s);
                for (ll j = 0; j < s; j++)
                {
                    if (v[j].first >= p[i].size())return -1;
                    c[j] = p[i][v[j].first];
                }
                if (conc(c))
                {
                    tmp = 0;
                    for (ll x : c)
                    {
                        tmp += x;
                    }
                    cout << tmp;
                    if (sum > tmp)sum = tmp;
                    cout << '\n';
                }
            }
        }
        return 1;
    }
}
namespace Q62//3,381 microseconds
{
    ll big(ll n)
    {
        vl d(10);
        while (n)
        {
            d[n % 10]++; n /= 10;
        }
        ll ans = 0;
        for (ll i = 9; ~i; i--)
        {
            for (ll j = 0; j < d[i]; j++)
            {
                ans *= 10; ans += i;
            }
        }
        return ans;
    }
    ll solve()
    {

        ll n = 1e4 - 1, cur = 1;//it works because we run only on cubes less than 10^12 so this are our only options (if there is at least 1)
        vl v, num;
        for (ll i = 1; i <= n; i++)
        {
            v.push_back(big(i * i * i));
        }
        sort(v.begin(), v.end());
        for (ll i = 0; i < v.size() - 1; i++)
        {
            if (v[i] == v[i + 1])cur++;
            else
            {
                if (cur == 5)//becuase we checked there only 2 such numbers there is the funny
                {
                    num.push_back(v[i]);
                }
                cur = 1;
            }
        }
        for (ll i = 1; i <= n; i++)
        {
            cur = big(i * i * i);
            for (ll j = 0; j < num.size(); j++)
            {
                if (cur == num[j])return i * i * i;
            }
        }
        return -1;
    }
}
namespace Q63
{
    //in python
}
namespace Q67//1,139 microseconds
{
    ll solve()
    {
        ifstream F("C:\\Users\\user\\source\\repos\\computer science\\project euler\\Q67 input.txt");
        string s;
        ll m = 0;
        vvl dp(100);
        vvl v(100);
        for (ll i = 0; i < 100; i++)
        {
            for (ll j = 0; j <= i; j++)
            {
                F >> s;
                v[i].push_back(string_to_int(s));
                dp[i].push_back(0);
            }
        }
        dp[0][0] = v[0][0];
        for (ll i = 1; i < 100; i++)
        {
            for (ll j = 0; j <= i; j++)
            {
                if (j == 0)dp[i][j] = dp[i - 1][j] + v[i][j];
                else if (j == i)dp[i][j] = dp[i - 1][j - 1] + v[i][j];
                else dp[i][j] = max(dp[i - 1][j], dp[i - 1][j - 1]) + v[i][j];
                m = max(m, dp[i][j]);
            }
        }
        return m;
    }
}
namespace Q68
{
    ll solve()//you want outer circle to be 10,9,8,7,6 and than case work (you want 6,5,3) gives this as the best (and only) solution
    {
        return 6531031914842725;
    }
}
ll Q69()//you want as many different small primes, which means multiple from 2 onward
{
    return 2 * 3 * 5 * 7 * 11 * 13 * 17;
}
namespace Q70//633,173 microseconds
{
    ll per(ll a, ll b)
    {
        ll c = log(a);
        if (log(a) != log(b))return 0;
        vl na, nb;
        for (ll i = 0; i < c; i++)
        {
            na.push_back(a % 10), nb.push_back(b % 10);
            a /= 10, b /= 10;
        }
        sort(na.begin(), na.end()), sort(nb.begin(), nb.end());
        if (na == nb)return 1;
        return 0;
    }
    ll solve()
    {
        ll n = 1e7, tmp, ans = -1;
        double m = 1e7;
        vector <ll> prime;
        vector<ll>div(n);
        vector<bool>c(n);
        vector<ll>f(n); f[1] = 1;
        sieve(prime, div, c, n);
        for (ll i = 2; i < n; i++)
        {
            tmp = i;
            tmp /= div[tmp];
            if (div[tmp] == div[i])f[i] = f[tmp] * div[i];
            else f[i] = f[tmp] * (div[i] - 1);
            if (m > double(i) / f[i] && per(i, f[i]))
            {
                m = double(i) / f[i], ans = i;
            }
        }
        return ans;
    }
}
namespace Q71
{
    pl big(pl a, pl b)
    {
        if (a.first * b.second > a.second * b.first)return a;
        return b;
    }
    ll solve()
    {
        ll n = 1e6;
        pl m = { -1,1 };
        vector<pl>v;
        for (ll i = 3; i <= n; i++)
        {
            v.push_back({ (3 * i - 1) / 7,i });
        }
        for (pl x : v)m = big(m, x);
        return m.first / gcd(m.first, m.second);
    }
}
namespace Q72
{
    ll solve()
    {
        ll n = 1e6 + 1, tmp, sum = 0;
        vector <ll> prime;
        vector<ll>div(n);
        vector<bool>c(n);
        vector<ll>f(n); f[1] = 1;
        sieve(prime, div, c, n);
        for (ll i = 2; i < n; i++)
        {
            tmp = i;
            tmp /= div[tmp];
            if (div[tmp] == div[i])f[i] = f[tmp] * div[i];
            else f[i] = f[tmp] * (div[i] - 1);
            sum += f[i];
        }
        return sum;
    }
}
namespace Q73
{
    ll solve()
    {
        ll n = 12e3, sum = 0;
        vl add(n + 1);
        for (ll i = 5; i <= n; i++)
        {
            add[i] = (i - 1) / 2 - i / 3 - add[i];
            if (add[i])for (ll j = 2 * i; j <= n; j += i)add[j] += add[i];
            sum += add[i];
        }
        return sum;
    }
}
namespace Q74//2,955,342 microseconds
{
    ll fd(ll n, vl& f)
    {
        ll sum = 0;
        while (n)
        {
            sum += f[n % 10];
            n /= 10;
        }
        return sum;
    }
    ll solve()//brute force
    {
        ll n = 1e6, k = 60, cnt = 0, tmp, flag;
        vl fact(10); fact[0] = 1; for (ll i = 1; i < 10; i++)fact[i] = i * fact[i - 1];
        for (ll i = 3; i < n; i++)
        {
            set<ll>s; s.insert(i);
            tmp = i, flag = 1;
            for (ll j = 1; j < k; j++)
            {
                tmp = fd(tmp, fact);
                s.insert(tmp);
                if (s.size() != j + 1)
                {
                    flag = 0; break;
                }
            }
            if (flag)cnt++;

        }
        return cnt;
    }

}
namespace Q75
{
    ll solve()
    {
        ll n = 15e5, sum, cnt = 0;
        vl v(n + 1);
        for (ll i = 1; 2 * i * (i + 1) <= n; i++)for (ll j = 1 + (i & 1); j < i; j += 2)
        {
            sum = 2 * i * (i + j); if (sum > n)break;
            if (gcd(i, j) != 1)continue;
            for (ll k = sum; k <= n; k += sum)v[k]++;
        }
        for (ll x : v)if (x == 1)cnt++;
        return cnt;
    }
}
namespace Q76
{
    ll solve()//imported from my cses cause im lazy :)
    {
        ll n = 99, x = 100, tmp;
        vector<ll>coin(n);
        for (ll i = 0; i < n; i++)coin[i] = i + 1;
        vector<ll>window(x + 1);
        sort(coin.begin(), coin.end());
        for (int i = 0; i < n; i++)
        {
            window[coin[i]]++;
            for (int j = coin[i] + 1; j <= x; j++)
            {
                window[j] += window[j - coin[i]];
            }
        }
        return window[x];
    }
}
namespace Q77
{
    ll c(vl coin, ll n, ll x)//basically 76 with binary search <0/ clearly f(x+2)>f(x) and so 2 binary searches
    {
        ll tmp;
        vector<ll>window(x + 1);
        for (int i = 0; i < n; i++)
        {
            if (coin[i] > x)break;
            window[coin[i]]++;
            for (int j = coin[i] + 1; j <= x; j++)
            {
                window[j] += window[j - coin[i]];
            }
        }
        return window[x];
    }
    ll solve()
    {
        ll n = 100, ans;
        vector <ll> prime;
        vector<ll>div(n);
        vector<bool>p(n);
        sieve(prime, div, p, n);
        ll low = 1, high = n / 2 - 1, m;
        while (low < high - 1)
        {
            m = (low + high) / 2;
            if (c(prime, prime.size(), 2 * m + 1) < 5000)low = m;
            else high = m - 1;
        }
        ans = 2 * low + 3;
        low = 1, high = n / 2 - 1, m;
        while (low < high - 1)
        {
            m = (low + high) / 2;
            if (c(prime, prime.size(), 2 * m) < 5000)low = m;
            else high = m - 1;
        }
        return min(2 * low + 2, ans);
    }
}
string Q79()
{
    string s = "";
    vector<vector<ll>>v(50, vector<ll>(3));
    vector<ll>in(50);
    ll x;
    for (ll i = 0; i < 50; i++)
    {
        cin >> in[i];
    }
    for (ll k = 0; k < 1000; k++)
    {
        vector<bool>c(10);
        vector<bool>seen(10);
        for (ll i = 0; i < 50; i++)
        {
            for (ll j = 0; j < log(in[i]) - 1; j++)
            {
                if (j == 0)
                {
                    seen[in[i] % 10] = c[in[i] % 10] = 1;
                }
                if (j == 1)
                {
                    seen[(in[i] / 10) % 10] = c[(in[i] / 10) % 10] = 1;
                }
            }
            seen[in[i] / power(10, log(in[i]) - 1)] = 1;
        }
        for (ll j = 0; j < 10; j++)
        {
            if (seen[j] && !c[j])
            {
                for (ll m = 0; m < 50; m++)
                {
                    if ((in[m] / power(10, log(in[m]) - 1)) % 10 != j)continue;
                    if (log(in[m]) == 3)in[m] %= 100;
                    else if (log(in[m]) == 2)in[m] %= 10;
                    else in[m] = 0;
                }
                s += to_string(j);
                break;
            }
        }
    }
    return s;
}
namespace Q80
{
    //in python
}
namespace Q81//900 microseconds
{
    ll solve()
    {
        ifstream F("C:\\Users\\user\\source\\repos\\computer science\\project euler\\Q81 input.txt");
        vvl v;
        string s, t;
        ll j = -1;
        while (getline(F, s))
        {
            j++;
            v.push_back({});
            ll i = 0; t = "";
            while (i < s.size())
            {
                if (s[i] != ',')t += s[i];
                else
                {
                    v[j].push_back(string_to_int(t));
                    t = "";
                }
                i++;
            }
            v[j].push_back(string_to_int(t));
        }
        ll n = v[0].size();
        vvl dp(n + 1, vl(n + 1, 1e9)); dp[n - 1][n - 1] = v[n - 1][n - 1];
        for (ll i = n - 1; ~i; i--)for (ll j = n - 1; ~j; j--)
        {
            if (i == n - 1 && j == n - 1)continue;
            dp[i][j] = min(dp[i + 1][j], dp[i][j + 1]) + v[i][j];
        }
        return dp[0][0];
    }
}
namespace Q85
{
    ll calc(ll n, ll m)
    {
        return (n * (n + 1)) * (m * m + m) / 4;
    }
    ll solve()//there are n(n+1)m(m+1)/4 options for m,n
    {
        ll n = 2e6, cl = n, tmp, ans;
        for (ll i = 1; i < 3e3; i++)for (ll j = i; j < 3e3; j++)
        {
            tmp = calc(i, j);
            if (abs(tmp - n) < cl)
            {
                cl = abs(tmp - n);
                ans = i * j;
                cout << i << "  " << j << '\n';
            }
        }
        return ans;
    }
}
namespace Q86//5,083 microseconds
{
    ll solve()//for a>=b>=c min length is sqrt(a^2+(b+c)^2)
    {
        ll sum = 0, bound = 1e6, a;
        for (a = 1; sum < bound; a++)
        {
            for (ll s = 2; s < 2 * a; s++)
            {
                ll d = sqrt(a * a + s * s);
                if (d * d == a * a + s * s)
                {
                    sum += min(a, s - 1) - (s - 1) / 2;
                }
            }
        }
        return a - 1;
    }
}
namespace Q87
{
    ll solve()
    {
        ll n = 5e7, m = sqrt(n), tmp, cnt = 0;
        vl p, div(m), ans(n + 1);
        vector<bool>c(m);
        sieve(p, div, c, m);
        for (ll i = 0; i < p.size(); i++)
        {
            if (p[i] * p[i] >= n)break;
            for (ll j = 0; j < p.size(); j++)
            {
                if (p[i] * p[i] + p[j] * p[j] * p[j] >= n)break;
                for (ll k = 0; k < p.size(); k++)
                {
                    tmp = p[i] * p[i] + p[j] * p[j] * p[j] + power(p[k], 4);
                    if (tmp >= n)break;
                    ans[tmp]++;
                }
            }
        }
        for (ll x : ans)if (x)cnt++;
        return cnt;
    }
}
namespace Q88//127,000 microseconds
{
    ll solve()
    {
        ll n = 12e3, x = 2 * n, m, sum = 0;
        vvl d(x + 1); d[1] = { 1 };
        vl pos(x + 1, 1e9); pos[1] = { 1 };
        set<ll>ans;
        for (ll i = 2; i <= x; i++)
        {
            d[i] = divisors::find_divisors(i);
        }
        for (ll i = 2; i <= x; i++)sort(d[i].begin(), d[i].end());
        for (ll i = 2; i <= x; i++)
        {
            stack<pl>s;
            ll prod = 1;
            ll summed = 0;
            ll upper_bound = i;
            //the divisors are montonic decreasing
            s.push({ i,1 });
            while (!s.empty()) {
                //current elements are the partition
                //prod is their product
                if (prod == i) {
                    if (prod - summed + s.size() - 1 == 65)
                    {
                        ll asd = 2;
                    }
                    pos[prod - summed + s.size() - 1] = min(i, pos[prod - summed + s.size() - 1]);
                    prod /= s.top().first;
                    summed -= s.top().first;
                    s.pop();
                    upper_bound = s.top().first;
                    continue;
                }
                int idx = s.top().second;
                if (idx < d[i / prod].size() && d[i / prod][idx] <= upper_bound) {
                    auto pp = s.top();
                    s.pop();
                    s.push({ pp.first,pp.second + 1 });
                    int div = d[i / prod][idx];
                    s.push({ div,1 });
                    upper_bound = div;
                    prod *= div;
                    summed += div;
                }
                else {
                    prod /= s.top().first;
                    summed -= s.top().first;
                    s.pop();
                    if (!s.empty())upper_bound = s.top().first;
                }
            }
        }
        for (ll i = 2; i <= n; i++)ans.insert(pos[i]);
        for (auto a : ans)if (a != 1e9)sum += a;
        return sum;
    }
}
namespace Q91//Time taken: 719 microseconds
{
    double norm(double x, double y)
    {
        return sqrt(x * x + y * y);
    }
    ll solve()
    {
        ll n = 50, cnt = (n) * (n);//this options are for when (0,0) is the right angle
        double a, b, c, d, s = 0.001;
        for (double x = 0; x <= n; x++)for (double y = 0; y <= n; y++)
        {
            if (x == 0 && y == 0)continue;
            a = norm(x / 2, y / 2) * norm(x / 2, y / 2);
            for (double z = 0; z <= n; z++)
            {
                b = a - (z - x / 2) * (z - x / 2);
                if (b < -s)continue;
                if (b <= 0)b = 0;
                else b = sqrt(b);
                c = b + y / 2; if (c <= n + s && c + -int(c + s * s) < s)cnt++;
                d = y / 2 - b; if (d >= -s && d + -int(d + s * s) < s)cnt++;
                if (abs(b) < s && c - int(c + s * s) < s)cnt--;
            }
            cnt -= 2;//for (0,0) and(x,y)
        }
        return cnt;
    }
}
namespace Q92//155,567 microseconds
{
    ll solve()
    {
        ll n = 1e7;
        vl parent(n), size(n);
        for (ll i = 1; i < n; i++)
        {
            dsu::make_set(parent, size, i);
        }
        for (ll i = 1; i < n; i++)
        {
            ll tmp = i, s = 0;
            while (tmp)
            {
                s += tmp % 10 * (tmp % 10); tmp /= 10;
            }
            dsu::union_sets(i, s, parent, size);
        }
        return size[dsu::find_set(89, parent)];
    }
}
namespace Q94//954,213 microseconds
{
    ll solve()//herons formula+brute force
    {
        ll n = 1e9, sum = 0, tmp, d;
        for (ll x = 2; 3 * x + 1 <= n; x++)
        {
            tmp = ((x - 1) * (3 * x + 1));
            d = sqrt(tmp);
            if (d * d == tmp)sum += 3 * x + 1;
        }
        for (ll x = 3; 3 * x - 1 <= n; x++)
        {
            tmp = ((x + 1) * (3 * x - 1));
            d = sqrt(tmp);
            if (d * d == tmp)sum += 3 * x - 1;
        }
        return sum;
    }
}
namespace Q95//1,064,254 microseconds
{
    ll solve()//yet another brute force
    {
        ll n = 1e6, chain = 0, cnt, ans = 1e6, tmp;
        vl div(n + 1), prime;
        vector<bool>p(n + 1);
        vl sum(n + 1);
        sieve(prime, div, p, n + 1);
        for (ll i = 2; i <= n; i++)
        {
            if (p[i])
            {
                sum[i] = 1;
                continue;
            }
            tmp = div[i], cnt = 0;
            while (i % tmp == 0)
            {
                tmp *= div[i];
                cnt++;
            }
            tmp /= div[i];
            sum[i] = (sum[i / tmp] + i / tmp) * (power(div[i], cnt + 1) - 1) / (div[i] - 1) - i;
        }
        for (ll i = 2; i < n; i++)
        {
            set<ll>s; s.insert(i);
            cnt = 1, tmp = i;
            while (sum[tmp] <= n && sum[tmp] != 1)
            {
                tmp = sum[tmp];
                cnt++;
                s.insert(tmp);
                if (s.size() != cnt)
                {
                    if (tmp != i)
                    {
                        break;
                    }
                    if (chain == cnt - 1)
                    {
                        chain = cnt - 1, ans = min(ans, i);
                        break;
                    }
                    if (chain < cnt - 1)
                    {
                        chain = cnt - 1, ans = i;
                        break;
                    }
                }
            }
        }
        return ans;
    }
}
namespace Q97//218 microseconds
{
    ll solve()
    {
        ll p = 1e10, ans = 0;
        ans += power(2, 7830457, p, [](ll l1, ll l2, ll m) {return power(l1, l2, m, [](ll x, ll y, ll z) {return (x + y) % z; }, 0); }, 1);
        return (ans * 28433 + 1) % p;
    }
}
namespace Q99//624 microseconds
{
    ll solve()
    {
        ifstream F("C:\\Users\\user\\source\\repos\\computer science\\project euler\\Q99 input.txt");
        vector <pl> v;
        string s, t;
        ll j = -1;
        while (getline(F, s))
        {
            j++;
            ll i = 0; t = "";
            while (i < s.size())
            {
                if (s[i] != ',')t += s[i];
                else
                {
                    v.push_back({ string_to_int(t),0 });
                    t = "";
                }
                i++;
            }
            v[v.size() - 1] = { v[v.size() - 1].first,string_to_int(t) };
        }
        ll idx = 0;
        double m = -1;
        for (ll i = 0; i < v.size(); i++)
        {
            if (v[i].second * log10(v[i].first) > m)
            {
                m = v[i].second * log10(v[i].first); idx = i + 1;
            }
        }
        return idx;
    }
}
namespace Q100
{
    //solved in python
}
namespace Q102//1,666 microseconds
{
    ll side(ll x1, ll y1, ll x2, ll y2, ll x3, ll y3)
    {
        ll tmp1, tmp2;
        tmp1 = x2 - x1; tmp2 = y2 - y1;
        ll c = (x3 - x1) * tmp2 - (y3 - y1) * tmp1;
        if (c == 0)
        {
            return 0;
        }
        if (c < 0)return -1;
        if (c > 0)return 1;
    }
    bool boundry(pl p, vector<pair<ll, ll>> v)
    {
        for (int i = 0; i < v.size() - 1; i++)
        {
            if (side(v[i].first, v[i].second, v[i + 1].first, v[i + 1].second, p.first, p.second) == 0)
            {
                if (min(v[i].first, v[i + 1].first) <= p.first && max(v[i].first, v[i + 1].first) >= p.first && min(v[i].second, v[i + 1].second) <= p.second && max(v[i].second, v[i + 1].second) >= p.second)return true;
            }
        }
        return false;
    }
    bool intersect(ll x1, ll y1, ll x2, ll y2, ll x3, ll y3, ll x4, ll y4) {
        ll a = side(x1, y1, x2, y2, x3, y3);
        ll b = side(x1, y1, x2, y2, x4, y4);
        ll c = side(x3, y3, x4, y4, x1, y1);
        ll d = side(x3, y3, x4, y4, x2, y2);
        if (a != b && c != d)return true;
        else
        {
            if (a != 0 || b != 0 || c != 0 || d != 0)return false;
            else {
                pl p1 = { x1,y1 };
                pl p2 = { x2,y2 };
                pl p3 = { x3,y3 };
                pl p4 = { x4,y4 };
                if (max(p1, p2) >= min(p3, p4) && min(p1, p2) <= max(p3, p4))return true;
                else if (max(p3, p4) >= min(p1, p2) && min(p3, p4) <= max(p1, p2))return false;
                else return false;
            }
        }

    }
    ll solve()
    {
        ifstream F("C:\\Users\\user\\source\\repos\\computer science\\project euler\\Q102 input.txt");
        string s, t;
        ll ans = 0;
        while (getline(F, s))
        {
            vector<ll> v1;
            vector<pl>v;
            ll i = 0; t = "";
            while (i < s.size())
            {
                if (s[i] != ',')t += s[i];
                else
                {
                    v1.push_back(string_to_int(t));
                    t = "";
                }
                i++;
            }
            v1.push_back(string_to_int(t));
            for (ll i = 0; i < 3; i++)v.push_back({ v1[2 * i],v1[2 * i + 1] });
            v.push_back(v[0]);
            ll sum = 0, a = 0, b = 0;
            if (boundry({ a,b }, v))
            {
                ans++;
                continue;
            }
            for (ll i = 0; i < v.size() - 1; i++)
            {
                if (intersect(v[i].first, v[i].second, v[i + 1].first, v[i + 1].second, a, b, a + 1e6, b + 1e6 + 1))sum++;
            }
            if (sum & 1)ans++;
        }
        return ans;
    }
}
namespace Q103// 4,301 microseconds
{
    bool check(vl a)
    {
        vl v(256);
        for (ll i = 1; i < 8; i++)for (ll j = i + 1; j < 8; j++)v[a[i] + a[j]]++;
        for (ll i = 1; i < 8; i++)for (ll j = i + 1; j < 8; j++)for (ll k = j + 1; k < 8; k++)v[a[i] + a[j] + a[k]]++;
        for (ll i : v)if (i > 1)return 0;
        return 1;
    }
    ll solve()
    {
        ll n = 7, k = 255, cnt = 0;//we know there is a solution for 255 (and apperently there isnt a smaller one)
        vector<ll>a(8);
        for (a[3] = 35; a[3] <= (k - 16) / 6; a[3]++)
        {
            for (a[4] = a[3] + 1; a[3] + 5 * a[4] + 14 < k; a[4]++)
            {
                for (a[5] = a[4] + 1; a[3] + a[4] + 4 * a[5] + 10 < k; a[5]++)
                {
                    for (a[6] = a[5] + 1; a[3] + a[4] + a[5] + 3 * a[6] + 7 < k; a[6]++)
                    {
                        for (a[7] = a[6] + 1; a[3] + a[4] + a[5] + a[6] + 2 * a[7] + 5 < k; a[7]++)
                        {
                            for (a[2] = 21; a[2] < a[3]; a[2]++)
                            {
                                for (a[1] = max(ll(7), a[7] + a[6] + a[5] - a[4] - a[3] - a[2] + 1); a[1] <= min(a[2], k - a[7] - a[6]
                                    - a[5] - a[4] - a[3] - a[2]); a[1]++)
                                {
                                    cnt++;
                                    if (check(a))
                                    {
                                        k = 0; for (ll i : a)k += i;
                                        for (ll i = 1; i < 8; i++)cout << a[i];
                                        cout << "\n\n";
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        //cout << cnt << '\n'; its around 18k so its very fast to compute :)
        return 0;
    }
}
namespace Q105
{
    bool check(vl a)
    {
        ll tmp, s = 0, cur;
        for (ll i : a)s += i;
        vl v(s + 1);
        for (ll i = 1, j; i < (1 << a.size()); i++)
        {
            tmp = i; cur = j = 0;
            while (tmp)
            {
                if (tmp & 1)cur += a[j];
                j++; tmp /= 2;
            }
            v[cur]++;
        }
        for (ll i : v)if (i > 1)return 0;
        return 1;
    }
    ll solve()
    {
        ifstream F("C:\\Users\\user\\source\\repos\\computer science\\project euler\\Q105 input.txt");
        ll ans = 0;
        string s, t;
        while (getline(F, s))
        {
            vl v; ll i = 0; t = "";
            while (i < s.size())
            {
                if (s[i] != ',')t += s[i];
                else
                {
                    v.push_back(string_to_int(t));
                    t = "";
                }
                i++;
            }
            v.push_back(string_to_int(t));
            sort(v.begin(), v.end());
            ll a1 = 0, a2 = 0;
            for (ll i = 0; i < (v.size() + 1) / 2; i++)
            {
                a1 += v[i];
            }
            for (ll i = 0; i < (v.size() + 1) / 2 - 1; i++)
            {
                a2 += v[v.size() - i - 1];
            }
            if (a2 >= a1)continue;
            if (!check(v))continue;
            for (ll i : v)ans += i;
        }
        return ans;
    }
}
namespace Q106
{
    ll factorial(ll n)
    {
        if (n == 0)return 1;
        return n * factorial(n - 1);
    }
    ll choose(ll a, ll b)
    {
        return factorial(a) / factorial(b) / factorial(a - b);
    }
    bool check(vl& v1, vl& v2)
    {
        for (ll i = 0; i < v1.size(); i++)
        {
            if (v1[i] > v2[i])return 0;
        }
        return 1;
    }
    ll solve()//1,686 microseconds
    {
        ll n = 12, sum = 0, tmp, cur, cnt = 0;
        for (ll i = 2; i <= n / 2; i++)
        {
            cnt = 0;
            for (ll j = 1; j < 1 << (2 * i - 1); j++)
            {
                vl v1 = { 1 }, v2;
                tmp = j, cur = 2;
                while (cur <= 2 * i)
                {
                    if (tmp & 1)v1.push_back(cur);
                    else v2.push_back(cur);
                    tmp /= 2; cur++;
                }
                if (v1.size() != i)continue;
                sort(v1.begin(), v1.end()); sort(v2.begin(), v2.end());
                if (check(v1, v2))
                {
                    cnt++;
                }
            }
            sum += choose(n, 2 * i) * (choose(2 * i, i) / 2 - cnt);
        }
        return sum;
    }
}
namespace Q108//31,679 microseconds
{
    ll solve()//for n there are (sigma(n^2)+1)/2 answers (num of div)
    {
        ll n = 1e6, k = 1000, tmp, cur, d;
        vl div(n), prime, s(n, 1);
        vector<bool>q(n);
        sieve(prime, div, q, n);
        for (ll i = 2; i < n; i++)
        {
            tmp = i; cur = 0; d = 0;
            while (tmp > 1)
            {
                if (d == div[tmp])
                {
                    cur++;
                }
                else
                {
                    s[i] *= (2 * cur + 1);
                    d = div[tmp]; cur = 1;
                }
                tmp /= div[tmp];
            }
            s[i] *= (2 * cur + 1);
            s[i] = (s[i] + 1) / 2;
        }
        for (ll i = 1; i < n; i++)
        {
            if (s[i] >= k)return i;
        }
        return 0;
    }
}
namespace Q112//16,912 microseconds
{
    bool bouncy(ll n)
    {
        bool d = 1, i = 1;
        ll cur = -1;
        while (n)
        {
            if (cur == -1)
            {
                cur = n % 10;
                continue;
            }
            if (cur < n % 10)i = 0;//or opposite
            if (cur > n % 10)d = 0;//or opposite
            cur = n % 10;
            n /= 10;
        }
        return i || d;
    }
    ll solve()
    {
        ll n = 99, cnt = 0;
        for (ll i = 1;; i++)
        {
            if (!bouncy(i))cnt++;
            if (cnt * 100 == n * i)return i;
        }
    }
}
namespace Q113//286 microseconds
{
    ll choose(ll a, ll b)
    {
        ll ans = 1;
        for (ll i = 0; i < b; i++)ans *= (a - i);
        for (ll i = 2; i <= b; i++)ans /= i;
        return ans;
    }
    ll solve()//maf
    {
        ll n = 100, ans = 0;
        for (ll k = 1; k <= n; k++)
        {
            ans += choose(k + 9, 9) + choose(k + 8, 8) - 10;
        }
        return ans;
    }
}
namespace Q114//436 microseconds
{
    ll solve()
    {
        ll n = 50;
        vl dp(n + 1); dp[0] = dp[1] = dp[2] = 1;
        for (ll i = 3; i <= n; i++)
        {
            dp[i] = 1 + dp[i - 1];//whole row and first one is gray
            for (ll j = 3; j < i; j++)
            {
                dp[i] += dp[i - j - 1];
            }
        }
        return dp[n];
    }
}
namespace Q115//143 microseconds
{
    ll solve()
    {
        ll m = 50, n = 1e6;
        vl dp(m, 1);
        for (ll i = m;; i++)
        {
            dp.push_back(1 + dp[i - 1]);//whole row and first one is gray
            for (ll j = m; j < i; j++)
            {
                dp[i] += dp[i - j - 1];
            }
            if (dp[i] > n)return i;
        }
        return -1;
    }
}
namespace Q116//385 microseconds
{
    ll solve()
    {
        ll n = 50, sum = 0;
        vl v(n + 1); v[0] = v[1] = 1;
        for (ll i = 2; i <= n; i++)v[i] = v[i - 1] + v[i - 2];
        sum += v[n];
        v[0] = v[1] = v[2] = 1;
        for (ll i = 3; i <= n; i++)v[i] = v[i - 1] + v[i - 3];
        sum += v[n];
        v[0] = v[1] = v[2] = v[3] = 1;
        for (ll i = 4; i <= n; i++)v[i] = v[i - 1] + v[i - 4];
        sum += v[n];
        return sum - 3;//all gray happens 3 times
    }
}
namespace Q117//247 microseconds
{
    ll solve()
    {
        ll n = 50, sum = 0;
        vl v(n + 1); v[0] = v[1] = 1; v[2] = 2; v[3] = 4;
        for (ll i = 4; i <= n; i++)v[i] = v[i - 1] + v[i - 2] + v[i - 3] + v[i - 4];
        return sum = v[n];
    }
}
namespace Q119//880759 microseconds
{
    ll ds(ll n)
    {
        ll sum = 0; n *= 10;
        while (n /= 10)sum += n % 10;
        return sum;
    }
    ll solve()
    {
        ll n = 1e15, m = sqrt(n), tmp;
        set<ll>s;
        for (ll i = 2; i <= m; i++)
        {
            tmp = i;
            while (1)
            {
                if (ds(tmp) == i && tmp >= 10)s.insert(tmp);
                if (tmp > n / i)break;
                tmp *= i;
            }
        }
        auto ans = s.begin();
        for (ll i = 0; i < 29; i++)ans++;
        return *ans;
    }
}
namespace Q120//293 microseconds
{
    ll solve()//the value is max 2*a*n mod a^2 for odd n (notice there is period 2*a when a odd for only odd n)
    {
        ll n = 1e3, sum = 0, m = 0;
        for (ll a = 3; a <= n; a++)
        {
            if (a % 2 == 1)m = a * a - a;
            if (a % 2 == 0)m = a * a - 2 * a;
            sum += m;
        }
        return sum;
    }
}
namespace Q123
{
    ll solve()//ive got no fucking idea why the answer is off by 1 so i wont submit it for now
    {
        ll n = 1e6, m = 1e10, p;
        vl prime; nsieve(prime, n);
        for (ll i = 1;; i++)
        {
            p = prime[i];
            if (2 * p * (i + 1) > m)
            {
                cout << p << '\n';
                return i + 1;
            }
        }
    }
}
namespace Q124//7,022 microseconds
{
    ll solve()
    {
        ll n = 1e5, m = 1e4;
        vl div(n + 1), prime;
        vector<bool>p(n + 1);
        sieve(prime, div, p, n + 1);
        vector<pl> rad(n + 1);
        rad[0] = rad[1] = { 1,1 };
        for (ll i = 2; i <= n; i++)
        {
            if (p[i])rad[i] = { i,i };
            else
            {
                if (i / div[i] % div[i])rad[i] = { div[i] * rad[i / div[i]].first,i };
                else rad[i] = { rad[i / div[i]].first,i };
            }
        }
        sort(rad.begin(), rad.end());
        return rad[m].second;
    }
}
namespace Q125
{
    ll solve()
    {
        ll n = 1e8, sum = 0;
        ll k = sqrt(n);
        vector<ll>v(k + 1), ans;
        for (ll i = 1; i <= k; i++)v[i] = v[i - 1] + i * i;
        for (ll i = 2; i <= k; i++)for (ll j = i - 2; j >= 0; j--)
        {
            if (v[i] - v[j] > n)break;
            if (palindrome(v[i] - v[j], 10))ans.push_back(v[i] - v[j]);
        }
        sort(ans.begin(), ans.end());
        for (ll i = 0; i < ans.size(); i++)
        {
            if (!i || ans[i] != ans[i - 1])sum += ans[i];
        }
        return sum;
    }
}
namespace Q131//9,933 microseconds
{
    ll solve()//need 3x^2+3x+1 to be prime
    {
        ll n = 1e6, cnt = 0;
        vl div(n), prime;
        vector<bool>p(n);
        sieve(prime, div, p, n);
        for (ll i = 0; 3 * i * (i + 1) + 1 <= n; i++)
        {
            if (p[3 * i * (i + 1) + 1])cnt++;
        }
        return cnt;
    }
}
namespace Q133
{
    ll solve()
    {
        ll n = 1e5, sum = 0, k, cur;
        vector <ll> prime, div(n);
        vector<bool>q(n);
        sieve(prime, div, q, n);
        for (ll p : prime)
        {
            if (p == 2 || p == 5)continue;
            k = 1;
            for (ll i = 1; i <= p; i++)
            {
                if (k % p)
                {
                    k = 10 * k + 1; k %= p;
                }
                else
                {
                    cur = i; break;
                }
            }
            while (cur > 1)
            {
                if (div[cur] != 2 && div[cur] != 5)
                {
                    sum += p; break;
                }
                cur /= div[cur];
            }
        }
        return sum + 7;//7 for 2,5
    }
}
namespace Q134//66,495 microseconds
{
    ll solve()//basically need inverse of 10^k
    {
        ll n = 2e6, sum = 0, p1, p2, k;
        vl prime; nsieve(prime, n);
        for (ll i = 2; prime[i] <= n / 2; i++)
        {
            p1 = prime[i], p2 = prime[i + 1];
            k = log(p1);
            sum += p1;
            sum += power(10, k) * (power(power(10, k), p2 - 2, p2) * (p2 - p1) % p2);
        }
        return sum;
    }
}
namespace Q135//5,645,778 microseconds
{
    ll solve()
    {
        ll n = 1e6, cnt = 0, cur, k = 10;
        for (ll i = 1; i < n; i++)
        {
            if (i % 4 == 1 || i % 4 == 2)continue;
            cur = 0;
            vl d = divisors::find_divisors(i);
            for (ll x : d)
            {
                if ((i + x * x) % (4 * x) == 0)
                {
                    if (3 * x * x > i)cur++;
                }
            }
            if (cur == k)
            {
                cnt++;
                cout << i << '\n';
            }
        }
        return cnt;
    }
}
namespace Q136//451,883 microseconds
{
    ll solve()
    {
        ll n = 5e7, cnt = 0, cur;
        vl div(n), prime;
        vector<bool>p(n);
        sieve(prime, div, p, n);
        for (ll i = 1; i < 1e3; i++)//small cases cant bother to check
        {
            if (i % 4 == 1 || i % 4 == 2)continue;
            cur = 0;
            vl d = divisors::find_divisors(i);
            for (ll x : d)
            {
                if ((i + x * x) % (4 * x) == 0)
                {
                    if (3 * x * x > i)cur++;
                }
            }
            if (cur == 1)
            {
                cnt++;
            }
        }
        for (ll x : prime)//we find that either the number is 4 times odd prime or 16 times prime, or prime that is 3 mod 4 
        {
            if (4 * x >= 1e3 && 4 * x < n)cnt++;
            if (16 * x >= 1e3 && 16 * x < n)cnt++;
            if (x % 4 == 3 && x >= 1e3 && x < n)cnt++;
        }
        return cnt;
    }
}
namespace Q137//2,301 microseconds
{
    ll solve()
    {
        ll m, k = 15, cnt = 0;
        vl ans;
        for (m = 2; cnt < k; m++)
        {
            ll d = sqrt(5 * m * m - 4);
            if (d * d - 5 * m * m + 4)continue;
            ans.push_back((-m + d) / 2 * m);
            cnt++;
        }
        return ans[k - 1];
    }
}
ll Q142()
{
    vector<ll> sq;
    for (ll i = 0; i < 1e6; i++)
    {
        sq.push_back(i * i);
    }
    ll r = 1e7, a, b, c, d, e, f, mn = 1e15;
    for (ll i = 2; i * i * i * i < r; i++) {
        for (ll j = 1; j < i; j++) {
            for (ll mlt = 1; mlt * i * i * i * i < r; mlt++) {
                a = mlt * (i * i + j * j), c = mlt * max(2 * i * j, i * i - j * j), f = mlt * min(2 * i * j, i * i - j * j);
                for (ll k = 1; k * k * k * k < r; k++) {
                    ll m = Find(sq, a - k * k);
                    if (m < k && m>0) {
                        e = 2 * m * k, d = k * k - m * m;
                        b = Find(sq, c * c - e * e);
                        if (b > 0) {
                            if ((~(a ^ b) & 1) && (~(c ^ d) & 1) && (a > 0 && b > 0 && c > 0 && d > 0 && e > 0 && f > 0)) {
                                mn = min(mn, a * a + c * c + e * e);
                                if (a * a + c * c + e * e == 22252)
                                    cout << a << ' ' << c << ' ' << e << ' ' << b << ' ' << d << ' ' << f << endl;
                            }
                        }
                        e = d ^ e; d = d ^ e; e = d ^ e;
                        b = Find(sq, c * c - e * e);
                        if (b > 0) {
                            if ((~(a ^ b) & 1) && (~(c ^ d) & 1) && (a > 0 && b > 0 && c > 0 && d > 0 && e > 0 && f > 0)) {
                                mn = min(mn, a * a + c * c + e * e);
                                if (a * a + c * c + e * e == 22252)
                                    cout << a << ' ' << c << ' ' << e << ' ' << b << ' ' << d << ' ' << f << endl;
                            }
                        }
                    }
                }
            }
        }
    }
    return mn / 2;
}
namespace Q145//251,687,492 microseconds ~ 4 minutes
{
    ll solve()//just brute force
    {
        ll n = 1e9, cnt = 0, b, m, tmp, flag;
        for (ll i = 1; i < n; i++)
        {
            flag = 1;
            if (i % 10 == 0)continue;
            if (i % 10000000 == 1)cout << i << '\n';
            tmp = b = i, m = 1;
            vl v;
            while (tmp)
            {
                v.push_back(tmp % 10); tmp /= 10;
            }
            for (ll i = v.size() - 1; ~i; i--)
            {
                b += m * v[i];
                m *= 10;
            }
            while (b)
            {
                if (b % 2 == 0)
                {
                    flag = 0; break;
                }
                b /= 10;
            }
            cnt += flag;
        }
        return cnt;
    }
}
namespace Q146//60,278,632 microseconds
{
    ll solve()// can probably automate the case work and do faster time
    {
        ll n = 15e7, sum = 0, a, tmp;
        for (ll i = 10; i <= n; i += 2)
        {
            bool flag = 1;
            if (i % 3 == 0 || i % 13 == 0 || i % 5)continue;
            if (i % 7 != 3 && i % 7 != 4)continue;
            tmp = i % 11;
            if (tmp == 2 || tmp == 3 || tmp == 8 || tmp == 9)continue;
            tmp = i % 13;
            if (tmp != 1 && tmp != 3 && tmp != 4 && tmp != 9 && tmp != 10 && tmp != 12)continue;
            a = i * i;
            if (!is_prime(a + 1) || !is_prime(a + 3) || !is_prime(a + 7) || !is_prime(a + 9) || !is_prime(a + 13)
                || !is_prime(a + 27) || is_prime(a + 21))
            {
                flag = 0;
            }
            if (flag)sum += i;
        }
        return sum;
    }
}
namespace Q149
{
    ll solve()
    {
        ll n = 2000, p = 1000000, m = 0, cur;
        vvl v(n, vl(n));
        for (ll i = 0; i < n; i++)for (ll j = 0; j < n; j++)
        {
            if (i == 0 && j < 55)v[i][j] = (100003 - 200003 * (j + 1) + 300007 * (j + 1) * (j + 1) * (j + 1)) % p - p / 2;
            else
            {
                v[i][j] = (v[(n * i + j - 24) / n][(n * i + j - 24) % n] +
                    v[(n * i + j - 55) / n][(n * i + j - 55) % n] + p) % p - p / 2;
            }
        }
        for (ll i = 0; i < n; i++)
        {
            cur = 0;
            for (ll j = 0; j < n; j++)
            {
                cur += v[i][j];
                m = max(m, cur);
                if (cur < 0)cur = 0;
            }
        }
        for (ll j = 0; j < n; j++)
        {
            cur = 0;
            for (ll i = 0; i < n; i++)
            {
                cur += v[i][j];
                m = max(m, cur);
                if (cur < 0)cur = 0;
            }
        }
        for (ll d = 0; d < 2 * n - 1; d++)
        {
            cur = 0;
            for (ll i = max(ll(0), d - n + 1); i <= min(d, n - 1); i++)
            {
                cur += v[i][d - i];
                m = max(m, cur);
                if (cur < 0)cur = 0;
            }
        }
        for (ll d = -n + 1; d <= n - 1; d++)
        {
            cur = 0;
            for (ll i = max(ll(0), -d); i <= min(n - 1, n - d - 1); i++)
            {
                cur += v[i + d][i];
                m = max(m, cur);
                if (cur < 0)cur = 0;
            }
        }
        return m;
    }
}
namespace Q150//5,522,644 microseconds
{
    ll solve()
    {
        ll t = 0, n = 1e3, m = 1e8;
        vector<vvl>v(n + 1);
        for (ll i = 0; i <= n; i++)
        {
            v[i] = vvl(i + 1, vl(n - i + 1));
        }
        for (ll k = 0; k < n; k++)
        {
            for (ll j = 0; j <= k; j++)
            {
                t = 615949 * t + 797807; t %= (1 << 20);
                v[k][j][1] = t - (1 << 19);
            }
        }
        for (ll i = 2; i <= n; i++)
        {
            for (ll k = 0; k <= n - i; k++)
            {
                for (ll j = 0; j <= k; j++)
                {
                    v[k][j][i] = v[k + 1][j][i - 1] + v[k + 1][j + 1][i - 1] - v[k + 2][j + 1][i - 2] + v[k][j][1];
                }
            }
        }
        for (ll i = 0; i <= n; i++)
        {
            for (ll k = 0; k <= n - i; k++)
            {
                for (ll j = 0; j <= k; j++)
                {
                    m = min(m, v[k][j][i]);
                }
            }
        }
        return m;
    }
}
ll Q152()
{
    ll p = 1e9 + 7;
    ll sum = 0;
    vector<long double> v = { 2,3,4,5,6,7,8,9,10,12,14,15,16,18,20,21,24,25,27,28,30,32,35,36,40,42,45,48,49,50,54,56,60,63,64,70,72,75,80,12 };
    vector<long double> ans(40);//v.size()=40
    for (ll i = 0; i < ans.size(); i++)
    {
        ans[i] = 1 / (v[i] * v[i]);
    }
    vector<long double>ans1(power(2, 20)), ans2(power(2, 20));
    for (ll i = 0; i < (1 << 20); i++)
    {
        ll tmp = i, cnt = 0;
        while (tmp)
        {
            if (tmp & 1)
            {
                ans1[i] += ans[cnt];
                ans2[i] += ans[cnt + 20];
            }
            tmp /= 2; cnt++;
        }
    }
    sort(ans2.begin(), ans2.end());
    for (long double x : ans1)
    {
        sum += upper_bound(ans2.begin(), ans2.end(), 0.5 - x + 1 / 1e15) - lower_bound(ans2.begin(), ans2.end(), 0.5 - x - 1 / 1e15);//need 1e15
    }
    return sum;
}
namespace Q159//1,592,252 microseconds
{
    ll solve()
    {
        ll n = 1e6, sum = 0;
        vl dp(n + 1); dp[1] = 0;
        for (ll i = 2; i < n; i++)
        {
            dp[i] = (i + 8) % 9 + 1;
            if (i == 24)
            {
                ll c = 2;
            }
            for (ll j = 2;; j++)
            {
                if (j * j > i)break;
                if (i % j)continue;
                dp[i] = max((j + 8) % 9 + 1 + dp[i / j], dp[i]);
                dp[i] = max(dp[i], (i / j + 8) % 9 + 1 + dp[j]);
            }
            sum += dp[i];
        }
        return sum;
    }
}
namespace Q164//943 microseconds
{
    ll solve()
    {
        ll n = 20, cnt = 0;
        vector<vvl>dp(n + 1, vvl(10, vl(10))); dp[0][0][0] = 1;//keeps how many numbers with ith,i+1th digits are a,b
        for (ll i = 0; i < n; i++)
        {
            for (ll a = 0; a < 10; a++)for (ll b = 0; b < 10; b++)
            {
                if (!dp[i][a][b])continue;
                for (ll c = 0; c <= 9 - a - b; c++)
                {
                    if (i == 0 && c == 0)continue;
                    dp[i + 1][b][c] += dp[i][a][b];
                }
            }
        }
        for (ll a = 0; a < 10; a++)for (ll b = 0; b < 10; b++)
        {
            cnt += dp[n][a][b];
        }
        return cnt;
    }
}
namespace Q166//435,108 microseconds
{
    int solve()
    {
        int x, a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, cnt = 0;
        for (a = 0; a < 10; a++)for (b = 0; b < 10; b++)for (c = 0; c < 10; c++)for (d = 0; d < 10; d++)
        {
            x = a + b + c + d;
            for (e = 0; e <= min(9, x - a); e++)for (f = 0; f <= min(9, x - b); f++)
                for (g = 0; g <= min(9, x - e - f); g++)
                {
                    h = x - e - f - g;
                    if (h > 9)continue;
                    for (i = 0; i <= min(9, x - a - e); i++)
                    {
                        j = a + e + i - g - d;
                        if (j < 0 || j > 9)continue;
                        k = x - i - j + d + h - a - f;
                        l = x - i - j - d - h + a + f;
                        if (k < 0 || l < 0)continue;
                        if (k & 1 || l & 1)continue;
                        k /= 2; l /= 2;
                        if (k > 9 || l > 9)continue;
                        m = x - a - e - i, n = x - b - f - j, o = x - c - g - k, p = x - d - h - l;
                        if (m > 9 || m < 0)continue;
                        if (n > 9 || n < 0)continue;
                        if (o > 9 || o < 0)continue;
                        if (p > 9 || p < 0)continue;
                        cnt++;
                    }
                }
        }
        return cnt;
    }
}
namespace Q169//254 microseconds
{
    ll solve()//classic dp
    {
        ll n = power(5, 25), tmp = n;
        vl v(25);
        while (tmp)
        {
            v.push_back(tmp & 1); tmp >>= 1;
        }
        vvl dp(v.size(), vl(2));
        dp[0][0] = 1;
        for (ll i = 0; i < v.size() - 1; i++)
        {
            if (v[i] == 0)
            {
                dp[i + 1][0] = dp[i][0];
                dp[i + 1][1] = dp[i][1] + dp[i][0];
            }
            else
            {
                dp[i + 1][0] = dp[i][0] + dp[i][1];
                dp[i + 1][1] = dp[i][1];
            }
        }
        return dp[v.size() - 1][1] + dp[v.size() - 1][0];
    }
}
namespace Q173//1,336 microseconds
{
    ll solve()
    {
        ll n = 1e6, sum = 0;
        for (ll a = n / 4 + 1; a > 1; a--)
        {
            for (ll b = a - 2; b > 0; b -= 2)
            {
                if (a * a - b * b > n)break;
                sum++;
            }
        }
        return sum;
    }
}
namespace Q174//6,434 microseconds
{
    ll solve()
    {
        ll n = 1e6, sum = 0;
        vl v(n + 1);
        for (ll a = n / 4 + 1; a > 1; a--)
        {
            for (ll b = a - 2; b > 0; b -= 2)
            {
                if (a * a - b * b > n)break;
                v[a * a - b * b]++;
            }
        }
        for (ll x : v)if (x <= 10 && x)sum += 1;
        return sum;
    }
}
namespace Q179//216,344 microseconds
{
    ll solve()
    {
        ll n = 1e7, cnt = 0, tmp, cur;
        vl div(n + 1), prime;
        vector<bool>p(n + 1);
        sieve(prime, div, p, n + 1);
        vl d(n + 1); d[1] = 1;
        for (ll i = 2; i < n; i++)
        {
            tmp = i, cur = 1;
            while (tmp % div[i] == 0)
            {
                tmp /= div[i]; cur++;
            }
            d[i] = d[tmp] * cur;
        }
        for (ll i = 2; i < n; i++)if (d[i] == d[i + 1])cnt++;
        return cnt;
    }
}
namespace Q183//976 microseconds
{
    ll solve()
    {
        double n = 1e4, m1, m2;
        ll sum = 0, k;
        for (double i = 5; i <= n; i++)
        {
            k = i / exp(1);
            m1 = k * (log(i) - std::log(k));
            m2 = (k + 1) * (log(i) - std::log(k + 1));
            if (m1 > m2)
            {
                k /= gcd(i, k);
                while (k % 2 == 0)k /= 2;
                while (k % 5 == 0)k /= 5;
                if (k != 1)sum += i;
                else sum -= i;
            }
            else
            {
                k++;
                k /= gcd(i, k);
                while (k % 2 == 0)k /= 2;
                while (k % 5 == 0)k /= 5;
                if (k != 1)sum += i;
                else sum -= i;
            }
        }
        return sum;
    }
}
namespace Q187//1,372,193 microseconds
{
    ll solve()
    {
        ll n = 1e8, cnt = 0;
        vl div(n), prime;
        vector<bool>q(n);
        sieve(prime, div, q, n);
        for (ll i = 4; i < n; i++)
        {
            if (q[i])continue;
            if (q[i / div[i]])cnt++;
        }
        return cnt;
    }
}
ll Q188(ll a, ll b, ll MAX)
{
    if (MAX == 1)return 0;
    ll tmp = MAX;
    if (MAX % 5 == 0)tmp = (MAX / 5) * 2;
    else tmp /= 2;
    return power(a, Q188(a, b - 1, tmp), MAX);
}
namespace Q190
{
    //in python
}
ll Q193()
{
    ll n = power(2, 25), sz, cnt = 0, mult = 1, k = n * n;
    vector <ll> prime;
    vector<ll>div(n);
    vector<bool>p(n);
    sieve(prime, div, p, n);
    sz = prime.size();
    for (ll i = 0; i < sz; i++)
    {
        cnt += (n * n) / (prime[i] * prime[i]);
    }
    for (ll i = 2; i < 9; i++)
    {
        vl v(i);
        for (ll j = 0; j < i; j++)
        {
            v[j] = j;
        }
        while (1)
        {
            if (i == 9)
            {
                ll c = 2;
            }
            mult = 1;
            for (ll j = 0; j < i; j++)
            {
                mult *= prime[v[j]] * prime[v[j]];
            }
            if (i == 9)cout << mult - k << "\n";
            if (mult <= n * n)
            {
                cnt += (2 * (i & 1) - 1) * ((n * n) / mult);
                v[i - 1]++;
            }
            else
            {
                ll tmp = i - 2;
                while (tmp >= 0 && v[tmp] + 1 == v[tmp + 1])tmp--;
                if (tmp == -1)break;
                else
                {
                    v[tmp]++;
                    for (ll r = tmp + 1; r < i; r++)
                    {
                        v[r] = v[r - 1] + 1;
                    }
                }
            }
        }
    }
    return n * n - cnt;
}
ll Q202()
{
    ll n = 6008819575, cnt = 0;
    for (ll i = 2; i < n; i += 3)
    {
        if (gcd(i, n) == 1)cnt++;
    }
    return cnt;
}
ll Q203()
{
    ll n = 100, tmp, sum = 0;
    vector <ll> prime;
    vector<ll>div(n);
    vector<bool>p(n);
    sieve(prime, div, p, n);
    vector<pl>v;
    vvl bin(52, vl(52));
    bin[1][0] = bin[1][1] = 1;
    v.push_back({ 1,1 });
    for (ll i = 2; i <= 51; i++)
    {
        bin[i][0] = bin[i][i] = 1;
        for (ll j = 1; j < i; j++)
        {
            bin[i][j] = bin[i - 1][j - 1] + bin[i - 1][j];
        }
    }
    for (ll i = 2; i <= 50; i++)
    {
        for (ll j = 1; j < i; j++)
        {
            ll flag = 1;
            for (ll k = 2; k <= i; k++)
            {
                ll cnt = 0, tmp = k;
                if (!p[k])continue;
                while (tmp <= i)
                {
                    cnt += (i / tmp - (i - j) / tmp - j / tmp);
                    tmp *= k;
                }
                if (cnt > 1)
                {
                    flag = 0;
                    break;
                }
            }
            if (flag)v.push_back({ i,j });
        }
    }
    set<ll> ans;
    for (pl x : v)
    {
        ans.insert(bin[x.first][x.second]);
    }
    for (auto& num : ans)sum += num;
    return sum;
}
ll Q204()
{
    ll cnt = 0, x, least;
    vector<ll>p = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97 };
    queue<ll>q;
    for (ll i : p)q.push(i);
    while (!q.empty())
    {
        x = q.front(); q.pop();
        cnt++;
        for (int i = 0; i < p.size(); i++) {
            if (x % p[i] == 0)
            {
                for (int j = i; j >= 0; j--)
                {
                    if (p[j] * x <= 1e9)
                    {
                        q.push(p[j] * x);
                    }
                }
                break;
            }
        }
    }
    return cnt + 1;//+1 for the number 1
}
namespace Q205
{
    double solve()
    {
        ll a = 6, b = 9, c = 0;
        vvl va(7, vl(37)), vb(10, vl(37)); va[0][0] = vb[0][0] = 1;
        for (ll i = 1; i <= a; i++)
        {
            for (ll j = 0; j <= 36; j++)
            {
                for (ll k = 1; k <= 6; k++)
                {
                    if (j + k > 36)break;
                    va[i][j + k] += va[i - 1][j];
                }
            }
        }
        for (ll i = 1; i <= 9; i++)
        {
            for (ll j = 0; j <= 36; j++)
            {
                for (ll k = 1; k <= 4; k++)
                {
                    if (j + k > 36)break;
                    vb[i][j + k] += vb[i - 1][j];
                }
            }
        }
        for (ll i = 0; i <= 36; i++)for (ll j = 0; j < i; j++)
        {
            c += va[6][j] * vb[9][i];
        }
        return double(c) / (power(6, 6) * power(4, 9));
    }
}
ll Q206()
{
    bool odd = 0, flag;
    ll x;
    for (ll i = 1e8 + 3; i < 15 * 1e7; odd ? i += 4 : i += 6)
    {
        flag = 1;
        odd = 1 - odd;
        x = i * i;
        for (int j = 0; j < 9; j++)
        {
            if (x % 10 != 9 - j)
            {
                flag = 0;
                break;
            }
            x /= 100;
        }
        if (flag)return 10 * i;
    }
    return 0;
}
ll Q209()
{
    vector<ll>v(64);
    vector<bool>seen(64);
    vector<ll>fib(64); fib[0] = 0; fib[1] = 1;
    for (ll i = 2; i < 64; i++)
    {
        fib[i] = fib[i - 1] + fib[i - 2];
    }
    ll a, b, c, j, cnt = 0, mlt = 1, tmp;
    for (ll i = 0; i < 64; i++)
    {
        a = i / 32; b = (i / 16) % 2; c = (i / 8) % 2;
        v[i] = (i % 32) * 2 + (a ^ (b & c));
    }
    for (ll i = 0; i < 64; i++)
    {
        cnt = 0;
        j = i;
        if (seen[j])continue;
        while (!seen[j])
        {
            seen[j] = 1;
            j = v[j];
            cnt++;
        }
        if (cnt > 1)
        {
            ll sum = fib[cnt - 1] + 2;
            for (ll i = 2; i < cnt; i++)
            {
                sum += fib[cnt - i + 1];
            }
            mlt *= sum;
        }
    }
    return mlt;
}
namespace Q211//2,481,742 microseconds
{
    ll solve()
    {
        ll n = 64e6, sum = 0, d, tmp, cur;
        vl div(n), prime, sd(n);
        vector<bool>p(n);
        sieve(prime, div, p, n);
        sd[1] = 1;
        for (ll i = 2; i < n; i++)
        {
            if (p[i])
            {
                sd[i] = i * i + 1; continue;
            }
            tmp = i, cur = 1;
            while (tmp % div[i] == 0)
            {
                tmp /= div[i]; cur++;
            }
            sd[i] = sd[tmp] * ((power(div[i], 2 * cur) - 1) / (div[i] * div[i] - 1));
        }
        for (ll i = 1; i < n; i++)
        {
            d = sqrt(sd[i]);
            if (d * d == sd[i])sum += i;
        }
        return sum;
    }
}
namespace Q214//1,135,533 microseconds
{
    ll solve()
    {
        ll n = 4e7, sum = 0;
        vl div(n), prime, phi(n), len(n); len[1] = phi[1] = 1;
        vector<bool>p(n);
        sieve(prime, div, p, n);
        for (ll i = 2; i < n; i++)
        {
            if (p[i])phi[i] = i - 1;
            else
            {
                if (i / div[i] % div[i] == 0)phi[i] = div[i] * phi[i / div[i]];
                else phi[i] = (div[i] - 1) * phi[i / div[i]];
            }
            len[i] = len[phi[i]] + 1;
        }
        for (ll x : prime)if (len[x] == 25)sum += x;
        return sum;
    }
}
namespace Q218//962,088 microseconds
{
    ll solve()
    {
        return 0;
    }
    ll solve1()// yes the answer is actually 0 :< disvisivility by 4,6 trivial so only 7 which comes from casework (r=1,2,3,k=1,2,3 r!=k)
    {
        ll n = 1e8, m = sqrt(n), cnt = 0;
        for (ll r = 1; r <= m; r++)
        {
            if (r % 200 == 0)cout << r << '\n';
            for (ll k = 1; k < r; k++)
            {
                if (r * r + k * k > n)break;
                if (r & 1 && k & 1)continue;
                if (r % 7 == 0 || k % 7 == 0)continue;
                if (gcd(r, k) - 1)continue;
                ll a = r * r - k * k, b = 2 * r * k;
                if (a % 7 == 0)continue;
                if ((a * a - b * b) % 7 == 0)continue;
                cnt++;
                cout << r << "  " << k << '\n';
            }
        }
        return cnt;
    }
}
namespace Q221//7,205,724 microseconds
{
    ll solve()
    {
        ll cnt = 0, k = 150000, n, q = 2, m = 1e18;
        vl ans;
        for (; q * q * q < 2 * m; q++)
        {
            if (ans.size() >= 3 * k)
            {
                sort(ans.begin(), ans.end());
                ans.resize(k);
                m = ans[k - 1];
            }
            n = q * q + 1;
            vl div = divisors::find_divisors(n);
            if (is_prime(n))
            {
                if (q * (q - 1) > m / (n - q))continue;
                ans.push_back(q * (q - 1) * (n - q));
                continue;
            }
            else for (ll d : div)
            {
                if (d > q / 2 || (q - d) * q > m / (n / d - q))continue;
                cnt++;
                ans.push_back((q - d) * q * (n / d - q));
            }
        }
        sort(ans.begin(), ans.end());
        return ans[149999];
    }
}
namespace Q231//1,198,916 microseconds
{
    ll vp(ll n, ll m, ll p)//vp of n C m
    {
        ll sum = 0, tmp = p;
        while (tmp <= n)
        {
            sum += n / tmp - m / tmp - (n - m) / tmp;
            tmp *= p;
        }
        return sum;
    }
    ll solve()
    {
        ll n = 2e7, k = 15e6, sum = 0;
        vl div(n + 1), prime;
        vector<bool>p(n + 1);
        sieve(prime, div, p, n + 1);
        for (ll x : prime)
        {
            sum += vp(n, k, x) * x;
        }
        return sum;
    }
}
long double Q232()
{
    ll a = 128;
    vector<vector<long double>> dp(101 + a, vector<long double>(101 + a, -1));
    for (ll i = 0; i <= 100 + a; i++)
    {
        for (ll j = 0; j <= 100 + a; j++)
        {
            if (j <= a)
            {
                dp[i][j] = 1;
                continue;
            }
            if (i <= a)
            {
                dp[i][j] = 0;
                continue;
            }
            long double m = 0;
            for (ll t = 1; t < 9; t++)
            {
                long double pw = pow(2, -t);
                m = max(m, (dp[i - 1][j] * (1 - pw) + dp[i - 1][j - (1 << (t - 1))] * pw + dp[i][j - (1 << (t - 1))] * pw) / (1 + pw));
            }
            dp[i][j] = m;
        }
    }
    return (0.5 * dp[100 + a][100 + a] + 0.5 * dp[99 + a][100 + a]);
}
namespace Q234//8,368 microseconds
{
    ll s(ll a, ll b, ll p)
    {
        if (a > b)return 0;
        ll start = (a - 1) + p - (a - 1) % p;
        ll end = b - b % p;
        ll n = (end - start) / p + 1;
        return n * (start + end) / 2;
    }
    ll solve()
    {
        ll n = 999966663333, k = 2 * sqrt(n), sum = 0, p, q;
        vl prime; nsieve(prime, k);
        for (ll i = 0;; i++)
        {
            p = prime[i], q = prime[i + 1];
            if (q * q > n)break;
            sum += s(p * p + 1, q * q - 1, p);
            sum += s(p * p + 1, q * q - 1, q);
            sum -= 2 * s(p * p + 1, q * q - 1, p * q);
        }
        sum += s(p * p + 1, n, p);
        sum += s(p * p + 1, n, q);
        sum -= 2 * s(p * p + 1, n, p * q);
        return sum;
    }
}
namespace Q237//3,745,647,533 microseconds ~ 1h 2min
{
    ll solve()//After a lot of work we get the following code (its too long to explain how but trust)
    {
        ll n = 1e12;
        int p = 1e8;
        int a0 = 1, a1 = 2, a2 = 6, a3 = 14, tmp;
        vl A = { 1,2,6,14 };
        for (ll i = 4; i <= n; i++)
        {
            if (i % 1000000000 == 0)cout << i / 1000000000 << '\n';
            tmp = ((p + a3 + a2 - a1) << 1) + a0; tmp %= p;
            a0 = a1; a1 = a2; a2 = a3; a3 = tmp;
        }
        return (p + a2 - a1) % p;
    }
}
namespace Q250//17,310,231 microseconds
{
    ll solve()
    {
        ll n = 250250, k = 250, sum = 0, p = 1e16, a, b, m = -1;
        vl v(k), ans(k), tmp(k); ans[0] = 1;
        for (ll i = 1; i <= n; i++)
        {
            v[power(i, i, k)]++;
        }
        for (ll x : v)m = max(m, x);
        vvl ch(m + 1);
        for (ll i = 0; i <= m; i++)
        {
            ch[i].resize(i + 1, 1);
            for (ll j = 1; j < i; j++)ch[i][j] = (ch[i - 1][j - 1] + ch[i - 1][j]) % p;
        }
        for (ll i = 0; i < k; i++)
        {
            cout << i << '\n';
            tmp = vl(k, 0);
            for (ll j = 0; j <= v[i]; j++)
            {
                b = ch[v[i]][j];
                a = (i * j) % k;
                for (ll r = 0; r < k; r++)
                {
                    tmp[(a + r) % k] += power(ans[r], b, p, [](ll a, ll b, ll c) {return (a + b) % c; }, 0);
                    tmp[(a + r) % k] %= p;
                }
            }
            ans = tmp;
        }
        return ans[0] - 1;
    }
}
namespace Q254//1,743,764 microseconds
{
    ll sm(ll n) {
        ll s = 0;
        while (n)s += n % 10, n /= 10;
        return s;
    }
    ll solve(ll n = 150) {
        vl fct = { 1,1,2,6,24,120,720,5040,40320,362880 };
        ll s = 0;
        for (ll i = 1; i <= n; i++) {
            ll x = i % 9, j = i;
            while (j > 8) {
                x *= 10;
                x += 9;
                j -= 9;
            }
            ll c = x / fct[9], mn = 1e18, l2 = 0, a3 = x % fct[9];
            for (int b = 8; b; b--) {
                l2 += a3 / fct[b];
                a3 %= fct[b];
            }
            vl mnv(10, 100);
            for (ll k = c; k <= c + l2; k++) {
                ll y = k * fct[9];
                for (int a2 = i % 9; a2 < fct[9]; a2 += 9) {
                    ll a = a2;
                    if (sm(a + y) == i) {
                        vl v(10);
                        ll l = v[9] = k;
                        for (int b = 8; b; b--) {
                            l += v[b] = a / fct[b];
                            a %= fct[b];
                        }
                        if (mn > l) {
                            mn = l, mnv = v;
                        }
                        if (mn == l) {
                            for (int b = 1; b < 10; b++) {
                                if (mnv[b] > v[b])break;
                                else if (mnv[b] < v[b]) { mnv = v; break; }
                            }
                        }
                    }
                }
            }
            for (int b = 1; b < 10; b++)
                s += mnv[b] * b;
        }
        return s;
    }
}
ll Q257() {
    set<tuple<ll, ll, ll>>v;
    ll a, b, c, cnt = 0, k, s = 100000000;
    vl p = { 3,5 }, q = { 4,8 };
    for (ll i = 0; i < 2; i++)
    {
        for (ll m = 1; m < s * q[i] / (p[i] - 1); m++)
        {
            for (ll n = m / p[i] + 1; n < m; n++)
            {
                if ((p[i] - 1) * n * (m + n) > q[i] * s)break;
                if (gcd(m, n) != 1)continue;
                b = 2 * p[i] * n * n - 2 * m * n;
                c = m * m - n * n;
                a = (m - n) * (p[i] * n - m);
                if (a > min(b, c))continue;
                k = (gcd(a, gcd(b, c)));
                if (b > c)swap(b, c);
                if (a + b <= c)continue;
                if ((a + b + c) / k > s) continue;
                v.insert({ a / k,b / k,c / k });
            }
        }
    }
    for (auto p : v)
        cnt += s / (get<0>(p) + get<1>(p) + get<2>(p));
    return cnt + s / 3;
}
namespace Q263
{
    vector<ll> v, pr;
    bool f(int m) {
        if (m + 9 > v.size())return 0;
        return v[m - 9] == m - 9 && v[m - 3] == m - 3 && v[m + 3] == m + 3 && v[m + 9] == m + 9
            && v[m - 7] != m - 7 && v[m - 5] != m - 5 && v[m - 1] != m - 1
            && v[m + 1] != m + 1 && v[m + 5] != m + 5 && v[m + 7] != m + 7;
    }
    bool practical(ll m) {
        ll s = 1;
        while (m != 1) {
            ll p = v[m], c = 1;
            while (v[m] == p)m /= p, c *= p;
            if (s < p - 1)return 0;
            s *= ~- (c * p) / ~- p;
        }
        return 1;
    }
    ll solve() {
        ll n = 1200000000;
        v = vector<ll>(n);
        for (ll i = 2; i < n; i++) {
            if (!v[i])pr.push_back(v[i] = i);
            for (ll j = 0; j < pr.size() && pr[j] <= v[i] && pr[j] * i < n; j++)
                v[pr[j] * i] = pr[j];
        }
        set<ll> prs;
        for (ll i = 0; i < n / 840; i++) {
            if (f(840 * i + 20))prs.insert(840 * i + 20);
            if (f(840 * i + 400))prs.insert(840 * i + 400);
            if (f(840 * i + 440))prs.insert(840 * i + 440);
            if (f(840 * i + 820))prs.insert(840 * i + 820);
            ll s = 0, c = 0;
        }
        ll s = 0, c = 0;
        for (ll i : prs) {
            if (practical(i - 8) && practical(i - 4) && practical(i) && practical(i + 4) && practical(i + 8))
                c++, s += i;
            if (c == 4)break;
        }
        return s;
    }
}
namespace Q265//1,529,166 microseconds
{
    ll solve()//we have somewhere 100001 so we start from this (000001......1)
    {
        int n = 5, flag;
        unsigned int a;
        ll sum = 0;
        for (unsigned int i = (1 << ((1 << n) - n - 1)) + 1; i < (1 << ((1 << n) - n)); i += 2)
        {
            vector<int>v(32);
            flag = 1;
            for (int j = (1 << n) - 1; j >= n - 1; j--)
            {
                a = (i << (31 - j)) >> (32 - n);
                if (v[a])
                {
                    flag = 0;
                    break;
                }
                v[a] = 1;
            }
            if (flag)for (int j = n - 2; ~j; j--)
            {
                a = (((i << (31 - j)) >> (31 - j)) << (n - 1 - j));//from the left are all zeroes so it works
                if (v[a])
                {
                    flag = 0;
                    break;
                }
                v[a] = 1;
            }
            if (flag)sum += i;
        }
        return sum;
    }
}
namespace Q267
{
    //python
}
ll Q268()
{
    vector<ll>p = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97 };
    ll sum = 0, tmp, mlt, j, cnt = 1;
    ll lim = 1e16;
    vector<ll>tri(26, 0);
    for (ll i = 4; i < 26; i++)
    {
        tri[i] = abs(tri[i - 1]) + cnt;
        if (i & 1)tri[i] *= -1;
        cnt += (i - 2);
    }
    for (ll i = 1; i < power(2, 25); i++)
    {
        j = 0; tmp = i;
        mlt = 1;
        while (tmp && mlt <= lim)
        {
            mlt *= p[log(tmp & (-tmp), 2) - 1];
            tmp = tmp ^ (tmp & (-tmp));
            j++;
        }
        sum += (lim / mlt) * tri[j];
    }
    return sum;
}
namespace Q269//3,995,711 microseconds
{//the template creates k dimensional vector and can access negative indexes
    template<typename T>
    struct Vector {
        vector<T> v;
        Vector(ll sz, T t) : v(sz * 2, t) {}
        T& operator[](ll i) { return v[i + v.size() / 2]; }
    };
    Vector<ll> init(ll sz) {
        return Vector<ll>(sz, 0);
    }
    template<typename... Args>
    auto init(ll sz, Args... a) {
        auto V = init(a...);
        return Vector<decltype(V)>(sz, V);
    }
    vl mult(vl& a, vl& b)
    {
        vl c(a.size() + b.size() - 1);
        for (ll i = 0; i < c.size(); i++)for (ll j = max(i - i, ll(i - a.size() + 1)); j <= min(i, ll(b.size()) - 1); j++)
        {
            c[i] += a[i - j] * b[j];
        }
        return c;
    }
    ll solve()
    {
        ll n = 1e16, sum = n / 10, m, k = 16;//for P(0)=0
        for (ll i = 1; i < 10; i++)
        {
            m = 1000;
            vl v = { i,1 };
            vector<vector<ll>>dp(k, vl(m));//dp[a][b] ~ (b-m/2)x^a (x+1)(C+C1x+...+Cax^a)
            dp[0][m / 2] = 1;
            for (ll a = 0; a < k - 1; a++)
            {
                for (ll b = 0; b < m; b++)
                {
                    if (dp[a][b] == 0)continue;
                    ll d = b - m / 2;//0<=d+i*x<=9
                    if (a == 10 && b == 5)
                    {
                        ll acs = 2;
                    }
                    for (ll x = -d / i; x <= (9 - d) / i; x++)
                    {
                        if (d + x * i >= 0 && d + x * i <= 9 && (d + x * i != 0 || a != 0))dp[a + 1][x + m / 2] += dp[a][b];
                    }
                }
            }
            for (ll a = m / 2; a < m / 2 + 10; a++)
            {
                sum += dp[k - 1][a];
            }
        }
        for (ll i = 1; i < 10; i++)for (ll j = i + 1; j < 10; j++)
        {
            m = 1000;
            if (i * j >= 10)break;
            vl a = { i,1 }, b = { j,1 };
            vector<vvl>dp(k, vvl(m, vl(m))); dp[0][m / 2][m / 2] = 1;
            for (ll a = 0; a < k - 2; a++)
            {
                for (ll b = 0; b < m; b++)
                {
                    for (ll c = 0; c < m; c++)
                    {
                        if (dp[a][b][c] == 0)continue;
                        ll d = b - m / 2, e = c - m / 2;
                        for (ll x = -e / (i * j); x <= (9 - e) / (i * j); x++)
                        {
                            if (e + x * i * j >= 0 && e + x * i * j <= 9 && (e + x * i * j != 0 || a != 0))dp[a + 1][x + m / 2][b + x * (i + j)] += dp[a][b][c];
                        }
                    }
                }
            }
            for (ll a = m / 2; a < m / 2 + 10; a++)
            {
                for (ll b = m / 2; b < m / 2 + 10; b++)
                {
                    sum -= dp[k - 2][a][b];
                }
            }
        }
        for (ll i = 1; i < 10; i++)for (ll j = i + 1; j < 10; j++)for (ll l = j + 1; l < 10; l++)
        {
            ll m = 80;
            if (i * j * l > 10)continue;
            auto dp = init(k, m, m, m); dp[0][0][0][0] = 1;
            for (ll a = 0; a < k - 3; a++)
            {
                for (ll b = -m; b < m; b++)
                {
                    for (ll c = -m; c < m; c++)
                    {
                        for (ll d = -m; d < m; d++)
                        {
                            if (dp[a][b][c][d] == 0)continue;
                            for (ll x = -d / (i * j * l); x <= (9 - d) / (i * j * l); x++)
                            {
                                if (d + x * i * j * l >= 0 && d + x * i * j * l <= 9
                                    && (d + x * i * j * l != 0 || a != 0))dp[a + 1][x][b + (i + j + l) * x][c + (i * j + i * l + j * l) * x] += dp[a][b][c][d];
                            }
                        }
                    }
                }
            }
            for (ll a = 0; a < 10; a++)
            {
                for (ll b = 0; b < 10; b++)
                {
                    for (ll c = 0; c < 10; c++)
                    {
                        sum += dp[k - 3][a][b][c];
                    }
                }
            }
        }
        return sum;
    }
}
ll Q271()
{
    ll n = 13082761331670030;
    ll sum = 0, b, q;
    vector<ll>p = { 7,13,19,31,37,43 };
    ll s = 2ll * 3ll * 5ll * 11ll * 17ll * 23ll * 29ll * 41ll;
    vector<vector<ll>>v(6);
    for (ll i = 0; i < 6; i++)
    {
        for (ll j = 1; j < p[i]; j++)
        {
            if (j * j * j % p[i] == 1)v[i].push_back(j);
        }
    }
    for (ll i = 0; i < 729; i++)
    {
        ll tmp = i;
        vector<ll>temp;
        while (temp.size() < 6)
        {
            temp.push_back(tmp % 3);
            tmp /= 3;
        }
        q = p[0]; b = v[0][temp[0]];
        for (ll j = 1; j < 6; j++)
        {
            b = crt(p[j], q, v[j][temp[j]], b);
            q *= p[j];
        }
        sum += crt(s, q, 1, b);
    }
    return sum - 1;
}
ll Q272()
{
    ll n = ll(1e11), flag = 0, j, k, a = 7 * 9 * 13 * 18;
    ll sum = 0;
    vector<ll>num;
    pl tmp;
    vector<ll>prime;
    vector<ll>div(n / a);
    vector<bool>p(n / a);
    sieve(prime, div, p, n / a);
    vector<ll>p1, p2;//primes 1 mod 3 and 2 mod 3
    for (ll i = 0; i < prime.size(); i++)
    {
        if (prime[i] % 3 == 1)p1.push_back(prime[i]);
        if (prime[i] % 3 == 2)p2.push_back(prime[i]);
    }
    queue<pl>s;
    queue<ll>ans;
    s.push({ 1,0 });
    tmp = { 9,1 };
    while (tmp.first <= n)
    {
        s.push(tmp);
        tmp.first *= 3;
    }
    for (ll i = 0; i < p1.size(); i++)
    {
        flag = 0;
        while (s.front().first != 1 || !flag)
        {
            tmp = s.front();
            s.pop();
            if (tmp.second == 5)
            {
                ans.push(tmp.first);
                continue;
            }
            if (tmp.first * p1[i] > n)continue;
            s.push(tmp);
            tmp.second++;
            while (tmp.first * p1[i] <= n)
            {
                tmp.first *= p1[i];
                s.push(tmp);
            }
            flag = 1;
        }
        if (i % 1000 == 0)cout << i << "\n";
    }
    cout << ans.size() << "\n";
    for (ll i = 0; i < p2.size(); i++)
    {
        k = ans.size();
        while (k > 0)
        {
            k--;
            j = ans.front();
            ans.pop();
            if (j * p2[i] > n)
            {
                if (j % 3 != 0 && 3 * j <= n)sum += 4 * j;
                else sum += j;
                continue;
            }
            ans.push(j);
            while (j * p2[i] <= n)
            {
                j *= p2[i];
                ans.push(j);
            }
        }
        if (i % 1000 == 0)cout << i << "\n";
    }
    while (!ans.empty())
    {
        j = ans.front();
        ans.pop();
        if (j % 3 != 0 && 3 * j <= n)sum += 4 * j;
        else sum += j;
    }
    return sum;
}
namespace Q288//1,251,438 microseconds
{
    ll solve()
    {
        ll n = power(10, 7), k = 61, p = power(k, 10), s, sum = 0;
        vl v(n + 11);
        s = 290797;
        for (ll i = 0; i <= n; i++)
        {
            v[i] = s % k;
            s = s * s % 50515093;
        }
        for (ll i = 1; i <= n; i++)
        {
            for (ll j = 0; j < 10; j++)
            {
                sum += v[i + j] * power(k, j);
            }
            sum %= p;
        }
        return sum;
    }
}
namespace Q291//629,360,125 microseconds so ~ 10.5 min
{
    ll solve()//mathing and only 2x^2+2x+1 can be prime
    {
        ll n = 5e15, sum = 0;
        for (ll i = 1; 2 * i * i + 2 * i + 1 < n; i++)
        {
            if (is_prime(2 * i * i + 2 * i + 1))
            {
                sum++;
            }
        }
        return sum;
    }
}
namespace Q297//358 microseconds
{
    pl s(ll n, vl& f, vl& c)
    {
        ll i;
        if (n == 1)return { 0,1 };
        for (i = 0; f[i] < n; i++); i--;
        pl ans = s(n - f[i], f, c);
        return { ans.first + c[i - 1] + ans.second,n };
    }
    ll solve()
    {
        ll n = 1e17;
        vl f = { 1,2 }, c = { 1,2 };
        for (ll i = 2;; i++)
        {
            f.push_back(f[i - 1] + f[i - 2]);
            c.push_back(c[i - 1] + c[i - 2] + f[i - 1]);
            if (f[i] >= n)break;
        }
        return s(n, f, c).first;
    }
}
namespace Q301//530,310 microseconds
{
    ll solve()
    {
        ll sum = 0;
        for (ll i = 1; i <= (1 << 30); i++)
        {
            if ((i ^ (i * 2) ^ (i * 3)) == 0)sum++;
        }
        return sum;
    }
}
ll Q302()//not working
{
    ll n = 1e8;
    ll cnt = 0, tmp, cur, curp, flag;
    vector<ll> prime;
    vector<ll> div(n);
    vector<bool> p(n);
    sieve(prime, div, p, n);
    for (ll i = 2; i < n; i++)
    {
        flag = 1;
        tmp = i;
        while (tmp > 1)
        {
            curp = div[tmp];
            cur = 0;
            while (div[tmp] == curp)
            {
                tmp /= curp;
                cur++;
            }
            if (cur == 1)
            {
                flag = 0;
                break;
            }
        }
        if (flag && cur > 2)cnt++;
    }
    return cnt;
}
namespace Q306
{
    ll mex(vector<ll>& v, ll i)
    {
        vector<ll>m(i + 1);
        for (ll j = 0; j <= (i - 2) / 2; j++)
        {
            m[v[j] ^ v[i - 2 - j]] = 1;
        }
        for (ll j = 0; j <= i; j++)if (!m[j])return j;
    }
    ll solve()
    {
        ll n = 1e6, cnt = 0, a = 0;
        vector<ll>v(105);
        for (ll i = 2; i <= 104; i++)
        {
            v[i] = mex(v, i);
        }
        vector<ll>period(34);
        for (ll i = 0; i < 34; i++)period[i] = v[68 + i];
        /*for (ll i = 0; i <= n; i++)
        {
            for (ll j = 0; j < 34&&i<=n; j++, i++)cout << v[i] << "  ";//after 3rd row its periodic
            cout << "\n\n";
        }*/
        for (ll i = 1; i < 68; i++)
        {
            if (!v[i])cnt++;
        }
        for (ll x : period)if (!x)a++;
        cnt += ((n / 34) - 2) * a;
        for (ll i = 0; i <= (n % 34); i++)if (!period[i])cnt++;
        return n - cnt;
    }
}
namespace Q310//3,234,762 microseconds
{
    ll solve() {
        ll n = 1e5;
        vl v(n + 1), h(n + 1);
        h[0]++;
        for (ll i = 1; i <= n; i++) {
            vl s;
            for (ll j = 1; j * j <= i; j++)
                s.push_back(v[i - j * j]);
            sort(s.begin(), s.end());
            s.push_back(s[s.size() - 1] + 2);
            ll m = -1;
            for (auto a : s) {
                if (m < a - 1) {
                    m++;
                    break;
                }
                m = a;
            }
            h[v[i] = m]++;
        }
        ll ans = h[0] * (-~n * 3 + 2);
        for (ll a = 0; a <= n; a++)
            for (ll b = 0; b <= n; b++)
                ans += h[v[a] ^ v[b]];
        return ans / 6;
    }
}
namespace Q323//219 microseconds
{
    long double solve()
    {
        ll n = 32;
        vector<long double> e(n + 1); e[0] = 0;
        vvl ch(n + 1);
        for (ll i = 0; i <= n; i++)
        {
            ch[i].resize(i + 1, 1);
            for (ll j = 1; j < i; j++)ch[i][j] = (ch[i - 1][j - 1] + ch[i - 1][j]);
        }
        for (ll i = 1; i <= n; i++)
        {
            e[i] = 1;
            for (ll j = 0; j < i; j++)
            {
                e[i] += e[j] / (1ll << i) * ch[i][j];
            }
            e[i] /= double((1ll << i) - 1) / (1ll << i);
        }
        return e[n];
    }
}
namespace Q346
{
    ll solve()
    {
        ll n = 1e12, tmp, sum = 0;
        set<ll>s;
        for (ll i = 2; i * i < n; i++)
        {
            tmp = i * i + i + 1;
            while (1)
            {
                if (tmp < n)
                {
                    s.insert(tmp);
                    tmp *= i; tmp += 1;
                }
                else break;
            }
        }
        for (auto x : s)sum += x;
        return sum + 1;//+1 for the number 1
    }
}
namespace Q347
{
    ll big(ll n, ll p, ll q)
    {
        ll a = n / p / q, cur, m = 1;
        for (ll i = 0;; i++)
        {
            cur = power(p, i);
            if (cur > a)break;
            m = max(m, cur);
            for (ll j = 1;; j++)
            {
                cur *= q;
                if (cur > a)break;
                m = max(m, cur);
            }
        }
        return m * p * q;
    }
    ll solve()
    {
        ll n = 1e7, p, q, sum = 0;
        vector <ll> prime;
        vector<ll>div(n);
        vector<bool>c(n);
        sieve(prime, div, c, n);
        for (ll i = 0; i < prime.size(); i++)for (ll j = i + 1; j < prime.size(); j++)
        {
            p = prime[i], q = prime[j];
            if (p * q > n)break;
            sum += big(n, p, q);
        }
        return sum;
    }
}
namespace Q348//123,077 microseconds
{
    ll pali(ll n, ll i)
    {
        ll ans = n * power(10, i / 2), m = 0;
        if (i % 2)n /= 10;
        while (n)
        {
            m *= 10; m += n % 10;
            n /= 10;
        }
        return ans + m;
    }
    ll solve()
    {
        ll n = 5, k = 0, sum = 0, cnt = 0, s, e, p, d;
        for (ll i = 1; k < n; i++)
        {
            s = power(10, (i - 1) / 2); e = 10 * s;
            for (ll j = s; j < e; j++)
            {
                p = pali(j, i);
                cnt = 0;
                for (ll a = pow(p, 0.333333333); a; a--)
                {
                    d = sqrt(p - a * a * a);
                    if (d * d + a * a * a == p)
                    {
                        cnt++;
                    }
                }
                if (cnt == 4)
                {
                    sum += p;
                    k++;
                }
            }
        }
        return sum;
    }
}
namespace Q351//1,903,322 microseconds
{
    ll solve()
    {
        ll n = 1e8, sum = 0, tmp;
        vl div(n + 1), prime;
        vector<bool>p(n + 1);
        sieve(prime, div, p, n + 1);
        vl phi(n + 1); phi[1] = 1;
        for (ll i = 2; i <= n; i++)
        {
            if (p[i])phi[i] = i - 1;
            else
            {
                tmp = i / div[i];
                if (div[i] == div[tmp])phi[i] = div[i] * phi[tmp];
                else phi[i] = (div[i] - 1) * phi[tmp];
            }
            sum += (i - phi[i]);
        }
        return 6 * sum;
    }
}
namespace Q357
{
    ll check(ll a, vector<ll>& div, vector<bool>& p)
    {
        vector<ll>d;
        ll tmp = a, b, j;
        if (!p[a + 1])return 0;
        if (!p[a / 2 + 2])return 0;
        while (tmp > 1)
        {
            if (!d.size() || d[d.size() - 1] != div[tmp])
            {
                d.push_back(div[tmp]);
                tmp /= div[tmp];
            }
            else return 0;
        }
        for (ll i = 1; i < (1 << d.size()) - 1; i += 2)
        {
            b = 1; tmp = i, j = 0;
            while (tmp)
            {
                if (tmp & 1)b *= d[j];
                tmp /= 2; j++;
            }
            if (!p[b + a / b])return 0;
        }
        return a;
    }
    ll solve()
    {
        ll n = 1e8 + 2, sum = 0;
        vector <ll> prime;
        vector<ll>div(n);
        vector<bool>p(n);
        sieve(prime, div, p, n);
        for (ll i = 2; i <= n - 2; i += 4)
        {
            sum += check(i, div, p);
        }
        return sum + 1;//+1 for i=1
    }
}
namespace Q365//2541214 microseconds
{
    ll factmod(ll n, ll p) {
        vl f(p);
        f[0] = 1;
        for (ll i = 1; i < p; i++)
            f[i] = power(f[i - 1], i, p, [](ll a, ll b, ll c) {return (a + b) % c; }, 0);
        ll res = 1;
        while (n > 1) {
            if ((n / p) % 2)
                res = p - res;
            res = power(res, f[n % p], p, [](ll a, ll b, ll c) {return (a + b) % c; }, 0);
            n /= p;
        }
        return res;
    }
    ll check(ll n, ll m, ll p)
    {
        ll ans = 0, tmp = p;
        while (tmp <= n)
        {
            ans += n / tmp - (m / tmp) - (n - m) / tmp;
            tmp *= p;
        }
        return ans;
    }
    ll chinese(ll p, ll r, ll a, ll b)//finds m=a mod p,m=b mod r,and p is prime
    {
        ll m = (p + (a - b) % p) * power(r, p - 2, p) % p;
        return b + m * r;
    }
    ll solve()
    {
        ll n = 1e18, x = 0, m = 1e9, s = 1e3, e = 5e3, a, b, c, sum = 0;
        vl prime; nsieve(prime, e);
        vl res(prime.size());
        for (; prime[x] <= s; x++);
        for (ll i = x; i < prime.size(); i++)
        {
            if (check(n, m, prime[i]))
            {
                res[i] = 0; continue;
            }
            a = factmod(n, prime[i]), b = factmod(m, prime[i]), c = factmod(n - m, prime[i]);
            res[i] = a * power(b, prime[i] - 2, prime[i]) * power(c, prime[i] - 2, prime[i]) % prime[i];
        }
        for (ll i = x; i < prime.size(); i++)
        {
            for (ll j = i + 1; j < prime.size(); j++)for (ll r = j + 1; r < prime.size(); r++)
            {
                sum += chinese(prime[i], prime[j] * prime[r], res[i], chinese(prime[j], prime[r], res[j], res[r]));
            }
        }
        return sum;
    }
}
ll Q375()
{
    ll a = 290797;
    ll num = a;
    ll mod = 50515093;
    ll n; cin >> n;
    pl tmp;
    stack<pl>sl;
    sl.push({ -1e18,0 });
    vector<ll>l(n + 1);
    vector<ll>r(n + 1);
    vector<ll>v(n + 1);
    for (ll i = 1; i <= n; i++)
    {
        num = (num * num) % mod;
        v[i] = num;
        while (sl.top().first > num)
        {
            sl.pop();
        }
        l[i] = i - sl.top().second;
        sl.push({ num,i });
    }
    stack<pl>sr;
    sr.push({ -1e18,n + 1 });
    for (ll i = n; i > 0; i--)
    {
        while (sr.top().first >= v[i])
        {
            sr.pop();
        }
        r[i] = sr.top().second - i;
        sr.push({ v[i],i });
    }
    ll sum = 0;
    for (ll i = 1; i <= n; i++)
    {
        sum += v[i] * l[i] * r[i];
    }
    return sum;
}
namespace Q381
{
    ll solve()
    {
        ll n = 1e8, tmp, sum = 0;
        vector <ll> prime, inv(5);
        vector<ll>div(n);
        vector<bool>q(n);
        sieve(prime, div, q, n);
        for (ll p : prime)
        {
            if (p < 5)continue;
            inv[0] = p - 1;
            for (ll i = 1; i < 5; i++)
            {
                inv[i] = power(p - i, p - 2, p) * inv[i - 1]; inv[i] %= p;
            }
            tmp = 0;
            for (ll j = 0; j < 5; j++)tmp += inv[j];
            sum += tmp % p;
        }
        return sum;
    }
}
namespace Q387//16,496 microseconds
{
    ll solve()
    {
        ll sz = 0, flag = 1, n = 1e14, i = 0, tmp, sum = 0;
        vector<pl> had;
        vl a;
        for (ll i = 0; i < 9; i++)had.push_back({ i + 1,i + 1 });
        while (flag)
        {
            for (i = sz; i < had.size(); i++)
            {
                if (had[i].first * 100 > n)
                {
                    flag = 0; break;
                }
                tmp = 10 * had[i].first;
                for (ll j = 0; j < 10; j++)
                {
                    if ((tmp + j) % (had[i].second + j))continue;
                    had.push_back({ tmp + j,had[i].second + j });
                }
            }
            sz = i;
        }
        for (auto x : had)
        {
            if (!is_prime(x.first / x.second))continue;
            for (ll j = 1; j < 10; j += 2)
            {
                if (10 * x.first + j > n)continue;
                if (is_prime(10 * x.first + j))
                {
                    a.push_back(10 * x.first + j);
                    sum += a[a.size() - 1];
                }
            }
        }
        return sum;
    }
}
namespace Q421
{
    ll generator(vector<ll>& div, ll p)
    {
        vector<ll>d;
        ll n = p - 1, cur = -1;
        while (n > 1)
        {
            if (cur != div[n])d.push_back(div[n]);
            cur = div[n]; n /= div[n];
        }
        for (ll i = 2; i < p; i++)
        {
            bool flag = 1;
            for (ll j = 0; j < d.size(); j++)
            {
                if (power(i, (p - 1) / d[j], p) == 1)
                {
                    flag = 0; break;
                }
            }
            if (flag)return i;
        }
        return -1;
    }
    ll bigger(vector<ll>& v, ll p, ll n)
    {
        sort(v.begin(), v.end());
        ll sum = (n / p) * v.size();
        n %= p;
        for (ll i = 0; i < v.size(); i++)
        {
            if (v[i] > n)return sum;
            sum++;
        }
        return sum;
    }
    ll solve()
    {
        ll n = 1e8, k = 1e11, p, g, cur, sum = 0;
        vector<ll> prime;
        vector<ll>div(n);
        vector<bool>q(n);
        sieve(prime, div, q, n);
        for (ll i = 0; i < prime.size(); i++)
        {
            vector<ll>v;
            p = prime[i];
            if (p % 3 != 1 && p % 5 != 1)
            {
                v.push_back(p - 1);
                sum += p * bigger(v, p, k);
                continue;
            }
            g = generator(div, p);
            if (p % 3 == 1 && p % 5 != 1)
            {
                g = power(g, (p - 1) / 3, p); cur = g;
                for (ll j = 0; j < 3; j++)
                {
                    v.push_back(p - cur); cur *= g; cur %= p;
                }
                sum += p * bigger(v, p, k);
                continue;
            }
            if (p % 3 != 1 && p % 5 == 1)
            {
                g = power(g, (p - 1) / 5, p); cur = g;
                for (ll j = 0; j < 5; j++)
                {
                    v.push_back(p - cur); cur *= g; cur %= p;
                }
                sum += p * bigger(v, p, k);
                continue;
            }
            if (p % 3 == 1 && p % 5 == 1)
            {
                g = power(g, (p - 1) / 15, p); cur = g;
                for (ll j = 0; j < 15; j++)
                {
                    v.push_back(p - cur); cur *= g; cur %= p;
                }
                sum += p * bigger(v, p, k);
                continue;
            }
        }
        return sum;
    }
}
ll Q425()//idk if works
{
    ll n = 1e3, tmp, len, sum = 0, sum1 = 0;
    vector<ll>prime;
    vector<ll>div(n);
    vector<bool>p(n);
    vector<bool>seen(n);
    sieve(prime, div, p, n);
    vector<vector<ll>>G(n);
    for (ll i : prime)
    {
        len = log(i, 10);
        if (len < 3) {
            for (ll j = 1; j < 10; j++)
            {
                tmp = power(10, len);
                if (p[i + j * tmp])G[i].push_back(i + j * tmp);
            }
        }
        for (ll j = 0; j < len; j++)
        {
            tmp = (i % power(10, j + 1)) / power(10, j);
            for (ll k = 0; k < 10; k++)
            {
                if (k == 0 && j == len - 1)continue;
                if (k == tmp)continue;
                if (p[i + power(10, j) * (k - tmp)])G[i].push_back(i + power(10, j) * (k - tmp));
            }
        }
    }
    queue<ll>q;
    q.push(2);
    while (!q.empty())
    {
        tmp = q.front(); q.pop();
        if (seen[tmp])continue;
        sum += tmp;
        seen[tmp] = 1;
        for (ll i : G[tmp])
        {
            if (!seen[i])q.push(i);
        }
    }
    for (ll i : prime)sum1 += i;
    return sum1 - sum;
}
ll Q440()
{
    ll L = 2000;
    ll p = 987898789;
    ll x = 445663912;
    ll invx = power(x, p - 2, p);
    ll inv2 = (p + 1) / 2;
    ll a = 5 + x, b = p + 5 - x;
    ll j;
    ll cnt = 0;
    vl poscnt(L + 1);
    vl negcnt(L + 1);
    vvl ppc(L + 1, vl(L + 1));//T(a^n+1)
    vvl npc(L + 1, vl(L + 1));//T(a^n-1)
    vl pc(L + 1), nc(L + 1);
    for (ll i = 1; i <= L; i++)
    {
        poscnt[i]++;
        for (j = i + 1; j <= L; j++)//gcd(c^i+1,c^j+1)
        {
            pl f = { i,1 };
            pl s = { j,1 };
            while (f.first != 0)
            {
                s.first %= 2 * f.first;
                if (s.first >= f.first)
                {
                    s.first -= f.first;
                    s.second *= -f.second;
                }
                swap(f, s);
            }
            if (f.second == 1)cnt += 2 * 11 * (L / 2) + 20 * (L % 2);
            else if (s.second == 1)poscnt[s.first] += 2;
            else negcnt[s.first] += 2;
        }
    }
    for (ll i = 1; i <= L; i++)
    {
        for (j = 1; j <= L; j++)
        {
            ll exp = power(i, j, p - 1);
            exp -= 1; if (exp < 0)exp += p - 1;
            npc[i][j] = power(a, exp, p) - power(b, exp, p);
            npc[i][j] %= p;
            if (npc[i][j] < 0)npc[i][j] += p;
            npc[i][j] *= inv2; npc[i][j] %= p;
            npc[i][j] *= invx; npc[i][j] %= p;
            exp += 2; exp %= (p - 1);
            ppc[i][j] = power(a, exp, p) - power(b, exp, p);
            ppc[i][j] %= p;
            if (ppc[i][j] < 0)ppc[i][j] += p;
            ppc[i][j] *= inv2; ppc[i][j] %= p;
            ppc[i][j] *= invx; ppc[i][j] %= p;
        }
    }
    for (ll i = 1; i <= L; i++)
    {
        for (j = 1; j <= L; j++)
        {
            pc[i] += ppc[j][i];
            pc[i] %= p;
            nc[i] += npc[j][i];
            nc[i] %= p;
        }
    }
    for (ll i = 1; i <= L; i++)
    {
        cnt += nc[i] * negcnt[i] + pc[i] * poscnt[i];
        cnt %= p;
    }
    return cnt;
}
/*namespace Q443//uses gmp
{
    ll solve()
    {
        ll n = 4, x = 13, m = 1000000000000000, p;
        mpz_t d;
        mpz_init(d);
        for (; ++n <= m;) {
            p = x - n;
            mpz_set_ui(d, p);
            if (mpz_probab_prime_p(d, 50) == 2) {
                x += 2 * p - n % p;
                n += p - n % p;
            }
            else {
                x += gcd(n, x);
            }
        }
        return p + m + 1;
    }
}*/
namespace Q457//448,026 microseconds
{
    ll solve()//the answer is like 2e18 lol
    {
        ll n = 1e7, r, a, b, m, t, sum = 0;
        vl div(n), prime;
        vector<bool>q(n);
        sieve(prime, div, q, n);
        for (ll p : prime)
        {
            m = n * n;
            if (p == 2 || p == 13)continue;
            r = stonelli(13, p);
            for (ll i = 0; i < 2; i++)
            {
                if (i == 0)t = r;
                if (i == 1)t = p - r;
                if (t == -1)
                {
                    m = 0; break;
                }
                a = (t * t - 13) / p; a += p; a %= p;
                b = (-a) * (power(2 * t, p - 2, p)) % p; b += p; b %= p;
                t += b * p;
                if (t % 2 == 0)t += p * p; m = min(m, (t + 3) / 2 % (p * p));
            }
            sum += m;
        }
        return sum;
    }
}
namespace Q461//10,759,854 microseconds
{
    ll solve()
    {
        ll n = 1e4, ans = 0, a, b, c, d;
        double p = 3.141592653589793, err = 1;
        ll lim = n * log(p + 1) + 2;//e^(lim/n)-1 will be bigger than pi later by enough to not be considered
        vector<tuple<double, ll, ll>>v((lim + 1) * (lim + 1));
        for (double i = 0; i <= lim; i++)for (double j = i; j <= lim; j++)
        {
            v[(n + 1) * i + j] = make_tuple(exp(i / n) + exp(j / n) - 2, i, j);
        }
        sort(v.begin(), v.end());
        ll s = 0, e = v.size() - 1;
        while (s <= e)
        {
            if (get<0>(v[s]) + get<0>(v[e]) > p)e--;
            else
            {
                if (abs(get<0>(v[s]) + get<0>(v[e]) - p) < err)
                {
                    err = abs(get<0>(v[s]) + get<0>(v[e]) - p);
                    a = get<1>(v[s]); b = get<2>(v[s]);
                    c = get<1>(v[e]); d = get<2>(v[e]);
                    ans = a * a + b * b + c * c + d * d;
                }
                e++;
                if (abs(get<0>(v[s]) + get<0>(v[e]) - p) < err)
                {
                    err = abs(get<0>(v[s]) + get<0>(v[e]) - p);
                    a = get<1>(v[s]); b = get<2>(v[s]);
                    c = get<1>(v[e]); d = get<2>(v[e]);
                    ans = a * a + b * b + c * c + d * d;
                }
                s++;
            }
        }
        return ans;
    }
}
namespace Q493//300 microseconds
{
    double solve()
    {
        ll n = 70, k = 20;
        double e = 1;
        for (ll i = 0; i < 10; i++)
        {
            e *= (5e1 - i) / (7e1 - i);
        }
        return 7 * (1 - e);//linearity of expectation
    }
}
namespace Q500
{
    ll solve()
    {
        ll k = 500500, n = 1e7, p = 500500507, ans = 1, cnt = 0, a;
        vector <ll> prime;
        vector<ll>div(n);
        vector<bool>q(n);
        sieve(prime, div, q, n);
        priority_queue < pair<double, ll>, vector<pair<double, ll>>, greater<pair<double, ll>>> s;
        for (ll i = 0; i < k; i++)
        {
            s.push({ log2(log2(prime[i])),prime[i] });
        }
        while (cnt++ < k)
        {
            pair<double, ll>tmp = s.top();
            s.pop(); s.push({ tmp.first + 1,tmp.second });
            a = (tmp.first + 0.1 - log2(log2(tmp.second)));//+0.1 so the answer wont get rounded from .9999 to below
            ans *= power(tmp.second, power(2, a, p - 1), p); ans %= p;
        }
        return ans;
    }
}
namespace Q504
{
    ll solve()
    {
        ll m = 100, cnt = 0;
        double n, x, a, b, c, d;
        for (a = 1; a <= m; a++)
        {
            for (b = 1; b <= m; b++)
            {
                for (c = 1; c <= m; c++)
                {
                    for (d = 1; d <= m; d++)
                    {
                        n = ((a + c) * (b + d) - (gcd(a, b) + gcd(a, d) + gcd(c, b) + gcd(c, d))) / 2 + 1;
                        x = sqrt(n);
                        if (ll(x) == x)cnt++;
                    }
                }
            }
        }
        return cnt;
    }
}
namespace Q506//212 microseconds
{
    ll solve()//has period 15 cause 123432 is sum 15 so just some sum calculation
    {
        ll n = 1e14, p = 123454321, ph = phi(p), sum = 0, inv, lim, cur;
        vector<ll>v = { 0 }, g = { 123432 };
        string s = "123432", tmp;
        for (ll i = 0; i < 200; i++)s += s[i];
        for (ll i = 1, j = 0; i < 15; i++)
        {
            cur = 0;
            tmp = "";
            while (cur < i)
            {
                tmp += s[j]; cur += s[j] - '0'; j++;
            }
            v.push_back(string_to_int(tmp));
            tmp = "";
            for (ll k = j; k < j + 6; k++)tmp += s[k];
            g.push_back(string_to_int(tmp));
        }
        for (ll i = 0; i < 15; i++)
        {
            if (i > n)break;
            lim = (n - i) / 15; inv = power(1e6 - 1, ph - 1, p);
            sum += v[i] * ((power(10, 6 * lim + 6, p) - 1) * inv % p) % p; sum %= p;
            sum += (((power(10, 6 * lim + 6, p) - 1) * inv % p - (lim % p + 1)) * inv % p + p) * g[i] % p; sum %= p;
        }
        return sum;
    }
}
namespace Q509
{
    //ill solve later(solved last year)
}
namespace Q511//1,847,077 microseconds
{
    vl mult(vl& a, vl& b, ll k, ll p)
    {
        vl ans(k);
        for (ll i = 0; i < k; i++)for (ll j = 0; j < k; j++)
        {
            ans[(i + j) % k] += (a[i] * b[j]);
            ans[(i + j) % k] %= p;
        }
        return ans;
    }
    vl bin(vl& v, ll n, ll k, ll p)
    {
        vl a, b, c;
        cout << n << '\n';
        if (n == 1)return v;
        if (n & 1)
        {
            a = bin(v, n / 2, k, p);
            b = mult(a, v, k, p); b = mult(a, b, k, p);
            return b;
        }
        a = bin(v, n / 2, k, p);
        b = mult(a, a, k, p);
        return b;
    }
    ll solve()//we caculate logrithmic by disecting a sum of n to sum of n/2 and (n+1)/2 ,multiplying them is linear so its log(n)*k*k
    {
        ll n = 1234567898765, k = 4321, d, tmp, cur, p = 1e9;//prime div are 5,41,25343,237631
        vector<ll>vp = { 5,41,25343,237631 };
        vector<ll>v(k), ans;
        for (ll i = 0; i < 1 << 4; i++)
        {
            d = 1; tmp = i; cur = 0;
            while (tmp)
            {
                if (tmp & 1)d *= vp[cur];
                cur++; tmp >>= 1;
            }
            v[d % k]++;
        }
        ans = bin(v, n, k, p);
        return ans[(-n) % k + k];
    }
}
namespace Q512//less than a minute
{
    ll solve()
    {
        ll n = 5e8, sum = 0, tmp;
        vector <ll> prime;
        vector<ll>div(n + 1);
        vector<bool>c(n + 1);
        vector<ll>f(n + 1); f[1] = 1;
        sieve(prime, div, c, n + 1);
        for (ll i = 3; i < n; i += 2)
        {
            tmp = i;
            tmp /= div[tmp];
            if (div[tmp] == div[i])f[i] = f[tmp] * div[i];
            else f[i] = f[tmp] * (div[i] - 1);
            sum += f[i];
        }
        return sum + 1;//1 for f(1)
    }
}
namespace Q513//need much optimization
{
    ll solve()
    {
        ll n = 1e5, cnt = 0, d, e;
        for (ll a = 1; a < n; a++)
        {
            if (a % 100 == 0)cout << a << '\n';
            for (ll b = a; b < n; b += 2)
            {
                for (ll c = b + b % 2; c < min(n + 1, a + b); c += 2)
                {
                    e = 2 * (a * a + b * b) - c * c;
                    d = sqrt(e);
                    if (d * d == e)cnt++;
                }
            }
        }
        return cnt;
    }
}
ll Q521()//check if works
{
    ll mod = power(10, 9), tmp, cnt;
    ll n = power(10, 6);
    vector<ll>prime;
    vector<ll>div(n);
    vector<bool>p(n);
    sieve(prime, div, p, n);
    vector<vector<ll>>S(prime.size() + 1);
    for (int i = 0; i < 1000; i++)S[0].push_back(0);
    for (int i = 0; i < prime.size(); i++)
    {
        tmp = 1; cnt = 0;
        for (int j = i + 1; j < prime.size(); j++)
        {
            S[i + 1].push_back((n * n) / tmp);
        }
    }
    return 0;
}
namespace Q560 {
    ll mod = 1000000007, inv;
    void fwht(vl& v) {
        for (ll h = 1; h < v.size(); h <<= 1) {
            for (ll i = 0; i < v.size(); i += h * 2) {
                for (ll j = i; j < i + h; j++) {
                    ll x = v[j], y = v[j + h];
                    v[j] = x + y;
                    v[j + h] = x - y;
                }
            }
        }
        for (ll& i : v)i %= mod;
    }
    vl mult(vl v1, vl v2) {
        fwht(v1), fwht(v2);
        vl v3(v1.size());
        for (ll i = 0; i < v1.size(); i++)
            v3[i] = v1[i] * v2[i] % mod;
        fwht(v3);
        for (ll& i : v3)i = i * inv % mod;
        return v3;
    }
    vvl dp;
    void f(ll n) {
        if (dp[n].size())
            return;
        f(n / 2), f(n + 1 >> 1);
        dp[n] = mult(dp[n / 2], dp[n + 1 >> 1]);
    }
    ll powmod(ll a, ll b) {
        if (!b)return 1;
        ll x = powmod(a, b / 2);
        return x * x % mod * (b & 1 ? a : 1) % mod;
    }
    ll solve(ll x, ll n) {//need to input x=n=1e7
        dp = vvl(n + 1);
        vl pr, v(x);
        for (int i = 2; i < x; i++) {
            if (!v[i])pr.push_back((v[i] = pr.size() + 1, i));
            for (int j = 0; j < pr.size() && pr[j] * i < x && pr[j] <= pr[v[i] - 1]; j++)
                v[pr[j] * i] = j + 1;
        }
        ll n2 = 1;
        while (n2 <= pr.size())n2 <<= 1;
        inv = powmod(n2, mod - 2);
        dp[1] = vl(n2);
        for (ll i = 2; i < x; i++)
            dp[1][v[i] - 1 ? v[i] : 0]++;
        dp[1][1]++;
        f(n);
        return dp[n][0];
    }
}
/*namespace Q591 {//yoad says around 22 min
    ld pi = 0.141592653589793238462643383279;
    inline ld frac(ld x) { ld y; mpf_floor(y.get_mpf_t(), x.get_mpf_t()); return x - y; }
    ll solved(ld d, ll n, ll n1, ll m) {
        vector<pair<ld, ll>> v(m + 1);
        for (ll i = 1; i <= m; i++) {
            v[i].first = frac(v[i - 1].first + d);
            v[i].second = i;
        }
        ld mx = v[m].first;
        sort(v.begin(), v.end());
        ll b = 0, ansb = 0;
        ld x = 0, ansx = 1;
        for (; b + m < n; b += m) {
            ld dif = pi - x + (x > pi);
            auto a = lower_bound(v.begin(), v.end(), pair<ld, ll>{ dif, 0 });
            a--;
            if (a->first != 0 && ansx > abs(pi - frac(x + a->first))) {
                ansx = abs(pi - frac(x + a->first));
                ansb = a->second + b;
            }
            if (++a != v.end() && a->first && ansx > frac(x + a->first) - pi) {
                ansx = frac(x + a->first) - pi;
                ansb = a->second + b;
            }
            x = frac(x + mx);
        }
        for (; b < n; b++) {
            if (ansx > abs(x - pi))
                ansx = abs(x - pi), ansb = b;
            x = frac(x + d);
        }
        b = -m, x = 1 - mx;
        for (; b > n1; b -= m) {
            ld dif = pi - x + (x > pi);
            auto a = lower_bound(v.begin(), v.end(), pair<ld, ll>{ dif, 0 });
            a--;
            if (ansx > abs(pi - frac(x + a->first))) {
                ansx = abs(pi - frac(x + a->first));
                ansb = a->second + b;
            }
            if (++a != v.end() && ansx > frac(x + a->first) - pi) {
                ansx = frac(x + a->first) - pi;
                ansb = a->second + b;
            }
            x = frac(1 + x - mx);
        }
        for (b += m, x = frac(x + mx); b > n1; b--) {
            if (ansx > abs(x - pi))
                ansx = abs(x - pi), ansb = b;
            x = frac(x - d);
        }
        return ansb;
    }
    ll solve(ll n = 1e13, ll m = 5e6) {
        mpf_set_default_prec(250);
        ll s = 0;
        for (ll d = 2; d < 100; d++) {
            ld sd, d2(d * 1.0);
            mpf_sqrt(sd.get_mpf_t(), d2.get_mpf_t());
            if (abs(sqrt(d) - round(sqrt(d))) < 1e-10)
                continue;
            ll b = solved(sd, mpf_get_d((ld((ld(n + 3.0) + pi) / sd + ld(1.0))).get_mpf_t()), mpf_get_d((ld((pi + ld(3.0 - n)) / sd)).get_mpf_t()), m);
            s += abs(round(mpf_get_d((ld(pi - ld(double(b)) * sd + ld(3.0))).get_mpf_t())));
        }
        return s;
    }
}*/
ll Q628()
{
    ll p = 1008691207, sum = -2, fact = 1, s = 1;
    ll n = 1e8;
    for (ll i = 1; i < n; i++)
    {
        fact *= i;
        fact %= p;
        s += fact;
        s %= p;
        sum -= s;
    }
    sum += 3 * s;
    ll ans = fact * n - sum; ans %= p;
    if (ans < 0)ans += p;
    return ans;
}
namespace Q650//3,848,715 microseconds
{
    ll leg(ll n, ll p)
    {
        ll sum = 0, tmp = p;
        while (n >= tmp)
        {
            sum += n / tmp; tmp *= p;
        }
        return sum;
    }
    ll solve()
    {
        ll n = 2e4, sum = 0, p = 1e9 + 7, cur;
        vl div(n + 1), prime;
        vector<bool>q(n + 1);
        sieve(prime, div, q, n + 1);
        vl v(n + 1);//will keep prime factorization of (2!*...*k!)^2
        for (ll i = 1; i <= n; i++)
        {
            vl pr(n + 1);
            if (i == 5)
            {
                ll asd = 2;
            }
            for (ll x : prime)
            {
                if (x > i)break;
                v[x] += 2 * leg(i, x);
                pr[x] = (i + 1) * leg(i, x) - v[x];
            }
            cur = 1;
            for (ll x : prime)
            {
                if (x > i)break;
                cur *= (power(x, pr[x] + 1, p) - 1) * power(x - 1, p - 2, p) % p; cur %= p;
            }
            sum += cur; sum %= p;
        }
        return sum;
    }
}
ll Q665()//check if works
{
    ll n = 1e2, flag, j, sum = 0, cnt = 0;
    vector<vector<bool>>v(1000, vector<bool>(1000));
    set<ll>a, b, c, d;//each will hold values of illegal tiles
    for (ll k = 1; k < n; k++)
    {
        sum++;
    }
    return 0;
}
ll Q684()
{
    ll p = 1000000007;
    ll a = 0, b = 1, sum = 0;
    for (ll i = 2; i <= 90; i++)
    {
        a += b;
        swap(a, b);
        sum += p - b % p - 1;
        sum += 45 * (((power(10, b / 9, p) - 1) * power(9, p - 2, p)) % p);
        sum %= p;
        sum += (((b % 9 + 2) * (b % 9 + 1)) / 2) * power(10, b / 9, p);
        sum %= p;
    }
    return sum;
}
namespace Q686//111,425 microseconds
{
    ll solve()
    {
        ll n = 678910, l = 123, ans = 0, k = 0;
        double a = log10(l), b = log10(l + 1), c = log10(2);
        while (ans != n)//because we look at powers of 2 ans grows by at most one
        {
            ans += ll((b + k) / c) - ll((k + a) / c);
            k++;
        }
        return ll((b + k - 1) / c);
    }
}
namespace Q692//461 microseconds
{
    ll solve()//the smallest to take is the lsb in fib base!
    {
        ll n = 23416728348467685, sum = 0;
        vl fib = { 1,1 };
        for (ll i = 2; fib[i - 1] < n; i++)//n is fib number
        {
            fib.push_back(fib[i - 1] + fib[i - 2]);
        }
        for (ll i = 1; i < fib.size() - 1; i++)
        {
            sum += fib[i] * fib[fib.size() - i - 2];
        }
        sum += fib[fib.size() - 1];
        return sum;
    }
}
namespace Q700//7,882,371 microseconds
{
    ll solve()
    {
        ll k = 1504170715041707, n = 4503599627370517, m = 1, sum = 0, pos = n, tmp;//n prime
        ll kinv = power(k, n - 2, n, [](ll l1, ll l2, ll m1) {return power(l1, l2, m1, [](ll x, ll y, ll z) {return (x + y) % z; }, 0); }, 1);
        while (pos > 1e9)
        {
            tmp = power(m, kinv, n, [](ll x, ll y, ll z) {return (x + y) % z; }, 0);
            if (tmp < pos)
            {
                sum += m; pos = tmp;
            }
            m++;
        }
        cout << 1 << '\n';
        m = n; tmp = 0;
        for (ll i = 1; i < pos; i++)
        {
            tmp += k; tmp %= n;
            if (tmp < m)
            {
                sum += tmp; m = tmp;
            }
        }
        return sum;
    }
}
namespace Q710//95,249 microseconds
{
    ll solve()
    {
        ll n = 1e8, p = 1e6;
        vl t = { 0,0,1,0 }, st = { 0,0,1,0 };
        for (ll i = 4; i <= n; i++)
        {
            t.push_back(st[i - 2] - t[i - 4] + power(2, i / 2 - 2, p));
            t[i] %= p; st.push_back(st[i - 2] + t[i]); st[i] %= p;
            if (t[i] == 0)
            {
                return i;
            }
        }
        return -1;
    }
}
ll Q719()
{
    ll n = 1e4;
    ll sum, ans = 0;
    bool flag;
    for (ll i = 4; i <= n; i++)
    {
        flag = 0;
        for (ll j = 1; j < power(2, log(i * i) - 1); j++)
        {
            sum = 0;
            ll tmp1 = i * i, tmp2 = j;
            vector <ll> v = { -1 };
            for (ll i = 0; i < log(tmp1); i++)
            {
                if (tmp2 % 2)v.push_back(i);
                tmp2 /= 2;
            }
            v.push_back(log(tmp1));
            for (ll i = 0; i < v.size(); i++)
            {
                ans += tmp1 % power(10, v[i + 1] - v[i]);
                tmp1 /= power(10, v[i + 1] - v[i]);
            }
            if (ans == i)
            {
                flag = 1; break;
            }
        }
        if (flag)sum += i * i;
    }
    return sum;
}
namespace Q743//works for k|n but easy to generalize
{
    ll choose(ll a, ll b, ll c, ll p)
    {
        ll ans = a * power(b, p - 2, p); ans %= p;
        return (ans * power(c, p - 2, p)) % p;
    }
    ll solve()//takes about 1min
    {
        ll n = 1e16, k = 1e8, p = 1e9 + 7, sum = 0, tmp;
        vector<ll>fact(k + 1); fact[0] = 1;
        for (ll i = 1; i <= k; i++)
        {
            fact[i] = fact[i - 1] * i; fact[i] %= p;
        }
        for (ll i = k; i >= 0; i -= 2)
        {
            if (i % 1000000 == 0)cout << i << '\n';
            tmp = power(2, (n / k) * i, p) * choose(fact[k], fact[i], fact[k - i], p); tmp %= p;
            tmp *= choose(fact[k - i], fact[(k - i) / 2], fact[(k - i) / 2], p); tmp %= p;
            sum += tmp; sum %= p;
        }
        return sum;
    }
}
namespace Q757
{
    ll solve()//around 2min
    {
        ll n = 1e14, a, d;
        set<ll>s;
        for (ll k = 1; k <= n; k++)
        {
            cout << k << '\n';
            if ((k * k + k) * (k * k + k) > n)return s.size();
            d = k * k;
            for (;; d += k)
            {
                a = (k + 1) * (d + k) * (d / k);
                if (a > n)break;
                s.insert(a);
            }
        }
    }
}
namespace Q769 {
    ll cop_in_pfx(ll k, vl& pr) {
        ll sum = 0;
        for (ll i = 1; i < 1ll << pr.size(); i++) {
            ll sz = 0, mlt = 1, j = i, l = 0;
            while (j) {
                j & 1 ? mlt *= pr[l], sz++ : 0;
                j >>= 1;
                l++;
            }
            sum += k / mlt * (sz % 2 * 2 - 1);
        }
        return k - sum;
    }
    ll cop_in_rng(ll l, ll r, ll div, vl& pr) {
        if (l > r)return 0;
        return cop_in_pfx(r / div, pr) - cop_in_pfx((l + div - 1) / div - 1, pr);
    }
    pl get_rng(ll z, ll m) {
        return { 13 * m * m > z ? ceil(sqrtl(13 * m * m - z)) : 1, (5 - 2 * sqrtl(3)) * m };
    }
    ll solve(ll z = 1e14) {
        ll sum = 0;
        ll n = 22200000;
        vector <ll> prime;
        vector<ll>div(n);
        vector<bool>p(n);
        sieve(prime, div, p, n + 1);
        vl lp = div;
        for (ll m = 2; m < n; m++) {
            vl pr = { lp[m] };
            ll m2 = m / lp[m];
            while (m2 != 1) {
                if (lp[m2] != pr[pr.size() - 1])
                    pr.push_back(lp[m2]);
                m2 /= lp[m2];
            }
            if (m % 13 == 0) {
                if (m % 2) {
                    pl p1 = get_rng(z, m), p4 = get_rng(z * 4, m);
                    ll cnt1_2 = cop_in_rng(p1.first, p1.second, 2, pr),
                        cnt4_1 = cop_in_rng(p4.first, p4.second, 1, pr),
                        cnt4_2 = cop_in_rng(p4.first, p4.second, 2, pr);
                    sum += cnt1_2 + cnt4_1 - cnt4_2;
                }
                else {
                    pl p = get_rng(z, m);
                    sum += cop_in_rng(p.first, p.second, 1, pr);
                }
            }
            else if (m % 2) {
                pl p1 = get_rng(z, m), p4 = get_rng(z * 4, m), p13 = get_rng(z * 13, m), p52 = get_rng(z * 52, m);
                ll cnt1_2 = cop_in_rng(p1.first, p1.second, 2, pr),
                    cnt1_26 = cop_in_rng(p1.first, p1.second, 26, pr),
                    cnt13_26 = cop_in_rng(p13.first, p13.second, 26, pr),
                    cnt4_1 = cop_in_rng(p4.first, p4.second, 1, pr),
                    cnt4_2 = cop_in_rng(p4.first, p4.second, 2, pr),
                    cnt4_13 = cop_in_rng(p4.first, p4.second, 13, pr),
                    cnt4_26 = cop_in_rng(p4.first, p4.second, 26, pr),
                    cnt52_13 = cop_in_rng(p52.first, p52.second, 13, pr),
                    cnt52_26 = cop_in_rng(p52.first, p52.second, 26, pr);
                sum += cnt1_2 - cnt1_26 + cnt13_26 + cnt4_1 - cnt4_2 - cnt4_13 + cnt4_26 + cnt52_13 - cnt52_26;
            }
            else {
                pl p1 = get_rng(z, m), p13 = get_rng(z * 13, m);
                ll cnt1_1 = cop_in_rng(p1.first, p1.second, 1, pr),
                    cnt1_13 = cop_in_rng(p1.first, p1.second, 13, pr),
                    cnt13_13 = cop_in_rng(p13.first, p13.second, 13, pr);
                sum += cnt1_1 - cnt1_13 + cnt13_13;
            }
        }
        return sum + 1;
    }
}
namespace Q788//531,639 microseconds
{
    ll choose(vl& f, ll a, ll b, ll p)
    {
        return f[a] * power(f[b], p - 2, p) % p * power(f[a - b], p - 2, p) % p;
    }
    ll solve()
    {
        ll n = 2022, cur, cnt = 9, p = 1e9 + 7;//9 for 9 first numbers
        vl f(n + 1); f[0] = 1;
        for (ll i = 1; i <= n; i++)
        {
            f[i] = f[i - 1] * i; f[i] %= p;
        }
        for (ll k = 2; k <= n; k++)
        {
            cur = 0;
            for (ll i = k / 2; i < k; i++)
            {
                cur += choose(f, k - 1, i, p) * power(9, k - 1 - i, p); cur %= p;
            }
            cnt += 9 * cur; cnt %= p;
            cur = 0;
            for (ll i = k / 2 + 1; i < k; i++)
            {
                cur += choose(f, k - 1, i, p) * power(9, k - 1 - i, p); cur %= p;
            }
            cnt += 81 * cur; cnt %= p;
        }
        return cnt;
    }
}
namespace Q800//69,499 microseconds
{
    ll solve()
    {
        ll k = 800800, n = k * log2(k), ans = 0, p, q;
        double d = k * log2(k);
        vl c;
        nsieve(c, n);
        for (ll i = 0; i < c.size(); i++)
        {
            p = c[i];
            if (c[i + 1] * log2(p) + p * log2(c[i + 1]) > d)break;
            ll s = i + 1, e = c.size() - 1;
            while (s < e - 1)
            {
                ll m = (s + e) / 2;
                if (c[m] * log2(p) + p * log2(c[m]) <= d)s = m;
                else e = m;
            }
            ans += s - i;
        }
        return ans;
    }
}
ll Q808()
{
    ll n = 1e8, sum = 0, cnt = 0;
    vector<ll>prime, div(n); vector<bool>pr(n);
    sieve(prime, div, pr, n);
    for (ll p : prime)
    {
        if (cnt == 50)break;
        string s = to_string(p * p);
        reverse(s.begin(), s.end());
        if (s == to_string(p * p))continue;
        if (power(ll(sqrt(string_to_int(s))), 2) == string_to_int(s) && pr[ll(sqrt(string_to_int(s)))])
        {
            sum += p * p;
            cnt++;
        }
    }
    return sum;
}
namespace Q820
{
    ll solve()
    {
        ll n = 1e7, sum = 0;
        for (ll i = 1; i <= n; i++)sum += power(10, n, 10 * i) / i;
        return sum;
    }
}
namespace Q834//1,297,778 microseconds
{
    ll solve()
    {
        ll n = 1234567, sum = 0;
        vl div(n + 1), prime;
        vector<bool>q(n + 1);
        sieve(prime, div, q, n + 1);
        for (ll i = 3; i <= n; i++)
        {
            vl d = { 1 };
            ll sz = 1;
            ll tmp = i, divisor = 0, cnt = 0, m;
            for (ll x : {0, 1})
            {
                tmp = i - x;
                while (tmp > 1)
                {
                    divisor = div[tmp]; cnt = 0;
                    while (divisor == div[tmp])
                    {
                        tmp /= divisor; cnt++;
                    }
                    m = divisor;
                    for (ll i = 0; i < cnt; i++)
                    {
                        for (ll j = 0; j < sz; j++)
                        {
                            d.push_back(d[j] * m);
                        }
                        m *= divisor;
                    }
                    sz = d.size();
                }
            }
            for (ll x : d)
            {
                if (x <= i)continue;//m+i|m*(m+1)/2
                x -= i;
                if (x & 1 && i * ((x >> 1) + 1) % (x + i) == 0)sum += x;//to avoid overflow
                if (x % 2 == 0 && (i - 1) * (x / 2) % (x + i) == 0)sum += x;//to avoid overflow
            }
        }
        return sum;
    }
}
string Q836()
{
    return "aprilfoolsjoke";
}
ll Q874()
{
    ll n = 100000, tmp, sum = 0;
    vector <ll> prime;
    vector<ll>div(n);
    vector<bool>p(n);
    sieve(prime, div, p, n);
    ll a = 7000;
    ll b = prime[a];
    ll x = (a - 1) * b;
    x %= a;
    vector <ll> v = { 0 };
    for (ll i = 1; i <= x; i++)
    {
        v.push_back(prime[a - 1] - prime[a - 1 - i]);
    }
    vector<ll>dp(x + 1); dp[1] = v[1];
    for (ll i = 2; i <= x; i++)
    {
        ll m = 1e18;
        for (ll j = 0; j < i; j++)
        {
            m = min(m, dp[j] + v[i - j]);
        }
        dp[i] = m;
    }
    return prime[a - 1] * b - dp[x];
}
namespace Q885//1,234,333 microseconds
{
    ll inc(vl& v, ll k = 18)//generates all incresaingly ordered numbers of k digits
    {
        if (v[1] == k)return 0;
        for (ll i = v.size() - 2; ~i; i--)
        {
            if (v[i] < k)
            {
                v[i]++;
                for (ll j = i + 1; j < v.size() - 1; j++)v[j] = v[i];
                return 1;
            }
        }
    }
    ll fact(ll n)
    {
        if (!n)return 1;
        return  n * fact(n - 1);
    }
    ll c(ll a, ll b)
    {
        return fact(a) / fact(b) / fact(a - b);
    }
    ll choose(vl& v, ll k = 18)
    {
        ll ans = 1;
        for (ll i = 1; i < v.size(); i++)
        {
            ans *= c(k, v[i] - v[i - 1]);
            k -= (v[i] - v[i - 1]);
        }
        return ans;
    }
    ll solve()
    {
        ll k = 18, sum = 0, cur, tmp, p = 1123455689;
        vl v = { 0 };
        for (ll i = 0; i < 9; i++)v.push_back(0);
        v.push_back(k);
        while (1)
        {
            cur = 0;
            for (ll i = 1; i < v.size(); i++)
            {
                tmp = v[i] - v[i - 1];
                while (tmp--)
                {
                    cur *= 10; cur += i - 1;
                }
            }
            cur %= p;
            sum += cur * (choose(v, k) % p); sum %= p;
            if (!inc(v, k))break;
        }
        return sum;
    }
}
ll Q832()
{
    ll p = 1e9 + 7, sum = 0;
    ll p2 = 4, cnt = 1;
    //ll lim = 1000000000000000000;
    ll lim = 1000;
    while (p2 < lim)
    {
        p2 *= 2; cnt++;
    }
    p2 /= 2; cnt--;
    sum += ((((p2 % p) * ((p2 + 1) % p)) % p) * ((p + 1) / 2)) % p;
    ll x = (p2 + 4 * (2 * (cnt % 2) - 1)) / 3;
    sum = sum + (((lim - p2) % p) * ((lim + p2 - x) % p)) % p;
    return sum %= p;
}
namespace Q328 {
    ll solve()
    {
        ll n = 1e3, m, sum = 0;
        vvl dp(n, vl(n));
        for (ll add = 1; add < n; add++)
        {
            for (ll k = 1; k < n - add; k++)
            {
                m = 1e18;
                for (ll i = k; i < k + add; i++)
                {
                    m = min(m, i + max(dp[k][i - 1], dp[i + 1][k + add]));
                }
                dp[k][k + add] = m;
            }
        }
        for (ll x = 2; x <= 100; x++)
        {
            sum += dp[1][x];
            if (x == 100)cout << dp[1][x] << '\n';
        }
        cout << sum << '\n';
        return 0;
    }
}
namespace Q398 {
    long double dbp(long double a, ll b)
    {
        long double ans = 1, tmp;
        while (b != 0)
        {
            tmp = (b & 1) ? a : 1;
            ans *= ans * tmp;
            b /= 2;
        }
        return ans;
    }
    ll solve()
    {
        ll n = 1e7, m = 100;
        vvl dp(n / m, vl(m));
        vector <long double> p(n);
        for (long double i = 0; i < n; i++)
        {
            p[i] = dbp((n - i) / n, m);
        }
        for (long double i = 1; i < n - 1; i++)
        {
            p[i] -= p[i + 1];
        }
        for (ll i = 1000; p[i] > 1e-11; i += 1000)cout << p[i] << " " << i << "\n";
        return 0;
    }
}
/*namespace Q240
{
    ll coin(ll n,ll s,ll x)
    {
        vector<vvl>v(x + 1, vvl(s + 1, vl(s + 1)));
        for (ll i = 0; i <= s; i++)v[0][i][0] = 1;

    }
    ll solve()
    {
        ll n = 20, k = 12, sum = 0;
        for (ll i = 1; i <= k; i++)
        {
            for (ll j = i; j <= k; j++)
            {
                ll a = 70; a -= j + i;
                a -= (n - 2) * j;
            }
        }
    }
}*/
namespace Q770//the code is supposed to work but the big numbers reqiure python
{
    ll solve()
    {
        ll n = 64;
        long double t = 1.8;
        vector<vector<long double>> v(n, vector<long double>(n));
        v[0][0] = 1;
        for (ll i = 1; i < n; i++)
        {
            v[i][0] = v[i - 1][0];
            v[i][0] *= 2;
            v[0][i] = 1;
        }
        for (ll i = 1; i < n; i++)
        {
            for (ll j = 1; j < n; j++)
            {
                long double a = v[i - 1][j], b = v[i][j - 1];
                long double x = (1 - (2 * a) / (a + b));
                v[i][j] = (1 + x) * a;
            }
        }
        for (ll i = 0; i < n; i++)
        {
            if (v[i][i] > t)
            {
                return i;
            }
        }
        return 0;
    }
}
namespace Q678
{
    //for e=2 we have criterion for a^2+b^2=c^odd (prime factor of c)
    //for e,f>=3 we can probably brute force (cause e!=f so itll be like 10^10.5)
}

namespace Q822
{
    ll solve()//log2log2 solves
    {
        return 0;
    }
}
namespace Q472
{
    ll solve()
    {
        ll n = 1e3;
        vector<ll>v = { 0,0,0 };
        for (ll i = 3; i <= n; i++)
        {
            v.push_back(1 + v[(i - 1) / 2] + v[i / 2]);
            cout << i << " " << v[i] << '\n';
        }
        return 0;
    }
}
namespace Q611//not done
{
    ll solve1()
    {
        ll n = 100, sum = 0;
        vector<bool>b(101);
        for (ll i = 1; i < 10; i++)for (ll j = i + 1; j < 10; j++)
        {
            if (i * i + j * j > n)break;
            b[i * i + j * j] = 1 - b[i * i + j * j];
        }
        for (ll i : b)sum += i;
        return sum;
    }
    ll solve()
    {
        ll n = 1000, tmp, cur = 4, p, ans = 0, cnt = 0, index = 0;
        ll k = sqrt(n);
        vector <ll> prime;
        vector<ll>div(k + 2);
        vector<bool>q(k + 2);
        sieve(prime, div, q, k + 2);
        vector<ll>c = { 1 }, p1, p3, a;//p1 primes that arent 3 mod 4, p3 is all numbers whose only factors are 3mod4
        for (ll i = 1; i <= k; i++)
        {
            bool flag = 1;
            tmp = i;
            while (tmp > 1)
            {
                if (div[tmp] % 4 == 1)
                {
                    flag = 0; break;
                }
                tmp /= div[tmp];
            }
            if (flag)p3.push_back(i);
        }
        for (ll i = p3.size() - 1; ~i; i--)
        {
            a.push_back(n / (p3[i] * p3[i]));
            a.push_back(n / (2 * p3[i] * p3[i]));
        }
        sort(a.begin(), a.end());
        for (ll i = 1; i < 1; i++)//sould be 6
        {
            if (i >= prime.size())break;
            tmp = c.size(), p = prime[i];
            for (ll j = 0; j < tmp; j++)
            {
                if (!(c[j] % p))
                {
                    for (ll k = 2; k < p; k++)
                    {
                        c.push_back(c[j] + k * cur);
                    }
                    c[j] = c[j] + cur;
                }
                else
                {
                    for (ll k = 1; k < p; k++)
                    {
                        if ((c[j] + k * cur) % p)c.push_back(c[j] + k * cur);
                    }
                }
            }
            cur *= p;
        }
        sort(c.begin(), c.end());
        for (ll i = 0; i <= k; i++)
        {
            p = prime[i];
            if (p > cur * (k / cur))break;
            else if (prime[i] % 4 == 1)
            {
                while (a[index] < p)
                {
                    ans += cnt; index++;
                }
                cnt++;
            }
        }
        for (ll j = cur * (k / cur); j <= n; j += cur)
        {
            for (ll i : c)
            {
                if (j + i > n)break;
                if (isPrime(j + i, 50))
                {
                    while (a[index] < j + i)
                    {
                        ans += cnt; index++;
                    }
                    cnt++;
                }
            }
        }
        ans += cnt;
        for (ll i = 1; i <= k; i++)
        {
            if (2 * i * i > n)break;
            bool flag = 0;
            tmp = i;
            while (tmp > 1)
            {
                if (div[tmp] % 4 == 1)
                {
                    flag = 1;
                }
                tmp /= div[tmp];
            }
            if (flag)ans++;
        }
        for (ll i = 1; i <= k; i++)
        {
            bool flag = 0;
            tmp = i;
            while (tmp > 1)
            {
                if (div[tmp] % 4 == 1)
                {
                    flag = 1;
                }
                tmp /= div[tmp];
            }
            if (flag)ans++;
        }
        return ans;
    }
}
namespace Q526//need gmp to replace isprime
{
    ll bigp(ll n)
    {
        for (ll i = 2;; i++)
        {
            if (n % i)continue;
            while (!(n % i))n /= i;
            if (n == 1)return i;
        }
    }
    ll sum(ll n, ll t)
    {
        ll g = 0;
        for (ll i = n; i < n + 10; i += 2)
        {
            if (i == n + 4)continue;
            if (!isPrime(i, 35))return -1;
        }
        if (t)
        {
            if (!isPrime((n + 3) / 4, 35))return -1;
            if (!isPrime((n + 5) / 2, 35))return -1;
            if (!isPrime((n + 1) / 6, 35))return -1;
            cout << "2  " << n << '\n';
            g = n + n + 2 + n + 6 + n + 8 + (n + 3) / 4 + (n + 5) / 2 + (n + 1) / 6 + bigp(n + 7) + bigp(n + 4);
            return g;
        }
        else
        {
            if (!isPrime((n + 3) / 2, 35))return -1;
            if (!isPrime((n + 5) / 4, 35))return -1;
            if (!isPrime((n + 7) / 6, 35))return -1;
            cout << "3  " << n << '\n';
            g = n + n + 2 + n + 6 + n + 8 + (n + 3) / 2 + (n + 5) / 4 + (n + 7) / 6 + bigp(n + 1) + bigp(n + 4);
            return g;
        }
    }
    ll solve()//after computations we need n%840=311,521 and assuming the ans is bigger than 4.9n (which makes sense based on smaller numbers))we need 2e11 calcs
    {
        ll n = 1e16, m = 49166659840957733, k = 9, g, t = 0, cnt = 0;//this is a value of m i already found
        ll n1 = n - n % 840 + 311;
        ll n2 = n1 + 210;
        double a = 4 + 1.0 / 2 + 1.0 / 4 + 1.0 / 6 + 1.0 / 315 + 1.0 / 24;
        ll lim = 9909620315061464;
        for (ll i = 9999975546198431 + 840, j = n2;;)
        {
            cnt++;
            if (cnt % 10000000 == 0)cout << cnt << '\n';
            if (t)j -= 840;
            else i -= 840;
            t = 1 - t;
            if (t)
            {
                if ((j + 7) % 9 && (j + 4) % 9)continue;
                g = sum(j, 1);
                if (m < g)
                {
                    cout << j << "  " << g << "\n";
                    lim = g / a;
                    m = g;
                }
            }
            else
            {
                if ((i + 1) % 9 && (i + 4) % 9)continue;
                g = sum(i, 0);
                if (m < g)
                {
                    cout << i << "  " << g << "\n";
                    lim = g / a;
                    cout << lim << "\n";
                    m = g;
                }
            }
        }
        return m;
    }
    ll solve1()
    {
        ll n = 1e2, g = 0, k = 9;
        pl m = { -1,-1 };
        vector <ll> prime;
        vector<ll>div(n + k);
        vector<bool>q(n + k);
        vl bp(n + k);
        sieve(prime, div, q, n + k);
        cout << 1;
        for (ll i = 2; i < n + k; i++)
        {
            if (q[i])bp[i] = i;
            else bp[i] = bp[i / div[i]];
        }
        cout << 2 << "\n";
        for (ll j = 2; j < 2 + k; j++)g += bp[j];
        cout << 3 << "\n\n\n";
        for (ll i = 3; i <= n; i++)
        {
            if (!(i % 1000000))cout << i << "\n";
            g -= bp[i - 1]; g += bp[i + k - 1];
            if (g > m.first)
            {
                m = make_pair(g, i);
            }
        }
        for (ll i = m.second; i < m.second + 9; i++)cout << i << "   " << bp[i] << "   " << i / bp[i] << "\n";
        cout << m.second << "  ";
        return m.first;
    }
}
namespace Q545
{
    ll solve()
    {
        ll sum = 0;
        for (ll i = 1; i < 1e4; i++)
        {
            for (ll j = 1; j < 1e4; j++)
            {
                sum += gcd(i, j);
            }
        }
        return sum;
    }
}
namespace Q110//didnt implement, but go over all sorted lists of powers of primes and calculate, while keeping the best answer.
{//its harder version of 108
    ll solve()
    {
        ll k = 1000, m = 1;
        vl prime = { 2,3,5,7,11,13,17,19,23,31,37,41,43 };//3^14 >4e6
        for (ll i : prime)m *= i;
        vl v(prime.size());
        return 1;
    }
}
/*namespace Q377
{
    ll solve()
    {
        vvl v = { {1} };
        v.push_back({ 0,1 });
        for (ll i = 2; i < 22; i++)
        {
            v.push_back(vl(i + 1));
            for (ll j = 1; j <= i; j++)
            {
                for (ll k = 1; k < min(ll(10), i - j + 2); k++)
                {
                    v[i][j] += v[i - k][j - 1];
                }
            }
        }
        for (vl a : v)
        {
            ll sum = 0;
            for (ll x : a)
            {
                cout << x << " ";
                sum += x;
            }
            cout << sum << '\n';
        }
        return 1;
    }
}*/
//Q853
/*namespace Q61
{
    vector <string> pos;
    void fill(string s, int i, int n)
    {
        string d = "";
        if (i == n - 1)
        {
            pos.push_back(s);
            return;
        }
        fill(s, i + 1, n);
        for (int j = i + 1; j < n; j++)
        {
            if (s[i] == s[j])continue;
            d = "";
            for (int m = 0; m < n; m++)
            {
                if (m < i)d += s[m];
                if (m > i && m != j)d += s[m];
                if (m == i)
                {
                    d += s[j];
                    d += s[i];
                }
            }
            fill(d, i + 1, n);
        }
    }
    vl range(ll s, ll type)
    {
        ll tmp, l = 100 * s, r = 100 * s + 100;
        if (type == 3)
        {
            for (ll i = sqrt(200 * s) - 1;; i++)
            {
                tmp = i * (i + 1) / 2;
                if (tmp >= r)break;
            }
        }
    }
    ll solve()
    {
        string s = "34567";
        vector <int>a(26);
        ll sum = 0;
        for (int i = 0; i < s.size(); i++)
        {
            int b = s[i] - '0';
            a[b]++;
        }
        s = "";
        for (int i = 0; i < 26; i++)
        {
            for (int j = 0; j < a[i]; j++)
            {
                char y = i + '0';
                s += y;
            }
        }
        fill(s, 0, s.size());
        //now pos has all permutations of 3-7
        for (string s : pos)
        {
            vvl v(5, vl(1e4));
            for (ll i = 19; i < 59; i++)
            {
                if (i * (3 * i - 2) % 100 / 10 == 0)continue;
                v[0].push_back(i * (3 * i - 2));
            }
        }
        return 1;
    }
}*/
namespace Q776
{
    ll inc(vl& v)//needs to generate all k digit numbers with first digit less then something
    {
        if (!v[0])return 0;
        for (ll i = v.size() - 1; ~i; i--)
        {
            if (v[i])
            {
                v[i]--;
                for (ll j = i + 1; j < v.size() - 1; j++)v[j] = v[i];
                return 1;
            }
        }
    }
    ll solve()
    {
        ll n = 10, add = 0, cnt = 0, num, k, tmp, cur;
        double avg = 0, sum;
        vl fact(19); fact[0] = 1;
        for (ll i = 1; i < 19; i++)fact[i] = i * fact[i - 1];
        while (n)
        {
            k = log(n);
            vl v(k, 9); v[k - 1] = 10;
            while (inc(v))
            {
                sum = add, num = 1, tmp = k, cur = 0;
                vl digits(10);
                for (ll x : v)
                {
                    sum += x;
                    digits[x]++;
                }
                for (ll i = 0; i < 9; i++)
                {
                    num *= fact[tmp] / fact[digits[i]] / fact[tmp - digits[i]];
                    tmp -= digits[i];
                }
                for (ll i = 0; i < k; i++)
                {
                    cur *= 10; cur++;
                }
                for (ll i = 1; i < 10; i++)
                {
                    if (!digits[i])continue;
                    avg += (cur * i) * (num * digits[i] / k) / sum;
                }
            }
            add += n / power(10, k - 1); n %= power(10, k - 1);
        }
        return avg;
    }
}
namespace Q714//the code works fine but there are bigger numbers than longlong that divide, so should write in python or implement with arrays
{
    pl minimal(ll a, ll b, ll m,ll lim)//xa+yb=0 mod m for 0<x<10 and 0<=y<10
    {
        for (ll x = 1; x <= lim; x++)for (ll y = 0; y < 10; y++)
        {
            if ((x * a + y * b) % m == 0)return { x,y };
        }
        return { 0,0 };
    }
    ll solve()
    {
        ll n = 18, k = 50000, tmp, mod, cur, sum = 0, lim, cnt = 0;
        pl a;
        vl duo(1ll << n);
        vl ans(k + 1);
        for (ll i = 1; i < 100; i++)ans[i] = i;
        for (ll i = 1; i < (1ll << n); i++)
        {
            tmp = i, cur = 1;
            while (tmp)
            {
                if (tmp & 1)duo[i] += cur;
                cur *= 10; tmp >>= 1;
            }
        }
        cout << duo[duo.size() - 1] << '\n';
        for (ll m = 100; m <= k; m++)
        {
            if (ans[m])
            {
                cnt++; continue;
            }
            ll flag = 1;
            lim = 9;
            for (ll i = 1; i < n && flag; i++)
            {
                cur = 1ll << (i + 1); cur--;
                for (ll j = (1ll << i); j <= cur; j++)
                {
                    tmp = duo[j] % m;
                    a = minimal(tmp, (duo[cur] - tmp) % m, m, lim);
                    if (a.first != 0)
                    {
                        if (ans[m])
                        {
                            ans[m] = min(ans[m], a.first * duo[j] + a.second * (duo[cur] - duo[j]));
                        }
                        else ans[m] = a.first * duo[j] + a.second * (duo[cur] - duo[j]);
                        lim = a.first;
                        flag = 0;
                    }
                }
            }
            if (ans[m] == 0)return -1;
            tmp = ans[m] / m;
            for (ll i = 2; i <= k / m; i++)
            {
                if (tmp % i == 0)ans[m * i] = ans[m];
            }
        }
        cout << cnt << '\n';
        for (ll x : ans)
        {
            sum += x;
            if (sum > 3e18)return -1;
        }
        return sum;
    }
}
namespace Q569
{
    ll solve()//just O(n^2) to see if everyone can see only a couple peaks
    {
        ll n = 5e7, k = 10000, sum = 0, cur, mdist = 0;
        vl prime;
        double tmp, m;
        vector<pl>peak(k); peak[0] = { 2,2 };
        nsieve(prime, n);
        for (ll i = 1; i < k; i++)
        {
            peak[i] = { peak[i - 1].first + prime[2 * i - 1] + prime[2 * i],peak[i - 1].second - prime[2 * i - 1] + prime[2 * i] };
        }
        for (ll i = 1; i < k; i++)
        {
            m = 1e9, cur = 0;
            for (ll j = i - 1; ~j; j--)
            {
                tmp = (peak[i].second - peak[j].second) / double(peak[i].first - peak[j].first);
                if (tmp < m)
                {
                    cur++; m = tmp;
                    mdist = max(i - j, mdist);
                }
            }
            sum += cur;
        }
        cout << mdist << '\n';
        return sum;
    }
}
namespace Q312
{
    ll solve()//a_n=b_(n-1)^3,and we find that b_n=3*b_(n-1)^3, where b_2=2
    {
        ll b = 2, k = 1e4, p = power(13, 8);
        for (ll i = 3; i < k; i++)
        {
            b = 3 * power(b, 3, p); b %= p;
        }
        return power(b, 3, p);
    }
}
namespace Q140
{
    ll solve()//same as finding all a such that 5a^2+14a+1 is square
    {
        ll sum = 0, cnt = 0, d, p;
        for (ll a = 1; cnt < 25; a++)
        {
            p = a * (5 * a + 14) + 1;
            d = sqrt(p);
            if (d * d == p)
            {
                cnt++; sum += a;
            }
        }
        return sum;
    }
}
namespace Q216
{
    ll solve()
    {
        ll n = 1e4, cnt = 0, p, flag;
        vl prime; nsieve(prime, 5e3);
        vvl v;
        for (ll p : prime)
        {
            for (ll i = 1; i <= p/2; i++)
            {
                if ((2 * i * i - 1) % p == 0)
                {
                    v.push_back({ p,i,p - i });
                    break;
                }
            }
        }
        for (ll i = 2; i <= n; i++)
        {
            flag = 0;
            for (ll j = 0; j < v.size(); j++)
            {
                p = v[j][0];
                if (i % p == v[j][1] || i % p == v[j][2])
                {
                    flag = 1; break;
                }
            }
            if (flag)continue;
            if (is_prime(2 * i * i - 1))cnt++;
        }
        return cnt;
    }
}
namespace Q243
{
    ll inc(vl& v, ll k = 10)
    {
        if (v[1] == k)return 0;
        for (ll i = v.size() - 2; ~i; i--)
        {
            if (v[i] < k)
            {
                v[i]++;
                for (ll j = i + 1; j < v.size() - 1; j++)v[j] = v[i];
                return 1;
            }
        }
    }
    ll solve()
    {
        ll n = 1e3,lim= 6469693230,m=1,phi=1;//found from multiplying the first primes
        vl prime;
        nsieve(prime, n);
        for (ll i = 0;; i++)
        {
            ll p = prime[i];
            m *= p; phi *= (p - 1);
            if (94744 * phi < (m - 1) * 15499)return m;
        }
    }
}
namespace Q303
{
    ll solve()
    {
        ll n = 1e4, m, tmp, k;
        vl v(n + 1);
        for (ll i = 0; i < power(3, 14) - 1; i++)
        {
            tmp = i, m = 1, k = 0;
            while (tmp)
            {
                k = tmp % 3 * m; m *= 10;
            }
        }
        return 1;
    }
}
namespace Q304//ive got no idea why it doesnt work
{
    pl fib(ll n,ll p=1e18) {//finds the fibonacci number in logrithmic time
        if (n == 0)return { 0, 1 };
        auto f = fib(n >> 1, p);
        ll c = power(f.first, (2 * f.second - f.first), p, [](ll a, ll b, ll c) {return (a + b) % c; }, 0);
        ll d = power(f.first, f.first, p, [](ll a, ll b, ll c) {return (a + b) % c; }, 0);
        d += power(f.second, f.second, p, [](ll a, ll b, ll c) {return (a + b) % c; }, 0); d %= p;
        if (n & 1)
            return { d, (c + d) % p };
        else
            return { c, d };
    }
    ll solve()
    {
        ll n = 1e14, k = 1e5, i = n + 1, mod = 1234567891011, sum = 0;
        while (i % 6 != 1 && i % 6 != 5)i++;
        vl p;
        while (p.size() < k)
        {
            if (is_prime(i))p.push_back(i);
            if (i % 6 == 1)i += 4;
            else i += 2;
        }
        for (ll x : p)
        {
            sum += fib(x, mod).first;
            sum %= mod;
        }
        return sum;
    }
}
namespace Q329
{
    ll solve()
    {
        ll n = 500, k = 15, a, b, sum = 0, cur, tmp, c, flag;
        string s = "PPPPNNPPPNPPNPN";
        vl div(n + 1), prime; vector<bool>p(n + 1); sieve(prime, div, p, n + 1);
        for (ll i = 1; i <= n; i++)
        {
            cur = i;
            for (ll j = 0; j < (1ll << k); j++)
            {
                c = 1;
                flag = 1;
                tmp = j;
                for (ll m = 0; m < k; m++)
                {
                    if (tmp & 1)
                    {
                        if (cur == n)
                        {
                            flag = 0; break;
                        }
                        if (cur == 1)c <<= 1;
                        cur++;
                        if (s[m] == 'N' && !p[cur] || s[m] == 'P' && p[cur])c <<= 1;
                    }
                    else
                    {
                        if (cur == 1)
                        {
                            flag = 0; break;
                        }
                        if (cur == n)c <<= 1;
                        cur--;
                        if (s[m] == 'N' && !p[cur] || s[m] == 'P' && p[cur])c <<= 1;
                    }
                    if (flag)sum += c;
                }
            }
        }
        return sum;
    }
}
namespace Q366//from another question we know that the minimum to take is the lsb in fib base and so we can calculate the answer from that
{
    ll fbase(ll n,vl &f)
    {
        for (ll i = f.size() - 1; ~i; i--)
        {
            if (n < f[i])continue;
            if (n == f[i])return f[i];
            n -= f[i];
        }
    }
    ll solve()
    {
        ll n = 1e2, sum=0;
        vl f = { 1,2 };
        for (ll i = 2;; i++)
        {
            f.push_back(f[i - 1] + f[i - 2]);
            if (f[i] >= n)break;
        }
        for (ll i = 1; i <= n; i++)
        {
            cout << i << " ";
            for (ll j = (i - 1) / 3; ~j; j--)
            {
                if (j == 0)
                {
                    cout << 0 << '\n';
                }
                else if (fbase(i - j, f) > 2 * j)
                {
                    cout << j << '\n';
                    sum += j;
                    break;
                }
            }
        }
        return sum;
    }
}
namespace Q374
{
    ll solve()
    {
        ll n = 1e2, sum = 1, fact = 2, e = 1, k;//p = 982451653
        for (k = 1; (k * k + 5 * k + 2) / 2 <= n; k++)
        {
            fact *= (k + 2);
            for (ll a = 2; a <= k + 2; a++)
            {
                sum += k * fact / a;
            }
            sum += k * fact * (k + 3) / 2 / (k + 2);
        }
        fact *= (k + 2);
        for (ll a = k + 2; a > (n + 2 - (k * k + 3 * k) / 2); a--)
        {
            sum += k * fact / a;
        }
        return sum;
    }
    ll solve1()
    {
        ll n = 1e1, sum = 1, p = 1e9 + 7, fact = 2, e = 1, k;//p = 982451653
        for (; (e * e + 5 * e + 2) / 2 <= n; e++);
        vl inv(e + 3);
        for (ll i = 1; i <= e + 2; i++)
        {
            inv[i] = inv[i - 1] + power(i, p - 2, p); inv[i] %= p;
        }
        for (k = 1; (k * k + 5 * k + 2) / 2 <= n; k++)
        {
            fact *= (k + 2); fact %= p;
            sum += (k * fact % p) * (inv[k + 2] - 1);sum %= p;
            sum += k * (k + 3) % p * fact % p * (inv[2] - 1) % p * (p + inv[k + 2] - inv[k + 1]); sum %= p;
        }
        fact *= (k + 2); fact %= p;
        sum += k * fact % p * (p + inv[k + 2] - inv[n + 2 - (k * k + 3 * k) / 2]); sum %= p;
        return sum;// -1683550844462 % p;
    }
}
namespace Q378
{
    /*so you can find the number of inversions in nlogn
    so you find for every index the number of invesrions, and for every index j the number of i<j such that A[i]>A[j] 
    (which you can find by running inverison algo on reverse array with minus values).
    after that you run over j and find the wanted i,k (which are independent of each other, so it runs linearly*/
}
namespace Q386
{
    ll solve()
    {
        ll n = 1e6, k = log2(n), sum = 0;
        vvl ch(k + 2);
        vl div(n + 1), prime;
        vector<bool>p(n + 1); sieve(prime, div, p, n + 1);
        vl ds(n + 1);
        for (ll i = 0; i <= k+1; i++)
        {
            ch[i].resize(i + 1, 1);
            for (ll j = 1; j < i; j++)ch[i][j] = (ch[i - 1][j - 1] + ch[i - 1][j]);
        }
        for (ll i = 2; i <= n; i++)ds[i] = ds[i / div[i]] + 1;
        for (ll i = 1; i <= n; i++)sum += ch[ds[i]][ds[i] / 2];
        return sum;
    }
}
using namespace Q237;
int main()
{
    auto start = chrono::high_resolution_clock::now();
    cout << solve() << "\n\n";
    //cout << setprecision(11) << solve() << "\n\n";
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(end - start).count();
    cout << "Time taken: " << duration << " microseconds" << std::endl;
}

