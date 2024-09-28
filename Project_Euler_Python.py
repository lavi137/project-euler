from ast import Continue
from base64 import b16decode
from http.client import CONTINUE
from math import sqrt
from math import log10
from decimal import Decimal, getcontext
def sieve(prime,div,p,n):
    for i in range(2,n):
        if(div[i]==0):
            prime.append(i)
            div[i]=i
            p[i]=1
        j=0
        while(j<len(prime) and prime[j]<=div[i] and prime[j]*i<n):
            div[prime[j]*i]=prime[j]
            j+=1
    return 1
def factorial(a):
    if(a==0):return 1
    ans=1
    for b in range(a,0,-1):ans*=b
    return ans
def choose(a,b):
    return (factorial(a)//factorial(b))//factorial(a-b)
def gcd(a,b):
    if(b==0):
        return a
    return gcd(b,a%b)
def Q16():
    n=1
    s=0
    for i in range(0,1000):
        n*=2
    while(n>0):
        s+=n%10
        n=n//10
    return s
def Q20():
    n=1
    s=0
    for i in range(1,101):
        n*=i
    while(n>0):
        s+=n%10
        n=n//10
    return s
def Q44():
    i=0
    while 1:
        i+=1
        a=(3*i-1)*i/2
        for j in range(1,i):
            b=(3*j-1)*j/2
            if(sqrt(1+24*a+48*b)%6==5 and sqrt(1+24*a+24*b)%6==5):
                return a
                break;
def Q53():
    s=0
    for i in range(0,101):
        for j in range(0,i+1):
            if(choose(i,j)>1e6):
                s+=1
    return s
def Q55():
    cnt=0
    for i in range(10001):
        flag=0
        tmp=i
        for j in range(50):
            tmp=tmp+int(str(tmp)[::-1])
            if(tmp==int(str(tmp)[::-1])):
                flag=1
                j=50
        if(not flag):
            cnt=cnt+1
    return cnt         
def Q56():
    sum=-1;
    for a in range(101):
        for b in range(101):
            cur=0;
            c=a**b;
            while(c!=0):
                cur+=(c%10)
                c=c//10
            sum=max(sum,cur)
    return sum;
def Q57():
    a=1
    b=1
    cnt=0
    for i in range(1000):
        a+=2*b
        b=a-b
        if(int(log10(a))!=int(log10(b))):cnt+=1
    return cnt        
def Q63():
    i=0
    cnt=0
    while 1:
        i+=1
        if(9**i<10**(i-1)):break
        for j in range(1,10):
            if(int(log10(j**i))==i-1):cnt+=1
    return cnt
def Q80():
    sum=0
    for i in range(2,101):
        getcontext().prec=200
        x=Decimal(i).sqrt()
        if(int(x)==x):continue
        while(x>=1):x=x/10
        for j in range(1,101):
            x*=10
            x-=10*(x//10)
            sum+=int(x)
    return sum
def Q190():#by induction on x_m we find the recursion of v
    v=[1];
    s=0
    for i in range(2,16):
        a=(i-1)/2
        v.append(v[i-2]*pow(pow(a,a)/pow(a+1,a+1),2*a+1))
    for i in range(1,15):
        s+=int(v[i]*pow((i+1),(i+1)*(i+2)/2))
    return s
def Q267():#for const k its just finding max for p(x)=(1+2x)^k*(1-x)^(1000-k)
    a=0
    for k in range(334,1000):
        x=(3*k-1000)/2000
        if(k*log10(1+2*x)+(1000-k)*log10(1-x)>9):
            break
    for i in range(k,1001):
        a+=choose(1000,i)
    return round(a/(2**1000),12)
def Q493():
    a=[0,0,0,0,0,0,0,0]
    b=0
    c=0
    for i in range(2,8):
        a[i]=choose(7,i)*choose(10*i,20)
        for j in range(i-1,1,-1):
            a[i]-=a[j]
        b+=i*a[i]
    return round(b/choose(70,20),9)
def Q686():
    cnt=0
    cur=0
    i=10
    y=2**10
    l=4
    while 1:
        if(cnt==678910):return cur
        i+=1
        if((y//(10**(l-1)))>4):l+=1
        y*=2
        x=y
        x=x//(10**(l-3))
        if(x==123):
            cnt+=1
            cur=i
            if(cnt%100==0):print(cnt)
    return 1
def Q853():#idk why dont work
    f=[0]*124
    f[1]=1
    s=0
    for i in range(2,124):
        f[i]=f[i-1]+f[i-2]
    a=gcd(f[120]-f[0],f[121]-f[1])
    b=set()
    for i in range(1,int(sqrt(a))):
        if(a%i==0):
            b.add(i)
            b.add(a//i)
    c=[60,40,24]
    for i in range(0,3):
        d=gcd(f[c[i]]-f[0],f[c[i]+1]-f[1])
        for j in range (1,int(sqrt(d))):
            if(d%j==0):
                b.discard(j)
                b.discard(d//j)
    for x in b:
        if(x<1e9):
            s+=x
    return s;
print(Q190())

            

