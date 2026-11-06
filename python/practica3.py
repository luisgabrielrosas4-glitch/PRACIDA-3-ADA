import sys
import time
from array import array

def memoria():
    if sys.platform == 'win32':
        import ctypes
        class PM(ctypes.Structure):
            _fields_ = [('cb',ctypes.c_ulong),('PageFaultCount',ctypes.c_ulong),('PeakWorkingSetSize',ctypes.c_size_t),('WorkingSetSize',ctypes.c_size_t),('QuotaPeakPagedPoolUsage',ctypes.c_size_t),('QuotaPagedPoolUsage',ctypes.c_size_t),('QuotaPeakNonPagedPoolUsage',ctypes.c_size_t),('QuotaNonPagedPoolUsage',ctypes.c_size_t),('PagefileUsage',ctypes.c_size_t),('PeakPagefileUsage',ctypes.c_size_t)]
        x=PM()
        x.cb=ctypes.sizeof(x)
        return x.PeakWorkingSetSize/1048576 if ctypes.windll.psapi.GetProcessMemoryInfo(ctypes.windll.kernel32.GetCurrentProcess(),ctypes.byref(x),x.cb) else -1
    import resource
    v=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return v/1048576 if sys.platform=='darwin' else v/1024

def datos(n,bits):
    a=array('i' if bits==32 else 'q',(1+(i*7919+1234567)%10000000 for i in range(n)))
    seed=734251
    for i in range(n-1,0,-1):
        seed^=(seed<<13)&0xffffffff
        seed^=seed>>17
        seed^=(seed<<5)&0xffffffff
        seed&=0xffffffff
        j=seed%(i+1)
        a[i],a[j]=a[j],a[i]
    return a

def quickselect(a,k):
    lo,hi=0,len(a)
    k-=1
    while True:
        if hi-lo==1:
            return a[lo]
        pivot=a[lo+(hi-lo)//2]
        lt,i,gt=lo,lo,hi
        while i<gt:
            if a[i]<pivot:
                a[lt],a[i]=a[i],a[lt]
                lt+=1
                i+=1
            elif a[i]>pivot:
                gt-=1
                a[i],a[gt]=a[gt],a[i]
            else:
                i+=1
        if k<lt:
            hi=lt
        elif k>=gt:
            lo=gt
        else:
            return a[k]

def suma(a,b):
    n=len(a)
    return [[a[i][j]+b[i][j] for j in range(n)] for i in range(n)]

def resta(a,b):
    n=len(a)
    return [[a[i][j]-b[i][j] for j in range(n)] for i in range(n)]

def cortar(a,fi,fj):
    h=len(a)//2
    return [a[i+fi*h][fj*h:fj*h+h] for i in range(h)]

def tradicional_2d(a,b):
    n=len(a)
    c=[[0]*n for _ in range(n)]
    for i in range(n):
        ci=c[i]
        ai=a[i]
        for k in range(n):
            v=ai[k]
            bk=b[k]
            for j in range(n):
                ci[j]+=v*bk[j]
    return c

def strassen_2d(a,b):
    n=len(a)
    if n<=64:
        return tradicional_2d(a,b)
    h=n//2
    a11,a12,a21,a22=(cortar(a,x,y) for x,y in ((0,0),(0,1),(1,0),(1,1)))
    b11,b12,b21,b22=(cortar(b,x,y) for x,y in ((0,0),(0,1),(1,0),(1,1)))
    p1=strassen_2d(suma(a11,a22),suma(b11,b22))
    p2=strassen_2d(suma(a21,a22),b11)
    p3=strassen_2d(a11,resta(b12,b22))
    p4=strassen_2d(a22,resta(b21,b11))
    p5=strassen_2d(suma(a11,a12),b22)
    p6=strassen_2d(resta(a21,a11),suma(b11,b12))
    p7=strassen_2d(resta(a12,a22),suma(b21,b22))
    q11=suma(resta(suma(p1,p4),p5),p7)
    q12=suma(p3,p5)
    q21=suma(p2,p4)
    q22=suma(resta(suma(p1,p3),p2),p6)
    return [q11[i]+q12[i] for i in range(h)]+[q21[i]+q22[i] for i in range(h)]

def sumar_flat(a,b):
    return array('q',(x+y for x,y in zip(a,b)))

def restar_flat(a,b):
    return array('q',(x-y for x,y in zip(a,b)))

def cortar_flat(a,n,fi,fj):
    h=n//2
    c=array('q')
    for i in range(h):
        s=(fi*h+i)*n+fj*h
        c.extend(a[s:s+h])
    return c

def tradicional_flat(a,b,n):
    c=array('q',[0])*(n*n)
    for i in range(n):
        base=i*n
        for k in range(n):
            v=a[base+k]
            bk=k*n
            for j in range(n):
                c[base+j]+=v*b[bk+j]
    return c

def strassen_flat(a,b,n):
    if n<=64:
        return tradicional_flat(a,b,n)
    h=n//2
    a11,a12,a21,a22=(cortar_flat(a,n,x,y) for x,y in ((0,0),(0,1),(1,0),(1,1)))
    b11,b12,b21,b22=(cortar_flat(b,n,x,y) for x,y in ((0,0),(0,1),(1,0),(1,1)))
    p1=strassen_flat(sumar_flat(a11,a22),sumar_flat(b11,b22),h)
    p2=strassen_flat(sumar_flat(a21,a22),b11,h)
    p3=strassen_flat(a11,restar_flat(b12,b22),h)
    p4=strassen_flat(a22,restar_flat(b21,b11),h)
    p5=strassen_flat(sumar_flat(a11,a12),b22,h)
    p6=strassen_flat(restar_flat(a21,a11),sumar_flat(b11,b12),h)
    p7=strassen_flat(restar_flat(a12,a22),sumar_flat(b21,b22),h)
    q11=sumar_flat(restar_flat(sumar_flat(p1,p4),p5),p7)
    q12=sumar_flat(p3,p5)
    q21=sumar_flat(p2,p4)
    q22=sumar_flat(restar_flat(sumar_flat(p1,p3),p2),p6)
    c=array('q',[0])*(n*n)
    for i in range(h):
        for j in range(h):
            v=i*h+j
            c[i*n+j]=q11[v]
            c[i*n+j+h]=q12[v]
            c[(i+h)*n+j]=q21[v]
            c[(i+h)*n+j+h]=q22[v]
    return c

def correr(args):
    modo=args[0]
    if modo=='select':
        bits=int(args[1]);n=int(args[2]);k=int(args[3]) if len(args)>3 else (n+1)//2
        if bits not in (32,64) or not 1<=k<=n<=10000000:
            raise ValueError('Parametros invalidos')
        a=datos(n,bits)
        s=time.perf_counter()
        v=quickselect(a,k)
        secs=time.perf_counter()-s
        menor=sum(x<v for x in a)
        igual=sum(x==v for x in a)
        ok='si' if menor<k<=menor+igual else 'no'
        print(f'python,quickselect,{bits},-,{n},{k},{secs:.6f},{memoria():.6f},{ok},{v}')
    elif modo=='matrix':
        alg,rep,n=args[1],args[2],int(args[3])
        if alg not in ('strassen','tradicional') or rep not in ('2d','flat') or n not in (16,32,64,128,256,512):
            raise ValueError('Parametros invalidos')
        if rep=='2d':
            a=[[( (i*n+j)*17+3)%11-5 for j in range(n)] for i in range(n)]
            b=[[( (i*n+j)*37+9)%11-5 for j in range(n)] for i in range(n)]
            s=time.perf_counter()
            c=strassen_2d(a,b) if alg=='strassen' else tradicional_2d(a,b)
            secs=time.perf_counter()-s
            ok=all(c[i][j]==sum(a[i][k]*b[k][j] for k in range(n)) for i in (0,n//2,n-1) for j in (0,n//2,n-1))
            checksum=sum(map(sum,c))
        else:
            a=array('q',(((i*17+3)%11)-5 for i in range(n*n)))
            b=array('q',(((i*37+9)%11)-5 for i in range(n*n)))
            s=time.perf_counter()
            c=strassen_flat(a,b,n) if alg=='strassen' else tradicional_flat(a,b,n)
            secs=time.perf_counter()-s
            ok=all(c[i*n+j]==sum(a[i*n+k]*b[k*n+j] for k in range(n)) for i in (0,n//2,n-1) for j in (0,n//2,n-1))
            checksum=sum(c)
        print(f'python,{alg},64,{rep},{n},-,{secs:.6f},{memoria():.6f},{"si" if ok else "no"},{checksum}')
    else:
        raise ValueError('Modo desconocido')

if __name__=='__main__':
    correr(sys.argv[1:])
