import numpy as np
Scal = 128
Count1=0
Count2=0
global Num, Flen, Srange, Startq, CurRow, Count
global Mat, Fbas, Rvec, Lvec, Sieve, Pivot, Uvec, Vvec, StackUV
def trial_factor(n):
    factors = []
    d = 2  # Start with the smallest prime factor

    while d * d <= n:
        if n % d == 0:
            factors.append(d)
            n //= d
        else:
            d += 1

    if n > 1:
        factors.append(n)

    return factors
def test_trial_factor():
    n = 274177 * 67280421310721
    print("n = ", n, "=", 274177, "*", 67280421310721)
    factors = trial_factor(n)
    print("factors = ", factors)
test_trial_factor()
def fermat_factor(N):
    """
    Benutzt das Fermat'sches Faktorisierungs-Verfahren um positive Integer zu faktorisieren

    Parameter:
    - N: Der zu faktorisierende Integerwert

    Rückgabe:
    - ein Tupel (p, q) der Faktoren von N oder None bei negativen/ nicht Integer Werten
    """
    # Prüft den Faktoren 2 und 5
    if ( not N%2 ):
        return (2, N/2)
    if mod(N,5)==0:
        return (5,N/5)

    
    # Try to find a nontrivial square root of n
    x0 = ceil(sqrt(N))
    y_square = x0**2 - N

    while not is_square(y_square):
        x0 += 1
        y_square = x0**2 - N
        if x0>=N :
            return (1,N)
    
    # Berechnet die Zerlegung in (x-y)(x+y)
    p = x0 - sqrt(y_square)
    q = x0 + sqrt(y_square)
    return (p, q)

def test_fermat_factor():
    n = 274177 * 274179
    print("n = ", n, "=", 274177, "*", 274179)
    factors = fermat_factor(n)
    print("factors = ", factors)
test_fermat_factor()
def make_factorbase(N, len):
    """
    Konstruiert eine Faktorbasis zur gegebenen Zahl N und Länge len

    Parameter:
    - N: Der zu faktorisierende Integerwert
    - len: Die gefordete Länger der Faktorbasis

    Rückgabe:
    - Eine Liste der Länge len von Primzahlen, sodass qudratische Reste mit N entstehen 
    """
    fb = [-1,2]
    p = 2
    for k in range(len-2):
        p = next_prime(p)
        while legendre_symbol(N, p) != 1:
            p = next_prime(p)
        fb.append(p)
    return fb

def test_make_factorbase() :
    N = 274177 * 67280421310721
    fb = make_factorbase(N, 10)
    print("fb = ", fb)
test_make_factorbase()
def rtpvec(N, pvec):
    """
    Berechnet die Wurzel von N im p-Köper der p aus der Faktorbasis 

    Parameter:
    - N: Integer dessen Wurzel berechnet wird
    - pvec: Liste der Zahlen, die die Faktorbasis darstellen

    Rückgabe:
    - Eine Liste der Wurzel von N auf Basis der Faktorbasis
    """
    rvec = [mod(N,p).sqrt() for p in pvec]
    rvec[0]=rvec[1]=1
    return rvec

def test_rtpvec():
    N = 274177 * 67280421310721
    fb = make_factorbase(N, 10)
    rvec = rtpvec(N, fb)
    print("rvec = ", rvec)
test_rtpvec()
class QuadPol:
    """
    Klasse zum speichern der Werte Parameter des Polynoms ax^2 + bx + c.
    """
    def __init__(self, a, b, c, q):
        self.a = a
        self.b = b
        self.c = c
        self.q = q
    
def Zp2_sqrt(p, n):
    """
    Berechnet die Wurzel von N modulo p^2 nach: 
    Satz:   Für ungerade Primzahl p, n>=2 und a mit p teilt nicht a quadratischer Rest mod p
            => Es gibt  eindeutige Zahl x mit
        x^2 ≡ a mod p^n UND x ≡ x0 mod p.

    Parameter:
    - p: Primzahl
    - n: Integer

    Rückgabe:
    - Die Wurzel von n modulo p^2
    """

    x0 = mod(n, p).sqrt()
    xi = inverse_mod(2*x0, p**2)
    x = (x0**2 + n)*xi
    return mod(x, p**2)
def make_quadpol(N, q):
    """
    Konstruiert zur Faktorisierung von N ein quadratisches Polynom ax^2 + bx + c 

    Parameter:
    - N: Der zu faktorisierende Integerwert
    - q: Ausgangswert der Suche nach der kleinsten Primzahl mit (N/q) = 1

    Rückgabe:
    - Ein QuadPol Objekt welches das Polynom ax^2 + bx + c repräsentiert
    """
    q = next_prime((q-1))
    while kronecker_symbol(N, q) != 1:
        q = next_prime(q+1)
    a = q**2
    b = mod(N,a).sqrt()
    c = (Integer(b)**2 - N) // a
    return QuadPol(a,b,c,q)

def test_make_quadpol():
    N = 445546739
    q = 2
    qp = make_quadpol(N, q)
    print("a = ", qp.a)
    print("b = ", qp.b)
    print("c = ", qp.c)
    print("q = ", qp.q)
test_make_quadpol()

def logpvec(pvec):
    """
    Berechnet den Logarithmus der p aus der Faktorbasis basierend auf der Skalierung scal

    Parameter:
    - pvec: Liste der Zahlen, die die Faktorbasis darstellen

    Rückgabe:
    - Eine Liste der Logarithmen von N auf Basis der Faktorbasis als Ganzzahlige Werte
    """
    global Scal

    lvec = [0] * len(pvec)
    
    for k in range(1, len(pvec)):
        lvec[k] = round(log(pvec[k]) * Scal)
    
    return lvec
def test_logpvec():
    N = 274177 * 67280421310721
    fb = make_factorbase(N, 10)
    lvec = logpvec(fb)
    print("lvec = ", lvec)
test_logpvec()
def QSinitialize(N):
    """
    Die Funktion QSfactorize(N) initialisiert zuerst die erforderlichen Variablen:

    Parameter:
    - N: Die zu faktorisierende Zahl N.

    Rückgabe:
    - Flen: Die Länge der verwendeten Faktorbasis.
    """

    global Num, Flen, Srange, Startq, CurRow, Count
    global Mat, Fbas, Rvec, Lvec, Sieve, Pivot, Uvec, Vvec, StackUV
    
    Num = N
    blen = int(N).bit_length()
    Flen = max(8, blen**2 // 32)
    Srange = min(blen * 256, (sys.maxsize - 1) // 2)
    Startq = isqrt(isqrt(2 * N) // Srange)
    Fbas = make_factorbase(N, Flen)
    Rvec = rtpvec(N, Fbas)
    Lvec = logpvec(Fbas)
    Mat = MatrixSpace(GF(2),Flen+1,2*Flen +1,sparse=True).matrix()
    Pivot = [i for i in range(Flen)]
    CurRow = 0
    Count = 0
    Vvec = Uvec = [0] * (Flen + 1)
    Sieve = [0] * (2 * Srange)
    StackUV = []
    
    return Flen

def QSfactorize(N):
    """
    Die Funktion QSfactorize(N) implementiert den Faktorisierungsprozess des Quadratischen Siebs (Quadratic Sieve) für eine gegebene Zahl N.

    Parameter:
    - N: Die zu faktorisierende Zahl N.

    Rückgabe:
    - none
    """    
    
    global Srange, Fbas, CurRow, Uvec, Vvec, Flen, Count1, Count2

    QSinitialize(N)
    print("quadratic sieve length", 2 * Srange,
          ", factorbase 2 ...", Fbas[Flen - 1], "of length", Flen - 1)
    print('working ', end='')
    nbfail = 0
    while nbfail < 32:
        (u, v) = get_qres(N)
        if not qr_trialdiv(v):
            continue
        Uvec[CurRow] = u
        Vvec[CurRow] = v
        if gausselim(CurRow) == 0:
            CurRow += 1
            print('.', end='')
        else:
            relvec = getrel(CurRow)
            print('!', end='')
            d = findfactor(N, relvec)
            if d > 0:
                print()
                print(Count1, "polynomials,", Count2,
                      "completely factorized quadratic residues")
                return d
            else:
                nbfail += 1
    print('\nfailed to often, aborting.... \n')
def get_qres(N):
    """
    Die Funktion get_qres(N) hat den Zweck, Paare von u und v zu generieren, für die u^2 ≡ v (mod N)

    Parameter:
    - N: Die zu faktorisierende Zahl N.

    Rückgabe:
    - UV (Tupel): Ein StackUV Element mit einem Paar (u, v).
    """
    global Startq, Count1, StackUV
    Q = None
    UV = []
    while not StackUV:
        Q = make_quadpol(N, Startq)
        Count1 += 1
        print('_', end='')
        dosieve(Q)
        sieveresults(Q)
        Startq = Q.q + 2
    UV = (StackUV.pop())
    return UV

def dosieve(Q):
    """ 
    dosieve implementiert das Sieben des Siebarrays für eine gegebene quadratische Polynomfunktion Q, indem sie die Siebwerte entsprechend der Faktorbasis und den Wurzeln von Q aktualisiert.

    Parameter:
    - Q (quadpol): Ein quadratisches Polynom, das verarbeitet wird.

    Rückgabe:
    - none
    """
    global Sieve, Srange, Fbas, Lvec, Rvec, Flen

    for i in range(2 * Srange):
        Sieve[i] = 0
    
    for k in range(2, Flen):
        p = Fbas[k]
        z = Lvec[k]
        r = Rvec[k]
        s = (-Srange) % p
        
        if Q.a % p != 0:
            a1 = inverse_mod(Q.a, p)
            b1 = mod(Q.b,p)
            r1 = (r - b1) * a1 % p
            if r1 >= s:
                i0 = r1 - s
            else:
                i0 = p + r1 - s
            
            for i in range(i0, 2 * Srange, p):
                Sieve[i] += z
            
            r2 = (p - r - b1) * a1 % p
            if r2 >= s:
                i0 = r2 - s
            else:
                i0 = p + r2 - s
            
            for i in range(i0, 2 * Srange, p):
                Sieve[i] += z

def sieveresults(Q):
    """
    Die Funktion "sieveresults" identifiziert potenzielle Faktorisierungskandidaten, indem sie die x-Werte überprüft, bei denen die Summe der Logarithmen einen bestimmten Schwellenwert überschreitet, und entsprechende (u, v)-Paare auf den Stack legt.

    Parameter:
    - Q (quadpol): Ein quadratisches Polynom, das verarbeitet wird.

    Rückgabe:
    - len(StackUV): aktuelle Anzahl der Elemente auf der Stack-Liste
    """
    global Sieve, StackUV, Flen, Lvec, Num, Scal, Srange

    target = round(log(Num/Q.a)*Scal) - Lvec[Flen-1]
    qinv = int(inverse_mod(Q.q, Num))
    for k in range(2*Srange):
        if Sieve[k] >= target:
            x = k - Srange
            u = Q.a * x + Q.b
            v = (u + Q.b) * x + Q.c
            u = int(mod(qinv * u , Num))
            StackUV.append((u, v))
    return len(StackUV)

def qr_trialdiv(v):
    """
    Die Funktion qr_trialdiv überprüft, ob eine gegebene Zahl v vollständig faktorisiert werden kann.

    Parameter:
    - v: die zu überprüfende Zahl

    Rückgabe:
    - true/false ob zutreffend

    """

    global Mat, Fbas, Flen, CurRow, Count2

    Mat[CurRow] += Mat[CurRow]
    
    if v < 0:
        v = abs(v)
        Mat[CurRow,0] = 1
    
    for i in range(1, Flen):
        p = Fbas[i]

        
        while (v > 1) and (mod(v,p) == 0):
            v = int(v) // int(p)
            Mat[CurRow,i] += 1    
    
        if v <= 1:
            Mat[CurRow, Flen + CurRow] = 1
            Count2 += 1
            return True
    
    return False

def gausselim(k):
    """
    Führt die Gauß-Elimination auf der Matrix Mat durch, um sie in eine obere Dreiecksform zu bringen.

    Parameter:
    - k: Zeile in der Matrix, ab der die Zeilenstufenform noch erreicht werden muss

    Rückgabe:
    - 1/0 ob erfolgreich
    """
    global Mat, Pivot

    for i in range(k):
        if Mat[k][Pivot[i]]:
            Mat[k] = Mat[k] + Mat[i]
    
    for i in range(k, len(Pivot)):
        j = Pivot[i]
        if Mat[k][j]:
            Pivot[i] = Pivot[k]
            Pivot[k] = j
            return 0
    
    return 1

def getrel(row):
    """
    Diese Funktion durchläuft die Bits in der Zeile Mat[row] und sucht nach gesetzten Bits ab dem Index Flen. 
    """

    global Mat, Flen

    st = []
    for i in range(row+1):
        if Mat[row][Flen+i]:
            st.append(i)
    return st

def findfactor(N, relvec):
    """
    Die Funktion findfactor sucht einen Faktor von N basierend auf den gegebenen Informationen in relvec.
    """
    global Uvec, VvecMat

    k = relvec[0]
    x = Integer(Uvec[k])
    v = Integer(Vvec[k])
    y = 1

    for i in range(1, len(relvec)):
        k = Integer(relvec[i])
        u = Integer(Uvec[k])
        v1 = Integer(Vvec[k])
        x = (x * u) % N
        d = gcd(v, v1)
        v = (v // d) * (v1 // d)
        y = (y * d) % N

    y = (y * isqrt(v)) % N
    d = gcd(x + y, N)
    
    if d <= 1 or d == N:
        d = 0
    
    return d
