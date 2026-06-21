# Esercizio Completo 1 — Floating Point (versione integrata)

**a)** Calcola in MATLAB `realmax * 10` e `realmin / 100`. Commenta i fenomeni osservati.

**b)** Calcola `1 + eps` e confrontalo con `1` usando l'operatore `==`. Commenta.

**c)** Date le espressioni g(s) = s − √(s²−1) e p(s) = 1/(s+√(s²−1)), calcola entrambe in MATLAB per s = 8.888×10⁸. Commenta quale soffre di cancellazione numerica e perché.

**d)** Calcola analiticamente (a+b)+c e a+(b+c) con a = 0.2×10⁻⁵, b = 0.45×10³, c = −0.45×10³. Poi verifica in MATLAB. Commenta il risultato in relazione alla proprietà associativa.

**e)** Calcola le radici dell'equazione x² − 0.6001x + 9×10⁻⁵ = 0 sia analiticamente con la formula quadratica standard sia con MATLAB. Una delle due radici soffre di cancellazione: identificala e proponi una formula alternativa numericamente stabile.

**f)** Siano a = 2.5×10¹⁸⁰, b = 4×10¹⁵⁰, c = 1×10⁻¹⁶⁰. Calcola (a·b)·c e a·(b·c) in MATLAB. I due risultati coincidono? Commenta in relazione a `realmax` spiegando da cosa dipende l'eventuale differenza.

**g)** Tramite un ciclo `for`, calcola la somma Σ (per k=1 a 3000) di 1/1000. Confronta il risultato con il valore esatto (che dovrebbe essere 3) e commenta l'origine dell'errore osservato.

---

# Esercizio Completo 2 — Sistemi Lineari (versione integrata)

Dato il sistema lineare Ax=b con:

```
A = | 4   1   0 |      b = | 5 |
    | -1  3  0.3|          | 6 |
    | 1   β   3 |          | 7 |
```

**a)** Scegliere un valore di β per cui il sistema abbia un'unica soluzione e il metodo di Jacobi converga alla soluzione esatta. Giustificare la scelta sul foglio.

**b)** Descrivere sul foglio l'algoritmo di Jacobi e scriverlo per componenti per questo sistema specifico (n=3).

**c)** Eseguire analiticamente (a mano) un passo del metodo di Jacobi partendo da x⁽⁰⁾=[1,1,1]ᵗ.

**d)** Implementare in MATLAB il metodo di Jacobi e approssimare la soluzione, fissata una tolleranza 10⁻⁴.

**e)** Calcolare il raggio spettrale della matrice di iterazione di Jacobi e commentare riguardo alla convergenza.

**f)** Risolvere lo stesso sistema con l'operatore `\` di MATLAB e confrontare con la soluzione di Jacobi calcolando l'errore.

**g)** Calcolare il numero di condizionamento di A in norma infinito. Il sistema è ben condizionato?

**h)** Implementare anche il metodo di Gauss-Seidel sullo stesso sistema e confrontare il numero di iterazioni necessarie per raggiungere la stessa tolleranza rispetto a Jacobi. Rappresentare graficamente l'andamento dell'errore per entrambi i metodi in funzione delle iterazioni.

**i)** Trascrivere in MATLAB il seguente codice di Gauss-Seidel e farlo girare sullo stesso sistema, verificando che dia lo stesso risultato del punto h):

```matlab
function [x, B_Gauss] = GaussS(a, b, x0, nmax)
    [n, n] = size(a);
    iter = 0;
    x = [];
    x_old = x0;
    while iter < nmax
        iter = iter + 1;
        D = diag(diag(a));
        E = tril(a, -1);
        F = triu(a, 1);
        B_Gauss = -inv(D+E)*F;
        x(:,iter) = -inv(D+E)*F*x_old + inv(D+E)*b;
        x_old = x(:,iter);
    end
end
```

**j)** Costruire, al variare della dimensione di una matrice tridiagonale A∈ℝⁿˣⁿ con 2 sulla diagonale e −1 sulle codiagonali (n=10,20,...,100), una tabella che riporti: errore commesso risolvendo con `\`, errore con sostituzione diretta (se A è triangolarizzata via LU), tempo computazionale (con `tic`/`toc`), e numero di condizionamento. Il termine noto b va costruito in modo che la soluzione esatta sia un vettore di tutti 1. Commentare l'andamento al crescere di n.

**k)** Ripetere il punto j) usando la matrice di Hilbert `hilb(n)` invece della tridiagonale, per n=2,...,15. Confrontare il condizionamento delle due matrici e commentare quale sistema è più difficile da risolvere numericamente.

---

# Esercizio Completo 3 — Radici di Equazioni Non Lineari (versione integrata)

Data la funzione f(s) = s³ − 2s − 5 nell'intervallo (1,3).

**a)** Descrivere sul foglio il problema della ricerca di radici di equazioni non lineari, sia dal punto di vista analitico che numerico.

**b)** Fare il grafico di f(s) nell'intervallo (1,3) in MATLAB.

**c)** Descrivere sul foglio l'algoritmo di bisezione e dimostrare che è globalmente convergente.

**d)** Implementare la bisezione in MATLAB con criterio d'arresto basato sulla tolleranza sull'ampiezza dell'intervallo, tolleranza 10⁻⁸. Quante iterazioni servono? Verificarlo anche con la formula a priori.

**e)** Descrivere sul foglio il metodo delle corde e il criterio d'arresto basato sul controllo del residuo.

**f)** Eseguire analiticamente (a mano) tre passi del metodo delle corde nell'intervallo (1,3) partendo da s⁽⁰⁾=1.

**g)** Implementare il metodo di Newton in MATLAB con punto di innesco s₀=2 e criterio d'arresto sull'incremento, tolleranza 10⁻⁸. Costruire una tabella con iterazione, approssimante e residuo.

**h)** Implementare anche il metodo delle secanti partendo da x₋₁=1, x₀=3.

**i)** Confrontare il numero di iterazioni di bisezione, corde, Newton e secanti per raggiungere la stessa precisione. Rappresentare sullo stesso grafico la convergenza dei quattro metodi (errore vs iterazioni, scala logaritmica).

**j)** Calcolare la radice anche con `fzero` e verificare la coerenza di tutti i risultati.

**k)** Trascrivere in MATLAB il seguente codice di bisezione dato dalla prof e farlo girare sulla stessa funzione, verificando che il risultato coincida con quello del punto d):

```matlab
function [xvect, nit, xdif] = bisect(a, b, toll, nmax, fun)
    err = abs(b-a);
    nit = 1;
    xvect = b;
    fx = [];
    xdif = [];
    while nit < nmax && err > toll
        nit = nit + 1;
        c = (a+b)/2;
        x = c;
        fc = fun(x);
        xvect = [xvect; x];
        fx = [fx; fc];
        x = a;
        if fc*fun(x) > 0
            a = c;
        else
            b = c;
        end
        err = abs(xvect(nit) - xvect(nit-1));
        xdif = [xdif; err];
    end
end
```

**l)** Ora si consideri la funzione g(x) = (x−2)² nell'intervallo (0,4), che ha una radice doppia in x=2. Applicare il metodo di Newton con punto di innesco x₀=0 e osservare quante iterazioni servono per raggiungere tolleranza 10⁻⁸ sull'incremento, confrontandolo con il numero di iterazioni necessario per la radice semplice del punto g). Applicare anche la bisezione su (0,4): funziona? Commentare perché.

---

# Esercizio Completo 4 — Integrazione Numerica (versione integrata)

Si consideri l'integrale:

I = ∫₀² (1 + x·e^(−x)) dx

**a)** Descrivere sul foglio la formula del punto medio semplice e il suo grado di esattezza.

**b)** Descrivere sul foglio la formula del trapezio composita su una decomposizione uniforme dell'intervallo.

**c)** Fare il grafico della funzione integranda nell'intervallo [0,2].

**d)** Implementare in MATLAB la formula del punto medio semplice per approssimare I.

**e)** Implementare in MATLAB la formula del trapezio composita per approssimare I, con m sottointervalli.

**f)** Approssimare I anche con `integral` (o `quad`) e usarlo come valore di riferimento.

**g)** Costruire una tabella con il decadimento dell'errore della formula del trapezio composita all'aumentare di m = 2,4,8,16,32,64,128 rispetto al riferimento ottenuto al punto f).

**h)** Verificare la correttezza del trapezio composito confrontando il risultato anche con il comando `trapz`.

**i)** Usare `tic`/`toc` per misurare il tempo computazionale al variare di m e commentare come cresce.

**j)** Implementare anche la formula di Cavalieri-Simpson composita e confrontare il decadimento dell'errore con quello del trapezio, a parità di m. Quale dei due algoritmi risulta più accurato?

**k)** Sapendo che max|f''(x)| sull'intervallo [0,2] vale circa 0.74, calcolare analiticamente quanti sottointervalli m servono per garantire un errore inferiore a 10⁻⁶ con la formula del trapezio composita, usando la stima a priori dell'errore. Verificare poi numericamente in MATLAB se il valore di m trovato rispetta effettivamente la tolleranza richiesta.

---

# Esercizio Completo 5 — Interpolazione e Minimi Quadrati (versione integrata)

Dati i punti di coordinate xᵢ, yᵢ per i=1,...,5:

```
xᵢ  | -2   -1   0    1    2
yᵢ  | 4.2  1.1  0.3  1.4  4.0
```

**a)** Scrivere sul foglio la definizione di polinomio interpolante n+1 punti assegnati. Esiste sempre? È unico?

**b)** Calcolare analiticamente (a mano) il polinomio interpolante di grado 2 usando solo i punti (−1,1.1), (0,0.3), (1,1.4), tramite la forma di Lagrange.

**c)** In MATLAB, costruire il polinomio interpolante di grado 4 di tutti i 5 punti tramite la matrice di Vandermonde.

**d)** Fare il grafico nello stesso sistema di riferimento della funzione interpolante di grado 4, dei 5 punti dati, e del polinomio di regressione lineare (ai minimi quadrati) calcolato sugli stessi punti.

**e)** Sul foglio, impostare i passaggi per calcolare i coefficienti m e a della retta di regressione lineare per questi 5 punti.

**f)** Disegnare un'altra retta a piacere (diversa da quella ottimale) e verificare matematicamente (calcolando la somma dei residui quadratici) che non soddisfa la definizione di retta ai minimi quadrati.

**g)** Ripetere l'interpolazione di grado 4 usando nodi di Chebyshev sull'intervallo [−2,2] invece dei nodi xᵢ originali (supponendo che i valori yᵢ provengano da una funzione nota, es. g(x)=x²−2cos(2x)). Confrontare il grafico con quello a nodi equispaziati.

**h)** Commentare la differenza di comportamento tra i due polinomi all'aumentare ipotetico del numero di nodi, collegandolo al fenomeno di Runge.

**i)** Usando solo i tre punti (−1,1.1), (0,0.3), (1,1.4) del punto b), costruire la tabella delle differenze divise e ricavare il polinomio interpolante in forma di Newton. Verificare che valutando questo polinomio nei tre nodi si ottengano esattamente i valori yᵢ, e confrontare il risultato con il polinomio di Lagrange calcolato al punto b) — devono coincidere.

**j)** Aggiungere alla tabella delle differenze divise del punto i) il quarto punto (2,4.0), e costruire il nuovo polinomio interpolante di grado 3 senza dover ricalcolare da zero le differenze divise già trovate, ma solo aggiungendo una riga alla tabella.

---

Stessa modalità di prima: provali con calma, teoria a mano prima e poi MATLAB, e mandami quello che fai così correggiamo insieme.
