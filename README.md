# RodHeatDiffusion.jl

Pakiet Julia do numerycznego modelowania jednowymiarowej dyfuzji ciepła
w pręcie złożonym z jednego lub kilku materiałów (np. stal–miedź).

Projekt grupowy na przedmiot *Obliczenia naukowe*, WMat 2026.

**Autorzy:** Maria Wilgosz, Maja Boruszek, Jakub Klimek, Aleksandra Rześniowiecka

## Opis problemu

Rozważamy równanie przewodnictwa ciepła dla materiału o właściwościach
zmiennych w przestrzeni:

```
ρ(x) c(x) ∂T/∂t = ∂/∂x( k(x) ∂T/∂x )
```

z warunkami brzegowymi Dirichleta na końcach pręta. Na granicach między
różnymi materiałami przewodność wyznaczana jest za pomocą średniej
harmonicznej, co poprawnie odwzorowuje efekt „dławienia" przepływu ciepła
przez materiał o mniejszej przewodności.

## Zaimplementowane metody

Pakiet udostępnia jeden uniwersalny solver schematu θ (theta-method),
z którego wynikają trzy klasyczne metody jako szczególne przypadki:

| Metoda           | θ   | Rząd dokładności (czas) | Stabilność                  |
|------------------|-----|--------------------------|------------------------------|
| Euler jawny      | 0   | O(Δt)                    | warunkowa (limit CFL)         |
| Euler niejawny   | 1   | O(Δt)                    | bezwarunkowa                  |
| Crank–Nicolson   | 0.5 | O(Δt²)                   | bezwarunkowa (może oscylować) |

Dyskretyzacja przestrzenna wykorzystuje rzadką macierz trójprzekątniową,
a w metodach niejawnych macierz układu jest faktoryzowana jednokrotnie
i wielokrotnie wykorzystywana w pętli czasowej.

## Instalacja

```julia
] add https://github.com/<twoja-nazwa-uzytkownika>/RodHeatDiffusion.jl
```

albo lokalnie, po sklonowaniu repozytorium:

```julia
] dev .
```

## Szybki start

```julia
using RodHeatDiffusion

# Definicja pręta złożonego ze stali i miedzi
rod = [
    Segment(name="Stal",  length=0.5, k=50.0,  rho=7800.0, c=450.0),
    Segment(name="Miedź", length=0.5, k=400.0, rho=8900.0, c=390.0),
]
print_segments(rod)

# Siatka obliczeniowa
x, k, rhoc, names = build_rod(rod; Nx=151)

# Warunek początkowy: gaussowski impuls w środku pręta
T0 = 100 .* exp.(-((x .- 0.5).^2) ./ (2 * 0.05^2))

# Symulacja metodą Cranka-Nicolson, T_left = 0°C, T_right = 0°C
times, T_history = solve_rod_crank_nicolson(
    x, k, rhoc;
    T_initial = T0,
    T_left  = t -> 0.0,
    T_right = t -> 0.0,
    dt = 1.0,
    t_end = 500.0,
    save_every = 10,
)

# T_history[i, :] to profil temperatury w chwili times[i]
```

Limit stabilności dla metody jawnej można sprawdzić przed symulacją:

```julia
dx = x[2] - x[1]
dt_max = explicit_dt_limit(k, rhoc, dx; safety=0.8)
```

## Dokumentacja

Wszystkie eksportowane funkcje i typy mają docstringi opisujące
parametry, wartości zwracane oraz wzory matematyczne — dostępne w
Julii bezpośrednio przez:

```julia
?solve_rod_theta
?build_rod
?Segment
```

Lista eksportowanych elementów: `Segment`, `print_segments`, `build_rod`,
`harmonic_mean`, `face_conductivity`, `build_operator`, `explicit_dt_limit`,
`solve_rod_theta`, `solve_rod_explicit`, `solve_rod_implicit`,
`solve_rod_crank_nicolson`.

## Testy

```julia
] test
```

Testy sprawdzają m.in. poprawność przypisania materiałów na siatce,
średnią harmoniczną na granicach komórek oraz zgodność stanu ustalonego
wszystkich trzech metod z rozwiązaniem analitycznym dla pręta
jednorodnego.

## Analiza, wykresy i wnioski

Szczegółowa analiza dokładności, wydajności, demonstracje pułapek
numerycznych (np. naruszenie warunku CFL, niespójne warunki brzegowe)
oraz dyskusja praktycznego znaczenia metod znajdują się w osobnym
notatniku prezentacyjnym (Jupyter/Pluto), nie wchodzącym w skład tego
pakietu.
