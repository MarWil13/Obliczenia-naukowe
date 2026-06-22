# RodHeatDiffusion.jl

Pakiet w języku Julia do numerycznego modelowania jednowymiarowej dyfuzji ciepła w pręcie złożonym z jednego lub kilku materiałów, np. stal–miedź.

Projekt grupowy na przedmiot **Obliczenia naukowe**, WMat 2026.

## Autorzy

- Maria Wilgosz
- Maja Boruszek
- Jakub Klimek
- Aleksandra Rześniowiecka

## Opis problemu

Rozważamy równanie przewodnictwa ciepła:

```text
rho(x) c(x) dT/dt = d/dx( k(x) dT/dx )
```

z warunkami brzegowymi Dirichleta na końcach pręta. Materiał może być jednorodny albo złożony z kilku segmentów. Na granicach między materiałami przewodność jest liczona za pomocą średniej harmonicznej.

## Zaimplementowane metody

Pakiet udostępnia uniwersalny solver schematu theta:

| Metoda | theta | Dokładność w czasie | Stabilność |
|---|---:|---:|---|
| Euler jawny | 0 | O(dt) | warunkowa, wymaga limitu CFL |
| Euler niejawny | 1 | O(dt) | bezwarunkowa dla dyfuzji |
| Crank-Nicolson | 0.5 | O(dt^2) | bezwarunkowa, ale może oscylować dla dużego dt |

## Struktura repozytorium

```text
src/                         kod pakietu Julia
test/                        testy jednostkowe
examples/                    krótkie przykłady uruchomieniowe
docs/                        opis numeryczny i dokumentacja API
notebooks/                   notebook prezentacyjno-analityczny
figures/                     miejsce na wygenerowane wykresy
.github/workflows/ci.yml     automatyczne testy na GitHub Actions
```

## Instalacja lokalna

Po sklonowaniu repozytorium:

```julia
]
dev .
```

albo z terminala:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## Szybki start

```julia
using RodHeatDiffusion

rod = [
    Segment(name="Stal", length=0.5, k=50.0, rho=7800.0, c=450.0),
    Segment(name="Miedź", length=0.5, k=400.0, rho=8900.0, c=390.0),
]

x, k, rhoc, names = build_rod(rod; Nx=151)
T0 = 100 .* exp.(-((x .- 0.5).^2) ./ (2 * 0.05^2))

times, T_history = solve_rod_crank_nicolson(
    x, k, rhoc;
    T_initial = T0,
    T_left = 0.0,
    T_right = 0.0,
    dt = 1.0,
    t_end = 500.0,
    save_every = 10,
)
```

## Testy

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
```

Testy sprawdzają m.in. średnią harmoniczną, budowę siatki, przewodności na granicach komórek, limit stabilności metody jawnej oraz zgodność stanu ustalonego z rozwiązaniem analitycznym.

## Dokumentacja

Opis eksportowanych funkcji znajduje się w pliku [`docs/api.md`](docs/api.md). Wszystkie eksportowane funkcje mają też docstringi dostępne z poziomu Julii, np.:

```julia
?solve_rod_theta
?build_rod
?Segment
```

## Analiza, wykresy i notebook

Opis założeń, algorytmów, złożoności, możliwych problemów, dokładności, wydajności i praktycznego znaczenia znajduje się w pliku [`docs/projekt_numeryczny.md`](docs/projekt_numeryczny.md).

Notebook z analizą i wykresami znajduje się w pliku [`notebooks/projekt_dyfuzja.ipynb`](notebooks/projekt_dyfuzja.ipynb).

Wygenerowane wykresy i animacje znajdują się w katalogu [`figures/`](figures/).

Przykładowe skrypty uruchomieniowe znajdują się w katalogu [`examples/`](examples/).
