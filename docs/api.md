# Dokumentacja eksportowanych elementów pakietu

Ten plik streszcza eksportowane typy i funkcje pakietu `RodHeatDiffusion`.

## `Segment`

Opisuje jeden jednorodny fragment pręta.

Pola:

- `name::String` — nazwa materiału,
- `length::Float64` — długość fragmentu,
- `k::Float64` — przewodność cieplna,
- `rho::Float64` — gęstość,
- `c::Float64` — ciepło właściwe.

## `print_segments(segments)`

Wypisuje strukturę pręta oraz przybliżony całkowity opór cieplny wynikający z segmentów.

## `build_rod(segments; Nx=151)`

Buduje jednorodną siatkę punktów na całej długości pręta i przypisuje każdemu punktowi właściwości materiału.

Zwraca:

- `x` — współrzędne punktów,
- `k` — przewodność cieplną w punktach,
- `rhoc` — iloczyn `rho*c`,
- `names` — nazwy materiałów przypisane do punktów.

## `harmonic_mean(a, b)`

Liczy średnią harmoniczną dwóch liczb. W projekcie jest używana do modelowania przewodności na styku dwóch materiałów, ponieważ odpowiada szeregowemu charakterowi oporów cieplnych.

## `face_conductivity(k)`

Wyznacza przewodności na granicach między sąsiednimi komórkami siatki.

## `build_operator(k, rhoc, dx)`

Tworzy rzadką macierz operatora dyfuzji dla punktów wewnętrznych. Zwraca macierz operatora oraz współczynniki związane z lewym i prawym warunkiem brzegowym.

## `explicit_dt_limit(k, rhoc, dx; safety=0.8)`

Wyznacza bezpieczny krok czasowy dla jawnej metody Eulera. Przekroczenie tego limitu może prowadzić do niestabilności numerycznej.

## `solve_rod_theta(...)`

Uniwersalny solver schematu theta:

- `theta = 0` — Euler jawny,
- `theta = 1` — Euler niejawny,
- `theta = 0.5` — Crank-Nicolson.

Zwraca wektor zapisanych czasów oraz macierz historii temperatury.

## `solve_rod_explicit(...)`

Wrapper na `solve_rod_theta` z `theta = 0`. Metoda jest prosta i szybka w jednym kroku, ale stabilna tylko przy odpowiednio małym `dt`.

## `solve_rod_implicit(...)`

Wrapper na `solve_rod_theta` z `theta = 1`. Metoda jest bezwarunkowo stabilna dla dyfuzji, ale dla dużych kroków czasowych może nadmiernie wygładzać rozwiązanie.

## `solve_rod_crank_nicolson(...)`

Wrapper na `solve_rod_theta` z `theta = 0.5`. Metoda ma drugi rząd dokładności w czasie, ale przy bardzo dużym kroku czasowym i ostrych danych początkowych może dawać oscylacje.

## `final_temperature(solver, x, k, rhoc; T_initial, T_left, T_right, dt, t_end)`

Uruchamia wybrany solver i zwraca końcowy profil temperatury. Funkcja jest używana do porównywania metod dla tych samych warunków początkowych, brzegowych i czasu końcowego.

## `error_norms(T, Tref, dx)`

Liczy błędy między profilem temperatury `T` a dokładniejszym wynikiem numerycznym `Tref`. Zwraca błąd `L2`, opisujący różnicę całego profilu, oraz błąd `L∞`, oznaczający największe lokalne odchylenie.
