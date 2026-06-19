"""
    RodHeatDiffusion

Pakiet do numerycznego modelowania dyfuzji ciepła w jednowymiarowym pręcie,
złożonym z jednego lub kilku materiałów (np. stal–miedź).

Implementuje rodzinę schematów θ dla równania:

    ρ(x) c(x) ∂T/∂t = ∂/∂x( k(x) ∂T/∂x )

obejmującą jako szczególne przypadki:
- Euler jawny      (θ = 0),
- Euler niejawny    (θ = 1),
- Crank–Nicolson    (θ = 1/2).

# Eksportowane elementy
- [`Segment`](@ref) — opis jednego segmentu materiałowego pręta,
- [`print_segments`](@ref) — wypisanie struktury pręta i jego oporu cieplnego,
- [`build_rod`](@ref) — budowa siatki obliczeniowej z przypisanymi własnościami materiałów,
- [`harmonic_mean`](@ref) — średnia harmoniczna dwóch przewodności,
- [`face_conductivity`](@ref) — przewodności na granicach komórek siatki,
- [`build_operator`](@ref) — rzadka macierz operatora przewodzenia ciepła,
- [`explicit_dt_limit`](@ref) — limit stabilności CFL dla metody jawnej,
- [`solve_rod_theta`](@ref) — uniwersalny solver schematu θ,
- [`solve_rod_explicit`](@ref), [`solve_rod_implicit`](@ref), [`solve_rod_crank_nicolson`](@ref) — konkretne metody.
"""
module RodHeatDiffusion

using LinearAlgebra
using SparseArrays
using Printf

export Segment, print_segments, build_rod,
       harmonic_mean, face_conductivity, build_operator, explicit_dt_limit,
       solve_rod_theta, solve_rod_explicit, solve_rod_implicit, solve_rod_crank_nicolson

#materialy
"""
    Segment(; name, length, k, rho, c)

Pojedynczy odcinek pręta o jednorodnych właściwościach materiałowych.

# Pola
- `name::String`   — nazwa materiału (np. `"Stal"`, `"Miedź"`),
- `length::Float64` — długość odcinka `[m]`,
- `k::Float64`      — przewodność cieplna `[W/(m·K)]`,
- `rho::Float64`    — gęstość `[kg/m³]`,
- `c::Float64`      — ciepło właściwe `[J/(kg·K)]`.

Pełny pręt opisuje się jako `Vector{Segment}` — kolejne segmenty są
„zespawane" jeden za drugim w porządku, w jakim występują w wektorze.

"""
Base.@kwdef struct Segment
    name::String
    length::Float64
    k::Float64
    rho::Float64
    c::Float64
end

"""
    print_segments(segments::Vector{Segment})

Wypisuje na konsolę strukturę pręta (nazwę, długość, przewodność,
pojemność cieplną na jednostkę objętości) oraz całkowity zastępczy
opór cieplny pręta

    R = Σᵢ length_i / k_i   [m²K/W]

przy założeniu jednostkowego przekroju poprzecznego.
"""
function print_segments(segments::Vector{Segment})
    println("Struktura pręta:")
    for s in segments
        @printf("%-15s długość = %.3f m, k = %6.2f W/(mK), rho*c = %.2e J/(m³K)\n",
                s.name, s.length, s.k, s.rho * s.c)
    end

    R = sum(s.length / s.k for s in segments)
    @printf("\nCałkowity opór cieplny R = %.5f m²K/W\n", R)
end


#siatka obliczeniowa
"""
    build_rod(segments::Vector{Segment}; Nx::Int=151) -> (x, k, rhoc, names)

Buduje jednorodną siatkę `Nx` punktów na całej długości pręta i przypisuje
każdemu punktowi właściwości materiału, do którego należy.

# Argumenty
- `segments` — wektor segmentów materiałowych (patrz [`Segment`](@ref)),
- `Nx`       — liczba punktów siatki (domyślnie 151).

# Zwraca
- `x::Vector{Float64}`     — współrzędne punktów siatki `[m]`,
- `k::Vector{Float64}`     — przewodność cieplna w każdym punkcie `[W/(mK)]`,
- `rhoc::Vector{Float64}`  — `ρ·c` w każdym punkcie `[J/(m³K)]`,
- `names::Vector{String}`  — nazwa materiału w każdym punkcie.

"""
function build_rod(segments::Vector{Segment}; Nx::Int=151)
    L = sum(s.length for s in segments)
    x = collect(range(0.0, L; length=Nx))
    edges = cumsum([s.length for s in segments])

    k = zeros(Float64, Nx)
    rhoc = zeros(Float64, Nx)
    names = Vector{String}(undef, Nx)

    for (i, xi) in enumerate(x)
        # +eps pomaga poprawnie przypisać punkt x=0 do pierwszego segmentu
        j = searchsortedfirst(edges, xi + eps(Float64))
        j = clamp(j, 1, length(segments))
        segment = segments[j]
        k[i] = segment.k
        rhoc[i] = segment.rho * segment.c
        names[i] = segment.name
    end

    return x, k, rhoc, names
end


#dyskretyzacja operatora przewodzenia
"""
    harmonic_mean(a, b)

Średnia harmoniczna dwóch liczb: `2ab/(a+b)`.

Używana do wyznaczania przewodności cieplnej na granicy dwóch komórek
o różnych materiałach. Jest fizycznie poprawniejsza niż średnia
arytmetyczna, ponieważ na styku dwóch metali głównym ograniczeniem dla
przepływu ciepła jest materiał o mniejszej przewodności (analogia do
oporów połączonych szeregowo).
"""
harmonic_mean(a, b) = 2a * b / (a + b)

"""
    face_conductivity(k::AbstractVector) -> Vector

Wyznacza przewodność cieplną na granicach (ścianach) między sąsiednimi
punktami siatki, korzystając ze średniej harmonicznej (patrz
[`harmonic_mean`](@ref)).

Dla wektora `k` o długości `N` zwraca wektor długości `N-1`.
"""
function face_conductivity(k)
    return [harmonic_mean(k[i], k[i+1]) for i in 1:length(k)-1]
end

"""
    build_operator(k, rhoc, dx) -> (L, leftcoef, rightcoef)

Buduje rzadką (trójprzekątniową) macierz operatora przewodzenia ciepła
dla punktów wewnętrznych siatki, w postaci konserwatywnej:

    ρᵢcᵢ dTᵢ/dt = (1/Δx²) [ k_{i+1/2}(T_{i+1} - Tᵢ) − k_{i−1/2}(Tᵢ − T_{i−1}) ]

# Argumenty
- `k`    — przewodność cieplna w każdym punkcie siatki,
- `rhoc` — `ρ·c` w każdym punkcie siatki,
- `dx`   — odstęp między punktami siatki `[m]`.

# Zwraca
- `L::SparseMatrixCSC` — macierz operatora dla punktów wewnętrznych (rozmiar `(N-2)×(N-2)`),
- `leftcoef`, `rightcoef` — wektory współczynników wiążących punkty wewnętrzne
  z warunkami brzegowymi (lewym i prawym końcem pręta).
"""
function build_operator(k, rhoc, dx)
    N = length(k)
    n = N - 2
    kh = face_conductivity(k)

    rows, cols, vals = Int[], Int[], Float64[]
    leftcoef, rightcoef = zeros(n), zeros(n)

    for j in 1:n
        i = j + 1
        aL = kh[i-1] / (rhoc[i] * dx^2)
        aR = kh[i]   / (rhoc[i] * dx^2)

        push!(rows, j); push!(cols, j); push!(vals, -(aL + aR))
        if j > 1
            push!(rows, j); push!(cols, j - 1); push!(vals, aL)
        else
            leftcoef[j] = aL
        end
        if j < n
            push!(rows, j); push!(cols, j + 1); push!(vals, aR)
        else
            rightcoef[j] = aR
        end
    end
    return sparse(rows, cols, vals, n, n), leftcoef, rightcoef
end

"""
    explicit_dt_limit(k, rhoc, dx; safety=0.8) -> Float64

Wyznacza bezpieczny krok czasowy dla jawnej metody Eulera, wynikający
z warunku stabilności CFL:

    Δt ≤ safety · min_i  ρᵢcᵢ(Δx)² / (k_{i-1/2} + k_{i+1/2})

# Argumenty
- `k`, `rhoc`, `dx` — jak w [`build_operator`](@ref),
- `safety` — współczynnik bezpieczeństwa `∈ (0, 1)` (domyślnie 0.8);
  wartości bliższe 1 są ryzykowne ze względu na błędy numeryczne.
"""
function explicit_dt_limit(k, rhoc, dx; safety=0.8)
    kh = face_conductivity(k)
    rates = [(kh[i-1] + kh[i]) / (rhoc[i] * dx^2) for i in 2:length(k)-1]
    return safety / maximum(rates)
end

#solvery czasowe - schemat θ
"""
    solve_rod_theta(x, k, rhoc; T_initial, T_left, T_right, dt, t_end, theta, save_every=1)
        -> (times, T_history)

Uniwersalny solver schematu θ dla równania dyfuzji ciepła w pręcie
z warunkami brzegowymi Dirichleta:

    (I − θΔt L) T^{n+1} = (I + (1−θ)Δt L) T^n

# Argumenty
- `x`, `k`, `rhoc`  — siatka i właściwości materiałowe (z [`build_rod`](@ref)),
- `T_initial`       — wektor temperatury początkowej,
- `T_left`, `T_right` — temperatura na końcach pręta; albo stała liczba,
  albo funkcja czasu `t -> wartość`,
- `dt`        — krok czasowy `[s]`,
- `t_end`     — czas końcowy symulacji `[s]`,
- `theta`     — parametr schematu θ ∈ [0, 1],
- `save_every` — co ile kroków zapisywać stan (domyślnie każdy krok).

# Zwraca
- `times::Vector{Float64}` — chwile czasu, dla których zapisano rozwiązanie,
- `T_history::Matrix{Float64}` — wiersz `i` to profil temperatury `T(x, times[i])`.

Dla `theta > 0` macierz układu `A = I - θΔt·L` jest faktoryzowana
(`factorize`) tylko raz, przed pętlą czasową, i wykorzystywana
ponownie w każdym kroku — to znacząco redukuje koszt metod niejawnych.
"""
function solve_rod_theta(x, k, rhoc; T_initial, T_left, T_right, dt, t_end, theta, save_every=1)
    dx = x[2] - x[1]
    Lop, leftcoef, rightcoef = build_operator(k, rhoc, dx)
    Id = spdiagm(0 => ones(length(x) - 2))

    A = Id - theta * dt * Lop
    B = Id + (1 - theta) * dt * Lop
    F = factorize(A)

    T = copy(T_initial)

    get_T_left(t)  = T_left isa Function ? T_left(t) : T_left
    get_T_right(t) = T_right isa Function ? T_right(t) : T_right

    T[1] = get_T_left(0.0)
    T[end] = get_T_right(0.0)

    nt = ceil(Int, t_end / dt)
    times = Float64[]
    saved = Vector{Vector{Float64}}()

    for n in 0:nt
        t = n * dt
        if n % save_every == 0
            push!(times, t)
            push!(saved, copy(T))
        end
        n == nt && break

        b_now = leftcoef .* get_T_left(t) .+ rightcoef .* get_T_right(t)
        b_next = leftcoef .* get_T_left(t + dt) .+ rightcoef .* get_T_right(t + dt)

        rhs = B * T[2:end-1] .+ dt .* ((1 - theta) .* b_now .+ theta .* b_next)
        T[2:end-1] .= F \ rhs
        T[1] = get_T_left(t + dt)
        T[end] = get_T_right(t + dt)
    end

    return times, permutedims(hcat(saved...))
end

"""
    solve_rod_explicit(args...; kwargs...)

Jawna metoda Eulera — wrapper na [`solve_rod_theta`](@ref) z `theta = 0`.
Wymaga `dt` poniżej limitu CFL (patrz [`explicit_dt_limit`](@ref)),
inaczej rozwiązanie jest niestabilne.
"""
solve_rod_explicit(args...; kwargs...) = solve_rod_theta(args...; theta=0.0, kwargs...)

"""
    solve_rod_crank_nicolson(args...; kwargs...)

Schemat Cranka–Nicolsona — wrapper na [`solve_rod_theta`](@ref) z `theta = 0.5`.
Bezwarunkowo stabilny, drugi rząd dokładności w czasie; dla bardzo dużego
`dt` i nieciągłych warunków początkowo-brzegowych może generować oscylacje.
"""
solve_rod_crank_nicolson(args...; kwargs...) = solve_rod_theta(args...; theta=0.5, kwargs...)

"""
    solve_rod_implicit(args...; kwargs...)

Niejawna metoda Eulera — wrapper na [`solve_rod_theta`](@ref) z `theta = 1`.
Bezwarunkowo stabilna dla równania dyfuzji; dla dużego `dt` silnie tłumi
(wygładza) rozwiązanie.
"""
solve_rod_implicit(args...; kwargs...) = solve_rod_theta(args...; theta=1.0, kwargs...)

end 
