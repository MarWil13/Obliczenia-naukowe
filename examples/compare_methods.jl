# Porównanie trzech metod czasowych: Euler jawny, Euler niejawny, Crank-Nicolson.
# Uruchomienie z katalogu repozytorium:
#     julia --project=. examples/compare_methods.jl

using RodHeatDiffusion
using Printf

rod = [Segment(name="Stal", length=1.0, k=50.0, rho=7800.0, c=450.0)]
x, k, rhoc, _ = build_rod(rod; Nx=101)
dx = x[2] - x[1]

T0 = zeros(length(x))
T_left = 0.0
T_right = 100.0
T_ref = collect(range(T_left, T_right; length=length(x)))

dt_exp = explicit_dt_limit(k, rhoc, dx; safety=0.8)
t_end = 6000.0

_, sol_exp = solve_rod_explicit(
    x, k, rhoc;
    T_initial=copy(T0), T_left=T_left, T_right=T_right,
    dt=dt_exp, t_end=t_end, save_every=10^9,
)

_, sol_imp = solve_rod_implicit(
    x, k, rhoc;
    T_initial=copy(T0), T_left=T_left, T_right=T_right,
    dt=10.0, t_end=t_end, save_every=10^9,
)

_, sol_cn = solve_rod_crank_nicolson(
    x, k, rhoc;
    T_initial=copy(T0), T_left=T_left, T_right=T_right,
    dt=10.0, t_end=t_end, save_every=10^9,
)

err_exp = maximum(abs.(sol_exp[end, :] .- T_ref))
err_imp = maximum(abs.(sol_imp[end, :] .- T_ref))
err_cn = maximum(abs.(sol_cn[end, :] .- T_ref))

@printf("Błąd względem profilu liniowego w stanie ustalonym:
")
@printf("  Euler jawny:        %.6f
", err_exp)
@printf("  Euler niejawny:     %.6f
", err_imp)
@printf("  Crank-Nicolson:     %.6f
", err_cn)
@printf("
Limit stabilności dla Eulera jawnego: dt <= %.6e s
", dt_exp)
