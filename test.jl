# Creates a function to run the vfi algorithm
function vfi_solver(;
    alpha=0.4,
    delta=0.1,
    A=1.0,
    beta=0.95,
    N=11,
    k_min=0.1,
    k_max=50.0,
    tolerance=1e-8,
    maxiter=10_000,

    # NEW: stopping rule switch
    stop_rule::Symbol=:policy,   # :policy or :value
    vtol=1e-10                   # only used if stop_rule = :value
)

    # Capital grid
    k = collect(range(k_min, k_max, length=N))

    V = zeros(N)
    V_new = similar(V)

    s = zeros(N)
    s_new = similar(s)

    max_index = zeros(Int, N)

    # Runs VFI algorithm
    for iter in 1:maxiter # for testing that loop converges, (to be replaced with while loop once code works)

        # maximization step
        for i in 1:N

            # For each element in the grid defining the functions (V(k), s(k)) computes the index of k
            # which maximizes the objective function

            # Computes the value of c for each choice of k and stores it in an array
            c = A * k[i]^alpha + (1 - delta) * k[i] .- k
            val = fill(-Inf, N)

            # Loops over the array of c and calculates the objective function evaluated at the i'th point on the grid
            # Using the j'th level of capital in the next period
            for j in 1:N
                if c[j] > 0
                    val[j] = log(c[j]) + beta * V[j]
                end
            end

            # Defines the index of k which maximizes the objective function
            max_index[i] = argmax(val)
        end

        # Evaluation step
        for i in 1:N

            # Evaluates V_{n+1} on the i-th step of the grid
            c = A * k[i]^alpha + (1 - delta) * k[i] - k[max_index[i]]
            if c > 0
                V_new[i] = log(c) + beta * V[max_index[i]]
            else
                V_new[i] = -Inf
            end

            s_new[i] = k[max_index[i]]
        end

        # Checks if stopping criteria is met
        if stop_rule == :policy
            # policy convergence (your original rule)
            if maximum(abs.(s_new .- s)) < tolerance
                V .= V_new
                s .= s_new
                return "Iterations until converge (policy): $(iter)", k, V, s
            end
        elseif stop_rule == :value
            # value convergence (THIS is the “other stopping rule” that makes V match closed form)
            if maximum(abs.(V_new .- V)) < vtol
                V .= V_new
                s .= s_new
                return "Iterations until converge (value): $(iter)", k, V, s
            end
        else
            error("stop_rule must be :policy or :value")
        end

        V .= V_new
        s .= s_new
    end

    # If convergence is not reached
    return "Did not converge: $(maxiter)", k, V, s
end


using Plots
using Measures


# c)

# setting delta = 1 and repeating part a)
# IMPORTANT: use stop_rule=:value so value function is actually converged
sol11 = vfi_solver(N=11, delta=1, stop_rule=:value, vtol=1e-10)
sol101 = vfi_solver(N=101, delta=1, stop_rule=:value, vtol=1e-10)
sol1001 = vfi_solver(N=1001, delta=1, stop_rule=:value, vtol=1e-10)

print(sol11[1])
print(sol101[1])
print(sol1001[1])


# plots with N = 11
p1 = plot(sol11[2], sol11[3], xlabel="Capital (k)",
    ylabel="Value Function V(k)",
    legend=false);

p2 = plot(sol11[2], sol11[4], xlabel="Capital (k)",
    ylabel="Policy Function g(k)",
    legend=false);

display(plot(p1, p2, layout=(1, 2),
    plot_title="Value and Policy Functions N = 11, Delta = 1 (value-stopping)",
    size=(1000, 500),
    left_margin=5mm,
    bottom_margin=5mm))


# plots with N = 101
p1 = plot(sol101[2], sol101[3], xlabel="Capital (k)",
    ylabel="Value Function V(k)",
    legend=false);

p2 = plot(sol101[2], sol101[4], xlabel="Capital (k)",
    ylabel="Policy Function g(k)",
    legend=false);

display(plot(p1, p2, layout=(1, 2),
    plot_title="Value and Policy Functions N = 101, Delta = 1 (value-stopping)",
    size=(1000, 500),
    left_margin=5mm,
    bottom_margin=5mm))


# plots with N = 1001
p1 = plot(sol1001[2], sol1001[3], xlabel="Capital (k)",
    ylabel="Value Function V(k)",
    legend=false);

p2 = plot(sol1001[2], sol1001[4], xlabel="Capital (k)",
    ylabel="Policy Function g(k)",
    legend=false);

display(plot(p1, p2, layout=(1, 2),
    plot_title="Value and Policy Functions N = 1001, Delta = 1 (value-stopping)",
    size=(1000, 500),
    left_margin=5mm,
    bottom_margin=5mm))


# Analytical solution
alpha = 0.4
A = 1.0
beta = 0.95

a = (1 / (1 - beta)) * (1 / (1 - alpha * beta)) * (log(A) + (1 - alpha * beta) * log(1 - alpha * beta) + alpha * beta * log(alpha * beta))
b = alpha / (1 - alpha * beta)

PolicyFunction = alpha * beta * A .* sol1001[2] .^ alpha
ValueFunction = a .+ b .* log.(sol1001[2])


# Plots closed form vs VFI value function solution
p1 = plot(sol1001[2], ValueFunction, xlabel="Capital (K)",
    ylabel="Value Function V(k)",
    legend=true,
    label="Closed Form Solution");

plot!(p1, sol1001[2], sol1001[3],
    label="VFI Solution");

# Plots closed form vs VFI policy function solution
p2 = plot(sol1001[2], PolicyFunction, xlabel="Capital (K)",
    ylabel="Policy Function g(k)",
    legend=true,
    label="Closed Form Solution");

plot!(p2, sol1001[2], sol1001[4],
    label="VFI Solution");

display(plot(p1, p2, layout=(1, 2),
    plot_title="Closed Form vs VFI Solutions N = 1001, Delta = 1 (value-stopping)",
    size=(1000, 500),
    left_margin=5mm,
    bottom_margin=5mm))

