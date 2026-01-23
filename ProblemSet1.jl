
# Question 3

#a) 

#i)

# Creats a function to run the vfi algorithum
function vfi_solver(;
    alpha=0.4,
    delta=0.1,
    A=1.0,
    beta=0.96,
    N=11,
    k_min=0.1,
    k_max=50.0,
    tolerance=1e-8,
    maxiter=10_000
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

        # Checks if stopping criteria is met (policy convergence)
        if maximum(abs.(s_new .- s)) < tolerance
            V .= V_new
            s .= s_new
            return "Iterations until converge: $(iter)", k, V, s
        else
            V .= V_new
            s .= s_new
        end
    end

    # If convergence is not reached
    return "Did not converge: $(maxiter)", k, V, s
end

# Call vfi_solver() to return number of iterations until convergance, V,s
print(vfi_solver(N=11)[1])
print(vfi_solver(N=101)[1])
print(vfi_solver(N=1001)[1])


using Plots
# plots with N = 11
plot(vfi_solver(N=11)[2], vfi_solver(N=11)[3], xlabel="Capital (k)",
    ylabel="Value Function V(k)",
    title="Value Function Iteration Result",
    legend=false)


plot(vfi_solver(N=11)[2], vfi_solver(N=11)[4], xlabel="Capital (k)",
    ylabel="Savings Function g(k)",
    title="Value Function Iteration Result",
    legend=false)


# plots with N = 101
plot(vfi_solver(N=101)[2], vfi_solver(N=101)[3], xlabel="Capital (k)",
    ylabel="Value Function V(k)",
    title="Value Function Iteration Result",
    legend=false)


plot(vfi_solver(N=101)[2], vfi_solver(N=101)[4], xlabel="Capital (k)",
    ylabel="Savings Function g(k)",
    title="Value Function Iteration Result",
    legend=false)

# plots with N = 1001
plot(vfi_solver(N=1001)[2], vfi_solver(N=1001)[3], xlabel="Capital (k)",
    ylabel="Value Function V(k)",
    title="Value Function Iteration Result",
    legend=false)


plot(vfi_solver(N=1001)[2], vfi_solver(N=1001)[4], xlabel="Capital (k)",
    ylabel="Savings Function g(k)",
    title="Value Function Iteration Result",
    legend=false)

# b) 

# changing tolerance to 1e-5
print(vfi_solver(N=11, tolerance=1e-5)[1])
print(vfi_solver(N=101, tolerance=1e-5)[1])
print(vfi_solver(N=1001, tolerance=1e-5)[1])


# plots with N = 11
plot(vfi_solver(N=11, tolerance=1e-5)[2], vfi_solver(N=11, tolerance=1e-5)[3], xlabel="Capital (k)",
    ylabel="Value Function V(k)",
    title="Value Function Iteration Result",
    legend=false)


plot(vfi_solver(N=11, tolerance=1e-5)[2], vfi_solver(N=11, tolerance=1e-5)[4], xlabel="Capital (k)",
    ylabel="Savings Function g(k)",
    title="Value Function Iteration Result",
    legend=false)


# plots with N = 101
plot(vfi_solver(N=101, tolerance=1e-5)[2], vfi_solver(N=101, tolerance=1e-5)[3], xlabel="Capital (k)",
    ylabel="Value Function V(k)",
    title="Value Function Iteration Result",
    legend=false)


plot(vfi_solver(N=101, tolerance=1e-5)[2], vfi_solver(N=101, tolerance=1e-5)[4], xlabel="Capital (k)",
    ylabel="Savings Function g(k)",
    title="Value Function Iteration Result",
    legend=false)

# plots with N = 1001
plot(vfi_solver(N=1001, tolerance=1e-5)[2], vfi_solver(N=1001, tolerance=1e-5)[3], xlabel="Capital (k)",
    ylabel="Value Function V(k)",
    title="Value Function Iteration Result",
    legend=false)


plot(vfi_solver(N=1001, tolerance=1e-5)[2], vfi_solver(N=1001, tolerance=1e-5)[4], xlabel="Capital (k)",
    ylabel="Savings Function g(k)",
    title="Value Function Iteration Result",
    legend=false)






# changing tolerance to 1e-6
print(vfi_solver(N=11, tolerance=1e-6)[1])
print(vfi_solver(N=101, tolerance=1e-6)[1])
print(vfi_solver(N=1001, tolerance=1e-6)[1])


# plots with N = 11
plot(vfi_solver(N=11, tolerance=1e-6)[2], vfi_solver(N=11, tolerance=1e-6)[3], xlabel="Capital (k)",
    ylabel="Value Function V(k)",
    title="Value Function Iteration Result",
    legend=false)


plot(vfi_solver(N=11, tolerance=1e-6)[2], vfi_solver(N=11, tolerance=1e-6)[4], xlabel="Capital (k)",
    ylabel="Savings Function g(k)",
    title="Value Function Iteration Result",
    legend=false)


# plots with N = 101
plot(vfi_solver(N=101, tolerance=1e-6)[2], vfi_solver(N=101, tolerance=1e-6)[3], xlabel="Capital (k)",
    ylabel="Value Function V(k)",
    title="Value Function Iteration Result",
    legend=false)


plot(vfi_solver(N=101, tolerance=1e-6)[2], vfi_solver(N=101, tolerance=1e-6)[4], xlabel="Capital (k)",
    ylabel="Savings Function g(k)",
    title="Value Function Iteration Result",
    legend=false)

# plots with N = 1001
plot(vfi_solver(N=1001, tolerance=1e-6)[2], vfi_solver(N=1001, tolerance=1e-6)[3], xlabel="Capital (k)",
    ylabel="Value Function V(k)",
    title="Value Function Iteration Result",
    legend=false)


plot(vfi_solver(N=1001, tolerance=1e-6)[2], vfi_solver(N=1001, tolerance=1e-6)[4], xlabel="Capital (k)",
    ylabel="Savings Function g(k)",
    title="Value Function Iteration Result",
    legend=false)
