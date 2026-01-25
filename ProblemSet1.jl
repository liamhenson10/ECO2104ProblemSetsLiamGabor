# Question 3

#a) 

#i)

# Creates a function to run the vfi algorithm
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
sol11 = vfi_solver(N=11);
sol101 = vfi_solver(N=101);
sol1001 = vfi_solver(N=1001);



using Plots
# plots with N = 11
plot(sol11[2], sol11[3], xlabel="Capital (k)",
    ylabel="Value Function V(k)",
    title="Value Function N = 11",
    legend=false)


plot(sol11[2], sol11[4], xlabel="Capital (k)",
    ylabel="Policy Function g(k)",
    title="Policy Function N = 11",
    legend=false)


# plots with N = 101
plot(sol101[2], sol101[3], xlabel="Capital (k)",
    ylabel="Value Function V(k)",
    title=" Value Function N = 101",
    legend=false)


plot(sol101[2], sol101[4], xlabel="Capital (k)",
    ylabel="Policy Function g(k)",
    title="Policy Function N = 101",
    legend=false)

# plots with N = 1001
plot(sol1001[2], sol1001[3], xlabel="Capital (k)",
    ylabel="Value Function V(k)",
    title="Value Function N = 1001",
    legend=false)


plot(sol1001[2], sol1001[4], xlabel="Capital (k)",
    ylabel="Policy Function g(k)",
    title="Policy Function N = 1001",
    legend=false)



# b) 

# changing tolerance to 1e-5
sol11 = vfi_solver(N=11, tolerance=1e-5);
sol101 = vfi_solver(N=101, tolerance=1e-5);
sol1001 = vfi_solver(N=1001, tolerance=1e-5);


# plots with N = 11
plot(sol11[2], sol11[3], xlabel="Capital (k)",
    ylabel="Value Function V(k)",
    title="Value Function N = 11, Tolerance = 1e-5",
    legend=false)


plot(sol11[2], sol11[4], xlabel="Capital (k)",
    ylabel="Policy Function g(k)",
    title="Policy Function N = 11, Tolerance = 1e-5",
    legend=false)


# plots with N = 101
plot(sol101[2], sol101[3], xlabel="Capital (k)",
    ylabel="Value Function V(k)",
    title="Value Function N = 101, Tolerance = 1e-5",
    legend=false)


plot(sol101[2], sol101[4], xlabel="Capital (k)",
    ylabel="Policy Function g(k)",
    title="Policy Function N = 101, Tolerance = 1e-5",
    legend=false)

# plots with N = 1001
plot(sol1001[2], sol1001[3], xlabel="Capital (k)",
    ylabel="Value Function V(k)",
    title="Value Function N = 1001, Tolerance = 1e-5",
    legend=false)


plot(sol1001[2], sol1001[4], xlabel="Capital (k)",
    ylabel="Policy Function g(k)",
    title="Policy Function N = 11, Tolerance = 1e-5",
    legend=false)








# changing tolerance to 1e-6
sol11 = vfi_solver(N=11, tolerance=1e-6)
sol101 = vfi_solver(N=101, tolerance=1e-6)
sol1001 = vfi_solver(N=1001, tolerance=1e-6)


# plots with N = 11
plot(sol11[2], sol11[3], xlabel="Capital (k)",
    ylabel="Value Function V(k)",
    title="Value Function N = 11, Tolerance = 1e-6",
    legend=false)


plot(sol11[2], sol11[4], xlabel="Capital (k)",
    ylabel="Policy Function g(k)",
    title="Policy Function N = 11, Tolerance = 1e-6",
    legend=false)


# plots with N = 101
plot(sol101[2], sol101[3], xlabel="Capital (k)",
    ylabel="Value Function V(k)",
    title="Value Function N = 101, Tolerance = 1e-6",
    legend=false)


plot(sol101[2], sol101[4], xlabel="Capital (k)",
    ylabel="Policy Function g(k)",
    title="Policy Function N = 101, Tolerance = 1e-6",
    legend=false)

# plots with N = 1001
plot(sol1001[2], sol1001[3], xlabel="Capital (k)",
    ylabel="Value Function V(k)",
    title="Value Function N = 1001, Tolerance = 1e-6",
    legend=false)


plot(sol1001[2], sol1001[4], xlabel="Capital (k)",
    ylabel="Policy Function g(k)",
    title="Policy Function N = 11, Tolerance = 1e-6",
    legend=false)



# c) 

# setting delta = 1 and repeating part a) 
sol11 = vfi_solver(N=11, delta=1)
sol101 = vfi_solver(N=101, delta=1)
sol1001 = vfi_solver(N=1001, delta=1)


# plots with N = 11
plot(sol11[2], sol11[3], xlabel="Capital (k)",
    ylabel="Value Function V(k)",
    title="Value Function N = 11, Delta = 1 ",
    legend=false)


plot(sol11[2], sol11[4], xlabel="Capital (k)",
    ylabel="Policy Function g(k)",
    title="Policy Function N = 11, Delta = 1 ",
    legend=false)


# plots with N = 101
plot(sol101[2], sol101[3], xlabel="Capital (k)",
    ylabel="Value Function V(k)",
    title="Value Function N = 101, Delta = 1 ",
    legend=false)


plot(sol101[2], sol101[4], xlabel="Capital (k)",
    ylabel="Policy Function g(k)",
    title="Policy Function N = 101, Delta = 1 ",
    legend=false)

# plots with N = 1001
plot(sol1001[2], sol1001[3], xlabel="Capital (k)",
    ylabel="Value Function V(k)",
    title="Value Function N = 1001, Delta = 1 ",
    legend=false)


plot(sol1001[2], sol1001[4], xlabel="Capital (k)",
    ylabel="Policy Function g(k)",
    title="Policy Function N = 11, Delta = 1 ",
    legend=false)