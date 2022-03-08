include("utils.jl")


# Función a evaluar 1
function ackley(X)
    return -20*ℯ^(-2*√0.5*(X[1]^2+X[2]^2))-ℯ^(0.5*(cos(2*π*X[1])+cos(2*π*X[2])))+ℯ+20
end
# Search domain -5 < x,y <5

n_variables = 2
i_sup_vec = [5,5]
i_inf_vec = [-5,-5]
precision = 4
m = 50

dimension_vec = []
poblacion_inicial_vec = []

for i in 1:n_variables
    length_cromosoma = length_variable(i_sup_vec[i],i_inf_vec[i],precision)
    push!(dimension_vec,length_cromosoma)
    push!(poblacion_inicial_vec,rand_population_binary(m, length_cromosoma))
end

feno = DECODE(n_variables,m,i_sup_vec,i_inf_vec,dimension_vec,poblacion_inicial_vec)
objv = OBJFUN(ackley,feno)
