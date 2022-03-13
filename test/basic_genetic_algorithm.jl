include("../src/utils.jl")
using Plots
using Printf

# Función a evaluar Ackley
function ackley(X)
    return -20*ℯ^(-2*√0.5*(X[1]^2+X[2]^2))-ℯ^(0.5*(cos(2*π*X[1])+cos(2*π*X[2])))+ℯ+20
end
# Search domain -5 < x,y <5

n_variables = 2
i_sup_vec = [5,5]
i_inf_vec = [-5,-5]
precision = 4
m = 200
maxiter = 50
dimension_vec = []
genotipo = []
length_total_cromosoma = 0

## Generamos población inicial
for i in 1:n_variables
    length_cromosoma = length_variable(i_sup_vec[i],i_inf_vec[i],precision)
    length_total_cromosoma += length_cromosoma
    push!(dimension_vec,length_cromosoma)
    push!(genotipo,rand_population_binary(m, length_cromosoma))
end

## Iniciamos el algoritmo genético
feno = DECODE(n_variables,m,i_sup_vec,i_inf_vec,dimension_vec,genotipo)
objv = OBJFUN(ackley,feno)

resultados = []
mejor_individuo = 0
for it in 1:maxiter
    println(it)
    aptitud = APTITUD(objv,"min")
    seleccion = SELECCION(aptitud,"ruleta",n_variables,genotipo)
    genotipo = CRUZA(seleccion,"unpunto",length_total_cromosoma)
    genotipo = MUTACION(genotipo,length_total_cromosoma,n_variables,dimension_vec)
    feno = DECODE(n_variables,m,i_sup_vec,i_inf_vec,dimension_vec,genotipo)
    objv = OBJFUN(ackley,feno)
    push!(resultados,minimum(objv))
    mejor_individuo = findmin(objv)[2]
    @printf "Mejor valor fun.obj ---> %s. Variables de decision ---> %s" objv[mejor_individuo] feno[mejor_individuo]
end

plot(resultados, title = "Ackley")


# Funcion a evaluar Rosenbrock
function rosenbrock(X)
    suma = 0
    for i in 1:length(X)-1
        suma += (100*(X[i+1]-X[i]^2)^2 + (1-X[i])^2)
    end
    return suma
end


n_variables = 3
i_sup_vec = [50,50,50]
i_inf_vec = [-50,-50,-50]
precision = 4
m = 40
maxiter = 100
dimension_vec = []
genotipo = []
length_total_cromosoma = 0

## Generamos población inicial
for i in 1:n_variables
    length_cromosoma = length_variable(i_sup_vec[i],i_inf_vec[i],precision)
    length_total_cromosoma += length_cromosoma
    push!(dimension_vec,length_cromosoma)
    push!(genotipo,rand_population_binary(m, length_cromosoma))
end

## Iniciamos el algoritmo genético
feno = DECODE(n_variables,m,i_sup_vec,i_inf_vec,dimension_vec,genotipo)
objv = OBJFUN(rosenbrock,feno)

resultados = []
mejor_individuo = 0
for it in 1:maxiter
    println(it)
    aptitud = APTITUD(objv,"min")
    seleccion = SELECCION(aptitud,"ruleta",n_variables,genotipo)
    genotipo = CRUZA(seleccion,"unpunto",length_total_cromosoma)
    genotipo = MUTACION(genotipo,length_total_cromosoma,n_variables,dimension_vec)
    feno = DECODE(n_variables,m,i_sup_vec,i_inf_vec,dimension_vec,genotipo)
    objv = OBJFUN(rosenbrock,feno)
    push!(resultados,minimum(objv))
    mejor_individuo = findmin(objv)[2]
    @printf "Mejor valor fun.obj ---> %s. Variables de decision ---> %s" objv[mejor_individuo] feno[mejor_individuo]
end

plot(resultados)

plot(resultados, title = "Rosenbrock")

# Funcion a evaluar Himmelblau
function himmelblau(X)
    return (X[1]^2 + X[2] - 11)^2 + (X[1] + X[2]^2 - 7)^2
end


n_variables = 2
i_sup_vec = [5,5]
i_inf_vec = [-5,-5]
precision = 2
m = 40
maxiter = 100
dimension_vec = []
genotipo = []
length_total_cromosoma = 0

## Generamos población inicial
for i in 1:n_variables
    length_cromosoma = length_variable(i_sup_vec[i],i_inf_vec[i],precision)
    length_total_cromosoma += length_cromosoma
    push!(dimension_vec,length_cromosoma)
    push!(genotipo,rand_population_binary(m, length_cromosoma))
end

## Iniciamos el algoritmo genético
feno = DECODE(n_variables,m,i_sup_vec,i_inf_vec,dimension_vec,genotipo)
objv = OBJFUN(himmelblau,feno)

resultados = []
mejor_individuo = 0
for it in 1:maxiter
    println(it)
    aptitud = APTITUD(objv,"min")
    seleccion = SELECCION(aptitud,"ruleta",n_variables,genotipo)
    genotipo = CRUZA(seleccion,"unpunto",length_total_cromosoma)
    genotipo = MUTACION(genotipo,length_total_cromosoma,n_variables,dimension_vec)
    feno = DECODE(n_variables,m,i_sup_vec,i_inf_vec,dimension_vec,genotipo)
    objv = OBJFUN(himmelblau,feno)
    push!(resultados,minimum(objv))
    mejor_individuo = findmin(objv)[2]
    @printf "Mejor valor fun.obj ---> %s. Variables de decision ---> %s" objv[mejor_individuo] feno[mejor_individuo]
end

plot(resultados)
plot(resultados, title = "Himmelblau")

# Funcion a evaluar Eggholder
function eggholder(X)
    return -(X[2]+47)*sin((√abs((X[1]/2)+(X[2]+47)))) - (X[1]*sin(√abs(X[1]-(X[2]+47))))
end


n_variables = 2
i_sup_vec = [512,512]
i_inf_vec = [-512,-512]
precision = 2
m = 600
maxiter = 100
dimension_vec = []
genotipo = []
length_total_cromosoma = 0

## Generamos población inicial
for i in 1:n_variables
    length_cromosoma = length_variable(i_sup_vec[i],i_inf_vec[i],precision)
    length_total_cromosoma += length_cromosoma
    push!(dimension_vec,length_cromosoma)
    push!(genotipo,rand_population_binary(m, length_cromosoma))
end

## Iniciamos el algoritmo genético
feno = DECODE(n_variables,m,i_sup_vec,i_inf_vec,dimension_vec,genotipo)
objv = OBJFUN(eggholder,feno)

resultados = []
mejor_individuo = 0
for it in 1:maxiter
    println(it)
    aptitud = APTITUD(objv,"min")
    seleccion = SELECCION(aptitud,"ruleta",n_variables,genotipo)
    genotipo = CRUZA(seleccion,"unpunto",length_total_cromosoma)
    genotipo = MUTACION(genotipo,length_total_cromosoma,n_variables,dimension_vec)
    feno = DECODE(n_variables,m,i_sup_vec,i_inf_vec,dimension_vec,genotipo)
    objv = OBJFUN(eggholder,feno)
    push!(resultados,minimum(objv))
    mejor_individuo = findmin(objv)[2]
    @printf "Mejor valor fun.obj ---> %s. Variables de decision ---> %s" objv[mejor_individuo] feno[mejor_individuo]
end

plot(resultados)
plot(resultados, title = "Eggholder")
