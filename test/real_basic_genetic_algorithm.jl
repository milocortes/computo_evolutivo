using Printf
using Random
using Distributions


#=
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Algoritmo genético
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
=#

function torneo_k_enario(y,k)
    p = randperm(length(y))
    return p[argmin(y[p[1:k]])]
end

function muestra_uniforme(m,a,b)
    return [a+rand(length(a)).*(b-a) for i in 1:m]
end

function cruza_SBX(p1,p2)
    u = rand()
    nc = 1

    if u <= 0.5
        β = (2*u)^(1/(nc+1))
    else
        β =(1/(2*(1-u)))^(1/(nc+1))
    end

    return [abs.(0.5 .* ((p1+p2) - β .*abs.(p2-p1))),abs.(0.5 .* ((p1+p2) + β .*abs.(p2-p1)))]

end

function mutacion_polinomial(hijo,σ)
    return hijo + randn(length(hijo))*σ
end

function algoritmo_genetico(f,m,i_inf_vec,i_sup_vec,k,it_max)
    # Generamos una población aleatoria de tamaño M
    poblacion = muestra_uniforme(m,i_inf_vec,i_sup_vec);

    for it in 1:it_max
        @printf("Iteración %d \n", it)
        y = f.(poblacion);
        @printf "Mejor valor fun.obj ---> %s. Variables de decision ---> %s\n"  minimum(f.(poblacion)) poblacion[argmin(f.(poblacion))]
        padres = [[torneo_k_enario(y,k),torneo_k_enario(y,k)] for i in y];
        hijos =  [cruza_SBX(poblacion[p[1]],poblacion[p[2]]) for p in padres];
        poblacion = [mutacion_polinomial(hijo,0.1) for hijo in vcat(hijos...)];
    end

    return poblacion[argmin(f.(poblacion))]
end


# Función a evaluar Ackley
function ackley(X)
    return -20*ℯ^(-2*√0.5*(X[1]^2+X[2]^2))-ℯ^(0.5*(cos(2*π*X[1])+cos(2*π*X[2])))+ℯ+20
end
# Search domain -5 < x,y <5
i_sup_vec = [5,5]
i_inf_vec = [-5,-5]
m = 40
it_max = 20
k = 2

algoritmo_genetico(ackley,m,i_inf_vec,i_sup_vec,k,it_max)

# Funcion a evaluar Rosenbrock
function rosenbrock(X)
    suma = 0
    for i in 1:length(X)-1
        suma += (100*(X[i+1]-X[i]^2)^2 + (1-X[i])^2)
    end
    return suma
end


i_sup_vec = [50,50,50]
i_inf_vec = [-50,-50,-50]
m = 40
it_max = 10
k = 2
algoritmo_genetico(rosenbrock,m,i_inf_vec,i_sup_vec,k,it_max)

# Funcion a evaluar Himmelblau
function himmelblau(X)
    return (X[1]^2 + X[2] - 11)^2 + (X[1] + X[2]^2 - 7)^2
end

i_sup_vec = [5,5]
i_inf_vec = [-5,-5]
m = 40
maxiter = 10
k = 2

algoritmo_genetico(himmelblau,m,i_inf_vec,i_sup_vec,k,it_max)


# Funcion a evaluar Eggholder
function eggholder(X)
    return -(X[2]+47)*sin((√abs((X[1]/2)+(X[2]+47)))) - (X[1]*sin(√abs(X[1]-(X[2]+47))))
end


i_sup_vec = [512,512]
i_inf_vec = [-512,-512]
m = 100
maxiter = 10
k = 6

algoritmo_genetico(eggholder,m,i_inf_vec,i_sup_vec,k,it_max)
