using Random

# Genera población aleatoria binaria de m bit-string y cromosomas de tamaño n
rand_population_binary(m, n) = [bitrand(n) for i in 1:m]

# Función que codifica las variables
length_variable(i_sup,i_inf,precision) = Int(ceil(log2((i_sup-i_inf)*10^(precision))))

# Función que obtiene las potencias en base dos de un vector de bits
to_decimal(dimension,v) = [2^(i-1) for i in 1:dimension]'*reverse(v)

# Función que codifica el vector de bits a un valor real
binary2real(i_sup,i_inf,dimension,pob) = [i_inf + (to_decimal(dimension,v)*(i_sup-i_inf)/(2^(dimension)-1)) for v in pob]

# Función que genera la estructura de datos Fenotipo
function DECODE(n_variables,m,i_sup_vec,i_inf_vec,dimension_vec,pob_vec)

    feno = Array{Int,2}(undef, m, 0)

    for i in 1:n_variables
        i_sup = i_sup_vec[i]
        i_inf = i_inf_vec[i]
        pob = pob_vec[i]
        dim = dimension_vec[i]
        feno = hcat(feno,binary2real(i_sup,i_inf,dim,pob))
    end

    return mapslices(x->[x],feno,dims=2)[:]
end

# Funcion que genera la estructura de datos de la función objetivo

function OBJFUN(f,feno)

    return f.(feno)

end

# Función que genera la aptitud de los individuos
function APTITUD(objv,operacion)

    val_max = maximum(objv)
    val_min = minimum(objv)

    if operacion == "min"
        objv_norm = [(((i-val_min)/(val_max-val_min))+0.1)^-1 for i in objv]
        suma = sum(objv_norm)
        key_objv = [(k,i/suma) for (k,i) in enumerate(objv_norm)]
        objv_sort = sort(key_objv, by = x -> x[2],rev=true)

    elseif operacion == "max"
        objv_norm = [(((i-val_min)/(val_max-val_min))+0.1) for i in objv]
        suma = sum(objv_norm)
        key_objv = [(k,i/suma) for (k,i) in enumerate(objv_norm)]
        objv_sort = sort(key_objv, by = x -> x[2],rev=true)

    end

    return objv_sort
end

# Función que selecciona a los mejores individuos

function SELECCION(aptitud,tipo,n_variables,población)
    if tipo == "ruleta"
        n = Int(length(aptitud)/2)
        suma_acumulada = cumsum([v for (k,v) in aptitud])

        individuos_dict = Dict(i => Dict() for i in 1:n)

        for pareja in 1:n
            for individuo in 1:2
                aleatorio = rand()
                index_ind = findall(x->x >= aleatorio , suma_acumulada)[1]
                cromosoma = []
                for gen in 1:n_variables
                    push!(cromosoma,poblacion_inicial_vec[gen][index_ind])
                end
                cromosoma = vcat(cromosoma...)
                individuos_dict[pareja][individuo] = cromosoma
            end
        end
    end

    return individuos_dict
end

function CRUZA(seleccion,tipo,length_total_cromosoma)
    if tipo == "unpunto"
        n = length(seleccion)

        nueva_poblacion = []

        for pareja in 1:n
            punto_cruza = Int(floor(rand(Uniform(1,length_total_cromosoma))))

            primer_nuevo_individuo = vcat([seleccion[pareja][1][1:punto_cruza],seleccion[pareja][2][punto_cruza+1:length_total_cromosoma]]...)
            segundo_nuevo_individuo = vcat([seleccion[pareja][2][1:punto_cruza],seleccion[pareja][1][punto_cruza+1:length_total_cromosoma]]...)

            push!(nueva_poblacion,primer_nuevo_individuo)
            push!(nueva_poblacion,segundo_nuevo_individuo)

        end

    end
    return nueva_poblacion

end

function MUTACION(nueva_poblacion,length_total_cromosoma)

    mutacion_param = 1/length_total_cromosoma
    n = length(nueva_poblacion)

    for individuo in 1:n
         muta_random = rand(length_total_cromosoma)
         muta_index = findall(x->x < mutacion_param ,muta_random)

         for i in muta_index
             nueva_poblacion[individuo][i] = !nueva_poblacion[individuo][i]
         end
    end

    return nueva_poblacion
end
