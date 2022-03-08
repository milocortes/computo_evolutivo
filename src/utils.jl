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
