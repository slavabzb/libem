function accum(m)
	tmp = 1
	for i in find(m)
		tmp *= abs(m[i])
	end
	tmp
end

function get_size(a, b, c)
	P = 1
	P *= accum(a)
	P *= accum(b)
	P *= accum(c)
	P = size(a)[1] * size(a)[2] + log10(abs(P))
	P = round(Int, P)
	println("size: ", P)
	P
end

function set_prec(n, L)
	println("prec: ", 2 * log2(n) + 2 * L)
end

function get_maxiter(b, c, n, R2, eps)
	p = sqrt(sum(x for x in b))
	m = sqrt(sum(x^2 for x in c))
	k = 2*(n+1)^2 * (abs(log(eps)) + log(1/p*m*R2))
	k = round(Int, k)
	k
end

a = Float64[1 0; 0 2; 3 2]
b = Float64[180, 150, 300]
c = Float64[-3, -5]

L = get_size(a, b, c)

n = size(a)[2]

set_prec(n, L)

x = zeros(2, 1)
R2 = Int128(n ^ 2) * (2 ^ (2 * L))

H = zeros(n, n)

for i in 1:n
	for j in 1:n
		if i == j
			H[i, j] = R2
		end
	end
end

eps = 1e-3
maxiter = get_maxiter(b, c, n, R2, eps)

function check(x, a, b)
	max = -Inf
	idx = 0
	for i in 1:size(a)[1]
		lhs = sum(a[i,:] .* x)
		rhs = b[i]
		#println(lhs, " >= ", rhs)
		if lhs >= rhs
			diff = lhs - rhs
			if diff > max
				max = diff
				idx = i
			end
		end
	end
	idx
end

for iter in 1:maxiter
	#println("iter: ", iter)
	
	g = c
	
	idx = check(x, a, b)
	if idx != 0
		g = a[idx,:]
	end
	
	#println("g: ", g)
	
	Hg = H * g
	#println("Hg: ", Hg)
	
	Hgg = sum(x for x in (Hg .* g))
	#println("Hgg: ", Hgg)
	
	if Hgg < 0
		print("iters: ", iter)
		break
	end
	
	x = x - 1 / ((n + 1) * sqrt(Hgg)) * Hg
	H = n^2/(n^2-1) * (H - 2/(n+1) * (Hg * (g' * H)) / Hgg)
	
	#println("x: ", x)
	#println("H: ", H)
	
	if iter == 2
		#break
	end
end

println("/", maxiter)

println("x: ", x)
println("f(x): ", c[1]*x[1] + c[2]*x[2])


