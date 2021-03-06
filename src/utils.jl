

"""
Bequest function after DeNardi (2004)
"""
function bequest(b::Float64, p::Param)
	@assert b >= 0.0  # no negative bequests
	p.ν * u(b + p.bbar,p)
end

function bequest(b::Vector{Float64}, p::Param)
	n = length(b)
	y = similar(b)
	for i in 1:n
		y[i] = bequest(b[i],p)
	end
	y
end

#utility functions

# utility without discrete choice
function u(x::Float64,p::Param)
	if p.gamma == 1.0
		log(x)
	else
		p.oneover_oneminusgamma * (x^p.oneminusgamma)
	end
end
function u(x::Array{T}, p::Param) where T
	n = length(x)
	y = zeros(T,n)
	for i in 1:n
		y[i] = u(x[i],p)
	end
	y
end

# utility with discrete choice
function u(x::Float64,working::Bool,p::Param)
	if p.gamma == 1.0
		log(x) - p.alpha*working
	else
		p.oneover_oneminusgamma * (x^p.oneminusgamma - 1.0) - p.alpha*working
	end
end
function u(x::Array{T}, working::Bool, p::Param) where T
	n = length(x)
	y = zeros(T,n)
	for i in 1:n
 		y[i] = u(x[i],working,p)
	end
	y
end

"bequest with discrete choice"
function bequest(b::Float64, working::Bool, p::Param)
	bequest(b,p) - p.alpha * working
end

function bequest(b::Vector{Float64}, working::Bool, p::Param)
	n = length(b)
	y = similar(b)
	for i in 1:n
		y[i] = bequest(b[i],working,p)
	end
	y
end


# partial derivative of utility wrt c
function up(c::Float64,p::Param)
	if p.gamma == 1.0
		1.0 / c
	else
		c ^ (p.neg_gamma)
	end
end
function up(c::Array{Float64,2},p::Param)
	n = length(c)
	x = zeros(size(c))
	for i in 1:n
		x[i] = up(c[i],p)
	end
	x
end
function up!(c::Array{Float64},p::Param)
	for i in eachindex(c)
		c[i] = up(c[i],p)
	end
end
function up(c::Array{Float64},p::Param)
	n = length(c)
	x = zeros(n)
	for i in 1:n
		x[i] = up(c[i],p)
	end
	x
end

# inverse of partial derivative
function iup(u::Float64,p::Param)
	if p.gamma == 1.0
		1.0 / u
	else
		u ^ p.neg_oneover_gamma
	end
end
function iup(u::Vector{Float64},p::Param)
	n = length(u)
	x = zeros(n)
	for i in 1:n
		x[i] = iup(u[i],p)
	end
	return x
end

NBLt(ylow,t::Int,p::Param) = sum(-p.R^(-(j-1) ) * income(j,p,ylow) for j in t:p.nT-1)

nbl(ylow,it,p::Param) = sum(-p.R^(j) * income(it+j,p,ylow) for j in 1:((p.nT)-1)-it )
nblv(ylow,it,p::Param) = [-p.R^(j) * income(it+j,p,ylow) for j in 1:((p.nT)-1)-it ]

# NBLt(ylow,t::Int,p::Param) = sum(-p.R^(-j ) * income(j,p,ylow) for j in (p.nT-1)-t:-1:0)
NBLs(ylow,p::Param) = [-p.R^(-j) * income(j,p,ylow) for j in (p.nT-1)-1:-1:0]


# nextbound(ylow,it,p::Param) = sum(-income(j,p,ylow)*p.R^(-(j-it)) for j in (it+1):(p.nT-1))
nextbound(ylow,it,p::Param) = -income(it+1,p,ylow)*p.R^(-1)
pnextbound(ylow,it,p::Param) = ["income($j,p,ylow)*p.R^($(-(j-it)))" for j in (it+1):(p.nT-1)]




"""
	abounds(ylow,p::Param)

returns a vector with the lowest possible asset level in each period.
This assumes the worst case of income draw in each period and zero consumption, while ``a_T=0``.
"""
function abounds(ylow,p::Param)
	x = [nextbound(ylow,it,p) for it in 1:p.nT-2]
	prepend!(x, -ylow + x[1]) # adds period 1 bound
	x
	# n = NBLs(ylow,p)  # vector of borrowing limit in period t
	# reverse(cumsum(reverse(n)))
end




allincomes(m::Model,p::Param) = [DCEGM.income(j,p,y) for y in m.yvec, j in 1:p.nT]

"expected borrowing limit"
function EBL(m::Model,p::Param)
	yy = [DCEGM.income(j,p,y) for y in m.yvec, j in 1:p.nT]  # income at each age and each state
	r = similar(yy)
	r[:,p.nT] = yy[:,p.nT]

	for it in p.nT-1:-1:1
		r[:,it] = m.ywgt * yy[:,it+1]
	end
	return r
end







"""
rouwenhorst AR1 approximation


This is taken from [http://karenkopecky.net/RouwenhorstPaperFinal.pdf](Karen Kopecky's paper)

"""
function rouwenhorst(rho::Float64,mu_eps,sigma_eps,n)
	q = (rho+1)/2
	nu = ((n-1)/(1-rho^2))^(1/2) * sigma_eps
	P = [q 1-q ; 1-q q]

	for i=2:n-1

		# P = q * vcat(hcat(P , zeros(i,1)),zeros(1,i+1)) .+ (1-q).* vcat( hcat(zeros(i,1),P), zeros(1,i+1)) .+
		# (1-q) .* vcat(zeros(1,i+1),hcat(P,zeros(i,1))) .+ q .*vcat(zeros(1,i+1),hcat(zeros(i,1),P))
		# P[2:i,:] = P[2:i,:] ./ 2

		P = q .* [P zeros(i,1); zeros(1,i+1)] .+ (1-q) .* [zeros(i,1) P; zeros(1,i+1)] .+
		(1-q) .* [zeros(1,i+1); P zeros(i,1)] .+   q   .* [zeros(1,i+1); zeros(i,1) P]

		P[2:i,:] = P[2:i,:] ./ 2
	end

	z = range(mu_eps/(1-rho)-nu,stop = mu_eps/(1-rho)+nu,length = n);
	return (collect(z),P)
end

# asset grid scaling
function scaleGrid(lb::Float64,ub::Float64,n::Int;logorder::Int=1,offset=0.05)
	if lb >= ub
		@warn("lower bound $lb is greater than upper bound $ub: forcing lb < ub.")
		lb = ub - abs(ub)*offset
	end
	out = zeros(n)
	if logorder==0
		out    = collect(range(lb,stop = ub,length = n))
	elseif logorder==1
		off = 1
		if lb<0
			off = 1 - lb #  adjust in case of neg bound
		end
		out[1] = log(lb + off)
		out[n] = log(ub + off)
		out    = collect(range(out[1],stop = out[n],length = n))
		out    = exp.(out) .- off
	elseif logorder==2
		off = 1
		if lb<0
			off = 1 - lb #  adjust in case of neg bound
		end
		out[1] = log( log(lb + off) + off )
		out[n] = log( log(ub + off) + off )
		out    = collect(range(out[1],stop = out[n],length = n))
		out    = exp.( exp.(out) .- off ) .- off
	elseif logorder==3
		off = 1
		if lb<0
			off = 1 - lb #  adjust in case of neg bound
		end
		out[1] = log( log( log(lb + off) + off ) + off )
		out[n] = log( log( log(ub + off) + off ) + off )
		out    = collect(range(out[1],stop = out[n],length = n))
		out    = exp.( exp.( exp.(out) .- off ) .- off ) .- off
	else
		error("supports only up to tripple log grid")
	end
end

"""
lifecycle profile in income
"""
function income(it::Int,p::Param,shock::Float64)
	age = it + 19
	if it > p.retage
		p.pension
	else
		exp( p.inc0 + p.inc1*age - p.inc2*(age^2) + shock)
	end
end
function income(it::Int,p::Param,shock::Array{Float64})
	x = similar(shock)
	for i=1:length(x)
		x[i] = income(it,p,shock[i])
	end
	return x
end



function quadpoints(n,lbnd,ubnd)

   x2  = ubnd
   x1  = lbnd
   x   = zeros(n)
   w   = zeros(n)
   EPS = 3.e-14
   m   = floor((n+1)/2)
   xm  = (x2+x1)/2
   xl  = (x2-x1)/2
   i = 1
   z1 = 1.e99
   pp = p1 = p2 = p3 = z = 0.0
   while (i <= m)
	   z  = cos(pi*(i-0.25)/(n+0.5))
	   while (abs(z-z1)>EPS)
	       p1 = 1
	       p2 = 0
	       j=1
	       while (j <= n)
				 p3 = copy(p2)
				 p2 = copy(p1)
				 p1 = ((2*j-1)*z*p2-(j-1)*p3)/j
		         j=j+1
	       end
	       pp = n*(z*p1-p2)/(z*z-1)
	       z1 = z
	       z  = z1 - p1/pp
    	end
     x[i]     = xm - xl*z
     x[n+1-i] = xm + xl*z
     w[i]     = 2*xl/((1-z*z)*pp*pp)
     w[n+1-i] =w[i]
     i = i+1
   end
   return (x,w)
end


function linearapprox(x::Vector{Float64},y::Vector{Float64},xi::Float64,lo::Int,hi::Int)
	r = 0.0
	n = length(x)
	@assert n==length(y)

	# determining bounds
	if xi == x[1]
		r = y[1]
		return r
	elseif xi < x[1]
		# get linear approx below
		@inbounds r = y[1] + (y[2] - y[1]) * (xi - x[1])  / (x[2] - x[1])
		return r
	end
	if xi == x[n]
		r = y[n]
		return (r,n)
	elseif xi > x[n]
		# get linear approx above
		@inbounds r = y[n] + (y[n] - y[n-1]) * (xi - x[n])  / (x[n] - x[n-1])
		return r
	end

	# if have to find interval
	if hi - lo > 1
		jinf = searchsortedlast(x,xi,lo,hi,Base.Forward)	# get rid
	# if not, lo is jinf
	else
		jinf = lo
	end
	@inbounds r = (y[jinf] * (x[jinf+1] - xi) + y[jinf+1] * (xi - x[jinf]) ) / (x[jinf+1] - x[jinf])
	return r
end
linearapprox(x::Vector{Float64},y::Vector{Float64},xi::Float64) = linearapprox(x,y,xi,1,length(x))

function linearapprox(x::Vector{Float64},y::Vector{Float64},xi::Vector{Float64})
	n = length(xi)
	z = zeros(n)
	for i in 1:n
		z[i] = linearapprox(x,y,xi[i])
	end
	return z
end
