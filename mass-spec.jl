### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ b86c2020-c1ab-4957-86d0-021ca3155b6b
begin
	using Plots; gr()
	using PlutoUI
	import Contour: contours, levels, level, lines, coordinates
	using LaTeXStrings
	using Nabla
	using LinearAlgebra
end

# ╔═╡ 0fab30c2-c295-11eb-3da9-6360f936fdba
md"
# Modelling Ion Traps
"

# ╔═╡ 91d0ed27-bc9b-4bd1-bf16-be7c05810209
md"We can work with the following potential function which can be generated using a ring electrode and two cap electrodes:

$\Phi(x,y,z) = \dfrac{U_0}{r_0^2+2z_0^2}(2z^2 - x^2 - y^2)$

"

# ╔═╡ da94bfd8-22a9-4383-9a94-e7b3d47103de
Φ(x,y,z;U₀=0.2,r₀=1,z₀=1) = (U₀/((r₀^2)+2*(z₀^2)))*(2*(z^2) -(x^2)- (y^2))

# ╔═╡ a8a0f282-97b7-4f96-a9ee-ddb11571ed70
function slice(F, c=0;
                       xrng=(-1,1), yrng=xrng, zrng=xrng,
                       nlevels=6,         # number of levels in a direction
                       slices=Dict(:x => :black,
                                   :y => :black,
                                   :z => :black), # which directions and color
                       kwargs...          # passed to initial `plot` call
                       )

    _linspace(rng, n=150) = range(rng[1], stop=rng[2], length=n)

    X1, Y1, Z1 = _linspace(xrng), _linspace(yrng), _linspace(zrng)

    p = Plots.plot(;legend=false,kwargs...)

    if :x ∈ keys(slices)
        for x in _linspace(xrng, nlevels)
            local X1 = [F(x,y,z) for y in Y1, z in Z1]
            cnt = contours(Y1,Z1,X1, [c])
            for line in lines(levels(cnt)[1])
                ys, zs = coordinates(line) # coordinates of this line segment
                plot!(p, x .+ 0 * ys, ys, zs, color=slices[:x])
          end
        end
    end

    if :y ∈ keys(slices)
        for y in _linspace(yrng, nlevels)
            local Y1 = [F(x,y,z) for x in X1, z in Z1]
            cnt = contours(Z1,X1,Y1, [c])
            for line in lines(levels(cnt)[1])
                xs, zs = coordinates(line) # coordinates of this line segment
                plot!(p, xs, y .+ 0 * xs, zs, color=slices[:y])
            end
        end
    end

    if :z ∈ keys(slices)
        for z in _linspace(zrng, nlevels)
            local Z1 = [F(x, y, z) for x in X1, y in Y1]
            cnt = contours(X1, Y1, Z1, [c])
            for line in lines(levels(cnt)[1])
                xs, ys = coordinates(line) # coordinates of this line segment
                plot!(p, xs, ys, z .+ 0 * xs, color=slices[:z])
            end
        end
    end


    p
end

# ╔═╡ a8febe07-406a-40f1-8a1e-09562cbb4274
md"Let's try to visualise equipotential surfaces in 3D. Toggle the value of $\Phi$ and observe the equipotential surface corresponding to it below."

# ╔═╡ 306ebda7-53bb-4cda-b04d-8c318716e9d0
md" Φ = : $@bind c Slider(-0.12:0.01:0.12, default=-0.01, show_value=true)"

# ╔═╡ d45c9c1b-a022-40a6-8aaf-9bba825de77a
begin
	slice(Φ,c; nlevels=7,slices=Dict(:x=>:indigo,:y=>:indigo,:z=>:indigo))
	plot!(xlabel = L"x", ylabel= L"y", zlabel= L"z",
	xlim = (-1.1, 1.1),
    ylim = (-1.1, 1.1),
    zlim = (-1.1, 1.1),
	camera=(40,30),size=(400,300))
end

# ╔═╡ 5f071554-2fa7-4fed-886a-5f3380cb5b18
md"Lower values of the potential are obtained as you move radially away from the center and towards the bottom or top cap electrodes. It's clear that the potential is radially symmetric since we can also write the potential in cylindrical coordinates as follows:

$\phi(z,r) = \dfrac{U_0}{r_0^2+2z_0^2}(2z^2 - r^2)$

"

# ╔═╡ cba5969c-396d-4b17-94d3-eeda3fc9cfae
ϕ(z,r;U₀=0.2,r₀=1,z₀=1) = (U₀/((r₀^2)+2*(z₀^2)))*(2*(z^2) -(r^2))

# ╔═╡ e6c2c4a1-18e5-4e98-b815-24512652497d
begin	
	r=range(0,stop=1,length=100)
	z=range(-1,stop=1,length=100)
	plot(z, r, ϕ, st = :surface, colorbar=false,
	xlabel = L"z", ylabel= L"r", zlabel= L"\phi",zlim=(-0.12,0.15),
	xlim=(-1.1,1.1),ylim=(-0.1,1.1),
	camera=(40,30),size=(400,300))
end

# ╔═╡ 26295eb6-730c-4dda-b6ae-717c003a384e
md"The farther you move from the origin in the radial direction, the lower the potential energy. However, one can trap the particle along the axial direction since it will experience a parabolic potential."

# ╔═╡ d2cb4496-71d0-48d2-a353-178cadccdbe5
md" ## Penning Trap

To model a Penning Trap we need to introduce a magnetic field along the $z$-direction.
"

# ╔═╡ 1ff809d0-412c-447b-96ac-5dac9547adb3
md" B₀ = : $@bind B₀ Slider(0:0.1:1, default=0.6, show_value=true)"

# ╔═╡ b9168984-5503-42da-b316-a8989cbd5e86
md" m/q = : $@bind mbq Slider(-0.9:0.2:0.9, default=0.5, show_value=true)"

# ╔═╡ fb26007d-58c3-475c-ac9d-37c755a9cbda
∇Φ = ∇(Φ)

# ╔═╡ f2750a03-baa0-49d5-b644-c88dd298864a
F(x,y,z,V;B₀=0) = -collect(∇Φ(x,y,z)) + cross(V,[0,0,B₀])

# ╔═╡ bd5964f7-0aed-488e-a5cc-5a4d14a15d4a
X₀ = [0.0,0.0,-1.0] #initial position

# ╔═╡ 1c393adb-aa7b-431b-b28c-c54d89e46b85
V₀ = [0.01,0.01,0.1] #initial velocity

# ╔═╡ 05ad4fc6-f1b7-4a97-a308-4ecc077a0537
T = 20 #total time of evolving the system

# ╔═╡ 7ffcbfe8-c097-4fe0-9c26-05ed44092d93
let
	dt = 0.1
	
	Xₜ = X₀
	Vₜ = V₀
	
	plt = plot3d(1,
    xlim = (-1.1, 1.1),
    ylim = (-1.1, 1.1),
    zlim = (-1.1, 1.1),
    legend=false,
	xlabel = L"x", ylabel= L"y", zlabel= L"z",
	camera=(40,30),size=(400,300))
	
	@gif for t in 0:dt:T
		Vₜ += F(Xₜ[1],Xₜ[2],Xₜ[3],Vₜ;B₀=B₀)/mbq*dt
		Xₜ += Vₜ*dt
		push!(plt, Xₜ[1],Xₜ[2],Xₜ[3])
	end
end

# ╔═╡ 6cc49635-7920-4200-8064-87c49adb0e28
md"The larger the magnetic field strength, the heavier the ion that can be _trapped_. For mass ratios greater than this special ratio, the ion will follow the electric potential gradient and fall into one of the caps while for mass ratios lower than it, the ions spiral out into the wall. To trap -ve ions, the potential configuration must be reversed.
"

# ╔═╡ 88d7df5b-7409-4775-a563-d45cd69d8b2f
md"## Paul Trap
In a Paul Trap, instead of applying a constant magnetic field in the $z$-direction, a sinusoidally varying electric potential is applied —

$\Phi_p(x,y,z,t) = \dfrac{U_0-U_1cos(\Omega t)}{r_0^2+2z_0^2}(2z^2 - x^2 - y^2)$

or equivalently

$\phi_p(z,r,t) = \dfrac{U_0-U_1cos(\Omega t)}{r_0^2+2z_0^2}(2z^2 - r^2)$

"

# ╔═╡ f0b5768d-23d1-4b8a-8c61-bf6122f28465
Ω = 5 #fixed

# ╔═╡ b3d7a1e7-03d9-450e-a223-2c8a8ab64393
md"Let's observe how the potential function changes in time."

# ╔═╡ 56626c23-8b67-40d4-934d-7a8246f2da8e
md" U₁ = : $@bind U₁ Slider(0:0.5:4, default=3, show_value=true)"

# ╔═╡ 80c8d186-f6de-4b30-8c05-bdd330c6b717
Φₚ(x,y,z;t,U₀=0.2,r₀=1,z₀=1) = ((U₀-U₁*cos(Ω*t))/((r₀^2)+2*(z₀^2)))*(2*(z^2) -(x^2)- (y^2))

# ╔═╡ baf9f7e2-1fc4-4254-acfe-a61877ccd8e6
ϕₚ(z,r;t,U₀=0.2,r₀=1,z₀=1) = ((U₀-U₁*cos(Ω*t))/((r₀^2)+2*(z₀^2)))*(2*(z^2) -(r^2))

# ╔═╡ 95bde903-01fb-460b-a368-fc96b0061abe
md" m/q = : $@bind MbQ Slider(-1:0.4:1, default=0.6, show_value=true)"

# ╔═╡ 77092aad-fedd-4325-ae6e-dccfef3b4fce
md"The potential flips in a way that the trap centre is a minima along the radial axis for half the time period and a maxima for the other half."

# ╔═╡ 17b25b04-9af4-4794-a01e-322b649d261e
md"Similar to a Penning Trap, the greater the magnitude of the AC potential, the heavier the ion that can be trapped. Ions with mass ratios lower than this special ratio oscillate wildly in the $z$-direction and get ejected out of the trap, while ions with with mass ratios greater than it follow the DC potential to their doom. 

Unlike a Penning Trap, a Paul Trap can stabilise ions irrespective of the sign of their charge without changing the potential configuration. However, for $U_0 \neq 0$ the effect on oppositely charged ions is not symmetric, i.e. for the same value of $U_1$ a -ve ion might be stable while a +ve one with same $|q|/m$ might not be.

It is also possiible to make the trap function as a mass spectrometer by sweeping across increasing values of $U_1$ for $U_0=0$. The lighter ions will get ejected in the $z$-direction first."

# ╔═╡ 19a67711-8435-4645-bfff-729d4754d3e0
md" U₀ = : $@bind U₀ Slider(0:0.2:0.2, default=0.2, show_value=true)"

# ╔═╡ 62cc4d48-93f2-41dd-864d-cf5beb924ee4
let
	dt = 0.1
	@gif for t in 0:dt:2*π/Ω
		r=range(0,stop=1,length=100)
		z=range(-1,stop=1,length=100)
		ϕₚₜ(z,r) = ϕₚ.(z,r;t=t,U₀=U₀)
		plot(z, r, ϕₚₜ, st = :surface, colorbar=false,
		xlabel = L"z", ylabel= L"r", zlabel= L"\phi",zlim=(-4,4),
		xlim=(-1.1,1.1),ylim=(-0.1,1.1),
		camera=(40,30),size=(400,300))
	end
end

# ╔═╡ 2da65d77-dbbf-44d6-bced-52df23a6aa77
let
	dt = 0.1
	
	Xₜ = X₀
	Vₜ = V₀
	
	plt = plot3d(1,
    xlim = (-1.1, 1.1),
    ylim = (-1.1, 1.1),
    zlim = (-1.1, 1.1),
    legend=false,
	xlabel = L"x", ylabel= L"y", zlabel= L"z",
	camera=(40,30),size=(400,300))
	
	@gif for t in 0:dt:T
		Φₚₜ(x,y,z) = Φₚ.(x,y,z;t=t,U₀=U₀)
		Fₚ = -collect(∇(Φₚₜ)(Xₜ[1],Xₜ[2],Xₜ[3]))
		Vₜ += Fₚ/MbQ*dt
		Xₜ += Vₜ*dt
		push!(plt, Xₜ[1],Xₜ[2],Xₜ[3])
	end
end

# ╔═╡ 106e2cc0-6ddf-4559-8ba2-70c9c72af879
md"## Quadrupole Potential
A slight simplification of the Paul Trap potential gives us the quadrupole potential which is commonly used in quadrupole mass analysers: 

$\phi_4(x,y) = \dfrac{U_0-U_1cos(\Omega t)}{2r_0^2}(x^2 - y^2)$

The quadrupole potential has no $z$-dependence."

# ╔═╡ 66feb652-b667-4ac3-a374-5feb2a75ef70
md"We can visualise the equipotential surfaces at $t=\pi/2\Omega$ (i.e. after half a cycle) to get a sense of what kind of shape the electrodes might have in this setup."

# ╔═╡ 5cf37cfe-ea52-46f6-b1ff-9f16db5fc842
md" ϕ₄(t=π/2Ω) = : $@bind p Slider(-0.2:0.05:0.2, default=0.05, show_value=true)"

# ╔═╡ 1de3f53e-8229-48b1-8c4d-c978c02069f7
md"We can also see how ions with different charge-to-mass ratios will behave under different magnitudes of alternating potential."

# ╔═╡ 9d61ea5b-898b-4529-9f66-3fdbeca5beed
md" U₁ = : $@bind U₁₄ Slider(2:2:6, default=4, show_value=true)"

# ╔═╡ 2a931bd7-b868-4b77-bc2e-2eb18d4d1ef2
ϕ₄(x,y;t,U₀=0.5,r₀=1) = ((U₀-U₁₄*cos(Ω*t))/(2*(r₀^2)))*((x^2) -(y^2))

# ╔═╡ 7cb416aa-8944-48f4-8f84-6dc316f1faea
begin
	Φ₄(x,y,z) = ϕ₄(x,y;t=π/(2*Ω))
	slice(Φ₄,p; nlevels=7,slices=Dict(:x=>:indigo,:y=>:indigo,:z=>:indigo))
	plot!(xlabel = L"x", ylabel= L"y", zlabel= L"z",
	xlim = (-1.1, 1.1),
    ylim = (-1.1, 1.1),
    zlim = (-1.1, 1.1),
	camera=(40,30),size=(400,300))
end

# ╔═╡ ea992085-7822-4ac2-8251-49687d6483cc
md" m/q = : $@bind MbQ₄ Slider(-0.9:0.2:0.9, default=0.5, show_value=true)"

# ╔═╡ 4aed0cd0-ba1c-452c-8ebd-1ab69961ae8b
let
	dt = 0.1
	@gif for t in 0:dt:2*π/Ω
		x=range(-1,stop=1,length=100)
		y=range(-1,stop=1,length=100)
		ϕ₄ₜ(x,y) = ϕ₄.(x,y;t=t)
		plot(x, y, ϕ₄ₜ, st = :surface, colorbar=false,
		xlabel = L"x", ylabel= L"y", zlabel= L"\phi",zlim=(-8,8),
		xlim=(-1.1,1.1),ylim=(-1.1,1.1),
		camera=(40,50),size=(400,300))
	end
end

# ╔═╡ a45cfe79-188b-437d-87c4-3fa73e47b5c5
md"There is no confinement along the $z$-direction since we don't want the ions to sit in the centre but pass through the quadrupole if they fit our charge-to-mass ratio criteria."

# ╔═╡ 33a93d3d-339e-4ac8-903a-c3abbd513944
let
	dt = 0.1
	
	Xₜ = X₀
	Vₜ = V₀
	
	plt = plot3d(1,
    xlim = (-1.1, 1.1),
    ylim = (-1.1, 1.1),
    zlim = (-1.1, 1.1),
    legend=false,
	xlabel = L"x", ylabel= L"y", zlabel= L"z",
	camera=(40,30),size=(400,300))
	
	@gif for t in 0:dt:T
		ϕ₄ₜ(x,y) = ϕ₄.(x,y;t=t)
		Fₚ = [-collect(∇(ϕ₄ₜ)(Xₜ[1],Xₜ[2])); 0]
		Vₜ += Fₚ/MbQ₄*dt
		Xₜ += Vₜ*dt
		push!(plt, Xₜ[1],Xₜ[2],Xₜ[3])
	end
end

# ╔═╡ 3ff933fd-b807-4797-a098-3f7bf3322503
md"Note how heavier ions follow the DC potential and deviate from their trajectories. Lighter ions have highly unstable trajectories and crash into the walls. Only for a small intermediate range of $m/q$ can ions pass through a quadrupole. A quadrupole can filter ions with particular mass-to-charge ratios irrespective of sign."

# ╔═╡ Cell order:
# ╟─0fab30c2-c295-11eb-3da9-6360f936fdba
# ╟─b86c2020-c1ab-4957-86d0-021ca3155b6b
# ╟─91d0ed27-bc9b-4bd1-bf16-be7c05810209
# ╟─da94bfd8-22a9-4383-9a94-e7b3d47103de
# ╟─a8a0f282-97b7-4f96-a9ee-ddb11571ed70
# ╟─a8febe07-406a-40f1-8a1e-09562cbb4274
# ╟─306ebda7-53bb-4cda-b04d-8c318716e9d0
# ╟─d45c9c1b-a022-40a6-8aaf-9bba825de77a
# ╟─5f071554-2fa7-4fed-886a-5f3380cb5b18
# ╟─cba5969c-396d-4b17-94d3-eeda3fc9cfae
# ╟─e6c2c4a1-18e5-4e98-b815-24512652497d
# ╟─26295eb6-730c-4dda-b6ae-717c003a384e
# ╟─d2cb4496-71d0-48d2-a353-178cadccdbe5
# ╟─1ff809d0-412c-447b-96ac-5dac9547adb3
# ╟─b9168984-5503-42da-b316-a8989cbd5e86
# ╟─fb26007d-58c3-475c-ac9d-37c755a9cbda
# ╟─f2750a03-baa0-49d5-b644-c88dd298864a
# ╟─bd5964f7-0aed-488e-a5cc-5a4d14a15d4a
# ╟─1c393adb-aa7b-431b-b28c-c54d89e46b85
# ╟─05ad4fc6-f1b7-4a97-a308-4ecc077a0537
# ╟─7ffcbfe8-c097-4fe0-9c26-05ed44092d93
# ╟─6cc49635-7920-4200-8064-87c49adb0e28
# ╟─88d7df5b-7409-4775-a563-d45cd69d8b2f
# ╟─80c8d186-f6de-4b30-8c05-bdd330c6b717
# ╟─baf9f7e2-1fc4-4254-acfe-a61877ccd8e6
# ╟─f0b5768d-23d1-4b8a-8c61-bf6122f28465
# ╟─b3d7a1e7-03d9-450e-a223-2c8a8ab64393
# ╟─56626c23-8b67-40d4-934d-7a8246f2da8e
# ╟─95bde903-01fb-460b-a368-fc96b0061abe
# ╟─62cc4d48-93f2-41dd-864d-cf5beb924ee4
# ╟─77092aad-fedd-4325-ae6e-dccfef3b4fce
# ╟─2da65d77-dbbf-44d6-bced-52df23a6aa77
# ╟─17b25b04-9af4-4794-a01e-322b649d261e
# ╟─19a67711-8435-4645-bfff-729d4754d3e0
# ╟─106e2cc0-6ddf-4559-8ba2-70c9c72af879
# ╟─2a931bd7-b868-4b77-bc2e-2eb18d4d1ef2
# ╟─66feb652-b667-4ac3-a374-5feb2a75ef70
# ╟─5cf37cfe-ea52-46f6-b1ff-9f16db5fc842
# ╟─7cb416aa-8944-48f4-8f84-6dc316f1faea
# ╟─1de3f53e-8229-48b1-8c4d-c978c02069f7
# ╟─9d61ea5b-898b-4529-9f66-3fdbeca5beed
# ╟─ea992085-7822-4ac2-8251-49687d6483cc
# ╟─4aed0cd0-ba1c-452c-8ebd-1ab69961ae8b
# ╟─a45cfe79-188b-437d-87c4-3fa73e47b5c5
# ╟─33a93d3d-339e-4ac8-903a-c3abbd513944
# ╟─3ff933fd-b807-4797-a098-3f7bf3322503
