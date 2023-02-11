### A Pluto.jl notebook ###
# v0.19.16

using Markdown
using InteractiveUtils

# ╔═╡ 516346f0-43c8-4c43-af5f-566d56fb0a7a
begin
	using LinearAlgebra
	using StarAlgebras
	using JuMP
	using Pkg
	# Pkg.resolve()
end

# ╔═╡ e59da2e3-f60f-4550-be8e-a9035e2f19eb
begin
	Pkg.build("MATLAB");
	using MATLAB  # A working MATLAB installation is required
	cd("SeDuMi_1_3/") do  # Download binaries from sedumi.ie.lehigh.edu/?page_id=58 and move them to the given folder. In this case it's the project folder, "SeDuMi_1_3/".
		   mat"install_sedumi"
	end;
	using SeDuMi
end

# ╔═╡ dbe625c4-d474-4c3d-bb01-e35912d5dad9
md"""
## Proving Kazhdan's Property (T) for ``SL(3, \mathbb{Z})``
"""

# ╔═╡ fb061d85-1a5f-439a-81d8-78606c005586
md"""
In this notebook, we reproduce the proof of Property (T) for ``SL(3, \mathbb{Z})`` as presented in the paper "Kazhdan's Property (T) via Semidefinite Optimization" by Tim Netzer and Andreas Thom, published [here](https://arxiv.org/abs/1411.2488).
"""

# ╔═╡ 617a09fb-c847-467d-93b9-30a6aed3e058
begin
	begin
		As = -[
		            0.0 0.0 0.0 1.0 0.0
		            1.0 0.0 1.0 0.0 0.0
		        ]
		        bs = [1.0, 0.0]
		        cs = [0.0, 1.0, 0.0, 0.0, 1.0]
		        primal, dual, info =
		            sedumi(As, bs, cs, SeDuMi.Cone(0, 1, [], [], [2]), fid = 0)
	end;
	md"""This cell tests whether MATLAB and SeDuMi are installed correctly."""
end

# ╔═╡ 7affa801-20e5-4d2b-9637-dab2e35b3269
begin
	begin
		ℤ = BigInt
		
		function generatingset_SLnZ(n::Integer = 3)
			S = Set{Matrix{ℤ}}()
			for i ∈ 1:n
				for j ∈ 1:n
					if i ≠ j
						Mᵢⱼ = Matrix{ℤ}(I, n, n)
						Mᵢⱼ[i, j] = 1
						push!(S, Mᵢⱼ)
						Mᵢⱼ = Matrix{ℤ}(I, n, n)
						Mᵢⱼ[i, j] = -1
						push!(S, Mᵢⱼ)
					end
				end
			end
			return S
		end

		Σ = sum

		is_positive_semi_definite(M::Matrix) = all(map(x -> x ≥ 0, eigvals(M)))
	end;
	md"""Set up auxilliary functions and definitions."""
end

# ╔═╡ f7e30187-e485-421f-a9d7-fa14a33c867c
S = generatingset_SLnZ(3)

# ╔═╡ d04b139a-08e9-438c-8393-5cf91424300d
@assert length(S) == 12

# ╔═╡ e562e9e4-a583-11ed-360d-2daaf7526351
md"""
``G = SL(3, \mathbb{Z}) = \langle S \rangle``.
"""

# ╔═╡ 7dae365e-c45e-48b0-8720-9eb6ecb8dc90
md"""
We define ``M`` to be the set of all products of elements of ``S`` of length ``\leq 4``.
"""

# ╔═╡ f40a2b8c-4dde-4fbb-aeab-4337bc2e7c53


# ╔═╡ 6351206c-8699-4a28-bf3c-6d4335412bd7
md"""
Let ``P = (p_{ij})``. The equation
```math
a = (a^*_1,\dots,a^*_m)P(a_1,\dots,a_m)^T
```
yields the following conditions

"""

# ╔═╡ 02514078-89bc-4c1d-9233-8442ece8a4b9
md"""
some conditions go here
"""

# ╔═╡ 003cc82c-bac7-4856-8051-f32bf3f97e49
begin
	abstract type Group end

	mutable struct MatrixGroup{T} <: Group
	    S::Vector{Matrix{T}}
	
	    function MatrixGroup(S::Vector{Matrix{T}}) where {T}
	        new{T}(S)
	    end
	end
	
	Base.one(G::MatrixGroup{T}) where {T} = one(first(G.S))
end

# ╔═╡ cfc39d8f-231e-4a7b-baf2-f94515e98b83
one(first(S));

# ╔═╡ cd42c36d-0f38-4169-9ada-34bb401ecd94
[[one(first(S))]; collect(S)]

# ╔═╡ e83d5691-34d7-4c9e-b94e-fd695e7ed731
M, sizes = begin
	ball = Set([[one(first(S))]; collect(S)])
	let S = S
		for Mᵢⱼ ∈ S
			for Nᵢⱼ ∈ S
				push!(ball, Mᵢⱼ * Nᵢⱼ)
				push!(ball, Nᵢⱼ * Mᵢⱼ)
			end
		end
	end
	sizes = [length(ball)]
	
	let S = copy(ball)
		for Mᵢⱼ ∈ S
			# push!(M, Mᵢⱼ)
			for Nᵢⱼ ∈ S
				push!(ball, Mᵢⱼ * Nᵢⱼ)
				push!(ball, Nᵢⱼ * Mᵢⱼ)
			end
		end
	end
	push!(sizes, length(ball))
	ball, sizes
end

# ╔═╡ ba85cb87-2a25-46e4-a42d-f2018e9df332
begin
	@assert length(M) == 5455
	md"""M has $(length(M)) elements"""
end

# ╔═╡ 04c3b592-8264-4ec7-a066-01ca33e6f1f7
begin
	G = MatrixGroup(collect(M))
	basis = StarAlgebras.Basis{Float64}(collect(M))
	ℝG = StarAlgebra(MGroup, basis)
	𝒜 = ℝG
end

# ╔═╡ 29cd69ae-4d9c-4533-b90b-38659c6ceee5
3.141 * one(ℝG)

# ╔═╡ 1a05667f-b42d-40b5-89b1-3340f2f25ca9
begin
	Δ = length(S) * one(𝒜)
	for s ∈ S
		Δ -= 𝒜(s)
	end
	Δ
end

# ╔═╡ b1443e53-69eb-4417-af10-eb71b8ebf837
ε = 0.2805

# ╔═╡ 7fbba67d-580a-4033-8425-f5a2e13b7307
md"""
Let's define the unnormalized Laplace operator:
"""

# ╔═╡ 50a7ee5b-2f09-4744-ae57-87affde137bf
parent(Δ).mstructure.basis[1]

# ╔═╡ 2a66e100-6fd5-4ef7-b132-17811d81e74b
parent(Δ).mstructure.basis[2]

# ╔═╡ 153766c0-2c4d-4ad1-b3a8-ca3827932a60


# ╔═╡ 1be46ea1-2d3c-45a3-8ddd-52e684f874c0
a = Δ * Δ

# ╔═╡ 2578cf1d-f9d0-404a-a37b-c893c995bd5f
begin
	model = Model(SeDuMi.Optimizer)
	@variable(model, P[1:121, 1:121], PSD)
end

# ╔═╡ cae53bd5-ba53-4591-9ca7-40278746729c


# ╔═╡ 0b62796b-824a-4c86-9396-175d901d0ba2


# ╔═╡ f066b06b-08a4-445b-9b37-a965edb4a6f0
md"""
Then we compute the root:
"""

# ╔═╡ d3b9af2a-8677-40fc-b7a5-2ff498ac3256
√(first(S))

# ╔═╡ 8a2bdbc4-7151-4a1a-b7a0-b57fa5f94903
md"""
We change it slightly to make all rows sum to 0
"""

# ╔═╡ fb357b59-f420-4505-8fdc-75e77201b843
md"""
Let ``P = 10^{-12}Q^TQ``
"""

# ╔═╡ c458d180-005f-4c6c-a0a7-8374d7c03b3c
Q = Matrix(I, 3, 3)

# ╔═╡ 0147d4c2-a7d9-45d4-bdb7-8219423b9ebc
P´ = 10^-12 * Q' * Q

# ╔═╡ bd554b57-57f9-42ec-9957-7efb6b3a3538
@assert is_positive_semi_definite(P)

# ╔═╡ 5ffd0d63-82ec-4514-8389-26e2317abc66
b = star(a) * P * transpose(a)

# ╔═╡ 96157e19-0938-4b69-98cc-e50492f6a518
c = a - b

# ╔═╡ 9f050d68-ee64-4eff-8f27-efd8faeca574
md"""
### Lemma 2.1.
Let G be a group with finite generating set ``S = S^{-1}``, and let
```math
	\Delta = |S| - \sum_{s \in S}s \in \mathbb{R}[G]
```
be the Laplace operator. Let ``c = \sum_gc_gg \in \omega[G]^h`` be such that whenever ``c_g \neq 0``, then g is a product of at most ``2^d`` elements from ``S``. Then
```math
	c + 2^{2d-1}||c||_1 \cdot \Delta \in \Sigma^2\mathbb{R},
```
where ``||c||_1 = \sum_g|c_g|``.
"""

# ╔═╡ 27e8136b-b62e-4969-824d-73c6059f6447
md"""
Sources and references:

	* Kazhdan's Property (T) via Semidefinite Optimization [arxiv](arxiv.org/abs/1411.2488)
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
JuMP = "4076af6c-e467-56ae-b986-b466b2749572"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
MATLAB = "10e44e05-a98a-55b3-a45b-ba969058deb6"
Pkg = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
SeDuMi = "a895aaad-f784-5544-9392-bb281339c1b2"
StarAlgebras = "0c0c59c1-dc5f-42e9-9a8b-b5dc384a6cd1"

[compat]
JuMP = "~1.7.0"
MATLAB = "~0.8.3"
SeDuMi = "~0.4.1"
StarAlgebras = "~0.1.8"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.5"
manifest_format = "2.0"
project_hash = "e61d8441d2f8daf01f887af25e54922e7c433e94"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BenchmarkTools]]
deps = ["JSON", "Logging", "Printf", "Profile", "Statistics", "UUIDs"]
git-tree-sha1 = "d9a9701b899b30332bbcb3e1679c41cce81fb0e8"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "1.3.2"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c6d890a52d2c4d55d326439580c3b8d0875a77d9"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.7"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "844b061c104c408b24537482469400af6075aae4"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.5"

[[deps.CodecBzip2]]
deps = ["Bzip2_jll", "Libdl", "TranscodingStreams"]
git-tree-sha1 = "2e62a725210ce3c3c2e1a3080190e7ca491f18d7"
uuid = "523fee87-0ab8-5b00-afb7-3ecf72e48cfd"
version = "0.7.2"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "9c209fb7536406834aa938fb149964b985de6c83"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.1"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "61fdd77467a5c3ad071ef8277ac6bd6af7dd4c04"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.6.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.1+0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "c5b6685d53f933c11404a3ae9822afe30d522494"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.12.2"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "a69dd6db8a809f78846ff259298678f0d6212180"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.34"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "49510dfcb407e572524ba94aeae2fced1f3feb0f"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.8"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.JuMP]]
deps = ["LinearAlgebra", "MathOptInterface", "MutableArithmetics", "OrderedCollections", "Printf", "SnoopPrecompile", "SparseArrays"]
git-tree-sha1 = "8ebcb407c28617ea075563c550ec766dddf25a2e"
uuid = "4076af6c-e467-56ae-b986-b466b2749572"
version = "1.7.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "680e733c3a0a9cea9e935c8c2184aea6a63fa0b5"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.21"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MATLAB]]
deps = ["Libdl", "SparseArrays"]
git-tree-sha1 = "e263657fe013cb02450c5d4210d2c50a354a5e08"
uuid = "10e44e05-a98a-55b3-a45b-ba969058deb6"
version = "0.8.3"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "42324d08725e200c23d4dfb549e0d5d89dede2d2"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.10"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MathOptInterface]]
deps = ["BenchmarkTools", "CodecBzip2", "CodecZlib", "DataStructures", "ForwardDiff", "JSON", "LinearAlgebra", "MutableArithmetics", "NaNMath", "OrderedCollections", "Printf", "SnoopPrecompile", "SparseArrays", "SpecialFunctions", "Test", "Unicode"]
git-tree-sha1 = "b577d6c6b484f35fc27c1e767dc32458815da0e5"
uuid = "b8f27783-ece8-5eb3-8dc8-9495eed66fee"
version = "1.11.5"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.0+0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.2.1"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "aa532179d4a643d4bd9f328589ca01fa20a0d197"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.1.0"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "a7c3d1da1189a1c2fe843a3bfa04d18d20eb3211"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.Parsers]]
deps = ["Dates", "SnoopPrecompile"]
git-tree-sha1 = "151d91d63d8d6c1a5789ecb7de51547e00480f1b"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.4"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Profile]]
deps = ["Printf"]
uuid = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SeDuMi]]
deps = ["LinearAlgebra", "MATLAB", "MathOptInterface", "SparseArrays"]
git-tree-sha1 = "a3d5e618812270facbc47410dbe424ccfdc8d0d2"
uuid = "a895aaad-f784-5544-9392-bb281339c1b2"
version = "0.4.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SnoopPrecompile]]
deps = ["Preferences"]
git-tree-sha1 = "e760a70afdcd461cf01a575947738d359234665c"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "d75bda01f8c31ebb72df80a46c88b25d1c79c56d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.7"

[[deps.StarAlgebras]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "b40512d1f692da224ce48c9e6722e3e22a71ca7e"
uuid = "0c0c59c1-dc5f-42e9-9a8b-b5dc384a6cd1"
version = "0.1.8"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "6954a456979f23d05085727adb17c4551c19ecd1"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.5.12"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6b7ba252635a5eff6a0b0664a41ee140a1c9e72a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "94f38103c984f89cf77c402f2a68dbd870f8165f"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.11"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.12+3"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"
"""

# ╔═╡ Cell order:
# ╟─dbe625c4-d474-4c3d-bb01-e35912d5dad9
# ╟─fb061d85-1a5f-439a-81d8-78606c005586
# ╟─516346f0-43c8-4c43-af5f-566d56fb0a7a
# ╟─e59da2e3-f60f-4550-be8e-a9035e2f19eb
# ╟─617a09fb-c847-467d-93b9-30a6aed3e058
# ╟─7affa801-20e5-4d2b-9637-dab2e35b3269
# ╠═f7e30187-e485-421f-a9d7-fa14a33c867c
# ╠═d04b139a-08e9-438c-8393-5cf91424300d
# ╟─e562e9e4-a583-11ed-360d-2daaf7526351
# ╟─7dae365e-c45e-48b0-8720-9eb6ecb8dc90
# ╟─cfc39d8f-231e-4a7b-baf2-f94515e98b83
# ╠═f40a2b8c-4dde-4fbb-aeab-4337bc2e7c53
# ╠═cd42c36d-0f38-4169-9ada-34bb401ecd94
# ╟─e83d5691-34d7-4c9e-b94e-fd695e7ed731
# ╠═ba85cb87-2a25-46e4-a42d-f2018e9df332
# ╟─6351206c-8699-4a28-bf3c-6d4335412bd7
# ╟─02514078-89bc-4c1d-9233-8442ece8a4b9
# ╟─003cc82c-bac7-4856-8051-f32bf3f97e49
# ╠═04c3b592-8264-4ec7-a066-01ca33e6f1f7
# ╠═29cd69ae-4d9c-4533-b90b-38659c6ceee5
# ╠═1a05667f-b42d-40b5-89b1-3340f2f25ca9
# ╠═b1443e53-69eb-4417-af10-eb71b8ebf837
# ╟─7fbba67d-580a-4033-8425-f5a2e13b7307
# ╠═50a7ee5b-2f09-4744-ae57-87affde137bf
# ╠═2a66e100-6fd5-4ef7-b132-17811d81e74b
# ╠═153766c0-2c4d-4ad1-b3a8-ca3827932a60
# ╠═1be46ea1-2d3c-45a3-8ddd-52e684f874c0
# ╠═2578cf1d-f9d0-404a-a37b-c893c995bd5f
# ╠═cae53bd5-ba53-4591-9ca7-40278746729c
# ╠═0b62796b-824a-4c86-9396-175d901d0ba2
# ╠═f066b06b-08a4-445b-9b37-a965edb4a6f0
# ╠═d3b9af2a-8677-40fc-b7a5-2ff498ac3256
# ╟─8a2bdbc4-7151-4a1a-b7a0-b57fa5f94903
# ╟─fb357b59-f420-4505-8fdc-75e77201b843
# ╠═c458d180-005f-4c6c-a0a7-8374d7c03b3c
# ╠═0147d4c2-a7d9-45d4-bdb7-8219423b9ebc
# ╠═bd554b57-57f9-42ec-9957-7efb6b3a3538
# ╠═5ffd0d63-82ec-4514-8389-26e2317abc66
# ╠═96157e19-0938-4b69-98cc-e50492f6a518
# ╟─9f050d68-ee64-4eff-8f27-efd8faeca574
# ╠═27e8136b-b62e-4969-824d-73c6059f6447
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
