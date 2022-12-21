abstract type AbstractFPGroup <: Group end
abstract type GroupElement end

struct FreeGroup <: AbstractFPGroup #(?)
	A::Alphabet
	gens::Vector{...}
	# ... ?
end

struct FPGroupElement{W<:AbstractWord, G<:AbstractFPGroup} <: GroupElement
	word::W
	parent::G
	# ....
end