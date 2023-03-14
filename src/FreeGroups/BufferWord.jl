struct Cell{T} where {T}
    elements::Vector{T}
    length::Integer  # number of _defined_ entries in elements
    previous::Cell{T}
    next::Cell{T}

    function Cell{T}(cellsize::Integer=42) where {T}
        new{T}(Vector{T}(undef, cellsize), 0)
    end
end

function length(cell::Cell)
    return cell.length
end

function push!(cell::Cell{T}, e::T) where {T}
    @assert cell.length + 1 <= length(cell.elements)
    cell.elements[cell.length + 1] = e
    cell.length += 1
end

function pop!(cell::Cell{T}) where {T}
    @assert cell.length - 1 >= 0
    cell.length -= 1
    return cell.elements[cell.length + 1]
end

function pushfirst!(cell::Cell{T}, e::T) where {T}
    @assert cell.length + 1 <= length(cell.elements)
    pushfirst!(cell.elements[start:end - 1], e)
end

function popfirst!(cell::Cell{T}) where {T}
    @assert cell.length - 1 >= 0
    first = cell.elements[1]
    append!(v[2:end], zero(T)) 
    return first
end

function Base.getindex(cell::Cell{T}, i::Integer)
    @assert i <= length(cell) - 1
    return cell.elements[i]
end

function Base.setindex!(cell::Cell{T}, e::T, i::Integer)
    @assert i <= length(cell) - 1
    cell.elements[i] = e
end

hasspace(cell::Cell) = cell.length < length(cell.elements)
isempty(cell::Cell) = cell.length == 0

struct BufferWord{T} <: AbstractWord{T}
    first::Cell{T}
    last::Cell{T}
    cellsize::Integer

    function BufferWord{T}(cellsize::Integer=42) where {T}
        first = Cell{T}(cellsize)
        first.next = first
        first.previous = first
        new{T}(first, first, cellsize)
    end
end

first(bufferWord::BufferWord) = bufferWord.first
last(bufferWord::BufferWord) = bufferWord.last
cellsize(bufferWord::BufferWord) = bufferWord.cellsize

function push!(bufferWord::BufferWord{T}, e::T) where {T}
    last = last(bufferWord)
    if !hasspace(last)
        newCell = Cell{T}(cellsize(bufferWord))
        bufferWord.last.next = newCell
        newCell.previous = last
        last = newCell
        bufferWord.last = last
    end
    push!(last(bufferWord), e)
end

function pop!(bufferWord::BufferWord{T}) where {T}
    if isempty(bufferWord.last)
        dump = bufferWord.last
        bufferWord.last = bufferWord.last.previous
        bufferWord.last.next = bufferWord.first
        dump = nothing
    end
    return pop!(last(bufferWord))
end

function pushfirst!(bufferWord::BufferWord{T}, e::T) where {T}
    first = first(bufferWord)
    if !hasspace(first)
        newCell = Cell{T}(cellsize(bufferWord))
        newCell.next = bufferWord.first
        newCell.previous = bufferWord.last
        bufferWord.first.previous = newCell
        bufferWord.last.next = newCell
        bufferWord.first = newCell
    end
    pushfirst!(first(bufferWord), e)
end

function popfirst!(bufferWord::BufferWord{T}) where {T}
    if isempty(bufferWord.first)
        dump = bufferWord.first
        bufferWord.first = bufferWord.first.next
        bufferWord.first.previous = bufferWord.last
        bufferWord.last.next = bufferWord.first
        dump = nothing
    end
    return popfirst!(first(bufferWord))
end
