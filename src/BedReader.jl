import Base: eof, close, position
export BedReader, close, value, eof, advance!, eachposition

type BedReader <: SequenceReader
    bedStream
    done::Bool
    position::Int64
    currRegionStart::Int64
    currRegionEnd::Int64
    currRegionValue::Float64
    contigs::ReferenceContigs
end

function BedReader(bedFileName::ASCIIString, contigs)
    f = open(bedFileName)
    r = BedReader(f, false, 0, 0, 0, 0, contigs)
    advance!(r)
    r
end
function BedReader(f::IO, contigs)
    r = BedReader(f, false, 0, 0, 0, 0, contigs)
    advance!(r)
    r
end
close(reader::BedReader) = close(reader.bedStream)
value(reader::BedReader) = r.currRegionValue
position(reader::BedReader) = reader.position
eof(reader::BedReader) = reader.position == -1

function advance!(r::BedReader)

    while !r.done
        if eof(r.bedStream)
            r.done = true
            r.position = -1
            return
        end

        r.position += 1

        # read until we have the right BED line
        while r.position > r.currRegionEnd && !eof(r.bedStream)
            line = readline(r.bedStream)
            parts = split(line, '\t', limit=6)
            r.currRegionValue = 1.0
            if length(parts) >= 5
                r.currRegionValue = parse(Float64, parts[5])
            end
            offset = contig_offset(r.contigs, parts[1])
            r.currRegionStart = offset + parse(Float64, parts[2]) + 1
            r.currRegionEnd = offset + parse(Float64, parts[3])
        end

        # skip blank areas
        if r.position < r.currRegionStart
            r.position = r.currRegionStart
        end

        # If we found a region then
        if r.position <= r.currRegionEnd
            return
        end
    end
end

# here we want to update the reader
eachposition(r::BedReader) = BedReaderIterator(r)
immutable BedReaderIterator
	reader::BedReader
end
Base.start(it::BedReaderIterator) = it.reader.position
Base.done(it::BedReaderIterator, position) = position == -1
function Base.next(it::BedReaderIterator, position)
	pos = it.reader.position
	advance!(it.reader)
	pos,it.reader.position
end
