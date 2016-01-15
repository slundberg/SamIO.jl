using SamIO
using Base.Test

## BamReader

# test forward reads
reader = BamReader("data/small.bam", :forward, ReferenceContigs_hg38)
@test position(reader) == 10544
advance!(reader)
@test position(reader) == 11247
close(reader)

# test reverse reads
reader = BamReader("data/small.bam", :reverse, ReferenceContigs_hg38)
@test position(reader) == 10056
advance!(reader)
@test position(reader) == 10102
close(reader)

# test reading through a whole file
reader = BamReader("data/small.bam", :any, ReferenceContigs_hg38)
lastPos = -1
while !eof(reader)
	lastPos = position(reader)
	advance!(reader)
end
close(reader)

# test reading through a whole file with an iterator
reader = BamReader("data/small.bam", :reverse, ReferenceContigs_hg38)
for pos in eachposition(reader)
end
close(reader)


## BedReader

# test forward reads
reader = BedReader("data/small.bed", ReferenceContigs_hg38)
@test position(reader) == 1001
advance!(reader)
@test position(reader) == 1002
advance!(reader)
@test position(reader) == 5001
close(reader)

# test reading through a whole file
reader = BedReader("data/small.bed", ReferenceContigs_hg38)
lastPos = -1
while !eof(reader)
	lastPos = position(reader)
	advance!(reader)
end
close(reader)

# test reading through a whole file with an iterator
reader = BedReader("data/small.bed", ReferenceContigs_hg38)
for pos in eachposition(reader)
end
close(reader)


## BinningMap

reader = BinningMap(BamReader("data/small.bam", :any, ReferenceContigs_hg38), 1000)
i = 1
while !eof(reader)
	if i == 1
		@test position(reader) == 11
		@test value(reader) == 4.0
	end
    advance!(reader)
	i += 1
end
close(reader)

reader = BinningMap(BedReader("data/small.bed", ReferenceContigs_hg38), 1000)
i = 1
while !eof(reader)
	if i == 1
		@test position(reader) == 2
		@test value(reader) == 2.0
	end
    advance!(reader)
	i += 1
end
close(reader)
