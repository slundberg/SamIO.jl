using SamIO
using Base.Test
using GZip

## BamReader

# test forward reads
reader = BamReader("data/small.bam", :forward, SamIO.assembly["GRCh38_alt"])
@test position(reader) == 10544
advance!(reader)
@test position(reader) == 11247
close(reader)

# test reverse reads
reader = BamReader("data/small.bam", :reverse, SamIO.assembly["GRCh38_alt"])
@test position(reader) == 10056
advance!(reader)
@test position(reader) == 10102
close(reader)

# test reading through a whole file
reader = BamReader("data/small.bam", :any, SamIO.assembly["GRCh38_alt"])
lastPos = -1
while !eof(reader)
	lastPos = position(reader)
	advance!(reader)
end
close(reader)

# test reading through a whole file with an iterator
reader = BamReader("data/small.bam", :reverse, SamIO.assembly["GRCh38_alt"])
for pos in eachposition(reader)
end
close(reader)

# test with a different assembly
reader = BamReader("data/small.bam", :forward, SamIO.assembly["GRCh38"])
@test position(reader) == 10544
advance!(reader)
@test position(reader) == 11247
close(reader)

# test with a different assembly
reader = BamReader("data/mm10.bam", :any, SamIO.assembly["mm10"])
@test position(reader) == 3007701
advance!(reader)
@test position(reader) == 3011073
close(reader)




## BedReader

# test first few positions
reader = BedReader("data/small.bed", SamIO.assembly["GRCh38_alt"])
@test position(reader) == 1001
advance!(reader)
@test position(reader) == 1002
advance!(reader)
@test position(reader) == 5001
close(reader)

# test reading through a whole file as a stream
f = open("data/small.bed")
reader = BedReader(f, SamIO.assembly["GRCh38_alt"])
lastPos = -1
while !eof(reader)
	lastPos = position(reader)
	advance!(reader)
end
close(reader)

# test reading through a whole file as a zipped stream
f = GZip.open("data/small.bed.gz")
reader = BedReader(f, SamIO.assembly["GRCh38_alt"])
lastPos = -1
while !eof(reader)
	lastPos = position(reader)

	advance!(reader)
end
close(reader)

# test reading through a whole file with an iterator
reader = BedReader("data/small.bed", SamIO.assembly["GRCh38_alt"])
for pos in eachposition(reader)
end
close(reader)


## BinningMap

reader = BinningMap(BamReader("data/small.bam", :any, SamIO.assembly["GRCh38_alt"]), 1000)
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

reader = BinningMap(BedReader("data/small.bed", SamIO.assembly["GRCh38_alt"]), 1000)
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
