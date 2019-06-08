function find_position_partition(position::Int64, partition_dict::Dict{String, DataFrame})
    foundPos = false
    position_partition_id = ""
    three_prime_extent = 0
    sample_strand = 0
    for (partition_id, partition) in partition_dict
        hitindex = findfirst(x->x>=position, partition.End)
        if hitindex != nothing
            if position >= partition.Start[hitindex]
                foundPos=true
                position_partition_id = partition_id
                three_prime_extent = partition.End[hitindex] - position
                if partition_id == "exon" || partition_id == "periexonic"
                    strand = partition.Strand[hitindex]
                end
            end
        end
    end

    if foundPos == false
        @error "Position not found among partition coordinates!"
    else
        return position_partition_id, three_prime_extent, strand
    end
end

function get_partition_code_dict(dict_forward::Bool=true)
    if dict_forward == true
        return partition_code_dict = Dict("intergenic"=>1,"periexonic"=>2,"exon"=>3)
    else
        return code_partition_dict = Dict(1=>"intergenic", 2=>"periexonic", 3=>"exon")
    end
end

function mask_sequence_by_partition(maskLength::Int64, seqStart::Int64, scaffold_coords_dict::Dict{String, DataFrame})
    partition_code_dict = get_partition_code_dict()
    seqMask = zeros(Int64, maskLength, 2)
    position = seqStart
    while position <= seqStart+maskLength
        position_partition, partition_extent, position_strand = find_position_partition(position, scaffold_coords_dict)

        partition_code = partition_code_dict[position_partition]
        mask_position = position - seqStart + 1
        seqMask[mask_position:min(maskLength,mask_position + partition_extent),1] .= partition_code
        if position_strand == '+'
            seqMask[mask_position:min(maskLength,mask_position + partition_extent),2] .= 1
        elseif position_strand == '-'
            seqMask[mask_position:min(maskLength,mask_position + partition_extent),2] .= -1
        else
            seqMask[mask_position:min(maskLength,mask_position + partition_extent),2] .= 0
        end

        position += partition_extent + 1
    end

    if any(iszero, seqMask)
        @error "Malformed seqmask!!"
    else
        return seqMask
    end
end

function divide_partitions_by_scaffold(partition_coords_dict::Dict{String, DataFrame})
    scaffold_coords_dict = Dict{Tuple{String, String}, DataFrame}()
    for (partition_id, partition_df) in partition_coords_dict
        for scaffold_subframe in groupby(partition_df, :SeqID)
            scaffold_id = scaffold_subframe.SeqID[1]
            scaffold_df = copy(scaffold_subframe)
            scaffold_coords_dict[(partition_id, scaffold_id)] = scaffold_df
        end
    end
    return scaffold_coords_dict
end

function add_partition_masks(position_df::DataFrame, gff3_path::String, perigenic_pad::Int64=500)
    partitions=["exon", "periexonic", "intergenic"]
    partition_coords_dict = BGHMM.partition_genome_coordinates(gff3_path, perigenic_pad)
    partitioned_scaffolds = divide_partitions_by_scaffold(partition_coords_dict)
    maskcol = Vector{Vector{Int64}}()
    masked_df = DataFrame(Sequence = DNASequence[], Mask = Matrix{Int64}[], Scaffold = String[], Start = Int64[], End = Int64[], Smt = Float64[], Fuzziness = Float64[])

    @showprogress 1 "Masking..." for entry in eachrow(position_df)
        scaffold = entry.Scaffold
        if scaffold != "MT"
            maskLength = length(entry.Sequence)
            seqStart = entry.Start

            scaffold_coords_dict = Dict{String,DataFrame}()
            for partition in partitions
                if haskey(partitioned_scaffolds, (partition, scaffold))
                    scaffold_coords_dict[partition] = partitioned_scaffolds[(partition, scaffold)]
                end
            end

            #@info "Masking from $seqStart on scaffold $scaffold, with $(length(scaffold_coords_dict)) partitions"

            push!(masked_df, [entry.Sequence, mask_sequence_by_partition(maskLength, seqStart, scaffold_coords_dict), entry.Scaffold, entry.Start, entry.End, entry.Smt, entry.Fuzziness])
        end
    end

    return masked_df
end

function get_mask_matrix(mask_matrices::Vector{Matrix{Int64}})
    output = zeros(Int64, (findmax(length.(mask_matrices))[1]+1), length(mask_matrices), 1) #leave 1 "missing value" 0 after the longest sequence to match code_seqs marking of sequence length
    @inbounds for (o, mask) in enumerate(mask_matrices)
        @inbounds for t in 1:(length(mask))
            output[t,o,1]=mask[t,1] #partition code
            output[t,o,2]=mask[t,2] #strand code
        end
    end
    return output
end
