using BGHMM,BioSequences,DataFrames,Distributed,ProgressMeter,Test

include("synthetic_sequence_gen.jl")

#JOB FILEPATHS
#GFF3 feature database, FASTA genome and index paths
Sys.islinux() ? genome =  (@__DIR__) * "/synthetic.fna" : genome = (@__DIR__) * "\\synthetic.fna"
Sys.islinux() ? index =  (@__DIR__) * "/synthetic.fna.fai" : index = (@__DIR__) * "\\synthetic.fna.fai"
Sys.islinux() ? gff =  (@__DIR__) * "/synthetic.gff3" : gff = (@__DIR__) * "\\synthetic.gff3"
Sys.islinux() ? posfasta =  (@__DIR__) * "/syntheticpos.fa" : posfasta = (@__DIR__) * "\\syntheticpos.fa"


!isfile(genome) && print_synthetic_fasta(genome)
!isfile(index) && print_synthetic_index(index)
!isfile(gff) && print_synthetic_gff(gff)
!isfile(posfasta) && print_synthetic_position(posfasta)

@testset "Order coding functions" begin
    test_seqs = [BioSequences.DNASequence("ACGTACGTACGTACGT"),BioSequences.DNASequence("TTTTTTT")]
    target0= [1 4
              2 4
              3 4
              4 4
              1 4
              2 4
              3 4
              4 0
              1 0
              2 0
              3 0
              4 0
              1 0
              2 0
              3 0
              4 0
              0 0]
    target2= [37 64
              58 64
              15 64
              20 64
              37 64
              58 0
              15 0
              20 0
              37 0
              58 0
              15 0
              20 0
              37 0
              58 0
              0 0]
    order0_seqs = BGHMM.get_order_n_seqs(test_seqs,0)
    code0_seqs = BGHMM.code_seqs(order0_seqs)
    @test target0 == code0_seqs
    order2_seqs = BGHMM.get_order_n_seqs(test_seqs,2)
    code2_seqs = BGHMM.code_seqs(order2_seqs)
    @test target2 == code2_seqs
end

@testset "Partition masker functions" begin
    synthetic_seq = BioSequences.DNASequence(generate_synthetic_seq())
    position_start = 501
    position_length=350
    perigenic_pad = 250

    position_df = BGHMM.make_padded_df(posfasta, gff, genome, index, position_length)
    BGHMM.add_partition_masks!(position_df, gff, perigenic_pad)

    @test position_df.SeqID[1]=="1"
    @test length(position_df.PadSeq[1]) == 2*position_length == length(position_df.PadStart[1]:position_df.PadStart[1]+position_length*2-1) == position_df.End[1]-position_df.Start[1]+1+position_length == size(position_df.MaskMatrix[1])[1]
    @test synthetic_seq[position_df.PadStart[1]:position_df.End[1]]==position_df.PadSeq[1]

    mm = position_df.MaskMatrix[1]
    idx=1#pos 151
    @test sum(mm[idx:idx+99,1] .== 1)==length(mm[idx:idx+99,1]); idx+=100#pos 251
    @test sum(mm[idx:idx+259,1] .== 2)==length(mm[idx:idx+259,1]); idx+=260 #pos 511
    @test sum(mm[idx:idx+59,1] .== 3)==length(mm[idx:idx+59,1]); idx+=60#pos 571
    @test sum(mm[idx:idx+29,1] .== 2)==length(mm[idx:idx+29,1]); idx+=30#pos601
    @test sum(mm[idx:idx+59,1] .== 3)==length(mm[idx:idx+59,1]); idx+=60#pos661
    @test sum(mm[idx:idx+39,1] .== 2)==length(mm[idx:idx+39,1]);
end

@testset "Sequence sampler functions" begin
    partitions = 3 #exonic, periexonic, intragenic
    sample_set_length=75
    min_sample_window=5
    max_sample_window=25
    perigenic_pad=250
    syn_intergenic_starts = [1]
    syn_intergenic_ends = [250]
    syn_periexonic_starts = [251,571,661,761,861,961]
    syn_periexonic_ends = [510,600,700,800,900,1000]
    syn_periexonic_strands = ['-','-','-','-','-','-']
    syn_exonic_starts = [511,601,701,801,901]
    syn_exonic_ends = [570,660,760,860,960]
    syn_exonic_strands = ['-','-','-','-','-']
    partition_lengths = Dict("exon"=>5*60,"intergenic"=>250,"periexonic"=>450)

    @info "Testing BGHMM sequence sampler fns..."
    synthetic_seq = BioSequences.DNASequence(generate_synthetic_seq())

    @info "Partitioning synthetic genome coordinates..."
    coordinate_partitions = BGHMM.partition_genome_coordinates(gff, perigenic_pad)

    @test coordinate_partitions["intergenic"].Start == syn_intergenic_starts
    @test coordinate_partitions["intergenic"].End == syn_intergenic_ends

    @test coordinate_partitions["periexonic"].Start == syn_periexonic_starts
    @test coordinate_partitions["periexonic"].End == syn_periexonic_ends
    @test coordinate_partitions["periexonic"].Strand == syn_periexonic_strands

    @test coordinate_partitions["exon"].Start == syn_exonic_starts
    @test coordinate_partitions["exon"].End == syn_exonic_ends
    @test coordinate_partitions["exon"].Strand == syn_exonic_strands

    @info "Checking sampling functions at all synthetic indices..."

    input_test_channel, completed_test_channel = BGHMM.setup_sample_jobs(genome, index, gff, sample_set_length, min_sample_window, max_sample_window, perigenic_pad; deterministic=true)
    while isready(input_test_channel)
        genome_path, genome_index_path, partition_df, partitionid, sample_set_length, sample_window_min, sample_window_max, deterministic = take!(input_test_channel)

        for feature in eachrow(partition_df)
            seqid, scaffold_start, scaffold_end = BGHMM.meta_to_feature_coord(feature.MetaStart, feature.MetaEnd, partition_df)
            @test scaffold_start == feature.Start && scaffold_end == feature.End
        end

        stranded = BGHMM.get_strand_dict()[partitionid]
        @test typeof(stranded) == Bool

        scaffold_sequence_record_dict = BGHMM.build_scaffold_seq_dict(genome, index)
        @test scaffold_sequence_record_dict["1"] == synthetic_seq

        partition_length = partition_df.MetaEnd[end]
        @test partition_length == partition_lengths[partitionid]

        metacoordinate_bitarray = trues(partition_df.MetaEnd[end])

        for bitindex in findall(metacoordinate_bitarray)
            feature_metaStart, feature_metaEnd, strand = BGHMM.get_feature_params_from_metacoord(bitindex, partition_df, stranded)
            @test 1 <= feature_metaStart < feature_metaEnd <= length(metacoordinate_bitarray)
            @test feature_metaStart in partition_df.MetaStart
            @test feature_metaEnd in partition_df.MetaEnd
            if partitionid == "exon" || partitionid == "periexonic"
                @test strand == '-'
            end

            feature_length = length(feature_metaStart:feature_metaEnd)
            window = BGHMM.determine_sample_window(feature_metaStart, feature_metaEnd, bitindex, metacoordinate_bitarray, sample_window_min, sample_window_max) #get an appropriate sampling window around the selected index, given the feature boundaries and params
            @test 1 <= feature_metaStart <= window[1] < window[1]+sample_window_min-1 < window[1]+sample_window_max-1<= window[2] <= feature_metaEnd <= length(metacoordinate_bitarray)
            sample_scaffoldid, sample_scaffold_start, sample_scaffold_end = BGHMM.meta_to_feature_coord(window[1],window[2],partition_df)
            @test sample_scaffoldid == "1"
            @test 1 <= sample_scaffold_start <= sample_scaffold_start+sample_window_min <= sample_scaffold_end <= min(sample_scaffold_start+sample_window_max,1000) <= 1000

            strand == '-' ? target_seq=reverse_complement(synthetic_seq[sample_scaffold_start:sample_scaffold_end]) :
                target_seq = synthetic_seq[sample_scaffold_start:sample_scaffold_end]

            proposal_sequence = BGHMM.fetch_sequence(sample_scaffoldid, scaffold_sequence_record_dict, sample_scaffold_start, sample_scaffold_end, strand; deterministic=deterministic) #get the sequence associated with the sample window

            @test proposal_sequence == target_seq
        end
    end

    @info "Verifying sampling channels..."

    input_sample_channel, completed_sample_channel = BGHMM.setup_sample_jobs(genome, index, gff, sample_set_length, min_sample_window, max_sample_window, perigenic_pad; deterministic=true)
    progress_channel = RemoteChannel(()->Channel{Tuple}(20))
    BGHMM.get_sample_set(input_sample_channel, completed_sample_channel, progress_channel)

    #collect sample dfs by partition id when ready
    sample_record_dfs = Dict{String,DataFrame}()
    collected_counter = 0
    while collected_counter < partitions
        wait(completed_sample_channel)
        partition_id, sample_df = take!(completed_sample_channel)
        sample_record_dfs[partition_id] = sample_df
        collected_counter += 1
    end

    
    synthetic_seq = BioSequences.DNASequence(generate_synthetic_seq())
    for (partid, df) in sample_record_dfs
        for sample in eachrow(df)
            target_seq = synthetic_seq[sample.SampleStart:sample.SampleEnd]
            strand = sample.Strand
            if sample.Strand == '-'
                target_seq = reverse_complement(target_seq)
            end
            @test sample.SampleSequence == target_seq
        end
    end
end
