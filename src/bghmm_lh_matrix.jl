#function to obtain positional likelihoods for a sequence under a given HMM.
function get_BGHMM_symbol_lh(seq::Matrix{Int64}, hmm::HMM)
    @assert size(seq)[2] == 1 # one sequence at a time only
    symbol_lhs = zeros(length(seq))
    length_mask = [length(seq)-1]
    
    lls = log_likelihoods(hmm, seq) #obtain log likelihoods for sequences and states
    log_α = messages_forwards_log(hmm.π0, hmm.π, lls, length_mask) #get forward messages
    log_β = messages_backwards_log(hmm.π, lls, length_mask) # get backwards messages

    #calculate observation probability and γ weights
    K,Tmaxplus1,Strand = size(lls) #the last T value is the 0 end marker of the longest T

    #transforms to cut down log_ξ, log_γ assignment times
    lls = permutedims(lls, [2,1,3]) # from (K,T) to (T,K)
    log_α = permutedims(log_α, [2,1,3])
    log_β = permutedims(log_β, [2,1,3])

    log_γ = fill(-Inf, Tmaxplus1,K)
    log_pobs = logsumexp(lps.(log_α[1,:], log_β[1,:]))

    @inbounds for i = 1:K, t = 1:length_mask[1]
            log_γ[t,i] = lps(log_α[t,i],log_β[t,i],-log_pobs)
    end

    for t in 1:length_mask[1]
        symbol_lh::Float64 = -Inf #ie log(p=0)
        for k = 1:K #iterate over states
                state_symbol_lh::Float64 = lps(log_γ[t,k], log(hmm.D[k].p[seq[t]])) #state symbol likelihood is the γ weight * the state symbol probability (log implementation)
                symbol_lh = logaddexp(symbol_lh, state_symbol_lh) #sum the probabilities over states
        end
        symbol_lhs[t] = symbol_lh
    end

    return symbol_lhs[1:end-1] #remove trailing index position
end
                #multisequence competent lls fn
                function log_likelihoods(hmm::AbstractHMM{Univariate}, observations::Matrix{Int64})
                    lls = zeros(length(hmm.D), size(observations)...)
                    Threads.@threads for d in 1:length(hmm.D)
                        lls[d,:,:] = logpdf.(hmm.D[d], observations)
                    end
                    return lls
                end

                #Multisequence competent log implementations of forward and backwards algos
                function messages_forwards_log(init_distn::AbstractVector{Float64}, trans_matrix::AbstractMatrix{Float64}, log_likelihoods::Array{Float64,3}, obs_lengths::Vector{Int64})
                    log_alphas = zeros(size(log_likelihoods))
                    log_trans_matrix = log.(trans_matrix)
                    log_alphas[:,1,:] = log.(init_distn) .+ log_likelihoods[:,1,:]
                    Threads.@threads for o in 1:size(log_likelihoods)[3]
                        @inbounds for t in 2:obs_lengths[o], j in 1:size(log_likelihoods)[1]
                            log_alphas[j,t,o]=logsumexp(log_alphas[:,t-1,o] .+ log_trans_matrix[:,j]) + log_likelihoods[j,t,o]
                        end   
                    end
                    return log_alphas
                end

                function messages_backwards_log(trans_matrix::AbstractMatrix{Float64}, log_likelihoods::Array{Float64,3}, obs_lengths::Vector{Int64})
                    log_betas = zeros(size(log_likelihoods))
                    log_trans_matrix = log.(trans_matrix)
                    Threads.@threads for o in 1:size(log_likelihoods)[3]
                        @inbounds for t in obs_lengths[o]-1:-1:1
                            tmp = view(log_betas, :, t+1, o) .+ view(log_likelihoods, :, t+1, o)
                            @inbounds for i in 1:size(log_likelihoods)[1]
                                log_betas[i,t,o] = logsumexp(view(log_trans_matrix, i,:) .+ tmp)
                            end
                        end
                    end
                    return log_betas
                end


#function to calculate BGHMM from an observation set and a dict of BGHMMs
function BGHMM_likelihood_calc(observations::DataFrame, BGHMM_dict::Dict{String,Tuple{HMM, Int64, Float64}}, code_partition_dict = BGHMM.get_partition_code_dict(false))
    lh_matrix_size = ((findmax(length.(collect(values(observations.PadSeq))))[1]), length(observations.PadSeq))
    BGHMM_lh_matrix = zeros(lh_matrix_size) #T, Strand, O

    BGHMM_fragments = fragment_observations_by_BGHMM(observations.PadSeq, observations.MaskMatrix, observations.SeqOffset)

    @showprogress 1 "Writing frags to matrix.." for (jobid, frag) in BGHMM_fragments
        (offset, frag_start, o, partition, strand) = jobid

        partition_BGHMM::HMM = BGHMM_dict[code_partition_dict[partition]][1]
        no_symbols = length(partition_BGHMM.D[1].p)
        order = Int(log(4,no_symbols) - 1)

        order_seq = BGHMM.get_order_n_seqs([frag], order)
        coded_seq = BGHMM.code_seqs(order_seq)
        
        subseq_symbol_lh = get_BGHMM_symbol_lh(coded_seq, partition_BGHMM)

        if strand == -1
            subseq_symbol_lh = reverse(subseq_symbol_lh)
        end #positive, unstranded frags  are inserted as-is
        BGHMM_lh_matrix[offset+frag_start:offset+frag_start+length(frag)-1,o] = subseq_symbol_lh
    end

    return BGHMM_lh_matrix
end

function fragment_observations_by_BGHMM(seqs::Vector{BioSequence{DNAAlphabet{4}}}, masks::Vector{Matrix{Int64}},offsets::Vector{Int64})
    likelihood_jobs = Vector{Tuple{Tuple,BioSequence{DNAAlphabet{4}}}}()
    @showprogress 1 "Fragmenting observations by partition..." for (o, obs_seq) in enumerate(seqs)
        mask = masks[o]
        frags = Vector{BioSequence{DNAAlphabet{4}}}() #container for all subsequences in observation
        frag = BioSequence{DNAAlphabet{4}}()

        frag_end=0
        frag_start = 1

        while frag_start < length(obs_seq) # while we're not at the sequence end
            curr_partition = mask[frag_start,1] #get the partition code of the frag start
            curr_strand = mask[frag_start,2] #get the strand of the frag start

            #JOBID COMPOSED HERE
            jobid = (offsets[o], frag_start, o, curr_partition, curr_strand) #compose an identifying index for this frag

            findnext(!isequal(curr_partition),mask[:,1],frag_start) != nothing ? frag_end = findnext(!isequal(curr_partition),mask[:,1],frag_start) -1 : frag_end = length(obs_seq) #find the next position in the frag that has a different partition mask value from hte current one and set that position-1 to frag end, alternately frag end is end of the overall  sequence 
            frag = obs_seq[frag_start:frag_end] #get the frag bases
            if curr_strand == -1 #if the fragment is reverse stranded
                frag = reverse_complement(frag) #use the reverse complement sequence
            end

            push!(likelihood_jobs,(jobid, frag)) #put the frag in the jobs vec
            frag_start = frag_end + 1 #move on
        end        
    end
    return likelihood_jobs
end