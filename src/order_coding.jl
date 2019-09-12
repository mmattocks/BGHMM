#KMER ORDER/SEQUENCE INTEGER CODING UTILITIES
#higher order DNA alphabet
struct CompoundAlphabet
  symbols::Dict{Kmer,Int64}
end

#holds a DNA sequence, a higher-order index/kmer sequence, the alphabet used to produce the Kmers, and the order number
struct N_Order_ntSequence
  alphabet::CompoundAlphabet
  seq_lengths::Vector{Int64}
  order_kmers::Vector{Vector{Kmer}}
end

#build a CompoundAlphabet for DNA of some order_no
function compound_DNA_alphabet(alphabet::Tuple, order_no::Int64)
  symbols = Array{Kmer}(undef, length(alphabet)^(order_no+1))
  tuples = Array{Tuple}
  if order_no > 0
      alphabet_list = [alphabet]
      for ord_no in 1:order_no
          push!(alphabet_list, alphabet)
      end
      tuples = collect(Iterators.product(alphabet_list...))
  else #0th order compound product of an alphabet is that alphabet
      tuples = collect(Iterators.product(alphabet))
  end

  @inbounds for index in eachindex(tuples)
      tuple_seq = Kmer(DNASequence(collect(tuples[index])))
      symbols[index] = tuple_seq
  end

  code_dict=Dict{Kmer,Int64}()
  @inbounds for i in 1:length(symbols)
      code_dict[symbols[i]]=i
  end

  return CompoundAlphabet(code_dict)
end

#from a vector of DNASequences, get
function get_order_n_seqs(seqs::Vector{DNASequence}, order_no::Int64, base_tuple::Tuple=ACGT)
    kmer_vecs = Vector{Vector{Kmer}}()
    length_vec = Vector{Int64}()
    window = order_no + 1

    for seq in seqs
        kmer_vec = Vector{Kmer}()
        @inbounds for (i, kmer) in collect(each(Kmer{DNA,window},seq))
            push!(kmer_vec, kmer)
        end

        push!(kmer_vecs, kmer_vec)
        push!(length_vec, length(seq))
    end

    return nordseqs = N_Order_ntSequence(compound_DNA_alphabet(base_tuple, order_no), length_vec, kmer_vecs)
end

#convert tuple kmers to symbol codes
function code_seqs(input::N_Order_ntSequence, offsets::Array{Int64}=[0 for i in 1:length(input.order_kmers)])
    alphabet = input.alphabet
    output = zeros(Int64, (maximum([length(seq) for seq in input.order_kmers])+1), length(input.order_kmers)) #leave 1 missing value after the longest sequence for indexing sequence length in CLHMM messages
    for (i, seq) in enumerate(input.order_kmers)
        for t in 1:length(seq)
            curr_kmer = input.order_kmers[i][t]
            curr_code = alphabet.symbols[curr_kmer]
            output[t+offsets[i],i]=curr_code
        end
    end
    return output
end
