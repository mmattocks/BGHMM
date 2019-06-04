#KMER ORDER/SEQUENCE INTEGER CODING UTILITIES
#higher order DNA alphabet
struct CompoundAlphabet
  symbols::Dict{Kmer,Int64}
end

#holds a DNA sequence, a higher-order index/kmer sequence, the alphabet used to produce the Kmers, and the order number
struct N_Order_ntSequence
  order_no::Int8
  alphabet::CompoundAlphabet
  seqs::Dict{Tuple{Int64,DNASequence},Array{Tuple{Int64,Kmer}}}
end

#build a CompoundAlphabet for DNA of some order_no
function compound_DNA_alphabet(alphabet::Tuple, order_no::Int64)
  symbols = Array{Kmer}(undef, length(alphabet)^(order_no+1))
  tuples = Array{Tuple}
  if order_no > 0
      alphabet_list = [alphabet]
      @inbounds for ord_no in 1:order_no
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
    n_order_dict = Dict{Tuple{Int64,DNASequence},Array{Tuple{Int64,Kmer}}}()
    window = order_no + 1

    @inbounds for (i, seq) in enumerate(seqs)
        kmers = collect(each(Kmer{DNA,window},seq))
        n_order_dict[(i, seq)] = kmers
    end

    return nordseqs = N_Order_ntSequence(order_no, compound_DNA_alphabet(base_tuple, order_no), n_order_dict)
end

#convert tuple kmers to symbol codes
function code_seqs(input::N_Order_ntSequence)
    seqs = collect(keys(input.seqs))
    alphabet = input.alphabet
    output = zeros(Int64, (findmax(length.(collect(values(input.seqs))))[1]+1), length(seqs)) #leave 1 missing value after the longest sequence for indexing sequence length in MS_HMMBase messages
    @inbounds for (i, seq) in enumerate(seqs)
        @inbounds for t in 1:(length(input.seqs[seq]))
            curr_kmer = input.seqs[seq][t][2]
            curr_code = alphabet.symbols[curr_kmer]
            output[t,i]=curr_code
        end
    end
    return output
end
