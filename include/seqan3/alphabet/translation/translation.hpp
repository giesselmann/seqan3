// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2017, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2017, Knut Reinert & MPI Molekulare Genetik
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ============================================================================

/*!\file alphabet/translation/translation.hpp
 * \author Sara Hetzel <sara.hetzel AT fu-berlin.de>
 */

#pragma once

#include <tuple>

#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/alphabet/translation/translation_details.hpp>
#include <seqan3/range/concept.hpp>

namespace seqan3
{

/*!\brief Translate one nucleotide triplet into single amino acid (single nucleotide interface)
 *
 * \details
 *
 * Translates single nucleotides into amino acid according to given genetic code.
 *
 * \par Complexity
 *
 * Constant.
 *
 * \par Exceptions
 *
 * Guaranteed not to throw.
*/
template <typename nucl_type, genetic_code gc = genetic_code::CANONICAL>
constexpr aa27 translate_triplet(nucl_type const & n1,
                       nucl_type const & n2,
                       nucl_type const & n3)
{
    return translation_table_dna5_to_aa27_<gc>::VALUE[n1.to_rank()][n2.to_rank()][n3.to_rank()];
}

/*!\brief Translate one nucleotide triplet into single amino acid (tuple interface)
 *
 * \details
 *
 * Translates std::tuple or std::array with 3 nucleotides into amino acid according to given genetic code.
 *
 * \par Complexity
 *
 * Constant.
 *
 * \par Exceptions
 *
 * Guaranteed not to throw.
*/
template <typename tuple_type, genetic_code gc = genetic_code::CANONICAL>
constexpr aa27 translate_triplet(tuple_type const & input_tuple)
{
    static_assert(std::tuple_size<tuple_type>::value == 3);
    return translate_triplet(std::get<0>(input_tuple), std::get<1>(input_tuple), std::get<2>(input_tuple));
}

/*!\brief Translate one nucleotide triplet into single amino acid (range interface)
 *
 * \details
 *
 * Translates range with 3 nucleotides into amino acid according to given genetic code.
 *
 * \par Complexity
 *
 * Constant.
 *
 * \par Exceptions
 *
 * Guaranteed not to throw.
*/
template <typename range_type, genetic_code gc = genetic_code::CANONICAL>
    requires random_access_range_concept<range_type>
constexpr aa27 translate_triplet(range_type && input_range)
{
    static_assert(ranges::size(input_range) == 3);
    return translate_triplet(input_range[0], input_range[1], input_range[2]);
}

}
