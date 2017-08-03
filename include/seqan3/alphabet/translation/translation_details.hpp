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

/*!\file translation_details.hpp
 * \ingroup translation
 * \author Sara Hetzel <sara.hetzel AT fu-berlin.de>
 * \brief Contains translation details for nucleotide to aminoacid conversion.
 */

#pragma once

#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/alphabet/aminoacid/all.hpp>

namespace seqan3
{
enum struct genetic_code : uint8_t
{
    CANONICAL=1,
    VERT_MITOCHONDRIAL,
    YEAST_MITOCHONDRIAL,
    MOLD_MITOCHONDRIAL,
    INVERT_MITOCHONDRIAL,
    CILIATE,
    FLATWORM_MITOCHONDRIAL = 9,
    EUPLOTID,
    PROKARYOTE,
    ALT_YEAST,
    ASCIDIAN_MITOCHONDRIAL,
    ALT_FLATWORM_MITOCHONDRIAL,
    BLEPHARISMA,
    CHLOROPHYCEAN_MITOCHONDRIAL,
    TREMATODE_MITOCHONDRIAL = 21,
    SCENEDESMUS_MITOCHONDRIAL,
    THRAUSTOCHYTRIUM_MITOCHONDRIAL,
    PTEROBRANCHIA_MITOCHONDRIAL,
    GRACILIBACTERIA
};

template <genetic_code, typename void_type = void>
struct translation_table_dna5_to_aa27_
{
    static aa27 const VALUE[5][5][5];
};


template <typename void_type>
struct translation_table_dna5_to_aa27_<genetic_code::CANONICAL, void_type>
{
    static aa27 const VALUE[5][5][5];
};

template <typename void_type>
aa27 const translation_table_dna5_to_aa27_<genetic_code::CANONICAL, void_type>::VALUE[5][5][5] =
{
    { // a??
        { aa27::K,          aa27::N, aa27::K,          aa27::N, aa27::X }, // aa?
        { aa27::T,          aa27::T, aa27::T,          aa27::T, aa27::T }, // ac?
        { aa27::R,          aa27::S, aa27::R,          aa27::S, aa27::X }, // ag?
        { aa27::I,          aa27::I, aa27::M,          aa27::I, aa27::X }, // au?
        { aa27::X,          aa27::X, aa27::X,          aa27::X, aa27::X }  // an?
    }, { // c??
        { aa27::Q,          aa27::H, aa27::Q,          aa27::H, aa27::X }, // ca?
        { aa27::P,          aa27::P, aa27::P,          aa27::P, aa27::P }, // cc?
        { aa27::R,          aa27::R, aa27::R,          aa27::R, aa27::R }, // cg?
        { aa27::L,          aa27::L, aa27::L,          aa27::L, aa27::L }, // cu?
        { aa27::X,          aa27::X, aa27::X,          aa27::X, aa27::X }  // cn?
    }, { // g??
        { aa27::E,          aa27::D, aa27::E,          aa27::D, aa27::X }, // ga?
        { aa27::A,          aa27::A, aa27::A,          aa27::A, aa27::A }, // gc?
        { aa27::G,          aa27::G, aa27::G,          aa27::G, aa27::G }, // gg?
        { aa27::V,          aa27::V, aa27::V,          aa27::V, aa27::V }, // gu?
        { aa27::X,          aa27::X, aa27::X,          aa27::X, aa27::X }  // gn?
    }, { // u??
        { aa27::TERMINATOR, aa27::Y, aa27::TERMINATOR, aa27::Y, aa27::X }, // ua?
        { aa27::S,          aa27::S, aa27::S,          aa27::S, aa27::S }, // uc?
        { aa27::TERMINATOR, aa27::C, aa27::W,          aa27::C, aa27::X }, // ug?
        { aa27::L,          aa27::F, aa27::L,          aa27::F, aa27::X }, // uu?
        { aa27::X,          aa27::X, aa27::X,          aa27::X, aa27::X }  // un?
    }, { // n??
        { aa27::X,          aa27::X, aa27::X,          aa27::X, aa27::X }, // na?
        { aa27::X,          aa27::X, aa27::X,          aa27::X, aa27::X }, // nc?
        { aa27::X,          aa27::X, aa27::X,          aa27::X, aa27::X }, // ng?
        { aa27::X,          aa27::X, aa27::X,          aa27::X, aa27::X }, // nu?
        { aa27::X,          aa27::X, aa27::X,          aa27::X, aa27::X }  // nn?
    }
};

/*
template <genetic_code, typename void_type = void>
struct translation_table_nucl16_to_aa27_
{
    static aa27 const VALUE[16][16][16];
};


template <typename void_type>
struct translation_table_nucl16_to_aa27_<genetic_code::CANONICAL, void_type>
{
    static aa27 const VALUE[16][16][16];
};

template <typename void_type>
char const translation_table_nucl16_to_aa27_<genetic_code<CANONICAL>, void_type>::VALUE[16][16][16] =
{
    { // a??
        // a,               b,       c,       d,       g,                h,       k,       m,       n,       r,                s,       t,       u,       v,       w,       y
        { aa27::K,          aa27::X, aa27::N, aa27::X, aa27::K,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::K,          aa27::X, aa27::N, aa27::N, aa27::X, aa27::X, aa27::N }, // aa?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // ab?
        { aa27::T,          aa27::T, aa27::T, aa27::T, aa27::T,          aa27::T, aa27::T, aa27::T, aa27::T, aa27::T,          aa27::T, aa27::T, aa27::T, aa27::T, aa27::T, aa27::T }, // ac?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // ad?
        { aa27::R,          aa27::X, aa27::S, aa27::X, aa27::R,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::R,          aa27::X, aa27::S, aa27::S, aa27::X, aa27::X, aa27::S }, // ag?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // ah?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // ak?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // am?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // an?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // ar?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // as?
        { aa27::I,          aa27::X, aa27::I, aa27::X, aa27::M,          aa27::I, aa27::X, aa27::I, aa27::X, aa27::X,          aa27::X, aa27::I, aa27::I, aa27::X, aa27::I, aa27::I }, // at?
        { aa27::I,          aa27::X, aa27::I, aa27::X, aa27::M,          aa27::I, aa27::X, aa27::I, aa27::X, aa27::X,          aa27::X, aa27::I, aa27::I, aa27::X, aa27::I, aa27::I }, // au?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // av?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // aw?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }  // ay?
   }, { // b??
        // a,               b,       c,       d,       g,                h,       k,       m,       n,       r,                s,       t,       u,       v,       w,       y
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // ba?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // bb?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // bc?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // bd?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // bg?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // bh?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // bk?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // bm?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // bn?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // br?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // bs?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // bt?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // bu?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // bv?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // bw?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }  // by?
    }, { // c??
        // a,               b,       c,       d,       g,                h,       k,       m,       n,       r,                s,       t,       u,       v,       w,       y
        { aa27::Q,          aa27::X, aa27::H, aa27::X, aa27::Q,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::Q,          aa27::X, aa27::H, aa27::H, aa27::X, aa27::X, aa27::H }, // ca?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // cb?
        { aa27::P,          aa27::P, aa27::P, aa27::P, aa27::P,          aa27::P, aa27::P, aa27::P, aa27::P, aa27::P,          aa27::P, aa27::P, aa27::P, aa27::P, aa27::P, aa27::P }, // cc?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // cd?
        { aa27::R,          aa27::R, aa27::R, aa27::R, aa27::R,          aa27::R, aa27::R, aa27::R, aa27::R, aa27::R,          aa27::R, aa27::R, aa27::R, aa27::R, aa27::R, aa27::R }, // cg?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // ch?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // ck?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // cm?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // cn?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // cr?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // cs?
        { aa27::L,          aa27::L, aa27::L, aa27::L, aa27::L,          aa27::L, aa27::L, aa27::L, aa27::L, aa27::L,          aa27::L, aa27::L, aa27::L, aa27::L, aa27::L, aa27::L }, // ct?
        { aa27::L,          aa27::L, aa27::L, aa27::L, aa27::L,          aa27::L, aa27::L, aa27::L, aa27::L, aa27::L,          aa27::L, aa27::L, aa27::L, aa27::L, aa27::L, aa27::L }, // cu?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // cv?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // cw?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }  // cy?
    }, { // d??
        // a,               b,       c,       d,       g,                h,       k,       m,       n,       r,                s,       t,       u,       v,       w,       y
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // da?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // db?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // dc?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // dd?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // dg?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // dh?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // dk?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // dm?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // dn?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // dr?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // ds?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // dt?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // du?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // dv?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // dw?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }  // dy?

    }, { // g??
        // a,               b,       c,       d,       g,                h,       k,       m,       n,       r,                s,       t,       u,       v,       w,       y
        { aa27::E,          aa27::X, aa27::D, aa27::X, aa27::E,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::E,          aa27::X, aa27::D, aa27::D, aa27::X, aa27::X, aa27::D }, // ga?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // gb?
        { aa27::A,          aa27::A, aa27::A, aa27::A, aa27::A,          aa27::A, aa27::A, aa27::A, aa27::A, aa27::A,          aa27::A, aa27::A, aa27::A, aa27::A, aa27::A, aa27::A }, // gc?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // gd?
        { aa27::G,          aa27::G, aa27::G, aa27::G, aa27::G,          aa27::G, aa27::G, aa27::G, aa27::G, aa27::G,          aa27::G, aa27::G, aa27::G, aa27::G, aa27::G, aa27::G }, // gg?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // gh?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // gk?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // gm?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // gn?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // gr?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // gs?
        { aa27::V,          aa27::V, aa27::V, aa27::V, aa27::V,          aa27::V, aa27::V, aa27::V, aa27::V, aa27::V,          aa27::V, aa27::V, aa27::V, aa27::V, aa27::V, aa27::V }, // gt?
        { aa27::V,          aa27::V, aa27::V, aa27::V, aa27::V,          aa27::V, aa27::V, aa27::V, aa27::V, aa27::V,          aa27::V, aa27::V, aa27::V, aa27::V, aa27::V, aa27::V }, // gu?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // gv?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // gw?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }  // gy?
    }, { // h??
        // a,               b,       c,       d,       g,                h,       k,       m,       n,       r,                s,       t,       u,       v,       w,       y
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // ha?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // hb?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // hc?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // hd?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // hg?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // hh?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // hk?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // hm?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // hn?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // hr?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // hs?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // ht?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // hu?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // hv?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // hw?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }  // hy?

    }, { // k??
        // a,               b,       c,       d,       g,                h,       k,       m,       n,       r,                s,       t,       u,       v,       w,       y
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // ka?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // kb?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // kc?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // kd?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // kg?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // kh?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // kk?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // km?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // kn?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // kr?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // ks?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // kt?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // ku?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // kv?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // kw?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }  // ky?

    }, { // m??
        // a,               b,       c,       d,       g,                h,       k,       m,       n,       r,                s,       t,       u,       v,       w,       y
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // ma?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // mb?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // mc?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // md?
        { aa27::R,          aa27::X, aa27::X, aa27::X, aa27::R,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::R,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // mg?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // mh?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // mk?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // mm?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // mn?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // mr?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // ms?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // mt?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // mu?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // mv?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // mw?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }  // my?


    }, { // n??
        // a,               b,       c,       d,       g,                h,       k,       m,       n,       r,                s,       t,       u,       v,       w,       y
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // nb?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // nc?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // nd?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // ng?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // nh?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // nk?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // nm?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // nn?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // nr?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // ns?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // nt?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // nu?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // nv?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // nw?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }  // ny?
    }, { // r??
        // a,               b,       c,       d,       g,                h,       k,       m,       n,       r,                s,       t,       u,       v,       w,       y
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // rb?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // rc?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // rd?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // rg?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // rh?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // rk?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // rm?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // rn?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // rr?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // rs?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // rt?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // ru?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // rv?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // rw?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }  // ry?
    }, { // s??
        // a,               b,       c,       d,       g,                h,       k,       m,       n,       r,                s,       t,       u,       v,       w,       y
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // sb?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // sc?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // sd?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // sg?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // sh?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // sk?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // sm?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // sn?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // sr?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // ss?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // st?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // su?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // sv?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // sw?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }  // sy?
    }, { // t??
        // a,               b,       c,       d,       g,                h,       k,       m,       n,       r,                s,       t,       u,       v,       w,       y
        { aa27::TERMINATOR, aa27::X, aa27::Y, aa27::X, aa27::TERMINATOR, aa27::X, aa27::X, aa27::X, aa27::X, aa27::TERMINATOR, aa27::X, aa27::Y, aa27::Y, aa27::X, aa27::X, aa27::Y }, // ta?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // tb?
        { aa27::S,          aa27::S, aa27::S, aa27::S, aa27::S,          aa27::S, aa27::S, aa27::S, aa27::S, aa27::S,          aa27::S, aa27::S, aa27::S, aa27::S, aa27::S, aa27::S }, // tc?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // td?
        { aa27::TERMINATOR, aa27::X, aa27::C, aa27::X, aa27::W,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::C, aa27::C, aa27::X, aa27::C }, // tg?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // th?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // tk?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // tm?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // tn?
        { aa27::TERMINATOR, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // tr?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // ts?
        { aa27::L,          aa27::X, aa27::F, aa27::X, aa27::X,          aa27::L, aa27::X, aa27::X, aa27::X, aa27::L,          aa27::X, aa27::F, aa27::F, aa27::X, aa27::X, aa27::F }, // tt?
        { aa27::L,          aa27::X, aa27::F, aa27::X, aa27::X,          aa27::L, aa27::X, aa27::X, aa27::X, aa27::L,          aa27::X, aa27::F, aa27::F, aa27::X, aa27::X, aa27::F }, // tu?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // tv?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // tw?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }  // ty?
    }, { // u??
        // a,               b,       c,       d,       g,                h,       k,       m,       n,       r,                s,       t,       u,       v,       w,       y
        { aa27::TERMINATOR, aa27::X, aa27::Y, aa27::X, aa27::TERMINATOR, aa27::X, aa27::X, aa27::X, aa27::X, aa27::TERMINATOR, aa27::X, aa27::Y, aa27::Y, aa27::X, aa27::X, aa27::Y }, // ua?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // ub?
        { aa27::S,          aa27::S, aa27::S, aa27::S, aa27::S,          aa27::S, aa27::S, aa27::S, aa27::S, aa27::S,          aa27::S, aa27::S, aa27::S, aa27::S, aa27::S, aa27::S }, // uc?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // ud?
        { aa27::TERMINATOR, aa27::X, aa27::C, aa27::X, aa27::W,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::C, aa27::C, aa27::X, aa27::C }, // ug?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // uh?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // uk?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // um?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // un?
        { aa27::TERMINATOR, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // ur?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // us?
        { aa27::L,          aa27::X, aa27::F, aa27::X, aa27::X,          aa27::L, aa27::X, aa27::X, aa27::X, aa27::L,          aa27::X, aa27::F, aa27::F, aa27::X, aa27::X, aa27::F }, // ut?
        { aa27::L,          aa27::X, aa27::F, aa27::X, aa27::X,          aa27::L, aa27::X, aa27::X, aa27::X, aa27::L,          aa27::X, aa27::F, aa27::F, aa27::X, aa27::X, aa27::F }, // uu?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // uv?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // uw?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }  // uy?
    }, { // v??
        // a,               b,       c,       d,       g,                h,       k,       m,       n,       r,                s,       t,       u,       v,       w,       y
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // va?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // vb?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // vc?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // vd?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // vg?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // vh?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // vk?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // vm?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // vn?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // vr?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // vs?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // vt?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // vu?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // vv?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // vw?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }  // vy?
    }, { // w??
        // a,               b,       c,       d,       g,                h,       k,       m,       n,       r,                s,       t,       u,       v,       w,       y
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // wa?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // wb?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // wc?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // wd?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // wg?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // wh?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // wk?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // wm?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // wn?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // wr?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // ws?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // wt?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // wu?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // wv?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // ww?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }  // wy?
    }, { // y??
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // ya?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // yb?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // yc?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // yd?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // yg?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // yh?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // yk?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // ym?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // yn?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // yr?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // ys?
        { aa27::L,          aa27::X, aa27::X, aa27::X, aa27::L,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::L,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // yt?
        { aa27::L,          aa27::X, aa27::X, aa27::X, aa27::L,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::L,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // yu?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // yv?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // yw?
        { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }  // yy?
    }
};*/

}
