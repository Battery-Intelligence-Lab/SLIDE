/*
 * Span.hpp
 *
 * A small class for returning non-owning view of some storages.
 *  Created on: 13 Feb 2022
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#pragma once

#include <array>
#include <vector>

namespace slide {
//!< std::span is decided to be used now.

//!< template <typename Tdata>
//!< class Span
//!< {
//!<     Tdata *T_begin{nullptr};
//!<     size_t n{};

//!< public:
//!<     Span() = default;

//!<     template <typename Tcontainer>
//!<     Span(Tcontainer &data) : T_begin(&data[0]), n(data.size()) {}

//!<     template <typename Tcontainer>
//!<     Span(Tcontainer &data, size_t n) : T_begin(&data[0]), n(n) {}

//!<     Span(Tcontainer &data, size_t n) : T_begin(&data[0]), n(n) {}

//!<     template <typename Titer>
//!<     Span(Titer begin, Titer end) : T_begin(&(*begin)), T_end(&(*end))
//!<     {
//!<     }

//!<     Tdata *begin() { return T_begin; }
//!<     Tdata *end() { return T_end; }
//!< };
} // namespace slide