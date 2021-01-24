// -*- coding:utf-8; mode:c++; mode:auto-fill; fill-column:80; -*-

/// @file      color.hpp
/// @brief     Output stream manipulators to add ANSI console colors.
/// @author    J. Arrieta <Juan.Arrieta@nablazerolabs.com>
/// @date      March 14, 2017
/// @copyright (c) 2017 Nabla Zero Labs
///
/// This code is released under The MIT License
///
/// Copyright 2017 Nabla Zero Labs
///
/// Permission is hereby granted, free of charge, to any person obtaining a copy
/// of this software and associated documentation files (the "Software"), to
/// deal in the Software without restriction, including without limitation the
/// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
/// sell copies of the Software, and to permit persons to whom the Software is
/// furnished to do so, subject to the following conditions:
///
/// The above copyright notice and this permission notice shall be included in
/// all copies or substantial portions of the Software.
///
/// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
/// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
/// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
/// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
/// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
/// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
/// IN THE SOFTWARE.

#pragma once

#include <ostream>
#include <cstdio>
#include <unistd.h>

namespace color {

// ref: https://rosettacode.org/wiki/Check_output_device_is_a_terminal#C.2B.2B
int stdout_isatty() {
	return (isatty(fileno(stdout)));
}

#define MAKE_COLOR_MANIPULATOR(name, code)                               \
	template < typename CharT, typename Traits = std::char_traits<CharT> > \
	inline std::basic_ostream< CharT, Traits >&                            \
	name(std::basic_ostream< CharT, Traits >& os)                          \
	{ return os << code; }

// These color definitions are based on the color scheme used by Git (the
// distributed version control system) as declared in the file
// `git/color.h`. You can add more manipulators as desired.

MAKE_COLOR_MANIPULATOR( normal       , ""           )
MAKE_COLOR_MANIPULATOR( reset        , (stdout_isatty() ? "\033[m"     : ""))
MAKE_COLOR_MANIPULATOR( bold         , (stdout_isatty() ? "\033[1m"    : ""))
MAKE_COLOR_MANIPULATOR( underline    , (stdout_isatty() ? "\033[4m"    : ""))
MAKE_COLOR_MANIPULATOR( ul           , (stdout_isatty() ? "\033[4m"    : ""))
MAKE_COLOR_MANIPULATOR( inverse      , (stdout_isatty() ? "\033[7m"    : ""))
MAKE_COLOR_MANIPULATOR( red          , (stdout_isatty() ? "\033[31m"   : ""))
MAKE_COLOR_MANIPULATOR( green        , (stdout_isatty() ? "\033[32m"   : ""))
MAKE_COLOR_MANIPULATOR( yellow       , (stdout_isatty() ? "\033[33m"   : ""))
MAKE_COLOR_MANIPULATOR( blue         , (stdout_isatty() ? "\033[34m"   : ""))
MAKE_COLOR_MANIPULATOR( magenta      , (stdout_isatty() ? "\033[35m"   : ""))
MAKE_COLOR_MANIPULATOR( cyan         , (stdout_isatty() ? "\033[36m"   : ""))
MAKE_COLOR_MANIPULATOR( bold_red     , (stdout_isatty() ? "\033[1;31m" : ""))
MAKE_COLOR_MANIPULATOR( bold_green   , (stdout_isatty() ? "\033[1;32m" : ""))
MAKE_COLOR_MANIPULATOR( bold_yellow  , (stdout_isatty() ? "\033[1;33m" : ""))
MAKE_COLOR_MANIPULATOR( bold_blue    , (stdout_isatty() ? "\033[1;34m" : ""))
MAKE_COLOR_MANIPULATOR( bold_magenta , (stdout_isatty() ? "\033[1;35m" : ""))
MAKE_COLOR_MANIPULATOR( bold_cyan    , (stdout_isatty() ? "\033[1;36m" : ""))
MAKE_COLOR_MANIPULATOR( bg_red       , (stdout_isatty() ? "\033[41m"   : ""))
MAKE_COLOR_MANIPULATOR( bg_green     , (stdout_isatty() ? "\033[42m"   : ""))
MAKE_COLOR_MANIPULATOR( bg_yellow    , (stdout_isatty() ? "\033[43m"   : ""))
MAKE_COLOR_MANIPULATOR( bg_blue      , (stdout_isatty() ? "\033[44m"   : ""))
MAKE_COLOR_MANIPULATOR( bg_magenta   , (stdout_isatty() ? "\033[45m"   : ""))
MAKE_COLOR_MANIPULATOR( bg_cyan      , (stdout_isatty() ? "\033[46m"   : ""))

} // namespace color
