/// \file htslib/hts_log.h
/// Configuration of log levels.
/* The MIT License
Copyright (C) 2017 Genome Research Ltd.

Author: Anders Kaplan

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#ifndef HTS_LOG_H
#define HTS_LOG_H

#ifdef __cplusplus
extern "C" {
#endif

/// Log levels.
enum htsLogLevel {
    HTS_LOG_OFF,        ///< All logging disabled.
    HTS_LOG_ERROR,      ///< Logging of errors only.
    HTS_LOG_WARNING,    ///< Logging of errors and warnings.
    HTS_LOG_INFO,       ///< Logging of errors, warnings, and normal but significant events.
    HTS_LOG_DEBUG,      ///< Detailed logging enabled. Intended for troubleshooting.
    HTS_LOG_ALL = 100   ///< All logging enabled.
};

/// Sets the selected log level.
void hts_set_log_level(enum htsLogLevel level);

/// Gets the selected log level.
enum htsLogLevel hts_get_log_level();

/// Selected log level.
/*!
 * One of the HTS_LOG_* values. The default is HTS_LOG_INFO.
 * \note Avoid direct use of this variable, use hts_set_log_level and hts_get_log_level instead.
 */
extern int hts_verbose;

#ifdef __cplusplus
}
#endif

#endif // #ifndef HTS_LOG_H
