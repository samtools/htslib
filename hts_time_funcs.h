/*  hts_time_funcs.h -- Implementations of non-standard time functions

    Copyright (C) 2022 Genome Research Ltd.

    Author: Rob Davies <rmd@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

/*
  This mainly exists because timegm() is not a standard function, and so
  Cannot be used in portable code.  Unfortunately the standard one (mktime)
  always takes the local timezone into accout so doing a UTC conversion
  with it involves changing the TZ environment variable, which is rather
  messy and not likely to go well with threaded code.

  The code here is a much simplified version of the BSD timegm() implementation.
  It currently rejects dates before 1970, avoiding problems with -ve time_t.
  It also works strictly in UTC, so doesn't have to worry about tm_isdst
  which makes the calculation much easier.

  Some of this is derived from BSD sources, for example
  https://github.com/NetBSD/src/blob/trunk/lib/libc/time/localtime.c
  which state:

  ** This file is in the public domain, so clarified as of
  ** 1996-06-05 by Arthur David Olson.

  Non-derived code is copyright as above.
*/

#include <stdint.h>
#include <limits.h>
#include <errno.h>
#include <time.h>

static inline int hts_time_normalise(int *tens, int *units, int base) {
    if (*units < 0 || *units >= base) {
        int delta = *units >= 0 ? *units / base : (-1 - (-1 - *units) / base);
        int64_t tmp = (int64_t) (*tens) + delta;
        if (tmp < INT_MIN || tmp > INT_MAX) return 1;
        *tens = tmp;
        *units -= delta * base;
    }
    return 0;
}

static inline int hts_year_is_leap(int64_t year) {
    return ((year % 4 == 0) && (year % 100 != 0)) || (year % 400 == 0);
}

// Number of leap years to start of year
// Only works for year >= 1.
static inline int64_t hts_leaps_to_year_start(int64_t year) {
    --year;
    return year / 4 - year / 100 + year / 400;
}

static inline int hts_time_normalise_tm(struct tm *t)
{
    const int days_per_mon[2][12] = {
        { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 },
        { 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 }
    };
    const int year_days[2] = { 365, 366 };
    int overflow = 0;
    int64_t year;

    if (t->tm_sec > 62) {
        overflow |= hts_time_normalise(&t->tm_min, &t->tm_sec, 60);
    }
    overflow |= hts_time_normalise(&t->tm_hour, &t->tm_min,  60);
    overflow |= hts_time_normalise(&t->tm_mday, &t->tm_hour, 24);
    overflow |= hts_time_normalise(&t->tm_year, &t->tm_mon,  12);
    if (overflow)
        return 1;

    year = (int64_t) t->tm_year + 1900LL;
    while (t->tm_mday <= 0) {
        --year;
        t->tm_mday += year_days[hts_year_is_leap(year + (1 < t->tm_mon))];
    }
    while (t->tm_mday > 366) {
        t->tm_mday -= year_days[hts_year_is_leap(year + (1 < t->tm_mon))];
        ++year;
    }
    for (;;) {
        int mdays = days_per_mon[hts_year_is_leap(year)][t->tm_mon];
        if (t->tm_mday <= mdays)
            break;
        t->tm_mday -= mdays;
        t->tm_mon++;
        if (t->tm_mon >= 12) {
            year++;
            t->tm_mon = 0;
        }
    }
    year -= 1900;
    if (year != t->tm_year) {
        if (year < INT_MIN || year > INT_MAX)
            return 1;
        t->tm_year = year;
    }
    return 0;
}

/**
 *  Convert broken-down time to an equivalent time_t value
 *  @param target  Target broken-down time structure
 *  @return Equivalent time_t value on success; -1 on failure
 *
 *  This function first normalises the time in @p target so that the
 *  structure members are in the valid range.  It then calculates the
 *  number of seconds (ignoring leap seconds) between midnight Jan 1st 1970
 *  and the target date.
 *
 *  If @p target is outside the range that can be represented in a time_t,
 *  or tm_year is less than 70 (which would return a negative value) then
 *  it returns -1 and sets errno to EOVERFLOW.
 */

static inline time_t hts_time_gm(struct tm *target)
{
    int month_start[2][12] = {
        { 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334 },
        { 0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335 }
    };
    int years_from_epoch, leaps, days;
    int64_t secs;

    if (hts_time_normalise_tm(target) != 0)
        goto overflow;

    if (target->tm_year < 70)
        goto overflow;

    years_from_epoch = target->tm_year - 70;
    leaps = (hts_leaps_to_year_start(target->tm_year + 1900)
        - hts_leaps_to_year_start(1970));
    days = ((365 * (years_from_epoch - leaps) + 366 * leaps)
        + month_start[hts_year_is_leap(target->tm_year + 1900)][target->tm_mon]
        + target->tm_mday - 1);
    secs = ((int64_t) days * 86400LL
        + target->tm_hour * 3600
        + target->tm_min * 60
        + target->tm_sec);
    if (sizeof(time_t) < 8 && secs > INT_MAX)
        goto overflow;

    return (time_t) secs;

 overflow:
    errno = EOVERFLOW;
    return (time_t) -1;
}
