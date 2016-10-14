//
// Created by Arthur Rand on 10/13/16.
//

#include "common.h"

void st_uglyf(const char *string, ...) {
    va_list ap;
    va_start(ap, string);
    vfprintf(stderr, string, ap);
    va_end(ap);
}