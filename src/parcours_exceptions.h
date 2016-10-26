//
// Created by Arthur Rand on 10/6/16.
//

#ifndef PARCOURS_HMM_GRAPH_EXCEPTIONS_H
#define PARCOURS_HMM_GRAPH_EXCEPTIONS_H


#include <exception>

class ParcoursException: public std::exception {
public:
    char text[1000];

    ParcoursException(char const* fmt, ...) __attribute__((format(printf,2,3))) {
        va_list ap;
        va_start(ap, fmt);
        vsnprintf(text, sizeof text, fmt, ap);
        va_end(ap);
    }


    virtual const char *what() const throw() {
        return text;
    }
};

#endif //PARCOURS_GRAPH_EXCEPTIONS_H
