#include "diagonal.h"

Diagonal::Diagonal(int64_t ay, int64_t L, int64_t R) {
    if ((ay + L) % 2 != 0 || (ay + R) % 2 != 0 || L > R) {
        throw ParcoursException("[Diagonal::Diagonal]: Trying to make Diagonal with invalid "
                                "coordinates xay=%" PRIi64 ", xmyL=%" PRIi64 ", xmyR=%" PRIi64 "\n", ay, L, R);
    }
    xay = ay;
    xmyL = L;
    xmyR = R;
    //assert(xmyL <= xmyR);
    //assert(xay >= 0);
}

int64_t Diagonal::Xay() { return xay; }

int64_t Diagonal::MinXmy() { return xmyL; }

int64_t Diagonal::MaxXmy() { return xmyR; }

int64_t Diagonal::Width() { return (xmyR - xmyL) / 2 + 1; }

void Diagonal::Set(int64_t XaY, int64_t xmyl, int64_t xmyr) {
    xay  = XaY;
    xmyL = xmyl;
    xmyR = xmyr;
}

std::string Diagonal::ToString() {
    char buff[100];
    std::sprintf(buff, "Diagonal, xay: %" PRIi64 " xmyL %" PRIi64 ", xmyR: %" PRIi64 ", width %" PRIi64 "", 
                 Xay(), MinXmy(), MaxXmy(), Width());
    std::string s(buff);
    return s;
}

bool Diagonal::operator == (Diagonal& other) const {
    return xay == other.xay && xmyL == other.xmyL && xmyR == other.xmyR;
}

int64_t Diagonal::XCoordinate(int64_t ay, int64_t my) {
    if ((ay + my) % 2 != 0) throw ParcoursException(
            "[Diagonal::XCoordinate]: Illegal input xay %" PRIi64 "xmy %" PRIi64 "\n", ay, my);
    return (ay + my) / 2;
}

int64_t Diagonal::YCoordinate(int64_t ay, int64_t my) {
    if ((ay + my) % 2 != 0) throw ParcoursException(
            "[Diagonal::YCoordinate]: Illegal input xay %" PRIi64 "xmy %" PRIi64 "\n", ay, my);
    return (ay - my) / 2;
}
