#include "band.h"
#include "common.h"

template<class T, size_t sn>
Band<T, sn>::Band(AnchorPairs anchors, int64_t lX, int64_t lY, int64_t expansion): index(0) {
    if (lX <= 0 || lY <= 0 || expansion % 2 != 0) {
        throw ParcoursException("[Band::Band] Illegal band construct input: \n"
                                "lX: %" PRIi64 " lY: %" PRIi64 ", expansion %" PRIi64 "\n", lX, lY, expansion);
    }

    diagonals.reserve(lX + lY + 1);
    maxLxLy = lX + lY;

    // Now initialise the diagonals
    uint64_t anchorPairIndex = 0;
    int64_t XaY = 0;
    int64_t pxay = 0, pxmy = 0;
    int64_t nxay = 0, nxmy = 0;
    int64_t xL = 0, yL = 0, xU = 0, yU = 0;
    
    auto set_current_diagonal = [&] (int64_t xL, int64_t yL, int64_t xU, int64_t yU) {
        auto avoid_off_by_one = [XaY] (int64_t xmy) -> int64_t { 
            return (XaY + xmy) % 2 == 0 ? xmy : xmy + 1; 
        };
        
        auto set_boundary = [] (int64_t xmy, int64_t i, int64_t j, int64_t k) -> int64_t {
            if (i < j) {
                return xmy += (2 * (j - i) *  k);
            } else {
                return xmy;
            }
        };
    
        int64_t xmyL = xL - yL;
        int64_t xmyR = xU - yU;
    
        assert(XaY >= xL + yU);
        assert(XaY <= xU + yL);
        
        xmyL = avoid_off_by_one(xmyL);
        xmyR = avoid_off_by_one(xmyR);

        xmyL = set_boundary(xmyL, Diagonal::XCoordinate(XaY, xmyL), xL, 1);
        xmyL = set_boundary(xmyL, yL, Diagonal::YCoordinate(XaY, xmyL), 1);
        xmyR = set_boundary(xmyR, xU, Diagonal::XCoordinate(XaY, xmyR), -1);
        xmyR = set_boundary(xmyR, Diagonal::YCoordinate(XaY, xmyR), yU, -1);
        
        Diagonal d(XaY, xmyL, xmyR);
        diagonals.insert(begin(diagonals) + XaY, d);
    };
    
    auto bound_coordinate = [] (int64_t z, int64_t lZ) -> int64_t { return z < 0 ? 0 : (z > lZ ? lZ : z); };

    while (XaY <= maxLxLy) {
        set_current_diagonal(xL, yL, xU, yU);
        if (nxay == XaY++) {
            //The previous diagonals become the next
            pxay = nxay;
            pxmy = nxmy;

            int64_t x = lX;
            int64_t y = lY;
            
            if (anchorPairIndex < anchors.size()) { 
                std::pair<int64_t, int64_t> anchor = anchors.at(anchorPairIndex);
                anchorPairIndex++;
                x = anchor.first + 1;
                y = anchor.second + 1;
                // Check the anchor pairs
                assert(x > Diagonal::XCoordinate(pxay, pxmy));
                assert(y > Diagonal::YCoordinate(pxay, pxmy));
                assert(x <= lX);
                assert(y <= lY);
                assert(x > 0);
                assert(y > 0);
            }

            nxay = x + y;
            nxmy = x - y;

            xL = bound_coordinate(Diagonal::XCoordinate(pxay, pxmy - expansion), lX);
            yL = bound_coordinate(Diagonal::YCoordinate(nxay, nxmy - expansion), lY);
            xU = bound_coordinate(Diagonal::XCoordinate(nxay, nxmy + expansion), lX);
            yU = bound_coordinate(Diagonal::YCoordinate(pxay, pxmy + expansion), lY);
        }
    }
}

template<class T, size_t sn>
Band<T, sn>::Band(int64_t lx, int64_t ly): Band(EmptyAnchors(), lx, ly, 0) {}

template<class T, size_t sn>
Diagonal Band<T, sn>::Next() {
    int64_t idx_check = index > maxLxLy ? maxLxLy : index;
    if (index <= maxLxLy) index++;
    return diagonals.at(idx_check);
}

template<class T, size_t sn>
Diagonal Band<T, sn>::Previous() {
    if (index > 0) index--;
    return diagonals.at(index);
}

template class Band<double, fiveState>;
