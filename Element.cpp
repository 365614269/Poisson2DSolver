#include "Element.h"

Element::Element(int elemIndex) {
    this->elemIndex = elemIndex;

    for (int i = 0; i < Nlb; i++) {
        this->nodes[i] = Node(elemIndex, i);
    }
}

Element::Element() {}