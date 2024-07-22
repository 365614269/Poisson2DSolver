#include "include/Element.h"

Element::Element(int elemIndex) {
    this->elemIndex = elemIndex;

    for (int i = 0; i < Nlb; i++) {
        this->nodes[i] = Node(elemIndex, i);
    }
}

Element::Element() {}

int Element::getElemIndex() {
    return this->elemIndex;
}

Node& Element::getNode(int localNodeIndex) {
    return this->nodes[localNodeIndex];
}