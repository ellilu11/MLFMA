#include "node.h"

using enum Dir;

/* getNeighborGeqSize(dir)
 * Find the unique neighbor node (if it exists) of equal or greater size
 * along a given direction
 * dir : direction
 */
std::shared_ptr<Node> Node::getNeighborGeqSize(const Dir dir) const {
    if (isRoot()) return nullptr;

    std::shared_ptr<Node> nbor;

    switch (dir) {
        case U:
            if (branchIdx < 4)
                return base->branches[branchIdx+4];
            nbor = base->getNeighborGeqSize(dir);
            if (nbor == nullptr) return nbor;
            if (nbor->isNodeType<Leaf>()) return nbor;
            return nbor->branches[branchIdx-4];
            break;

        case N:
            if (branchIdx == 0 || branchIdx == 1 || branchIdx == 4 || branchIdx == 5)
                return base->branches[branchIdx+2];
            nbor = base->getNeighborGeqSize(dir);
            if (nbor == nullptr) return nbor;
            if (nbor->isNodeType<Leaf>()) return nbor;
            return nbor->branches[branchIdx-2];
            break;

        case E:
            if (branchIdx == 0 || branchIdx == 2 || branchIdx == 4 || branchIdx == 6)
                return base->branches[branchIdx+1];
            nbor = base->getNeighborGeqSize(dir);
            if (nbor == nullptr) return nbor;
            if (nbor->isNodeType<Leaf>()) return nbor;
            return nbor->branches[branchIdx-1];
            break;

        case W:
            if (branchIdx == 1 || branchIdx == 3 || branchIdx == 5 || branchIdx == 7)
                return base->branches[branchIdx-1];
            nbor = base->getNeighborGeqSize(dir);
            if (nbor == nullptr) return nbor;
            if (nbor->isNodeType<Leaf>()) return nbor;
            return nbor->branches[branchIdx+1];
            break;

        case S:
            if (branchIdx == 2 || branchIdx == 3 || branchIdx == 6 || branchIdx == 7)
                return base->branches[branchIdx-2];
            nbor = base->getNeighborGeqSize(dir);
            if (nbor == nullptr) return nbor;
            if (nbor->isNodeType<Leaf>()) return nbor;
            return nbor->branches[branchIdx+2];
            break;

        case D:
            if (branchIdx >= 4)
                return base->branches[branchIdx-4];
            nbor = base->getNeighborGeqSize(dir);
            if (nbor == nullptr) return nbor;
            if (nbor->isNodeType<Leaf>()) return nbor;
            return nbor->branches[branchIdx+4];
            break;

        case UN:
            if (branchIdx == 0 || branchIdx == 1)
                return base->branches[branchIdx+6];
            if (branchIdx == 6 || branchIdx == 7) {
                nbor = base->getNeighborGeqSize(dir);
                if (nbor == nullptr) return nbor;
                if (nbor->isNodeType<Leaf>()) return nbor;
                return nbor->branches[branchIdx-6];
            } 
            nbor = (branchIdx == 4 || branchIdx == 5) ?
                base->getNeighborGeqSize(U) :
                base->getNeighborGeqSize(N);
            if (nbor == nullptr) return nbor;
            // do not double count neighbor that will be found along a cardinal direction
            if (nbor->isNodeType<Leaf>()) return nullptr;
            return (branchIdx == 2 || branchIdx == 4) ?
                nbor->branches[6-branchIdx] :
                nbor->branches[8-branchIdx];
            break;

        case DS:
            if (branchIdx == 6 || branchIdx == 7)
                return base->branches[branchIdx-6];
            if (branchIdx == 0 || branchIdx == 1) {
                nbor = base->getNeighborGeqSize(dir);
                if (nbor == nullptr) return nbor;
                if (nbor->isNodeType<Leaf>()) return nbor;
                return nbor->branches[branchIdx+6];
            } 
            nbor = (branchIdx == 2 || branchIdx == 3) ?
                base->getNeighborGeqSize(D) :
                base->getNeighborGeqSize(S);
            if (nbor == nullptr) return nbor;
            if (nbor->isNodeType<Leaf>()) return nullptr;
            return (branchIdx == 2 || branchIdx == 4) ?
                nbor->branches[6-branchIdx] :
                nbor->branches[8-branchIdx];
            break;

        case US:
            if (branchIdx == 2 || branchIdx == 3)
                return base->branches[branchIdx+2];
            if (branchIdx == 4 || branchIdx == 5) {
                nbor = base->getNeighborGeqSize(dir);
                if (nbor == nullptr) return nbor;
                if (nbor->isNodeType<Leaf>()) return nbor;
                return nbor->branches[branchIdx-2];
            }
            nbor = (branchIdx == 6 || branchIdx == 7) ?
                base->getNeighborGeqSize(U) :
                base->getNeighborGeqSize(S);
            if (nbor == nullptr) return nbor;
            if (nbor->isNodeType<Leaf>()) return nullptr;
            return (branchIdx == 0 || branchIdx == 6) ?
                nbor->branches[6-branchIdx] :
                nbor->branches[8-branchIdx];
            break;

        case DN:
            if (branchIdx == 4 || branchIdx == 5)
                return base->branches[branchIdx-2];
            if (branchIdx == 2 || branchIdx == 3) {
            nbor = base->getNeighborGeqSize(dir);
                if (nbor == nullptr) return nbor;
                if (nbor->isNodeType<Leaf>()) return nbor;
            return nbor->branches[branchIdx+2];
            }
            nbor = (branchIdx == 0 || branchIdx == 1) ?
                base->getNeighborGeqSize(D) :
                base->getNeighborGeqSize(N);
            if (nbor == nullptr) return nbor;
            if (nbor->isNodeType<Leaf>()) return nullptr;
            return (branchIdx == 0 || branchIdx == 6) ?
                nbor->branches[6-branchIdx] :
                nbor->branches[8-branchIdx];
            break;

        case UE:
            if (branchIdx == 0 || branchIdx == 2)
                return base->branches[branchIdx+5];
            if (branchIdx == 5 || branchIdx == 7) {
                nbor = base->getNeighborGeqSize(dir);
                if (nbor == nullptr) return nbor;
                if (nbor->isNodeType<Leaf>()) return nbor;
                return nbor->branches[branchIdx-5];
            }
            nbor = (branchIdx == 4 || branchIdx == 6) ?
                base->getNeighborGeqSize(U) :
                base->getNeighborGeqSize(E);
            if (nbor == nullptr) return nbor;
            if (nbor->isNodeType<Leaf>()) return nullptr;
            return (branchIdx == 1 || branchIdx == 4) ?
                nbor->branches[5-branchIdx] :
                nbor->branches[9-branchIdx];
            break;

        case DW:
            if (branchIdx == 5 || branchIdx == 7)
                return base->branches[branchIdx-5];
            if (branchIdx == 0 || branchIdx == 2) {
                nbor = base->getNeighborGeqSize(dir);
                if (nbor == nullptr) return nbor;
                if (nbor->isNodeType<Leaf>()) return nbor;
                return nbor->branches[branchIdx+5];
            }
            nbor = (branchIdx == 1 || branchIdx == 3) ?
                base->getNeighborGeqSize(D) :
                base->getNeighborGeqSize(W);
            if (nbor == nullptr) return nbor;
            if (nbor->isNodeType<Leaf>()) return nullptr;
            return (branchIdx == 1 || branchIdx == 4) ?
                nbor->branches[5-branchIdx] :
                nbor->branches[9-branchIdx];
            break;

        case UW:
            if (branchIdx == 1 || branchIdx == 3)
                return base->branches[branchIdx+3];
            if (branchIdx == 4 || branchIdx == 6) {
                nbor = base->getNeighborGeqSize(dir);
                if (nbor == nullptr) return nbor;
                if (nbor->isNodeType<Leaf>()) return nbor;
                return nbor->branches[branchIdx-3];
            }
            nbor = (branchIdx == 5 || branchIdx == 7) ?
                base->getNeighborGeqSize(U) :
                base->getNeighborGeqSize(W);
            if (nbor == nullptr) return nbor;
            if (nbor->isNodeType<Leaf>()) return nullptr;
            return (branchIdx == 0 || branchIdx == 5) ?
                nbor->branches[5-branchIdx] :
                nbor->branches[9-branchIdx];
            break;

        case DE:
            if (branchIdx == 4 || branchIdx == 6)
                return base->branches[branchIdx-3];
            if (branchIdx == 1 || branchIdx == 3) {
                nbor = base->getNeighborGeqSize(dir);
                if (nbor == nullptr) return nbor;
                if (nbor->isNodeType<Leaf>()) return nbor;
                return nbor->branches[branchIdx+3];
            }
            nbor = (branchIdx == 0 || branchIdx == 2) ?
                base->getNeighborGeqSize(D) :
                base->getNeighborGeqSize(E);
            if (nbor == nullptr) return nbor;
            if (nbor->isNodeType<Leaf>()) return nullptr;
            return (branchIdx == 0 || branchIdx == 5) ?
                nbor->branches[5-branchIdx] :
                nbor->branches[9-branchIdx];
            break;

        case NE:
            if (branchIdx == 0 || branchIdx == 4)
                return base->branches[branchIdx+3];
            if (branchIdx == 3 || branchIdx == 7) {
                nbor = base->getNeighborGeqSize(dir);
                if (nbor == nullptr) return nbor;
                if (nbor->isNodeType<Leaf>()) return nbor;
                return nbor->branches[branchIdx-3];
            }
            nbor = (branchIdx == 2 || branchIdx == 6) ?
                base->getNeighborGeqSize(N) :
                base->getNeighborGeqSize(E);
            if (nbor == nullptr) return nbor;
            if (nbor->isNodeType<Leaf>()) return nullptr;
            return (branchIdx == 1 || branchIdx == 2) ?
                nbor->branches[3-branchIdx] :
                nbor->branches[11-branchIdx];
            break;

        case SW:
            if (branchIdx == 3 || branchIdx == 7)
                return base->branches[branchIdx-3];
            if (branchIdx == 0 || branchIdx == 4) {
                nbor = base->getNeighborGeqSize(dir);
                if (nbor == nullptr) return nbor;
                if (nbor->isNodeType<Leaf>()) return nbor;
                return nbor->branches[branchIdx+3];
            }
            nbor = (branchIdx == 1 || branchIdx == 5) ?
                base->getNeighborGeqSize(S) :
                base->getNeighborGeqSize(W);
            if (nbor == nullptr) return nbor;
            if (nbor->isNodeType<Leaf>()) return nullptr;
            return (branchIdx == 1 || branchIdx == 2) ?
                nbor->branches[3-branchIdx] :
                nbor->branches[11-branchIdx];
            break;

        case NW:
            if (branchIdx == 1 || branchIdx == 5)
                return base->branches[branchIdx+1];
            if (branchIdx == 2 || branchIdx == 6) {
                nbor = base->getNeighborGeqSize(dir);
                if (nbor == nullptr) return nbor;
                if (nbor->isNodeType<Leaf>()) return nbor;
                return nbor->branches[branchIdx-1];
            } 
            nbor = (branchIdx == 3 || branchIdx == 7) ?
                base->getNeighborGeqSize(N) :
                base->getNeighborGeqSize(W);
            if (nbor == nullptr) return nbor;
            if (nbor->isNodeType<Leaf>()) return nullptr;
            return (branchIdx == 0 || branchIdx == 3) ?
                nbor->branches[3-branchIdx] :
                nbor->branches[11-branchIdx];
            break;

        case SE:
            if (branchIdx == 2 || branchIdx == 6)
                return base->branches[branchIdx-1];
            if (branchIdx == 1 || branchIdx == 5) {
                nbor = base->getNeighborGeqSize(dir);
                if (nbor == nullptr) return nbor;
                if (nbor->isNodeType<Leaf>()) return nbor;
                return nbor->branches[branchIdx+1];
            }
            nbor = (branchIdx == 0 || branchIdx == 4) ?
                base->getNeighborGeqSize(S) :
                base->getNeighborGeqSize(E);
            if (nbor == nullptr) return nbor;
            if (nbor->isNodeType<Leaf>()) return nullptr;
            return (branchIdx == 0 || branchIdx == 3) ?
                nbor->branches[3-branchIdx] :
                nbor->branches[11-branchIdx];
            break;

        case UNE:
            if (branchIdx == 0)
                return base->branches[7];
            if (branchIdx == 7) {
                nbor = base->getNeighborGeqSize(dir);
                if (nbor == nullptr) return nbor;
                if (nbor->isNodeType<Leaf>()) return nbor;
                return nbor->branches[0];
            }
            nbor = [this] {
                switch (branchIdx) {
                    case 1: return base->getNeighborGeqSize(E);
                    case 2: return base->getNeighborGeqSize(N);
                    case 3: return base->getNeighborGeqSize(NE);
                    case 4: return base->getNeighborGeqSize(U);
                    case 5: return base->getNeighborGeqSize(UE);
                    case 6: return base->getNeighborGeqSize(UN);
                } } ();
            if (nbor == nullptr) return nbor;
            if (nbor->isNodeType<Leaf>()) return nullptr;
            return nbor->branches[7-branchIdx];
            break;

        case DSW:
            if (branchIdx == 7)
                return base->branches[0];
            if (branchIdx == 0) {
                nbor = base->getNeighborGeqSize(dir);
                if (nbor == nullptr) return nbor;
                if (nbor->isNodeType<Leaf>()) return nbor;
                return nbor->branches[7];
            } 
            nbor = [this] {
                switch (branchIdx) {
                    case 1: return base->getNeighborGeqSize(DS);
                    case 2: return base->getNeighborGeqSize(DW);
                    case 3: return base->getNeighborGeqSize(D);
                    case 4: return base->getNeighborGeqSize(SW);
                    case 5: return base->getNeighborGeqSize(S);
                    case 6: return base->getNeighborGeqSize(W);
                } } ();
            if (nbor == nullptr) return nbor;
            if (nbor->isNodeType<Leaf>()) return nullptr;
            return nbor->branches[7-branchIdx];
            break;

        case UNW:
            if (branchIdx == 1)
                return base->branches[6];
            if (branchIdx == 6) {
                nbor = base->getNeighborGeqSize(dir);
                if (nbor == nullptr) return nbor;
                if (nbor->isNodeType<Leaf>()) return nbor;
                return nbor->branches[1];
            } 
            nbor = [this] {
                switch (branchIdx) {
                    case 0: return base->getNeighborGeqSize(W);
                    case 2: return base->getNeighborGeqSize(NW);
                    case 3: return base->getNeighborGeqSize(N);
                    case 4: return base->getNeighborGeqSize(UW);
                    case 5: return base->getNeighborGeqSize(U);
                    case 7: return base->getNeighborGeqSize(UN);
                } } ();
            if (nbor == nullptr) return nbor;
            if (nbor->isNodeType<Leaf>()) return nullptr;
            return nbor->branches[7-branchIdx];
            break;

        case DSE:
            if (branchIdx == 6)
                return base->branches[1];
            if (branchIdx == 1) {
                nbor = base->getNeighborGeqSize(dir);
                if (nbor == nullptr) return nbor;
                if (nbor->isNodeType<Leaf>()) return nbor;
                return nbor->branches[6];
            }
            nbor = [this] {
                switch (branchIdx) {
                    case 0: return base->getNeighborGeqSize(DS);
                    case 2: return base->getNeighborGeqSize(D);
                    case 3: return base->getNeighborGeqSize(DE);
                    case 4: return base->getNeighborGeqSize(S);
                    case 5: return base->getNeighborGeqSize(SE);
                    case 7: return base->getNeighborGeqSize(E);
                } } ();
                if (nbor == nullptr) return nbor;
                if (nbor->isNodeType<Leaf>()) return nullptr;
                return nbor->branches[7-branchIdx];
            break;

        case USE:
            if (branchIdx == 2)
                return base->branches[5];
            if (branchIdx == 5) {
                nbor = base->getNeighborGeqSize(dir);
                if (nbor == nullptr) return nbor;
                if (nbor->isNodeType<Leaf>()) return nbor;
                return nbor->branches[2];
            }
            nbor = [this] {
                switch (branchIdx) {
                    case 0: return base->getNeighborGeqSize(S);
                    case 1: return base->getNeighborGeqSize(SE);
                    case 3: return base->getNeighborGeqSize(E);
                    case 4: return base->getNeighborGeqSize(US);
                    case 6: return base->getNeighborGeqSize(U);
                    case 7: return base->getNeighborGeqSize(UE);
                } } ();
            if (nbor == nullptr) return nbor;
            if (nbor->isNodeType<Leaf>()) return nullptr;
            return nbor->branches[7-branchIdx];
            break;

        case DNW:
            if (branchIdx == 5)
                return base->branches[2];
            if (branchIdx == 2) {
                nbor = base->getNeighborGeqSize(dir);
                if (nbor == nullptr) return nbor;
                if (nbor->isNodeType<Leaf>()) return nbor;
                return nbor->branches[5];
            }
            nbor = [this] {
                switch (branchIdx) {
                    case 0: return base->getNeighborGeqSize(DW);
                    case 1: return base->getNeighborGeqSize(D);
                    case 3: return base->getNeighborGeqSize(DN);
                    case 4: return base->getNeighborGeqSize(W);
                    case 6: return base->getNeighborGeqSize(NW);
                    case 7: return base->getNeighborGeqSize(N);
                } } ();
            if (nbor == nullptr) return nbor;
            if (nbor->isNodeType<Leaf>()) return nullptr;
            return nbor->branches[7-branchIdx];
            break;

        case USW:
            if (branchIdx == 3)
                return base->branches[4];
            if (branchIdx == 4) {
                nbor = base->getNeighborGeqSize(dir);
                if (nbor == nullptr) return nbor;
                if (nbor->isNodeType<Leaf>()) return nbor;
                return nbor->branches[3];
            }
            nbor = [this] {
                switch (branchIdx) {
                    case 0: return base->getNeighborGeqSize(SW);
                    case 1: return base->getNeighborGeqSize(S);
                    case 2: return base->getNeighborGeqSize(W);
                    case 5: return base->getNeighborGeqSize(US);
                    case 6: return base->getNeighborGeqSize(UW);
                    case 7: return base->getNeighborGeqSize(U);
                } } ();
            if (nbor == nullptr) return nbor;
            if (nbor->isNodeType<Leaf>()) return nullptr;
            return nbor->branches[7-branchIdx];
            break;

        case DNE:
            if (branchIdx == 4)
                return base->branches[3];
            if (branchIdx == 3) {
                nbor = base->getNeighborGeqSize(dir);
                if (nbor == nullptr) return nbor;
                if (nbor->isNodeType<Leaf>()) return nbor;
                return nbor->branches[4];
            } 
            nbor = [this] {
                switch (branchIdx) {
                    case 0: return base->getNeighborGeqSize(D);
                    case 1: return base->getNeighborGeqSize(DE);
                    case 2: return base->getNeighborGeqSize(DN);
                    case 5: return base->getNeighborGeqSize(E);
                    case 6: return base->getNeighborGeqSize(N);
                    case 7: return base->getNeighborGeqSize(NE);
                } } ();
            if (nbor == nullptr) return nbor;
            if (nbor->isNodeType<Leaf>()) return nullptr;
            return nbor->branches[7-branchIdx];
            break;
    }
}

/* getNeighborsLeqSize(nborGeqSize,dir)
 * Find all neighbor leaves contained in nborGeqSize
 * nborGeqSize : neighbor node of greater or equal size
 * dir         : direction (must be consistent with direction of nborGeqSize)
 */
NodeVec Node::getNeighborsLeqSize(
    const std::shared_ptr<Node> nborGeqSize, const Dir dir) const
{
    NodeVec nbors{};
    std::queue<std::shared_ptr<Node>> queue{ { nborGeqSize } };
    if (isRoot()) return nbors;

    while (!queue.empty()) {
        auto nbor = queue.front();

        if (nbor->isNodeType<Leaf>())
            nbors.push_back(nbor);
        else {
            switch (dir) {
                case U:
                    queue.push(nbor->branches[0]);
                    queue.push(nbor->branches[1]);
                    queue.push(nbor->branches[2]);
                    queue.push(nbor->branches[3]);
                    break;

                case D:
                    queue.push(nbor->branches[4]);
                    queue.push(nbor->branches[5]);
                    queue.push(nbor->branches[6]);
                    queue.push(nbor->branches[7]);
                    break;

                case N:
                    queue.push(nbor->branches[0]);
                    queue.push(nbor->branches[1]);
                    queue.push(nbor->branches[4]);
                    queue.push(nbor->branches[5]);
                    break;

                case S:
                    queue.push(nbor->branches[2]);
                    queue.push(nbor->branches[3]);
                    queue.push(nbor->branches[6]);
                    queue.push(nbor->branches[7]);
                    break;

                case E:
                    queue.push(nbor->branches[0]);
                    queue.push(nbor->branches[2]);
                    queue.push(nbor->branches[4]);
                    queue.push(nbor->branches[6]);
                    break;

                case W:
                    queue.push(nbor->branches[1]);
                    queue.push(nbor->branches[3]);
                    queue.push(nbor->branches[5]);
                    queue.push(nbor->branches[7]);
                    break;

                case UN:
                    queue.push(nbor->branches[0]);
                    queue.push(nbor->branches[1]);
                    break;

                case DS:
                    queue.push(nbor->branches[6]);
                    queue.push(nbor->branches[7]);
                    break;

                case US:
                    queue.push(nbor->branches[2]);
                    queue.push(nbor->branches[3]);
                    break;

                case DN:
                    queue.push(nbor->branches[4]);
                    queue.push(nbor->branches[5]);
                    break;

                case UE:
                    queue.push(nbor->branches[0]);
                    queue.push(nbor->branches[2]);
                    break;

                case DW:
                    queue.push(nbor->branches[5]);
                    queue.push(nbor->branches[7]);
                    break;

                case UW:
                    queue.push(nbor->branches[1]);
                    queue.push(nbor->branches[3]);
                    break;

                case DE:
                    queue.push(nbor->branches[4]);
                    queue.push(nbor->branches[6]);
                    break;

                case NE:
                    queue.push(nbor->branches[0]);
                    queue.push(nbor->branches[4]);
                    break;

                case SW:
                    queue.push(nbor->branches[3]);
                    queue.push(nbor->branches[7]);
                    break;

                case NW:
                    queue.push(nbor->branches[1]);
                    queue.push(nbor->branches[5]);
                    break;

                case SE:
                    queue.push(nbor->branches[2]);
                    queue.push(nbor->branches[6]);
                    break;

                case UNE:
                    queue.push(nbor->branches[0]);
                    break;

                case DSW:
                    queue.push(nbor->branches[7]);
                    break;

                case UNW:
                    queue.push(nbor->branches[1]);
                    break;

                case DSE:
                    queue.push(nbor->branches[6]);
                    break;

                case USE:
                    queue.push(nbor->branches[2]);
                    break;

                case DNW:
                    queue.push(nbor->branches[5]);
                    break;

                case USW:
                    queue.push(nbor->branches[3]);
                    break;

                case DNE:
                    queue.push(nbor->branches[4]);
                    break;
            }
        }

        queue.pop();
    }

    return nbors;
}